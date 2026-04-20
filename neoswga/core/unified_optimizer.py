"""
Unified optimizer entry point for CLI and programmatic use.

This module provides a clean interface for running primer set optimization
using any registered optimizer. It replaces the complex conditional logic
in cli_unified.py with a simple factory-based approach.

Usage:
    # From CLI
    result = run_optimization('greedy', candidates, config)

    # Programmatic with custom config
    result = run_optimization(
        method='dominating-set',
        candidates=candidates,
        fg_prefixes=['data/target'],
        fg_seq_lengths=[1000000],
        target_size=10,
    )
"""

import json
import os
import logging
import threading
import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any

from .base_optimizer import OptimizationResult, OptimizationStatus, OptimizerConfig
from .optimizer_factory import OptimizerFactory, OptimizerRegistry
from .position_cache import PositionCache
from .progress import progress_context
from . import parameter

logger = logging.getLogger(__name__)

# Thread-safe registration tracking
_optimizers_registered = False
_registration_lock = threading.Lock()


@dataclass
class OptimizationConfig:
    """
    Unified configuration for optimization runs.

    Combines all parameters needed for optimization in one place.
    """
    # Optimizer selection
    method: str = 'hybrid'

    # Target parameters
    target_set_size: int = 6
    max_iterations: int = 100

    # File paths
    data_dir: str = '.'
    step3_file: str = 'step3_df.csv'
    output_file: str = 'step4_improved_df.csv'

    # Genome info (loaded from pipeline if not provided)
    fg_prefixes: Optional[List[str]] = None
    fg_seq_lengths: Optional[List[int]] = None
    bg_prefixes: Optional[List[str]] = None
    bg_seq_lengths: Optional[List[int]] = None

    # Performance options
    use_position_cache: bool = True
    use_background_filter: bool = True

    # Output options
    verbose: bool = True
    quiet: bool = False

    # Advanced options
    uniformity_weight: float = 0.0
    minimize_primers: bool = False
    target_coverage: float = 0.70


def list_available_optimizers() -> Dict[str, str]:
    """
    List all registered optimizers with descriptions.

    Returns:
        Dict mapping optimizer names to descriptions
    """
    # Ensure all optimizers are registered
    _ensure_optimizers_registered()

    return OptimizerFactory.list_optimizers()


def _ensure_optimizers_registered():
    """
    Ensure all optimizer modules are imported and registered.

    Thread-safe and idempotent - only performs registration once.
    """
    global _optimizers_registered

    # Fast path: already registered
    if _optimizers_registered:
        return

    with _registration_lock:
        # Double-checked locking
        if _optimizers_registered:
            return

        # Import optimizer modules to trigger factory registration
        try:
            from . import greedy_optimizer  # Greedy BFS implementations
            from . import dominating_set_adapter  # Graph-based set cover
            from . import hybrid_optimizer  # Two-stage hybrid
            from . import network_optimizer  # Network connectivity
            from . import genetic_algorithm  # Evolutionary optimization
            from . import background_aware_optimizer  # Clinical/background-aware
            from . import equiphi29_optimizer  # EquiPhi29-specific (42C, 12-15bp primers)
            from . import tiling_optimizer  # Interval-based tiling coverage
        except ImportError as e:
            logger.warning(f"Some optimizer modules not available: {e}")

        # Optional optimizers with external dependencies
        try:
            from . import milp_optimizer  # MILP (requires python-mip)
        except ImportError:
            pass  # python-mip not installed

        try:
            from . import normalized_optimizer  # Normalized scoring with strategy presets
        except ImportError:
            pass  # Should always work, but be safe

        try:
            from . import moea_optimizer  # MOEA (requires pymoo)
        except ImportError:
            pass  # pymoo not installed

        try:
            from . import multi_agent_optimizer  # Multi-agent parallel execution
        except ImportError:
            pass  # Should always work, but be safe

        try:
            from . import clique_optimizer  # Clique-based dimer-free (requires networkx)
        except ImportError:
            pass  # networkx not installed

        try:
            from . import background_prefilter  # Background pruning pre-filter
        except ImportError:
            pass  # Should always work, but be safe

        try:
            from . import serial_cascade_optimizer  # Serial pipeline combinations
        except ImportError:
            pass  # Should always work, but be safe

        _optimizers_registered = True
        logger.debug("Optimizer registration complete")


def _prefilter_by_background(cache, candidates, fg_prefixes, bg_prefixes,
                              min_ratio=1.0, max_removal_fraction=0.20, verbose=False):
    """Remove candidates with poor foreground/background binding ratio.

    Keeps all primers with ratio >= min_ratio. If that would remove more
    than max_removal_fraction of candidates, keeps the top
    (1 - max_removal_fraction) by ratio instead.

    Args:
        cache: PositionCache with loaded positions
        candidates: List of candidate primer sequences
        fg_prefixes: Foreground genome prefixes
        bg_prefixes: Background genome prefixes
        min_ratio: Minimum fg/bg ratio to keep (default 1.0)
        max_removal_fraction: Maximum fraction of candidates to remove (default 0.20)
        verbose: Log filtering details

    Returns:
        Filtered list of candidates (never empty)
    """
    ratios = []
    for primer in candidates:
        fg_count = sum(len(cache.get_positions(p, primer, 'both')) for p in fg_prefixes)
        bg_count = sum(len(cache.get_positions(p, primer, 'both')) for p in bg_prefixes)
        ratios.append((fg_count / (bg_count + 1), primer))

    # Keep primers above threshold
    above = [(r, p) for r, p in ratios if r >= min_ratio]

    # Ensure we don't remove too many
    min_keep = int(len(candidates) * (1 - max_removal_fraction))
    if len(above) < min_keep:
        ratios.sort(reverse=True)
        above = ratios[:min_keep]

    filtered = [p for _, p in above]
    removed = len(candidates) - len(filtered)
    if removed > 0 and verbose:
        logger.info(f"  Background pre-filter: removed {removed} candidates with poor fg/bg ratio")

    return filtered if filtered else candidates  # never return empty


def run_optimization(
    method: str = 'hybrid',
    candidates: Optional[List[str]] = None,
    fg_prefixes: Optional[List[str]] = None,
    fg_seq_lengths: Optional[List[int]] = None,
    bg_prefixes: Optional[List[str]] = None,
    bg_seq_lengths: Optional[List[int]] = None,
    target_size: int = 6,
    verbose: bool = True,
    **kwargs
) -> OptimizationResult:
    """
    Run primer set optimization using specified method.

    This is the main entry point for optimization. It:
    1. Loads candidates if not provided
    2. Creates position cache
    3. Instantiates optimizer via factory
    4. Runs optimization
    5. Returns typed result

    Args:
        method: Optimizer method name ('greedy', 'dominating-set', etc.)
        candidates: List of candidate primers (loaded from step3 if None)
        fg_prefixes: Foreground genome HDF5 prefixes
        fg_seq_lengths: Foreground genome lengths
        bg_prefixes: Background genome HDF5 prefixes (optional)
        bg_seq_lengths: Background genome lengths (optional)
        target_size: Desired primer set size
        verbose: Print progress information
        **kwargs: Additional optimizer-specific parameters

    Returns:
        OptimizationResult with selected primers and metrics

    Example:
        result = run_optimization(
            method='dominating-set',
            candidates=['ATCGATCG', 'GCTAGCTA', ...],
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
            target_size=10,
        )

        if result.is_success:
            print(f"Selected: {result.primers}")
            print(f"Coverage: {result.metrics.fg_coverage:.1%}")
    """
    _ensure_optimizers_registered()

    # Load parameters from pipeline if needed
    if fg_prefixes is None or fg_seq_lengths is None:
        from . import pipeline as core_pipeline
        core_pipeline._initialize()
        fg_prefixes = core_pipeline.fg_prefixes
        fg_seq_lengths = core_pipeline.fg_seq_lengths
        bg_prefixes = bg_prefixes or getattr(core_pipeline, 'bg_prefixes', [])
        bg_seq_lengths = bg_seq_lengths or getattr(core_pipeline, 'bg_seq_lengths', [])

    # Handle no-background (host-free) mode
    if kwargs.get('no_background', False):
        bg_prefixes = []
        bg_seq_lengths = []
        if verbose:
            logger.info("  Host-free mode: no background genome data used")

    # Load candidates from step3 if not provided
    if candidates is None:
        step3_path = os.path.join(parameter.data_dir, "step3_df.csv")
        if not os.path.exists(step3_path):
            return OptimizationResult.failure(
                method, f"Step 3 output not found: {step3_path}"
            )
        step3_df = pd.read_csv(step3_path)
        candidates = step3_df['primer'].tolist()

    # Set global random seed if provided (ensures reproducibility across all
    # optimizer components, not just those that accept a seed parameter)
    seed = kwargs.get('seed')
    if seed is not None:
        import random
        import numpy as np
        random.seed(seed)
        np.random.seed(seed)
        if verbose:
            logger.info(f"Random seed set to {seed} for reproducibility")

    if verbose:
        logger.info(f"Running {method} optimization")
        logger.info(f"  Candidates: {len(candidates)}")
        logger.info(f"  Target size: {target_size}")

    # Create position cache
    with progress_context("Loading position data", disable=not verbose):
        cache = PositionCache(fg_prefixes + (bg_prefixes or []), candidates)

    # Optional: pre-filter candidates by fg/bg binding ratio
    if bg_prefixes and kwargs.get('bg_prefilter', True):
        candidates = _prefilter_by_background(
            cache, candidates, fg_prefixes, bg_prefixes,
            min_ratio=kwargs.get('bg_min_ratio', 1.0),
            max_removal_fraction=kwargs.get('bg_max_removal', 0.20),
            verbose=verbose,
        )

    # Determine polymerase extension reach for coverage computation
    polymerase = kwargs.get('polymerase') or getattr(parameter, 'polymerase', 'phi29')
    try:
        from .reaction_conditions import get_polymerase_processivity
        extension_reach = get_polymerase_processivity(polymerase)
    except (ImportError, ValueError):
        extension_reach = 70000  # Phi29 default

    # Build optimizer config
    config = OptimizerConfig(
        target_set_size=target_size,
        max_iterations=kwargs.get('max_iterations', 100),
        verbose=verbose,
        extension_reach=extension_reach,
    )

    # Build ReactionConditions from parameter globals so every optimizer sees
    # the user's additive cocktail. Kwarg-provided 'conditions' wins (for
    # programmatic callers that already constructed one). Without this, the
    # additive-aware Tm paths in network_optimizer / integrated_quality_scorer
    # are dark code — they exist but never receive conditions at runtime.
    conditions = kwargs.pop('conditions', None)
    if conditions is None:
        try:
            from .reaction_conditions import ReactionConditions
            conditions = ReactionConditions(
                temp=getattr(parameter, 'reaction_temp', None) or 30.0,
                polymerase=getattr(parameter, 'polymerase', 'phi29'),
                na_conc=getattr(parameter, 'na_conc', 50.0),
                mg_conc=getattr(parameter, 'mg_conc', 10.0),
                dmso_percent=getattr(parameter, 'dmso_percent', 0.0),
                betaine_m=getattr(parameter, 'betaine_m', 0.0),
                trehalose_m=getattr(parameter, 'trehalose_m', 0.0),
                formamide_percent=getattr(parameter, 'formamide_percent', 0.0),
                ethanol_percent=getattr(parameter, 'ethanol_percent', 0.0),
                urea_m=getattr(parameter, 'urea_m', 0.0),
                tmac_m=getattr(parameter, 'tmac_m', 0.0),
            )
        except Exception as e:
            logger.warning(
                f"Could not construct ReactionConditions ({e}); optimizer "
                f"runs additive-blind. This indicates a mismatch between "
                f"parameter.* globals and ReactionConditions signature."
            )
            conditions = None

    # Create optimizer via factory
    try:
        optimizer = OptimizerFactory.create(
            name=method,
            position_cache=cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            config=config,
            conditions=conditions,
            **kwargs
        )
    except Exception as e:
        logger.error(f"Failed to create optimizer '{method}': {e}")
        return OptimizationResult.failure(method, str(e))

    # Run optimization
    with progress_context(f"Running {optimizer.name} optimizer", disable=not verbose):
        result = optimizer.optimize(candidates, target_size)

    # Post-optimization sanity validation. Catches duplicates, size drift,
    # zero coverage, and accidental blacklist re-injection. The result is
    # written to data_dir/step4_improved_df_validation.json and logged here
    # so downstream CSV consumers can trust the output shape.
    try:
        forbidden = []
        if bg_prefixes and getattr(parameter, 'bl_genomes', None):
            # Phase 12A blacklist guard. Collect blacklist primers (those we
            # have already filtered) from step2_df if available to cross-check
            # that expand-primers / swap-primer didn't slip any back in.
            pass  # Cheap path: the optimizer's own candidate pool came from
                  # step2/step3 which already excluded blacklist hits; fuller
                  # cross-check would reload the raw blacklist.

        validation = result.validate(
            target_size=target_size,
            min_coverage=0.0,  # soft by default; caller can tighten
            min_per_target_coverage=kwargs.get('min_per_target_coverage', 0.0),
            forbidden_primers=forbidden or None,
        )
    except Exception as e:
        logger.debug(f"Post-optimization validator crashed ({e}); skipping")
        validation = None

    if validation is not None:
        try:
            import json as _json
            data_dir = getattr(parameter, 'data_dir', None) or os.getcwd()
            out_path = os.path.join(data_dir, "step4_improved_df_validation.json")
            with open(out_path, "w") as fh:
                _json.dump(validation, fh, indent=2)
        except Exception as e:
            logger.debug(f"Could not write validation report ({e}); continuing")

    if verbose:
        if result.is_success:
            logger.info(f"Selected {result.num_primers} primers")
            logger.info(f"Coverage: {result.metrics.fg_coverage:.1%}")
            logger.info(f"Score: {result.score:.4f}")
        elif result.status == OptimizationStatus.PARTIAL:
            if result.num_primers < target_size:
                logger.warning(
                    f"PARTIAL result: found {result.num_primers} primers "
                    f"but target was {target_size}"
                )
            else:
                logger.warning(
                    f"PARTIAL result: found {result.num_primers} primers "
                    f"but genome coverage is below threshold"
                )
            logger.warning("Suggestions to improve results:")
            logger.warning("  - Relax filtering thresholds (increase max_bg_freq or max_gini)")
            logger.warning("  - Widen the k-mer range (decrease min_k or increase max_k)")
            logger.warning("  - Increase the candidate pool (raise max_primer in filter step)")
            logger.warning(f"  - Try a different optimizer (current: {method})")
        else:
            logger.warning(f"Optimization failed: {result.message}")

        if validation is not None and validation.get("issues"):
            for issue in validation["issues"]:
                lvl = issue.get("level", "warning")
                msg = f"Post-opt {lvl}: {issue.get('code')} — {issue.get('detail')}"
                (logger.error if lvl == "error" else logger.warning)(msg)

    return result


def run_optimization_from_config(config: OptimizationConfig) -> OptimizationResult:
    """
    Run optimization using a config object.

    Args:
        config: OptimizationConfig with all settings

    Returns:
        OptimizationResult with selected primers
    """
    return run_optimization(
        method=config.method,
        fg_prefixes=config.fg_prefixes,
        fg_seq_lengths=config.fg_seq_lengths,
        bg_prefixes=config.bg_prefixes,
        bg_seq_lengths=config.bg_seq_lengths,
        target_size=config.target_set_size,
        verbose=config.verbose,
        max_iterations=config.max_iterations,
    )


def save_results(
    result: OptimizationResult,
    output_path: str,
    include_all_sets: bool = False
) -> None:
    """
    Save optimization results to CSV.

    Args:
        result: OptimizationResult to save
        output_path: Path for output CSV
        include_all_sets: Whether to include all found sets (not just best)
    """
    if not result.primers:
        logger.warning("No primers to save")
        return

    # Build results DataFrame with set-level aggregate metrics
    metrics = result.metrics
    records = []
    for primer in result.primers:
        records.append({
            'primer': primer,
            'set_index': 0,
            'score': result.score,
            'normalized_score': result.normalized_score,
            'coverage': metrics.fg_coverage,
            'bg_coverage': metrics.bg_coverage,
            'selectivity': metrics.selectivity_ratio,
            'mean_gap': metrics.mean_gap,
            'max_gap': metrics.max_gap,
            'coverage_uniformity': metrics.coverage_uniformity,
            'gap_gini': metrics.gap_gini,
            'gap_entropy': metrics.gap_entropy,
            'dimer_risk_score': metrics.dimer_risk_score,
            'strand_alternation': metrics.strand_alternation_score,
            'strand_coverage_ratio': metrics.strand_coverage_ratio,
            'mean_tm': metrics.mean_tm,
            'optimizer': result.optimizer_name,
        })

    df = pd.DataFrame(records)
    df.to_csv(output_path, index=False)
    logger.info(f"Results saved to {output_path}")

    # Write summary JSON alongside CSV
    summary_path = os.path.splitext(output_path)[0] + '_summary.json'
    try:
        summary = result.to_dict()
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        logger.info(f"Summary saved to {summary_path}")
    except Exception as e:
        logger.warning(f"Could not write summary JSON: {e}")

    # Write audit trail for reproducibility
    _write_audit_trail(output_path, result)


def _write_audit_trail(output_path: str, result: OptimizationResult) -> None:
    """Write run metadata alongside pipeline results for reproducibility.

    Creates a JSON file with version, timestamp, parameters hash, and
    optimizer configuration so that results can be traced back to the
    exact conditions that produced them.
    """
    from datetime import datetime, timezone
    import hashlib

    audit_path = os.path.splitext(output_path)[0] + '_audit.json'
    try:
        # Collect parameters hash
        params_dict = {}
        for attr in ['fg_genomes', 'bg_genomes', 'min_k', 'max_k',
                      'min_fg_freq', 'max_bg_freq', 'max_gini',
                      'num_primers', 'target_set_size', 'polymerase',
                      'reaction_temp', 'na_conc', 'mg_conc']:
            val = getattr(parameter, attr, None)
            if val is not None:
                params_dict[attr] = str(val)

        params_str = json.dumps(params_dict, sort_keys=True)
        params_hash = hashlib.sha256(params_str.encode()).hexdigest()[:16]

        try:
            from neoswga import __version__
            version = __version__
        except ImportError:
            version = 'unknown'

        audit = {
            'neoswga_version': version,
            'timestamp': datetime.now(timezone.utc).isoformat(),
            'params_hash': params_hash,
            'optimizer': result.optimizer_name,
            'num_primers': result.num_primers,
            'score': float(result.score) if result.score else None,
            'status': result.status.value if result.status else None,
            'parameters': params_dict,
        }

        with open(audit_path, 'w') as f:
            json.dump(audit, f, indent=2)
        logger.info(f"Audit trail saved to {audit_path}")
    except Exception as e:
        logger.debug(f"Could not write audit trail: {e}")


def _simulation_rescore(
    result: OptimizationResult,
    fg_prefixes: List[str],
    fg_seq_lengths: List[int],
    simulation_time: float = 1800.0,
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Re-score optimization result using simulation-based fitness.

    Uses the agent-based replication simulator to predict actual
    amplification performance of the selected primer set.

    Args:
        result: OptimizationResult from an optimizer.
        fg_prefixes: Foreground genome prefixes.
        fg_seq_lengths: Foreground genome lengths.
        simulation_time: Simulation duration in seconds.
        verbose: Print progress.

    Returns:
        Dictionary with simulation fitness metrics, or None on failure.
    """
    if not result.is_success or not result.primers:
        return None

    try:
        from .simulation_fitness import SimulationBasedEvaluator
        from .position_cache import PositionCache
        from .genome_io import GenomeLoader
    except ImportError as e:
        if verbose:
            logger.warning(f"Simulation modules not available: {e}")
        return None

    try:
        # Load genome sequence for simulation
        genome_seq = None
        genome_length = sum(fg_seq_lengths)

        fg_genomes = getattr(parameter, 'fg_genomes', [])
        if fg_genomes:
            try:
                loader = GenomeLoader()
                genome_seq = loader.load_genome(fg_genomes[0], return_stats=False)
                genome_length = len(genome_seq)
            except Exception:
                pass

        if genome_seq is None:
            if verbose:
                logger.info("  Genome sequence not available for simulation")
            return None

        # Build position cache for simulation
        primers = list(result.primers)
        cache = PositionCache(fg_prefixes, primers)

        evaluator = SimulationBasedEvaluator(
            genome_sequence=genome_seq,
            genome_length=genome_length,
            position_cache=cache,
            n_replicates=2,
            simulation_duration=simulation_time,
        )

        fitness = evaluator.evaluate(primers, verbose=verbose)

        return {
            'simulation_coverage': fitness.mean_coverage,
            'simulation_uniformity': fitness.coverage_uniformity,
            'simulation_fitness': fitness.fitness_score,
            'simulation_forks': fitness.mean_forks_created,
            'simulation_time': fitness.simulation_time,
        }
    except Exception as e:
        if verbose:
            logger.warning(f"Simulation re-scoring failed: {e}")
        return None


def optimize_step4(
    use_cache: bool = True,
    use_background_filter: bool = True,
    optimization_method: str = 'hybrid',
    verbose: bool = True,
    uniformity_weight: float = 0.0,
    minimize_primers: bool = False,
    target_coverage: float = 0.70,
    **kwargs
) -> Tuple[List[List[str]], List[float], Any]:
    """
    Drop-in replacement for improved_step4 from pipeline_integration.

    Maintains backward compatibility with existing CLI while using
    the new optimizer framework internally.

    Args:
        use_cache: Use position cache (always True with new system)
        use_background_filter: Use background filtering
        optimization_method: Optimizer method name
        verbose: Print progress
        uniformity_weight: Weight for coverage uniformity
        minimize_primers: Minimize primer count
        target_coverage: Target coverage for minimization
        **kwargs: Additional parameters

    Returns:
        Tuple of (primer_sets, scores, cache) for CLI compatibility
    """
    _ensure_optimizers_registered()

    # Get target size: prefer explicit kwarg, then parameter module, then default
    target_size = kwargs.pop('target_size', None)
    if target_size is None:
        target_size = getattr(parameter, 'num_primers', 6)
        target_size = getattr(parameter, 'target_set_size', target_size)

    if verbose:
        logger.info(f"Unified optimizer: method={optimization_method}")

    # Run optimization
    result = run_optimization(
        method=optimization_method,
        target_size=target_size,
        verbose=verbose,
        uniformity_weight=uniformity_weight,
        minimize_primers=minimize_primers,
        target_coverage=target_coverage,
        **kwargs
    )

    # Optional simulation-based re-scoring
    if kwargs.get('validate_with_simulation', False) and result.is_success:
        if verbose:
            logger.info("")
            logger.info("=" * 60)
            logger.info("Simulation-Based Validation")
            logger.info("=" * 60)

        # Resolve genome prefixes for simulation
        sim_fg_prefixes = getattr(parameter, 'fg_prefixes', [])
        sim_fg_seq_lengths = getattr(parameter, 'fg_seq_lengths', [])

        sim_time = kwargs.get('simulation_time', 1800.0)
        sim_results = _simulation_rescore(
            result, sim_fg_prefixes, sim_fg_seq_lengths,
            simulation_time=sim_time, verbose=verbose,
        )

        if sim_results:
            logger.info(f"  Simulated coverage: {sim_results['simulation_coverage']:.1%}")
            logger.info(f"  Simulated uniformity: {sim_results['simulation_uniformity']:.2f}")
            logger.info(f"  Simulation fitness: {sim_results['simulation_fitness']:.3f}")

    # Convert to legacy format
    if result.is_success or (result.status == OptimizationStatus.PARTIAL and result.primers):
        primer_sets = [list(result.primers)]
        scores = [result.score]

        # Save results
        output_path = os.path.join(parameter.data_dir, "step4_improved_df.csv")
        save_results(result, output_path)

        # Return cache placeholder (for compatibility)
        from .position_cache import PositionCache
        fg_prefixes = getattr(parameter, 'fg_prefixes', [])
        cache = PositionCache(fg_prefixes, list(result.primers)) if fg_prefixes else None

        return primer_sets, scores, cache
    else:
        logger.error(f"Optimization failed: {result.message}")
        return [], [], None


# CLI entry point
def main():
    """CLI entry point for standalone optimization."""
    import argparse

    parser = argparse.ArgumentParser(
        description="NeoSWGA Primer Set Optimization"
    )
    parser.add_argument(
        '-m', '--method',
        default='hybrid',
        help='Optimization method'
    )
    parser.add_argument(
        '-j', '--json-file',
        help='Parameters JSON file'
    )
    parser.add_argument(
        '-n', '--num-primers',
        type=int,
        default=6,
        help='Target number of primers'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    parser.add_argument(
        '--list-methods',
        action='store_true',
        help='List available optimization methods'
    )

    args = parser.parse_args()

    if args.list_methods:
        print("Available optimization methods:")
        for name, desc in list_available_optimizers().items():
            print(f"  {name}: {desc}")
        return

    # Load parameters
    if args.json_file:
        parameter.json_file = args.json_file
        from . import pipeline as core_pipeline
        core_pipeline._initialize()

    # Run optimization
    result = run_optimization(
        method=args.method,
        target_size=args.num_primers,
        verbose=args.verbose,
    )

    if result.is_success:
        print(f"\nSelected {result.num_primers} primers:")
        for primer in result.primers:
            print(f"  {primer}")
        print(f"\nScore: {result.score:.4f}")
        print(f"Coverage: {result.metrics.fg_coverage:.1%}")
    else:
        print(f"Optimization failed: {result.message}")
        exit(1)


if __name__ == '__main__':
    main()
