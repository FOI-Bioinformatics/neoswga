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

import os
import logging
import threading
import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional, Any

from .base_optimizer import OptimizationResult, OptimizerConfig
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

        _optimizers_registered = True
        logger.debug("Optimizer registration complete")


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

    # Load candidates from step3 if not provided
    if candidates is None:
        step3_path = os.path.join(parameter.data_dir, "step3_df.csv")
        if not os.path.exists(step3_path):
            return OptimizationResult.failure(
                method, f"Step 3 output not found: {step3_path}"
            )
        step3_df = pd.read_csv(step3_path)
        candidates = step3_df['primer'].tolist()

    if verbose:
        logger.info(f"Running {method} optimization")
        logger.info(f"  Candidates: {len(candidates)}")
        logger.info(f"  Target size: {target_size}")

    # Create position cache
    with progress_context("Loading position data", disable=not verbose):
        cache = PositionCache(fg_prefixes, candidates)

    # Build optimizer config
    config = OptimizerConfig(
        target_set_size=target_size,
        max_iterations=kwargs.get('max_iterations', 100),
        verbose=verbose,
    )

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
            **kwargs
        )
    except Exception as e:
        logger.error(f"Failed to create optimizer '{method}': {e}")
        return OptimizationResult.failure(method, str(e))

    # Run optimization
    with progress_context(f"Running {optimizer.name} optimizer", disable=not verbose):
        result = optimizer.optimize(candidates, target_size)

    if verbose:
        if result.is_success:
            logger.info(f"Selected {result.num_primers} primers")
            logger.info(f"Coverage: {result.metrics.fg_coverage:.1%}")
            logger.info(f"Score: {result.score:.4f}")
        else:
            logger.warning(f"Optimization failed: {result.message}")

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

    # Build results DataFrame
    records = []
    for primer in result.primers:
        records.append({
            'primer': primer,
            'set_index': 0,
            'score': result.score,
            'coverage': result.metrics.fg_coverage,
            'bg_coverage': result.metrics.bg_coverage,
            'selectivity': result.metrics.selectivity_ratio,
            'mean_gap': result.metrics.mean_gap,
            'optimizer': result.optimizer_name,
        })

    df = pd.DataFrame(records)
    df.to_csv(output_path, index=False)

    logger.info(f"Results saved to {output_path}")


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

    # Get target size from parameters
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

    # Convert to legacy format
    if result.is_success:
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
