#!/usr/bin/env python3
"""
Two-stage hybrid primer set optimization for SWGA.

Combines the strengths of two complementary approaches:
1. Dominating-set (Stage 1): Maximize genome coverage
2. Network-based (Stage 2): Maximize amplification connectivity

This hybrid approach addresses the key limitation of each method:
- Dominating-set ignores amplification network structure
- Network-based can have poor coverage in sparse regions

By combining them, we get both broad coverage AND efficient amplification.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.0 - Phase 2.1
"""

import logging
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from dataclasses import replace as _dc_replace
from typing import Dict, List, Optional, Tuple

import networkx as nx
import numpy as np

from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.network_optimizer import AmplificationNetwork, NetworkOptimizer
from neoswga.core.registry import POLYMERASES as _POLYMERASES

logger = logging.getLogger(__name__)

# How much Stage-2 removal weighs the coverage a removal costs against the
# amplification-network score it leaves behind. Both terms are normalised
# within each removal step, so this is a genuine 0-1 trade-off rather than a
# coefficient between incommensurable units. Half and half: Stage 2 exists to
# improve connectivity, so coverage must not simply dominate it and turn the
# stage into a second set-cover pass.
_STAGE2_COVERAGE_WEIGHT = 0.5


def _removal_network_score(connectivity: float, largest_component: int) -> float:
    """How good the amplification network left behind by a removal is.

    This was `connectivity + predicted_fold / 100`, and the predicted fold is
    capped at 2**20. Any component past roughly 200 sites reaches that cap --
    which is every real bacterial genome -- so every candidate removal scored
    the identical constant and the term ranked nothing. On a M. tuberculosis
    run the largest component was 410 sites and algebraic connectivity read
    0.00 throughout, because the site graph is disconnected; between them
    Stage 2 had no signal at all, spent 150 s, and produced a set covering
    53.2% against Stage 1's 54.8%, which the guard then threw away.

    Ranking on the component size instead costs nothing: fold is a monotone
    function of it, so the ordering is identical everywhere the fold is
    informative, and it stays an ordering where the fold has saturated. The
    fold itself is unchanged and still what gets reported -- capping a
    predicted yield is reasonable, ranking on a capped value is not.
    """
    return connectivity + largest_component


# =========================================================================
# Polymerase presets for polymerase-aware optimization
# =========================================================================


@dataclass
class PolymeraseConfig:
    """Per-polymerase configuration defaults for hybrid optimization."""

    max_extension: int = 70000
    thermo_filter: bool = False
    primer_multiplier: float = 1.0
    reaction_temp: float = 30.0
    min_primer_tm: float = 20.0
    max_primer_tm: float = 50.0
    min_gc: float = 0.25
    max_gc: float = 0.75


# Standard polymerase presets, derived from the registry.
#
# NOTE `max_extension` is the Stage-2 amplification-NETWORK reach -- "could two
# primers connect via one uninterrupted extension?" -- so single-molecule
# processivity is the right quantity for it. The Stage-1 set-cover COVERAGE
# objective uses `coverage_reach` (the realistic per-primer reach), threaded
# separately by unified_optimizer.
#
# These were previously hand-maintained numbers that disagreed with processivity
# for bst (10000 vs 2000, physically impossible) and klenow (5000 vs 40). They
# now derive from processivity_bp, so the two cannot diverge again.
POLYMERASE_PRESETS: Dict[str, PolymeraseConfig] = {
    key: PolymeraseConfig(
        max_extension=spec.processivity_bp,
        thermo_filter=spec.thermo_filter,
        primer_multiplier=spec.primer_multiplier,
        reaction_temp=spec.preset_reaction_temp,
        min_primer_tm=spec.primer_tm_range[0],
        max_primer_tm=spec.primer_tm_range[1],
        min_gc=spec.gc_range[0],
        max_gc=spec.gc_range[1],
    )
    for key, spec in _POLYMERASES.items()
}


def _get_polymerase_config(polymerase: str) -> PolymeraseConfig:
    """Get a polymerase config, falling back to phi29 defaults.

    Returns a COPY. `PolymeraseConfig` is a mutable dataclass and
    `HybridOptimizer._adjust_for_gc` writes to `min_gc`/`max_gc` on the instance it
    holds. Handing out the shared preset let one genome's GC adaptation leak into
    every optimizer constructed later in the same process -- so an AT-rich target
    would silently widen the GC window for an unrelated GC-rich target in an
    ensemble, multi-genome, or long-running session.
    """
    preset = POLYMERASE_PRESETS.get(polymerase, POLYMERASE_PRESETS["phi29"])
    return _dc_replace(preset)


@dataclass
class HybridResult:
    """Result from hybrid optimization"""

    # Final primer set
    primers: List[str]

    # Stage 1 (Coverage) results
    stage1_primers: List[str]
    stage1_coverage: float
    stage1_regions_covered: int

    # Stage 2 (Network) results
    stage2_primers: List[str]
    stage2_connectivity: float
    stage2_predicted_amplification: float
    stage2_largest_component: int

    # Overall metrics
    final_coverage: float
    final_connectivity: float
    final_predicted_amplification: float

    # Stage-1 primers in greedy selection order. `stage1_primers` comes from a
    # set and has no meaningful order; this does. Set-cover picks the largest
    # marginal gain at each step, so its first N picks are the same N whatever
    # `stage1_count` was, and coverage over that prefix is non-decreasing in N.
    # Stage 2 is floored against it, which is what makes final coverage
    # monotone in the requested count.
    stage1_ordered_primers: List[str] = field(default_factory=list)

    # Simulation validation (optional)
    simulation_fitness: Optional[object] = None  # SimulationFitness if validated

    # Metadata
    runtime_stage1: float = 0.0
    runtime_stage2: float = 0.0
    runtime_simulation: float = 0.0
    total_runtime: float = 0.0

    def __str__(self):
        base = f"""Hybrid Optimization Result:
  Final primers: {len(self.primers)}
  Coverage: {self.final_coverage:.1%}
  Connectivity: {self.final_connectivity:.2f}
  Amplification: {self.final_predicted_amplification:.1f}×"""

        if self.simulation_fitness:
            base += f"""
  Simulation validation:
    Coverage: {self.simulation_fitness.mean_coverage:.1%} ± {self.simulation_fitness.std_coverage:.1%}
    Uniformity: {self.simulation_fitness.coverage_uniformity:.2f}
    Fitness: {self.simulation_fitness.fitness_score:.3f}"""

        base += f"\n  Runtime: {self.total_runtime:.2f}s"
        return base


class HybridOptimizer:
    """
    Two-stage hybrid optimizer combining coverage and network approaches.

    Stage 1 (Dominating Set):
    - Select N primers (typically 20-25)
    - Goal: Maximize genome coverage
    - Fast (0.1s for typical genomes)

    Stage 2 (Network Refinement):
    - From Stage 1 primers, select M primers (typically 10-15)
    - Goal: Maximize amplification network connectivity
    - Considers strand pairing and extension distances
    - Moderate runtime (~1s for typical sets)

    Expected improvement: 20-30% better amplification vs dominating-set alone
    while maintaining good coverage.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        bin_size: int = 10000,
        max_extension: Optional[int] = None,
        coverage_reach: Optional[int] = None,
        uniformity_weight: float = 0.0,
        min_tm: Optional[float] = None,
        max_tm: Optional[float] = None,
        polymerase: str = "phi29",
        genome_gc_content: Optional[float] = None,
        background_pruning: bool = False,
        background_weight: float = 2.0,
        min_coverage_threshold: float = 0.95,
        conditions=None,
        mechanistic_weight: float = 0.0,
    ):
        """
        Initialize hybrid optimizer.

        Args:
            position_cache: Cache with primer binding positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers (optional)
            bg_seq_lengths: Background genome lengths (optional)
            bin_size: Bin size for coverage analysis (bp)
            max_extension: Maximum extension distance (bp), overridden by
                polymerase preset if not explicitly provided
            uniformity_weight: Weight for coverage uniformity (0.0-1.0)
            polymerase: Polymerase type ('phi29', 'equiphi29', 'bst', 'klenow').
                Applies preset config for thermo-filtering and extension distance.
            genome_gc_content: Target genome GC content (0-1). Used for
                GC-adaptive adjustments when polymerase requires it.
            background_pruning: Enable background-pruning stage between coverage
                and network refinement. Removes primers with high background
                binding while maintaining coverage above min_coverage_threshold.
            background_weight: Weight for background sites in pruning score.
                Higher values favor more aggressive background removal.
            min_coverage_threshold: Minimum coverage to maintain during
                background pruning (0.0-1.0).
        """
        self.position_cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.bin_size = bin_size
        self.uniformity_weight = uniformity_weight

        # Polymerase-aware configuration
        self.polymerase = polymerase
        self.poly_config = _get_polymerase_config(polymerase)
        self.genome_gc_content = genome_gc_content

        # Apply polymerase preset for max_extension (caller can override).
        # max_extension is the single-molecule processivity used for the Stage-2
        # amplification-NETWORK connectivity question ("can two primers connect
        # via one extension?").
        #
        # `None` means "not chosen" rather than the literal 70000 this used to
        # test for. That was a sentinel-by-value: 70000 is phi29's processivity
        # and a perfectly reasonable explicit argument, and a caller who passed
        # it for bst had it silently replaced by 2000 -- a 35-fold change to the
        # network reach, on the strength of a collision between a default and a
        # real value.
        self.max_extension = (
            self.poly_config.max_extension if max_extension is None else max_extension
        )

        # Explicit Tm window, when the caller configured one. Stage-0 screening
        # built its criteria from the polymerase preset alone, so a user who
        # widened max_tm in params.json had the filter step honour it and the
        # optimizer silently narrow it back.
        self.min_tm = min_tm
        self.max_tm = max_tm

        # coverage_reach is the realistic per-primer reach used for the Stage-1
        # set-cover COVERAGE objective, so selection optimizes the same coverage
        # definition that base_optimizer.compute_metrics scores the result on
        # (~3 kb for phi29). Defaults to max_extension only for backward
        # compatibility; unified_optimizer passes the resolved realistic reach.
        self.coverage_reach = coverage_reach if coverage_reach is not None else self.max_extension

        # GC-adaptive adjustments for polymerases that benefit from it
        if genome_gc_content is not None and self.poly_config.thermo_filter:
            self._adjust_for_gc(genome_gc_content)

        # Background pruning configuration
        self.background_pruning = background_pruning
        self.background_weight = background_weight
        self.min_coverage_threshold = min_coverage_threshold

        # Initialize both optimizers
        self.dominating_optimizer = DominatingSetOptimizer(
            cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bin_size=bin_size,
            # Stage-1 coverage uses the realistic per-primer reach so the
            # selection objective matches how the result is scored.
            extension_reach=self.coverage_reach,
        )

        self.network_optimizer = NetworkOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            bg_prefixes=self.bg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_seq_lengths=self.bg_seq_lengths,
            max_extension=self.max_extension,
            uniformity_weight=uniformity_weight,
            # Propagate ReactionConditions so the inner NetworkOptimizer's
            # _get_primer_tm applies additive corrections (DMSO / betaine / etc.)
            # when computing Tm-weighted edges and Tm scores.
            conditions=conditions,
            # Phase 13B: forward --use-mechanistic-model weight so the
            # NetworkOptimizer's scoring includes a mechanistic term.
            mechanistic_weight=mechanistic_weight,
        )
        # Retain for introspection / rescoring hooks.
        self.conditions = conditions

        logger.info("Hybrid optimizer initialized")
        logger.info(f"  Polymerase: {polymerase}")
        logger.info(f"  Bin size: {bin_size:,} bp")
        logger.info(f"  Max extension: {self.max_extension:,} bp")
        if self.poly_config.thermo_filter:
            logger.info(f"  Thermo-filtering: enabled ({self.poly_config.reaction_temp}C)")
        if background_pruning:
            logger.info(f"  Background pruning: enabled (weight={background_weight})")

    def _adjust_for_gc(self, gc: float):
        """
        Adjust polymerase config based on genome GC content.

        GC-rich genomes benefit from higher betaine, longer primers, and
        wider GC acceptance ranges. AT-rich genomes use shorter primers.
        """
        if gc > 0.65:
            self.poly_config.max_gc = 0.80
            logger.info("  GC-rich genome detected - widening GC acceptance")
        elif gc < 0.35:
            self.poly_config.min_gc = 0.20
            logger.info("  AT-rich genome detected - widening GC acceptance")

    def _rescale_for_polymerase(self, final_count: int, verbose: bool) -> int:
        """Apply the polymerase preset's primer-count multiplier and floor.

        Two faults lived at the call site. The assignment sat inside an
        `and verbose` branch, so the rescale took effect only when logging was
        on: equiphi29 asked for 10 returned 10 quiet and 8 loud, while phi29
        (multiplier 1.0) returned 10 either way. A design must not depend on
        how chatty the run is.

        And it overrode counts a caller had named. At the default target of 6
        the multiplier cannot bite -- `max(6, int(6 * 0.85))` is 6 -- so the
        only counts it ever changed were deliberate ones. A params.json asking
        for 32 primers got 27, and the completion check then blamed the
        candidate pool, advice that would degrade a design that was not
        deficient. Hence the call is now opt-in; see `optimize`.
        """
        adjusted = max(6, int(final_count * self.poly_config.primer_multiplier))
        if adjusted != final_count and verbose:
            logger.info(
                f"Adjusted target from {final_count} to {adjusted} "
                f"({self.polymerase} multiplier: {self.poly_config.primer_multiplier})"
            )
        return adjusted

    def optimize(
        self,
        candidates: List[str],
        final_count: int = 12,
        apply_polymerase_multiplier: bool = False,
        stage1_count: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        verbose: bool = True,
        validate_with_simulation: bool = False,
        genome_sequence: Optional[str] = None,
        simulation_replicates: int = 3,
    ) -> HybridResult:
        """
        Two-stage hybrid optimization.

        Args:
            candidates: List of candidate primers (already filtered)
            final_count: Final number of primers to select (Stage 2 output)
            apply_polymerase_multiplier: Let the preset rescale `final_count`
                (see `_rescale_for_polymerase`). Defaults False: a caller who
                names a count means it.
            stage1_count: Number of primers for Stage 1 (None = auto)
            fixed_primers: Optional list of primers that must be included in
                the final set. The optimizer will select additional primers
                to complement these fixed primers. Useful for iterative design
                where some primers have been experimentally validated.
            verbose: Print progress
            validate_with_simulation: Run simulation to validate final set
            genome_sequence: Genome sequence (required if validate_with_simulation=True)
            simulation_replicates: Number of simulation replicates

        Returns:
            HybridResult with complete optimization details
        """
        total_start = time.time()

        # Process fixed primers
        fixed_primers = fixed_primers or []
        fixed_primers = [p.upper() for p in fixed_primers]
        n_fixed = len(fixed_primers)

        # Remove fixed primers from candidates (they're already selected)
        candidates_filtered = [c for c in candidates if c.upper() not in fixed_primers]

        if verbose:
            stages = "Two-Stage"
            if self.poly_config.thermo_filter:
                stages = "Multi-Stage (thermo-filter + coverage + network)"
            if self.background_pruning:
                stages = "Multi-Stage (coverage + bg-pruning + network)"
            if self.poly_config.thermo_filter and self.background_pruning:
                stages = "Multi-Stage (thermo-filter + coverage + bg-pruning + network)"
            logger.info("=" * 80)
            logger.info(f"HYBRID OPTIMIZATION ({stages})")
            logger.info("=" * 80)
            logger.info(f"Input: {len(candidates)} candidates")
            logger.info(f"Polymerase: {self.polymerase}")
            if n_fixed > 0:
                logger.info(f"Fixed primers: {n_fixed} (pre-selected)")
                logger.info(f"Candidates after exclusion: {len(candidates_filtered)}")
            logger.info(f"Target: {final_count} final primers")

        # =================================================================
        # PRE-STAGE: Thermodynamic Filtering (polymerase-dependent)
        # =================================================================
        if self.poly_config.thermo_filter and candidates_filtered:
            candidates_filtered = self._thermo_filter_candidates(
                candidates_filtered, verbose=verbose
            )

        # Adjust target count for fixed primers
        target_new = final_count - n_fixed
        if target_new <= 0:
            if verbose:
                logger.info(f"Fixed primers ({n_fixed}) meet or exceed target ({final_count})")
                logger.info("Returning fixed primers only")

            # Calculate metrics for fixed primers
            network = self._build_network(fixed_primers)
            stats = network.get_statistics()
            coverage = self._calculate_coverage(fixed_primers)

            return HybridResult(
                primers=fixed_primers,
                stage1_primers=fixed_primers,
                stage1_coverage=coverage,
                stage1_regions_covered=0,
                stage2_primers=fixed_primers,
                stage2_connectivity=stats["connectivity"],
                stage2_predicted_amplification=stats["predicted_amplification"],
                stage2_largest_component=stats["largest_component"],
                final_coverage=coverage,
                final_connectivity=stats["connectivity"],
                final_predicted_amplification=stats["predicted_amplification"],
                runtime_stage1=0.0,
                runtime_stage2=0.0,
                runtime_simulation=0.0,
                total_runtime=time.time() - total_start,
            )

        if apply_polymerase_multiplier:
            final_count = self._rescale_for_polymerase(final_count, verbose)

        # Auto-determine Stage 1 count if not specified
        if stage1_count is None:
            # Stage 1 should select more primers than final
            # Rule of thumb: 1.5-2x the final count
            stage1_count = max(final_count + 8, int(final_count * 1.67))
            # When background pruning is enabled, select more for pruning headroom
            if self.background_pruning:
                stage1_count = max(stage1_count, final_count * 2)
            stage1_count = min(stage1_count, len(candidates))

        if verbose:
            if n_fixed > 0:
                logger.info(
                    f"Stage 1 target: {stage1_count - n_fixed} new primers + {n_fixed} fixed"
                )
                logger.info(
                    f"Stage 2 target: {target_new} new primers + {n_fixed} fixed = {final_count} total"
                )
            else:
                logger.info(f"Stage 1 target: {stage1_count} primers (coverage)")
                logger.info(f"Stage 2 target: {final_count} primers (network)")

        # ===================================================================
        # STAGE 1: Dominating Set (Coverage Optimization)
        # ===================================================================

        if verbose:
            logger.info("\n" + "-" * 80)
            logger.info("STAGE 1: Coverage Optimization (Dominating Set)")
            logger.info("-" * 80)
            if n_fixed > 0:
                logger.info(f"Pre-selecting {n_fixed} fixed primers")

        stage1_start = time.time()

        # Adjust stage1 count to account for fixed primers
        stage1_new_count = max(1, stage1_count - n_fixed)

        stage1_result = self.dominating_optimizer.optimize_greedy(
            candidates=candidates_filtered,
            max_primers=stage1_new_count,
            fixed_primers=fixed_primers,
            verbose=verbose,
        )

        stage1_runtime = time.time() - stage1_start

        # Combine fixed primers with newly selected primers
        stage1_new_primers = stage1_result["primers"]
        stage1_primers = fixed_primers + [p for p in stage1_new_primers if p not in fixed_primers]
        # Greedy selection order, which `stage1_primers` (built from a set) loses.
        stage1_ordered = list(stage1_result.get("ordered_primers") or stage1_primers)
        stage1_coverage = stage1_result["coverage"]
        stage1_regions = stage1_result["covered_regions"]

        if verbose:
            logger.info(f"\nStage 1 complete:")
            if n_fixed > 0:
                logger.info(f"  Fixed primers: {n_fixed}")
                logger.info(f"  New primers: {len(stage1_new_primers)}")
                logger.info(f"  Total: {len(stage1_primers)} primers")
            else:
                logger.info(f"  Selected: {len(stage1_primers)} primers")
            logger.info(f"  Coverage: {stage1_coverage:.1%}")
            logger.info(f"  Regions: {stage1_regions}")
            logger.info(f"  Runtime: {stage1_runtime:.2f}s")

        # If Stage 1 gave us fewer primers than target, use them all
        if len(stage1_primers) <= final_count:
            if verbose:
                logger.info(f"\nStage 1 selected <= {final_count} primers, skipping Stage 2")

            # Calculate network metrics for these primers
            network = self._build_network(stage1_primers)
            stats = network.get_statistics()

            total_runtime = time.time() - total_start

            return HybridResult(
                primers=stage1_primers,
                stage1_primers=stage1_primers,
                stage1_ordered_primers=stage1_ordered,
                stage1_coverage=stage1_coverage,
                stage1_regions_covered=stage1_regions,
                stage2_primers=stage1_primers,
                stage2_connectivity=stats["connectivity"],
                stage2_predicted_amplification=stats["predicted_amplification"],
                stage2_largest_component=stats["largest_component"],
                final_coverage=stage1_coverage,
                final_connectivity=stats["connectivity"],
                final_predicted_amplification=stats["predicted_amplification"],
                simulation_fitness=None,  # Skip simulation for early return
                runtime_stage1=stage1_runtime,
                runtime_stage2=0.0,
                runtime_simulation=0.0,
                total_runtime=total_runtime,
            )

        # ===================================================================
        # STAGE 1.5 (OPTIONAL): Background Pruning
        # ===================================================================

        stage1_5_runtime = 0.0
        if self.background_pruning and self.bg_prefixes:
            if verbose:
                logger.info("\n" + "-" * 80)
                logger.info("STAGE 1.5: Background Pruning")
                logger.info("-" * 80)

            stage1_5_start = time.time()

            # Target: keep enough primers for network refinement, prune the rest
            bg_prune_target = int(final_count * 1.5)

            stage1_primers, prune_coverage, prune_bg_sites = self._prune_background(
                stage1_primers,
                target_size=bg_prune_target,
                verbose=verbose,
            )
            stage1_coverage = prune_coverage
            stage1_5_runtime = time.time() - stage1_5_start

            if verbose:
                logger.info(f"\nBackground pruning complete:")
                logger.info(f"  Primers: {len(stage1_primers)}")
                logger.info(f"  Coverage: {prune_coverage:.1%}")
                logger.info(f"  Background sites: {prune_bg_sites}")
                logger.info(f"  Runtime: {stage1_5_runtime:.2f}s")

        # ===================================================================
        # STAGE 2: Network Refinement (Amplification Optimization)
        # ===================================================================

        if verbose:
            logger.info("\n" + "-" * 80)
            logger.info("STAGE 2: Network Refinement (Amplification)")
            logger.info("-" * 80)

        stage2_start = time.time()

        # Use network-based selection from Stage 1 primers
        # Fixed primers will never be removed during refinement
        stage2_primers = self._network_refine(
            stage1_primers,
            target_count=final_count,
            fixed_primers=fixed_primers,
            verbose=verbose,
            coverage_floor_set=stage1_ordered,
        )

        stage2_runtime = time.time() - stage2_start

        # Calculate final metrics
        final_network = self._build_network(stage2_primers)
        final_stats = final_network.get_statistics()

        # Calculate final coverage
        final_coverage_result = self._calculate_coverage(stage2_primers)

        if verbose:
            logger.info(f"\nStage 2 complete:")
            logger.info(f"  Selected: {len(stage2_primers)} primers")
            logger.info(f"  Connectivity: {final_stats['connectivity']:.2f}")
            logger.info(f"  Amplification: {final_stats['predicted_amplification']:.1f}×")
            logger.info(f"  Runtime: {stage2_runtime:.2f}s")

        # ===================================================================
        # STAGE 3 (OPTIONAL): Simulation Validation
        # ===================================================================

        simulation_fitness = None
        simulation_runtime = 0.0

        if validate_with_simulation:
            if genome_sequence is None:
                logger.warning(
                    "Simulation validation requested but no genome sequence provided - skipping"
                )
            elif len(stage2_primers) == 0:
                logger.warning("No primers selected - skipping simulation validation")
            else:
                if verbose:
                    logger.info("\n" + "-" * 80)
                    logger.info("STAGE 3: Simulation Validation")
                    logger.info("-" * 80)

                simulation_start = time.time()

                try:
                    from neoswga.core.simulation_fitness import SimulationBasedEvaluator

                    evaluator = SimulationBasedEvaluator(
                        genome_sequence=genome_sequence,
                        genome_length=self.fg_seq_lengths[0],
                        position_cache=self.position_cache,
                        n_replicates=simulation_replicates,
                    )

                    simulation_fitness = evaluator.evaluate(stage2_primers, verbose=verbose)

                except Exception as e:
                    logger.warning(f"Simulation validation failed: {e}")

                simulation_runtime = time.time() - simulation_start

                if verbose and simulation_fitness:
                    logger.info(f"\nStage 3 complete:")
                    logger.info(f"  Simulated coverage: {simulation_fitness.mean_coverage:.1%}")
                    logger.info(f"  Fitness score: {simulation_fitness.fitness_score:.3f}")
                    logger.info(f"  Runtime: {simulation_runtime:.2f}s")

        total_runtime = time.time() - total_start

        if verbose:
            logger.info("\n" + "=" * 80)
            logger.info("HYBRID OPTIMIZATION COMPLETE")
            logger.info("=" * 80)
            logger.info(f"Final primers: {len(stage2_primers)}")
            logger.info(f"Final coverage: {final_coverage_result:.1%}")
            logger.info(f"Final connectivity: {final_stats['connectivity']:.2f}")
            logger.info(f"Final amplification: {final_stats['predicted_amplification']:.1f}×")
            if simulation_fitness:
                logger.info(f"Simulation fitness: {simulation_fitness.fitness_score:.3f}")
            logger.info(f"Total runtime: {total_runtime:.2f}s")
            logger.info("=" * 80)

        return HybridResult(
            primers=stage2_primers,
            stage1_primers=stage1_primers,
            stage1_ordered_primers=stage1_ordered,
            stage1_coverage=stage1_coverage,
            stage1_regions_covered=stage1_regions,
            stage2_primers=stage2_primers,
            stage2_connectivity=final_stats["connectivity"],
            stage2_predicted_amplification=final_stats["predicted_amplification"],
            stage2_largest_component=final_stats["largest_component"],
            final_coverage=final_coverage_result,
            final_connectivity=final_stats["connectivity"],
            final_predicted_amplification=final_stats["predicted_amplification"],
            simulation_fitness=simulation_fitness,
            runtime_stage1=stage1_runtime,
            runtime_stage2=stage2_runtime,
            runtime_simulation=simulation_runtime,
            total_runtime=total_runtime,
        )

    def _network_refine(
        self,
        primers: List[str],
        target_count: int,
        fixed_primers: Optional[List[str]] = None,
        verbose: bool = True,
        coverage_floor_set: Optional[List[str]] = None,
    ) -> List[str]:
        """
        Refine primer set using network analysis.

        From a set of primers with good coverage, select subset that
        maximizes amplification network connectivity.

        Uses O(1) subgraph views instead of rebuilding the full network
        for each removal candidate, providing 10-50x speedup.

        Args:
            primers: Input primer set (from Stage 1)
            target_count: Number of primers to select
            fixed_primers: Primers that must not be removed
            verbose: Print progress

        Returns:
            Refined primer set
        """
        if len(primers) <= target_count:
            return primers

        fixed_set = set(fixed_primers) if fixed_primers else set()

        if verbose:
            logger.info(f"Refining {len(primers)} primers → {target_count}")
            if fixed_set:
                logger.info(f"  {len(fixed_set)} primers are fixed and will not be removed")

        # Build network once (instead of per-candidate rebuilds)
        full_network = self._build_network(primers)
        initial_stats = full_network.get_statistics()

        if verbose:
            logger.info(
                f"Initial network: {initial_stats['largest_component']} sites in largest component"
            )

        # Pre-index nodes by primer for O(1) lookup
        nodes_by_primer = defaultdict(set)
        for node in full_network.graph.nodes():
            nodes_by_primer[node.primer].add(node)

        # Track current node set for efficient subgraph views
        current_primers = primers.copy()
        current_node_set = set(full_network.graph.nodes())

        # Coverage bookkeeping for the removal criterion.
        #
        # The criterion used to be `connectivity + pred_amp / 100`, with no
        # coverage term at all. That is not a neutral omission: the primer
        # whose removal leaves the best network score is the most peripheral
        # one in the amplification graph, and a primer is peripheral precisely
        # when its sites sit far from everything else -- which is to say when
        # it is the only thing covering its region. Stage 2 was therefore
        # biased towards dropping exactly the primers carrying unique
        # coverage, and took a 94.0%-covering Stage-1 set down to 64.9%.
        #
        # `_calculate_coverage` per candidate per step would undo the O(1)
        # subgraph optimisation above, so instead count how many primers cover
        # each bin once, and get the cost of a removal from the bins where that
        # count is 1.
        bins_by_primer = self._coverage_bins_by_primer(current_primers)
        bin_counts = Counter()
        for owned in bins_by_primer.values():
            bin_counts.update(owned)
        total_bins = self._total_coverage_bins()

        while len(current_primers) > target_count:
            best_to_remove = None
            best_score_after_removal = -float("inf")
            candidates = []

            # Try removing each primer using subgraph views (O(1) each)
            for primer in current_primers:
                if primer in fixed_set:
                    continue

                # Create subgraph view without this primer's nodes
                remaining_nodes = current_node_set - nodes_by_primer.get(primer, set())
                if not remaining_nodes:
                    continue

                subgraph = full_network.graph.subgraph(remaining_nodes)

                # Compute connectivity on subgraph view
                if len(subgraph) < 2:
                    connectivity = 0.0
                else:
                    try:
                        connectivity = nx.algebraic_connectivity(subgraph, method="tracemin_lu")
                    except (nx.NetworkXError, ValueError, np.linalg.LinAlgError):
                        connectivity = 0.0

                components = list(nx.connected_components(subgraph))
                largest = max(len(c) for c in components) if components else 0
                network_score = _removal_network_score(connectivity, largest)

                # Bins only this primer covers: what removing it costs.
                unique_bins = sum(1 for b in bins_by_primer.get(primer, ()) if bin_counts[b] == 1)
                candidates.append((primer, network_score, unique_bins))

            if not candidates:
                best_to_remove = None
            else:
                # Rank on both axes within this step rather than adding raw
                # values: algebraic connectivity and a bin count share no
                # scale, and a fixed coefficient between them would be
                # arbitrary and input-dependent. Normalising per step keeps the
                # trade-off meaningful whatever the magnitudes happen to be.
                net_values = [c[1] for c in candidates]
                cost_values = [c[2] for c in candidates]
                net_lo, net_hi = min(net_values), max(net_values)
                cost_lo, cost_hi = min(cost_values), max(cost_values)
                net_span = (net_hi - net_lo) or 1.0
                cost_span = (cost_hi - cost_lo) or 1.0

                best_combined = -float("inf")
                for primer, network_score, unique_bins in candidates:
                    norm_net = (network_score - net_lo) / net_span
                    # Cheap to remove == loses few unique bins == score 1.
                    norm_keep = 1.0 - (unique_bins - cost_lo) / cost_span
                    combined = (
                        1.0 - _STAGE2_COVERAGE_WEIGHT
                    ) * norm_net + _STAGE2_COVERAGE_WEIGHT * norm_keep
                    if combined > best_combined:
                        best_combined = combined
                        best_to_remove = primer
                best_score_after_removal = best_combined

            if best_to_remove:
                current_primers.remove(best_to_remove)
                current_node_set -= nodes_by_primer.get(best_to_remove, set())
                for b in bins_by_primer.get(best_to_remove, ()):
                    bin_counts[b] -= 1
                    if bin_counts[b] <= 0:
                        del bin_counts[b]
                if verbose and len(current_primers) % 5 == 0:
                    logger.info(f"  Reduced to {len(current_primers)} primers...")
            else:
                if verbose:
                    logger.info(f"  Cannot reduce further - remaining primers are fixed")
                break

        # Compute final stats using subgraph view
        final_subgraph = full_network.graph.subgraph(current_node_set)
        final_components = list(nx.connected_components(final_subgraph))
        final_largest = max(len(c) for c in final_components) if final_components else 0

        if verbose:
            try:
                final_connectivity = (
                    nx.algebraic_connectivity(final_subgraph, method="tracemin_lu")
                    if len(final_subgraph) >= 2
                    else 0.0
                )
            except (nx.NetworkXError, ValueError, np.linalg.LinAlgError):
                final_connectivity = 0.0
            logger.info(f"Final network: {final_largest} sites in largest component")
            logger.info(
                f"Connectivity improved: {initial_stats['connectivity']:.2f} → {final_connectivity:.2f}"
            )

        # Floor the result at the Stage-1 greedy prefix of the same size.
        #
        # Coverage-aware removal makes Stage 2 much less destructive, but "much
        # less destructive" is a tendency, not a guarantee -- it holds on the
        # inputs someone happened to try. Greedy set cover selects in a
        # deterministic order, so its first N picks are the same N whatever
        # stage1_count was, and coverage over that prefix is non-decreasing in
        # N. Refusing to return anything worse than it turns the tendency into
        # a floor that is itself monotone in the requested count.
        if coverage_floor_set:
            baseline = list(coverage_floor_set)[:target_count]
            if len(baseline) == target_count:
                refined_coverage = self._calculate_coverage(current_primers)
                baseline_coverage = self._calculate_coverage(baseline)
                if baseline_coverage > refined_coverage + 1e-12:
                    if verbose:
                        logger.info(
                            f"  Refinement covered {refined_coverage:.1%} against "
                            f"{baseline_coverage:.1%} for the stage-1 prefix; "
                            f"keeping the prefix."
                        )
                    return baseline

        return current_primers

    def _coverage_bins_by_primer(self, primers: List[str]) -> Dict[str, set]:
        """Bins each primer covers, at the granularity `_calculate_coverage` uses.

        Computed once per refinement so the removal loop can price a candidate
        from a bin-occupancy counter instead of recomputing coverage, which
        would undo the O(1) subgraph views the loop is built around.
        """
        from neoswga.core.dominating_set_optimizer import BipartiteGraph, coverage_bin_size

        bin_size = coverage_bin_size(self.bin_size, self.coverage_reach)
        owned: Dict[str, set] = {}

        for primer in primers:
            graph = BipartiteGraph(bin_size=bin_size)
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                fw = self.position_cache.get_positions(prefix, primer, "forward")
                rv = self.position_cache.get_positions(prefix, primer, "reverse")
                positions = np.concatenate([fw, rv])
                if len(positions) > 0:
                    graph.add_primer_coverage(
                        primer,
                        positions,
                        prefix,
                        length,
                        extension_reach=self.coverage_reach,
                    )
            owned[primer] = set(graph.regions)
        return owned

    def _total_coverage_bins(self) -> int:
        from neoswga.core.dominating_set_optimizer import coverage_bin_size

        bin_size = coverage_bin_size(self.bin_size, self.coverage_reach)
        return sum((length + bin_size - 1) // bin_size for length in self.fg_seq_lengths)

    def _build_network(self, primers: List[str]) -> AmplificationNetwork:
        """Build amplification network for primer set"""
        network = AmplificationNetwork(max_extension=self.max_extension)

        for primer in primers:
            for prefix in self.fg_prefixes:
                # Use 'forward'/'reverse' (not '+'/'-') for PositionCache API
                positions_fwd = self.position_cache.get_positions(prefix, primer, "forward")
                positions_rev = self.position_cache.get_positions(prefix, primer, "reverse")

                if len(positions_fwd) > 0:
                    network.add_primer_sites(primer, positions_fwd, "+")
                if len(positions_rev) > 0:
                    network.add_primer_sites(primer, positions_rev, "-")

        network.build_edges()
        return network

    def _calculate_coverage(self, primers: List[str]) -> float:
        """Binned genome coverage for a primer set, at the realistic reach.

        This is an APPROXIMATION used for progress reporting and for the
        coverage floor in background pruning. It works in bins, so a bin counts
        as covered when any part of it falls within reach of a site -- which is
        only sound while a bin is no larger than the reach. `coverage_bin_size`
        enforces that: the configured 10 kb default is cut to the ~3 kb reach,
        because a 10 kb bin let a primer covering 30% of it claim the whole bin
        and report 1.000 for a set that covers 0.433.

        The authoritative number is `PrimerSetMetrics.fg_coverage`, computed
        base-by-base by `BaseOptimizer._compute_coverage` and written to
        `step4_improved_df_summary.json`. Compare like with like: this method
        and that metric will not agree exactly, by construction.
        """
        from neoswga.core.dominating_set_optimizer import BipartiteGraph, coverage_bin_size

        bin_size = coverage_bin_size(self.bin_size, self.coverage_reach)
        graph = BipartiteGraph(bin_size=bin_size)

        for primer in primers:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                # Use 'forward'/'reverse' (not '+'/'-') for PositionCache API
                positions_fwd = self.position_cache.get_positions(prefix, primer, "forward")
                positions_rev = self.position_cache.get_positions(prefix, primer, "reverse")
                positions = np.concatenate([positions_fwd, positions_rev])

                if len(positions) > 0:
                    # Thread the realistic per-primer reach, as Stage-1 set
                    # cover does. Omitting it left `extension_reach` at its
                    # default of 0, so this measured bin OCCUPANCY -- "does any
                    # site fall in this bin" -- and the answer depended on
                    # `bin_size` rather than on the polymerase.
                    #
                    # At the default 10 kb bin that reported 100% coverage for a
                    # set the shipped metric scores at 39%, and at a 1 kb bin the
                    # same set scored 13%. It also meant `_prune_background`
                    # checked its `min_coverage_threshold` floor against a
                    # different quantity from the one finally reported, so the
                    # guard against over-pruning was not measuring what it
                    # protected.
                    graph.add_primer_coverage(
                        primer,
                        positions,
                        prefix,
                        length,
                        extension_reach=self.coverage_reach,
                    )

        if len(graph.regions) == 0:
            return 0.0

        # Total bins at the SAME granularity the graph was built with. Dividing
        # bins counted at one size by a total computed at another reported
        # coverage above 1.0.
        total_bins = sum((length + bin_size - 1) // bin_size for length in self.fg_seq_lengths)

        coverage = len(graph.regions) / total_bins if total_bins > 0 else 0.0
        return coverage

    def _thermo_criteria(self):
        """Criteria for the Stage-0 thermodynamic screen.

        Two things used to be wrong here and both were invisible.

        The Tm bounds came from the polymerase preset alone, so a window
        widened in params.json was honoured by the filter step and then
        narrowed back by the optimizer. An explicitly configured bound now
        wins; the preset remains the default when nothing is configured.

        `na_conc=50.0, mg_conc=0.0` were hardcoded. Every hairpin, homodimer
        and heterodimer energy this screen rejects primers on was therefore
        computed for a reaction containing no magnesium -- the same
        zero-magnesium fault schema v2 corrected in the presets, and Mg2+ is
        the dominant term in the salt correction. The configured buffer is
        used instead.
        """
        from neoswga.core.thermodynamic_filter import ThermodynamicCriteria

        conditions = getattr(self, "conditions", None)
        if conditions is not None:
            na_conc = conditions.na_conc
            mg_conc = conditions.mg_conc
        else:
            from neoswga.core.parameter import default_mg_conc

            na_conc = 50.0
            mg_conc = default_mg_conc(self.polymerase)

        return ThermodynamicCriteria(
            min_tm=self.min_tm if self.min_tm is not None else self.poly_config.min_primer_tm,
            max_tm=self.max_tm if self.max_tm is not None else self.poly_config.max_primer_tm,
            target_tm=self.poly_config.reaction_temp + 5,
            na_conc=na_conc,
            mg_conc=mg_conc,
            max_homodimer_dg=-10.0,
            max_heterodimer_dg=-10.0,
            max_hairpin_dg=-3.0,
            min_gc=self.poly_config.min_gc,
            max_gc=self.poly_config.max_gc,
            reaction_temp=self.poly_config.reaction_temp,
            polymerase=self.polymerase,
        )

    def _thermo_filter_candidates(self, candidates: List[str], verbose: bool = True) -> List[str]:
        """
        Apply thermodynamic filtering based on polymerase requirements.

        Filters candidates by Tm range and GC content appropriate for
        the configured polymerase.
        """
        if verbose:
            logger.info("\n" + "-" * 80)
            logger.info(
                f"PRE-STAGE: Thermodynamic Filtering ({self.polymerase}, "
                f"{self.poly_config.reaction_temp}C)"
            )
            logger.info("-" * 80)

        try:
            from neoswga.core.thermodynamic_filter import ThermodynamicCriteria, ThermodynamicFilter

            criteria = self._thermo_criteria()

            thermo_filter = ThermodynamicFilter(criteria)
            filtered, stats = thermo_filter.filter_candidates(
                candidates,
                check_heterodimers=True,
                max_heterodimer_fraction=0.3,
            )

            if verbose:
                logger.info(f"Filtered: {len(filtered)}/{len(candidates)} passed")
                if stats.get("mean_tm") is not None:
                    logger.info(f"  Mean Tm: {stats['mean_tm']:.1f}C")
                if stats.get("mean_gc") is not None:
                    logger.info(f"  Mean GC: {stats['mean_gc']:.1%}")

            if len(filtered) == 0:
                logger.warning(
                    "No primers passed thermodynamic filtering, " "using unfiltered candidates"
                )
                return candidates

            return filtered

        except ImportError:
            logger.warning("Thermodynamic filter not available, skipping")
            return candidates

    def _prune_background(
        self, primers: List[str], target_size: int, verbose: bool = False
    ) -> Tuple[List[str], float, int]:
        """
        Greedy background pruning: remove primers with worst background/coverage ratio.

        Iteratively removes the primer whose removal causes the largest
        reduction in background binding relative to coverage loss. Stops
        when coverage drops below min_coverage_threshold or target_size
        is reached.

        Args:
            primers: Initial primer set from coverage stage
            target_size: Target number of primers after pruning
            verbose: Print removal details

        Returns:
            (pruned_primers, final_coverage, final_background_sites)
        """
        current_primers = list(primers)
        current_coverage = self._calculate_coverage(current_primers)

        if verbose:
            logger.info(f"  Background pruning: {len(current_primers)} -> {target_size} primers")
            logger.info(
                f"  Initial: coverage={current_coverage:.1%}, "
                f"background={self._count_background_sites(current_primers)} sites"
            )

        removed_count = 0

        while len(current_primers) > target_size:
            best_removal = None
            best_score = -np.inf

            for primer in current_primers:
                test_set = [p for p in current_primers if p != primer]
                test_coverage = self._calculate_coverage(test_set)
                coverage_loss = current_coverage - test_coverage

                if test_coverage < self.min_coverage_threshold:
                    continue

                primer_bg_sites = self._count_background_sites([primer])

                if coverage_loss > 0:
                    score = primer_bg_sites / coverage_loss * self.background_weight
                else:
                    score = primer_bg_sites * 1000

                if score > best_score:
                    best_score = score
                    best_removal = primer

            if best_removal is None:
                if verbose:
                    logger.info(f"  Stopping: coverage threshold reached")
                break

            current_primers.remove(best_removal)
            current_coverage = self._calculate_coverage(current_primers)
            removed_count += 1

            if verbose and removed_count % 2 == 0:
                logger.info(
                    f"    Removed {removed_count} primers, "
                    f"coverage={current_coverage:.1%}, "
                    f"background={self._count_background_sites(current_primers)} sites"
                )

        final_coverage = self._calculate_coverage(current_primers)
        final_bg_sites = self._count_background_sites(current_primers)

        return current_primers, final_coverage, final_bg_sites

    def _count_background_sites(self, primers: List[str]) -> int:
        """Count total background binding sites for primer set."""
        if not primers or not self.bg_prefixes:
            return 0

        total_sites = 0
        for primer in primers:
            for bg_prefix in self.bg_prefixes:
                fwd = self.position_cache.get_positions(bg_prefix, primer, "forward")
                rev = self.position_cache.get_positions(bg_prefix, primer, "reverse")
                total_sites += len(fwd) + len(rev)

        return total_sites


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

from neoswga.core.base_optimizer import (
    BaseOptimizer,
    OptimizationResult,
    OptimizationStatus,
    OptimizerConfig,
    PrimerSetMetrics,
)
from neoswga.core.optimizer_factory import OptimizerFactory


@OptimizerFactory.register("hybrid", aliases=["hybrid-optimizer", "two-stage"])
class HybridBaseOptimizer(BaseOptimizer):
    """
    Hybrid optimizer implementing BaseOptimizer interface.

    Combines dominating-set coverage (Stage 1) with network connectivity
    optimization (Stage 2). Stage 2 is conditions-aware: the inner
    NetworkOptimizer receives ReactionConditions and applies additive Tm
    corrections during refinement.
    """

    # Stage 2 (network refinement) applies ReactionConditions; this wrapper
    # is therefore additive-aware end-to-end even though Stage 1 is
    # coverage-only by design.
    ADDITIVE_AWARE = True

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        conditions=None,
        **kwargs,
    ):
        super().__init__(
            position_cache,
            fg_prefixes,
            fg_seq_lengths,
            bg_prefixes,
            bg_seq_lengths,
            config,
            conditions=conditions,
            # Forward the rest so background_profile / aggregate reach
            # BaseOptimizer. Dropping **kwargs here made a compositional
            # background silently inert on every optimizer.
            **kwargs,
        )
        self._hybrid = HybridOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            bin_size=kwargs.get("bin_size", 10000),
            max_extension=kwargs.get("max_extension"),
            # Stage-1 coverage selection uses the realistic reach resolved by
            # unified_optimizer (OptimizerConfig.extension_reach), so selection
            # and the scored metrics.fg_coverage share one coverage definition.
            coverage_reach=getattr(self.config, "extension_reach", None),
            polymerase=kwargs.get("polymerase", "phi29"),
            min_tm=getattr(self.config, "min_tm", None),
            max_tm=getattr(self.config, "max_tm", None),
            genome_gc_content=kwargs.get("genome_gc_content"),
            background_pruning=kwargs.get("background_pruning", False),
            background_weight=kwargs.get("background_weight", 2.0),
            min_coverage_threshold=kwargs.get("min_coverage_threshold", 0.95),
            conditions=conditions,
            mechanistic_weight=kwargs.get("mechanistic_weight", 0.0),
        )

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "hybrid"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return "Two-stage hybrid optimizer (coverage + network connectivity)"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs,
    ) -> OptimizationResult:
        """Run hybrid optimization."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        if self.config.verbose:
            logger.info(f"Running hybrid optimization: {len(candidates)} candidates")
            if fixed_primers:
                logger.info(f"  Fixed primers: {len(fixed_primers)}")

        try:
            result = self._hybrid.optimize(
                candidates=candidates,
                final_count=target,
                fixed_primers=fixed_primers,
                verbose=self.config.verbose,
                # `target` came from params.json or --auto-size, both
                # deliberate. --auto-size already accounts for the polymerase
                # through the mechanistic model, so rescaling here would
                # double-count it; an explicit count it would simply override.
                apply_polymerase_multiplier=False,
            )

            primers = result.primers
            metrics = self.compute_metrics(primers)

            return OptimizationResult(
                primers=tuple(primers),
                score=result.final_predicted_amplification,
                status=OptimizationStatus.SUCCESS if primers else OptimizationStatus.NO_CONVERGENCE,
                metrics=metrics,
                iterations=1,
                optimizer_name=self.name,
                message=f"Coverage: {result.final_coverage:.1%}, Connectivity: {result.final_connectivity:.2f}",
            )

        except Exception as e:
            logger.error(f"Hybrid optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))
