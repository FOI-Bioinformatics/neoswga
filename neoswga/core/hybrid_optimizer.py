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

from neoswga.core.coverage_counter import CoverageCounter
from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.hybrid_thermo_screen import ThermoScreenMixin
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

# Share of the Stage-2 ranking given to host binding, and applied ONLY when the
# optimizer was built background-aware (`background_pruning=True`, which is what
# `--optimization-method background-aware` and `expand-primers` against a host
# genome both set). At 0.0 the expression below reduces exactly to the two-term
# one, so plain `hybrid` is untouched.
#
# Stage 1.5 pruning alone does not make a panel specific, because Stage 2 is
# what chooses the panel. Measured on a 40-candidate expansion over a 300 kb
# synthetic target and host, enabling pruning moved delivered host binding from
# 32 sites to 45: pruning only shrank the pool Stage 2 drew from, and a
# background-blind choice over a smaller pool is not a better choice. The
# background has to be one of the axes the deciding stage ranks on.
#
# The remaining 0.75 is split between connectivity and coverage in the existing
# 50/50 ratio, so adding this term re-weights the two rather than displacing
# either.
#
# Measured on the three GC-tier designs against hg38, host sites in the
# delivered panel before and after this term, with `hybrid` unchanged in every
# case (identical panels):
#
#     design            bg-aware before   after   coverage before -> after
#     low_saureus n=24             1674    1361      0.4014 -> 0.3830
#     low_saureus n=36             2381    2003      0.5348 -> 0.5144
#     mid_ecoli   n=24              746     669      0.4867 -> 0.4809
#     mid_ecoli   n=36              928     787      0.5965 -> 0.5951
#     high_mtb    n=24              215     199      0.7411 -> 0.7246
#     high_mtb    n=36              336     218      0.8358 -> 0.8052
#
# Host binding falls 7-35% and coverage falls 0.1-3.1 points. That trade is the
# reason the term is gated: it is what `background-aware` is for, and it is not
# what a `hybrid` user asked for.
#
# At n=12 nothing moves on any of the three pools. Stage 1 returns 18-19
# primers there, so almost every one is carrying coverage no other primer
# supplies, and the coverage term decides every removal on its own. The
# background only gets a say once the set has enough redundancy for the
# coverage differences to narrow.
_STAGE2_BACKGROUND_WEIGHT = 0.25


def _rank_removal_candidates(
    candidates: List[Tuple[str, float, int, int]], background_aware: bool
) -> Tuple[Optional[str], float]:
    """Pick which primer Stage 2 should drop this step.

    Each candidate is `(primer, network_score, unique_bins, bg_sites)`. Rank
    within the step rather than adding raw values: algebraic connectivity, a bin
    count and a host-site count share no scale, and a fixed coefficient between
    them would be arbitrary and input-dependent. Normalising per step keeps the
    trade-off meaningful whatever the magnitudes happen to be.

    Returns `(None, -inf)` for an empty candidate list, which the caller reads
    as "nothing left that may be removed".
    """
    if not candidates:
        return None, -float("inf")

    net_values = [c[1] for c in candidates]
    cost_values = [c[2] for c in candidates]
    bg_values = [c[3] for c in candidates]
    net_lo, net_hi = min(net_values), max(net_values)
    cost_lo, cost_hi = min(cost_values), max(cost_values)
    bg_lo, bg_hi = min(bg_values), max(bg_values)
    net_span = (net_hi - net_lo) or 1.0
    cost_span = (cost_hi - cost_lo) or 1.0
    bg_span = (bg_hi - bg_lo) or 1.0

    # The two original weights keep their ratio to each other and share whatever
    # the background term leaves. With `bg_weight == 0` this is the previous
    # two-term expression exactly, which is what keeps plain `hybrid` still.
    bg_weight = _STAGE2_BACKGROUND_WEIGHT if background_aware else 0.0
    net_weight = (1.0 - _STAGE2_COVERAGE_WEIGHT) * (1.0 - bg_weight)
    cov_weight = _STAGE2_COVERAGE_WEIGHT * (1.0 - bg_weight)

    best_primer = None
    best_combined = -float("inf")
    for primer, network_score, unique_bins, bg_sites in candidates:
        norm_net = (network_score - net_lo) / net_span
        # Cheap to remove == loses few unique bins == score 1.
        norm_keep = 1.0 - (unique_bins - cost_lo) / cost_span
        # Binds the host most == best to remove == score 1.
        norm_bg = (bg_sites - bg_lo) / bg_span
        combined = net_weight * norm_net + cov_weight * norm_keep + bg_weight * norm_bg
        if combined > best_combined:
            best_combined = combined
            best_primer = primer

    return best_primer, best_combined


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
  Coverage: {self.final_coverage:.1%} (estimated, binned)
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


class HybridOptimizer(ThermoScreenMixin):
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
        # How much coverage background pruning may give up, measured from the
        # coverage stage 1 handed it. See _prune_background for why the floor
        # has to be relative: an absolute 0.95 is unreachable on a real target,
        # so the stage removed nothing at all.
        min_coverage_drop: float = 0.05,
        # Opt back in to the absolute reading of `min_coverage_threshold`, for
        # a caller that wants a hard bound rather than a budget.
        absolute_coverage_floor: bool = False,
        conditions=None,
        mechanistic_weight: float = 0.0,
        # Stage-2 scoring parameters. These are handed straight to the inner
        # NetworkOptimizer, which is the object that actually scores candidates
        # during refinement, so the defaults here are ITS defaults
        # (network_optimizer.py) rather than anything new: a caller that omits
        # them gets exactly the previous behaviour.
        #
        # They were absent from this signature entirely, which is why
        # `--application` could report "selection weights applied: tm=0.30,
        # uniformity=0.05, dimer_penalty=0.45" while the object doing the
        # selecting used 0.0, 0.0 and 0.0, and why a params.json max_dimer_bp of
        # 3 became 4 -- the looser threshold -- inside the default method.
        reaction_temp: Optional[float] = None,
        tm_weight: float = 0.0,
        dimer_penalty: float = 0.0,
        max_dimer_bp: Optional[int] = None,
        template_gc: float = 0.5,
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
            reaction_temp: Reaction temperature (C) for Stage-2 Tm scoring.
            tm_weight: Weight for the Tm term in Stage-2 selection.
            dimer_penalty: Weight for the dimer term in Stage-2 selection. This
                is the only dimer consideration in Stage-2, so a value of 0.0
                means refinement weighs dimers not at all.
            max_dimer_bp: Longest complementary stretch tolerated between two
                primers before the dimer term applies. Resolved from this
                argument, then `parameter.max_dimer_bp`, then 3 -- the same
                threshold is then used to gate the Stage-0 thermodynamic
                pre-screen, so the pipeline no longer computes one criterion
                and reports the delivered pool against another.
            template_gc: Template GC fraction used by the mechanistic term.
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
        self.min_coverage_drop = min_coverage_drop
        self._absolute_coverage_floor = absolute_coverage_floor

        # The configured dimer limit, resolved once and passed to all three
        # stages below. Before this it reached dimer.is_dimer and nothing else:
        # the Stage-0 screen ran on a hardcoded free-energy threshold no user
        # could set, and the delivered pool was reported against max_dimer_bp,
        # which nothing had enforced.
        self.max_dimer_bp = self._resolve_max_dimer_bp(max_dimer_bp)

        # Initialize both optimizers
        self.dominating_optimizer = DominatingSetOptimizer(
            cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bin_size=bin_size,
            # Stage-1 coverage uses the realistic per-primer reach so the
            # selection objective matches how the result is scored.
            extension_reach=self.coverage_reach,
            # Omitting this let Stage 1 re-resolve its own threshold, so a
            # config supplying 4 selected under 3 here and reported against 4.
            max_dimer_bp=self.max_dimer_bp,
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
            # The Stage-2 scoring parameters. Omitting these left the object
            # that performs refinement at its own defaults, so `--application`
            # and a configured `max_dimer_bp` reached the adapter and stopped
            # there. `network` has always passed them; the default method did
            # not. See tests/test_optimizer_config_reaches_optimizers.py.
            reaction_temp=reaction_temp,
            tm_weight=tm_weight,
            dimer_penalty=dimer_penalty,
            max_dimer_bp=self.max_dimer_bp,
            template_gc=template_gc,
        )
        # Retain for introspection / rescoring hooks.
        self.conditions = conditions

        # collect_alternative_sets re-enters optimize() once per alternative
        # set over what is nearly the same pool, so without this the screen ran
        # up to nine times. Measured at 1372 us a pair before Task 2's
        # pre-screen, that was the dominant cost of the pipeline. The reuse is
        # an approximation, not an identity: see
        # ThermoScreenMixin._thermo_filter_with_cache.
        self._thermo_filter_cache = None

        self._log_configuration(polymerase, bin_size, background_pruning, background_weight)

    def _log_configuration(self, polymerase, bin_size, background_pruning, background_weight):
        """What this optimizer was constructed with, at INFO."""
        logger.info("Hybrid optimizer initialized")
        logger.info(f"  Polymerase: {polymerase}")
        logger.info(f"  Bin size: {bin_size:,} bp")
        logger.info(f"  Max extension: {self.max_extension:,} bp")
        if self.poly_config.thermo_filter:
            logger.info(f"  Thermo-filtering: enabled ({self.poly_config.reaction_temp}C)")
        if background_pruning:
            logger.info(f"  Background pruning: enabled (weight={background_weight})")

    @staticmethod
    def _resolve_max_dimer_bp(max_dimer_bp: Optional[int]) -> int:
        """Resolve the dimer limit: constructor argument, then params.json, then 3.

        `isinstance(..., int) and not isinstance(..., bool)` rather than a
        truthiness check, because several test modules replace the `parameter`
        module with a mock whose every attribute is truthy -- `getattr(...) or
        default` would silently reconfigure the threshold from such a mock.
        """
        if max_dimer_bp is not None:
            return int(max_dimer_bp)

        from neoswga.core import parameter as _parameter

        configured = getattr(_parameter, "max_dimer_bp", None)
        if isinstance(configured, int) and not isinstance(configured, bool):
            return int(configured)
        return 3

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
            candidates_filtered = self._thermo_filter_with_cache(
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
                fixed_primers=fixed_primers,
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
            logger.info(f"Final coverage: {final_coverage_result:.1%} (estimated, binned)")
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
        bin_counts = self._bin_occupancy(bins_by_primer)

        # Host binding, the third axis, and the reason this stage stopped being
        # background-blind. Counted once per primer here rather than inside the
        # removal loop: a primer's host load does not depend on which other
        # primers are still in the set, so recounting it every step would be the
        # same number at O(steps) times the cost.
        background_aware = bool(self.background_pruning and self.bg_prefixes)
        bg_by_primer = (
            {p: self._count_background_sites([p]) for p in current_primers}
            if background_aware
            else {}
        )

        while len(current_primers) > target_count:
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
                unique_bins = sum(
                    1 for b in bins_by_primer.get(primer, ()) if bin_counts[self._bin_key(b)] == 1
                )
                candidates.append((primer, network_score, unique_bins, bg_by_primer.get(primer, 0)))

            best_to_remove, _best_score = _rank_removal_candidates(candidates, background_aware)

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

        ONE graph for all the primers, not one per primer. `BipartiteGraph`
        already dedupes regions by coordinate through its `_region_lookup`, and
        already records which primers cover each region, so a shared graph gives
        both this mapping and the occupancy counts with no extra structure --
        the same structure the greedy path uses, which is why greedy was never
        affected by the bug this fixes.

        Per-primer graphs looked equivalent and were not. `CoverageRegion` hashes
        on `(chromosome, start, end)` but the dataclass `__eq__` also compares
        `covered_by`, so the same bin arriving from two different graphs hashed
        equal and compared unequal. Every bin then landed in the occupancy
        counter as its own key with a count of 1, "uniquely covered" became
        "covered at all", and the removal ranking preferred discarding the primer
        with the FEWEST bins -- the sole coverer of a region.
        """
        from neoswga.core.dominating_set_optimizer import BipartiteGraph, coverage_bin_size

        bin_size = coverage_bin_size(self.bin_size, self.coverage_reach)
        graph = BipartiteGraph(bin_size=bin_size)

        for primer in primers:
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

        return {primer: set(graph.primer_to_regions.get(primer, ())) for primer in primers}

    @staticmethod
    def _bin_key(region):
        """The identity of a coverage bin: where it is, not who covers it."""
        return (region.chromosome, region.start, region.end)

    @staticmethod
    def _bin_occupancy(bins_by_primer: Dict[str, set]) -> Dict[object, int]:
        """How many of the given primers cover each bin.

        Keyed on the region's coordinates rather than on the region object, so
        the count cannot be split by an unequal-but-equally-hashing key. That is
        belt and braces given `_coverage_bins_by_primer` now shares one graph,
        and it is what makes the count correct for any caller that does not.
        """
        counts: Dict[object, int] = Counter()
        for owned in bins_by_primer.values():
            counts.update({HybridOptimizer._bin_key(r) for r in owned})
        return counts

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

    def _build_coverage_counter(self, primers: List[str]) -> CoverageCounter:
        """A `CoverageCounter` over the same bins `_calculate_coverage` counts.

        Delegates the graph-building to `_coverage_bins_by_primer`, which
        already builds ONE shared `BipartiteGraph` for the whole set rather
        than one per primer -- that method's docstring records what a
        per-primer graph does to `CoverageRegion` equality. Bins are keyed here
        by `_bin_key` (coordinates) rather than by the region object, matching
        `_bin_occupancy`'s convention, for the same reason.

        The denominator is `_total_coverage_bins`, the same genome-bin total
        `_calculate_coverage` divides by -- NOT the graph's own covered-region
        count, which would always equal the numerator (every region the graph
        holds came from these primers) and read a constant 1.0.
        """
        bins_by_primer = self._coverage_bins_by_primer(primers)
        counter = CoverageCounter(total_bins=self._total_coverage_bins())
        for primer, bins in bins_by_primer.items():
            counter.add(primer, {self._bin_key(b) for b in bins})
        return counter

    def _prune_background(
        self,
        primers: List[str],
        target_size: int,
        verbose: bool = False,
        fixed_primers: Optional[List[str]] = None,
    ) -> Tuple[List[str], float, int]:
        """
        Greedy background pruning: remove primers with worst background/coverage ratio.

        Iteratively removes the primer whose removal causes the largest
        reduction in background binding relative to coverage loss. Stops when
        coverage would fall more than `min_coverage_drop` below the coverage
        this stage STARTED with, or when target_size is reached.

        The floor used to be absolute: a removal was rejected whenever the
        remaining coverage fell below `min_coverage_threshold`, default 0.95.
        The quantity compared is binned coverage at realistic reach, which on a
        real target is well under 0.95 -- the file's own note on
        `_calculate_coverage` records a set reading 39% under the corrected
        measure where the older bin-occupancy measure read 100%. So on any
        normal run every candidate removal was rejected on the first pass,
        `best_removal` stayed None, and the loop broke immediately:

            stage-1 coverage 0.99 -> 12 primers to 9, background 7800 -> 4500
            stage-1 coverage 0.94 -> 12 primers to 12, background unchanged
            stage-1 coverage 0.39 -> 12 primers to 12, background unchanged

        which made `--optimization-method background-aware` degenerate to plain
        hybrid on every run whose coverage sat at or below the floor, while the
        documentation advertised a 10-20x background reduction.

        A RELATIVE floor is what the stage needs: it exists to stop pruning
        eating the coverage stage 1 just bought, and "how much did we give up"
        is the question that asks. The absolute reading is still available by
        passing `min_coverage_threshold`, for a caller that wants a hard bound.

        Args:
            primers: Initial primer set from coverage stage
            target_size: Target number of primers after pruning
            verbose: Print removal details

        Returns:
            (pruned_primers, final_coverage, final_background_sites)
        """
        # Deduped up front: the counter is keyed by primer sequence, so a
        # duplicate entry would desync it from `current_primers` the moment
        # one copy was removed -- `list.remove` drops a single occurrence
        # while `counter.remove` drops that primer's bins entirely. The old
        # rebuild-per-step code had no such gap (it always saw the true
        # remaining set), so the two structures must agree from here on.
        current_primers = list(dict.fromkeys(primers))
        counter = self._build_coverage_counter(current_primers)
        current_coverage = counter.covered_fraction()
        # Both floors apply; the relative one is what normally binds.
        floor = max(
            self.min_coverage_threshold if self._absolute_coverage_floor else 0.0,
            current_coverage - self.min_coverage_drop,
        )

        if verbose:
            logger.info(f"  Background pruning: {len(current_primers)} -> {target_size} primers")
            logger.info(
                f"  Initial: coverage={current_coverage:.1%}, "
                f"background={self._count_background_sites(current_primers)} sites"
            )

        # Fixed primers are the caller's existing panel. `optimize` hands this
        # stage `fixed_primers + newly_selected`, and the ranking below is
        # purely background per unit of coverage, so a fixed primer that binds
        # the host heavily is the first thing it would remove. That is worse
        # than ignoring the background: `expand-primers` would silently drop
        # part of the panel it was asked to extend. `_network_refine` already
        # takes the same argument for the same reason.
        protected = set(fixed_primers or [])

        removed_count = 0

        while len(current_primers) > target_size:
            if len(current_primers) <= len(protected):
                break
            best_removal = None
            best_score = -np.inf

            for primer in current_primers:
                if primer in protected:
                    continue
                lost_bins = counter.loss_if_removed(primer)
                test_coverage = (
                    (counter.covered_count() - lost_bins) / counter.total_bins
                    if counter.total_bins
                    else 0.0
                )
                coverage_loss = current_coverage - test_coverage

                if test_coverage < floor:
                    continue

                primer_bg_sites = self._count_background_sites([primer])

                # Rank on background removed per unit of coverage given up, and
                # break ties on the cheaper removal.
                #
                # This was two disjoint branches -- `bg / loss * weight` when a
                # removal cost coverage and `bg * 1000` when it did not -- with
                # two consequences. `background_weight` was a positive constant
                # multiplier INSIDE one branch, so it could not change that
                # branch's argmax and could not make pruning "more aggressive"
                # as documented; all it did was move the crossover against the
                # other branch, favouring removals that cost coverage. And when
                # every candidate carried the same host load, which includes
                # every candidate carrying none, all scores tied and the strict
                # `>` kept the first primer in list order, so the choice was
                # made by input order with the coverage cost never consulted.
                #
                # One expression removes both. `background_weight` now scales
                # the background term against the coverage term, which is the
                # trade-off it is named for, and the coverage cost is always
                # part of the score, so a tie on host load is broken by keeping
                # the primer that covers more.
                score = (self.background_weight * primer_bg_sites) - (coverage_loss * 1000.0)

                if score > best_score:
                    best_score = score
                    best_removal = primer

            if best_removal is None:
                if verbose:
                    logger.info(f"  Stopping: coverage threshold reached")
                break

            current_primers.remove(best_removal)
            counter.remove(best_removal)
            current_coverage = counter.covered_fraction()
            removed_count += 1

            if verbose and removed_count % 2 == 0:
                logger.info(
                    f"    Removed {removed_count} primers, "
                    f"coverage={current_coverage:.1%}, "
                    f"background={self._count_background_sites(current_primers)} sites"
                )

        # current_coverage already reflects this exact set: it is updated from
        # the counter after every removal (and, if the loop never ran, from the
        # counter built at entry), so a fresh rebuild here would just recompute
        # the same number -- the thing this method exists to stop doing.
        final_coverage = current_coverage
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

    @staticmethod
    def _coverage_message(final_coverage: float, connectivity: float) -> str:
        """The `message` stored in step4_improved_df_summary.json.

        `final_coverage` is the binned progress figure from
        `HybridOptimizer._calculate_coverage`, not `metrics.fg_coverage`, which
        the same summary carries and which is computed base by base. The two
        disagree by construction, so the label says which this one is rather
        than leaving a reader to reconcile 55.2% against 50.9%.
        """
        return (
            f"Coverage: {final_coverage:.1%} (estimated, binned; "
            f"fg_coverage is the measured value), "
            f"Connectivity: {connectivity:.2f}"
        )

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
            min_coverage_drop=kwargs.get("min_coverage_drop", 0.05),
            absolute_coverage_floor=kwargs.get("absolute_coverage_floor", False),
            conditions=conditions,
            mechanistic_weight=kwargs.get("mechanistic_weight", 0.0),
            # Selection weights and the dimer threshold. `uniformity_weight`
            # was accepted by HybridOptimizer and simply not passed here; the
            # other four were not on its signature at all until this change.
            # The CLI always sends the weights, resolved from --application by
            # unified_optimizer._resolve_selection_weights, and max_dimer_bp
            # comes off the config the way `network` already reads it.
            uniformity_weight=kwargs.get("uniformity_weight", 0.0),
            reaction_temp=kwargs.get("reaction_temp"),
            tm_weight=kwargs.get("tm_weight", 0.0),
            dimer_penalty=kwargs.get("dimer_penalty", 0.0),
            max_dimer_bp=getattr(self.config, "max_dimer_bp", 4),
            template_gc=kwargs.get("template_gc", 0.5),
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
                message=self._coverage_message(result.final_coverage, result.final_connectivity),
            )

        except Exception as e:
            logger.error(f"Hybrid optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))
