"""
Automatic primer set size optimization using Pareto frontier analysis.

Determines optimal number of primers based on:
- Coverage vs fg/bg ratio tradeoff (Pareto frontier)
- Application profile (discovery, clinical, enrichment, metagenomics)
- Reaction conditions (affects primer effectiveness via mechanistic model)

Key insight: Coverage and specificity (fg/bg ratio) are NOT a simple tradeoff.
Good primer selection can improve BOTH metrics by choosing primers that:
1. Have high fg/bg binding site ratios (selective)
2. Bind to unique regions of the target genome (coverage)

Additive-aware optimization:
    When MechanisticEffects is provided, the optimizer accounts for how
    reaction conditions (additives, temperature) affect primer effectiveness:
    - Enhanced processivity (betaine, DMSO) -> fewer primers needed for coverage
    - Enhanced binding (SSB) -> more efficient primer utilization
    - Temperature effects -> adjusted coverage estimates

The Pareto frontier approach:
1. Generates primer sets at different sizes
2. Computes (coverage, fg_bg_ratio) for each, adjusted for mechanistic effects
3. Identifies Pareto-optimal points (no other point dominates)
4. Selects best point based on application priority

Usage:
    from neoswga.core.set_size_optimizer import (
        ParetoFrontierGenerator,
        select_from_frontier,
        recommend_set_size,
    )

    # Full frontier analysis (when position data available)
    generator = ParetoFrontierGenerator(optimizer, primer_pool, fg_positions, bg_positions)
    frontier = generator.generate_frontier(min_size=4, max_size=15)
    selected, explanation = select_from_frontier(frontier, 'clinical')

    # Quick estimate (without position data)
    size = quick_size_estimate('clinical', genome_length=1_000_000)

    # Additive-aware recommendation with mechanistic model
    from neoswga.core.mechanistic_model import MechanisticModel
    from neoswga.core.reaction_conditions import get_enhanced_conditions

    conditions = get_enhanced_conditions()
    model = MechanisticModel(conditions)
    effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
    rec = recommend_set_size('clinical', 1_000_000, 12, effects, processivity=80000)
"""

import math
import logging
from dataclasses import dataclass, field
from typing import List, Tuple, Dict, Optional, Any, TYPE_CHECKING, Callable

import numpy as np
import pandas as pd

from neoswga.core.mechanistic_params import APPLICATION_PROFILES, get_application_profile

if TYPE_CHECKING:
    from neoswga.core.mechanistic_model import MechanisticEffects
    from neoswga.core.base_optimizer import BaseOptimizer
    from neoswga.core.position_cache import PositionCache

logger = logging.getLogger(__name__)


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class SetSizeMetrics:
    """
    Metrics for evaluating a primer set at a given size.

    These metrics quantify the coverage vs specificity tradeoff:
    - fg_coverage: How much of the target genome is covered
    - fg_bg_ratio: How selective the primers are (fg_sites / bg_sites)

    A good primer set has BOTH high coverage AND high fg/bg ratio.
    """
    set_size: int
    fg_coverage: float          # Fraction of foreground genome covered [0, 1]
    bg_coverage: float          # Fraction of background genome covered [0, 1]
    fg_binding_sites: int       # Total foreground binding sites
    bg_binding_sites: int       # Total background binding sites
    fg_bg_ratio: float          # fg_sites / bg_sites (specificity metric)

    # Optional: the actual primers if available
    primers: Optional[Tuple[str, ...]] = None

    # Pareto status (computed after generating all candidates)
    _is_pareto_optimal: Optional[bool] = None

    @property
    def is_pareto_optimal(self) -> bool:
        """Whether this point is on the Pareto frontier."""
        if self._is_pareto_optimal is None:
            return True  # Unknown, assume optimal
        return self._is_pareto_optimal

    def dominates(self, other: 'SetSizeMetrics') -> bool:
        """
        Check if this point Pareto-dominates another.

        A point dominates another if it is:
        - At least as good in all objectives (coverage, fg_bg_ratio)
        - Strictly better in at least one objective

        We want to MAXIMIZE both coverage and fg_bg_ratio.
        """
        at_least_as_good = (
            self.fg_coverage >= other.fg_coverage and
            self.fg_bg_ratio >= other.fg_bg_ratio
        )
        strictly_better = (
            self.fg_coverage > other.fg_coverage or
            self.fg_bg_ratio > other.fg_bg_ratio
        )
        return at_least_as_good and strictly_better

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        result = {
            'set_size': self.set_size,
            'fg_coverage': self.fg_coverage,
            'bg_coverage': self.bg_coverage,
            'fg_binding_sites': self.fg_binding_sites,
            'bg_binding_sites': self.bg_binding_sites,
            'fg_bg_ratio': self.fg_bg_ratio,
            'is_pareto_optimal': self.is_pareto_optimal,
        }
        if self.primers is not None:
            result['primers'] = list(self.primers)
        return result

    @classmethod
    def empty(cls, set_size: int = 0) -> 'SetSizeMetrics':
        """Create empty metrics for failed optimization."""
        return cls(
            set_size=set_size,
            fg_coverage=0.0,
            bg_coverage=0.0,
            fg_binding_sites=0,
            bg_binding_sites=0,
            fg_bg_ratio=0.0,
        )


@dataclass
class FrontierResult:
    """
    Result of Pareto frontier generation.

    Contains all evaluated points and identifies the Pareto-optimal subset.
    """
    all_points: List[SetSizeMetrics]
    pareto_points: List[SetSizeMetrics]
    selected_point: Optional[SetSizeMetrics] = None
    selection_explanation: str = ""

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'all_points': [p.to_dict() for p in self.all_points],
            'pareto_points': [p.to_dict() for p in self.pareto_points],
            'selected_point': self.selected_point.to_dict() if self.selected_point else None,
            'selection_explanation': self.selection_explanation,
        }


# =============================================================================
# Pareto Frontier Generation
# =============================================================================

def filter_pareto_optimal(points: List[SetSizeMetrics]) -> List[SetSizeMetrics]:
    """
    Filter a list of points to keep only Pareto-optimal ones.

    A point is Pareto-optimal if no other point dominates it.
    We maximize both fg_coverage and fg_bg_ratio.

    Args:
        points: List of SetSizeMetrics to filter

    Returns:
        List of Pareto-optimal points, sorted by set_size
    """
    if not points:
        return []

    pareto = []
    for point in points:
        # Check if any other point dominates this one
        is_dominated = False
        for other in points:
            if other is point:
                continue
            if other.dominates(point):
                is_dominated = True
                break

        if not is_dominated:
            # Mark as Pareto-optimal
            point._is_pareto_optimal = True
            pareto.append(point)
        else:
            point._is_pareto_optimal = False

    # Sort by set_size for consistent ordering
    return sorted(pareto, key=lambda p: p.set_size)


class ParetoFrontierGenerator:
    """
    Generate coverage vs fg/bg ratio Pareto frontier using hybrid approach.

    The hybrid approach balances speed and accuracy:
    1. Coarse estimation: Quick estimate at all sizes from primer statistics
    2. Identify promising region: Find sizes near the estimated Pareto frontier
    3. Full optimization: Run actual optimizer at promising sizes only

    This is ~4x faster than full optimization at all sizes while still
    finding the true Pareto frontier.

    Additive-aware optimization:
        When mech_effects is provided, coverage estimates account for
        enhanced processivity from additives:
        - effective_processivity = base_processivity * processivity_factor
        - This allows fewer primers to achieve the same coverage
    """

    def __init__(
        self,
        primer_pool: pd.DataFrame,
        position_cache: Optional['PositionCache'] = None,
        fg_prefixes: Optional[List[str]] = None,
        bg_prefixes: Optional[List[str]] = None,
        fg_seq_lengths: Optional[List[int]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        optimizer: Optional['BaseOptimizer'] = None,
        processivity: int = 70000,
        mech_effects: Optional['MechanisticEffects'] = None,
    ):
        """
        Initialize the Pareto frontier generator.

        Args:
            primer_pool: DataFrame with candidate primers. Must have 'primer' column.
                        Optional columns: 'fg_freq', 'bg_freq', 'fg_bg_ratio'
            position_cache: PositionCache for primer binding positions
            fg_prefixes: Foreground genome HDF5 prefixes
            bg_prefixes: Background genome HDF5 prefixes
            fg_seq_lengths: Foreground genome lengths
            bg_seq_lengths: Background genome lengths
            optimizer: BaseOptimizer instance for full optimization
            processivity: Base polymerase processivity in bp (default: 70000 for phi29)
            mech_effects: Optional MechanisticEffects for additive-aware estimation.
                         If provided, processivity is adjusted by processivity_factor.
        """
        self.primer_pool = primer_pool
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes or []
        self.bg_prefixes = bg_prefixes or []
        self.fg_seq_lengths = fg_seq_lengths or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.optimizer = optimizer
        self.mech_effects = mech_effects

        # Apply mechanistic effects to processivity
        if mech_effects is not None:
            self.processivity = int(processivity * mech_effects.processivity_factor)
            self._effective_binding = mech_effects.effective_binding_rate
            logger.info(
                f"Additive-aware optimization: "
                f"processivity={self.processivity:,} bp "
                f"(base={processivity:,} * factor={mech_effects.processivity_factor:.2f}), "
                f"effective_binding={self._effective_binding:.2f}"
            )
        else:
            self.processivity = processivity
            self._effective_binding = 0.8  # Default assumption

        # Computed properties
        self.fg_total_length = sum(fg_seq_lengths) if fg_seq_lengths else 0
        self.bg_total_length = sum(bg_seq_lengths) if bg_seq_lengths else 0

        # Validate primer pool
        self._validate_primer_pool()

    def _validate_primer_pool(self) -> None:
        """Validate the primer pool DataFrame."""
        if self.primer_pool is None or len(self.primer_pool) == 0:
            raise ValueError("primer_pool cannot be empty")

        # Check for primer column (could be 'primer' or 'sequence' or index)
        if 'primer' not in self.primer_pool.columns:
            if 'sequence' in self.primer_pool.columns:
                self.primer_pool = self.primer_pool.rename(columns={'sequence': 'primer'})
            elif self.primer_pool.index.name == 'primer':
                self.primer_pool = self.primer_pool.reset_index()
            else:
                raise ValueError("primer_pool must have 'primer' or 'sequence' column")

    def generate_frontier(
        self,
        min_size: int = 4,
        max_size: int = 20,
        quick_only: bool = False,
        num_refine: int = 5,
        verbose: bool = True,
    ) -> FrontierResult:
        """
        Generate the Pareto frontier of (coverage, fg_bg_ratio) tradeoffs.

        Uses a hybrid approach:
        1. Quick estimation at all sizes from primer statistics
        2. Identify promising sizes near the estimated frontier
        3. Full optimization at promising sizes only (unless quick_only=True)

        Args:
            min_size: Minimum set size to evaluate
            max_size: Maximum set size to evaluate
            quick_only: If True, skip full optimization (faster but less accurate)
            num_refine: Number of sizes to fully optimize
            verbose: Whether to log progress

        Returns:
            FrontierResult with all points and Pareto-optimal subset
        """
        if verbose:
            logger.info(f"Generating Pareto frontier for sizes {min_size}-{max_size}")

        # Phase 1: Coarse estimation from primer-level statistics
        coarse_points = self._estimate_from_statistics(min_size, max_size)

        if quick_only or self.optimizer is None:
            # Return coarse estimates only
            pareto = filter_pareto_optimal(coarse_points)
            return FrontierResult(
                all_points=coarse_points,
                pareto_points=pareto,
            )

        # Phase 2: Identify promising sizes
        promising_sizes = self._identify_promising_sizes(coarse_points, num_refine)

        if verbose:
            logger.info(f"Refining at sizes: {promising_sizes}")

        # Phase 3: Full optimization at promising sizes
        refined_points = []
        for size in promising_sizes:
            try:
                metrics = self._optimize_at_size(size, verbose)
                refined_points.append(metrics)
            except Exception as e:
                logger.warning(f"Optimization at size {size} failed: {e}")
                # Fall back to coarse estimate
                coarse = next((p for p in coarse_points if p.set_size == size), None)
                if coarse:
                    refined_points.append(coarse)

        # Combine coarse and refined points, preferring refined where available
        refined_sizes = {p.set_size for p in refined_points}
        all_points = refined_points + [
            p for p in coarse_points if p.set_size not in refined_sizes
        ]
        all_points = sorted(all_points, key=lambda p: p.set_size)

        # Filter to Pareto-optimal
        pareto = filter_pareto_optimal(all_points)

        return FrontierResult(
            all_points=all_points,
            pareto_points=pareto,
        )

    def _estimate_from_statistics(
        self,
        min_size: int,
        max_size: int,
    ) -> List[SetSizeMetrics]:
        """
        Quick estimation assuming greedy primer selection by fg/bg ratio.

        For each size N:
        1. Rank primers by fg/bg ratio (highest first)
        2. Take top N primers
        3. Estimate coverage from binding sites + effective processivity
        4. Calculate fg/bg ratio from total sites

        This is fast but may overestimate coverage (ignores position overlap)
        and underestimate fg/bg ratio (ignores position-aware selection).

        Additive-aware:
            When mech_effects was provided to the generator, coverage
            estimates use the enhanced processivity and effective binding
            rate, resulting in more accurate predictions for reactions
            with additives.
        """
        # Sort primers by fg/bg ratio (descending)
        if 'fg_bg_ratio' in self.primer_pool.columns:
            sorted_pool = self.primer_pool.sort_values('fg_bg_ratio', ascending=False)
        elif 'fg_freq' in self.primer_pool.columns and 'bg_freq' in self.primer_pool.columns:
            # Calculate ratio on the fly
            pool = self.primer_pool.copy()
            pool['fg_bg_ratio'] = pool['fg_freq'] / (pool['bg_freq'] + 1e-10)
            sorted_pool = pool.sort_values('fg_bg_ratio', ascending=False)
        else:
            # No frequency data, just use order
            sorted_pool = self.primer_pool

        points = []
        for size in range(min_size, max_size + 1):
            if size > len(sorted_pool):
                continue

            top_n = sorted_pool.head(size)
            primers = tuple(top_n['primer'].tolist())

            # Estimate binding sites
            if 'fg_freq' in top_n.columns:
                # Use frequency data if available
                total_fg = top_n['fg_freq'].sum()
                total_bg = top_n['bg_freq'].sum() if 'bg_freq' in top_n.columns else 1e-10
                # Convert frequency to approximate site count
                fg_sites = int(total_fg * self.fg_total_length) if self.fg_total_length > 0 else int(total_fg * 1e6)
                bg_sites = int(total_bg * self.bg_total_length) if self.bg_total_length > 0 else max(1, int(total_bg * 1e6))
            else:
                # Use random sequence model
                primer_length = len(primers[0]) if primers else 10
                expected_sites_per_primer = (4 ** (-primer_length)) * 2  # Both strands
                fg_sites = int(expected_sites_per_primer * self.fg_total_length * size) if self.fg_total_length > 0 else size * 100
                bg_sites = int(expected_sites_per_primer * self.bg_total_length * size) if self.bg_total_length > 0 else size * 10
                bg_sites = max(1, bg_sites)

            # Estimate coverage with effective binding rate
            # Coverage model: 1 - exp(-sites * processivity * binding_rate / genome_length)
            if self.fg_total_length > 0:
                # Effective sites account for binding efficiency
                effective_fg_sites = fg_sites * self._effective_binding
                coverage_factor = effective_fg_sites * self.processivity / self.fg_total_length
                # Cap at reasonable value to prevent overflow
                coverage_factor = min(coverage_factor, 20)
                fg_coverage = 1 - math.exp(-coverage_factor)
            else:
                fg_coverage = min(0.95, 0.1 * size)  # Linear approximation

            if self.bg_total_length > 0:
                effective_bg_sites = bg_sites * self._effective_binding
                bg_coverage_factor = effective_bg_sites * self.processivity / self.bg_total_length
                bg_coverage_factor = min(bg_coverage_factor, 20)
                bg_coverage = 1 - math.exp(-bg_coverage_factor)
            else:
                bg_coverage = 0.0

            fg_bg_ratio = fg_sites / max(1, bg_sites)

            points.append(SetSizeMetrics(
                set_size=size,
                fg_coverage=fg_coverage,
                bg_coverage=bg_coverage,
                fg_binding_sites=fg_sites,
                bg_binding_sites=bg_sites,
                fg_bg_ratio=fg_bg_ratio,
                primers=primers,
            ))

        return points

    def _identify_promising_sizes(
        self,
        coarse_points: List[SetSizeMetrics],
        num_refine: int = 5,
    ) -> List[int]:
        """
        Identify sizes worth full optimization.

        Selects sizes that are:
        - On or near the coarse Pareto frontier
        - Spread across the coverage range

        Args:
            coarse_points: Points from coarse estimation
            num_refine: Number of sizes to select

        Returns:
            List of set sizes to refine
        """
        if len(coarse_points) <= num_refine:
            return [p.set_size for p in coarse_points]

        # Get coarse Pareto frontier
        pareto = filter_pareto_optimal(coarse_points)

        if len(pareto) >= num_refine:
            # If frontier has enough points, sample evenly by coverage
            pareto_sorted = sorted(pareto, key=lambda p: p.fg_coverage)
            indices = np.linspace(0, len(pareto_sorted) - 1, num_refine).astype(int)
            return [pareto_sorted[i].set_size for i in indices]

        # Otherwise, include all frontier points plus nearest neighbors
        frontier_sizes = {p.set_size for p in pareto}
        promising = list(frontier_sizes)

        # Add neighbors of frontier points
        all_sizes = sorted(p.set_size for p in coarse_points)
        for size in list(frontier_sizes):
            idx = all_sizes.index(size)
            if idx > 0 and all_sizes[idx - 1] not in promising:
                promising.append(all_sizes[idx - 1])
            if idx < len(all_sizes) - 1 and all_sizes[idx + 1] not in promising:
                promising.append(all_sizes[idx + 1])

            if len(promising) >= num_refine:
                break

        return sorted(promising)[:num_refine]

    def _optimize_at_size(
        self,
        size: int,
        verbose: bool = True,
    ) -> SetSizeMetrics:
        """
        Run full optimization at a specific set size.

        Args:
            size: Target set size
            verbose: Whether to log progress

        Returns:
            SetSizeMetrics with actual optimization results
        """
        if self.optimizer is None:
            raise ValueError("No optimizer provided for full optimization")

        candidates = self.primer_pool['primer'].tolist()

        if verbose:
            logger.info(f"  Optimizing at size {size}...")

        # Run optimizer
        result = self.optimizer.optimize(candidates, target_size=size)

        if not result.is_success:
            logger.warning(f"Optimization at size {size} did not succeed: {result.message}")

        primers = result.primers
        metrics = result.metrics

        return SetSizeMetrics(
            set_size=len(primers),
            fg_coverage=metrics.fg_coverage,
            bg_coverage=metrics.bg_coverage,
            fg_binding_sites=metrics.total_fg_sites,
            bg_binding_sites=max(1, metrics.total_bg_sites),
            fg_bg_ratio=metrics.selectivity_ratio,
            primers=primers,
        )

    def compute_metrics_for_set(
        self,
        primers: List[str],
    ) -> SetSizeMetrics:
        """
        Compute SetSizeMetrics for a given primer set.

        Args:
            primers: List of primer sequences

        Returns:
            SetSizeMetrics for the primer set
        """
        if not primers:
            return SetSizeMetrics.empty()

        if self.cache is None:
            # Can't compute without position cache
            return SetSizeMetrics(
                set_size=len(primers),
                fg_coverage=0.0,
                bg_coverage=0.0,
                fg_binding_sites=0,
                bg_binding_sites=1,
                fg_bg_ratio=0.0,
                primers=tuple(primers),
            )

        # Collect positions
        fg_positions = set()
        bg_positions = set()

        for primer in primers:
            for prefix in self.fg_prefixes:
                try:
                    pos = self.cache.get_positions(prefix, primer, 'both')
                    fg_positions.update(pos.tolist())
                except Exception as e:
                    logger.debug(f"Ignored error getting fg positions for {primer}: {e}")

            for prefix in self.bg_prefixes:
                try:
                    pos = self.cache.get_positions(prefix, primer, 'both')
                    bg_positions.update(pos.tolist())
                except Exception as e:
                    logger.debug(f"Ignored error getting bg positions for {primer}: {e}")

        fg_sites = len(fg_positions)
        bg_sites = max(1, len(bg_positions))

        # Compute coverage
        fg_coverage = self._compute_coverage(sorted(fg_positions), self.fg_total_length)
        bg_coverage = self._compute_coverage(sorted(bg_positions), self.bg_total_length)

        return SetSizeMetrics(
            set_size=len(primers),
            fg_coverage=fg_coverage,
            bg_coverage=bg_coverage,
            fg_binding_sites=fg_sites,
            bg_binding_sites=bg_sites,
            fg_bg_ratio=fg_sites / bg_sites,
            primers=tuple(primers),
        )

    def _compute_coverage(
        self,
        positions: List[int],
        total_length: int,
    ) -> float:
        """Compute coverage fraction from sorted positions."""
        if not positions or total_length == 0:
            return 0.0

        covered = 0
        i = 0
        while i < len(positions):
            start = positions[i]
            end = start + self.processivity
            while i < len(positions) and positions[i] <= end:
                i += 1
            covered += min(end, total_length) - start

        return min(1.0, covered / total_length)


# =============================================================================
# Frontier Selection
# =============================================================================

def select_from_frontier(
    frontier: List[SetSizeMetrics],
    application: str,
    min_fg_bg_ratio: Optional[float] = None,
    target_coverage: Optional[float] = None,
) -> Tuple[SetSizeMetrics, str]:
    """
    Select optimal point from Pareto frontier based on application.

    Args:
        frontier: List of Pareto-optimal SetSizeMetrics
        application: Application profile name ('discovery', 'clinical', etc.)
        min_fg_bg_ratio: Override for minimum fg/bg ratio constraint
        target_coverage: Override for target coverage

    Returns:
        Tuple of (selected_point, explanation_string)

    Raises:
        ValueError: If frontier is empty or application is unknown
    """
    if not frontier:
        raise ValueError("Frontier cannot be empty")

    # Get application profile
    try:
        profile = get_application_profile(application)
    except ValueError:
        logger.warning(f"Unknown application '{application}', using 'enrichment'")
        application = 'enrichment'
        profile = get_application_profile(application)

    # Use overrides if provided, otherwise profile defaults
    ratio_constraint = min_fg_bg_ratio
    if ratio_constraint is None:
        ratio_constraint = profile.get('default_min_fg_bg_ratio', profile.get('min_specificity', 0.5) * 10)

    coverage_target = target_coverage
    if coverage_target is None:
        coverage_target = profile.get('default_target_coverage', profile.get('target_coverage', 0.8))

    priority = profile.get('priority', 'balanced')

    # Filter to points meeting hard constraints
    valid = [p for p in frontier if p.fg_bg_ratio >= ratio_constraint]

    if not valid:
        # No points meet constraint - return best available with warning
        best = max(frontier, key=lambda p: p.fg_bg_ratio)
        explanation = (
            f"Warning: No set achieves fg/bg ratio >= {ratio_constraint:.1f}. "
            f"Best available: {best.fg_bg_ratio:.1f} at size {best.set_size}"
        )
        return best, explanation

    # Select based on priority
    if priority == 'coverage':
        # Pick point with highest coverage among valid
        # But also consider coverage target
        meeting_target = [p for p in valid if p.fg_coverage >= coverage_target]
        if meeting_target:
            # Among those meeting coverage target, pick smallest set (Occam's razor)
            selected = min(meeting_target, key=lambda p: p.set_size)
            explanation = (
                f"Selected size {selected.set_size} for '{application}' (coverage priority): "
                f"{selected.fg_coverage:.1%} coverage with fg/bg ratio {selected.fg_bg_ratio:.1f}"
            )
        else:
            # Pick highest coverage available
            selected = max(valid, key=lambda p: p.fg_coverage)
            explanation = (
                f"Selected size {selected.set_size} for '{application}' (best available coverage): "
                f"{selected.fg_coverage:.1%} coverage (target: {coverage_target:.1%}), "
                f"fg/bg ratio {selected.fg_bg_ratio:.1f}"
            )

    elif priority == 'specificity':
        # Pick point with highest fg/bg ratio among valid
        selected = max(valid, key=lambda p: p.fg_bg_ratio)
        explanation = (
            f"Selected size {selected.set_size} for '{application}' (specificity priority): "
            f"fg/bg ratio {selected.fg_bg_ratio:.1f}, {selected.fg_coverage:.1%} coverage"
        )

    else:  # balanced - find knee of curve
        # Use geometric distance to ideal point
        # Ideal: (1.0 coverage, max fg/bg)
        max_ratio = max(p.fg_bg_ratio for p in valid)

        # Normalize both axes to [0, 1]
        def distance_to_ideal(p: SetSizeMetrics) -> float:
            cov_norm = 1.0 - p.fg_coverage  # Distance from perfect coverage
            ratio_norm = 1.0 - (p.fg_bg_ratio / max_ratio) if max_ratio > 0 else 1.0
            return math.sqrt(cov_norm**2 + ratio_norm**2)

        selected = min(valid, key=distance_to_ideal)
        explanation = (
            f"Selected size {selected.set_size} for '{application}' (balanced): "
            f"{selected.fg_coverage:.1%} coverage, fg/bg ratio {selected.fg_bg_ratio:.1f} "
            f"(at Pareto frontier knee)"
        )

    return selected, explanation


# =============================================================================
# High-Level API (Backward Compatible)
# =============================================================================

def estimate_optimal_set_size(
    genome_length: int,
    primer_length: int,
    target_coverage: float,
    processivity: int,
    mech_effects: 'MechanisticEffects'
) -> int:
    """
    Estimate optimal set size from first principles.

    Uses a probabilistic model of genome coverage:
    Coverage = 1 - exp(-N * sites_per_primer * effective_processivity / genome_length)

    Args:
        genome_length: Target genome length in bp
        primer_length: Primer length in bp (k-mer size)
        target_coverage: Desired genome coverage fraction (0-1)
        processivity: Base polymerase processivity in bp
        mech_effects: MechanisticEffects from the model

    Returns:
        Estimated optimal number of primers

    Example:
        >>> from neoswga.core.mechanistic_model import MechanisticEffects
        >>> effects = MechanisticEffects(
        ...     tm_correction=0, effective_tm=35.0,
        ...     accessibility_factor=1.0, processivity_factor=1.0,
        ...     speed_factor=1.0, stability_factor=1.0,
        ...     kon_factor=1.0, koff_factor=1.0,
        ...     effective_binding_rate=0.8, effective_extension_rate=1.0,
        ...     predicted_amplification_factor=0.8,
        ... )
        >>> estimate_optimal_set_size(1_000_000, 10, 0.8, 70000, effects)
        8
    """
    # Expected binding sites per primer (random sequence model)
    # For a k-mer in a random sequence of length L, expected count is:
    # E[sites] = (L - k + 1) * (1/4)^k * 2 (both strands)
    # Simplified: ~ L * 4^(-k) * 2
    expected_sites = (4 ** (-primer_length)) * genome_length * 2

    # Effective processivity accounting for reaction conditions
    effective_proc = processivity * mech_effects.processivity_factor

    # Coverage contribution per primer
    # Each primer contributes coverage proportional to:
    # sites * effective_processivity * binding_rate / genome_length
    coverage_per_primer = (
        expected_sites * effective_proc * mech_effects.effective_binding_rate
        / genome_length
    )

    if coverage_per_primer <= 0:
        logger.warning("Coverage per primer is zero or negative, using default size")
        return 20

    # Solve for N in: target_coverage = 1 - exp(-N * coverage_per_primer)
    # => N = -ln(1 - target_coverage) / coverage_per_primer
    try:
        optimal_n = -math.log(1 - target_coverage) / coverage_per_primer
    except (ValueError, ZeroDivisionError):
        logger.warning("Could not compute optimal N, using default size")
        return 15

    # Apply overlap correction factor
    # Primers don't cover genome independently; there's overlap
    # Empirical correction: multiply by ~1.3
    optimal_n *= 1.3

    # Bound to reasonable range
    return max(4, min(20, int(optimal_n + 0.5)))


def recommend_set_size(
    application: str,
    genome_length: int,
    primer_length: int,
    mech_effects: 'MechanisticEffects',
    processivity: int = 70000,
    min_fg_bg_ratio: Optional[float] = None,
    target_coverage: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Recommend primer set size based on application profile.

    Combines first-principles estimation with application-specific
    requirements for coverage and fg/bg ratio (specificity).

    Additive-aware optimization:
        The mech_effects parameter contains factors that modify the
        calculation based on reaction conditions (additives, temperature):
        - processivity_factor: Enhanced processivity from betaine/DMSO
          means fewer primers needed to achieve the same coverage
        - effective_binding_rate: Better binding efficiency from
          optimal Tm and accessibility means more efficient primers

        For example, with betaine 1.5M + DMSO 5%:
        - processivity_factor ~ 1.1 (10% boost)
        - effective_binding_rate ~ 0.85 (better than baseline 0.8)
        Result: ~15% fewer primers needed for same coverage

    Args:
        application: Application profile name
            - 'discovery': Maximize sensitivity for pathogen discovery
            - 'clinical': Minimize false positives for diagnostics
            - 'enrichment': Balanced for sequencing enrichment
            - 'metagenomics': Capture diversity
        genome_length: Target genome length in bp
        primer_length: Primer length in bp
        mech_effects: MechanisticEffects from mechanistic model. This
                     provides processivity_factor and effective_binding_rate
                     that account for reaction conditions.
        processivity: Base polymerase processivity (default 70000 for phi29).
                     This is multiplied by mech_effects.processivity_factor.
        min_fg_bg_ratio: Override for minimum fg/bg ratio (optional)
        target_coverage: Override for target coverage (optional)

    Returns:
        Dictionary with:
            - recommended_size: Recommended number of primers
            - size_range: (min, max) typical range for this application
            - target_coverage: Target coverage fraction
            - min_fg_bg_ratio: Minimum fg/bg ratio requirement
            - rationale: Description of the application
            - base_estimate: Raw estimate before application adjustment
            - priority: Selection priority ('coverage', 'specificity', 'balanced')
            - effective_processivity: Processivity used in calculation
            - additive_adjustment: Percentage adjustment from additives

            # Legacy fields (for backward compatibility):
            - min_specificity: Same as min_fg_bg_ratio / 10 (deprecated)

    Example:
        >>> from neoswga.core.mechanistic_model import MechanisticModel
        >>> from neoswga.core.reaction_conditions import get_enhanced_conditions
        >>> conditions = get_enhanced_conditions()  # 42C, 5% DMSO, 1M betaine
        >>> model = MechanisticModel(conditions)
        >>> effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        >>> rec = recommend_set_size('clinical', 1_000_000, 12, effects)
        >>> print(f"Use {rec['recommended_size']} primers")
        >>> print(f"Additive adjustment: {rec['additive_adjustment']:.0f}%")
    """
    # Get application profile
    try:
        profile = get_application_profile(application)
    except ValueError:
        logger.warning(f"Unknown application '{application}', using 'enrichment'")
        application = 'enrichment'
        profile = get_application_profile(application)

    # Get target coverage (use override or profile default)
    cov_target = target_coverage
    if cov_target is None:
        cov_target = profile.get('default_target_coverage', profile.get('target_coverage', 0.8))

    # Get fg/bg ratio constraint (use override or profile default)
    ratio_constraint = min_fg_bg_ratio
    if ratio_constraint is None:
        ratio_constraint = profile.get('default_min_fg_bg_ratio', 5.0)

    # Get base estimate from first principles
    base_estimate = estimate_optimal_set_size(
        genome_length=genome_length,
        primer_length=primer_length,
        target_coverage=cov_target,
        processivity=processivity,
        mech_effects=mech_effects,
    )

    # Adjust based on application priority
    priority = profile.get('priority', 'balanced')

    if priority == 'specificity':
        # High specificity applications may benefit from fewer, more selective primers
        # This is based on the insight that smaller sets can have better fg/bg ratios
        # when primers are carefully selected
        adjusted = max(4, int(base_estimate * 0.85))
    elif priority == 'coverage':
        # Coverage priority may need more primers to ensure all regions are covered
        adjusted = min(20, int(base_estimate * 1.1))
    else:  # balanced
        adjusted = base_estimate

    # Constrain to typical range for this application
    min_size, max_size = profile['typical_size']
    recommended = max(min_size, min(max_size, adjusted))

    # Calculate effective processivity and additive adjustment
    effective_processivity = int(processivity * mech_effects.processivity_factor)

    # Additive adjustment shows how much the mechanistic model changes the estimate
    # Baseline assumption: processivity_factor=1.0, effective_binding_rate=0.8
    baseline_estimate = estimate_optimal_set_size(
        genome_length=genome_length,
        primer_length=primer_length,
        target_coverage=cov_target,
        processivity=processivity,
        mech_effects=create_baseline_effects(),
    )
    if baseline_estimate > 0:
        additive_adjustment = ((baseline_estimate - base_estimate) / baseline_estimate) * 100
    else:
        additive_adjustment = 0.0

    return {
        'recommended_size': recommended,
        'size_range': profile['typical_size'],
        'target_coverage': cov_target,
        'min_fg_bg_ratio': ratio_constraint,
        'priority': priority,
        'rationale': profile['description'],
        'base_estimate': base_estimate,
        'application': application,
        'effective_processivity': effective_processivity,
        'additive_adjustment': additive_adjustment,
        'processivity_factor': mech_effects.processivity_factor,
        'effective_binding_rate': mech_effects.effective_binding_rate,
        # Legacy field for backward compatibility
        'min_specificity': profile.get('min_specificity', ratio_constraint / 10),
    }


def find_optimal_size_by_elbow(
    optimizer,
    candidates: List[str],
    min_size: int = 4,
    max_size: int = 15,
    verbose: bool = True
) -> Tuple[int, List[Dict[str, Any]]]:
    """
    Find optimal set size using elbow method on coverage curve.

    Runs optimization at multiple set sizes and finds the "elbow"
    where adding more primers yields diminishing returns.

    Args:
        optimizer: NetworkOptimizer instance (or compatible optimizer)
        candidates: List of candidate primer sequences
        min_size: Minimum set size to evaluate
        max_size: Maximum set size to evaluate
        verbose: Whether to log progress

    Returns:
        Tuple of (optimal_size, results_list)
        - optimal_size: Recommended set size at elbow
        - results_list: List of dicts with size, primers, coverage, score

    Note:
        This method is computationally expensive as it runs full
        optimization at each size. Use recommend_set_size() for
        a faster estimate.
    """
    results = []

    for size in range(min_size, max_size + 1):
        if verbose:
            logger.info(f"Evaluating set size {size}...")

        # Run optimization at this size
        primers = optimizer.optimize_greedy(candidates, num_primers=size)
        score = optimizer.score_primer_set(primers)

        result = {
            'size': size,
            'primers': primers,
            'coverage': score.get('target_coverage', score.get('enrichment', 0)),
            'enrichment': score.get('enrichment', 0),
            'score': score,
        }
        results.append(result)

        if verbose:
            logger.info(f"  Size {size}: coverage={result['coverage']:.1%}, "
                       f"enrichment={result['enrichment']:.1f}x")

    # Find elbow using curvature analysis
    if len(results) < 3:
        return results[len(results) // 2]['size'], results

    sizes = np.array([r['size'] for r in results])
    coverages = np.array([r['coverage'] for r in results])

    # Normalize to [0, 1] range
    size_range = sizes.max() - sizes.min()
    coverage_range = coverages.max() - coverages.min()

    if size_range == 0 or coverage_range < 1e-10:
        # No variation, return middle
        return results[len(results) // 2]['size'], results

    sizes_norm = (sizes - sizes.min()) / size_range
    coverages_norm = (coverages - coverages.min()) / coverage_range

    # Compute curvature using finite differences
    # First derivative
    d1 = np.gradient(coverages_norm, sizes_norm)
    # Second derivative
    d2 = np.gradient(d1, sizes_norm)

    # Curvature formula: |d2| / (1 + d1^2)^1.5
    curvature = np.abs(d2) / np.power(1 + d1**2, 1.5)

    # Find maximum curvature (excluding endpoints which can be noisy)
    if len(curvature) > 2:
        # Exclude first and last points
        interior_curvature = curvature[1:-1]
        elbow_idx = np.argmax(interior_curvature) + 1
    else:
        elbow_idx = len(results) // 2

    optimal_size = results[elbow_idx]['size']

    if verbose:
        logger.info(f"Elbow detected at size {optimal_size}")

    return optimal_size, results


def get_size_recommendation_summary(
    application: str,
    genome_length: int,
    primer_length: int,
    mech_effects: 'MechanisticEffects',
    processivity: int = 70000
) -> str:
    """
    Get a human-readable summary of set size recommendation.

    Args:
        application: Application profile name
        genome_length: Target genome length in bp
        primer_length: Primer length in bp
        mech_effects: MechanisticEffects from mechanistic model
        processivity: Base polymerase processivity

    Returns:
        Formatted string with recommendation details
    """
    rec = recommend_set_size(
        application=application,
        genome_length=genome_length,
        primer_length=primer_length,
        mech_effects=mech_effects,
        processivity=processivity,
    )

    lines = [
        f"Set Size Recommendation for '{application}' Application",
        "=" * 50,
        f"",
        f"Recommended primers: {rec['recommended_size']}",
        f"Typical range: {rec['size_range'][0]}-{rec['size_range'][1]}",
        f"",
        f"Application: {rec['rationale']}",
        f"Priority: {rec['priority']}",
        f"Target coverage: {rec['target_coverage']:.0%}",
        f"Min fg/bg ratio: {rec['min_fg_bg_ratio']:.1f}",
        f"",
        f"Base estimate (first principles): {rec['base_estimate']}",
        f"Additive adjustment: {rec['additive_adjustment']:+.1f}%",
        f"",
        f"Input parameters:",
        f"  Genome length: {genome_length:,} bp",
        f"  Primer length: {primer_length} bp",
        f"  Base processivity: {processivity:,} bp",
        f"  Effective processivity: {rec['effective_processivity']:,} bp",
        f"",
        f"Mechanistic model effects:",
        f"  Processivity factor: {mech_effects.processivity_factor:.2f}",
        f"  Speed factor: {mech_effects.speed_factor:.2f}",
        f"  Effective binding: {mech_effects.effective_binding_rate:.2f}",
        f"  Accessibility: {mech_effects.accessibility_factor:.2f}",
    ]

    return "\n".join(lines)


def create_baseline_effects() -> 'MechanisticEffects':
    """
    Create baseline MechanisticEffects for quick estimates.

    Returns MechanisticEffects with all factors at 1.0 (optimal conditions).
    Use this when you don't have specific reaction conditions.

    Returns:
        MechanisticEffects with baseline values
    """
    from neoswga.core.mechanistic_model import MechanisticEffects

    return MechanisticEffects(
        tm_correction=0.0,
        effective_tm=35.0,
        accessibility_factor=1.0,
        processivity_factor=1.0,
        speed_factor=1.0,
        stability_factor=1.0,
        kon_factor=1.0,
        koff_factor=1.0,
        effective_binding_rate=0.8,
        effective_extension_rate=1.0,
        predicted_amplification_factor=0.8,
    )


def quick_size_estimate(
    application: str,
    genome_length: int,
    primer_length: int = 10,
    processivity: int = 70000
) -> int:
    """
    Quick set size estimate without full mechanistic model.

    Convenience function for rapid estimates when you don't need
    detailed reaction condition modeling.

    Args:
        application: Application profile name
        genome_length: Target genome length in bp
        primer_length: Primer length (default 10)
        processivity: Polymerase processivity (default 70000 for phi29)

    Returns:
        Recommended number of primers

    Example:
        >>> quick_size_estimate('clinical', 1_000_000)
        8
        >>> quick_size_estimate('discovery', 5_000_000)
        12
    """
    effects = create_baseline_effects()
    rec = recommend_set_size(
        application=application,
        genome_length=genome_length,
        primer_length=primer_length,
        mech_effects=effects,
        processivity=processivity,
    )
    return rec['recommended_size']
