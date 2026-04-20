"""
Interval-based tiling optimizer for SWGA primer set selection.

Models each primer binding site as covering a fragment-length interval
(position to position + fragment_length), tracks coverage per base pair
using sorted interval merging, and uses strand-separate optimization
with adaptive declining thresholds.

This approach is inspired by the CoatSWGA interval tiling strategy
(Bailey Lab, 2025) and provides improved coverage uniformity compared
to fixed-bin models, particularly for longer primers.

Key features:
- Fragment-length-aware coverage intervals (biologically relevant for MDA)
- Strand-separate optimization (forward then reverse)
- Adaptive declining thresholds (avoids over-selecting for easy regions)
- Specificity-integrated scoring during selection
- Multi-seed parallel runs (escape local optima)

Usage:
    optimizer = TilingOptimizer(cache, fg_prefixes, fg_lengths)
    result = optimizer.optimize(candidates, target_size=10)
"""

import bisect
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .base_optimizer import (
    BaseOptimizer,
    OptimizationResult,
    OptimizationStatus,
    OptimizerConfig,
)
from .optimizer_factory import OptimizerFactory

logger = logging.getLogger(__name__)

# Default fragment lengths per polymerase (bp)
FRAGMENT_LENGTHS = {
    'phi29': 70000,
    'equiphi29': 50000,
    'bst': 10000,
    'klenow': 5000,
}


class IntervalCoverage:
    """Tracks genomic interval coverage using sorted interval merging.

    Maintains a list of non-overlapping, sorted intervals representing
    covered regions of a genome. Supports efficient addition of new
    intervals and marginal gain queries.
    """

    def __init__(self, genome_length: int):
        self.genome_length = genome_length
        self._intervals: List[Tuple[int, int]] = []
        self._covered_bases = 0

    @property
    def coverage_fraction(self) -> float:
        """Fraction of genome covered."""
        if self.genome_length == 0:
            return 0.0
        return self._covered_bases / self.genome_length

    @property
    def uncovered_bases(self) -> int:
        """Number of uncovered bases."""
        return self.genome_length - self._covered_bases

    def add_intervals(self, positions: np.ndarray, fragment_length: int) -> None:
        """Add coverage intervals from binding positions.

        Each position generates an interval [pos, pos + fragment_length),
        clamped to [0, genome_length).

        Args:
            positions: Array of binding positions.
            fragment_length: Length of fragment produced from each site.
        """
        if len(positions) == 0:
            return

        new_intervals = []
        for pos in positions:
            start = max(0, int(pos))
            end = min(self.genome_length, int(pos) + fragment_length)
            if start < end:
                new_intervals.append((start, end))

        if not new_intervals:
            return

        combined = self._intervals + new_intervals
        self._intervals = _merge_intervals(combined)
        self._covered_bases = sum(end - start for start, end in self._intervals)

    def marginal_gain(self, positions: np.ndarray, fragment_length: int) -> int:
        """Count new bases that would be covered without modifying state.

        Args:
            positions: Array of candidate binding positions.
            fragment_length: Length of fragment produced from each site.

        Returns:
            Number of previously uncovered bases that would be newly covered.
        """
        if len(positions) == 0:
            return 0

        new_intervals = []
        for pos in positions:
            start = max(0, int(pos))
            end = min(self.genome_length, int(pos) + fragment_length)
            if start < end:
                new_intervals.append((start, end))

        if not new_intervals:
            return 0

        merged_new = _merge_intervals(new_intervals)
        gain = 0
        for start, end in merged_new:
            gain += _uncovered_length(self._intervals, start, end)
        return gain


def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping intervals.

    Args:
        intervals: List of (start, end) tuples (unsorted ok).

    Returns:
        Sorted list of non-overlapping (start, end) tuples.
    """
    if not intervals:
        return []
    sorted_ivs = sorted(intervals)
    merged = [sorted_ivs[0]]
    for start, end in sorted_ivs[1:]:
        if start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _uncovered_length(
    existing: List[Tuple[int, int]], query_start: int, query_end: int
) -> int:
    """Count bases in [query_start, query_end) not covered by existing intervals.

    existing must be sorted and non-overlapping.
    """
    if not existing:
        return query_end - query_start

    uncovered = 0
    cursor = query_start

    # Find the first interval that could overlap with query_start
    idx = bisect.bisect_right([iv[0] for iv in existing], query_start) - 1
    if idx < 0:
        idx = 0

    for i in range(idx, len(existing)):
        iv_start, iv_end = existing[i]
        if iv_start >= query_end:
            break
        if iv_end <= cursor:
            continue
        # Gap before this interval
        if iv_start > cursor:
            uncovered += min(iv_start, query_end) - cursor
        cursor = max(cursor, iv_end)

    # Gap after last overlapping interval
    if cursor < query_end:
        uncovered += query_end - cursor

    return uncovered


@dataclass
class TilingConfig(OptimizerConfig):
    """Configuration for tiling optimizer."""
    fragment_length: int = 70000
    adaptive_thresholds: Tuple[float, ...] = (0.50, 0.25, 0.15, 0.10, 0.05, 0.0)
    strand_separate: bool = True
    n_seeds: int = 3
    specificity_weight: float = 1.0
    min_coverage: float = 0.95


@OptimizerFactory.register('tiling', aliases=['coverage-tiling', 'interval-tiling'])
class TilingOptimizer(BaseOptimizer):
    """Interval-based tiling optimizer with strand-separate coverage.

    Models each primer binding site as covering a fragment-length interval,
    then greedily selects primers that maximize new base-pair coverage while
    penalizing background binding. Adaptive declining thresholds relax the
    minimum gain requirement when high-gain primers are exhausted.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[TilingConfig] = None,
        conditions=None,
        **kwargs,
    ):
        super().__init__(
            position_cache,
            fg_prefixes,
            fg_seq_lengths,
            bg_prefixes,
            bg_seq_lengths,
            config or TilingConfig(),
            conditions=conditions,
        )
        self.tiling_config = (
            config if isinstance(config, TilingConfig) else TilingConfig()
        )

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "tiling"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return "Interval-based tiling optimizer with strand-separate coverage"

    @property
    def supports_background(self) -> bool:
        """Indicates this optimizer uses background genome data."""
        return True

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs,
    ) -> OptimizationResult:
        """Find optimal primer set using interval-based tiling.

        Args:
            candidates: Pool of candidate primer sequences.
            target_size: Desired number of primers (overrides config).
            fixed_primers: Primers that must be included in result.
            **kwargs: Additional parameters (unused).

        Returns:
            OptimizationResult with selected primers and metrics.
        """
        candidates = self._validate_candidates(candidates)
        max_primers = target_size or self.config.target_set_size
        fragment_length = self.tiling_config.fragment_length

        fixed = list(fixed_primers) if fixed_primers else []
        # Remove fixed primers from candidate pool
        fixed_set = set(fixed)
        pool = [p for p in candidates if p not in fixed_set]

        if self.config.verbose:
            logger.info(
                "Tiling optimization: %d candidates, target=%d, fragment=%d bp",
                len(candidates),
                max_primers,
                fragment_length,
            )

        # Pre-compute per-primer position data
        fg_positions = {}  # primer -> {strand -> np.ndarray of all fg positions}
        bg_counts = {}     # primer -> total bg binding sites

        all_primers = list(set(pool + fixed))
        for primer in all_primers:
            fg_positions[primer] = {}
            for strand in ('forward', 'reverse', 'both'):
                all_pos = []
                for prefix in self.fg_prefixes:
                    pos = self.get_primer_positions(primer, prefix, strand)
                    all_pos.append(pos)
                fg_positions[primer][strand] = (
                    np.concatenate(all_pos) if all_pos else np.array([], dtype=np.int64)
                )

            bg_total = 0
            for prefix in self.bg_prefixes:
                bg_pos = self.get_primer_positions(primer, prefix, 'both')
                bg_total += len(bg_pos)
            bg_counts[primer] = bg_total

        # Count fg sites for seed ranking
        fg_counts = {
            p: len(fg_positions[p].get('both', np.array([]))) for p in all_primers
        }

        # Run optimization
        if self.tiling_config.strand_separate:
            fwd_selected, fwd_coverage = self._tiling_single_strand(
                pool,
                'forward',
                fg_positions,
                bg_counts,
                fg_counts,
                fragment_length,
                max_primers,
                fixed,
            )
            rev_selected, rev_coverage = self._tiling_single_strand(
                pool,
                'reverse',
                fg_positions,
                bg_counts,
                fg_counts,
                fragment_length,
                max_primers,
                fixed,
            )
            # Merge: union of both strand results, deduplicated
            merged = list(dict.fromkeys(fixed + fwd_selected + rev_selected))
            # Trim to target size if over (keep fixed, then interleave fwd/rev)
            if len(merged) > max_primers:
                merged = self._trim_merged(
                    merged, fixed, fwd_selected, rev_selected, max_primers
                )
            selected = merged
            # Recompute coverage with 'both' strand for final metric
            final_cov = IntervalCoverage(self.fg_total_length)
            for primer in selected:
                final_cov.add_intervals(
                    fg_positions[primer]['both'], fragment_length
                )
            coverage_frac = final_cov.coverage_fraction
        else:
            selected, cov_obj = self._tiling_single_strand(
                pool,
                'both',
                fg_positions,
                bg_counts,
                fg_counts,
                fragment_length,
                max_primers,
                fixed,
            )
            selected = fixed + selected
            # Deduplicate preserving order
            seen = set()
            deduped = []
            for p in selected:
                if p not in seen:
                    seen.add(p)
                    deduped.append(p)
            selected = deduped[:max_primers]
            coverage_frac = cov_obj.coverage_fraction

        # Compute metrics via base class
        metrics = self.compute_metrics(selected)

        status = OptimizationStatus.SUCCESS
        if coverage_frac < self.tiling_config.min_coverage:
            status = OptimizationStatus.PARTIAL

        return OptimizationResult(
            primers=tuple(selected),
            score=coverage_frac,
            status=status,
            metrics=metrics,
            iterations=self.tiling_config.n_seeds,
            optimizer_name=self.name,
            message=(
                f"Coverage: {coverage_frac:.1%}, "
                f"{len(selected)} primers selected"
            ),
        )

    def _tiling_single_strand(
        self,
        pool: List[str],
        strand: str,
        fg_positions: Dict[str, Dict[str, np.ndarray]],
        bg_counts: Dict[str, int],
        fg_counts: Dict[str, int],
        fragment_length: int,
        max_primers: int,
        fixed: List[str],
    ) -> Tuple[List[str], IntervalCoverage]:
        """Run tiling optimization for a single strand direction.

        Returns:
            Tuple of (selected primers excluding fixed, IntervalCoverage).
        """
        n_seeds = self.tiling_config.n_seeds
        thresholds = self.tiling_config.adaptive_thresholds
        spec_weight = self.tiling_config.specificity_weight
        min_cov = self.tiling_config.min_coverage

        # Sort pool by fg binding count for seed selection
        ranked = sorted(pool, key=lambda p: fg_counts.get(p, 0), reverse=True)

        best_selected = None
        best_coverage = None

        for seed_idx in range(min(n_seeds, max(1, len(ranked)))):
            coverage = IntervalCoverage(self.fg_total_length)
            selected = []

            # Seed with fixed primers first
            for primer in fixed:
                pos = fg_positions[primer].get(strand, np.array([], dtype=np.int64))
                coverage.add_intervals(pos, fragment_length)

            # Seed primer from ranked candidates
            if seed_idx < len(ranked):
                seed_primer = ranked[seed_idx]
                selected.append(seed_primer)
                pos = fg_positions[seed_primer].get(
                    strand, np.array([], dtype=np.int64)
                )
                coverage.add_intervals(pos, fragment_length)

            budget = max_primers - len(fixed)
            threshold_idx = 0

            while (
                len(selected) < budget
                and coverage.coverage_fraction < min_cov
            ):
                remaining = coverage.uncovered_bases
                if remaining == 0:
                    break
                min_gain = int(
                    remaining * thresholds[threshold_idx]
                ) if threshold_idx < len(thresholds) else 0

                best_primer = None
                best_score = float('-inf')

                for primer in pool:
                    if primer in selected:
                        continue

                    pos = fg_positions[primer].get(
                        strand, np.array([], dtype=np.int64)
                    )
                    gain = coverage.marginal_gain(pos, fragment_length)
                    if gain < min_gain:
                        continue

                    fg_count = fg_counts.get(primer, 0)
                    bg_count = bg_counts.get(primer, 0)
                    specificity_penalty = bg_count / (fg_count + 1)
                    score = gain / (1.0 + spec_weight * specificity_penalty)

                    if score > best_score:
                        best_score = score
                        best_primer = primer

                if best_primer is None:
                    threshold_idx += 1
                    if threshold_idx >= len(thresholds):
                        break
                    continue

                selected.append(best_primer)
                pos = fg_positions[best_primer].get(
                    strand, np.array([], dtype=np.int64)
                )
                coverage.add_intervals(pos, fragment_length)
                threshold_idx = 0  # Reset after successful addition

            if (
                best_coverage is None
                or coverage.coverage_fraction > best_coverage.coverage_fraction
            ):
                best_selected = selected
                best_coverage = coverage

        if best_selected is None:
            return [], IntervalCoverage(self.fg_total_length)

        return best_selected, best_coverage

    @staticmethod
    def _trim_merged(
        merged: List[str],
        fixed: List[str],
        fwd: List[str],
        rev: List[str],
        max_primers: int,
    ) -> List[str]:
        """Trim merged primer list to max_primers, keeping fixed primers."""
        result = list(fixed)
        seen = set(fixed)
        # Interleave forward and reverse picks
        fi, ri = 0, 0
        while len(result) < max_primers:
            added = False
            while fi < len(fwd) and fwd[fi] in seen:
                fi += 1
            if fi < len(fwd) and len(result) < max_primers:
                result.append(fwd[fi])
                seen.add(fwd[fi])
                fi += 1
                added = True
            while ri < len(rev) and rev[ri] in seen:
                ri += 1
            if ri < len(rev) and len(result) < max_primers:
                result.append(rev[ri])
                seen.add(rev[ri])
                ri += 1
                added = True
            if not added:
                break
        return result
