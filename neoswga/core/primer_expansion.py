"""
Primer set expansion for iterative SWGA design.

Allows users to provide existing primers that have been experimentally
validated and design additional primers to fill coverage gaps.

This supports iterative wet-lab workflows where:
1. Initial primer set is synthesized and tested
2. Some primers work well, others fail
3. User wants to design additional primers to improve coverage

Usage:
    from neoswga.core.primer_expansion import PrimerExpander

    expander = PrimerExpander(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
    )

    # Identify coverage gaps in current primer set
    gaps = expander.identify_gaps(
        primers=['ATCGATCG', 'GCTAGCTA'],
        min_gap_size=10000
    )

    # Expand primer set with additional primers
    result = expander.expand(
        candidates=candidates_df['primer'].tolist(),
        fixed_primers=['ATCGATCG', 'GCTAGCTA'],
        failed_primers=['BADPRIMER'],  # Exclude these
        target_new=6,
        optimization_method='hybrid'
    )
"""

import dataclasses
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from neoswga.core.deficit_objective import (
    DeficitObjective,
    dilate_intervals,
    recovered_deficit,
)
from neoswga.core.design_context import design_context_from_params
from neoswga.core.dominating_set_optimizer import coverage_bin_size as _coverage_bin_size

logger = logging.getLogger(__name__)


@dataclass
class CoverageGap:
    """A region of the genome not covered by current primers."""

    chromosome: str
    start: int
    end: int
    size: int

    def __post_init__(self):
        if self.size == 0:
            self.size = self.end - self.start


def merge_gap_intervals(gaps: List["CoverageGap"]) -> List["CoverageGap"]:
    """Union overlapping/adjacent gap intervals per chromosome.

    Used to combine in-silico gaps (from binding positions) with real
    sequencing-depth gaps (from a BAM) into one non-overlapping set. Gaps are
    grouped by ``chromosome``; within each group, intervals are sorted by start
    and merged when they overlap or touch. Returns the merged gaps sorted
    largest-first.
    """
    by_chrom: dict = {}
    for g in gaps:
        by_chrom.setdefault(g.chromosome, []).append((g.start, g.end))

    merged: List[CoverageGap] = []
    for chrom, intervals in by_chrom.items():
        intervals.sort()
        cur_start, cur_end = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_end:  # overlap or touch
                cur_end = max(cur_end, e)
            else:
                merged.append(CoverageGap(chrom, cur_start, cur_end, cur_end - cur_start))
                cur_start, cur_end = s, e
        merged.append(CoverageGap(chrom, cur_start, cur_end, cur_end - cur_start))

    merged.sort(key=lambda g: g.size, reverse=True)
    return merged


@dataclass
class ExpansionInput:
    """
    Input for primer set expansion.

    Attributes:
        fixed_primers: Primers that must be kept (already validated)
        failed_primers: Primers that failed and should be excluded
        coverage_gaps: Optional list of specific regions needing coverage
    """

    fixed_primers: List[str] = field(default_factory=list)
    failed_primers: List[str] = field(default_factory=list)
    coverage_gaps: List[CoverageGap] = field(default_factory=list)


@dataclass
class ExpansionResult:
    """
    Result of primer set expansion.

    Attributes:
        new_primers: Newly designed primers
        combined_set: Fixed + new primers
        fixed_primers: Original fixed primers
        coverage_before: Coverage with fixed primers only
        coverage_after: Coverage with combined set
        gap_coverage: Fraction of gap regions now covered
        predicted_improvement: Estimated improvement in amplification
        gaps_remaining: Number of significant gaps still present
        optimization_method: Method used for optimization
    """

    new_primers: List[str]
    combined_set: List[str]
    fixed_primers: List[str]
    coverage_before: float
    coverage_after: float
    gap_coverage: float
    predicted_improvement: float
    gaps_remaining: int
    optimization_method: str
    message: str = ""
    # The read-selection rules any BAM depth behind this result was read
    # under. `None` when no BAM was used. A breadth figure means nothing
    # without the rule that produced it, and two runs under different rules
    # are not comparable. See `core/depth_policy.py`.
    depth_policy: Optional[Dict] = None
    stage_history: tuple = ()

    @property
    def n_new(self) -> int:
        """Number of new primers designed."""
        return len(self.new_primers)

    @property
    def n_fixed(self) -> int:
        """Number of fixed primers."""
        return len(self.fixed_primers)

    @property
    def n_total(self) -> int:
        """Total primers in combined set."""
        return len(self.combined_set)

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "new_primers": self.new_primers,
            "combined_set": self.combined_set,
            "fixed_primers": self.fixed_primers,
            "n_new": self.n_new,
            "n_fixed": self.n_fixed,
            "n_total": self.n_total,
            "coverage_before": self.coverage_before,
            "coverage_after": self.coverage_after,
            "gap_coverage": self.gap_coverage,
            "predicted_improvement": self.predicted_improvement,
            "gaps_remaining": self.gaps_remaining,
            "optimization_method": self.optimization_method,
            "stage_history": list(self.stage_history),
            "message": self.message,
            "depth_policy": self.depth_policy,
        }

    def __str__(self) -> str:
        lines = [
            "Primer Expansion Result:",
            f"  Fixed primers: {self.n_fixed}",
            f"  New primers: {self.n_new}",
            f"  Total set: {self.n_total}",
            f"  Coverage: {self.coverage_before:.1%} -> {self.coverage_after:.1%}",
            f"  Gap coverage: {self.gap_coverage:.1%}",
            f"  Predicted improvement: {self.predicted_improvement:.1f}x",
            f"  Gaps remaining: {self.gaps_remaining}",
        ]
        if self.message:
            lines.append(f"  Note: {self.message}")
        return "\n".join(lines)


class PrimerExpander:
    """
    Expands existing primer sets with additional primers.

    Supports iterative design workflows where users have validated
    some primers experimentally and want to add more.
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
        context=None,
        candidate_source=None,
    ):
        """
        Initialize primer expander.

        Args:
            position_cache: PositionCache with primer binding positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers (optional)
            bg_seq_lengths: Background genome lengths (optional)
            bin_size: Bin size for coverage analysis (bp)
            max_extension: Deprecated alias for `coverage_reach`. It always
                meant the per-primer coverage reach despite the name, and was
                then handed to `HybridOptimizer(max_extension=...)`, which is
                the amplification-NETWORK reach -- a different quantity,
                answered by single-molecule processivity.
            coverage_reach: Per-primer extension reach in bp for the coverage
                figure and for Stage-1 set cover. Defaults to 3000, the
                realistic reach for phi29 in a dense SWGA design, matching
                base_optimizer.compute_metrics.
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.bin_size = bin_size
        # A `DesignContext` resolves the reach, chemistry, dimer limits and
        # circularity from one params file, so `expand-primers` designs under
        # the same conditions `plan-pool` does. Without one the old defaults
        # stand, because three callers and any library caller pass none; that
        # default of 3 kb whatever the polymerase is finding F7's
        # configuration half.
        self.context = context
        self.candidate_source = candidate_source
        self.conditions = getattr(context, "conditions", None)
        self.max_dimer_bp = getattr(context, "max_dimer_bp", None)
        self.fg_circular = bool(getattr(context, "fg_circular", False))
        self.coverage_reach = (
            coverage_reach or getattr(context, "coverage_reach", None) or max_extension or 3000
        )
        # Retained so existing callers reading the attribute still see the
        # value they set; it has only ever been the coverage reach.
        self.max_extension = self.coverage_reach
        # A bin larger than the reach lets a primer covering part of it claim
        # all of it, which is what `coverage_bin_size` exists to stop.
        self.coverage_bin_size = _coverage_bin_size(bin_size, self.coverage_reach)
        self.total_length = sum(fg_seq_lengths)

    def identify_gaps(
        self,
        primers: List[str],
        min_gap_size: int = 10000,
        extra_gaps: Optional[List[CoverageGap]] = None,
        merge: bool = True,
    ) -> List[CoverageGap]:
        """
        Identify coverage gaps in a primer set.

        Args:
            primers: Current primer sequences
            min_gap_size: Minimum gap size to report (bp)
            extra_gaps: Optional externally-derived gaps (e.g. low-depth
                regions from a mapped BAM, via ``bam_coverage.bam_gaps``) in
                the SAME coordinate space (chromosome == fg prefix). When
                provided and ``merge`` is True they are unioned with the
                in-silico gaps.
            merge: If True (default), union ``extra_gaps`` with the in-silico
                gaps via :func:`merge_gap_intervals`.

        Returns:
            List of CoverageGap objects sorted by size (largest first)
        """
        if not primers:
            # Entire genome is a gap
            gaps = []
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths, strict=True):
                gaps.append(CoverageGap(chromosome=prefix, start=0, end=length, size=length))
            if extra_gaps and merge:
                gaps = merge_gap_intervals(gaps + list(extra_gaps))
            return gaps

        # Collect all binding positions
        all_positions = {}  # prefix -> sorted positions

        for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths, strict=True):
            positions = []
            for primer in primers:
                pos = self.cache.get_positions(prefix, primer, "both")
                positions.extend(pos.tolist())
            all_positions[prefix] = sorted(set(positions))

        # Find gaps
        gaps = []

        for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths, strict=True):
            positions = all_positions.get(prefix, [])

            if not positions:
                # Entire chromosome is a gap
                gaps.append(CoverageGap(chromosome=prefix, start=0, end=length, size=length))
                continue

            # Check gap at start
            if positions[0] > min_gap_size:
                gaps.append(
                    CoverageGap(chromosome=prefix, start=0, end=positions[0], size=positions[0])
                )

            # Check internal gaps
            for i in range(1, len(positions)):
                gap_size = positions[i] - positions[i - 1]
                if gap_size > min_gap_size:
                    gaps.append(
                        CoverageGap(
                            chromosome=prefix,
                            start=positions[i - 1],
                            end=positions[i],
                            size=gap_size,
                        )
                    )

            # Check gap at end
            if length - positions[-1] > min_gap_size:
                gaps.append(
                    CoverageGap(
                        chromosome=prefix,
                        start=positions[-1],
                        end=length,
                        size=length - positions[-1],
                    )
                )

        # Merge in externally-derived gaps (e.g. BAM low-depth regions).
        if extra_gaps and merge:
            gaps = merge_gap_intervals(gaps + list(extra_gaps))
            gaps = [g for g in gaps if g.size >= min_gap_size]

        # Sort by size (largest first)
        gaps.sort(key=lambda g: g.size, reverse=True)

        return gaps

    def _calculate_coverage(self, primers: List[str]) -> float:
        """Calculate genome coverage fraction for a primer set."""
        if not primers:
            return 0.0

        from neoswga.core.dominating_set_optimizer import BipartiteGraph

        bin_size = self.coverage_bin_size
        graph = BipartiteGraph(bin_size=bin_size)

        for primer in primers:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths, strict=True):
                positions = self.cache.get_positions(prefix, primer, "both")
                if len(positions) > 0:
                    # Without `extension_reach` this marked only the bin each
                    # site sits in, so the result was bin OCCUPANCY and its
                    # value came from `bin_size` rather than the polymerase: a
                    # 10 kb bin let one site claim 10 kb of coverage.
                    graph.add_primer_coverage(
                        primer,
                        positions,
                        prefix,
                        length,
                        extension_reach=self.coverage_reach,
                    )

        if not graph.regions:
            return 0.0

        total_bins = sum((length + bin_size - 1) // bin_size for length in self.fg_seq_lengths)

        return min(1.0, len(graph.regions) / total_bins) if total_bins > 0 else 0.0

    def expand(
        self,
        candidates: List[str],
        fixed_primers: List[str],
        failed_primers: Optional[List[str]] = None,
        target_new: int = 6,
        optimization_method: str = "hybrid",
        verbose: bool = True,
        target_gaps: Optional[List[CoverageGap]] = None,
    ) -> ExpansionResult:
        """
        Expand primer set with additional primers.

        Args:
            candidates: Pool of candidate primers to choose from
            fixed_primers: Primers that must be included (validated)
            failed_primers: Primers to exclude (failed in wet lab)
            target_new: Number of new primers to select
            optimization_method: 'hybrid', 'dominating-set', or 'network'
            verbose: Print progress
            target_gaps: Optional list of gaps (in-silico and/or BAM-derived)
                to focus on. When provided, the candidate pool is HARD-filtered
                to primers with at least one binding site inside a gap, with a
                fallback to the full pool if too few remain to fill ``target_new``.

        Returns:
            ExpansionResult with new and combined primer sets
        """
        failed_primers = failed_primers or []
        self._excluded_primers = tuple(str(p).upper() for p in failed_primers)
        fixed_set = set(p.upper() for p in fixed_primers)
        failed_set = set(p.upper() for p in failed_primers)

        # Filter candidates
        candidates_filtered = [
            c for c in candidates if c.upper() not in fixed_set and c.upper() not in failed_set
        ]

        # Budget prescreen: candidates whose WINDOW can reach a target gap,
        # which is the gap intervals dilated by one reach. It used to require
        # the binding SITE inside the gap and never consulted the reach, so a
        # candidate binding just outside a 50 kb gap whose window would blanket
        # it was discarded while one binding at the gap's last base and
        # extending away was kept.
        #
        # There is no fallback. Abandoning the prescreen when it left fewer
        # than `target_new` candidates is what made a gap-targeted run silently
        # stop targeting gaps; a short list is a finding, not a reason to
        # search somewhere else.
        if target_gaps:
            in_gap = self._filter_candidates_to_gaps(candidates_filtered, target_gaps)
            if verbose:
                logger.info(
                    f"Gap-reachable pool: {len(in_gap)}/{len(candidates_filtered)} "
                    f"candidates can reach one of {len(target_gaps)} target gap(s) "
                    f"within the {self.coverage_reach} bp reach"
                )
            if len(in_gap) < target_new:
                logger.warning(
                    "Only %d candidate(s) can reach a target gap, fewer than the "
                    "%d requested. Selecting from those rather than widening to "
                    "the full pool: a pool that cannot reach the gaps is a "
                    "finding about the pool, and searching elsewhere would "
                    "deliver primers that do not address what was asked.",
                    len(in_gap),
                    target_new,
                )
            candidates_filtered = in_gap

        if verbose:
            logger.info("=" * 60)
            logger.info("PRIMER SET EXPANSION")
            logger.info("=" * 60)
            logger.info(f"Fixed primers: {len(fixed_primers)}")
            logger.info(f"Failed primers (excluded): {len(failed_primers)}")
            logger.info(f"Candidate pool: {len(candidates_filtered)}")
            logger.info(f"Target new primers: {target_new}")
            logger.info(f"Optimization method: {optimization_method}")

        # Calculate initial coverage
        coverage_before = self._calculate_coverage(list(fixed_set))

        if verbose:
            logger.info(f"\nCurrent coverage: {coverage_before:.1%}")

        # Identify gaps before expansion
        gaps_before = self.identify_gaps(list(fixed_set), min_gap_size=self.bin_size)

        if verbose:
            logger.info(f"Significant gaps: {len(gaps_before)}")
            if gaps_before and len(gaps_before) <= 5:
                for gap in gaps_before[:5]:
                    logger.info(f"  {gap.chromosome}: {gap.start:,}-{gap.end:,} ({gap.size:,} bp)")

        # Select optimization method.
        #
        # `background-aware` used to fall through to the `else` below and run
        # plain hybrid with a warning, so a user who asked for it got exactly
        # the run they were trying to avoid. In this codebase background-aware
        # IS hybrid with pruning switched on -- that is what the registered
        # optimizer of that name does -- so honouring it is one argument.
        #
        # The remaining methods the CLI offers are not implemented here, and an
        # unsupported method now RAISES. Substituting a different algorithm and
        # logging about it is how a flag comes to mean nothing; the same class
        # of defect as an inert params.json key.
        if optimization_method in ("hybrid", "two-stage"):
            result = self._expand_hybrid(
                candidates_filtered,
                fixed_primers,
                target_new,
                verbose,
                target_gaps=target_gaps,
            )
        elif optimization_method in ("background-aware", "bg-aware", "clinical"):
            result = self._expand_hybrid(
                candidates_filtered,
                fixed_primers,
                target_new,
                verbose,
                background_pruning=True,
                target_gaps=target_gaps,
            )
        elif optimization_method in ("dominating-set", "dominating_set", "ds"):
            result = self._expand_dominating_set(
                candidates_filtered, fixed_primers, target_new, verbose, target_gaps=target_gaps
            )
        else:
            raise ValueError(
                f"expand-primers does not implement optimization method "
                f"{optimization_method!r}. Supported: hybrid (the default), "
                f"background-aware, dominating-set. It used to run hybrid here "
                f"and log a warning, which meant asking for a different method "
                f"changed nothing."
            )

        # Calculate coverage improvement
        combined_set = list(fixed_set) + result["new_primers"]
        coverage_after = self._calculate_coverage(combined_set)

        # Count remaining gaps
        gaps_after = self.identify_gaps(combined_set, min_gap_size=self.bin_size)

        # Fraction of the targeted deficit BASES recovered. It used to be
        # `1 - len(after)/len(before)`, a count of gaps, which goes NEGATIVE
        # when one long gap splits into two -- progress reported as
        # regression. See `core/deficit_objective.py`.
        gap_coverage = self._recovered_deficit_fraction(target_gaps or gaps_before, combined_set)

        # Estimate improvement factor
        if coverage_before > 0:
            predicted_improvement = coverage_after / coverage_before
        else:
            predicted_improvement = float("inf") if coverage_after > 0 else 1.0

        if verbose:
            logger.info("\n" + "=" * 60)
            logger.info("EXPANSION COMPLETE")
            logger.info("=" * 60)
            logger.info(f"New primers selected: {len(result['new_primers'])}")
            logger.info(f"Coverage: {coverage_before:.1%} -> {coverage_after:.1%}")
            logger.info(f"Improvement: {predicted_improvement:.2f}x")
            logger.info(f"Gaps remaining: {len(gaps_after)}")

        return ExpansionResult(
            new_primers=result["new_primers"],
            combined_set=combined_set,
            fixed_primers=fixed_primers,
            coverage_before=coverage_before,
            coverage_after=coverage_after,
            gap_coverage=gap_coverage,
            predicted_improvement=predicted_improvement,
            gaps_remaining=len(gaps_after),
            optimization_method=optimization_method,
            message=result.get("message", ""),
            stage_history=tuple(result.get("stage_history", ())),
        )

    def _recovered_deficit_fraction(self, target_gaps, combined_set) -> float:
        """Fraction of the targeted deficit BASES the delivered panel reaches.

        Weight 1 inside a target gap and 0 outside, which is the in-silico
        case: a gap is a total deficit. A BAM-derived weight per base slots in
        here unchanged once selection reads it.

        In [0, 1] by construction, so it cannot report progress as regression
        the way a count of gaps did.
        """
        length_by_prefix = dict(zip(self.fg_prefixes, self.fg_seq_lengths, strict=True))
        total = 0.0
        recovered = 0.0
        for prefix, length in length_by_prefix.items():
            spans = [
                (g.start, min(g.end, length))
                for g in target_gaps
                if g.chromosome == prefix and g.start < length
            ]
            if not spans:
                continue
            weights = np.zeros(int(length), dtype=np.float64)
            for start, end in spans:
                weights[int(start) : int(end)] = 1.0
            total += float(weights.sum())

            positions = []
            for primer in combined_set:
                positions.extend(int(p) for p in self.cache.get_positions(prefix, primer, "both"))
            if positions:
                recovered += recovered_deficit(
                    positions,
                    weights,
                    extension=self.coverage_reach,
                    length=int(length),
                    circular=bool(getattr(self, "fg_circular", False)),
                    record_starts=self.cache.get_record_starts(prefix) or None,
                )
        if total <= 0:
            return 0.0
        return float(min(1.0, recovered / total))

    def _filter_candidates_to_gaps(
        self,
        candidates: List[str],
        target_gaps: List[CoverageGap],
    ) -> List[str]:
        """Keep candidates whose WINDOW can reach any target gap.

        The intervals are the gaps dilated by one coverage reach, so a
        candidate binding within a reach of a gap is kept: its modelled window
        blankets part of the gap, and the old membership test threw away
        exactly those. Finding F7.

        This is a budget prescreen. It narrows what is scored and decides
        nothing; `deficit_objective` ranks what survives by how much missing
        depth it actually recovers.

        Gaps are grouped by chromosome (== fg prefix). For each candidate, its
        binding positions on that prefix are tested against the dilated
        intervals with a vectorized ``np.searchsorted``. A gap whose ``end``
        exceeds the genome length (a circular wrap gap) is split into
        ``[start, length)`` and ``[0, end-length)`` before dilation.
        """
        # Build per-prefix sorted interval bounds, expanding wrap gaps.
        length_by_prefix = dict(zip(self.fg_prefixes, self.fg_seq_lengths, strict=True))
        bounds_by_prefix: dict = {}
        for g in target_gaps:
            length = length_by_prefix.get(g.chromosome)
            intervals = bounds_by_prefix.setdefault(g.chromosome, [])
            if length is not None and g.end > length:
                intervals.append((g.start, length))
                intervals.append((0, g.end - length))
            else:
                intervals.append((g.start, g.end))

        # Dilate by one reach and merge, then pre-sort for searchsorted.
        prepared = {}
        for prefix, intervals in bounds_by_prefix.items():
            widened = dilate_intervals(
                intervals,
                reach=self.coverage_reach,
                length=length_by_prefix.get(prefix, max(e for _, e in intervals)),
            )
            starts = np.array([s for s, _ in widened], dtype=np.int64)
            ends = np.array([e for _, e in widened], dtype=np.int64)
            prepared[prefix] = (starts, ends)

        kept = []
        for cand in candidates:
            hit = False
            for prefix, (starts, ends) in prepared.items():
                positions = self.cache.get_positions(prefix, cand, "both")
                if len(positions) == 0:
                    continue
                pos = np.asarray(positions, dtype=np.int64)
                # For each position, the candidate interval is the one whose
                # start is the largest <= pos; check pos < that interval's end.
                idx = np.searchsorted(starts, pos, side="right") - 1
                valid = idx >= 0
                if np.any(valid):
                    if np.any(pos[valid] < ends[idx[valid]]):
                        hit = True
                        break
            if hit:
                kept.append(cand)
        return kept

    def _deficit_weights(self, target_gaps):
        """Weight 1 inside a target gap, 0 outside, per prefix.

        The in-silico case, where a gap IS a total deficit. A BAM-derived
        weight per base slots in here unchanged.
        """
        weights_by_prefix = {}
        lengths_by_prefix = dict(zip(self.fg_prefixes, self.fg_seq_lengths, strict=True))
        for prefix, length in lengths_by_prefix.items():
            spans = [
                (g.start, min(g.end, length))
                for g in target_gaps
                if g.chromosome == prefix and g.start < length
            ]
            if not spans:
                continue
            weights = np.zeros(int(length), dtype=np.float64)
            for start, end in spans:
                weights[int(start) : int(end)] = 1.0
            weights_by_prefix[prefix] = weights
        return weights_by_prefix, lengths_by_prefix

    def _attach_deficit_objective(self, optimizer, target_gaps, verbose):
        """Make the search rank by recovered deficit rather than by breadth.

        `attach_search_config` sets it on the wrapper AND the delegate: every
        command-line path is handed a wrapper that delegates the search to an
        inner `HybridOptimizer`, and assigning to the wrapper alone left the
        refinement reading an attribute nobody had set, undetected for months.
        """
        weights_by_prefix, lengths_by_prefix = self._deficit_weights(target_gaps)
        if not weights_by_prefix:
            return
        from neoswga.core.pool_objective import PoolConstraints, PoolObjective
        from neoswga.core.swap_refinement import attach_search_config

        from .panel_refinement import objective_for_optimizer

        inner = objective_for_optimizer(optimizer) or PoolObjective(
            optimizer.compute_metrics, PoolConstraints()
        )
        objective = DeficitObjective(
            inner,
            cache=self.cache,
            weights_by_prefix=weights_by_prefix,
            lengths_by_prefix=lengths_by_prefix,
            extension=self.coverage_reach,
            circular=bool(getattr(self, "fg_circular", False)),
        )
        attach_search_config(optimizer, "pool_objective", objective)
        if verbose:
            total = sum(float(w.sum()) for w in weights_by_prefix.values())
            logger.info(
                "Selection ranks by recovered deficit over %.0f targeted bases, "
                "not by genome breadth.",
                total,
            )

    def _build_hybrid_optimizer(self, background_pruning=None, refinement_method="network"):
        """Hybrid optimizer configured with both reaches, kept distinct.

        `coverage_reach` governs Stage-1 set cover and the reported coverage;
        `max_extension` governs the Stage-2 amplification-network question
        "could two primers connect via one uninterrupted extension?", which is
        single-molecule processivity. This class passed its 3 kb coverage reach
        as `max_extension`, so the network stage ran at a reach 23x too short
        for phi29. Leaving `max_extension` unset takes the polymerase preset.
        """
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        return HybridOptimizer(
            # `network` does not read `pool_objective` at all, so a deficit
            # objective attached to it would be a measurement that reached the
            # code and not the user. `swap` is the Stage 2 that scores
            # `(-shortfall, coverage, -background)`. The caller chooses, because
            # the two Stage 2s carry different things: the host term Known
            # Issue 14 added lives in `network`.
            refinement_method=refinement_method,
            position_cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bg_prefixes=self.bg_prefixes,
            bg_seq_lengths=self.bg_seq_lengths,
            bin_size=self.bin_size,
            coverage_reach=self.coverage_reach,
            # Background pruning is the ONLY stage that reads `bg_prefixes`.
            # It defaults to False, and this builder used to leave it there, so
            # expansion was handed a host genome and never looked at it:
            # measured on the plasmid example, a real run queried the target
            # prefix 136 times and the background prefix zero times. Enabled
            # whenever there is a background to read, because specificity is
            # the property `expand-primers` exists to preserve.
            background_pruning=(
                background_pruning if background_pruning is not None else bool(self.bg_prefixes)
            ),
        )

    def _expand_configured(
        self,
        candidates,
        fixed_primers,
        target_new,
        verbose,
        method="hybrid",
        target_gaps=None,
        background_pruning=None,
    ):
        """Configured expansion uses the same panel stages as a new design."""
        from .optimization_service import OptimizationRequest, panel_violations, run_panel_search
        from .optimizer_factory import OptimizerFactory
        from .panel_refinement import objective_for_optimizer
        from .unified_optimizer import _ensure_optimizers_registered

        _ensure_optimizers_registered()
        optimizer = OptimizerFactory.create(
            name=method,
            position_cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bg_prefixes=self.bg_prefixes,
            bg_seq_lengths=self.bg_seq_lengths,
            conditions=self.conditions,
            config=self.context.optimizer_config(verbose=verbose, refinement_method="swap"),
            polymerase=self.context.polymerase,
            background_pruning=(
                bool(self.bg_prefixes) if background_pruning is None else background_pruning
            ),
        )
        objective_for_optimizer(optimizer, self.context.constraints)
        if target_gaps:
            self._attach_deficit_objective(optimizer, target_gaps, verbose)
        result = run_panel_search(
            OptimizationRequest(
                optimizer,
                tuple(candidates),
                len(fixed_primers) + target_new,
                fixed_primers=tuple(fixed_primers),
                excluded_primers=getattr(self, "_excluded_primers", ()),
                candidate_source=self.candidate_source,
                prepare_candidates=(
                    (lambda pool: self._filter_candidates_to_gaps(pool, target_gaps))
                    if target_gaps
                    else None
                ),
            )
        )
        violations = (
            panel_violations(optimizer, result.primers) if result.primers else ("no panel found",)
        )
        if violations:
            raise ValueError("Expansion did not find a qualifying pool: " + "; ".join(violations))
        return {
            "new_primers": [p for p in result.primers if p not in set(fixed_primers)],
            "message": result.message,
            "stage_history": list(result.stage_history),
        }

    def _expand_hybrid(
        self,
        candidates: List[str],
        fixed_primers: List[str],
        target_new: int,
        verbose: bool,
        background_pruning=None,
        target_gaps=None,
    ) -> Dict:
        """Use hybrid optimizer for expansion.

        The Stage 2 is chosen rather than fixed, because the two carry
        different things and neither carries both. `network` holds the host
        term Known Issue 14 added, which is the whole of what
        `background-aware` buys. `swap` is the one that reads a
        `pool_objective`, which is how a deficit objective can steer anything.

        So a host-aware expansion does not yet target the deficit, and a
        deficit-targeted one is not host-aware. That is stated here and in the
        log rather than left for a reader to discover from a panel that quietly
        ignored one of the two.
        """
        if self.context is not None:
            return self._expand_configured(
                candidates,
                fixed_primers,
                target_new,
                verbose,
                target_gaps=target_gaps,
                background_pruning=background_pruning,
            )
        wants_deficit = bool(target_gaps) and not background_pruning
        optimizer = self._build_hybrid_optimizer(
            background_pruning=background_pruning,
            refinement_method="swap" if wants_deficit else "network",
        )
        if wants_deficit:
            self._attach_deficit_objective(optimizer, target_gaps, verbose)
        elif target_gaps and background_pruning and verbose:
            logger.warning(
                "Host-aware expansion keeps the network refinement, which "
                "carries the host term but cannot read a deficit objective, so "
                "the %d target gap(s) narrow the candidate pool without steering "
                "selection. Run without background pruning to rank by recovered "
                "deficit instead.",
                len(target_gaps),
            )

        # Target total = fixed + new
        total_target = len(fixed_primers) + target_new

        result = optimizer.optimize(
            candidates=candidates,
            final_count=total_target,
            fixed_primers=fixed_primers,
            verbose=verbose,
            # total_target is len(fixed) + requested new, so rescaling it could
            # ask for fewer primers than the caller is already keeping.
            apply_polymerase_multiplier=False,
        )

        # Extract new primers (those not in fixed)
        fixed_set = set(fixed_primers)
        new_primers = [p for p in result.primers if p not in fixed_set]

        return {
            "new_primers": new_primers,
            "message": f"Connectivity: {result.final_connectivity:.2f}",
        }

    def _expand_dominating_set(
        self,
        candidates: List[str],
        fixed_primers: List[str],
        target_new: int,
        verbose: bool,
        target_gaps=None,
    ) -> Dict:
        """Use dominating set optimizer for expansion."""
        if self.context is not None:
            return self._expand_configured(
                candidates,
                fixed_primers,
                target_new,
                verbose,
                method="dominating-set",
                target_gaps=target_gaps,
            )
        from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

        optimizer = DominatingSetOptimizer(
            cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bin_size=self.bin_size,
            extension_reach=self.coverage_reach,
        )

        result = optimizer.optimize_greedy(
            candidates=candidates,
            max_primers=target_new,
            fixed_primers=fixed_primers,
            verbose=verbose,
        )

        return {
            "new_primers": result["new_primers"],
            "message": f"Coverage: {result['coverage']:.1%}",
        }


def expand_primers(
    params_path: str,
    fixed_primers: List[str],
    failed_primers: Optional[List[str]] = None,
    target_new: int = 6,
    optimization_method: str = "hybrid",
    output_dir: Optional[str] = None,
    verbose: bool = True,
    bam_path: Optional[str] = None,
    min_depth: int = 5,
    min_gap_size: int = 10000,
    contig_aliases: Optional[Dict[str, str]] = None,
) -> ExpansionResult:
    """
    Convenience function to expand primer set using params.json.

    Args:
        params_path: Path to params.json
        fixed_primers: Primers to keep (validated)
        failed_primers: Primers to exclude (failed)
        target_new: Number of new primers to add
        optimization_method: Optimization method to use
        output_dir: Directory to save results (optional)
        verbose: Print progress
        bam_path: Optional mapped BAM. When given, low-depth regions are
            computed and merged with the in-silico gaps; the candidate pool is
            then focused on primers that bind inside those gaps.
        min_depth: Depth below which a base counts as a BAM gap (default 5).
        min_gap_size: Minimum gap length to act on (bp, default 10000).
        contig_aliases: Optional {fg_prefix_or_basename: bam_contig} overrides.

    Returns:
        ExpansionResult with new and combined primer sets
    """
    import json
    import os

    import pandas as pd

    from neoswga.core import parameter
    from neoswga.core.position_cache import PositionCache

    # Load parameters
    with open(params_path) as f:
        params = json.load(f)

    # Get prefixes and lengths
    fg_prefixes = params.get("fg_prefixes", [])
    bg_prefixes = params.get("bg_prefixes", [])
    data_dir = params.get("data_dir", "./")

    # Load genome lengths (from step2 or calculate)
    # This is a simplified version - in practice, get from genome files
    fg_seq_lengths = params.get("fg_seq_lengths", [])
    bg_seq_lengths = params.get("bg_seq_lengths", [])

    # Load candidates from step3
    step3_path = os.path.join(data_dir, "step3_df.csv")
    if not os.path.exists(step3_path):
        raise FileNotFoundError(
            f"Step 3 output not found: {step3_path}. Run 'neoswga prepare-candidates' first."
        )

    # Through the shared source. `expand-primers` exists to add primers to an
    # existing panel, so the pool it may draw from is the whole point, and the
    # CSV shortlist is a fraction of what the inventory retains.
    from neoswga.core.candidate_source import open_design_source

    context = design_context_from_params(params)
    conditions = context.conditions
    shortlist = pd.read_csv(step3_path)["primer"].astype(str).tolist()
    candidate_source = open_design_source(
        data_dir,
        conditions.fingerprint() if hasattr(conditions, "fingerprint") else "",
        sorted({len(p) for p in shortlist}),
        fallback=shortlist,
    )
    candidates = candidate_source.initial()

    # Initialize position cache
    # Built over the background prefixes too. `get_positions` answers a
    # prefix the cache was not built over with an empty array, silently, so an
    # fg-only cache made every background lookup below read zero -- which is
    # indistinguishable downstream from a perfectly specific panel.
    cache = PositionCache(fg_prefixes + bg_prefixes, candidates + fixed_primers)

    # The same index check a new design makes. Expansion exists to ADD to a
    # delivered panel, so it is scored against the same references the panel
    # was, and an index that has since gone stale would put the additions on a
    # different genome from the primers they are joining.
    cache.require_record_metadata(
        fg_prefixes + bg_prefixes,
        genomes={
            str(prefix): str(genome)
            for prefixes, genomes in (
                (fg_prefixes, params.get("fg_genomes") or []),
                (bg_prefixes, params.get("bg_genomes") or []),
            )
            # strict=False DELIBERATELY, the same decision as
            # `DesignRequest.reference_manifest` and `pipeline.py`: a prefix
            # with no genome is ABSENT from the manifest rather than paired
            # with a guess, so the identity check does not run for it. Pairing
            # it with whatever was left in a mutable global is what made a
            # design refuse its own index under `pytest -n 8`.
            for prefix, genome in zip(list(prefixes), list(genomes), strict=False)
            if genome
        },
    )

    # Create expander
    # Resolved from the params this function already held. It used to build the
    # expander with none of it, so expansion ran at 3 kb with no chemistry
    # while the panel it was extending had been designed with both.
    context = design_context_from_params(params)
    expander = PrimerExpander(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bg_prefixes=bg_prefixes,
        bg_seq_lengths=bg_seq_lengths,
        context=context,
        candidate_source=candidate_source,
    )

    # Build target gaps: in-silico gaps from the fixed set, optionally merged
    # with low-depth regions from a mapped BAM.
    fg_circular = bool(params.get("fg_circular", False))
    bam_derived_gaps = None
    bam_depth_policy = None
    if bam_path:
        from neoswga.core.bam_coverage import bam_gaps
        from neoswga.core.depth_policy import DepthPolicy

        bam_derived_gaps = bam_gaps(
            bam_path,
            fg_prefixes,
            fg_seq_lengths,
            min_depth=min_depth,
            min_gap_size=min_gap_size,
            circular=fg_circular,
            contig_aliases=contig_aliases,
            # Finding F8: without the FASTA layout a prefix is matched to one
            # BAM contig, so a multi-record reference matches nothing and the
            # gap list is silently empty.
            fg_genomes=params.get("fg_genomes"),
        )
        bam_depth_policy = DepthPolicy().to_dict()
        if verbose:
            logger.info(f"BAM low-depth gaps: {len(bam_derived_gaps)}")
            logger.info(DepthPolicy().describe())

    target_gaps = expander.identify_gaps(
        fixed_primers,
        min_gap_size=min_gap_size,
        extra_gaps=bam_derived_gaps,
        merge=True,
    )

    # Run expansion (focus candidates on the gaps when we have any)
    result = expander.expand(
        candidates=candidates,
        fixed_primers=fixed_primers,
        failed_primers=failed_primers,
        target_new=target_new,
        optimization_method=optimization_method,
        verbose=verbose,
        target_gaps=target_gaps or None,
    )
    if bam_depth_policy is not None:
        # Attached here rather than threaded through `expand`, which does not
        # read a BAM and should not carry a BAM concern. Recorded because a
        # breadth figure means nothing without the rule that produced it.
        result = dataclasses.replace(result, depth_policy=bam_depth_policy)

    # Save results if output_dir specified
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

        # Save result as JSON
        result_path = os.path.join(output_dir, "expansion_result.json")
        with open(result_path, "w") as f:
            json.dump(result.to_dict(), f, indent=2)

        # Save combined primers as CSV
        primers_path = os.path.join(output_dir, "expanded_primers.csv")
        pd.DataFrame(
            {
                "primer": result.combined_set,
                "type": ["fixed"] * result.n_fixed + ["new"] * result.n_new,
            }
        ).to_csv(primers_path, index=False)

        if verbose:
            logger.info(f"\nResults saved to {output_dir}")

    return result
