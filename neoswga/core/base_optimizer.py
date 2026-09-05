"""
Base optimizer interface for SWGA primer set optimization.

All optimizer implementations should inherit from BaseOptimizer to ensure
consistent interfaces, enable polymorphism, and simplify testing.

Design principles:
- Strategy pattern: swap optimizers without changing client code
- Dependency injection: pass caches and conditions, don't create internally
- Immutable results: OptimizationResult is a frozen dataclass
"""

import logging
import math
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Set, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Selectivity ratio that scores a full 1.0 on the normalised scale. Not a
# threshold for a usable design -- it is the point past which more selectivity
# stops earning score, chosen so the term discriminates across the range real
# designs occupy rather than saturating inside it.
SELECTIVITY_REFERENCE = 100.0

# Stands in for unbounded selectivity when no background binding is
# detected at all, so the value stays finite and serialisable. Same
# convention as multi_genome_filter.MAX_ENRICHMENT.
MAX_SELECTIVITY = 1e6


class OptimizationStatus(Enum):
    """Status of optimization run."""

    SUCCESS = "success"
    PARTIAL = "partial"  # Converged but below target
    NO_CONVERGENCE = "no_convergence"
    ERROR = "error"


def json_safe(value: Any) -> Any:
    """
    Replace non-finite floats with None so the result is valid JSON.

    A failed optimization carries `score=-inf` and `mean_gap`/`max_gap` of `inf`
    -- deliberately, since "no coverage" really is unboundedly bad and any finite
    sentinel would rank it as merely mediocre. But `json.dump` writes those as
    the bare tokens `Infinity` / `-Infinity` / `NaN`, which RFC 8259 does not
    allow. Python reads them back happily, so the damage is invisible from
    inside: `step4_improved_df_summary.json` is the file the report and any
    downstream tooling read, and JavaScript, Go and Rust parsers all reject it.

    `null` is the honest encoding -- the quantity is not representable -- and it
    round-trips to None rather than to a number that would be silently wrong.
    """
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _selectivity_from_loads(fg_load: float, bg_load: float) -> float:
    """Selectivity from foreground and background binding loads.

    The integer version guarded with `max(bg, 1)`, which is right for counts --
    you cannot have a fraction of a site. Occupancy-weighted loads are floats,
    and carrying that guard over silently divided by 1.0 whenever the
    background load fell below it: a set with 0.3 effective background sites
    reported a third of its true selectivity, and one with 0.0 reported the
    foreground load as though it were a ratio.

    A genuinely zero background load is unbounded selectivity, not a large
    number. `MAX_SELECTIVITY` stands in so the value stays finite and
    JSON-serialisable, following the same convention and for the same reason as
    `multi_genome_filter.MAX_ENRICHMENT`: it means "no background binding was
    detected", not "measured this well".
    """
    if bg_load <= 0:
        return MAX_SELECTIVITY if fg_load > 0 else 0.0
    return fg_load / bg_load


def _selectivity_density_from_loads(
    fg_load: float, fg_length: float, bg_load: float, bg_length: float
) -> float:
    """Binding load per base of target against per base of background.

    `_selectivity_from_loads` divides two counts and carries no genome length,
    so its value moves with how much background sequence the caller supplied.
    Substituting a whole human genome (3.1 Gb) for the single chromosome often
    used as a stand-in (chr21, 46.7 Mb) multiplies the background load by
    roughly the length ratio and divides the reported selectivity by it, with
    nothing about the primers changed.

    Enrichment does not work that way. Priming is proportional to sites per
    genome copy, and a host contributes its sites spread over its own length,
    so the quantity that survives the substitution -- and the one comparable
    between two designs scored against different backgrounds -- is the ratio of
    densities.

    Both are reported. The count ratio remains what the objective is scored on:
    within a single run the background is fixed, so the two rank sets
    identically and selection is unaffected.
    """
    if fg_length <= 0 or bg_length <= 0:
        # No length is an absent measurement, not a clean background.
        return 0.0

    fg_density = fg_load / fg_length
    bg_density = bg_load / bg_length

    if bg_density <= 0:
        # Same convention as the count ratio: unbounded, not merely large.
        return MAX_SELECTIVITY if fg_density > 0 else 0.0
    return fg_density / bg_density


def _union_coverage(positions, total_length: int, reach: int, circular: bool) -> float:
    """Fraction of a sequence within `reach` of any binding site.

    Module-level and free of instance state because three call sites want it --
    the selection reach, each reporting reach, and the circular-wrap tests that
    duck-type an optimizer with only a `config`. Three implementations of one
    quantity is how this codebase has produced disagreeing coverage numbers
    before; one audit found three different semantics for "coverage" at once.
    """
    import numpy as np

    from .coverage import _mark_window

    if not positions or total_length == 0:
        return 0.0

    # On a circular target a single site whose window spans the whole sequence
    # covers everything, and the marking loop cannot represent that.
    if circular and 2 * reach >= total_length:
        return 1.0

    occupied = np.zeros(total_length, dtype=bool)
    for pos in positions:
        _mark_window(occupied, int(pos), reach, total_length, circular)
    return float(occupied.sum()) / total_length


@dataclass(frozen=True)
class PrimerSetMetrics:
    """
    Metrics for evaluating a primer set.

    Immutable to prevent accidental modification after optimization.
    """

    # Coverage metrics
    fg_coverage: float  # Fraction of foreground genome covered
    bg_coverage: float  # Fraction of background genome covered (lower is better)
    coverage_uniformity: float  # Gini coefficient of gap sizes (lower is better)

    # Binding metrics
    total_fg_sites: int  # Total exact-match binding sites in foreground
    total_bg_sites: int  # Total exact-match binding sites in background
    # Occupancy-weighted sites: each mismatch class weighted by how much of the
    # time a duplex of that stability is actually bound under the configured
    # conditions. Reported beside the raw counts so a reader can see both what
    # was counted and what it was worth.
    selectivity_ratio: float  # effective_fg / effective_bg (higher is better)

    # Thermodynamic metrics
    mean_tm: float  # Mean melting temperature
    tm_range: Tuple[float, float]  # (min_tm, max_tm)
    dimer_risk_score: float  # 0-1, fraction of primer pairs with dimer risk

    # Gap statistics
    mean_gap: float  # Mean gap between binding sites
    max_gap: float  # Maximum gap (coverage hole)
    gap_gini: float  # Gini coefficient of gaps
    gap_entropy: float  # Shannon entropy of gap distribution (bits)

    # Strand alternation metrics
    strand_alternation_score: float  # Fraction of adjacent pairs alternating strands (0-1)
    strand_coverage_ratio: float  # Balance between forward and reverse strand sites (0-1)

    # Per-target coverage (Phase 11D). Maps fg_prefix -> coverage fraction
    # so multi-genome runs can surface "target A 95% / target B 40%"
    # instead of a single aggregate. Empty dict in single-genome mode.
    per_target_coverage: Dict[str, float] = field(default_factory=dict)
    effective_fg_sites: float = 0.0
    effective_bg_sites: float = 0.0
    # "occupancy" when conditions were applied, "exact" when the jellyfish
    # count files needed for mismatch classes were unavailable and the metric
    # fell back to counting perfect matches. Recorded rather than inferred:
    # switching silently between two definitions of one number is exactly the
    # failure this audit has met most often.
    selectivity_mode: str = "exact"

    # Binding load per base of target against per base of background, beside
    # the count ratio in `selectivity_ratio`. The count ratio carries no genome
    # length, so it moves with how much background sequence was supplied: a
    # whole human genome in place of the single chromosome commonly used as a
    # stand-in divides it by roughly 66 with nothing about the primers changed.
    # Enrichment follows sites per genome copy, so the density ratio is the
    # comparable quantity. The lengths it was taken over are reported with it,
    # because a density cannot be read without them.
    selectivity_density: float = 0.0
    fg_total_length: int = 0
    bg_total_length: int = 0

    # Coverage weighted by how much of the time each site is actually bound
    # under the configured conditions, beside the raw `fg_coverage` that counts
    # every exact match as fully occupied.
    #
    # Phase 3 made selectivity occupancy-weighted and left coverage a count, so
    # the two halves of the objective sat on different footings: a set could be
    # rewarded for binding sites it does not occupy at the reaction temperature.
    # The bias is not small and it is not uniform. Measured on Prevotella at
    # equiphi29 42 C, a 10-mer set reports 52% coverage at mean occupancy 0.36
    # while a 12-mer set reports 40% at 0.86 -- so the raw numbers rank them in
    # the opposite order to the effective ones.
    #
    # This is also what connects additives to coverage. DMSO and betaine lower
    # effective Tm, which moves every site down its occupancy curve; against a
    # count-based coverage that cost is invisible, and a condition search would
    # spend coverage it could not see it was spending.
    effective_fg_coverage: float = 0.0

    # Raw coverage at several reaches, {reach_bp: fraction}, so the figure can
    # be compared across tools that do not share a convention. swga 2.0 reports
    # at phi29's ~70 kb processivity; NeoSWGA selects at the ~3 kb realistic
    # per-primer reach. The same set measures 0.418 / 0.836 / ~1.0 at 3 / 10 /
    # 70 kb, so a bare "coverage" number is not comparable to a published one.
    # Reporting the curve costs ~1 ms per reach and removes the ambiguity.
    coverage_by_reach: Dict[int, float] = field(default_factory=dict)

    # Coverage at the MEASURED product reach (~10 kb phi29), which is what a
    # reader should be shown as "the" coverage of a design, and the reach it was
    # computed at.
    #
    # Selection deliberately uses a shorter reach -- the design density that
    # published successful sets have (2-5 kb spacing) -- because selecting at
    # 10 kb would stop the greedy once the genome was covered at that reach and
    # yield designs sparser than any set with measured wet-lab success. But
    # reporting at the selection reach understates what is amplified, since
    # strand displacement means a downstream primer does not truncate extension.
    # Two questions, two numbers, both reported rather than one standing in for
    # the other.
    headline_coverage: float = 0.0
    headline_reach: int = 0

    # Per-primer extension reach (bp) used to compute fg_coverage. Recorded
    # so the report can label "95% coverage at 3 kb reach" rather than
    # leaving the granularity ambiguous. Defaults to the Phase 16 realistic
    # per-primer reach for phi29.
    extension_reach: int = 3000

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "fg_coverage": self.fg_coverage,
            # Same union of binding windows, weighted by how much of the time
            # each site is bound under the configured conditions.
            "effective_fg_coverage": self.effective_fg_coverage,
            # Keyed by reach in bp. Coverage is not comparable across reaches,
            # and the published tools do not share a convention.
            "coverage_by_reach": {str(k): v for k, v in self.coverage_by_reach.items()},
            # Coverage at the measured product reach; see headline_reach.
            "headline_coverage": self.headline_coverage,
            "headline_reach": self.headline_reach,
            "bg_coverage": self.bg_coverage,
            "coverage_uniformity": self.coverage_uniformity,
            "total_fg_sites": self.total_fg_sites,
            "total_bg_sites": self.total_bg_sites,
            "selectivity_ratio": self.selectivity_ratio,
            # Beside the raw counts, so a reader can see both what was counted
            # and what it was worth under the configured conditions.
            "effective_fg_sites": self.effective_fg_sites,
            "effective_bg_sites": self.effective_bg_sites,
            "selectivity_mode": self.selectivity_mode,
            "selectivity_density": self.selectivity_density,
            "fg_total_length": self.fg_total_length,
            "bg_total_length": self.bg_total_length,
            "mean_tm": self.mean_tm,
            "tm_range": list(self.tm_range),
            "dimer_risk_score": self.dimer_risk_score,
            "mean_gap": self.mean_gap,
            "max_gap": self.max_gap,
            "gap_gini": self.gap_gini,
            "gap_entropy": self.gap_entropy,
            "strand_alternation_score": self.strand_alternation_score,
            "strand_coverage_ratio": self.strand_coverage_ratio,
            "per_target_coverage": dict(self.per_target_coverage),
            "extension_reach": self.extension_reach,
        }

    # Application-profile weights for POST-HOC scoring (Phase 14C).
    # Different use cases weigh the component metrics differently:
    # clinical applications favour selectivity and dimer safety over raw
    # coverage, while metagenomics favours coverage above all.
    APPLICATION_WEIGHTS = {
        "balanced": {
            "coverage_w": 0.35,
            "selectivity_w": 0.30,
            "dimer_w": 0.15,
            "evenness_w": 0.10,
            "tm_w": 0.10,
        },
        "discovery": {
            "coverage_w": 0.50,
            "selectivity_w": 0.15,
            "dimer_w": 0.10,
            "evenness_w": 0.15,
            "tm_w": 0.10,
        },
        "clinical": {
            "coverage_w": 0.20,
            "selectivity_w": 0.45,
            "dimer_w": 0.20,
            "evenness_w": 0.10,
            "tm_w": 0.05,
        },
        "enrichment": {
            "coverage_w": 0.35,
            "selectivity_w": 0.30,
            "dimer_w": 0.15,
            "evenness_w": 0.10,
            "tm_w": 0.10,
        },
        "metagenomics": {
            "coverage_w": 0.60,
            "selectivity_w": 0.10,
            "dimer_w": 0.10,
            "evenness_w": 0.15,
            "tm_w": 0.05,
        },
    }

    def normalized_score(
        self,
        coverage_w: float = 0.35,
        selectivity_w: float = 0.30,
        dimer_w: float = 0.15,
        evenness_w: float = 0.10,
        tm_w: float = 0.10,
        application: Optional[str] = None,
    ) -> float:
        """Compute a [0,1] composite score from metrics.

        Provides a single comparable value regardless of which optimizer
        produced the result.  All component metrics are clamped to [0,1]
        before weighting.

        Args:
            coverage_w: Weight for foreground coverage (default 0.35).
            selectivity_w: Weight for selectivity (default 0.30).
            dimer_w: Weight for dimer safety (default 0.15).
            evenness_w: Weight for gap evenness (default 0.10).
            tm_w: Weight for Tm tightness (default 0.10).
            application: If provided (one of 'balanced', 'discovery',
                'clinical', 'enrichment', 'metagenomics'), overrides the
                individual weight arguments with the preset for that use
                case. See APPLICATION_WEIGHTS.
        """
        if application is not None:
            weights = self.APPLICATION_WEIGHTS.get(application.lower())
            if weights is not None:
                coverage_w = weights["coverage_w"]
                selectivity_w = weights["selectivity_w"]
                dimer_w = weights["dimer_w"]
                evenness_w = weights["evenness_w"]
                tm_w = weights["tm_w"]
        # Coverage: already [0,1].
        #
        # Prefer the occupancy-weighted figure when conditions were available to
        # compute it. Selecting on the raw count rewards a set for reaching
        # sites it does not occupy at the reaction temperature, and the bias is
        # not uniform across candidates -- on Prevotella at equiphi29 42 C a
        # 10-mer set reports 0.520 raw against 0.245 effective while a 12-mer
        # set reports 0.403 against 0.387, so the two figures rank them in
        # opposite orders. Optimizing the raw one also hides the coverage that
        # additives cost, which is the lever this tool exists to use.
        #
        # Falls back to the raw figure when no conditions are attached, matching
        # how `selectivity_mode` degrades: 0.0 means "not computed" here, not
        # "nothing is covered".
        cov = min(self.fg_coverage, 1.0)
        if self.effective_fg_coverage > 0.0:
            cov = min(self.effective_fg_coverage, 1.0)

        # Selectivity, on a log scale.
        #
        # This was `min(ratio / 100, 1.0)`, calibrated when selectivity was a
        # count of exact matches -- ratios routinely exceeded 100 and saturated
        # the cap, so the term behaved as a near-binary "good enough".
        # Occupancy weighting counts near-match background load too, so ratios
        # now land around 5-15 on the bundled example, where the same linear
        # form contributes 0.05-0.15 of its weight and the term is effectively
        # ignored instead. Saturating and vanishing are the same failure at
        # opposite ends.
        #
        # A ratio spanning orders of magnitude wants a log scale: a design at
        # SELECTIVITY_REFERENCE scores 1.0, and the term stays responsive
        # across everything below it rather than being decided in the first or
        # last few percent of the range.
        sel = 0.0
        if self.selectivity_ratio > 0:
            sel = math.log10(1.0 + self.selectivity_ratio) / math.log10(1.0 + SELECTIVITY_REFERENCE)
            sel = max(0.0, min(sel, 1.0))

        # Dimer safety: 1 - risk (lower risk is better)
        dimer_safe = 1.0 - min(self.dimer_risk_score, 1.0)

        # Evenness: 1 - Gini (lower Gini is more uniform)
        even = 1.0 - min(self.gap_gini, 1.0)

        # Tm tightness: narrow Tm spread is better
        tm_spread = abs(self.tm_range[1] - self.tm_range[0]) if self.tm_range else 0.0
        tm_tight = max(0.0, 1.0 - tm_spread / 20.0)  # 0C spread -> 1.0, 20C+ -> 0.0

        return (
            coverage_w * cov
            + selectivity_w * sel
            + dimer_w * dimer_safe
            + evenness_w * even
            + tm_w * tm_tight
        )

    @classmethod
    def empty(cls) -> "PrimerSetMetrics":
        """Create empty metrics for failed optimization."""
        return cls(
            fg_coverage=0.0,
            bg_coverage=0.0,
            coverage_uniformity=1.0,
            total_fg_sites=0,
            total_bg_sites=0,
            selectivity_ratio=0.0,
            mean_tm=0.0,
            tm_range=(0.0, 0.0),
            dimer_risk_score=1.0,
            mean_gap=float("inf"),
            max_gap=float("inf"),
            gap_gini=1.0,
            gap_entropy=0.0,
            strand_alternation_score=0.0,
            strand_coverage_ratio=0.0,
        )


# Phase 15B: application weights that affect SELECTION. Mapped onto
# NetworkOptimizer's internal knobs (tm_weight, uniformity_weight,
# dimer_penalty) so clinical users get selectivity-first primer sets,
# not just a different post-hoc ranking of the same set.
#
# The triple shapes the greedy graph selection:
#   - tm_weight: preference for primers whose effective Tm sits in the
#     reaction-appropriate window.
#   - uniformity_weight: preference for sets with even coverage.
#   - dimer_penalty: how aggressively to reject dimer-prone primer pairs.
#
# Re-projected from APPLICATION_WEIGHTS so the post-hoc normalized_score
# and the selection-time weighting push in the same direction per use case.
OPTIMIZER_APPLICATION_WEIGHTS = {
    "balanced": {"tm_weight": 0.25, "uniformity_weight": 0.10, "dimer_penalty": 0.30},
    "discovery": {"tm_weight": 0.15, "uniformity_weight": 0.20, "dimer_penalty": 0.20},
    "clinical": {"tm_weight": 0.30, "uniformity_weight": 0.05, "dimer_penalty": 0.45},
    "enrichment": {"tm_weight": 0.25, "uniformity_weight": 0.10, "dimer_penalty": 0.30},
    "metagenomics": {"tm_weight": 0.10, "uniformity_weight": 0.25, "dimer_penalty": 0.20},
}


@dataclass(frozen=True)
class OptimizationResult:
    """
    Result of primer set optimization.

    Immutable to ensure results are not accidentally modified.
    Contains the selected primers, score, and detailed metrics.
    """

    primers: Tuple[str, ...]  # Selected primer sequences (tuple for immutability)
    score: float  # Overall optimization score
    status: OptimizationStatus
    metrics: PrimerSetMetrics
    iterations: int  # Number of iterations/generations used
    optimizer_name: str  # Name of optimizer that produced this result

    # Optional detailed results
    all_scores: Optional[Tuple[float, ...]] = None  # Score history
    message: str = ""  # Status message or error description

    # Optional Pareto front for multi-objective optimizers (Phase 14A).
    # Each entry is (tuple-of-primers, metrics). MOEA populates this; other
    # optimizers leave it None. `pareto_metrics` carries the per-solution
    # metrics when the Pareto front is present.
    pareto_front: Optional[Tuple[Tuple[str, ...], ...]] = None
    pareto_metrics: Optional[Tuple[Dict[str, Any], ...]] = None

    # Per-method comparison rows when this result was produced by an ensemble
    # run (method='ensemble'). Each row is a dict:
    #   {method, normalized_score, score, n_primers, fg_coverage,
    #    bg_coverage, status, selected}
    # None for single-method runs.
    ensemble_comparison: Optional[Tuple[Dict[str, Any], ...]] = None

    @property
    def num_primers(self) -> int:
        """Number of primers in the set."""
        return len(self.primers)

    @property
    def is_success(self) -> bool:
        """Whether optimization succeeded."""
        return self.status == OptimizationStatus.SUCCESS

    @property
    def normalized_score(self) -> float:
        """Comparable [0,1] score derived from metrics.

        Unlike ``score`` (which varies in scale per optimizer), this value
        is always in [0,1] and comparable across different optimizer types.
        Returns 0.0 only for ERROR or NO_CONVERGENCE status; PARTIAL
        results with valid metrics are scored normally.
        """
        if self.status in (OptimizationStatus.ERROR, OptimizationStatus.NO_CONVERGENCE):
            return 0.0
        return self.metrics.normalized_score()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        d = {
            "primers": list(self.primers),
            "score": json_safe(self.score),
            "status": self.status.value,
            "metrics": {k: json_safe(v) for k, v in self.metrics.to_dict().items()},
            "iterations": self.iterations,
            "optimizer_name": self.optimizer_name,
            "num_primers": self.num_primers,
            "message": self.message,
        }
        if self.pareto_front is not None:
            d["pareto_front"] = [list(p) for p in self.pareto_front]
            if self.pareto_metrics is not None:
                d["pareto_metrics"] = list(self.pareto_metrics)
        if self.ensemble_comparison is not None:
            d["ensemble_comparison"] = list(self.ensemble_comparison)
        return d

    def validate(
        self,
        target_size: Optional[int] = None,
        min_coverage: float = 0.0,
        min_per_target_coverage: float = 0.0,
        forbidden_primers: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """Post-optimization sanity validation.

        Catches the broad class of silent failures an optimizer can produce
        without explicit assertions: duplicates, unexpected set size, zero
        coverage, blacklist re-injection. Returns a report dict that
        `unified_optimizer.run_optimization` stores on the result and writes
        alongside ``step4_improved_df.csv``.

        This never mutates the result. It emits warnings (``level='warning'``)
        for degenerate-but-acceptable outcomes (e.g., target_size underfilled
        when status is already PARTIAL) and errors (``level='error'``) for
        genuine correctness bugs (e.g., duplicate primers).
        """
        issues: List[Dict[str, Any]] = []

        primer_set = list(self.primers)
        n = len(primer_set)

        # Duplicate check — a hard error regardless of status.
        unique = set(primer_set)
        if len(unique) != n:
            dup_counts = {p: primer_set.count(p) for p in unique}
            duplicates = [p for p, c in dup_counts.items() if c > 1]
            issues.append(
                {
                    "level": "error",
                    "code": "duplicate_primers",
                    "detail": f"{len(duplicates)} duplicate primer(s): {duplicates[:5]}",
                }
            )

        # Set size check — warning unless the optimizer already flagged PARTIAL.
        if target_size is not None and n != target_size:
            issues.append(
                {
                    "level": ("warning" if self.status.value == "partial" else "error"),
                    "code": "set_size_mismatch",
                    "detail": f"Requested target_size={target_size}, got {n}",
                }
            )

        # Non-empty foreground coverage. ERROR status already carries 0.0
        # coverage implicitly, so skip that case.
        if self.status.value != "error":
            fg_cov = getattr(self.metrics, "fg_coverage", 0.0)
            if fg_cov < min_coverage:
                issues.append(
                    {
                        "level": ("warning" if self.status.value == "partial" else "error"),
                        "code": "coverage_below_threshold",
                        "detail": (f"fg_coverage={fg_cov:.3f} < min_coverage={min_coverage:.3f}"),
                    }
                )

        # Per-target coverage (multi-genome mode). Optimizers that populate
        # `metrics.per_target_coverage` (Phase 11D) have their minimum
        # checked here. Missing attribute means single-genome mode or the
        # optimizer did not report per-target numbers yet.
        per_target = getattr(self.metrics, "per_target_coverage", None)
        if per_target and min_per_target_coverage > 0.0:
            below = {k: v for k, v in per_target.items() if v < min_per_target_coverage}
            if below:
                issues.append(
                    {
                        "level": "warning",
                        "code": "per_target_coverage_below_threshold",
                        "detail": (
                            f"{len(below)} target(s) below "
                            f"{min_per_target_coverage:.2f}: {below}"
                        ),
                    }
                )

        # Blacklist re-injection guard. Caller passes forbidden_primers when
        # a blacklist is configured; catches the case where expand-primers
        # or swap-primer fed a non-filtered candidate pool into the optimizer.
        if forbidden_primers:
            hits = set(primer_set) & set(forbidden_primers)
            if hits:
                issues.append(
                    {
                        "level": "error",
                        "code": "blacklist_primer_in_set",
                        "detail": f"Primers in blacklist: {sorted(hits)[:5]}",
                    }
                )

        has_error = any(i["level"] == "error" for i in issues)
        return {
            "optimizer": self.optimizer_name,
            "num_primers": n,
            "ok": not has_error,
            "issues": issues,
        }

    @classmethod
    def failure(cls, optimizer_name: str, message: str) -> "OptimizationResult":
        """Create a failure result."""
        return cls(
            primers=tuple(),
            score=float("-inf"),
            status=OptimizationStatus.ERROR,
            metrics=PrimerSetMetrics.empty(),
            iterations=0,
            optimizer_name=optimizer_name,
            message=message,
        )


@dataclass
class OptimizerConfig:
    """
    Base configuration for optimizers.

    Subclasses can extend with optimizer-specific parameters.
    """

    target_set_size: int = 6
    max_iterations: int = 100
    max_dimer_bp: int = 4
    # Self-dimer threshold. The clique optimizer reached for this with
    # `getattr(self.config, "max_self_dimer_bp", max_dimer_bp + 1)` and the
    # fallback fired every time, because the field did not exist -- so a
    # params.json setting it had no effect anywhere.
    max_self_dimer_bp: int = 5
    # Highest mismatch class counted when weighting background load. 0
    # reproduces exact-match counting. 2 is affordable (405 lookups for a
    # 10-mer) but rarely changes the ranking, since a two-mismatch duplex is
    # mostly unbound at any usable stringency.
    max_mismatches: int = 1
    min_tm: float = 20.0
    max_tm: float = 50.0
    verbose: bool = True

    # Polymerase-specific extension reach (bp) for coverage computation.
    # Default 3000 is the realistic per-primer reach in a dense SWGA design
    # (phi29 ~3 kb, equiphi29 ~4 kb; see coverage.polymerase_extension_reach
    # and reaction_conditions.get_typical_amplicon_length). Per Phase 16
    # critical gap #2, this value controls the coverage objective that
    # optimizers maximize; using processivity here inflates fg_coverage
    # 5-20x over what the set can actually amplify.
    extension_reach: int = 3000

    # Whether the target sequence is closed-circular. Controls whether the
    # coverage window around a binding site wraps past the sequence ends.
    # Populated from params.json `fg_circular` via unified_optimizer.
    fg_circular: bool = False

    # Convergence criteria
    convergence_threshold: float = 0.001
    patience: int = 10  # Iterations without improvement before stopping

    def validate(self) -> None:
        """Validate configuration parameters."""
        if self.target_set_size < 1:
            raise ValueError(f"target_set_size must be >= 1, got {self.target_set_size}")
        if self.max_iterations < 1:
            raise ValueError(f"max_iterations must be >= 1, got {self.max_iterations}")
        if self.min_tm >= self.max_tm:
            raise ValueError(f"min_tm ({self.min_tm}) must be < max_tm ({self.max_tm})")


class BaseOptimizer(ABC):
    """
    Abstract base class for primer set optimizers.

    All optimizer implementations must inherit from this class and implement
    the optimize() method. This ensures consistent interfaces across different
    optimization strategies (greedy, genetic, network-based, etc.).

    Usage:
        optimizer = SomeOptimizer(cache, fg_prefixes, fg_lengths, ...)
        result = optimizer.optimize(candidates, target_size=10)

        if result.is_success:
            print(f"Selected primers: {result.primers}")
            print(f"Coverage: {result.metrics.fg_coverage:.1%}")

    Subclassing:
        class MyOptimizer(BaseOptimizer):
            @property
            def name(self) -> str:
                return "my-optimizer"

            def optimize(self, candidates, target_size=10, **kwargs):
                # Implementation
                return OptimizationResult(...)
    """

    # Structural flag: subclasses whose selection path actually uses
    # `self.conditions` (Tm correction, mechanistic weighting, polymerase-
    # specific penalties) set this to True. Consumed by `neoswga doctor` to
    # report an honest capability matrix without resorting to source-grep.
    ADDITIVE_AWARE: bool = False

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        conditions=None,
        background_profile=None,
        background_profile_lengths=None,
        background_aggregate: str = "worst-case",
        **_unused_kwargs,
    ):
        """
        Initialize optimizer with position data.

        Args:
            position_cache: PositionCache instance for primer position lookups
            fg_prefixes: Foreground (target) genome HDF5 prefixes
            fg_seq_lengths: Foreground genome lengths
            bg_prefixes: Background (off-target) genome HDF5 prefixes (optional)
            bg_seq_lengths: Background genome lengths (optional)
            config: Optimizer configuration
            conditions: Optional ReactionConditions. When provided, optimizers
                that do additive-aware Tm scoring (network, integrated quality
                scorer, genetic, greedy) pick up DMSO / betaine / etc. from
                this object. Optimizers that do not score on Tm ignore it,
                but we accept it here so every subclass has a consistent API
                and the unified_optimizer can forward it blindly.
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.config = config or OptimizerConfig()
        # Attach conditions to every optimizer for uniformity. Consumers that
        # use additive-aware Tm (network / integrated_quality_scorer) can read
        # self.conditions; others simply leave it as attached metadata.
        self.conditions = conditions

        # Compositional background, used only when there are no background
        # prefixes. A `HostProfile` (or a panel of them) supplies expected
        # k-mer counts so metagenomic and unknown-host designs still have a
        # selectivity term; see `_modelled_site_load`. Held as plain attributes
        # rather than threaded through every subclass constructor, since only
        # the metrics path reads them.
        self.background_profile = background_profile
        self.background_profile_lengths = background_profile_lengths
        self.background_aggregate = background_aggregate

        # Computed properties
        self.fg_total_length = sum(fg_seq_lengths)
        self.bg_total_length = sum(bg_seq_lengths) if bg_seq_lengths else 0

        # Validate
        self.config.validate()
        self._validate_inputs()

    def _validate_inputs(self) -> None:
        """Validate constructor inputs."""
        if len(self.fg_prefixes) != len(self.fg_seq_lengths):
            raise ValueError(
                f"fg_prefixes ({len(self.fg_prefixes)}) and fg_seq_lengths "
                f"({len(self.fg_seq_lengths)}) must have same length"
            )
        if self.bg_prefixes and len(self.bg_prefixes) != len(self.bg_seq_lengths):
            raise ValueError(
                f"bg_prefixes ({len(self.bg_prefixes)}) and bg_seq_lengths "
                f"({len(self.bg_seq_lengths)}) must have same length"
            )

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Human-readable optimizer name.

        Used for logging, result attribution, and factory registration.
        Should be kebab-case (e.g., 'greedy-bfs', 'dominating-set').
        """
        pass

    @property
    def description(self) -> str:
        """
        One-line description of the optimizer.

        Override in subclasses to provide helpful description.
        """
        return f"{self.name} optimizer"

    @property
    def supports_background(self) -> bool:
        """
        Whether this optimizer uses background genome data.

        Some optimizers (e.g., dominating-set) only consider foreground
        coverage and don't use background data for optimization.
        """
        return True

    @abstractmethod
    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs,
    ) -> OptimizationResult:
        """
        Find optimal primer set from candidates.

        This is the main entry point for optimization. Implementations should:
        1. Validate inputs
        2. Run optimization algorithm
        3. Compute metrics for result
        4. Return OptimizationResult

        Args:
            candidates: Pool of candidate primer sequences
            target_size: Desired number of primers (overrides config if provided)
            fixed_primers: Optional list of primers that must be included in the
                result. These primers are pre-selected and the optimizer will
                select additional primers from candidates to complement them.
                Useful for iterative design where some primers have already been
                validated experimentally.
            **kwargs: Optimizer-specific parameters

        Returns:
            OptimizationResult with selected primers and metrics

        Raises:
            ValueError: If candidates is empty or invalid
        """
        pass

    def _validate_candidates(self, candidates: List[str]) -> List[str]:
        """
        Validate and deduplicate candidate primers.

        Args:
            candidates: Raw candidate list

        Returns:
            Validated, deduplicated candidate list

        Raises:
            ValueError: If no valid candidates
        """
        if not candidates:
            raise ValueError("candidates list cannot be empty")

        # Deduplicate while preserving order.
        #
        # Normalised to uppercase, not merely checked in uppercase. This used to
        # validate `primer.upper()` and then append `primer` unchanged, so a
        # lowercase sequence entered the pool as-is -- and `dimer.py` compares
        # translated characters directly, so 'a' never matches 'A'. A perfect
        # 10 bp duplex between one lowercase primer and one uppercase primer
        # measured a complementary run of ZERO, in the clique optimizer's
        # compatibility graph and in the network penalty alike. That is the one
        # input that breaks the dimer-free guarantee outright.
        #
        # Normalising here rather than inside the dimer check keeps every other
        # consumer of the pool -- Tm, GC, k-mer lookup, HDF5 keys -- from having
        # to make the same allowance separately.
        seen = set()
        valid = []
        for primer in candidates:
            normalised = primer.upper() if isinstance(primer, str) else primer
            if normalised not in seen:
                seen.add(normalised)
                # Basic validation
                if normalised and all(c in "ATCG" for c in normalised):
                    valid.append(normalised)
                else:
                    logger.warning(f"Skipping invalid primer: {primer}")

        if not valid:
            raise ValueError("No valid primers in candidates list")

        return valid

    def get_primer_positions(self, primer: str, prefix: str, strand: str = "both") -> np.ndarray:
        """
        Get binding positions for a primer in a genome.

        Wrapper around position cache with consistent interface.

        Args:
            primer: Primer sequence
            prefix: Genome HDF5 prefix
            strand: 'forward', 'reverse', or 'both'

        Returns:
            Array of binding positions (empty array if cache unavailable)
        """
        if self.cache is None:
            logger.warning("Position cache not available")
            return np.array([], dtype=np.int64)

        try:
            return self.cache.get_positions(prefix, primer, strand)
        except Exception as e:
            logger.warning(f"Failed to get positions for {primer}: {e}")
            return np.array([], dtype=np.int64)

    def _effective_site_load(self, primers) -> Tuple[float, float, str]:
        """Occupancy-weighted binding load for the foreground and background.

        Each mismatch class is weighted by how much of the time a duplex of
        that stability is actually bound under the configured conditions:

            load = SUM_j n_j * theta(dH, Tm - j*penalty, T)

        This is what lets an additive change specificity. Selectivity used to be
        a ratio of exact-match counts, which contains no temperature and no free
        energy, so no chemistry could move it for a fixed primer set. Raising
        stringency drops occupancy faster for the mismatched background sites
        than for perfect foreground ones, and the ratio rises.

        Returns `(fg_load, bg_load, mode)`. `mode` is "exact" when the jellyfish
        count files needed for mismatch classes are unavailable -- an oligo set
        scanned straight from a FASTA, for instance -- in which case the caller
        falls back to exact counting. That is reported rather than inferred,
        because silently swapping between two definitions of one number is the
        failure this codebase has produced most often.
        """
        from neoswga.core.occupancy import weighted_site_load

        conditions = getattr(self, "conditions", None)
        if conditions is None:
            return 0.0, 0.0, "exact"

        max_mismatches = int(getattr(self.config, "max_mismatches", 1) or 0)

        if not self.bg_prefixes:
            return self._modelled_site_load(primers, conditions, max_mismatches)

        try:
            fg_load = weighted_site_load(primers, self.fg_prefixes, conditions, max_mismatches)
            bg_load = weighted_site_load(primers, self.bg_prefixes, conditions, max_mismatches)
        except (FileNotFoundError, OSError):
            return 0.0, 0.0, "exact"

        return fg_load, bg_load, "occupancy"

    def _modelled_site_load(self, primers, conditions, max_mismatches):
        """Background load from a compositional profile, when there is no genome.

        Metagenomic SWGA has no host to design against, and dropping the
        background makes every design report `MAX_SELECTIVITY` -- specificity
        switched off rather than absent. A Markov profile of the expected
        background supplies the missing term: expected k-mer counts follow from
        a transition table rather than from a genome, and feed the same
        occupancy weighting, so a modelled background is the same quantity as a
        real one measured a different way.

        Predicting real chr21 occupancy load for random 12-mers, the profile
        reaches rho 0.888 fitted to the actual host and 0.701 fitted to an
        unrelated organism -- so even a generic profile ranks background load
        far better than the nothing it replaces.

        Reports mode "modelled" rather than "occupancy", because a predicted
        background must not be mistaken for a measured one.
        """
        profile = getattr(self, "background_profile", None)
        if profile is None:
            return 0.0, 0.0, "exact"

        from neoswga.core.host_profile import aggregate_loads, expected_site_load
        from neoswga.core.occupancy import weighted_site_load

        try:
            fg_load = weighted_site_load(primers, self.fg_prefixes, conditions, max_mismatches)
        except (FileNotFoundError, OSError):
            return 0.0, 0.0, "exact"

        # A panel of candidate hosts is aggregated worst-case by default: one
        # unknown host is drawn from it, so a set that is clean against one and
        # filthy against another is not generic, and a mean would hide that.
        profiles = profile if isinstance(profile, (list, tuple)) else [profile]
        lengths = getattr(self, "background_profile_lengths", None) or [
            getattr(p, "fitted_length", 0) or 3_000_000_000 for p in profiles
        ]
        mode = getattr(self, "background_aggregate", "worst-case")
        loads = [
            expected_site_load(primers, prof, length, conditions, max_mismatches)
            for prof, length in zip(profiles, lengths)
        ]
        return fg_load, aggregate_loads(loads, mode), "modelled"

    def compute_metrics(
        self,
        primers: List[str],
    ) -> PrimerSetMetrics:
        """
        Compute detailed metrics for a primer set.

        Args:
            primers: List of primer sequences

        Returns:
            PrimerSetMetrics with coverage, binding, and thermodynamic stats
        """
        if not primers:
            return PrimerSetMetrics.empty()

        # Collect all positions. Kept per primer as well as pooled: the pooled
        # list cannot say which primer put a site there, and occupancy is a
        # per-primer property.
        # Positions are kept PER PREFIX, because each prefix is a separate
        # sequence with its own 0-based coordinate space. Pooling them into one
        # flat list made position 5000 in target A and position 5000 in target B
        # the same integer, so `set()` collapsed them, no position could exceed
        # the longest single genome, and the denominator was the sum of them
        # all. Two fully covered 100 kb targets reported fg_coverage 0.5145 and
        # a max_gap of 102,000 bp on a 100,000 bp target.
        #
        # It is not only a reporting fault: `--minimize-primers` gates on
        # `fg_coverage`, so a ceiling near 1/N put the 0.70 default out of reach
        # and the trim quietly stopped working on multi-target runs.
        fg_by_prefix = {prefix: set() for prefix in self.fg_prefixes}
        bg_by_prefix = {prefix: set() for prefix in self.bg_prefixes}
        fg_positions_by_primer = {}

        for primer in primers:
            for prefix in self.fg_prefixes:
                positions = self.get_primer_positions(primer, prefix, "both")
                fg_by_prefix[prefix].update(positions.tolist())
                fg_positions_by_primer.setdefault(primer, []).extend(positions.tolist())

            for prefix in self.bg_prefixes:
                positions = self.get_primer_positions(primer, prefix, "both")
                bg_by_prefix[prefix].update(positions.tolist())

        fg_by_prefix = {p: sorted(v) for p, v in fg_by_prefix.items()}
        bg_by_prefix = {p: sorted(v) for p, v in bg_by_prefix.items()}

        # Flat lists remain for the site COUNTS, where summing across sequences
        # is exactly right and de-duplication across them would be wrong.
        fg_positions = [pos for prefix in self.fg_prefixes for pos in fg_by_prefix[prefix]]
        bg_positions = [pos for prefix in self.bg_prefixes for pos in bg_by_prefix[prefix]]

        # Coverage
        fg_coverage = self._coverage_over_prefixes(
            fg_by_prefix, self.fg_prefixes, self.fg_seq_lengths
        )
        effective_fg_coverage = self._compute_effective_coverage(
            fg_positions_by_primer, self.fg_total_length
        )
        coverage_by_reach = self._coverage_by_reach_over_prefixes(
            fg_by_prefix, self.fg_prefixes, self.fg_seq_lengths
        )
        headline_reach = self._product_reach()
        headline_coverage = coverage_by_reach.get(headline_reach) or self._coverage_over_prefixes(
            fg_by_prefix, self.fg_prefixes, self.fg_seq_lengths, reach=headline_reach
        )
        bg_coverage = (
            self._coverage_over_prefixes(bg_by_prefix, self.bg_prefixes, self.bg_seq_lengths)
            if self.bg_total_length > 0
            else 0.0
        )

        # Gap statistics, per sequence then pooled. A "gap" spanning the end of
        # one molecule and the start of another is not a distance a polymerase
        # could ever travel.
        gaps = []
        for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
            gaps.extend(self._compute_gaps(fg_by_prefix.get(prefix, []), length))
        mean_gap = np.mean(gaps) if gaps else float("inf")
        max_gap = max(gaps) if gaps else float("inf")
        gap_gini = self._gini(gaps) if gaps else 1.0

        # Shannon entropy of gap distribution
        if len(gaps) >= 2:
            from scipy.stats import entropy as scipy_entropy

            counts, _ = np.histogram(gaps, bins=min(50, len(gaps)))
            counts = counts[counts > 0]
            probs = counts / counts.sum()
            gap_entropy = float(scipy_entropy(probs, base=2))
        else:
            gap_entropy = 0.0

        # Selectivity
        total_fg = len(fg_positions)
        total_bg = len(bg_positions)
        effective_fg, effective_bg, selectivity_mode = self._effective_site_load(primers)
        if selectivity_mode == "exact":
            effective_fg, effective_bg = float(total_fg), float(total_bg)
        selectivity = _selectivity_from_loads(effective_fg, effective_bg)
        fg_total_length = int(sum(self.fg_seq_lengths or []))
        bg_total_length = int(sum(self.bg_seq_lengths or []))
        selectivity_density = _selectivity_density_from_loads(
            effective_fg, fg_total_length, effective_bg, bg_total_length
        )

        # Tm calculation (simplified)
        tms = [self._estimate_tm(p) for p in primers]
        mean_tm = np.mean(tms)
        tm_range = (min(tms), max(tms))

        # Dimer risk (simplified - would use dimer module in full implementation)
        dimer_risk = self._estimate_dimer_risk(primers)

        # Strand alternation metrics
        strand_alt_score = 0.0
        strand_cov_ratio = 0.0
        if self.cache is not None and hasattr(self.cache, "compute_strand_alternation_stats"):
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                try:
                    strand_stats = self.cache.compute_strand_alternation_stats(
                        prefix, primers, length
                    )
                    strand_alt_score = strand_stats["strand_alternation_score"]
                    strand_cov_ratio = strand_stats["strand_coverage_ratio"]
                    break  # Use first fg genome stats
                except Exception:
                    pass

        return PrimerSetMetrics(
            fg_coverage=fg_coverage,
            effective_fg_coverage=effective_fg_coverage,
            coverage_by_reach=coverage_by_reach,
            headline_coverage=headline_coverage,
            headline_reach=headline_reach,
            bg_coverage=bg_coverage,
            coverage_uniformity=gap_gini,
            total_fg_sites=total_fg,
            total_bg_sites=total_bg,
            selectivity_ratio=selectivity,
            effective_fg_sites=effective_fg,
            effective_bg_sites=effective_bg,
            selectivity_mode=selectivity_mode,
            selectivity_density=selectivity_density,
            fg_total_length=fg_total_length,
            bg_total_length=bg_total_length,
            mean_tm=mean_tm,
            tm_range=tm_range,
            dimer_risk_score=dimer_risk,
            mean_gap=mean_gap,
            max_gap=max_gap,
            gap_gini=gap_gini,
            gap_entropy=gap_entropy,
            strand_alternation_score=strand_alt_score,
            strand_coverage_ratio=strand_cov_ratio,
            extension_reach=self.config.extension_reach,
        )

    def _coverage_over_prefixes(self, positions_by_prefix, prefixes, seq_lengths, reach=None):
        """Coverage across several sequences, weighted by their lengths.

        Each sequence is measured in its own coordinate space and the covered
        BASES are summed, which is the only aggregate that means anything when
        the parts are separate molecules. Equivalent to `_compute_coverage` when
        there is one prefix, so single-target behaviour is unchanged.
        """
        total = sum(int(length) for length in seq_lengths)
        if total <= 0:
            return 0.0

        covered = 0.0
        for prefix, length in zip(prefixes, seq_lengths):
            length = int(length)
            if length <= 0:
                continue
            positions = positions_by_prefix.get(prefix, [])
            fraction = (
                self._compute_coverage(positions, length)
                if reach is None
                else self._compute_coverage_at(positions, length, reach)
            )
            covered += fraction * length
        return covered / total

    def _coverage_by_reach_over_prefixes(self, positions_by_prefix, prefixes, seq_lengths):
        """`_compute_coverage_by_reach`, aggregated the same way."""
        reaches = set()
        for prefix, length in zip(prefixes, seq_lengths):
            reaches.update(
                self._compute_coverage_by_reach(positions_by_prefix.get(prefix, []), int(length))
            )
        return {
            reach: self._coverage_over_prefixes(
                positions_by_prefix, prefixes, seq_lengths, reach=reach
            )
            for reach in sorted(reaches)
        }

    def _compute_coverage(self, positions: List[int], total_length: int) -> float:
        """Compute coverage fraction from positions.

        Uses ``config.extension_reach`` (per-primer reach in bp) to
        determine how far each binding site extends along the genome.
        Each binding site marks ``[pos - reach, pos + reach]`` occupied
        and coverage is the union across all positions.

        When ``config.fg_circular`` is True, the window wraps past the
        sequence ends instead of being clipped. Delegates the per-site
        marking to :func:`coverage._mark_window` so this implementation
        stays in sync with :func:`coverage.compute_per_prefix_coverage`.
        """
        return _union_coverage(
            positions,
            total_length,
            self.config.extension_reach,
            getattr(self.config, "fg_circular", False),
        )

    def _compute_coverage_by_reach(self, positions, total_length: int) -> Dict[int, float]:
        """Raw coverage at each reporting reach, plus the one used for selection.

        Coverage is meaningless without the reach it was computed at, and the
        published tools do not share a convention: swga 2.0 reports at phi29's
        ~70 kb processivity, NeoSWGA selects at the ~3 kb realistic per-primer
        reach. The same set measures 0.418 at 3 kb and ~1.0 at 70 kb, so a bare
        coverage figure cannot be compared to a published one in either
        direction. Reporting the curve settles that without either convention
        having to move.

        The selection reach is always included even when it is not one of the
        standard reporting points, so `fg_coverage` can always be located on
        this curve rather than appearing to come from nowhere.
        """
        from .coverage import REPORTING_REACHES

        if not positions or total_length == 0:
            return {}

        reaches = sorted(
            set(REPORTING_REACHES) | {int(self.config.extension_reach), self._product_reach()}
        )
        return {
            reach: self._compute_coverage_at(positions, total_length, reach) for reach in reaches
        }

    def _product_reach(self) -> int:
        """Measured product reach for this run's polymerase; see coverage.product_reach."""
        from .coverage import product_reach

        polymerase = getattr(self.conditions, "polymerase", None) if self.conditions else None
        return product_reach(polymerase or "phi29")

    def _compute_coverage_at(self, positions, total_length: int, reach: int) -> float:
        """Raw union coverage at an explicit reach."""
        return _union_coverage(
            positions, total_length, reach, getattr(self.config, "fg_circular", False)
        )

    def _compute_effective_coverage(self, positions_by_primer, total_length: int) -> float:
        """Coverage weighted by how much of the time each site is occupied.

        `_compute_coverage` unions binding windows as booleans: a site either
        covers its window or does not. That is the right question for "could
        this primer reach here", and the wrong one for "will it". A primer whose
        effective Tm sits well below the reaction temperature is mostly not
        bound, and its sites contribute far less amplification than a count
        implies.

        Here a site covers its window with probability theta -- the same
        two-state occupancy the selectivity metric uses -- so a base is covered
        unless every site reaching it fails to bind:

            P(covered at x) = 1 - PRODUCT over sites s reaching x of (1 - theta_s)

        Independence across sites is an approximation. Sites on one template
        molecule compete for polymerase and are not independent, so this reads
        as an upper bound on what a single molecule does; across the many
        molecules in a reaction it is the right shape.

        Returns 0.0 when no reaction conditions are attached, since without them
        there is no temperature at which to evaluate occupancy and a fabricated
        number here would be indistinguishable from a measured one.
        """
        import numpy as np

        from .coverage import _mark_window
        from .occupancy import site_occupancy
        from .thermodynamics import calculate_enthalpy_entropy

        if not positions_by_primer or total_length == 0 or self.conditions is None:
            return 0.0

        extension_reach = self.config.extension_reach
        circular = getattr(self.config, "fg_circular", False)
        temp = self.conditions.temp

        # Accumulate the probability that a base is NOT reached by anything.
        not_covered = np.ones(total_length, dtype=np.float32)
        window = np.zeros(total_length, dtype=bool)

        for primer, positions in positions_by_primer.items():
            if not positions:
                continue
            tm = self.conditions.calculate_effective_tm(primer)
            dh, _ = calculate_enthalpy_entropy(primer)
            theta = site_occupancy(dh, tm, temp)
            if theta <= 0.0:
                continue

            # One primer's sites are marked together, then applied once. Marking
            # per site would multiply (1 - theta) in for every overlapping window
            # of the SAME primer, which understates coverage where a primer binds
            # densely -- the windows overlap, the primer does not stack with
            # itself.
            window[:] = False
            for pos in positions:
                _mark_window(window, int(pos), extension_reach, total_length, circular)
            not_covered[window] *= 1.0 - theta

        return float((1.0 - not_covered).sum()) / total_length

    def _compute_gaps(self, positions: List[int], total_length: int) -> List[float]:
        """Compute gaps between adjacent binding sites."""
        if len(positions) < 2:
            return [float(total_length)]
        positions = sorted(positions)
        gaps = []
        for i in range(1, len(positions)):
            gaps.append(positions[i] - positions[i - 1])
        # Add wrap-around gap for circular genomes
        gaps.append(total_length - positions[-1] + positions[0])
        return gaps

    def _gini(self, values: List[float]) -> float:
        """Compute Gini coefficient."""
        if not values:
            return 1.0
        values = sorted(values)
        n = len(values)
        total = sum(values)
        if total == 0:
            return 0.0
        cumsum = 0
        gini_sum = 0
        for i, v in enumerate(values):
            cumsum += v
            gini_sum += (2 * (i + 1) - n - 1) * v
        return gini_sum / (n * total)

    def _estimate_tm(self, primer: str) -> float:
        """Melting temperature under the configured reaction conditions.

        This was the Wallace rule, `4*GC + 2*AT` -- a 1979 approximation for
        14-20mers at 1 M salt, with no term for magnesium, additives or the
        nearest-neighbour stacking that sets the real value. Reported for
        equiphi29 12-mer designs it reads about 11 C low: one validated
        Prevotella set whose salt-corrected Tm is 45.9 C was published as
        34.9 C, below its own 42 C reaction temperature.

        That is the single check an equiphi29 design exists to pass -- do the
        primers melt above the reaction -- so the wrong number did not merely
        mislead, it inverted the answer. It also disagreed with the Tm the
        occupancy model was using on the same primers at the same moment
        (`conditions.calculate_effective_tm`, via `occupancy.site_occupancy`),
        which is the real one. One quantity, two definitions, is the failure
        this codebase has produced most often; this removes an instance of it.

        Without conditions there is no buffer to correct for, but that is not a
        reason to fall back to something worse than nearest-neighbour.
        """
        from neoswga.core.thermodynamics import calculate_tm_basic

        if self.conditions is not None:
            return self.conditions.calculate_effective_tm(primer)
        return calculate_tm_basic(primer)

    def _estimate_dimer_risk(self, primers: List[str]) -> float:
        """
        Estimate dimer risk for primer set.

        Returns fraction of primer pairs with potential dimer formation.
        Uses the centralized DimerValidator for consistent results across
        all optimizers.
        """
        if len(primers) < 2:
            return 0.0

        try:
            from .dimer_validator import DimerValidator

            validator = DimerValidator(max_dimer_bp=self.config.max_dimer_bp)
            return validator.dimer_risk(primers)
        except ImportError:
            logger.warning("Dimer module not available, using estimate")
            return 0.1  # Fallback estimate


class CompositeOptimizer(BaseOptimizer):
    """
    Optimizer that combines results from multiple sub-optimizers.

    Useful for ensemble approaches or multi-stage optimization.
    """

    def __init__(
        self,
        optimizers: List[BaseOptimizer],
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths, bg_prefixes, bg_seq_lengths, config
        )
        self.optimizers = optimizers

    @property
    def name(self) -> str:
        return "composite"

    @property
    def description(self) -> str:
        names = [opt.name for opt in self.optimizers]
        return f"Composite optimizer combining: {', '.join(names)}"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        application: str = "balanced",
        **kwargs,
    ) -> OptimizationResult:
        """Run all sub-optimizers and return the best result.

        Selection is by ``normalized_score`` (a [0,1] value comparable across
        optimizer types), NOT raw ``score`` — raw scores are on
        per-optimizer scales and are not comparable.
        """
        from dataclasses import replace as _dc_replace

        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        best_result = None
        best_norm = -1.0

        for optimizer in self.optimizers:
            try:
                result = optimizer.optimize(
                    candidates, target, fixed_primers=fixed_primers, **kwargs
                )
                norm = result.metrics.normalized_score(application=application)
                if best_result is None or norm > best_norm:
                    best_result = result
                    best_norm = norm
            except Exception as e:
                logger.warning(f"Optimizer {optimizer.name} failed: {e}")

        if best_result is None:
            return OptimizationResult.failure(self.name, "All sub-optimizers failed")

        return best_result
