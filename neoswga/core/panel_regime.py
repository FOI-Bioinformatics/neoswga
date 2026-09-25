"""Which criterion limits a delivered panel, and which ones nobody constrained.

Two wet-lab benchmarks in `docs/validation/published_primer_sets.md` disagree
about which spacing property predicts success: worst gap and evenness separate
the winners on the *Prevotella* sets and do not on Clarke's *M. tuberculosis*
sets, where binding density separates them instead. Clarke pre-filtered the
candidate pool for even binding, so evenness had no residual variance left to
explain the outcome. The conclusion recorded there is that **what predicts is
whichever property is currently limiting**, and that compressing the properties
into one weighted number is therefore the wrong move.

This module is the diagnostic form of that conclusion. It reports each criterion
separately, says which one is closest to binding, and -- the part that makes it
honest -- says which ones have no reference to be measured against at all.

Two things it deliberately does not do.

**It invents no references.** A criterion is ranked only against a line the user
drew (a coverage target, a ratio floor, a requested panel size, an optional
density or host-site limit) or one physics drew. `max_gap` gets neither: measured
on the 18 published sets with wet-lab outcomes, every one carries a hole wider
than twice the calibrated reach, the winners included, so no reach-derived
threshold separates anything
(`docs/validation/getting_ahead_on_spacing_2026-09-18.md`). Such a criterion is
reported with readable units and never ranked. Ranking it would name the worst
hole as limiting on essentially every real design, which is a slogan rather than
a diagnostic.

**It produces no composite score.** `normalized_score` already exists for
cross-optimizer comparison and is a different thing. Every criterion here stays
separately addressable, because the two benchmarks say a fixed weighting is
wrong for one of them whichever way it is set.

Not to be merged with `PoolObjective.shortfall`, which looks similar and answers
a different question. That sums distances over CONFIGURED CONSTRAINTS ONLY, to
rank infeasible panels during a search. This reports per-criterion distance over
configured and physical references, for one delivered panel, and never sums.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

# Comparison sense of a reference.
AT_LEAST = "at_least"
AT_MOST = "at_most"

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Criterion:
    """One property of a delivered panel, with its reference if it has one.

    `slack` is relative to the reference, so criteria measured on different
    scales are comparable and the tightest one is not merely the one with the
    biggest units. Negative means the reference is not met. It is `None` exactly
    when `reference` is `None`, which is what keeps an unreferenced criterion
    out of the ranking.
    """

    name: str
    value: float | None
    reference: float | None
    direction: str | None
    slack: float | None
    units: str
    note: str


@dataclass(frozen=True)
class PanelRegime:
    """What limits this panel, and what was never constrained."""

    criteria: tuple[Criterion, ...]
    limiting: str | None
    failing: tuple[str, ...]
    unreferenced: tuple[str, ...]


def _slack(value: float | None, reference: float, direction: str) -> float | None:
    """Relative distance from a reference. Negative means unmet.

    Scaled by the reference itself, and guarded for a zero reference, so that
    missing a coverage target of 0.8 by half scores the same as missing a ratio
    floor of 5.0 by half.
    """
    if value is None or not math.isfinite(value):
        return None
    scale = abs(reference) or 1.0
    if direction == AT_LEAST:
        return (value - reference) / scale
    return (reference - value) / scale


def _referenced(
    name: str,
    value: float | None,
    reference: float,
    direction: str,
    units: str,
    note: str,
) -> Criterion:
    return Criterion(
        name=name,
        value=value,
        reference=reference,
        direction=direction,
        slack=_slack(value, reference, direction),
        units=units,
        note=note,
    )


def _reported(name: str, value: float | None, units: str, note: str) -> Criterion:
    """A criterion with no reference: reported, never ranked."""
    return Criterion(
        name=name,
        value=value,
        reference=None,
        direction=None,
        slack=None,
        units=units,
        note=note,
    )


# `bg_coverage` is 0.0 both when measured zero and when nothing measured it, so
# its zero is indistinguishable from a measurement -- the shape of Known Issues
# 5, 6 and 13. The diagnostic says so rather than laundering it. The two strand
# figures used to share this defect and no longer do: `core/strand_metrics.py`
# makes them None when uncomputed, so a 0.0 there is a measurement.
_AMBIGUOUS_AT_ZERO = "; a zero here may mean not computed rather than measured"


def _maybe_ambiguous(value: float | None, note: str) -> str:
    if value is not None and value == 0.0:
        return note + _AMBIGUOUS_AT_ZERO
    return note


def _coverage_criterion(metrics: Any, target: float) -> Criterion:
    """The coverage figure this panel is judged on, and which one it was.

    `None` from the evaluator means occupancy could not be computed, which is
    not the same as zero, so the fallback to the raw figure is recorded in the
    note rather than left for a reader to infer.
    """
    effective = getattr(metrics, "effective_fg_coverage", None)
    if effective is not None:
        return _referenced(
            "coverage",
            float(effective),
            target,
            AT_LEAST,
            "fraction",
            "occupancy-weighted",
        )
    raw = getattr(metrics, "fg_coverage", None)
    return _referenced(
        "coverage",
        None if raw is None else float(raw),
        target,
        AT_LEAST,
        "fraction",
        "raw, because occupancy could not be computed",
    )


def _hole_note(max_gap: float | None, reach: int, genome_length: int | None) -> str:
    """The worst hole in units a reader can act on.

    A bare base count says nothing without the chemistry and the genome, and
    those are exactly what decide whether a hole matters.
    """
    if max_gap is None or not math.isfinite(max_gap):
        return "not measurable on this panel"
    parts = []
    if reach:
        parts.append(f"{max_gap / reach:.1f}x the {reach} bp reach")
    if genome_length:
        parts.append(f"{100.0 * max_gap / genome_length:.2f}% of the target")
    parts.append("no reference: no reach-derived threshold separates published winners")
    return ", ".join(parts)


def _occupancy_criterion(metrics: Any) -> list[Criterion]:
    """How weakly the weakest primer in the panel is bound.

    Reported, never ranked: no floor on occupancy can be validated against the
    outcome data available here. The note carries the median as well, because
    the minimum alone cannot say whether a panel is uniformly weak or carries
    one outlier, and that difference is what a reader acts on.

    An occupancy-WEIGHTED gap statistic is deliberately not offered; see
    `strand_metrics.panel_occupancy` for why.
    """
    import statistics

    by_primer = getattr(metrics, "primer_occupancy", None) or {}
    if not by_primer:
        return []
    values = sorted(by_primer.values())
    weakest, median, strongest = values[0], statistics.median(values), values[-1]
    spread = (strongest / weakest) if weakest > 0 else float("inf")
    return [
        _reported(
            "weakest_occupancy",
            weakest,
            "fraction",
            f"median {median:.3g}, strongest {strongest:.3g}, spread {spread:.1f}x; "
            "no reference, since no validated floor on occupancy exists",
        )
    ]


def _convergent_criteria(
    metrics: Any,
    fg_prefixes: Sequence[str],
    bg_prefixes: Sequence[str],
) -> list[Criterion]:
    """The widest gap between opposite-strand sites, per genome set.

    Exponential amplification needs two sites in convergent orientation within
    the polymerase's reach; one site primes linearly at best. So this is the
    closest quantity in this codebase to the mechanism, and on the HOST it is
    what swga 2.0 approximates with `within_mean_gap_ratio` and fits against
    measured sequencing breadth.

    A genome nobody measured yields NO criterion rather than a zero, because an
    unmeasured host must not read as one whose sites are all adjacent.
    """
    from .strand_metrics import worst_convergent_gap

    stats = getattr(metrics, "strand_stats", None) or {}
    out: list[Criterion] = []
    for name, prefixes, note in (
        (
            "convergent_gap",
            fg_prefixes,
            "widest target stretch with no convergent pair, no reference",
        ),
        (
            "host_convergent_gap",
            bg_prefixes,
            "widest host stretch with no convergent pair; LARGER is better "
            "here, since a host whose sites cannot face each other amplifies "
            "little. No reference",
        ),
    ):
        value = worst_convergent_gap(stats, list(prefixes or []))
        if value is not None:
            out.append(_reported(name, value, "bp", note))
    return out


def _spacing_criteria(metrics: Any, reach: int, genome_length: int | None) -> list[Criterion]:
    """The properties with no line to compare against.

    These are the ones the two benchmarks disagree about, and the ones a future
    `PoolConstraints` would make constrainable.
    """
    max_gap = getattr(metrics, "max_gap", None)
    mean_gap = getattr(metrics, "mean_gap", None)
    return [
        _reported(
            "worst_hole",
            None if max_gap is None else float(max_gap),
            "bp",
            _hole_note(max_gap, reach, genome_length),
        ),
        _reported(
            "mean_gap",
            None if mean_gap is None else float(mean_gap),
            "bp",
            (
                f"{mean_gap / reach:.2f}x the {reach} bp reach, no reference"
                if mean_gap is not None and math.isfinite(mean_gap) and reach
                else "no reference"
            ),
        ),
        _reported(
            "evenness",
            _as_float(getattr(metrics, "gap_gini", None)),
            "Gini",
            "0 is uniform spacing, no reference: it separates winners on one "
            "benchmark and runs backwards on the other",
        ),
        _reported(
            "host_coverage",
            _as_float(getattr(metrics, "bg_coverage", None)),
            "fraction",
            _maybe_ambiguous(
                _as_float(getattr(metrics, "bg_coverage", None)),
                "the only quantity here that sees host site position, no reference",
            ),
        ),
        _reported(
            "strand_balance",
            _as_float(getattr(metrics, "strand_coverage_ratio", None)),
            "ratio",
            "1 is balanced between strands, no reference",
        ),
        _reported(
            "strand_alternation",
            _as_float(getattr(metrics, "strand_alternation_score", None)),
            "fraction",
            "adjacent sites on opposite strands, no reference",
        ),
    ]


def _as_float(value: Any) -> float | None:
    return None if value is None else float(value)


def assess_panel(
    metrics: Any,
    *,
    coverage_target: float,
    min_fg_bg_ratio: float,
    requested_size: int,
    delivered_size: int,
    coverage_reach: int,
    genome_length: int | None = None,
    min_selectivity_density: float | None = None,
    max_background_sites: int | None = None,
    fg_prefixes: Sequence[str] = (),
    bg_prefixes: Sequence[str] = (),
) -> PanelRegime:
    """Assess one delivered panel against its references.

    Args:
        metrics: Anything carrying the `PrimerSetMetrics` fields read here.
        coverage_target: From the `--application` profile.
        min_fg_bg_ratio: From the same profile.
        requested_size: What the run asked for; `num_primers` is a request.
        delivered_size: What it returned.
        coverage_reach: The reach the panel was selected and scored on.
        genome_length: Target length, for expressing a hole as a fraction.
        min_selectivity_density: Optional floor; a criterion only when set.
        max_background_sites: Optional ceiling; a criterion only when set.

    An unset optional limit creates NO criterion. Reporting it as satisfied
    would read as specificity where there is only an absent constraint, which is
    the shape of Known Issues 5, 6 and 13.
    """
    criteria: list[Criterion] = [
        _coverage_criterion(metrics, coverage_target),
        _referenced(
            "fg_bg_ratio",
            _as_float(getattr(metrics, "selectivity_ratio", None)),
            min_fg_bg_ratio,
            AT_LEAST,
            "ratio",
            "from the application profile",
        ),
        _referenced(
            "panel_size",
            float(delivered_size),
            float(requested_size),
            AT_LEAST,
            "primers",
            "a request, not a guarantee: selection stops rather than admit a " "dimerising pair",
        ),
    ]

    if min_selectivity_density is not None:
        criteria.append(
            _referenced(
                "selectivity_density",
                _as_float(getattr(metrics, "selectivity_density", None)),
                float(min_selectivity_density),
                AT_LEAST,
                "per-base ratio",
                "configured floor; comparable across background sizes",
            )
        )
    if max_background_sites is not None:
        criteria.append(
            _referenced(
                "background_sites",
                _as_float(getattr(metrics, "total_bg_sites", None)),
                float(max_background_sites),
                AT_MOST,
                "sites",
                "configured ceiling",
            )
        )

    criteria.extend(_spacing_criteria(metrics, coverage_reach, genome_length))
    criteria.extend(_convergent_criteria(metrics, fg_prefixes, bg_prefixes))
    criteria.extend(_occupancy_criterion(metrics))

    ranked = [c for c in criteria if c.slack is not None]
    limiting = min(ranked, key=lambda c: c.slack).name if ranked else None
    failing = tuple(c.name for c in ranked if c.slack is not None and c.slack < 0)
    unreferenced = tuple(c.name for c in criteria if c.reference is None)

    return PanelRegime(
        criteria=tuple(criteria),
        limiting=limiting,
        failing=failing,
        unreferenced=unreferenced,
    )


def _value_text(criterion: Criterion) -> str:
    if criterion.value is None:
        return "not measured"
    if not math.isfinite(criterion.value):
        return "unmeasurable"
    if abs(criterion.value) >= 1000:
        return f"{criterion.value:,.0f}"
    return f"{criterion.value:.3f}".rstrip("0").rstrip(".")


def format_regime(regime: PanelRegime) -> list[str]:
    """Lines for the CLI. One criterion per line, no composite figure."""
    if not regime.criteria:
        return []

    lines = [
        "",
        "=" * 72,
        "What limits this panel",
        "=" * 72,
        f"{'criterion':<20} {'value':>12} {'reference':>12} {'slack':>8}",
        "-" * 54,
    ]
    for criterion in regime.criteria:
        reference = (
            "none"
            if criterion.reference is None
            else _value_text(
                Criterion(
                    criterion.name,
                    criterion.reference,
                    None,
                    None,
                    None,
                    criterion.units,
                    "",
                )
            )
        )
        slack = "  --  " if criterion.slack is None else f"{criterion.slack:+.3f}"
        lines.append(
            f"{criterion.name:<20} {_value_text(criterion):>12} {reference:>12} {slack:>8}"
        )

    lines.append("")
    if regime.limiting:
        # `limiting` is the minimum slack, which means the widest relative MISS
        # on a panel that fails something and the least HEADROOM on one that
        # clears everything. Calling both "closest to binding" is wrong in the
        # first case, which is the case a user sees when something needs fixing.
        if regime.failing:
            lines.append(
                f"Furthest from its reference: {regime.limiting}, in relative "
                "terms. That is the criterion to move first."
            )
        else:
            lines.append(
                f"Closest to binding: {regime.limiting}. Everything is met, so "
                "this is the one with the least headroom."
            )
    if regime.failing:
        lines.append(f"Not met: {', '.join(regime.failing)}.")
    ambiguous = [c.name for c in regime.criteria if _AMBIGUOUS_AT_ZERO.strip("; ") in c.note]
    if ambiguous:
        # Otherwise the table shows a bare 0 and a reader takes it for a
        # measurement, which is the shape of Known Issues 5, 6 and 13. The note
        # already says so in the summary JSON; the terminal needs it too.
        lines.append(
            f"Reported as zero and may not have been computed: "
            f"{', '.join(ambiguous)}. Check the background index before reading "
            "either as specificity."
        )
    if regime.unreferenced:
        lines.append(
            f"No reference for: {', '.join(regime.unreferenced)}. These are "
            "reported and never ranked, because no threshold for them "
            "separates the published wet-lab winners. Read them against each "
            "other across designs rather than against a line."
        )
    lines.append(
        "Criteria are kept separate deliberately: two published benchmarks "
        "disagree about which spacing property predicts success, so a fixed "
        "weighting is wrong for one of them either way."
    )
    lines.append("=" * 72)
    return lines


def criteria_for_summary(regime: PanelRegime) -> dict:
    """The regime as plain data for `step4_improved_df_summary.json`."""
    return {
        "limiting": regime.limiting,
        "failing": list(regime.failing),
        "unreferenced": list(regime.unreferenced),
        "criteria": [
            {
                "name": c.name,
                "value": c.value,
                "reference": c.reference,
                "direction": c.direction,
                "slack": c.slack,
                "units": c.units,
                "note": c.note,
            }
            for c in regime.criteria
        ],
    }


def assess_from_parameter(
    metrics: Any,
    parameter: Any,
    *,
    delivered_size: int,
    application: str | None = None,
) -> PanelRegime | None:
    """Assess a panel using the references a pipeline run already resolved.

    Returns `None` rather than a fabricated assessment when the reach or the
    application profile cannot be resolved, so a caller never prints a limit
    that nothing established.
    """
    from neoswga.core.coverage import resolve_coverage_reach
    from neoswga.core.mechanistic_params import get_application_profile

    try:
        profile = get_application_profile(str(application or "enrichment").lower())
    except ValueError:
        return None

    polymerase = getattr(parameter, "polymerase", "phi29") or "phi29"
    try:
        reach = resolve_coverage_reach(
            polymerase, override=getattr(parameter, "coverage_reach", None)
        )
    except Exception:
        return None

    lengths: Sequence[int] = getattr(parameter, "fg_seq_lengths", []) or []
    requested = getattr(parameter, "num_primers", None)
    requested = getattr(parameter, "target_set_size", requested) or delivered_size

    return assess_panel(
        metrics,
        coverage_target=float(profile["default_target_coverage"]),
        min_fg_bg_ratio=float(profile["default_min_fg_bg_ratio"]),
        requested_size=int(requested),
        delivered_size=int(delivered_size),
        coverage_reach=int(reach),
        genome_length=int(sum(lengths)) or None,
        min_selectivity_density=getattr(parameter, "min_selectivity_density", None),
        max_background_sites=getattr(parameter, "max_background_sites", None),
        fg_prefixes=getattr(parameter, "fg_prefixes", []) or [],
        bg_prefixes=getattr(parameter, "bg_prefixes", []) or [],
    )


def log_regime(regime: PanelRegime | None) -> None:
    """Print the assessment, or nothing at all.

    `None` means a reference could not be resolved, and a report of nothing has
    to stay silent: printing a limit that nothing established is the failure
    this module's docstring is about.
    """
    if regime is None:
        return
    for line in format_regime(regime):
        logger.info(line)
