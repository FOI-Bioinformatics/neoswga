"""Stage-2's removal ranking: which primer `hybrid` drops at each step.

Extracted from `hybrid_optimizer.py` on 2026-10-26, when the Tm axis took it
past its module-size budget. The three objects here form one unit -- the two
weight constants and the two functions that spend them -- and nothing else in
the optimizer reads them, so they move together. `hybrid_optimizer` re-exports
the same objects, which is what keeps `from .hybrid_optimizer import
_rank_removal_candidates` working for the tests that reach for it.

The ranking is rank-normalised per step rather than a weighted sum of raw
values, because algebraic connectivity, a bin count, a host-site count and a
Gaussian in degrees Celsius share no scale.
"""

from __future__ import annotations

__all__ = [
    "_STAGE2_BACKGROUND_WEIGHT",
    "_STAGE2_COVERAGE_WEIGHT",
    "_rank_removal_candidates",
    "_removal_network_score",
    "prefix_floors_refinement",
    "resolve_tm_inputs",
    "tm_axis_for",
]

import logging

logger = logging.getLogger(__name__)


def tm_axis_for(primers, conditions, reaction_temp, weight, verbose=False):
    """The Tm axis's per-primer input, or nothing when the term is off.

    Returns `(tm_by_primer, effective_weight)`. An empty mapping and a zero
    weight mean the removal ranking takes the identical branch it took before
    the Tm axis existed, so a run that does not ask for the term pays nothing.

    Counted once per primer rather than inside the removal loop, for the same
    reason as the host load: how well a primer suits the reaction does not
    depend on which other primers are still in the set, so recomputing it every
    step would buy the same number at O(steps) times the cost.

    A weight with no reaction temperature is reported rather than silently
    dropped. There is no optimum to score a Tm against without one, so the term
    cannot rank, and a configured option that quietly does nothing is the defect
    Known Issue 8 tracks.

    `effective_tm` is the canonical Tm door and applies the configured
    additives; `tm_fitness` is the one definition of a good Tm. Both are shared
    with `NetworkOptimizer` and with `BaseOptimizer._estimate_tm`, which
    `HybridOptimizer` does not inherit.
    """
    from neoswga.core.selection_weights import (
        OPTIMAL_TM_OFFSET_C,
        effective_tm,
        tm_fitness,
    )

    weight = float(weight or 0.0)
    if weight <= 0.0:
        return {}, 0.0
    if reaction_temp is None:
        logger.warning(
            "stage2_tm_weight=%.2f was configured but this run has no reaction "
            "temperature, so there is no optimum to score a Tm against and the "
            "term is not applied. Set reaction_temp or a polymerase preset.",
            weight,
        )
        return {}, 0.0
    if verbose:
        logger.info(
            "  Stage-2 Tm term active: weight %.2f, optimum %.1f C",
            weight,
            reaction_temp + OPTIMAL_TM_OFFSET_C,
        )
    return {p: tm_fitness(effective_tm(p, conditions), reaction_temp) for p in primers}, weight


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
    candidates: list[tuple[str, float, int, int, float]],
    background_aware: bool,
    tm_weight: float = 0.0,
) -> tuple[str | None, float]:
    """Pick which primer Stage 2 should drop this step.

    Each candidate is `(primer, network_score, unique_bins, bg_sites,
    tm_fitness)`. Rank within the step rather than adding raw values: algebraic
    connectivity, a bin count, a host-site count and a Gaussian in degrees
    Celsius share no scale, and a fixed coefficient between them would be
    arbitrary and input-dependent. Normalising per step keeps the trade-off
    meaningful whatever the magnitudes happen to be.

    `tm_weight` arrives from the `--application` profile and was inert until
    2026-10-26: it was handed to a `NetworkOptimizer` nothing read back, so the
    Tm preference never reached a delivered panel. The axis below is where it
    now lands. At `tm_weight == 0` every weight is exactly what it was before
    that date, so a run that does not ask for the term is bit-identical.

    Returns `(None, -inf)` for an empty candidate list, which the caller reads
    as "nothing left that may be removed".
    """
    if not candidates:
        return None, -float("inf")

    net_values = [c[1] for c in candidates]
    cost_values = [c[2] for c in candidates]
    bg_values = [c[3] for c in candidates]
    tm_values = [c[4] for c in candidates]
    net_lo, net_hi = min(net_values), max(net_values)
    cost_lo, cost_hi = min(cost_values), max(cost_values)
    bg_lo, bg_hi = min(bg_values), max(bg_values)
    tm_lo, tm_hi = min(tm_values), max(tm_values)
    net_span = (net_hi - net_lo) or 1.0
    cost_span = (cost_hi - cost_lo) or 1.0
    bg_span = (bg_hi - bg_lo) or 1.0
    tm_span = (tm_hi - tm_lo) or 1.0

    # The two original weights keep their ratio to each other and share whatever
    # the background and Tm terms leave. With `bg_weight == 0` and
    # `tm_share == 0` this is the previous two-term expression exactly, which is
    # what keeps a plain `hybrid` run unchanged.
    bg_weight = _STAGE2_BACKGROUND_WEIGHT if background_aware else 0.0
    tm_share = min(max(float(tm_weight or 0.0), 0.0), 1.0)
    remaining = (1.0 - bg_weight) * (1.0 - tm_share)
    net_weight = (1.0 - _STAGE2_COVERAGE_WEIGHT) * remaining
    cov_weight = _STAGE2_COVERAGE_WEIGHT * remaining
    tm_term_weight = (1.0 - bg_weight) * tm_share

    best_primer = None
    best_combined = -float("inf")
    for primer, network_score, unique_bins, bg_sites, tm_fit in candidates:
        norm_net = (network_score - net_lo) / net_span
        # Cheap to remove == loses few unique bins == score 1.
        norm_keep = 1.0 - (unique_bins - cost_lo) / cost_span
        # Binds the host most == best to remove == score 1.
        norm_bg = (bg_sites - bg_lo) / bg_span
        # Melts furthest from the reaction's optimum == best to remove == 1.
        # Inverted relative to the others because `tm_fitness` is "higher is a
        # better primer" while every axis here is "higher is better to drop".
        norm_tm = 1.0 - (tm_fit - tm_lo) / tm_span
        combined = (
            net_weight * norm_net
            + cov_weight * norm_keep
            + bg_weight * norm_bg
            + tm_term_weight * norm_tm
        )
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


def resolve_tm_inputs(stage2_tm_weight, reaction_temp, poly_config):
    """The Stage-2 Tm term's two inputs, resolved once at construction.

    Returns `(weight, reaction_temp)`, either of which may be None, in which
    case `tm_axis_for` leaves the removal ranking exactly as it was.

    Neither input was retained on the optimizer before 2026-10-26, which is why
    `reaction_temp`'s constructor docstring has read "for Stage-2 Tm scoring"
    since long before any Stage-2 code could consult it: both ends of the
    option existed and the path did not, which is Known Issue 8's usual shape.

    An explicit temperature wins and the polymerase preset supplies the
    default, the same precedence `_thermo_criteria` uses for the Tm window.

    `stage2_tm_weight` is deliberately a different key from the `tm_weight` the
    same constructor takes. All five `--application` profiles set `tm_weight`
    between 0.10 and 0.30 and `--application` defaults to `enrichment`, so
    reading the profile here would switch the term on for every run; the
    measurement in `docs/validation/2026-10-26-wiring-the-tm-term.md` shows a
    gain on a single-length pool and none on a mixed-length one, which is the
    shape the shipped k range produces.
    """
    weight = float(stage2_tm_weight) if stage2_tm_weight is not None else None
    if reaction_temp is not None:
        return weight, float(reaction_temp)
    return weight, getattr(poly_config, "reaction_temp", None)


def prefix_floors_refinement(prefix_axes, refined_axes, eps: float = 1e-9) -> bool:
    """Whether the Stage-1 prefix should replace the refined set.

    Each axis is oriented higher-is-better, and the two tuples carry the same
    axes in the same order. The prefix floors the refinement only when it
    PARETO-DOMINATES it: no worse on every axis and strictly better on at least
    one. A refinement that is better on one axis and worse on another is
    incomparable, so it survives.

    That is the whole of the 2026-10-26 change. The floor used to compare raw
    covered bases alone, so a Stage-2 refinement that traded a little coverage
    for a better melting-temperature profile -- which is exactly what a
    `stage2_tm_weight` run asks for -- was discarded for covering less. With
    more than one axis the single-number comparison cannot express "better
    overall", and dominance is the comparison that can. See
    `docs/validation/2026-10-26-wiring-the-tm-term.md`.

    On a run that declares neither an objective nor the Tm term there is ONE
    axis, raw coverage, and dominance on one axis is `prefix > refined`, which
    is bit-for-bit the comparison this replaced. Equal on every axis is not
    domination, so a tie keeps the refinement, matching the old `else` branch.
    """
    no_worse = all(p >= r - eps for p, r in zip(prefix_axes, refined_axes, strict=True))
    strictly_better = any(p > r + eps for p, r in zip(prefix_axes, refined_axes, strict=True))
    return no_worse and strictly_better
