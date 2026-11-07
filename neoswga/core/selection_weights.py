"""Selection weights, and the one definition of the Tm preference.

`--application` sets `tm_weight` and `uniformity_weight` from a profile.

A Tm term **exists** on `hybrid` as of 2026-10-26: Stage 2
(`hybrid_optimizer._network_refine`) carries a Tm axis in its removal ranking,
and `tm_fitness` below is the one function that says what a good Tm is.

It is driven by `stage2_tm_weight`, which is a different key from the
`tm_weight` a profile supplies, and it is **off by default**. The separation is
the decision, not an accident. Every profile sets `tm_weight` between 0.10 and
0.30 and `--application` defaults to `enrichment`, so feeding it to Stage 2
would switch the term on for every run -- and the measurement does not support
that. On a single-length pool the term improves mean Tm deviation from 5.86 to
4.86 C for two points of coverage; on a mixed-length pool it costs coverage,
trips the Stage-1 coverage floor, and the floor discards the Tm gain along with
the refinement, leaving a panel worse on both axes. Mixed-length pools are what
the shipped k range produces. See
`docs/validation/2026-10-26-wiring-the-tm-term.md`.

So `tm_weight` and `uniformity_weight` both remain inert on `hybrid` and both
are still warned about. What changed is that there is now a supported route for
the Tm preference, named in that warning.

See `tests/test_selection_weights_reach_the_method.py`.
"""

from __future__ import annotations

import logging
import math

logger = logging.getLogger(__name__)

__all__ = [
    "OPTIMAL_TM_OFFSET_C",
    "build_unread_network_optimizer",
    "effective_tm",
    "tm_fitness",
    "warn_about_inert_selection_weights",
]

# A primer anneals best a little above the reaction temperature. The offset and
# the Gaussian width are the shape `NetworkOptimizer._calculate_tm_score` has
# always used; they are named here so the Tm preference has exactly one
# definition, now that two stages consult it.
OPTIMAL_TM_OFFSET_C = 5.0

# exp(-0.05 * d^2) halves at d ~= 3.72 C. One step in k moves an effective Tm by
# 4-7 C, so this term discriminates strongly between oligo lengths. That is the
# property that makes it worth measuring before trusting it, and the reason
# `docs/validation/selection_weights_are_inert_2026-09-21.md` would not simply
# reconnect the wire.
_TM_GAUSSIAN_COEFFICIENT = 0.05


def effective_tm(primer: str, conditions) -> float:
    """The one effective-Tm lookup the selection path uses.

    `conditions.calculate_effective_tm` is the canonical door: it is the Tm the
    occupancy model scores sites with, so a selection term that used anything
    else would rank primers against a different number from the one the metrics
    report. `base_optimizer._estimate_tm` carries the long version of why that
    matters, including the published Prevotella set whose Tm was reported 11 C
    below its own reaction temperature by a second definition.

    Without conditions there is no buffer to correct for, which is not a reason
    to fall back to something worse than nearest-neighbour.

    Lives here because `HybridOptimizer` does not inherit `BaseOptimizer` --
    its MRO is `(HybridOptimizer, ThermoScreenMixin, object)` -- so the two
    stages that now need this number have no common base to put it on. A copy
    in each would be the "one quantity, two definitions" failure again.
    """
    from neoswga.core.thermodynamics import calculate_tm_basic

    if conditions is not None:
        return conditions.calculate_effective_tm(primer)
    return calculate_tm_basic(primer)


def tm_fitness(effective_tm: float | None, reaction_temp: float | None) -> float:
    """How well one primer's effective Tm suits the reaction. 1 is best.

    The single definition of the Tm preference, consulted by
    `NetworkOptimizer._calculate_tm_score` and by `hybrid`'s Stage-2 removal
    ranking. Two copies of this shape is the "one quantity, two definitions"
    failure this codebase has produced most often, so there is one.

    Returns 1.0 when there is no reaction temperature to judge against, which
    makes the term uniform across candidates and therefore unable to rank.
    That is deliberate: an absent temperature is not evidence that a primer
    melts badly, and a Tm term with nothing to compare against must not
    fabricate a preference. It is NOT the "unknown is not zero" case, because
    the value is not reported as a measurement of anything; it is the
    multiplicative identity for a ranking that is switched off.

    `effective_tm` of None is the same situation for the same reason.
    """
    if reaction_temp is None or effective_tm is None:
        return 1.0
    delta = abs(float(effective_tm) - (float(reaction_temp) + OPTIMAL_TM_OFFSET_C))
    return math.exp(-_TM_GAUSSIAN_COEFFICIENT * delta**2)


def warn_about_inert_selection_weights(tm_weight, uniformity_weight) -> None:
    """Say that a configured weight does not reach this method's selection.

    The precedent is `--use-gpu`, which says it is not implemented rather than
    logging that acceleration is enabled. A weight that is accepted,
    documented and silently ignored is the defect Known Issue 8 tracks, and
    saying so is the minimum owed until it is either wired or removed.

    Both weights still qualify. Stage 2 gained a Tm axis on 2026-10-26, but it
    reads `stage2_tm_weight` rather than the profile's `tm_weight`, and it is
    off unless that key is set. Naming the route is the useful half of the
    warning now: a user who wants the Tm preference has somewhere to go.
    """
    inert = [
        name
        for name, value in (
            ("tm_weight", tm_weight),
            ("uniformity_weight", uniformity_weight),
        )
        if value
    ]
    if not inert:
        return
    logger.warning(
        "%s set, but %s not reach the panel this method delivers: hybrid's "
        "Stage 2 is _network_refine, whose removal ranking carries network, "
        "coverage and background axes, plus a Tm axis driven by the separate "
        "stage2_tm_weight key and off by default. Set stage2_tm_weight in "
        "params.json for a Tm-weighted Stage 2, or use "
        "--optimization-method=network, which applies both profile weights. "
        "Measured trade-off: docs/validation/2026-10-26-wiring-the-tm-term.md",
        " and ".join(inert),
        "it does" if len(inert) == 1 else "they do",
    )


def build_unread_network_optimizer(
    owner,
    position_cache,
    fg_prefixes,
    fg_seq_lengths,
    conditions,
    mechanistic_weight,
    reaction_temp,
    tm_weight,
    uniformity_weight,
    dimer_penalty,
    allow_dimer_relaxation,
    template_gc,
):
    """The `NetworkOptimizer` a `HybridOptimizer` builds and never calls.

    It lives here rather than in the constructor because the weights it is
    handed are this module's subject: they are accepted, forwarded, and
    never read back. `self.network_optimizer` stays a public attribute, so
    an external caller may still drive it deliberately.

    Imported inside the function because `network_optimizer` imports from
    `hybrid_optimizer`, and this module is imported by it.
    """
    from neoswga.core.network_optimizer import NetworkOptimizer

    return NetworkOptimizer(
        position_cache=position_cache,
        fg_prefixes=fg_prefixes,
        bg_prefixes=owner.bg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bg_seq_lengths=owner.bg_seq_lengths,
        max_extension=owner.max_extension,
        uniformity_weight=uniformity_weight,
        max_dimer_dg=owner.max_dimer_dg,
        # Propagate ReactionConditions so the inner NetworkOptimizer's
        # _get_primer_tm applies additive corrections (DMSO / betaine etc.)
        conditions=conditions,
        # Phase 13B: forward --use-mechanistic-model weight so the
        # NetworkOptimizer's scoring includes a mechanistic term.
        mechanistic_weight=mechanistic_weight,
        # CAUTION, corrected 2026-09-21. An earlier comment called this
        # "the object that performs refinement". It is not: Stage 2 is
        # `_network_refine` on THIS class, and nothing reads
        # `self.network_optimizer` back, so `tm_weight` and
        # `uniformity_weight` are inert here. Still passed, because the
        # object is public. See `core/selection_weights.py`.
        reaction_temp=reaction_temp,
        tm_weight=tm_weight,
        dimer_penalty=dimer_penalty,
        max_dimer_bp=owner.max_dimer_bp,
        allow_dimer_relaxation=allow_dimer_relaxation,
        template_gc=template_gc,
    )
