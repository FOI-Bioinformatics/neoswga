"""Selection weights that a method accepts and does not apply.

`--application` sets `tm_weight` and `uniformity_weight` from a profile, and on
the default `hybrid` method neither reaches the delivered panel. `HybridOptimizer`
hands both to a `NetworkOptimizer` that nothing ever reads back; Stage 2 is
`_network_refine`, a method on the class itself, with no Tm or uniformity term.
Measured 2026-09-21: all four profiles deliver an identical panel.

This module exists so the warning is not buried in a 200-line constructor, and
so the concept has a name to grep for when the weights are eventually wired.

See `docs/validation/selection_weights_are_inert_2026-09-21.md` and
`tests/test_selection_weights_reach_the_method.py`.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

__all__ = ["build_unread_network_optimizer", "warn_about_inert_selection_weights"]


def warn_about_inert_selection_weights(tm_weight, uniformity_weight) -> None:
    """Say that a configured weight does not reach this method's selection.

    The precedent is `--use-gpu`, which says it is not implemented rather than
    logging that acceleration is enabled. A weight that is accepted,
    documented and silently ignored is the defect Known Issue 8 tracks, and
    saying so is the minimum owed until it is either wired or removed.

    A user asking for `clinical` over `discovery` gets an identical panel from
    `hybrid` and has no way to tell from the output.
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
        "Stage 2 is _network_refine, which has no Tm or uniformity term. "
        "Measured identical panels across --application profiles. Use "
        "--optimization-method=network for a Tm-weighted selection.",
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
