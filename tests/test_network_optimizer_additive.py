"""Verify NetworkOptimizer Tm cache honours ReactionConditions additives.

Before this fix, _calculate_tm_score called calculate_primer_tm() which ignored
DMSO/betaine/etc. so Tm-weighted edges in the network optimizer couldn't
distinguish a primer at 30 C pure buffer from the same primer at 42 C + 1 M
betaine.
"""

import pytest

from neoswga.core.reaction_conditions import ReactionConditions


def _make_optimizer(conditions):
    from neoswga.core.network_optimizer import NetworkOptimizer

    # NetworkOptimizer only needs a callable .cache with get_positions for
    # most paths, but _get_primer_tm does not touch it. Pass None for brevity.
    return NetworkOptimizer(
        position_cache=None,
        fg_prefixes=[],
        bg_prefixes=[],
        fg_seq_lengths=[],
        bg_seq_lengths=[],
        reaction_temp=conditions.temp if conditions else 30.0,
        tm_weight=1.0,
        conditions=conditions,
    )


def test_additives_change_network_tm_cache():
    primer = "GCGCGCGCGCGC"

    opt_plain = _make_optimizer(ReactionConditions(temp=30.0, polymerase='phi29'))
    opt_betaine = _make_optimizer(ReactionConditions(
        temp=30.0, polymerase='phi29', betaine_m=1.0,
    ))

    tm_plain = opt_plain._get_primer_tm(primer)
    tm_betaine = opt_betaine._get_primer_tm(primer)

    assert tm_plain != tm_betaine, (
        "Betaine should shift the network optimizer's cached Tm; got identical values"
    )
    assert tm_plain > tm_betaine, (
        "Betaine is Tm-lowering; cached value should drop"
    )


def test_no_conditions_uses_legacy_melting_temp():
    """Without conditions, fall back to the legacy path (must not crash)."""
    from neoswga.core.network_optimizer import NetworkOptimizer

    opt = NetworkOptimizer(
        position_cache=None,
        fg_prefixes=[], bg_prefixes=[],
        fg_seq_lengths=[], bg_seq_lengths=[],
        reaction_temp=30.0,
        tm_weight=1.0,
        conditions=None,
    )
    tm = opt._get_primer_tm("ATCGATCGATCG")
    assert 10 < tm < 100


def test_tm_score_differs_when_additives_shift_tm():
    """Score at a boundary should differ between plain and additive-heavy runs."""
    primer = "GCATCGATCGAT"
    opt_plain = _make_optimizer(ReactionConditions(temp=30.0, polymerase='phi29'))
    opt_heavy = _make_optimizer(ReactionConditions(
        temp=30.0, polymerase='phi29', dmso_percent=8.0, betaine_m=2.0,
    ))
    s_plain = opt_plain._calculate_tm_score(primer)
    s_heavy = opt_heavy._calculate_tm_score(primer)
    assert s_plain != s_heavy, (
        f"Heavy additive cocktail should change Tm score; got {s_plain} vs {s_heavy}"
    )
