"""The heterodimer screen examines the pairs a cheap test flags, once.

Finding A1. The screen materialised all n(n-1)/2 pairs before checking any of
them (2.0 GB measured for 2789 primers) and ran a 627 us thermodynamic
calculation on every one of them. It now streams, and only over the pairs the
exact substring test flags.
"""

import random

import pytest

from neoswga.core.thermodynamic_filter import ThermodynamicCriteria, ThermodynamicFilter


def _pool(n, seed=0):
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(12)) for _ in range(n)]


def test_screen_calls_the_thermodynamic_check_only_for_flagged_pairs(monkeypatch):
    primers = _pool(60, seed=3)
    criteria = ThermodynamicCriteria(min_tm=0.0, max_tm=100.0, min_gc=0.0, max_gc=1.0)
    filt = ThermodynamicFilter(criteria)

    calls = []
    import neoswga.core.thermodynamic_filter as tf

    real = tf.check_heterodimer
    monkeypatch.setattr(
        tf, "check_heterodimer", lambda a, b, c: calls.append((a, b)) or real(a, b, c)
    )

    filt.filter_candidates(primers, check_heterodimers=True, max_dimer_bp=3)

    from neoswga.core import dimer_matrix

    flagged = int(dimer_matrix.build(primers, 3).pairs.sum()) // 2
    total = len(primers) * (len(primers) - 1) // 2
    assert len(calls) <= flagged
    assert len(calls) < total, "the screen must not visit every pair"


def test_screen_does_not_build_a_full_pair_list(monkeypatch):
    """A list comprehension over all pairs is the allocation this removes."""
    import inspect

    import neoswga.core.thermodynamic_filter as tf

    source = inspect.getsource(tf.ThermodynamicFilter.filter_candidates)
    assert "for j in range(i + 1, len(passing))" not in source, (
        "the eager pair comprehension is back"
    )


def test_screen_falls_back_to_the_substring_test_when_the_matrix_would_be_too_big(
    monkeypatch, caplog
):
    """dimer_matrix.build raises ValueError when max_dimer_bp needs more t-mer
    codes than it will allocate (params.schema.json permits max_dimer_bp up to
    15 and max_k up to 30, so a threshold/length combination that trips this is
    reachable with legal configuration, e.g. a 12-base primer pool with
    max_dimer_bp=10). The screen must not crash; it must fall back to flagging
    pairs with dimer.is_dimer_fast pairwise rather than screening every pair
    thermodynamically, and it must say so at warning level.
    """
    primers = _pool(30, seed=7)
    criteria = ThermodynamicCriteria(min_tm=0.0, max_tm=100.0, min_gc=0.0, max_gc=1.0)
    filt = ThermodynamicFilter(criteria)

    import neoswga.core.thermodynamic_filter as tf

    calls = []
    real = tf.check_heterodimer
    monkeypatch.setattr(
        tf, "check_heterodimer", lambda a, b, c: calls.append((a, b)) or real(a, b, c)
    )

    total = len(primers) * (len(primers) - 1) // 2

    with caplog.at_level("WARNING"):
        # max_dimer_bp=10 on 12-mers needs 4**11 t-mer codes, above
        # dimer_matrix.MAX_CODES (4**8), so dimer_matrix.build raises.
        result, stats = filt.filter_candidates(
            primers, check_heterodimers=True, max_dimer_bp=10
        )

    assert result is not None
    assert len(calls) <= total
    warning_text = " ".join(r.message for r in caplog.records if r.levelno >= 30)
    assert "10" in warning_text
    assert "substring" in warning_text.lower()
