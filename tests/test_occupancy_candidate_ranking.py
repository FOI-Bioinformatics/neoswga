"""Ranking candidates by what binds, not by what matches exactly.

`step2` keeps the best `max_primer` survivors ranked by
`ratio` = bg_count / fg_count -- background load per unit of foreground
binding, ascending. That is the right question. Measured with exact k-mer
counts it is barely answerable: a primer long enough to have no perfect
background match scores 0/fg_count = 0, and against a distant background that
is nearly every candidate. At k=12 on Prevotella vs human chr21, all 369,431
survivors scored exactly 0.0.

Those candidates are not equivalent. Weighted by occupancy -- counting the
near-matches the model already computes -- their background load spans 664x
(0.54 to 358.68) and their fg/bg ratio spans 200x (0.057 to 11.6), none of it
visible to a count.

Discarding that capped the design. Coverage-greedy at 10 kb reach over the
abundance-ranked pool against the occupancy-ranked one:

    n     abundance-ranked      occupancy-ranked
     8    0.836 / sel 2.01      0.782 / sel 4.38
    16    0.929 / sel 1.44      0.860 / sel 4.59
    24    0.975 / sel 1.06      0.915 / sel 4.56
    32    0.998 / sel 0.95      0.952 / sel 4.55

Selectivity used to collapse as primers were added, because each new primer was
chosen for abundance and brought its background load with it. It now holds near
4.5, because every primer in the pool is clean. The frontier moved outward
rather than along: 16 occupancy-ranked primers beat 8 abundance-ranked ones on
coverage *and* on selectivity.

This is the same key the original design chose, measured properly -- not a
different objective.
"""

import numpy as np
import pandas as pd
import pytest

from neoswga.core import pipeline


@pytest.fixture
def restore_parameter():
    """These are module globals on `pipeline.parameter`; leave them as found."""
    names = (
        "fg_prefixes",
        "bg_prefixes",
        "occupancy_ranking",
        "occupancy_shortlist",
        "max_mismatches",
    )
    snapshot = {n: getattr(pipeline.parameter, n, "<absent>") for n in names}
    yield pipeline.parameter
    for name, value in snapshot.items():
        if value == "<absent>":
            if hasattr(pipeline.parameter, name):
                delattr(pipeline.parameter, name)
        else:
            setattr(pipeline.parameter, name, value)


def make_pool(rows):
    """rows: (primer, fg_count, bg_count)."""
    df = pd.DataFrame(rows, columns=["primer", "fg_count", "bg_count"])
    df["ratio"] = df["bg_count"] / df["fg_count"].replace(0, np.nan)
    df["ratio"] = df["ratio"].fillna(float("inf"))
    return df


# ----------------------------------------------------------------------
# It falls back rather than inventing a ranking
# ----------------------------------------------------------------------


def test_without_a_background_it_uses_the_exact_count_key(restore_parameter, monkeypatch):
    """Host-free designs have nothing to rank against. Returning some number
    anyway would be indistinguishable from a measured one."""
    monkeypatch.setattr(pipeline.parameter, "bg_prefixes", [], raising=False)
    pool = make_pool([("AAA", 1, 0), ("BBB", 9, 0)])

    kept = pipeline._rank_and_cut_candidates(pool, 1)

    # Falls through to ratio-then-abundance, so the abundant primer wins.
    assert list(kept["primer"]) == ["BBB"]
    assert "occupancy_ratio" not in kept.columns


def test_missing_count_files_fall_back_and_say_so(restore_parameter, monkeypatch, caplog):
    """`weighted_site_load` raises FileNotFoundError when the jellyfish files
    for a mismatch class are absent. A silent fallback would leave a user
    believing the pool was ranked on occupancy when it was not."""
    import logging

    monkeypatch.setattr(pipeline.parameter, "fg_prefixes", ["/nonexistent/fg"], raising=False)
    monkeypatch.setattr(pipeline.parameter, "bg_prefixes", ["/nonexistent/bg"], raising=False)
    pool = make_pool([("AAAACCCCGGGG", 4, 0), ("TTTTGGGGAAAA", 9, 0)])

    with caplog.at_level(logging.INFO):
        kept = pipeline._rank_and_cut_candidates(pool, 2)

    assert "occupancy_ratio" not in kept.columns
    assert "exact counts" in caplog.text


def test_it_can_be_turned_off(restore_parameter, monkeypatch, caplog):
    """The ranking changes every design against a background. A user
    reproducing an older result needs a way back."""
    import logging

    monkeypatch.setattr(pipeline.parameter, "bg_prefixes", ["bg"], raising=False)
    monkeypatch.setattr(pipeline.parameter, "occupancy_ranking", False, raising=False)
    pool = make_pool([("AAA", 1, 0), ("BBB", 9, 0)])

    with caplog.at_level(logging.INFO):
        kept = pipeline._rank_and_cut_candidates(pool, 2)

    assert "occupancy_ratio" not in kept.columns
    assert "disabled" in caplog.text


def test_an_empty_pool_does_not_raise(restore_parameter):
    assert len(pipeline._rank_and_cut_candidates(make_pool([]), 10)) == 0


# ----------------------------------------------------------------------
# The ranking itself
# ----------------------------------------------------------------------


def test_it_ranks_by_background_load_per_unit_of_foreground_binding():
    """Direction check. `ratio` is ascending-is-better, and the occupancy
    column has to match, or a reader comparing the two columns is comparing an
    ascending key with a descending one."""
    calls = {}

    class FakeConditions:
        temp = 42.0

    def fake_load(primers, prefixes, conditions, max_mismatches):
        primer = primers[0]
        calls.setdefault(primer, []).append(prefixes[0])
        table = {"CLEAN": (100.0, 1.0), "DIRTY": (100.0, 50.0)}
        fg, bg = table[primer]
        return fg if prefixes[0] == "fg" else bg

    import neoswga.core.occupancy as occ

    original = occ.weighted_site_load
    occ.weighted_site_load = fake_load
    try:
        ratios = pipeline._occupancy_ratio_column(
            ["CLEAN", "DIRTY"], ["fg"], ["bg"], FakeConditions()
        )
    finally:
        occ.weighted_site_load = original

    assert ratios == [0.01, 0.5], "lower must mean cleaner, matching `ratio`"


def test_a_primer_that_never_binds_the_target_ranks_last():
    """Mirrors the exact-count path, which assigns float('inf') when fg_count
    is 0. Dividing by it instead would raise or produce a spuriously good
    score."""

    class FakeConditions:
        temp = 42.0

    import neoswga.core.occupancy as occ

    original = occ.weighted_site_load
    occ.weighted_site_load = lambda primers, prefixes, c, m: (0.0 if prefixes[0] == "fg" else 5.0)
    try:
        ratios = pipeline._occupancy_ratio_column(["NEVER"], ["fg"], ["bg"], FakeConditions())
    finally:
        occ.weighted_site_load = original

    assert ratios == [float("inf")]


# ----------------------------------------------------------------------
# Configuration is type-checked
# ----------------------------------------------------------------------


def test_a_mocked_parameter_module_does_not_shrink_the_shortlist(monkeypatch):
    """The trap this codebase has hit twice. Several test modules replace the
    `parameter` module with a MagicMock, whose every attribute is a truthy stub
    with `int(...) == 1`. A `getattr(...) or default` would silently rank a
    one-primer shortlist."""
    from unittest.mock import MagicMock

    monkeypatch.setattr(pipeline, "parameter", MagicMock())

    assert pipeline._configured_positive_int("occupancy_shortlist", 50_000) == 50_000
    assert pipeline._configured_positive_int("max_mismatches", 1) == 1


@pytest.mark.parametrize("bad", ["many", -1, 0, None, 2.5, True])
def test_a_nonsense_shortlist_falls_back(restore_parameter, monkeypatch, bad):
    """A hand-edited params.json should not quietly produce a different pool.
    `True` is in the list because bool is an int in Python."""
    monkeypatch.setattr(pipeline.parameter, "occupancy_shortlist", bad, raising=False)

    assert pipeline._configured_positive_int("occupancy_shortlist", 50_000) == 50_000


def test_the_shortlist_is_never_smaller_than_the_cut(restore_parameter, monkeypatch):
    """A shortlist below `max_primer` would discard candidates before the
    ranking that is supposed to choose between them, making the cut arbitrary
    again in exactly the way this change exists to fix."""
    monkeypatch.setattr(pipeline.parameter, "occupancy_shortlist", 5, raising=False)
    monkeypatch.setattr(pipeline.parameter, "bg_prefixes", ["bg"], raising=False)
    monkeypatch.setattr(pipeline.parameter, "fg_prefixes", ["fg"], raising=False)

    seen = {}

    def fake_ratio(primers, fg, bg, conditions):
        seen["n"] = len(primers)
        return [1.0] * len(primers)

    monkeypatch.setattr(pipeline, "_occupancy_ratio_column", fake_ratio)
    monkeypatch.setattr(pipeline, "_occupancy_ranking_inputs", lambda: (object(), ["bg"]))

    pool = make_pool([(f"P{i:03d}", i + 1, 0) for i in range(40)])
    pipeline._rank_and_cut_candidates(pool, 20)

    assert seen["n"] >= 20
