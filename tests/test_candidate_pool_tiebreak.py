"""The candidate cut has to be a selection, not a slice.

`step2` ranks surviving primers by `ratio` = bg_count / fg_count and keeps the
best `max_primer` of them. The key is right but it is not a total order, and in
the regime this tool is built for it is barely an order at all: a primer long
enough to have no exact background match has ratio 0/fg_count = 0, and at
k >= 11 against a distant background that is essentially every candidate.

Measured on Prevotella (3.2 Mb) against human chr21 at k=12, all 369,431
survivors tied at exactly 0.0. `sort_values` is stable, so the cut returned
whichever 10,000 happened to come first -- mean foreground count 1.29, against
9.5 for the 10,000 most abundant, and 10 primers in common between them.

Coverage is set by how many sites a set binds, so this capped coverage before
optimization started. Breaking the tie on foreground abundance, end to end at
k=12:

    coverage      7.6%  ->  40.3%
    fg sites        43  ->    528
    selectivity   0.69  ->   2.54
    max gap      326 kb -> 108 kb

Both halves improve because this was a defect rather than a frontier: the extra
sites are foreground sites, and the background count stayed at zero.
"""

import numpy as np
import pandas as pd
import pytest

from neoswga.core.pipeline import _rank_and_cut_candidates


def make_pool(rows):
    """rows: (primer, fg_count, bg_count)."""
    df = pd.DataFrame(rows, columns=["primer", "fg_count", "bg_count"])
    df["ratio"] = df["bg_count"] / df["fg_count"].replace(0, np.nan)
    df["ratio"] = df["ratio"].fillna(float("inf"))
    return df


# ----------------------------------------------------------------------
# The regression
# ----------------------------------------------------------------------


def test_a_pool_with_no_background_hits_is_ranked_by_abundance():
    """The exact shape of the measured case: every ratio is 0.0, so the primary
    key carries no information and the tie-break decides everything."""
    pool = make_pool([(f"P{i:03d}", i, 0) for i in range(1, 51)])
    assert pool["ratio"].nunique() == 1, "the premise is that the primary key ties"

    kept = _rank_and_cut_candidates(pool, 5)

    assert list(kept["fg_count"]) == [50, 49, 48, 47, 46]


def test_the_cut_is_not_the_arrival_order():
    """What the bug looked like: stable sort on an all-equal key returns the
    first N rows, whatever order they arrived in."""
    pool = make_pool([(f"P{i:03d}", 1, 0) for i in range(20)] + [("BEST", 99, 0)])

    kept = _rank_and_cut_candidates(pool, 3)

    assert "BEST" in set(kept["primer"]), "the most abundant primer was cut"


def test_selectivity_still_outranks_abundance():
    """The tie-break must not become the primary key. A primer with background
    load loses to a clean one however abundant it is -- specificity is what the
    ratio is protecting and coverage is only the tie-breaker."""
    pool = make_pool([("ABUNDANT_BUT_DIRTY", 500, 50), ("SPARSE_BUT_CLEAN", 2, 0)])

    kept = _rank_and_cut_candidates(pool, 1)

    assert list(kept["primer"]) == ["SPARSE_BUT_CLEAN"]


def test_ordering_within_distinct_ratios_is_unchanged():
    """Where the primary key does discriminate, nothing moves."""
    pool = make_pool([("A", 10, 1), ("B", 10, 5), ("C", 10, 3)])

    kept = _rank_and_cut_candidates(pool, 3)

    assert list(kept["primer"]) == ["A", "C", "B"]


# ----------------------------------------------------------------------
# Edges
# ----------------------------------------------------------------------


def test_a_primer_that_never_binds_the_target_ranks_last():
    """fg_count 0 gives ratio inf by construction upstream; the tie-break must
    not resurrect it."""
    pool = make_pool([("NEVER_BINDS", 0, 0), ("BINDS", 1, 0)])

    kept = _rank_and_cut_candidates(pool, 1)

    assert list(kept["primer"]) == ["BINDS"]


def test_a_pool_smaller_than_the_cut_is_returned_whole():
    pool = make_pool([("A", 3, 0), ("B", 1, 0)])

    kept = _rank_and_cut_candidates(pool, 100)

    assert len(kept) == 2


def test_an_empty_pool_does_not_raise():
    """step2 can reach here with everything filtered out, and an exception at
    the cut would mask the real problem upstream."""
    kept = _rank_and_cut_candidates(make_pool([]), 10)

    assert len(kept) == 0


@pytest.mark.parametrize("max_primer", [1, 2, 7])
def test_the_cut_length_is_respected(max_primer):
    pool = make_pool([(f"P{i}", i, 0) for i in range(1, 21)])

    assert len(_rank_and_cut_candidates(pool, max_primer)) == max_primer
