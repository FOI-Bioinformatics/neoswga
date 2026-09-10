"""Evenness of spacing is not measurable from one or two binding sites.

Audit finding B4. `filter.get_gini` keeps a primer when
`gini.notna() & (gini < max_gini)`. The `.notna()` half was written for exactly
the case where evenness cannot be measured, but a primer with a single site
produced 0.0 -- the BEST score available -- so the guard never fired and the
evenness gate ranked that primer first.

How much this bit depends on the pool. On the shipped Prevotella-against-chr21
pool 86% of the 10,000 kept primers scored exactly 0.0, and every one of those
had two or fewer foreground sites. On the plasmid example the figure is 96.2%
of 500 rows. On the three whole-genome GC-tier pools, which reach step 3
already cut to a few hundred abundant primers, no primer scores 0.0 at all.
So this is a defect of small targets and of large sparse pools, and it is
invisible in the whole-genome runs this project usually inspects.
"""

import math

import pytest

from neoswga.core import primer_attributes as pa

PRIMER = "TTGACCATGA"
RC = "TCATGGTCAA"  # reverse complement of PRIMER


def _gini(forward, reverse=(), circular=True, seq_length=10_000, min_sites=None):
    cache = {("x", PRIMER): list(forward), ("x", RC): list(reverse)}
    return pa.get_gini_from_txt_for_one_k(
        [PRIMER], "x", None, seq_length, circular, cache, min_sites=min_sites
    )[PRIMER]


def test_the_default_minimum_is_three_sites():
    """Pinned so the default is a stated decision rather than a magic number."""
    assert pa.DEFAULT_MIN_GINI_SITES == 3


def test_the_threshold_is_a_parameter_not_a_literal():
    """Task 1b routes this from params.json and --min-gini-sites.

    Two sites are unmeasurable at the default and measurable at a threshold of
    2, so the argument has to be reaching the comparison.
    """
    forward, _ = _gini([500, 6_000], min_sites=2)
    assert not math.isnan(forward)
    forward, _ = _gini([500, 6_000], min_sites=3)
    assert math.isnan(forward)


def test_a_single_site_is_not_measurable():
    """The regression. One site scored 0.0, the best value, and passed the gate."""
    forward, reverse = _gini([500])
    assert math.isnan(forward) and math.isnan(reverse)


def test_two_sites_are_not_measurable():
    """Two sites give one gap per strand, whose Gini is identically 0.0."""
    forward, reverse = _gini([500, 6_000])
    assert math.isnan(forward) and math.isnan(reverse)


def test_three_sites_are_measurable():
    forward, reverse = _gini([500, 4_000, 8_000])
    assert not math.isnan(forward)
    assert 0.0 <= forward <= 1.0


def test_sites_are_counted_across_both_strands():
    """Two forward and one reverse is three sites, so the primer is measurable.

    Counting per strand instead would reject every primer that binds one strand
    only, which is a different and larger change.
    """
    forward, reverse = _gini([500, 4_000], reverse=[6_000])
    assert not math.isnan(forward)


def test_no_positions_at_all_is_still_nan():
    """The narrower guard this replaces must not be lost."""
    forward, reverse = _gini([])
    assert math.isnan(forward) and math.isnan(reverse)


def test_the_filters_notna_guard_rejects_an_unmeasurable_primer():
    """The point of the change: the existing guard now fires."""
    import pandas as pd

    forward, reverse = _gini([500])
    gini = pd.Series([max(forward, reverse)])
    kept = gini.notna() & (gini < 0.6)
    assert not kept.iloc[0]


def test_an_uneven_primer_still_scores_above_an_even_one():
    """The metric has to discriminate among the primers it does measure."""
    even, _ = _gini([1_000, 2_000, 3_000, 4_000, 5_000])
    lumpy, _ = _gini([10, 20, 30, 40, 9_500])
    assert lumpy > even
