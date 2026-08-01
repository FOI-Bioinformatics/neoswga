"""Behavioural tests for `primer_attributes`.

This module computes the per-primer numbers the filter actually decides on:
melting temperature, binding-site frequency, and the Gini index of binding
spacing. It sat at 32%.

The Gini path matters most. `filter.get_gini` keeps a primer when
`gini.notna() & (gini < max_gini)` -- a guard written for exactly the case where
evenness cannot be measured. But a primer with no recorded positions produced
0.0, the BEST possible score, so it sailed through the evenness filter as if it
were ideally spaced, and the guard never fired.
"""

import math

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core import primer_attributes as pa
from neoswga.core.thermodynamics import reverse_complement

PRIMER = "TTGACCATGA"


@pytest.fixture
def h5_prefix(tmp_path):
    """A position file in the layout the pipeline writes."""
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_{len(PRIMER)}mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array([100, 5_000, 9_000], dtype=np.int32))
        f.create_dataset(reverse_complement(PRIMER), data=np.array([2_000, 7_000], dtype=np.int32))
    return prefix


# ----------------------------------------------------------------------
# Unmeasurable evenness must not read as perfect evenness
# ----------------------------------------------------------------------


def test_no_positions_gives_nan_not_zero():
    """0.0 means "perfectly even" and is the best score a primer can get.

    Returning it for a primer nothing is known about let it pass the evenness
    filter on an absence of evidence.
    """
    result = pa.get_gini_from_txt_for_one_k(
        [PRIMER], "nonexistent", None, 10_000, True, position_cache={}
    )
    forward, reverse = result[PRIMER]
    assert math.isnan(forward) and math.isnan(reverse)


def test_the_filters_notna_guard_now_fires():
    """The guard was already there; it was the value that was wrong."""
    import pandas as pd

    result = pa.get_gini_from_txt_for_one_k(
        [PRIMER], "nonexistent", None, 10_000, True, position_cache={}
    )
    gini = pd.Series([max(result[PRIMER])])

    kept = gini.notna() & (gini < 0.6)
    assert not kept.iloc[0], "an unmeasurable primer still passes the evenness filter"


def test_positions_without_gaps_still_score_zero():
    """One binding site is data, not an absence of it.

    There is simply no spacing to be uneven about, so 0 is the right answer and
    must not be swept up by the NaN change.
    """
    cache = {("x", PRIMER): [500]}
    result = pa.get_gini_from_txt_for_one_k([PRIMER], "x", None, 10_000, True, cache)
    forward, _reverse = result[PRIMER]
    assert forward == 0 and not math.isnan(forward)


def test_uneven_spacing_scores_higher_than_even_spacing():
    """The metric has to actually discriminate, or the threshold means nothing."""
    even = {("x", PRIMER): [1_000, 2_000, 3_000, 4_000, 5_000]}
    lumpy = {("x", PRIMER): [10, 20, 30, 40, 9_500]}

    even_gini = pa.get_gini_from_txt_for_one_k([PRIMER], "x", None, 10_000, True, even)[PRIMER][0]
    lumpy_gini = pa.get_gini_from_txt_for_one_k([PRIMER], "x", None, 10_000, True, lumpy)[PRIMER][0]

    assert lumpy_gini > even_gini


def test_gini_is_bounded():
    cache = {("x", PRIMER): [10, 20, 30, 9_900]}
    forward, reverse = pa.get_gini_from_txt_for_one_k([PRIMER], "x", None, 10_000, True, cache)[
        PRIMER
    ]
    for value in (forward, reverse):
        assert 0.0 <= value <= 1.0 or math.isnan(value)


def test_forward_and_reverse_are_scored_separately():
    """They are reported as a pair because a primer can be even on one strand
    and clustered on the other."""
    cache = {
        ("x", PRIMER): [1_000, 2_000, 3_000],
        ("x", reverse_complement(PRIMER)): [10, 20, 9_500],
    }
    forward, reverse = pa.get_gini_from_txt_for_one_k([PRIMER], "x", None, 10_000, True, cache)[
        PRIMER
    ]
    assert forward != reverse


# ----------------------------------------------------------------------
# Reading positions from HDF5
# ----------------------------------------------------------------------


def test_rate_counts_both_strands(h5_prefix):
    """A primer binds double-stranded DNA; counting one strand halves the rate.

    Forward hits are stored under the primer, reverse hits under its reverse
    complement -- the distinction that has been the source of several bugs in
    this codebase.
    """
    assert pa.get_rate_from_h5py(PRIMER, [h5_prefix]) == 5  # 3 forward + 2 reverse


def test_rate_is_zero_for_an_absent_primer(h5_prefix):
    assert pa.get_rate_from_h5py("CGCGCGCGCG", [h5_prefix]) == 0


def test_rate_sums_across_prefixes(h5_prefix, tmp_path):
    """Multi-contig and multi-genome runs pass several prefixes."""
    second = str(tmp_path / "other")
    with h5py.File(f"{second}_{len(PRIMER)}mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array([42], dtype=np.int32))

    assert pa.get_rate_from_h5py(PRIMER, [h5_prefix, second]) == 6


def test_missing_prefix_warns_rather_than_raising(caplog):
    """A missing position file must be visible, not silently counted as zero."""
    import logging

    with caplog.at_level(logging.WARNING):
        count = pa.get_rate_from_h5py(PRIMER, ["definitely_not_a_prefix"])

    assert count == 0
    assert "Cannot find HDF5" in caplog.text


def test_gini_from_h5_runs_at_all(h5_prefix):
    """It called `_utility.gini`, which does not exist -- the module exposes
    `gini_exact` -- so every invocation raised AttributeError.

    Nothing inside the package calls it (`pipeline` uses the same-named helper
    in `filter`), which is why a public function being entirely broken went
    unnoticed.
    """
    value = pa.get_gini(PRIMER, [h5_prefix])
    assert 0.0 <= value <= 1.0


def test_gini_from_h5_reads_both_strands(h5_prefix, tmp_path):
    """Forward hits live under the primer, reverse under its revcomp."""
    only_forward = str(tmp_path / "fwd")
    with h5py.File(f"{only_forward}_{len(PRIMER)}mer_positions.h5", "w") as f:
        f.create_dataset(PRIMER, data=np.array([100, 5_000, 9_000], dtype=np.int32))

    both = pa.get_gini(PRIMER, [h5_prefix])
    forward_only = pa.get_gini(PRIMER, [only_forward])
    assert both != forward_only, "reverse-strand positions made no difference"


# ----------------------------------------------------------------------
# Melting temperature
# ----------------------------------------------------------------------


def test_melting_temperature_is_plausible_for_a_short_primer():
    """SWGA primers are short and run cold; a PCR-range Tm would mean the wrong
    model is being used."""
    tm = pa.get_melting_tm(PRIMER)
    assert 0.0 < tm < 60.0


def test_gc_rich_primer_melts_higher_than_at_rich():
    at_rich = pa.get_melting_tm("ATATATATAT")
    gc_rich = pa.get_melting_tm("GCGCGCGCGC")
    assert gc_rich > at_rich


def test_longer_primer_melts_higher():
    assert pa.get_melting_tm("ACGTACGTACGTACGT") > pa.get_melting_tm("ACGTACGT")
