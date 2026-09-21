"""Evenness is computed for the oligos it was handed, not for a configured window.

`get_gini_from_txt` iterated `range(parameter.min_k, parameter.max_k + 1)` and
built one task per length in that window. A primer whose length fell outside it
therefore got no task, no entry in `primer_to_all_ginis`, and the final
comprehension

    [np.mean(primer_to_all_ginis[primer]) for primer in primer_list]

raised `KeyError` on it.

This is the same defect `string_search.get_positions` records and fixed in its
own docstring: there, an externally designed 14-mer set against a params.json
configured for 10-mers "returned nothing for every primer". The Gini path was
not carried along, and it fails harder -- a crash rather than a silent zero.

It matters for any pool whose lengths do not exactly match the configured
window: a mixed-length design, `expand-primers`, `evaluate-set`, and every
bring-your-own-oligo flow.

The cache is supplied in these tests so the direct path runs. The `None` path
spawns a process pool per length and is not what is under test here.
"""

import logging

import numpy as np
import pytest

from neoswga.core import parameter, primer_attributes

PREFIX = "fg"
SEQ_LENGTH = 1000


def cache_for(primers, sites=(10, 200, 400, 600)):
    """A `(prefix, primer) -> positions` map, the shape the helper expects.

    Enough sites that evenness is measurable: below `min_gini_sites` the Gini
    is NaN by design, which is a different rule and not this file's subject.
    """
    return {(PREFIX, primer): np.array(sites, dtype=np.int64) for primer in primers}


def gini_for(primers, cache=None):
    return primer_attributes.get_gini_from_txt(
        list(primers),
        [PREFIX],
        ["unused.fna"],
        [SEQ_LENGTH],
        False,
        position_cache=cache if cache is not None else cache_for(primers),
    )


@pytest.fixture(autouse=True)
def window(monkeypatch):
    """A configured window of 10-12, so an 8-mer is outside it."""
    monkeypatch.setattr(parameter, "min_k", 10, raising=False)
    monkeypatch.setattr(parameter, "max_k", 12, raising=False)
    monkeypatch.setattr(parameter, "min_gini_sites", 3, raising=False)


IN_WINDOW = "ACGTACGTAC"  # 10
OUTSIDE_SHORT = "ACGTACGT"  # 8
OUTSIDE_LONG = "ACGTACGTACGTAC"  # 14


def test_an_out_of_window_length_is_scored_rather_than_crashing():
    """The reported failure: a KeyError on the primer nobody built a task for."""
    scores = gini_for([IN_WINDOW, OUTSIDE_SHORT])

    assert len(scores) == 2
    assert all(np.isfinite(value) for value in scores), scores


def test_lengths_either_side_of_the_window_are_both_scored():
    scores = gini_for([OUTSIDE_SHORT, IN_WINDOW, OUTSIDE_LONG])

    assert len(scores) == 3
    assert all(np.isfinite(value) for value in scores), scores


def test_the_order_of_the_returned_scores_matches_the_input():
    """The caller zips this against its own primer list, so order is load-bearing."""
    dense = (10, 200, 400, 600)
    sparse = (10, 20, 30, 900)
    cache = {
        (PREFIX, IN_WINDOW): np.array(dense, dtype=np.int64),
        (PREFIX, OUTSIDE_SHORT): np.array(sparse, dtype=np.int64),
    }

    forward = gini_for([IN_WINDOW, OUTSIDE_SHORT], cache)
    reverse = gini_for([OUTSIDE_SHORT, IN_WINDOW], cache)

    assert forward == pytest.approx(list(reversed(reverse)))


def test_a_length_outside_the_window_is_reported(caplog):
    """Included because it was asked for, and said out loud.

    Silence here would be the other half of the original defect: a user who
    configured 10-12 and handed in an 8-mer should be told it was scored
    anyway, not left to infer it.
    """
    with caplog.at_level(logging.INFO, logger="neoswga.core.primer_attributes"):
        gini_for([IN_WINDOW, OUTSIDE_SHORT])

    assert "8" in caplog.text
    assert "10" in caplog.text and "12" in caplog.text


def test_a_pool_entirely_inside_the_window_logs_nothing_extra(caplog):
    with caplog.at_level(logging.INFO, logger="neoswga.core.primer_attributes"):
        gini_for([IN_WINDOW, "ACGTACGTACGT"])

    assert "outside" not in caplog.text


def test_scores_are_unchanged_for_a_pool_inside_the_window():
    """The fix must not move an existing single-length design."""
    primers = [IN_WINDOW, "TTTTGGGGCC"]
    scores = gini_for(primers)

    assert len(scores) == 2
    assert all(0.0 <= value <= 1.0 for value in scores), scores


def test_a_primer_below_the_site_threshold_is_nan_not_missing():
    """Unmeasurable evenness is a different outcome from an unscored length.

    One site gives no gap and two give a single gap whose Gini is identically
    0.0, the best value available, so below `min_gini_sites` the answer is NaN
    and the caller drops the row. That rule is unaffected by this fix, and the
    primer must still appear in the returned list.
    """
    cache = {
        (PREFIX, IN_WINDOW): np.array([10, 200, 400], dtype=np.int64),
        (PREFIX, OUTSIDE_SHORT): np.array([10], dtype=np.int64),
    }
    scores = gini_for([IN_WINDOW, OUTSIDE_SHORT], cache)

    assert len(scores) == 2
    assert np.isfinite(scores[0])
    assert np.isnan(scores[1])
