"""One contract for where a design's candidates come from.

Phase 4 increment 1 of the plan for `docs/validation/pipeline_audit_2026-09-16/`,
audit finding F1.

`plan-pool`, `optimize` and `expand-primers` each read `step3_df.csv` and build a
list. The SQLite inventory beside it holds every candidate that cleared hard QC
-- 491,836 on the Wolbachia design against the CSV's 2,000 -- and nothing in
production reads it. `iter_eligible` has exactly one caller, `CandidateProvider`,
whose only construction site is inside `design_sweep`, which no command calls.

So the retained candidates are stored, indexed, and unable to affect any
delivered panel. That is what this phase closes, and a shared contract is the
first step: one thing a design asks for candidates, whether they come from a
file the user named or from the inventory the filter wrote.

This increment deliberately changes no delivered panel. The source starts at
exactly the shortlist the CSV held, in the order the inventory records, and
expansion is a later increment. What it adds is the seam, the counts, and the
honesty about what a run actually searched: a design that examined 2,000 of
491,836 candidates should say so rather than leaving a reader to assume the
shortlist was the universe.
"""

import pandas as pd
import pytest

from neoswga.core.candidate_inventory import record_stage2_inventory
from neoswga.core.candidate_source import (
    ListCandidateSource,
    as_candidate_source,
    open_candidate_source,
)

COND = "tm-2026-09-14:abc123"
A, C, G, T = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"


def _frame(sequences):
    return pd.DataFrame(
        {
            "primer": list(sequences),
            "fg_count": [10] * len(sequences),
            "bg_count": [1] * len(sequences),
        }
    )


def _inventory(tmp_path, cleared, shortlisted):
    record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=_frame(cleared),
        after_gini=_frame(cleared),
        shortlisted=_frame(shortlisted),
        indexed=list(cleared),
        enumerated=10_000,
    )


# -- the list source, which every existing caller gets ---------------------


def test_a_list_source_hands_back_what_it_was_given():
    source = ListCandidateSource([A, C, G])

    assert source.initial(10) == [A, C, G]
    assert source.universe_size() == 3
    assert source.exhausted() is True


def test_a_list_source_respects_the_frontier_limit():
    source = ListCandidateSource([A, C, G, T])

    assert source.initial(2) == [A, C]
    assert source.examined() == 2


def test_a_list_source_cannot_expand():
    """It has no universe beyond what it was handed, and says so."""
    source = ListCandidateSource([A, C])
    source.initial(2)

    assert source.advance() is False
    assert source.exhausted() is True


def test_a_bare_list_becomes_a_source():
    """So every current caller and test keeps working unchanged."""
    source = as_candidate_source([A, C])

    assert isinstance(source, ListCandidateSource)
    assert as_candidate_source(source) is source


# -- the inventory source --------------------------------------------------


def test_the_inventory_source_starts_at_the_shortlist(tmp_path):
    """Increment 1's whole promise: the same starting pool as the CSV."""
    _inventory(tmp_path, [A, C, G, T], shortlisted=[A, C])
    source = open_candidate_source(tmp_path, COND, [12], frontier=2)

    assert sorted(source.initial(2)) == sorted([A, C])


def test_the_inventory_source_knows_the_universe_it_did_not_search(tmp_path):
    """The number a reader needs to tell a shortlist from a universe."""
    _inventory(tmp_path, [A, C, G, T], shortlisted=[A, C])
    source = open_candidate_source(tmp_path, COND, [12], frontier=2)
    source.initial(2)

    assert source.universe_size() == 4
    assert source.examined() == 2
    assert source.exhausted() is False


def test_the_inventory_source_reports_where_it_came_from(tmp_path):
    """A run should say which of the two it used rather than leave it inferred."""
    _inventory(tmp_path, [A, C, G, T], shortlisted=[A, C])
    source = open_candidate_source(tmp_path, COND, [12], frontier=2)

    described = source.describe()
    assert described["kind"] == "inventory"
    assert described["universe"] == 4
    assert described["frontier"] == 2


def test_a_list_source_describes_itself_differently():
    described = ListCandidateSource([A, C]).describe()

    assert described["kind"] == "list"
    assert described["universe"] == 2


def test_the_factory_falls_back_to_the_candidates_it_is_given(tmp_path):
    """An explicit --candidates file wins, and no inventory is required."""
    source = open_candidate_source(tmp_path, COND, [12], candidates=[A, C])

    assert isinstance(source, ListCandidateSource)
    assert source.initial(10) == [A, C]


def test_the_factory_falls_back_when_there_is_no_inventory(tmp_path):
    """Older run directories have none; they must still design."""
    source = open_candidate_source(tmp_path, COND, [12], candidates=[A, C, G])

    assert source.describe()["kind"] == "list"


def test_an_inventory_with_nothing_eligible_is_not_silently_empty(tmp_path):
    """An empty universe is a condition to report, not a pool to search."""
    _inventory(tmp_path, [A, C], shortlisted=[A])

    with pytest.raises(ValueError, match="no eligible candidates"):
        open_candidate_source(tmp_path, "a-condition-nobody-recorded", [12])


def test_the_source_filters_to_the_requested_length(tmp_path):
    """A 10-mer in the inventory must not reach a 12-mer design."""
    cleared = [A, C, "AAAAAAAAAA"]
    record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=_frame(cleared),
        after_gini=_frame(cleared),
        shortlisted=_frame([A]),
        indexed=list(cleared),
    )
    source = open_candidate_source(tmp_path, COND, [12], frontier=10)

    assert all(len(p) == 12 for p in source.initial(10))
