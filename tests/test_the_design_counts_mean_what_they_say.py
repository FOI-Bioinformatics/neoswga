"""Three counters reported one thing under the name of another.

Phase 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, the last
of finding F3's neighbours.

`design_counts` exists because the filter funnel conflated numbers and a 90%
ranking loss went unnoticed. Three of the six it replaced them with had the same
fault:

`counted` was `COUNT(*) FROM candidates`, and candidates are only written after
hard QC. So it reported the post-QC survivors under a name that promises the
enumerated universe, and `counts()`'s own docstring said the enumerated number
"was never stated anywhere" while purporting to state it.

`hard_qc_passed` counted every passing assessment in the database, across every
condition, every policy and every generation. After Phase 3 gave verdicts a
generation, that became a running total over the run's own history.

`examined` counted candidates carrying a stage-3 assessment. Stage 3 records
what an efficacy filter carried forward, which is a real quantity and is not
what a search evaluated.

The fix is not better arithmetic. Two of them now measure what their name says,
and the third is reported as unknown, because the number that would answer it
does not exist until the search reports it. An absent count says "nobody
measured this"; a wrong one says "3,412", and only one of those can be checked.
"""

import pandas as pd
import pytest

from neoswga.core.candidate_inventory import (
    STAGE2_INVENTORY_NAME,
    CandidateInventory,
    design_counts,
    qc_policy_fingerprint,
    record_stage2_inventory,
    record_stage3_policy,
)

A, C, G, T = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"
COND = "tm-2026-09-14:abc123"
OTHER_COND = "tm-2026-09-14:def456"
POLICY = {"max_gini": 0.7, "min_tm": 25.0}


def _frame(sequences):
    return pd.DataFrame(
        {
            "primer": list(sequences),
            "fg_count": [10] * len(sequences),
            "bg_count": [1] * len(sequences),
        }
    )


def _record(tmp_path, cleared, shortlisted=None, condition=COND, enumerated=None, policy=POLICY):
    frame = _frame(cleared)
    keep = frame if shortlisted is None else _frame(shortlisted)
    return record_stage2_inventory(
        tmp_path,
        condition_id=condition,
        cleared_hard_gates=frame,
        after_gini=frame,
        shortlisted=keep,
        indexed=list(cleared),
        qc_policy=policy,
        enumerated=enumerated,
    )


# -- counted ---------------------------------------------------------------


def test_counted_is_the_enumerated_universe_not_the_survivors(tmp_path):
    """The defect: three survivors out of a thousand read as three counted."""
    path = _record(tmp_path, [A, C, G], enumerated=1000)

    assert design_counts(path, COND)["counted"] == 1000


def test_counted_is_absent_when_nobody_recorded_it(tmp_path):
    """Unknown is reported as unknown, not as the nearest number to hand.

    Substituting the post-QC count would give a plausible figure that cannot be
    told apart from a measured one, which is how the funnel hid the loss this
    whole structure exists to expose.
    """
    path = _record(tmp_path, [A, C, G])

    assert design_counts(path, COND)["counted"] is None


def test_the_survivors_are_still_reported_under_their_own_name(tmp_path):
    """Losing the old number would trade one gap for another."""
    path = _record(tmp_path, [A, C, G], enumerated=1000)

    assert design_counts(path, COND)["assessed"] == 3


# -- hard_qc_passed --------------------------------------------------------


def test_hard_qc_passed_counts_only_this_condition(tmp_path):
    """Another chemistry's verdicts are not this design's eligible pool."""
    _record(tmp_path, [A, C, G], condition=COND)
    path = _record(tmp_path, [T], condition=OTHER_COND)

    assert design_counts(path, COND)["hard_qc_passed"] == 3
    assert design_counts(path, OTHER_COND)["hard_qc_passed"] == 1


def test_hard_qc_passed_counts_only_the_newest_generation(tmp_path):
    """A re-run supersedes; the count must follow, or it is a running total."""
    _record(tmp_path, [A, C, G])
    path = _record(tmp_path, [A])

    assert design_counts(path, COND)["hard_qc_passed"] == 1


def test_hard_qc_passed_counts_only_the_current_policy(tmp_path):
    """A verdict under other thresholds is not eligibility under these."""
    _record(tmp_path, [A, C, G], policy={"max_gini": 0.7})
    path = _record(tmp_path, [A, C], policy={"max_gini": 0.6})

    assert design_counts(path, COND)["hard_qc_passed"] == 2


def test_hard_qc_passed_agrees_with_what_is_eligible(tmp_path):
    """Guard the guard: the count and the traversal must not drift apart."""
    _record(tmp_path, [A, C, G, T])
    path = _record(tmp_path, [A, C])

    with CandidateInventory(path) as inventory:
        policy = inventory.current_policy(COND)
        eligible = list(inventory.iter_eligible(COND, [12], policy))

    assert design_counts(path, COND)["hard_qc_passed"] == len(eligible)


# -- examined --------------------------------------------------------------


def test_examined_is_unknown_until_a_search_reports_it(tmp_path):
    """Nothing counts search evaluations yet, so nothing may claim to."""
    path = _record(tmp_path, [A, C, G])

    assert design_counts(path, COND)["examined"] is None


def test_what_stage_three_carried_is_reported_under_its_own_name(tmp_path):
    """The number is real and worth having; it was only mislabelled."""
    path = _record(tmp_path, [A, C, G])
    record_stage3_policy(path, COND, carried=[A, C], rejected={}, policy="carry_forward")

    counts = design_counts(path, COND)
    assert counts["carried_to_stage3"] == 2
    assert counts["examined"] is None, "stage-3 membership is not a search evaluation"


def test_a_stage_three_rejection_is_not_carried(tmp_path):
    path = _record(tmp_path, [A, C, G])
    record_stage3_policy(
        path, COND, carried=[A], rejected={C: "amp_pred below threshold"}, policy="amp_model"
    )

    assert design_counts(path, COND)["carried_to_stage3"] == 1


# -- the line a user reads -------------------------------------------------


def test_the_reported_line_survives_an_unknown_count(tmp_path, caplog):
    """A None must not crash the log line that exists to show these."""
    import logging

    from neoswga.core.candidate_inventory import report_design_counts

    path = _record(tmp_path, [A, C, G])

    with caplog.at_level(logging.INFO):
        report_design_counts(path, COND)

    line = next(r.getMessage() for r in caplog.records if "Design counts" in r.getMessage())
    assert "counted=unknown" in line
    assert "examined=unknown" in line


def test_the_reported_line_shows_a_known_count(tmp_path, caplog):
    import logging

    from neoswga.core.candidate_inventory import report_design_counts

    path = _record(tmp_path, [A, C, G], enumerated=874596)

    with caplog.at_level(logging.INFO):
        report_design_counts(path, COND)

    line = next(r.getMessage() for r in caplog.records if "Design counts" in r.getMessage())
    assert "counted=874596" in line


@pytest.mark.parametrize("key", ["counted", "assessed", "hard_qc_passed", "shortlisted", "indexed"])
def test_every_documented_count_is_still_present(tmp_path, key):
    """The six distinctions are the point; none may quietly disappear."""
    path = _record(tmp_path, [A, C, G], enumerated=10)

    assert key in design_counts(path, COND)
