"""Under `post_gini`, eligibility and the index have to be the same set.

Phase 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, audit
finding F5 part B, and the "record assessed failures explicitly" item with it.

`record_stage2_inventory` marked every hard-QC survivor eligible regardless of
retention mode. Under `post_gini` only the candidates clearing the evenness gate
get a background position index, so the inventory said 491,836 candidates were
selectable while 18,905 of them had any background data. A design reaching one
of the others would score it against an absent index, which reads as perfect
specificity rather than as a missing measurement -- the shape of Known Issues 5,
6 and 13, and the very thing retention was introduced to close.

Two readings were available and they are not interchangeable. `post_gini` could
be an INDEXING BUDGET, in which case the provider must index a candidate before
evaluating it; or a HARD ADMISSION POLICY, in which case those candidates are
not eligible at all. It was documented as the second and behaved as neither.

It is now the second. The eligible set and the indexed set are the same set, and
a candidate outside it carries an explicit failed assessment naming the gate,
so "assessed and failed" is distinguishable from "not yet assessed" -- the
second of which renders as "no eligible candidates" and means something quite
different.

Switching mode changes who qualifies, so the mode is part of the admission
policy digest. A run under one mode must not hand its verdicts to a run under
the other.
"""

import pandas as pd

from neoswga.core.candidate_inventory import (
    CandidateInventory,
    design_counts,
    record_stage2_inventory,
)
from neoswga.core.qc_policy import resolved_qc_policy

A, C, G, T = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"
COND = "tm-2026-09-14:abc123"

# A, C clear the evenness gate; G, T do not.
CLEARED = [A, C, G, T]
AFTER_GINI = [A, C]


def _frame(sequences):
    return pd.DataFrame(
        {
            "primer": list(sequences),
            "fg_count": [10] * len(sequences),
            "bg_count": [1] * len(sequences),
        }
    )


def _record(tmp_path, retention):
    """Record as the filter does, indexing exactly what the mode indexes."""
    indexed = AFTER_GINI if retention == "post_gini" else CLEARED
    return record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=_frame(CLEARED),
        after_gini=_frame(AFTER_GINI),
        shortlisted=_frame([A]),
        indexed=list(indexed),
        qc_policy={"max_gini": 0.7, "candidate_retention": retention},
        retention=retention,
    )


def _eligible(path):
    with CandidateInventory(path) as inventory:
        policy = inventory.current_policy(COND)
        return sorted(inventory.iter_eligible(COND, [12], policy))


def test_post_gini_admits_only_the_evenness_survivors(tmp_path):
    """The defect: four eligible, two indexed."""
    path = _record(tmp_path, "post_gini")

    assert _eligible(path) == sorted(AFTER_GINI)


def test_all_qc_still_admits_every_hard_qc_survivor(tmp_path):
    """The other mode is unchanged; that is its whole purpose."""
    path = _record(tmp_path, "all_qc")

    assert _eligible(path) == sorted(CLEARED)


def test_what_is_eligible_is_what_is_indexed(tmp_path):
    """The property, stated directly rather than through either mode.

    A candidate the search may reach and the index does not hold is a zero that
    reads as a measurement.
    """
    for retention in ("all_qc", "post_gini"):
        directory = tmp_path / retention
        directory.mkdir()
        path = _record(directory, retention)

        with CandidateInventory(path) as inventory:
            eligible = set(inventory.iter_eligible(COND, [12], inventory.current_policy(COND)))
            indexed = {
                sequence
                for sequence in CLEARED
                if inventory.metrics(sequence).get("indexed")
            }

        assert eligible == indexed, f"{retention}: eligible and indexed disagree"


def test_a_rejected_candidate_says_why(tmp_path):
    """Assessed-and-failed must be distinguishable from never-assessed.

    The second renders as "no eligible candidates", which is a different claim
    and sends a reader looking for a different problem.
    """
    path = _record(tmp_path, "post_gini")

    with CandidateInventory(path) as inventory:
        history = inventory.assessment_history(G, COND)

    assert history, "the rejected candidate has no assessment at all"
    assert history[0]["passed"] is False
    assert any("evenness" in reason for reason in history[0]["reasons"])


def test_a_rejected_candidate_is_still_recorded_as_enumerated(tmp_path):
    """It failed a gate; it did not stop existing."""
    path = _record(tmp_path, "post_gini")

    counts = design_counts(path, COND)
    assert counts["assessed"] == 4
    assert counts["hard_qc_passed"] == 2


def test_switching_mode_does_not_reuse_the_other_mode_s_verdicts(tmp_path):
    """The mode decides admission, so it belongs in the policy digest."""
    _record(tmp_path, "all_qc")
    path = _record(tmp_path, "post_gini")

    assert _eligible(path) == sorted(AFTER_GINI)


def test_the_resolved_policy_carries_the_retention_mode():
    """Guard the guard: without it the two modes share a digest."""
    from types import SimpleNamespace

    strict = resolved_qc_policy(SimpleNamespace(candidate_retention="post_gini"))
    loose = resolved_qc_policy(SimpleNamespace(candidate_retention="all_qc"))

    assert strict["candidate_retention"] == "post_gini"
    assert strict != loose


def test_an_unknown_mode_does_not_quietly_admit_everything(tmp_path):
    """A typo in the config must not widen the eligible set in silence."""
    import pytest

    with pytest.raises(ValueError, match="candidate_retention"):
        record_stage2_inventory(
            tmp_path,
            condition_id=COND,
            cleared_hard_gates=_frame(CLEARED),
            after_gini=_frame(AFTER_GINI),
            shortlisted=_frame([A]),
            indexed=list(AFTER_GINI),
            retention="whatever",
        )
