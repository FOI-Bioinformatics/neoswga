"""Six numbers, named, instead of one labelled as if it were the pool.

Task 3 of the condition-aware pool design plan, final steps.

The funnel reported a handful of counts and none of them answered "how many
candidates could the optimizer actually have chosen from". The plan names the
six distinctions that matter, and this pins them:

counted
    Enumerated at all.
assessed
    Judged against the declared hard gates under some reaction.
hard_qc_passed
    Cleared those gates. Eligible, whatever the ranking said.
shortlisted
    Written to `step2_df.csv`, the working pool `max_primer` chose.
indexed
    Given a background position index, so a specificity number computed for it
    is a measurement rather than a zero.
examined
    Carried into stage 3 and therefore reachable by the optimizer.

Stage 3's opt-in efficacy filter is a policy, not a gate. When it drops a
candidate the reason is recorded against that candidate rather than left as an
absent row, so a later reader can tell "the model scored it low" from "it never
cleared QC".
"""

import pytest

from neoswga.core.candidate_inventory import (
    CandidateInventory,
    design_counts,
    record_stage2_inventory,
    record_stage3_policy,
)

A, C, G = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"
COND = "tm-2026-09-14:abc123"


@pytest.fixture
def populated(tmp_path):
    import pandas as pd

    cleared = pd.DataFrame({"primer": [A, C, G], "fg_count": [9, 5, 1]})
    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared[cleared["primer"] != G],
        shortlisted=cleared[cleared["primer"] == A],
        indexed=[A, C, G],
    )
    return path


def test_the_counts_are_reported_by_name(populated):
    """`counted` and `examined` are unknown here, and say so.

    This fixture records no enumerated total and runs no search, so both are
    `None`. They used to report 3 and 0: the first was the post-QC survivors
    under a name promising the universe they came from, and the second was
    stage-3 membership under a name promising search evaluations. A number that
    cannot be told from a measurement is worse than an absent one.
    """
    counts = design_counts(populated, COND)

    assert counts["counted"] is None
    assert counts["assessed"] == 3
    assert counts["hard_qc_passed"] == 3
    assert counts["shortlisted"] == 1
    assert counts["indexed"] == 3
    assert counts["carried_to_stage3"] == 0, "nothing has reached stage 3 yet"
    assert counts["examined"] is None, "no search reports its evaluations yet"


def test_a_partly_indexed_inventory_still_records_the_hard_qc_verdict(tmp_path):
    """Indexing is a separate fact from passing the gates.

    Named for the shortlist-only `legacy` retention mode when that existed. The
    mode has gone; the invariant it exercised has not, because a caller can
    still pass a narrower `indexed` list than the candidates it cleared.
    """
    import pandas as pd

    cleared = pd.DataFrame({"primer": [A, C, G]})
    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared,
        after_gini=cleared,
        shortlisted=cleared[cleared["primer"] == A],
        indexed=[A],
    )

    counts = design_counts(path, COND)
    assert counts["indexed"] == 1
    assert (
        counts["hard_qc_passed"] == 3
    ), "not indexing a candidate does not make it fail the hard gates"


def test_stage_three_records_what_it_carried(populated):
    """Under its own name. What stage 3 kept is not what a search evaluated."""
    record_stage3_policy(populated, COND, carried=[A, C], rejected={}, policy="carry_forward")

    assert design_counts(populated, COND)["carried_to_stage3"] == 2


def test_an_efficacy_rejection_records_its_reason(populated):
    record_stage3_policy(
        populated,
        COND,
        carried=[A],
        rejected={C: "amp_pred 4.1 below min_amp_pred 10.0"},
        policy="amp_model",
    )

    with CandidateInventory(populated) as inv:
        reasons = inv.reasons(C, COND, policy_version="stage3:amp_model")
        assert reasons == ["amp_pred 4.1 below min_amp_pred 10.0"]
        # The hard-QC verdict is a separate record and is untouched.
        assert list(inv.iter_eligible(COND, [12])) == sorted([A, C, G])


def test_a_stage_three_rejection_is_not_a_hard_qc_failure(populated):
    """The distinction the funnel could not express.

    An opt-in model scoring a candidate low is a policy choice. Reading it as a
    QC failure is how a ranking came to look like a requirement.
    """
    record_stage3_policy(
        populated, COND, carried=[A], rejected={C: "low score"}, policy="amp_model"
    )

    counts = design_counts(populated, COND)
    assert counts["hard_qc_passed"] == 3
    assert counts["carried_to_stage3"] == 1
