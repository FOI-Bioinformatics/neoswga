"""Under `all_qc`, an unmeasurable Gini must not lose the inventory.

Phase 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, audit
finding F5, part A.

`check_gini_stage_kept_something` raises when hard-QC survivors exist and none
of them has a measurable Gini. That guard is right for the CSV: `step2_df.csv`
is the working shortlist, and writing an empty one leaves the design with
nothing to optimise over and no explanation.

It ran before the inventory was recorded, and regardless of retention mode. So
under `all_qc` -- whose entire purpose is to keep every hard-QC survivor
addressable -- a run where evenness happened to be unmeasurable aborted before a
single inventory row was written. The candidates had cleared every declared
requirement. What they had not done was bind often enough for their spacing to
be measured, which `min_gini_sites` exists to say is not a judgement about them.

This is most likely on exactly the targets the retention work was for: a small
genome where most candidates bind once or twice. On the bundled plasmid example
10,158 of 10,532 indexed k-mers bind exactly once.

So the inventory is written first and the guard still raises afterwards. The
run still fails, with the same message; it just no longer throws away the record
of what it enumerated on the way out.
"""

import pandas as pd
import pytest

from neoswga.core.candidate_inventory import STAGE2_INVENTORY_NAME, CandidateInventory

A, C, G = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"


def test_the_guard_still_refuses_an_empty_evenness_stage():
    """The behaviour that must not change: the run fails, and says why."""
    from neoswga.core.filter import check_gini_stage_kept_something

    with pytest.raises(Exception) as excinfo:
        check_gini_stage_kept_something(
            pd.DataFrame({"primer": [A, C]}), pd.DataFrame({"primer": []})
        )

    assert "min_gini_sites" in str(excinfo.value)


def test_the_guard_is_quiet_when_something_survived():
    from neoswga.core.filter import check_gini_stage_kept_something

    check_gini_stage_kept_something(pd.DataFrame({"primer": [A, C]}), pd.DataFrame({"primer": [A]}))


def test_the_inventory_is_recorded_before_the_evenness_guard():
    """The ordering, as a property of the source rather than of one run.

    Reproducing this end to end needs a genome whose candidates all bind once,
    which is a slow fixture for a question about statement order. The AST says
    it directly: the inventory write must come before the evenness guard.
    """
    import ast
    import inspect
    import textwrap

    from neoswga.core import stage2_recording

    tree = ast.parse(textwrap.dedent(inspect.getsource(stage2_recording._index_and_record)))
    calls = [
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    ]

    assert "_record_run_inventory" in calls, "the inventory is no longer recorded here"
    assert "check_gini_stage_kept_something" in calls, "the evenness guard has gone"
    assert calls.index("_record_run_inventory") < calls.index("check_gini_stage_kept_something"), (
        "the evenness guard runs before the inventory is written, so a run where "
        "no candidate has a measurable Gini loses the record of everything it "
        "enumerated"
    )


def test_step_two_still_reaches_that_ordering():
    """Guard the guard: the block above must still be on the filter's path."""
    import ast
    import inspect
    import textwrap

    from neoswga.core import pipeline

    tree = ast.parse(textwrap.dedent(inspect.getsource(pipeline.step2)))
    calls = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "_index_and_record" in calls


def test_an_inventory_written_before_the_guard_holds_the_survivors(tmp_path):
    """What the ordering buys: the record survives the failure.

    Recorded directly rather than through `step2`, because the point is what the
    inventory holds once written, not how the filter reached it.
    """
    from neoswga.core.candidate_inventory import record_stage2_inventory

    cleared = pd.DataFrame({"primer": [A, C, G], "fg_count": [1, 1, 1], "bg_count": [0, 0, 0]})
    empty = cleared.iloc[0:0]

    record_stage2_inventory(
        tmp_path,
        condition_id="cond",
        cleared_hard_gates=cleared,
        after_gini=empty,
        shortlisted=empty,
        indexed=list(cleared["primer"]),
    )

    with CandidateInventory(tmp_path / STAGE2_INVENTORY_NAME) as inventory:
        assert sorted(inventory.iter_eligible("cond", [12])) == sorted([A, C, G])
        assert all(not inventory.metrics(s)["passed_gini"] for s in (A, C, G))
