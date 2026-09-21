"""A failed design calculation must fail the run, not produce a plausible panel.

Task 1 of the 2026-09-21 valid-design plan. Three separable claims:

1. Output eligibility is a function of run state AND panel qualification. A
   panel that qualified on a run that did not finish is not recommendable, and
   neither is an unqualified panel on a run that did.
2. A required thermodynamic calculation that fails raises a `DesignError`
   rather than returning a substitute value. A numerical failure is never
   recorded as a candidate QC rejection, because the two mean opposite things:
   one says the candidate is bad, the other says we do not know.
3. The CLI turns a `DesignError` into a structured failure record and a
   nonzero exit, and no recommended FASTA is written for it.
"""

import json
import subprocess
import sys

import pytest

from neoswga.core.design_result import (
    RunState,
    TerminationReason,
    describe_failure,
    recommendation_allowed,
)
from neoswga.core.exceptions import (
    DesignError,
    InvalidDesignRequest,
    ModelEvaluationError,
    NeoSWGAError,
    ReferenceDataError,
    UnsupportedModelError,
)

# ---------------------------------------------------------------------------
# 1. Output eligibility
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "state,qualified,expected",
    [
        ("finished", True, True),
        ("finished", False, False),
        ("failed", True, False),
        ("interrupted", True, False),
    ],
)
def test_recommendation_requires_finished_qualified_run(state, qualified, expected):
    assert recommendation_allowed(state, qualified) is expected


def test_run_state_is_not_overloaded_with_qualification():
    """`finished` describes the run; it says nothing about the panel.

    The states the plan forbids collapsing: a run can finish having found
    nothing acceptable, and that is not the same event as a run that crashed.
    """
    assert RunState.FINISHED != RunState.FAILED
    assert recommendation_allowed(RunState.FINISHED, qualified=False) is False
    assert recommendation_allowed(RunState.FINISHED, qualified=True) is True


def test_budget_exhaustion_is_a_termination_reason_not_a_failed_run():
    """A search that spent its allowance still finished. See contract item 11."""
    assert TerminationReason.BUDGET_EXHAUSTED != TerminationReason.ERROR
    assert recommendation_allowed(RunState.FINISHED, qualified=True) is True


def test_unknown_run_state_is_refused_rather_than_treated_as_finished():
    with pytest.raises(ValueError, match="run state"):
        recommendation_allowed("done", True)


# ---------------------------------------------------------------------------
# 2. The error hierarchy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "error",
    [InvalidDesignRequest, ReferenceDataError, UnsupportedModelError, ModelEvaluationError],
)
def test_design_errors_are_catchable_as_one_family(error):
    assert issubclass(error, DesignError)
    assert issubclass(error, NeoSWGAError)


def test_design_error_is_dependency_light():
    """`exceptions.py` is imported at CLI startup; it must not reach sklearn.

    Known Issue 12 is the reason this module has a dependency budget at all.
    """
    import ast
    import pathlib

    source = pathlib.Path(
        __import__("neoswga.core.exceptions", fromlist=["x"]).__file__
    ).read_text()
    imported = {
        node.module.split(".")[0]
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.ImportFrom) and node.module
    }
    imported |= {
        alias.name.split(".")[0]
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    assert imported <= {"dataclasses", "typing", "enum"}, imported


def test_a_numerical_failure_is_not_a_qc_rejection():
    """`ModelEvaluationError` must not be reachable through a QC reason code.

    A candidate rejected by QC is a measurement. A candidate whose Tm could not
    be computed is an absence of one. Recording the second as the first is how
    an unmeasured primer becomes an apparently screened one.
    """
    error = ModelEvaluationError("tm", "ACGTACGTACGT", "no nearest-neighbour parameters")
    assert error.quantity == "tm"
    assert not getattr(error, "qc_reason", None)


def test_design_failures_preserve_their_cause():
    try:
        try:
            raise ZeroDivisionError("float division")
        except ZeroDivisionError as exc:
            raise ModelEvaluationError("tm", "ACGT", "non-finite") from exc
    except ModelEvaluationError as exc:
        assert isinstance(exc.__cause__, ZeroDivisionError)


def test_describe_failure_names_stage_input_and_model():
    record = describe_failure(
        ModelEvaluationError("tm", "ACGTACGTACGT", "non-finite result"),
        stage="panel_evaluation",
        request_hash="abc123",
    )
    assert record["run_state"] == RunState.FAILED
    assert record["termination"] == TerminationReason.ERROR
    assert record["stage"] == "panel_evaluation"
    assert record["request_hash"] == "abc123"
    assert record["error_type"] == "ModelEvaluationError"
    assert "tm" in record["quantity"]
    assert record["recommendation_written"] is False
    # Must survive a JSON round trip: this is written to a failure artifact.
    assert json.loads(json.dumps(record)) == record


# ---------------------------------------------------------------------------
# 3. A thermodynamic failure on a production design call
# ---------------------------------------------------------------------------


def test_tm_batch_raises_rather_than_returning_nan(monkeypatch):
    """`calculate_tm_batch` used to write NaN for an unexpected failure.

    NaN propagates into a Tm window comparison as False, so the candidate is
    dropped and the run reports a smaller pool. That is indistinguishable
    downstream from a candidate that genuinely missed the window.
    """
    from neoswga.core import thermodynamics

    def explode(seq, *args, **kwargs):
        raise RuntimeError("nearest-neighbour table unavailable")

    monkeypatch.setattr(thermodynamics, "calculate_tm_with_salt", explode)
    with pytest.raises(ModelEvaluationError) as caught:
        thermodynamics.calculate_tm_batch(["ACGTACGTACGT"])
    assert "ACGTACGTACGT" in str(caught.value)
    assert isinstance(caught.value.__cause__, RuntimeError)


def test_invalid_bases_are_a_named_rejection_not_a_numerical_failure():
    """A user sequence with a bad base is a QC rejection, and stays one."""
    from neoswga.core import thermodynamics

    with pytest.raises(thermodynamics.InvalidSequenceError) as caught:
        thermodynamics.calculate_tm_batch(["ACGTXXGTACGT"])
    assert caught.value.qc_reason == "invalid_base"


def test_heterodimer_failure_does_not_become_zero_free_energy(monkeypatch):
    """A dimer check that fails must not report a pair as dimer-free.

    Zero free energy is the most permissive answer available, so swallowing
    the error admits exactly the pairs the screen exists to reject.
    """
    from neoswga.core import thermodynamic_filter

    def explode(*args, **kwargs):
        raise RuntimeError("structure model failed")

    monkeypatch.setattr(thermodynamic_filter, "check_heterodimer", explode)
    with pytest.raises(ModelEvaluationError):
        thermodynamic_filter._check_heterodimer_pair(
            ("ACGTACGTACGT", "TGCATGCATGCA", 0, 1, {"temp": 30.0})
        )


def test_coverage_reach_failure_does_not_fall_back_to_a_default(monkeypatch):
    """A default reach silently rescales every coverage figure in the report."""
    from neoswga.core import coverage

    monkeypatch.setattr(coverage, "_reach_from_registry", None, raising=False)
    with pytest.raises(UnsupportedModelError):
        coverage.polymerase_extension_reach("not-a-polymerase")


def test_record_starts_getter_failure_propagates(monkeypatch):
    """Losing record boundaries lets a coverage window cross a record edge."""
    from neoswga.core import coverage

    class BrokenCache:
        def get_record_starts(self, prefix):
            raise OSError("index truncated")

    with pytest.raises(ReferenceDataError):
        coverage._record_starts_for(BrokenCache(), "fg")


def test_a_cache_without_record_starts_is_still_allowed():
    """Absence of the getter is a different thing from the getter failing."""
    from neoswga.core import coverage

    assert coverage._record_starts_for(object(), "fg") is None


def test_discrimination_profile_reports_failures_rather_than_skipping(monkeypatch):
    """An aggregate that drops its failures reports a mean over an unknown set."""
    from neoswga.core import occupancy

    def explode(seq):
        raise RuntimeError("no parameters")

    monkeypatch.setattr(
        "neoswga.core.thermodynamics.calculate_enthalpy_entropy", explode, raising=True
    )

    class Conditions:
        temp = 30.0

        def calculate_effective_tm(self, seq):
            return 35.0

    with pytest.raises(ModelEvaluationError):
        occupancy.discrimination_profile(["ACGTACGTACGT"], Conditions())


# ---------------------------------------------------------------------------
# 4. The CLI boundary writes a failure and no recommendation
# ---------------------------------------------------------------------------


def _run_cli(args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True,
        text=True,
        cwd=str(cwd),
    )


def test_a_design_error_exits_nonzero_and_writes_no_fasta(tmp_path):
    """An unresolvable request fails the command rather than designing anyway."""
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"coverage_reach": 0}))
    result = _run_cli(["optimize", "-j", str(params), "--data-dir", str(tmp_path)], tmp_path)
    assert result.returncode != 0
    assert not list(tmp_path.glob("*.fasta"))


def test_the_failure_record_lands_beside_the_run(tmp_path):
    """A stale output directory must not look like a current result.

    A directory holding last week's `step4_improved_df.csv` reads exactly like
    one holding this morning's. This file is what says otherwise, so it has to
    be written where the run's outputs are -- which several commands resolve
    into the `parameter` module rather than carrying on the argparse namespace.
    """
    import argparse

    from neoswga.cli._failure import report_design_failure

    args = argparse.Namespace(command="optimize", data_dir=str(tmp_path))
    report_design_failure(args, ReferenceDataError("positions", "digest mismatch", "re-run filter"))

    record = json.loads((tmp_path / "design_failure.json").read_text())
    assert record["run_state"] == RunState.FAILED
    assert record["stage"] == "optimize"
    assert record["artifact"] == "positions"
    assert record["remediation"] == "re-run filter"
    assert record["recommendation_written"] is False
    assert recommendation_allowed(record["run_state"], qualified=True) is False


def test_an_unexpected_exception_is_recorded_as_unexpected(tmp_path):
    """`expected: false` is what separates a contract refusal from a defect."""
    import argparse

    from neoswga.cli._failure import report_design_failure

    args = argparse.Namespace(command="filter", data_dir=str(tmp_path))
    report_design_failure(args, RuntimeError("something broke"), expected=False)

    record = json.loads((tmp_path / "design_failure.json").read_text())
    assert record["expected"] is False
    assert record["error_type"] == "RuntimeError"
    assert record["traceback"] is not None


def test_a_failure_record_that_cannot_be_written_does_not_mask_the_failure(tmp_path):
    """Losing the record is worse than the original error only if it hides it."""
    import argparse

    from neoswga.cli._failure import write_failure_artifact

    unwritable = tmp_path / "not-a-directory"
    unwritable.write_text("")
    args = argparse.Namespace(command="optimize", data_dir=str(unwritable))

    assert write_failure_artifact(args, ReferenceDataError("positions", "gone")) is None


def test_export_refuses_a_run_that_did_not_finish(tmp_path):
    """Task 1's gate on the recommendation export."""
    from neoswga.core.design_result import recommendation_allowed

    failure = json.loads(
        json.dumps(
            describe_failure(
                ReferenceDataError("positions", "digest mismatch"),
                stage="reference_check",
                request_hash="deadbeef",
            )
        )
    )
    assert recommendation_allowed(failure["run_state"], qualified=True) is False
