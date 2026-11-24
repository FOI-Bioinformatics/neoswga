"""Every pipeline-step failure leaves a record that names the step.

Before 2026-11-23 only a `DesignError` reaching the command boundary wrote
`design_failure.json`. `design` caught that error itself and exited, and a step
that failed any other way (`exit_on_step_failure`, a refused prerequisite, a
`ValueError` from the filter) wrote nothing. The previous run's
`step4_improved_df.csv` then exported as current
(2026-11-23 sweep, tier 1 items 1 and 4).

The record now carries the step, a downstream step refuses while an upstream
step's record is present, and a step's success clears only its own record.
"""

import json
import sys

import pytest

from neoswga import cli_unified
from neoswga.cli import _failure
from neoswga.core.exceptions import UnsupportedModelError


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    results = tmp_path / "results"
    results.mkdir()
    (tmp_path / "params.json").write_text(json.dumps({"data_dir": "results"}))
    monkeypatch.chdir(tmp_path)
    # The boundary reads the run directory from `parameter.data_dir` first.
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", None, raising=False)
    return results


def _main(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["neoswga", *argv])
    with pytest.raises(SystemExit) as info:
        cli_unified.main()
    return info.value.code


def _main_ok(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["neoswga", *argv])
    cli_unified.main()


def _record(results):
    path = results / "design_failure.json"
    return json.loads(path.read_text()) if path.exists() else None


def _write_record(results, stage):
    (results / "design_failure.json").write_text(
        json.dumps({"run_state": "failed", "stage": stage, "qualified": False})
    )


def test_design_failure_writes_a_record_naming_the_step(workspace, monkeypatch):
    from neoswga.cli import commands

    monkeypatch.setattr(commands, "check_kmer_counter_available", lambda *a, **k: None)
    for name in ("run_step1", "run_step2", "run_step3"):
        monkeypatch.setattr(commands, name, lambda args: None)

    def _refuse(args):
        raise UnsupportedModelError("oligo length for equiphi29", [6])

    monkeypatch.setattr(commands, "run_step4", _refuse)
    assert _main(monkeypatch, "design", "-j", "params.json") == 1
    record = _record(workspace)
    assert record is not None
    assert record["stage"] == "optimize"


def test_a_step_that_exits_writes_a_record(workspace, monkeypatch):
    def _exits(args):
        try:
            raise ValueError("Gini stage removed all 95 candidates")
        except ValueError as error:
            _failure.exit_on_step_failure("Step 2", error, _failure.logger)

    monkeypatch.setattr(cli_unified, "run_step2", _exits)
    assert _main(monkeypatch, "filter", "-j", "params.json") == 1
    record = _record(workspace)
    assert record["stage"] == "filter"
    assert "Gini stage" in record["message"]


def test_a_value_error_at_the_boundary_records_a_pipeline_command(workspace, monkeypatch):
    def _raises(args):
        raise ValueError("bad threshold")

    monkeypatch.setattr(cli_unified, "run_step2", _raises)
    assert _main(monkeypatch, "filter", "-j", "params.json") == 1
    assert _record(workspace)["stage"] == "filter"


def test_a_non_pipeline_command_writes_no_record(workspace, monkeypatch):
    def _raises(args):
        raise ValueError("not a primer")

    monkeypatch.setattr(cli_unified, "analyze_primer_set", _raises)
    _main(monkeypatch, "analyze-set", "--primers", "ACGT", "-j", "params.json")
    assert _record(workspace) is None


@pytest.mark.parametrize(
    "command,handler", [("optimize", "run_step4"), ("prepare-candidates", "run_step3")]
)
def test_a_downstream_step_refuses_after_a_failed_filter(workspace, monkeypatch, command, handler):
    _write_record(workspace, "filter")
    ran = []
    monkeypatch.setattr(cli_unified, handler, lambda args: ran.append(1))
    assert _main(monkeypatch, command, "-j", "params.json") == 1
    assert ran == []
    assert _record(workspace)["stage"] == "filter", "the upstream record was overwritten"


def test_optimize_runs_after_its_own_failure(workspace, monkeypatch):
    _write_record(workspace, "optimize")
    ran = []
    monkeypatch.setattr(cli_unified, "run_step4", lambda args: ran.append(1))
    _main_ok(monkeypatch, "optimize", "-j", "params.json")
    assert ran == [1]


def test_filter_success_clears_only_its_own_stage(workspace, monkeypatch):
    monkeypatch.setattr(cli_unified, "run_step2", lambda args: None)
    _write_record(workspace, "filter")
    _main_ok(monkeypatch, "filter", "-j", "params.json")
    assert _record(workspace) is None

    _write_record(workspace, "optimize")
    _main_ok(monkeypatch, "filter", "-j", "params.json")
    assert _record(workspace)["stage"] == "optimize"


@pytest.mark.parametrize(
    "stage,step,blocks",
    [
        ("filter", "optimize", True),
        ("count-kmers", "prepare-candidates", True),
        ("design", "optimize", True),
        ("optimize", "optimize", False),
        ("prepare-candidates", "prepare-candidates", False),
        ("simulate", "optimize", False),
    ],
)
def test_upstream_failure(tmp_path, stage, step, blocks):
    _write_record(tmp_path, stage)
    assert (_failure.upstream_failure(str(tmp_path), step) is not None) is blocks


def test_an_unreadable_record_blocks_downstream(tmp_path):
    (tmp_path / "design_failure.json").write_text("{ truncated")
    assert _failure.upstream_failure(str(tmp_path), "optimize") is not None


def test_design_from_optimize_refuses_after_a_failed_filter(workspace, monkeypatch):
    from neoswga.cli import commands

    monkeypatch.setattr(commands, "check_kmer_counter_available", lambda *a, **k: None)
    ran = []
    monkeypatch.setattr(commands, "run_step4", lambda args: ran.append(1))
    _write_record(workspace, "filter")
    assert _main(monkeypatch, "design", "-j", "params.json", "--start-from", "4") == 1
    assert ran == []
    assert _record(workspace)["stage"] == "filter", "the upstream record was overwritten"


# --- Final-review fixes (2026-11-23): records nothing can clear, and the
# --- directory the checks read against the directory the run uses.


def test_design_that_runs_no_step_writes_no_record(workspace, monkeypatch):
    from neoswga.cli import commands

    monkeypatch.setattr(commands, "check_kmer_counter_available", lambda *a, **k: None)
    assert (
        _main(monkeypatch, "design", "-j", "params.json", "--start-from", "4", "--stop-at", "3")
        == 1
    )
    assert _record(workspace) is None


def test_design_without_a_counter_writes_no_record(workspace, monkeypatch):
    from neoswga.cli import commands

    def _missing(*_a, **_k):
        sys.exit(1)

    monkeypatch.setattr(commands, "check_kmer_counter_available", _missing)
    assert _main(monkeypatch, "design", "-j", "params.json") == 1
    assert _record(workspace) is None


@pytest.mark.parametrize("content", ['{"run_state": "failed", "stage": "design"}', "not json"])
def test_count_kmers_success_clears_a_record_no_step_owns(workspace, monkeypatch, content):
    (workspace / "design_failure.json").write_text(content)
    monkeypatch.setattr(cli_unified, "run_step1", lambda args: None)
    _main_ok(monkeypatch, "count-kmers", "-j", "params.json")
    assert not (workspace / "design_failure.json").exists()


def test_a_design_rerun_recovers_from_a_legacy_design_record(workspace, monkeypatch):
    from neoswga.cli import commands

    monkeypatch.setattr(commands, "check_kmer_counter_available", lambda *a, **k: None)
    ran = []
    for name in ("run_step1", "run_step2", "run_step3", "run_step4"):
        monkeypatch.setattr(commands, name, lambda args, n=name: ran.append(n))
    _write_record(workspace, "design")
    _main_ok(monkeypatch, "design", "-j", "params.json")
    assert ran == ["run_step1", "run_step2", "run_step3", "run_step4"]


def test_checks_read_the_directory_the_run_uses(tmp_path, monkeypatch):
    """data_dir is resolved against the working directory, as `get_params` does,
    not against the params file's directory."""
    (tmp_path / "cfg").mkdir()
    (tmp_path / "cfg" / "params.json").write_text(json.dumps({"data_dir": "results"}))
    (tmp_path / "results").mkdir()
    _write_record(tmp_path / "results", "filter")
    monkeypatch.chdir(tmp_path)
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", None, raising=False)
    ran = []
    monkeypatch.setattr(cli_unified, "run_step4", lambda args: ran.append(1))
    assert _main(monkeypatch, "optimize", "-j", "cfg/params.json") == 1
    assert ran == []


def test_a_record_goes_to_the_run_directory_not_data_dir_flag(workspace, monkeypatch):
    """`-z` is not read by `optimize`, so since 2026-11-24 it refuses at
    parsing: nothing runs and no record is written under the flag's path."""

    def _exits(args):
        sys.exit(1)

    monkeypatch.setattr(cli_unified, "run_step4", _exits)
    assert _main(monkeypatch, "optimize", "-j", "params.json", "-z", "elsewhere") == 2
    assert not (workspace.parent / "elsewhere" / "design_failure.json").exists()
