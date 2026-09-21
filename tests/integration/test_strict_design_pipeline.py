"""The whole pipeline, and what each refusal does to the output directory.

Task 10 of the 2026-09-21 valid-design plan: an end-to-end run over a small
reference, "plus error injections at each required step. The failure cases
must return nonzero status and no qualifying recommendation artifact."

The individual contracts are unit-tested elsewhere. What this file adds is
that they hold through a real process: a refusal has to survive argument
parsing, parameter resolution, the step's own try/except, and the command
boundary, and each of those has swallowed one before. Known Issue 8's class
is precisely a check that exists and is not reached.

Every injection below is a configuration a user could plausibly write, not a
monkeypatched internal. What is asserted of each is the same three things:

1. the command exits nonzero;
2. it leaves `design_failure.json` saying what failed;
3. no recommendation can be exported from the directory afterwards.

The third matters most. The first two are about this run; the third is about
what the directory looks like to the next command, which is the only thing a
later reader sees.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("jellyfish") is None,
    reason="the pipeline's first step shells out to jellyfish",
)

PACKAGED = Path(__file__).resolve().parent.parent.parent / "neoswga" / "core" / "smoke"


def run(args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True,
        text=True,
        cwd=str(cwd),
        timeout=600,
    )


def write_params(directory, **overrides):
    params = {
        "fg_genomes": ["pcDNA.fasta"],
        "bg_genomes": [],
        "fg_prefixes": ["pcDNA"],
        "bg_prefixes": [],
        "data_dir": "results",
        "min_k": 10,
        "max_k": 10,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_fg_freq": 1e-6,
        "max_bg_freq": 1.0,
        "max_gini": 1.0,
        "min_gini_sites": 1,
        "max_primer": 60,
        "min_tm": 0,
        "max_tm": 100,
        "gc_min": 0.0,
        "gc_max": 1.0,
        "num_primers": 4,
        "target_set_size": 4,
        "max_sets": 1,
        "iterations": 2,
        "cpus": 1,
        "schema_version": 2,
    }
    params.update(overrides)
    path = directory / "params.json"
    path.write_text(json.dumps(params))
    return path


def stage_workspace(tmp_path):
    workspace = tmp_path / "work"
    workspace.mkdir()
    shutil.copy(PACKAGED / "pcDNA.fasta", workspace / "pcDNA.fasta")
    return workspace


def prepared(tmp_path, **overrides):
    """A workspace taken through count, filter and candidate preparation."""
    workspace = stage_workspace(tmp_path)
    write_params(workspace, **overrides)
    for step in ("count-kmers", "filter", "prepare-candidates"):
        result = run([step, "-j", "params.json"], workspace)
        assert result.returncode == 0, f"{step} failed:\n{result.stdout}\n{result.stderr}"
    return workspace


def failure_record(workspace):
    path = workspace / "results" / "design_failure.json"
    if not path.exists():
        return None
    return json.loads(path.read_text())


def assert_nothing_recommendable(workspace, expect_record=True):
    """The three things every refusal owes the next command."""
    results = workspace / "results"

    if expect_record:
        record = failure_record(workspace)
        assert record is not None, "a refused run left no failure record"
        assert record["run_state"] in {"failed", "interrupted"}
        assert record["recommendation_written"] is False

    export = run(["export", "-d", "results", "-o", "order", "--format", "fasta"], workspace)
    assert export.returncode != 0, (
        "a directory whose last run was refused still exported a recommendation:\n" + export.stdout
    )
    orders = workspace / "order"
    assert not (orders.exists() and list(orders.glob("*.fasta")))


# ---------------------------------------------------------------------------
# The happy path
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def completed(tmp_path_factory):
    workspace = prepared(tmp_path_factory.mktemp("strict"))
    result = run(["optimize", "-j", "params.json", "--seed", "1"], workspace)
    assert result.returncode == 0, result.stdout + result.stderr
    return workspace


def test_the_four_steps_produce_the_expected_artifacts(completed):
    results = completed / "results"

    assert (results / "step2_df.csv").exists()
    assert (results / "step3_df.csv").exists()
    assert (results / "step4_improved_df.csv").exists()
    assert (results / "step4_improved_df_summary.json").exists()


def test_a_finished_run_leaves_no_failure_record(completed):
    """`optimize` clears one, so an earlier failure cannot block this export."""
    assert failure_record(completed) is None


def test_the_finished_run_exports_a_recommendation(completed):
    export = run(["export", "-d", "results", "-o", "order", "--format", "fasta"], completed)

    assert export.returncode == 0, export.stdout + export.stderr
    assert list((completed / "order").glob("*.fasta"))


def test_the_run_manifest_records_the_request_hash(completed):
    manifest = json.loads((completed / "results" / "run_manifest.json").read_text())
    optimize_steps = [s for s in manifest["steps"] if s["step"] == "optimize"]

    assert optimize_steps
    assert optimize_steps[-1]["extra"]["request_hash"]


# ---------------------------------------------------------------------------
# Error injections, one per contract
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "overrides,expected",
    [
        pytest.param({"glycerol_percent": 10.0}, "glycerol", id="additive-with-no-model"),
        pytest.param(
            {"min_selectivity_density": 60.0},
            "min_selectivity_density",
            id="background-limit-with-no-background",
        ),
    ],
)
def test_a_refusal_on_the_design_path_leaves_nothing_recommendable(tmp_path, overrides, expected):
    """Each of these is a params file somebody could write.

    These reach the design path, so all three obligations apply: exit
    nonzero, leave a record saying what failed, and leave nothing the next
    command will export.

    The record is the one that had to be fixed here. Both refusals fire inside
    `resolve_design_request`, which runs BEFORE `get_params` populates the
    `parameter` module, so the failure writer looked for an output directory
    the module did not yet know about and wrote nowhere. The run exited
    nonzero and the directory still looked like its previous success. Same
    ordering trap as `warn_on_condition_drift`, reached from the other side.
    """
    workspace = prepared(tmp_path)
    write_params(workspace, **overrides)

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, f"accepted {overrides}:\n{output}"
    assert expected in output, output[-1500:]
    assert_nothing_recommendable(workspace)


@pytest.mark.parametrize(
    "overrides,expected",
    [
        pytest.param({"coverage_reach": 0}, "minimum of 1", id="explicit-zero-reach"),
        pytest.param({"polymerase": "taq"}, "not one of", id="unsupported-polymerase"),
    ],
)
def test_a_refusal_at_schema_validation_never_enters_the_design_path(tmp_path, overrides, expected):
    """Caught earlier, and better, than the design-request contract catches it.

    The schema refuses these before any handler runs, naming the permitted
    values. No failure record is written and none is wanted: the design path
    was never entered, so there is no stage to name, and the directory has no
    step-4 output for anyone to mistake for current.

    Asserted separately from the design-path refusals rather than folded in,
    because "exits nonzero" is the same observation from two different
    mechanisms and only one of them owes a record.
    """
    workspace = prepared(tmp_path)
    write_params(workspace, **overrides)

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, f"accepted {overrides}:\n{output}"
    assert expected in output, output[-1500:]
    assert failure_record(workspace) is None
    assert_nothing_recommendable(workspace, expect_record=False)


def test_a_missing_position_index_is_refused(tmp_path):
    """The index is deleted after preparation, so the pool is real and unscored."""
    workspace = prepared(tmp_path)
    for path in workspace.glob("*_positions.h5"):
        path.unlink()
    for path in (workspace / "results").glob("*_positions.h5"):
        path.unlink()

    result = run(["optimize", "-j", "params.json"], workspace)

    assert result.returncode != 0, result.stdout


def test_a_stale_index_from_another_reference_is_refused(tmp_path):
    """The structurally perfect failure: every dataset present, wrong genome.

    The FASTA is replaced after the index was built, so the index is complete,
    internally consistent, and describes a sequence this run is not targeting.
    """
    workspace = prepared(tmp_path)
    shutil.copy(PACKAGED / "pLTR.fasta", workspace / "pcDNA.fasta")

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, output
    assert "reference" in output.lower() or "count-kmers" in output
    assert_nothing_recommendable(workspace)


def test_a_reaction_the_filter_never_recorded_is_refused(tmp_path):
    """Candidates selected under one chemistry, scored under another.

    This used to warn and proceed, searching the shortlist while everything
    the inventory held went unreachable.
    """
    workspace = stage_workspace(tmp_path)
    write_params(workspace)
    assert run(["count-kmers", "-j", "params.json"], workspace).returncode == 0
    assert (
        run(["filter", "-j", "params.json", "--preset", "enhanced_equiphi29"], workspace).returncode
        == 0
    )
    assert run(["prepare-candidates", "-j", "params.json"], workspace).returncode == 0

    result = run(["optimize", "-j", "params.json"], workspace)
    output = result.stdout + result.stderr

    assert result.returncode != 0, output
    assert "filter" in output
    assert_nothing_recommendable(workspace)


def test_the_retired_command_name_is_refused(tmp_path):
    workspace = stage_workspace(tmp_path)
    write_params(workspace)

    result = run(["score", "-j", "params.json"], workspace)

    assert result.returncode != 0
    assert "prepare-candidates" in result.stdout + result.stderr
