"""The manifest must be able to answer "what produced this file".

Before 2026-09-10 an entry carried no timing of any kind, mixed the step's
output into `input_checksums`, and recorded `resolved_params` verbatim from
params.json -- so a run invoked with `-n 160` was recorded as num_primers 24.
`runs/gc_tiers/mid_ecoli/run_manifest.json` held three optimize entries under
two git SHAs, all hash-matching the single step4_improved_df.csv beside them.
"""

import json

import pytest

from neoswga.core import run_manifest as rm


@pytest.fixture
def workspace(tmp_path):
    (tmp_path / "in.csv").write_text("primer\nACGTACGTACGT\n")
    (tmp_path / "out.csv").write_text("primer,set_index\nACGTACGTACGT,0\n")
    return tmp_path


def _entries(data_dir):
    return json.loads((data_dir / rm.MANIFEST_FILENAME).read_text())["steps"]


def test_output_files_are_checksummed_separately_from_inputs(workspace):
    rm.write_manifest(
        step="optimize",
        data_dir=str(workspace),
        input_files=[str(workspace / "in.csv")],
        output_files=[str(workspace / "out.csv")],
    )

    entry = _entries(workspace)[0]

    assert list(entry["input_checksums"]) == [str(workspace / "in.csv")]
    assert list(entry["output_checksums"]) == [str(workspace / "out.csv")]
    assert entry["output_checksums"][str(workspace / "out.csv")] is not None


def test_a_rerun_that_changes_the_output_gets_a_different_output_checksum(workspace):
    rm.write_manifest(
        step="optimize",
        data_dir=str(workspace),
        output_files=[str(workspace / "out.csv")],
    )
    (workspace / "out.csv").write_text("primer,set_index\nTTTTTTTTTTTT,0\n")
    rm.write_manifest(
        step="optimize",
        data_dir=str(workspace),
        output_files=[str(workspace / "out.csv")],
    )

    first, second = _entries(workspace)
    key = str(workspace / "out.csv")
    assert (
        first["output_checksums"][key] != second["output_checksums"][key]
    ), "two optimize entries carry the same output checksum for different files"


def test_missing_output_files_are_skipped_not_recorded_as_none(workspace):
    rm.write_manifest(
        step="filter",
        data_dir=str(workspace),
        output_files=[str(workspace / "never_written.csv")],
    )
    assert _entries(workspace)[0]["output_checksums"] == {}


def test_extra_carries_the_elapsed_time(workspace):
    rm.write_manifest(
        step="filter",
        data_dir=str(workspace),
        extra={"elapsed_seconds": 12.5},
    )
    assert _entries(workspace)[0]["extra"]["elapsed_seconds"] == 12.5


def test_the_cli_wrapper_records_elapsed_time_and_the_effective_set_size(tmp_path):
    """`_record_run_manifest` is what the four step handlers call."""
    from neoswga.cli._common import _record_run_manifest

    class _Parameter:
        data_dir = str(tmp_path)
        num_primers = 160

    class _Args:
        json_file = None
        seed = None

    _record_run_manifest(
        "optimize",
        _Args(),
        _Parameter(),
        extra={"elapsed_seconds": 3.5, "effective_set_size": 160},
    )

    entry = _entries(tmp_path)[0]
    assert entry["extra"]["elapsed_seconds"] == 3.5
    assert entry["extra"]["effective_set_size"] == 160, (
        "resolved_params is a copy of params.json, so `-n 160` is recorded "
        "there as whatever the file said; the effective value must be its own "
        "field"
    )


def test_the_cli_wrapper_passes_outputs_through_as_outputs(tmp_path):
    """A wrapper that accepted output_files and dropped them would pass every
    other test here: the ones that check the separation call write_manifest
    directly."""
    from neoswga.cli._common import _record_run_manifest

    (tmp_path / "step4_improved_df.csv").write_text("primer\nACGTACGTACGT\n")

    class _Parameter:
        data_dir = str(tmp_path)

    class _Args:
        json_file = None
        seed = None

    _record_run_manifest(
        "optimize",
        _Args(),
        _Parameter(),
        output_files=[str(tmp_path / "step4_improved_df.csv")],
    )

    entry = _entries(tmp_path)[0]
    assert list(entry["output_checksums"]) == [str(tmp_path / "step4_improved_df.csv")]
    assert entry["input_checksums"] == {}


def test_every_step_handler_measures_before_it_records():
    """The four call sites had the manifest write above the elapsed
    measurement, so the timing could not have been passed in."""
    from pathlib import Path

    import neoswga.cli.pipeline as cli_pipeline

    lines = Path(cli_pipeline.__file__).read_text().splitlines()
    record_lines = [i for i, ln in enumerate(lines) if "_record_run_manifest(" in ln]
    assert len(record_lines) == 4, record_lines

    for i in record_lines:
        window = "\n".join(lines[max(0, i - 12) : i])
        assert "_elapsed = _time.time() - _t0" in window, (
            f"_record_run_manifest at line {i + 1} is not preceded by the "
            "elapsed measurement; it cannot be recording the step's time"
        )


def test_conditions_can_be_read_for_a_named_step(tmp_path):
    """Two score entries sitting after an optimize entry made a report describe
    the optimize result under the score step's reaction."""
    rm.write_manifest(
        step="optimize",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "phi29", "reaction_temp": 30.0},
    )
    rm.write_manifest(
        step="score",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "bst", "reaction_temp": 63.0},
    )

    assert rm.read_effective_conditions(str(tmp_path), step="optimize") == {
        "polymerase": "phi29",
        "reaction_temp": 30.0,
    }
    assert rm.read_effective_conditions(str(tmp_path), step="score") == {
        "polymerase": "bst",
        "reaction_temp": 63.0,
    }


def test_a_named_step_reads_its_latest_entry_not_its_first(tmp_path):
    """The manifest is append-only, so a step that ran twice has two entries.
    An implementation that scanned forward would return the stale one."""
    rm.write_manifest(
        step="optimize",
        data_dir=str(tmp_path),
        effective_conditions={"reaction_temp": 30.0},
    )
    rm.write_manifest(
        step="optimize",
        data_dir=str(tmp_path),
        effective_conditions={"reaction_temp": 42.0},
    )

    assert rm.read_effective_conditions(str(tmp_path), step="optimize") == {"reaction_temp": 42.0}


def test_reading_without_a_step_keeps_returning_the_latest(tmp_path):
    """export and report call it with no step; that behaviour is unchanged."""
    rm.write_manifest(step="filter", data_dir=str(tmp_path), effective_conditions={"na_conc": 50.0})
    rm.write_manifest(
        step="optimize", data_dir=str(tmp_path), effective_conditions={"na_conc": 75.0}
    )
    assert rm.read_effective_conditions(str(tmp_path)) == {"na_conc": 75.0}


def test_reading_a_step_that_never_ran_is_none(tmp_path):
    rm.write_manifest(step="filter", data_dir=str(tmp_path), effective_conditions={"na_conc": 50.0})
    assert rm.read_effective_conditions(str(tmp_path), step="optimize") is None
