"""`optimize` refuses a `step3_df.csv` prepared from a different `step2_df.csv`.

`prepare-candidates` records the checksum of the `step2_df.csv` it read and of
the `step3_df.csv` it wrote, and nothing compared them. After `filter` was
re-run with new thresholds and `prepare-candidates` skipped, `optimize` ran on
the old pool and logged its size as current
(2026-11-23 sweep, tier 1 item 4).
"""

import json
import sys

import pytest

from neoswga import cli_unified
from neoswga.core.run_manifest import step3_matches_step2, write_manifest


def _prepared(directory):
    (directory / "step2_df.csv").write_text("primer\nACGTACGTAC\nTTGACCATGA\n")
    (directory / "step3_df.csv").write_text("primer\nACGTACGTAC\nTTGACCATGA\n")
    write_manifest(
        "prepare-candidates",
        str(directory),
        input_files=[str(directory / "step2_df.csv")],
        output_files=[str(directory / "step3_df.csv")],
    )
    return directory


def test_matching_files_match(tmp_path):
    assert step3_matches_step2(str(_prepared(tmp_path))) is True


def test_a_refiltered_step2_does_not_match(tmp_path):
    _prepared(tmp_path)
    (tmp_path / "step2_df.csv").write_text("primer\nACGTACGTAC\n")
    assert step3_matches_step2(str(tmp_path)) is False


def test_no_manifest_is_unknown(tmp_path):
    (tmp_path / "step2_df.csv").write_text("primer\n")
    (tmp_path / "step3_df.csv").write_text("primer\n")
    assert step3_matches_step2(str(tmp_path)) is None


def test_a_step3_no_entry_describes_is_unknown(tmp_path):
    _prepared(tmp_path)
    (tmp_path / "step3_df.csv").write_text("primer\nGGCATTACGA\n")
    assert step3_matches_step2(str(tmp_path)) is None


def test_a_relative_data_dir_still_matches(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "step2_df.csv").write_text("primer\nA\n")
    (tmp_path / "results" / "step3_df.csv").write_text("primer\nA\n")
    write_manifest(
        "prepare-candidates",
        "results",
        input_files=["results/step2_df.csv"],
        output_files=["results/step3_df.csv"],
    )
    assert step3_matches_step2(str(tmp_path / "results")) is True


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    results = tmp_path / "results"
    results.mkdir()
    (tmp_path / "params.json").write_text(json.dumps({"data_dir": "results"}))
    monkeypatch.chdir(tmp_path)
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", None, raising=False)
    return results


def test_optimize_refuses_a_stale_step3(workspace, monkeypatch):
    _prepared(workspace)
    (workspace / "step2_df.csv").write_text("primer\nACGTACGTAC\n")
    ran = []
    monkeypatch.setattr(cli_unified, "run_step4", lambda args: ran.append(1))
    monkeypatch.setattr(sys, "argv", ["neoswga", "optimize", "-j", "params.json"])
    with pytest.raises(SystemExit):
        cli_unified.main()
    assert ran == []


def test_optimize_runs_without_provenance(workspace, monkeypatch):
    (workspace / "step2_df.csv").write_text("primer\n")
    (workspace / "step3_df.csv").write_text("primer\n")
    ran = []
    monkeypatch.setattr(cli_unified, "run_step4", lambda args: ran.append(1))
    monkeypatch.setattr(sys, "argv", ["neoswga", "optimize", "-j", "params.json"])
    cli_unified.main()
    assert ran == [1]
