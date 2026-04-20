"""Tests for scripts/csv_to_experimental_results_json.py."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "csv_to_experimental_results_json.py"


def test_script_exists():
    assert SCRIPT.is_file()


def test_converts_minimal_csv(tmp_path):
    csv_path = tmp_path / "lab.csv"
    csv_path.write_text(
        "primer_set_id,primers,enrichment_fold,uniformity,notes\n"
        "set_A,ATCGATCGATCG;GCGTAGCATAGC;TATACGCATGGA,12.4,0.18,phi29+betaine\n"
        "set_B,ATCGATCGATCG;GCGTAGCATAGC;CAGCTAGTCATT,6.1,0.31,no additive\n"
    )
    out = tmp_path / "results.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(csv_path), "-o", str(out)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr

    payload = json.loads(out.read_text())
    assert len(payload) == 2
    first = payload[0]
    # Required fields
    assert first["primer_set"] == ["ATCGATCGATCG", "GCGTAGCATAGC", "TATACGCATGGA"]
    assert abs(first["enrichment_fold"] - 12.4) < 1e-9
    assert abs(first["uniformity_score"] - 0.18) < 1e-9
    # Optional fields default in
    assert first["temperature"] == 30.0
    assert first["time_hours"] == 4.0
    assert first["primer_concentration"] == 1.0
    assert first["passed_qc"] is True
    assert first["notes"] == "phi29+betaine"


def test_respects_optional_temperature(tmp_path):
    csv_path = tmp_path / "lab.csv"
    csv_path.write_text(
        "primer_set_id,primers,enrichment_fold,uniformity,temperature,passed_qc\n"
        "set_A,ATCGATCGATCG;GCGTAGCATAGC,15.2,0.22,43.0,true\n"
    )
    out = tmp_path / "results.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(csv_path), "-o", str(out)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(out.read_text())
    assert payload[0]["temperature"] == 43.0
    assert payload[0]["passed_qc"] is True


def test_output_is_activelearn_compatible(tmp_path):
    """JSON emitted should be loadable by ExperimentalResult.from_dict."""
    from neoswga.core.active_learning import ExperimentalResult

    csv_path = tmp_path / "lab.csv"
    csv_path.write_text(
        "primer_set_id,primers,enrichment_fold,uniformity\n"
        "set_A,ATCGATCGATCG;GCGTAGCATAGC,15.2,0.22\n"
    )
    out = tmp_path / "results.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(csv_path), "-o", str(out)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0
    payload = json.loads(out.read_text())

    # Should round-trip through ExperimentalResult
    for record in payload:
        r = ExperimentalResult.from_dict(record)
        assert r.primer_set == record["primer_set"]
        assert r.enrichment_fold == record["enrichment_fold"]


def test_missing_required_column_errors(tmp_path):
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("primers,enrichment_fold\nACGT,1.0\n")  # missing uniformity
    out = tmp_path / "out.json"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "-i", str(csv_path), "-o", str(out)],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode != 0
    assert "required columns" in result.stderr.lower()
