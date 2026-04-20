"""Smoke tests for `neoswga doctor`."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _run(args, cwd=None, timeout=30):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True, text=True, timeout=timeout, cwd=cwd,
    )


def test_doctor_help():
    r = _run(["doctor", "--help"])
    assert r.returncode == 0, r.stderr
    assert "--json-file" in r.stdout or "-j" in r.stdout


def test_doctor_runs_without_params():
    r = _run(["doctor"])
    # Status depends on jellyfish availability — either way stdout
    # should have the header.
    assert "neoswga doctor" in r.stdout


def test_doctor_json_output_is_valid():
    r = _run(["doctor", "--json"])
    # Must be parseable JSON regardless of exit code
    payload = json.loads(r.stdout)
    assert "platform" in payload
    assert "optimizers" in payload
    assert "ok" in payload
    # At least one optimizer is registered
    assert len(payload["optimizers"]) > 0
    # Each optimizer carries the additive-aware flag
    for opt in payload["optimizers"]:
        assert "name" in opt and "additive_aware" in opt


def test_doctor_surfaces_params_summary():
    example = ROOT / "examples" / "plasmid_example" / "params.json"
    if not example.is_file():
        pytest.skip("plasmid example not available")
    r = _run(["doctor", "-j", str(example), "--json"])
    payload = json.loads(r.stdout)
    assert payload["params"] is not None
    assert payload["params"]["num_targets"] >= 1


def test_doctor_fails_on_missing_params_file(tmp_path):
    bogus = tmp_path / "does_not_exist.json"
    r = _run(["doctor", "-j", str(bogus)])
    assert r.returncode != 0
    # Blocking message should appear in stdout (not stderr, since we print)
    combined = r.stdout + r.stderr
    assert "not found" in combined.lower() or "blocking" in combined.lower()


def test_doctor_reports_at_least_three_additive_aware_optimizers():
    """After Phase 7 plumbing, at minimum network / hybrid / genetic are
    additive-aware. Guards against accidental regression of that work."""
    r = _run(["doctor", "--json"])
    payload = json.loads(r.stdout)
    aware_names = {o["name"] for o in payload["optimizers"] if o["additive_aware"]}
    # Must include the three known additive-aware optimizers
    assert {"network", "hybrid", "genetic"}.issubset(aware_names), (
        f"Expected network/hybrid/genetic to show as additive-aware; got {aware_names}"
    )
