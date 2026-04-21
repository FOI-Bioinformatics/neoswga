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


def test_doctor_aware_flag_uses_class_attribute_not_source_grep():
    """Phase 15D — the flag must be set on the optimizer class itself so it
    survives refactors that rename self.conditions references. Verifies
    ADDITIVE_AWARE attribute is True on the five known selection-aware
    wrappers and False on the coverage-only ones."""
    from neoswga.core.optimizer_factory import OptimizerRegistry
    from neoswga.core import unified_optimizer as _uo
    _uo._ensure_optimizers_registered()

    expected_aware = {"network", "hybrid", "background-aware",
                      "genetic", "equiphi29"}
    expected_coverage_only = {"greedy", "dominating-set", "weighted-set-cover",
                              "tiling", "clique"}

    for name in expected_aware:
        cls = OptimizerRegistry.get(name)
        assert getattr(cls, "ADDITIVE_AWARE", False) is True, (
            f"{name} ({cls.__name__}) should have ADDITIVE_AWARE = True"
        )

    for name in expected_coverage_only:
        cls = OptimizerRegistry.get(name)
        assert getattr(cls, "ADDITIVE_AWARE", False) is False, (
            f"{name} ({cls.__name__}) is coverage-only and should have "
            f"ADDITIVE_AWARE = False"
        )


def test_doctor_json_distinguishes_aware_vs_coverage_only():
    """End-to-end assertion through the CLI: the doctor JSON output
    partitions known-aware vs known-coverage-only optimizers correctly."""
    r = _run(["doctor", "--json"])
    payload = json.loads(r.stdout)
    by_name = {o["name"]: o["additive_aware"] for o in payload["optimizers"]}
    for name in ("network", "hybrid", "background-aware", "genetic", "equiphi29"):
        assert by_name.get(name) is True, (
            f"{name} should be additive-aware in doctor JSON; got {by_name.get(name)}"
        )
    for name in ("greedy", "dominating-set", "tiling", "clique"):
        assert by_name.get(name) is False, (
            f"{name} should be coverage-only in doctor JSON; got {by_name.get(name)}"
        )
