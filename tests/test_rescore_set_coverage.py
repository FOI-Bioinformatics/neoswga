"""Phase 12D: rescore-set includes coverage + per-target + bg under the
chosen reaction conditions.

The previous rescore-set was a half-measure — it computed per-primer
quality scores under new conditions but left users blind to whether the
set's coverage / specificity would actually hold under those conditions.
The added `coverage` block closes that gap.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = ROOT / "examples" / "plasmid_example"


def _run(args, cwd=None, timeout=120):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True, text=True, timeout=timeout, cwd=cwd,
    )


@pytest.fixture
def rescore_output():
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")
    result = _run(
        ["rescore-set", "-j", "params.json",
         "--primers", "AAACGCT", "CATCCGTAAG", "AGGAAAGGAC",
         "--betaine-m", "1.0", "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert result.returncode == 0, result.stderr
    return json.loads(result.stdout)


def test_rescore_output_has_coverage_block(rescore_output):
    cov = rescore_output.get("coverage")
    assert cov is not None, "rescore-set output missing 'coverage' block"
    for key in ("fg_coverage", "per_target_coverage",
                "bg_coverage", "per_background_coverage",
                "selectivity_ratio", "extension_reach_bp"):
        assert key in cov, f"coverage block missing key {key}"


def test_rescore_per_target_is_dict(rescore_output):
    per_target = rescore_output["coverage"]["per_target_coverage"]
    assert isinstance(per_target, dict)
    # Plasmid example has one target (pcDNA), so one key.
    assert len(per_target) >= 1
    for val in per_target.values():
        assert 0.0 <= val <= 1.0


def test_rescore_extension_reach_reflects_polymerase():
    if not EXAMPLE_DIR.is_dir():
        pytest.skip()
    # phi29 -> 70000
    r = _run(
        ["rescore-set", "-j", "params.json",
         "--primers", "AAACGCT", "CATCCGTAAG",
         "--polymerase", "phi29", "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert r.returncode == 0
    phi = json.loads(r.stdout)
    # bst -> 2000
    r = _run(
        ["rescore-set", "-j", "params.json",
         "--primers", "AAACGCT", "CATCCGTAAG",
         "--polymerase", "bst", "--reaction-temp", "63", "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert r.returncode == 0
    bst = json.loads(r.stdout)

    assert phi["coverage"]["extension_reach_bp"] > bst["coverage"]["extension_reach_bp"], (
        "phi29 should have larger extension reach than bst"
    )


def test_rescore_selectivity_ratio_positive(rescore_output):
    sel = rescore_output["coverage"]["selectivity_ratio"]
    assert sel > 0.0
