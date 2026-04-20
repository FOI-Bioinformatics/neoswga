"""Smoke tests for Phase 8 CLIs: swap-primer, contract-set, rescore-set.

Focuses on surface-level behaviour (help parses, flags register, end-to-end
runs produce valid JSON) rather than deep algorithm coverage (the underlying
modules already have dedicated tests).
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = ROOT / "examples" / "plasmid_example"


def _run(args, cwd=None, timeout=120):
    return subprocess.run(
        [sys.executable, '-m', 'neoswga.cli_unified', *args],
        capture_output=True, text=True, timeout=timeout, cwd=cwd,
    )


@pytest.mark.parametrize("cmd", ["swap-primer", "contract-set", "rescore-set"])
def test_help_is_available(cmd):
    result = _run([cmd, "--help"])
    assert result.returncode == 0, result.stderr
    assert "--primers" in result.stdout


def test_rescore_set_additives_shift_tm():
    """Under 1.5 M betaine the effective Tm should drop vs no additive."""
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")

    primer_args = ["AAACGCT", "CATCCGTAAG", "AGGAAAGGAC"]

    r_plain = _run(
        ["rescore-set", "-j", "params.json",
         "--primers", *primer_args,
         "--polymerase", "phi29",
         "--reaction-temp", "30",
         "--betaine-m", "0.0",
         "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert r_plain.returncode == 0, r_plain.stderr
    plain = json.loads(r_plain.stdout)

    r_betaine = _run(
        ["rescore-set", "-j", "params.json",
         "--primers", *primer_args,
         "--polymerase", "phi29",
         "--reaction-temp", "30",
         "--betaine-m", "1.5",
         "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert r_betaine.returncode == 0, r_betaine.stderr
    betaine = json.loads(r_betaine.stdout)

    tm_plain_avg = sum(p["tm_effective_c"] for p in plain["primers"]) / len(plain["primers"])
    tm_betaine_avg = sum(p["tm_effective_c"] for p in betaine["primers"]) / len(betaine["primers"])

    delta = tm_plain_avg - tm_betaine_avg
    assert delta > 0.5, (
        f"1.5 M betaine should lower Tm by >0.5 C on average; got delta={delta:.2f}"
    )


def test_contract_set_smoke():
    """contract-set runs end-to-end and emits valid JSON with expected keys."""
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")
    result = _run(
        ["contract-set", "-j", "params.json",
         "--primers", "AAACGCT", "CATCCGTAAG", "AGGAAAGGAC",
         "--min-coverage", "0.3", "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    for key in ["original_set", "contracted_set", "removed_primers",
                "baseline_coverage", "final_coverage"]:
        assert key in payload, f"missing {key} in contract-set output"
    assert payload["final_coverage"] >= payload["min_coverage_threshold"] - 1e-9


def test_swap_primer_smoke():
    """swap-primer runs end-to-end on the plasmid example."""
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")
    result = _run(
        ["swap-primer", "-j", "params.json",
         "--primers", "AAACGCT", "CATCCGTAAG", "AGGAAAGGAC",
         "--max-swaps", "1", "--quiet"],
        cwd=str(EXAMPLE_DIR),
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads(result.stdout)
    for key in ["before_set", "after_set", "num_swaps",
                "metrics_before", "metrics_after"]:
        assert key in payload
    assert len(payload["after_set"]) == len(payload["before_set"])
