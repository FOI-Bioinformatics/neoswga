"""CLI surface-consistency smoke tests.

These confirm that the --blacklist flag exists on the pipeline subcommands
where blacklist filtering is actually applied (and that --kmer-range on
`suggest` prints a structured recommendation).

Blacklist is a count/filter-stage concept: blacklist k-mers are counted at
`count-kmers` and the frequency penalty is applied at `filter`. `score` and
`optimize` consume the already-filtered candidate pool, so they do not surface
a --blacklist flag (a flag there would be inert and misleading). The optimizer
still cross-checks the candidate pool against any configured blacklist via the
library-level guard in unified_optimizer.
"""

import subprocess
import sys

import pytest


def _run(args, timeout=30):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True,
        text=True,
        timeout=timeout,
    )


@pytest.mark.parametrize("cmd", ["count-kmers", "filter"])
def test_blacklist_flag_present_on_pipeline_commands(cmd):
    result = _run([cmd, "--help"])
    assert result.returncode == 0, result.stderr
    assert "--blacklist" in result.stdout, (
        f"'{cmd}' should accept --blacklist (blacklist is applied at this "
        f"stage); help text does not mention it."
    )


def test_suggest_kmer_range_phi29():
    result = _run(
        [
            "suggest",
            "--kmer-range",
            "--genome-size",
            "5000",
            "--genome-gc",
            "0.5",
            "--polymerase",
            "phi29",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert "min_k = 6" in result.stdout
    assert "max_k" in result.stdout


def test_suggest_kmer_range_equiphi29():
    result = _run(
        [
            "suggest",
            "--kmer-range",
            "--genome-size",
            "2000000",
            "--genome-gc",
            "0.5",
            "--polymerase",
            "equiphi29",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert "min_k = 10" in result.stdout
    assert "max_k = 18" in result.stdout


def test_suggest_kmer_range_bst_bacterial():
    result = _run(
        [
            "suggest",
            "--kmer-range",
            "--genome-size",
            "5000000",
            "--polymerase",
            "bst",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert "min_k = 15" in result.stdout
    assert "max_k = 25" in result.stdout


def test_suggest_kmer_range_requires_size_or_genome():
    result = _run(
        [
            "suggest",
            "--kmer-range",
            "--polymerase",
            "phi29",
        ]
    )
    assert result.returncode != 0
    assert "requires" in (result.stderr + result.stdout).lower()
