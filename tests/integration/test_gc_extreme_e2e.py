"""GC-extreme scenario integration tests.

Uses synthetic AT-rich and GC-rich plasmid-scale FASTA files to exercise
the adaptive GC filter and extreme-GC-aware parameter path through
get_params() and step2. Marked `scale` because the count-kmers step needs
jellyfish; when jellyfish is unavailable the tests skip.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


def _has_jellyfish():
    return shutil.which("jellyfish") is not None


def _write_extreme_gc_fasta(path: Path, length: int, gc: float, seed: int = 0):
    """Write a synthetic FASTA with approximately the target GC fraction."""
    import random
    rng = random.Random(seed)
    with path.open("w") as fh:
        fh.write(f">synthetic_gc{int(gc*100)}\n")
        line = []
        for _ in range(length):
            if rng.random() < gc:
                line.append('G' if rng.random() < 0.5 else 'C')
            else:
                line.append('A' if rng.random() < 0.5 else 'T')
            if len(line) >= 80:
                fh.write("".join(line) + "\n")
                line = []
        if line:
            fh.write("".join(line) + "\n")


def _run_cli(cmd, cwd):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *cmd],
        capture_output=True, text=True, timeout=300, cwd=cwd,
    )


@pytest.mark.scale
@pytest.mark.slow
def test_at_rich_scenario_auto_adapts_gc_window(tmp_path):
    """AT-rich (~25% GC) target triggers adaptive GC window."""
    if not _has_jellyfish():
        pytest.skip("jellyfish not installed")

    target = tmp_path / "target_at_rich.fna"
    _write_extreme_gc_fasta(target, length=50_000, gc=0.25, seed=1)

    params = {
        "schema_version": 1,
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(target)],
        "fg_prefixes": [str(tmp_path / "target_at_rich")],
        "polymerase": "phi29",
        "min_k": 8,
        "max_k": 10,
        "min_fg_freq": 1e-4,
        "max_bg_freq": 1e-3,
        "min_tm": 10,
        "max_tm": 45,
        "max_gini": 0.8,
        "max_primer": 200,
        "cpus": 1,
        "fg_circular": False,
    }
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(params))

    result = _run_cli(["count-kmers", "-j", str(params_file)], cwd=tmp_path)
    assert result.returncode == 0, result.stderr[-1000:]

    result = _run_cli(["filter", "-j", str(params_file)], cwd=tmp_path)
    assert result.returncode == 0, result.stderr[-1000:]

    # Log must announce adaptive engagement (either captured by CLI or step2).
    combined = result.stdout + result.stderr
    assert "Extreme GC" in combined or "Adaptive GC filtering" in combined, (
        "Expected an adaptive-GC log line for a 25% GC genome"
    )

    step2_csv = tmp_path / "step2_df.csv"
    assert step2_csv.is_file()
    # Some primers should pass even with narrow GC window thanks to adaptive
    assert step2_csv.stat().st_size > 100, "No primers survived filter on AT-rich target"


@pytest.mark.scale
@pytest.mark.slow
def test_gc_rich_scenario_auto_adapts_gc_window(tmp_path):
    """GC-rich (~70% GC) target triggers adaptive GC window."""
    if not _has_jellyfish():
        pytest.skip("jellyfish not installed")

    target = tmp_path / "target_gc_rich.fna"
    _write_extreme_gc_fasta(target, length=50_000, gc=0.70, seed=2)

    params = {
        "schema_version": 1,
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(target)],
        "fg_prefixes": [str(tmp_path / "target_gc_rich")],
        "polymerase": "equiphi29",
        "reaction_temp": 43.0,
        "mg_conc": 10.0,
        "dmso_percent": 5.0,
        "betaine_m": 1.0,
        "min_k": 10,
        "max_k": 12,
        "min_fg_freq": 1e-4,
        "max_bg_freq": 1e-3,
        "min_tm": 30,
        "max_tm": 75,
        "max_gini": 0.8,
        "max_primer": 200,
        "cpus": 1,
        "fg_circular": False,
    }
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(params))

    result = _run_cli(["count-kmers", "-j", str(params_file)], cwd=tmp_path)
    assert result.returncode == 0, result.stderr[-1000:]

    result = _run_cli(["filter", "-j", str(params_file)], cwd=tmp_path)
    assert result.returncode == 0, result.stderr[-1000:]

    combined = result.stdout + result.stderr
    assert "Extreme GC" in combined or "Adaptive GC filtering" in combined, (
        "Expected an adaptive-GC log line for a 70% GC genome"
    )

    step2_csv = tmp_path / "step2_df.csv"
    assert step2_csv.is_file()
