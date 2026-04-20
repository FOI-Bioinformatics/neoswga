"""Opt-in scale test: end-to-end pipeline on a synthetic bacterial-scale
target, with a memory-usage guardrail. Not run in default `pytest`; enable
with `pytest -m scale`. Used by the nightly CI workflow.

The synthetic genome is ~2 Mb, small enough to run in ~30-60 s on CI while
still exercising the code path that large real bacteria trigger.
"""

import os
import random
import resource
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


def _has_jellyfish():
    return shutil.which("jellyfish") is not None


def _write_random_fasta(path: Path, length: int, seed: int = 0, gc: float = 0.5):
    rng = random.Random(seed)
    alphabet = []
    # Build an alphabet weighted by the target GC fraction.
    for _ in range(int(gc * 100)):
        alphabet.append('G' if rng.random() < 0.5 else 'C')
    for _ in range(int((1 - gc) * 100)):
        alphabet.append('A' if rng.random() < 0.5 else 'T')
    with path.open("w") as fh:
        fh.write(f">synthetic_{length}bp\n")
        # Write 80 bp per line for valid FASTA
        chunk = []
        for _ in range(length):
            chunk.append(rng.choice(alphabet))
            if len(chunk) >= 80:
                fh.write("".join(chunk) + "\n")
                chunk = []
        if chunk:
            fh.write("".join(chunk) + "\n")


@pytest.mark.scale
@pytest.mark.slow
def test_bacterial_scale_pipeline_memory(tmp_path):
    """End-to-end filter on a 2 Mb synthetic target completes and stays
    below a 2 GB RSS budget. Guards against accidental O(N^2) regressions
    in position handling or k-mer scanning."""
    if not _has_jellyfish():
        pytest.skip("jellyfish not installed; skipping scale test")

    target = tmp_path / "target.fna"
    _write_random_fasta(target, length=2_000_000, seed=42, gc=0.5)

    params = {
        "schema_version": 1,
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(target)],
        "fg_prefixes": [str(tmp_path / "target")],
        "polymerase": "phi29",
        "min_k": 10,
        "max_k": 12,
        "min_fg_freq": 1e-5,
        "max_bg_freq": 5e-6,
        "min_tm": 15,
        "max_tm": 45,
        "max_gini": 0.6,
        "max_primer": 200,
        "cpus": 1,
        "fg_circular": False,
    }
    import json
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps(params))

    original_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        # Run count-kmers via the CLI (jellyfish required)
        result = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified",
             "count-kmers", "-j", str(params_file)],
            capture_output=True, text=True, timeout=600,
        )
        assert result.returncode == 0, (
            f"count-kmers failed: {result.stderr[-2000:]}"
        )

        result = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified",
             "filter", "-j", str(params_file)],
            capture_output=True, text=True, timeout=600,
        )
        assert result.returncode == 0, (
            f"filter failed: {result.stderr[-2000:]}"
        )
    finally:
        os.chdir(original_cwd)

    # Memory guardrail: ~2 GB RSS max. ru_maxrss is kilobytes on Linux, bytes on macOS.
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    if sys.platform == "darwin":
        max_rss_kb = usage.ru_maxrss / 1024
    else:
        max_rss_kb = usage.ru_maxrss
    max_rss_mb = max_rss_kb / 1024
    assert max_rss_mb < 2048, (
        f"Pipeline exceeded 2 GB RSS: {max_rss_mb:.0f} MB. "
        f"Check for O(N^2) regression in position handling or k-mer scanning."
    )

    step2_csv = tmp_path / "step2_df.csv"
    assert step2_csv.is_file(), "step2_df.csv should exist after scale run"
