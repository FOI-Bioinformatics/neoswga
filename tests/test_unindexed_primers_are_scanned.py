"""Primers the position index does not hold are scanned, not read as absent.

`analyze-coverage` and `expand-primers` built their `PositionCache` with the
default `on_missing="warn"`. A primer outside the filtered shortlist, which is
what a user's own primer usually is, then read as binding nowhere:
`analyze-coverage --primers CTGCTCTGAT` on the plasmid example wrote one
0-6,157 bp gap though the primer binds pcDNA at 47 and 2950, and
`expand-primers` proposed a primer for that gap (2026-11-23 sweep, tier 2).
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.cli._common import design_position_cache
from neoswga.core.position_cache import MissingPositionsError

INDEXED, UNINDEXED = "ACGTACGTAC", "GGCATTACGA"


@pytest.fixture
def reference(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">t\n" + "T" * 100 + UNINDEXED + "T" * 100 + INDEXED + "T" * 100 + "\n")
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
        f.create_dataset(INDEXED, data=np.array([210], dtype=np.int64))
    return prefix, str(fasta)


def test_an_unindexed_primer_is_found_by_scanning(reference):
    prefix, fasta = reference
    cache = design_position_cache([prefix], [INDEXED, UNINDEXED], [fasta])
    assert 100 in cache.get_positions(prefix, UNINDEXED, "both").tolist()


def test_without_genomes_an_unindexed_primer_is_refused(reference):
    prefix, _fasta = reference
    with pytest.raises(MissingPositionsError):
        design_position_cache([prefix], [INDEXED, UNINDEXED], None)


EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "plasmid_example"


@pytest.mark.usefixtures("primed_plasmid_example")
def test_analyze_coverage_does_not_report_a_bound_primer_as_one_gap(tmp_path):
    for entry in EXAMPLE.iterdir():
        if entry.is_file():
            shutil.copy2(entry, tmp_path)
    result = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "analyze-coverage", "-j", "params.json",
         "--primers", "CTGCTCTGAT", "-o", "out"],
        cwd=tmp_path, capture_output=True, text=True, timeout=300,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    gaps = json.loads((tmp_path / "out" / "coverage_gaps.json").read_text())["gaps"]
    spans = [(g.get("start"), g.get("end")) for g in gaps]
    assert (0, 6157) not in spans and len(gaps) != 1, spans
