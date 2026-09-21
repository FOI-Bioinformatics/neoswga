"""A design may mix oligo lengths, and this proves it through a real run.

An external review claimed neoswga is "restricted to fixed k-mer lengths
(6-12 bp)". Both halves are wrong: `params.schema.json` allows k of 4 to 30,
and 6-12 is phi29's default window, not a limit (bst is 15-25, equiphi29
10-18). More to the point, nothing in the pipeline requires ONE k. The scan
writes `{prefix}_{k}mer_positions.h5` per length, `PositionCache.load` groups
its primers by length and opens the matching file, and the dimer screen codes
t-mers so unequal lengths compare without special handling.

What was missing was a demonstration. Every saved run in this repository used a
single k, so "mixed length works" was an argument rather than a measurement,
and an argument is what the review was entitled to disbelieve.

The reach is shortened deliberately. At phi29's realistic 3 kb reach one primer
covers the packaged 6 kb plasmid completely, Stage 1 returns after one pick,
and no panel of any composition is exercised -- the same limitation Known Issue
14 records for `expand-primers`. At 300 bp selection has to choose, and it
chooses across lengths.

This is an end-to-end test for the reason `test_strict_design_pipeline.py`
gives: the per-length handling has to survive argument parsing, parameter
resolution, four steps and three reporting commands, and a unit test that
constructs a mixed pool directly would pass whether or not a real run can
produce one.
"""

import collections
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("jellyfish") is None,
    reason="the pipeline's first step shells out to jellyfish",
)

PACKAGED = Path(__file__).resolve().parent.parent.parent / "neoswga" / "core" / "smoke"

MIN_K = 7
MAX_K = 11


def run(args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True,
        text=True,
        cwd=str(cwd),
        timeout=900,
    )


@pytest.fixture(scope="module")
def designed(tmp_path_factory):
    """One mixed-length design, taken through all four steps."""
    workspace = tmp_path_factory.mktemp("mixed_length") / "work"
    workspace.mkdir()
    for name in ("pcDNA.fasta", "pLTR.fasta"):
        shutil.copy(PACKAGED / name, workspace / name)

    params = {
        "schema_version": 2,
        "fg_genomes": ["pcDNA.fasta"],
        "bg_genomes": ["pLTR.fasta"],
        "fg_prefixes": ["pcDNA"],
        "bg_prefixes": ["pLTR"],
        "data_dir": "./",
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_k": MIN_K,
        "max_k": MAX_K,
        "min_fg_freq": 1e-5,
        "max_bg_freq": 5e-6,
        "min_tm": 15,
        "max_tm": 45,
        "max_gini": 0.6,
        "max_primer": 500,
        "max_dimer_bp": 3,
        "max_self_dimer_bp": 4,
        # Short enough that one primer cannot saturate the target, so Stage 1
        # must assemble a panel. See the module docstring.
        "coverage_reach": 300,
        "num_primers": 8,
        "target_set_size": 8,
        "iterations": 8,
        "max_sets": 5,
        "cpus": 1,
        "fg_circular": True,
        "bg_circular": False,
    }
    (workspace / "params.json").write_text(json.dumps(params, indent=2))

    for step in ("count-kmers", "filter", "prepare-candidates", "optimize"):
        result = run([step, "-j", "params.json"], workspace)
        assert result.returncode == 0, (
            f"{step} failed on a mixed-length configuration "
            f"(k {MIN_K}-{MAX_K}):\n{result.stdout}\n{result.stderr}"
        )
    return workspace


def lengths_of(sequences):
    return collections.Counter(len(str(s)) for s in sequences)


def delivered_panel(workspace):
    from neoswga.core.delivered_set import read_delivered_set

    return list(read_delivered_set(workspace / "step4_improved_df.csv").primers)


# ---------------------------------------------------------------------------
# The pipeline produces, and keeps, more than one length
# ---------------------------------------------------------------------------


def test_the_scan_writes_a_position_index_per_length(designed):
    """`PositionCache.load` opens `{prefix}_{k}mer_positions.h5`, so each k needs one."""
    written = sorted(
        int(path.name.split("_")[1].removesuffix("mer"))
        for path in designed.glob("pcDNA_*mer_positions.h5")
    )

    assert len(written) > 1, f"only one length was indexed: {written}"
    assert set(written) <= set(range(MIN_K, MAX_K + 1)), written


def test_the_candidate_pool_carries_several_lengths(designed):
    import pandas as pd

    pool = pd.read_csv(designed / "step3_df.csv")["primer"].astype(str)
    counts = lengths_of(pool)

    assert len(counts) > 1, f"the pool collapsed to one length: {dict(counts)}"


def test_the_delivered_panel_mixes_lengths(designed):
    """The claim under test, at the only place it matters.

    A mixed CANDIDATE pool proves nothing on its own: selection could still
    pick one length throughout, and then a user would be right to say the tool
    does not deliver mixed-length panels.
    """
    counts = lengths_of(delivered_panel(designed))

    assert len(counts) > 1, (
        f"selection delivered a single-length panel: {dict(counts)}. "
        "Mixed-length support is unproven if no run produces one."
    )


# ---------------------------------------------------------------------------
# Everything downstream survives it
# ---------------------------------------------------------------------------


def test_the_delivered_panel_honours_the_dimer_limit_across_lengths(designed):
    """The screen codes t-mers, so unequal lengths need no special handling."""
    from neoswga.core.dimer_matrix import build

    panel = delivered_panel(designed)
    violating = list(build(panel, 3).flagged_pairs())

    assert not violating, [(panel[i], panel[j]) for i, j in violating]


def test_the_summary_reports_measured_coverage_for_a_mixed_panel(designed):
    summary = json.loads((designed / "step4_improved_df_summary.json").read_text())
    metrics = summary.get("metrics", summary)

    assert 0.0 < float(metrics["fg_coverage"]) <= 1.0
    assert summary.get("unindexed_candidates", 0) == 0, (
        "a candidate the foreground index could not place covers nothing and "
        "is invisible to selection; on a mixed run that would most likely mean "
        "one length was never indexed"
    )


@pytest.mark.parametrize("command", [["export", "--format", "fasta"], ["interpret"], ["report"]])
def test_the_reporting_commands_accept_a_mixed_panel(designed, command):
    result = run([*command, "-d", "."], designed)

    assert result.returncode == 0, f"{command[0]} failed:\n{result.stdout}\n{result.stderr}"


def test_the_exported_file_is_the_mixed_panel_and_nothing_else(designed):
    """Lengths survive export, and no alternative set is pooled into it.

    The second half is `core/delivered_set.py`'s subject and is asserted here
    too, because export is where a length mix would most plausibly be lost and
    where a pooled file would most plausibly be missed.
    """
    result = run(["export", "-d", ".", "--format", "fasta", "-o", "out"], designed)
    assert result.returncode == 0, result.stderr

    exported = [
        line.strip()
        for line in (designed / "out" / "SWGA_primers.fasta").read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ]

    assert exported == delivered_panel(designed)
    assert len(lengths_of(exported)) > 1, dict(lengths_of(exported))


def test_the_fasta_header_records_each_oligos_own_length(designed):
    """A mixed panel is only orderable if the file says which is which."""
    run(["export", "-d", ".", "--format", "fasta", "-o", "out"], designed)

    headers = [
        line
        for line in (designed / "out" / "SWGA_primers.fasta").read_text().splitlines()
        if line.startswith(">")
    ]

    assert headers, "no FASTA records were written"
    assert all("len=" in header for header in headers), headers[:3]
