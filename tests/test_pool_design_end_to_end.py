"""The whole design path, on a fixture small enough to check by hand.

Task 8 of the condition-aware pool design plan.

Every earlier task was verified against its own unit. This runs counting, QC,
the inventory, scoring and the design report through the real entry points on a
two-record synthetic reference with known planted positions, and checks the
numbers that the separate tasks each own:

- record geometry, because concatenating records without separators produced
  k-mers present in neither and coverage windows crossing contigs;
- the six counts, because the funnel reported one number and labelled it as if
  it were the pool;
- the candidate inventory, because `step2_df.csv` is a shortlist.

Deliberately small and tracked. The plan's acceptance is that a clean checkout
can exercise the complete path without downloading large references; the
Wolbachia benchmark is a separate step that needs its declared references.
"""

import json
import shutil
from pathlib import Path

import pytest

from neoswga.core.kmer_counter import check_jellyfish_available

pytestmark = pytest.mark.skipif(
    not check_jellyfish_available(), reason="jellyfish is required to count k-mers"
)

FILLER = "TGTGTGTGTGTG"
# The junction. Record one ends in six A, record two begins with six C, and
# neither run occurs anywhere else, so this 12-mer exists ONLY across the join.
# A first attempt used a junction that also appeared inside record two, which
# made the test assert something untrue rather than catching anything.
JUNCTION = "AAAAAACCCCCC"


def _fasta(path: Path) -> None:
    """Two records, so record geometry is exercised rather than assumed."""
    first = ("ACGTTGCATTACG" + FILLER) * 12 + JUNCTION[:6]
    second = JUNCTION[6:] + ("GGCATTACGGCA" + FILLER) * 12
    path.write_text(f">contig_one\n{first}\n>contig_two\n{second}\n")


def _background(path: Path) -> None:
    path.write_text(">bg_one\n" + ("CCCCGGGGAAAA" + FILLER) * 12 + "\n")


@pytest.fixture
def design_dir(tmp_path):
    work = tmp_path / "work"
    work.mkdir()
    _fasta(work / "target.fna")
    _background(work / "background.fna")
    params = {
        "schema_version": 2,
        "fg_genomes": [str(work / "target.fna")],
        "bg_genomes": [str(work / "background.fna")],
        "fg_prefixes": [str(work / "target")],
        "bg_prefixes": [str(work / "background")],
        "data_dir": str(work),
        "min_k": 12,
        "max_k": 12,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_tm": 10.0,
        "max_tm": 80.0,
        "min_gini_sites": 1,
        "max_gini": 1.0,
        "min_fg_freq": 0.0,
        "max_bg_freq": 1.0,
        "max_primer": 5,
        "num_primers": 3,
        "cpus": 1,
    }
    (work / "params.json").write_text(json.dumps(params, indent=2))
    return work


def _run(step, work):
    import subprocess
    import sys

    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", step, "-j", "params.json"],
        cwd=work,
        capture_output=True,
        text=True,
    )


def test_the_design_path_runs_end_to_end_on_a_small_fixture(design_dir):
    for step in ("count-kmers", "filter", "score"):
        result = _run(step, design_dir)
        assert result.returncode == 0, f"{step} failed:\n{result.stderr[-2000:]}"

    assert (design_dir / "step2_df.csv").is_file()
    assert (design_dir / "step3_df.csv").is_file()
    assert (design_dir / "candidate_inventory.sqlite").is_file()


def test_plan_pool_forwards_search_controls_to_real_optimizer(design_dir, monkeypatch):
    import argparse

    from Bio import SeqIO

    from neoswga.cli.plan_pool import add_parsers, run_plan_pool
    from neoswga.core.optimizer_factory import OptimizerFactory

    for step in ("count-kmers", "filter", "score"):
        completed = _run(step, design_dir)
        assert completed.returncode == 0, completed.stderr
    path = design_dir / "params.json"
    params = json.loads(path.read_text())
    params.update(stage1_objective_width=8, swap_max_evaluations=17)
    for label, filename in (("fg", "target.fna"), ("bg", "background.fna")):
        params[f"{label}_seq_lengths"] = [
            sum(len(record.seq) for record in SeqIO.parse(design_dir / filename, "fasta"))
        ]
    path.write_text(json.dumps(params))
    seen = []
    create = OptimizerFactory.create

    def capture(*args, **kwargs):
        optimizer = create(*args, **kwargs)
        seen.append(optimizer.config)
        return optimizer

    monkeypatch.setattr(OptimizerFactory, "create", capture)
    parser = argparse.ArgumentParser()
    add_parsers(parser.add_subparsers())
    output = design_dir / "plan"
    args = parser.parse_args(
        [
            "plan-pool",
            "-j",
            str(path),
            "--candidates",
            str(design_dir / "step3_df.csv"),
            "--method",
            "dominating-set",
            "--max-size",
            "1",
            "--no-background",
            "-o",
            str(output),
        ]
    )
    run_plan_pool(args)
    assert seen and seen[0].stage1_objective_width == 8
    assert seen[0].swap_max_evaluations == 17
    assert (output / "pool_plan.json").is_file()


def test_the_inventory_holds_more_than_the_shortlist(design_dir):
    for step in ("count-kmers", "filter"):
        assert _run(step, design_dir).returncode == 0

    from neoswga.core.candidate_inventory import CandidateInventory

    shortlist = len((design_dir / "step2_df.csv").read_text().strip().split("\n")) - 1
    with CandidateInventory(design_dir / "candidate_inventory.sqlite") as inventory:
        counts = inventory.counts()

    assert counts["hard_qc_passed"] >= shortlist, (
        "the inventory holds fewer candidates than the shortlist it was meant to "
        "be a superset of"
    )


def test_the_index_records_the_record_geometry(design_dir):
    for step in ("count-kmers", "filter"):
        assert _run(step, design_dir).returncode == 0

    from neoswga.core.position_cache import PositionCache

    prefix = str(design_dir / "target")
    primers = [
        line.split(",")[1]
        for line in (design_dir / "step2_df.csv").read_text().strip().split("\n")[1:]
    ]
    cache = PositionCache([prefix], primers)

    starts = cache.get_record_starts(prefix)
    assert starts, "a new index carries no record geometry"
    assert starts[0] == 0
    assert len(starts) == 2, "the two records should both be recorded"
    # And the requirement a new design enforces is satisfied by what we just built.
    cache.require_record_metadata([prefix])


def test_no_reported_position_spans_the_record_join(design_dir):
    """The planted junction must not produce matches that exist in neither record."""
    for step in ("count-kmers", "filter"):
        assert _run(step, design_dir).returncode == 0

    from neoswga.core import string_search
    from neoswga.core.position_cache import PositionCache

    prefix = str(design_dir / "target")
    # A cache is built over primers, and the record geometry is read while their
    # index is opened, so it needs at least one primer to load anything.
    primers = [
        line.split(",")[1]
        for line in (design_dir / "step2_df.csv").read_text().strip().split("\n")[1:]
    ]
    starts = PositionCache([prefix], primers).get_record_starts(prefix)
    boundary = starts[1]

    sequence = string_search.get_cached_genome_sequence(str(design_dir / "target.fna"))
    assert sequence[boundary - 6 : boundary + 6] == JUNCTION
    # Guard the guard: the junction must exist ONLY across the join, or the
    # assertion below would pass for the wrong reason.
    assert JUNCTION not in sequence[:boundary]
    assert JUNCTION not in sequence[boundary:]

    found = string_search.get_all_positions_multi_k(
        {12: [JUNCTION]}, str(design_dir / "target.fna"), circular=False
    )

    assert (
        found[JUNCTION] == []
    ), "a 12-mer formed only by the join between two records was reported"
