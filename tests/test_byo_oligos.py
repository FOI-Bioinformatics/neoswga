"""Guard: an externally-designed oligo set must be evaluable, and never silently zero.

The failure this prevents: `PositionCache` skipped primers absent from the HDF5
index without a log line, and `get_positions` returned an empty array rather
than raising. A primer set designed elsewhere therefore scored 0% coverage with
nothing to indicate the number was meaningless rather than bad - which silently
corrupted expand-primers, contract-set, rescore-set and predict-efficiency.

The fix has two halves, both tested here:

* `on_missing` distinguishes "never indexed, so coverage is UNKNOWN" from
  "scanned and genuinely has no sites, so coverage is ZERO". Those produced
  identical output before.
* `evaluate-set` scans the FASTA directly, so no prior count-kmers run is
  needed at all.
"""

import json
import os
import subprocess
import sys

import pytest

from neoswga.core.position_cache import MissingPositionsError, PositionCache


def _make_genome(n=2000, seed=12345):
    """Deterministic pseudo-random sequence.

    Deliberately NOT a repeated motif: with a repeating genome two probes taken
    at different offsets can be the identical k-mer, which collect_primers dedups
    and makes the test assert the wrong thing.
    """
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


_SEQ = _make_genome()
GENOME = ">chr1\n" + _SEQ + "\n"


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "target.fasta"
    path.write_text(GENOME)
    return path


@pytest.fixture
def seq():
    return _SEQ


def test_on_missing_rejects_unknown_mode():
    with pytest.raises(ValueError, match="on_missing"):
        PositionCache([], [], on_missing="nonsense")


def test_warn_mode_preserves_historical_behaviour(tmp_path, caplog):
    """Default must stay non-raising so existing callers are unaffected."""
    cache = PositionCache([str(tmp_path / "absent")], ["ATCGATCGATCG"])
    assert len(cache.get_positions(str(tmp_path / "absent"), "ATCGATCGATCG")) == 0
    assert cache.missing_primers, "the primer should be recorded as unresolved"


def test_error_mode_raises_instead_of_scoring_zero(tmp_path):
    """The silent-zero case must be loud when the caller asks for it."""
    with pytest.raises(MissingPositionsError) as excinfo:
        PositionCache([str(tmp_path / "absent")], ["ATCGATCGATCG"], on_missing="error")
    assert "ATCGATCGATCG" in excinfo.value.missing
    # The message must say why the number would be wrong, not just that it failed.
    assert "would be 0" in str(excinfo.value)


def test_scan_mode_finds_primers_absent_from_the_index(genome, seq, tmp_path):
    """The motivating case: a real primer that was simply never indexed."""
    present = seq[100:112]
    cache = PositionCache(
        [str(tmp_path / "target")],
        [present],
        genome_paths=[str(genome)],
        circular=True,
        on_missing="scan",
    )
    hits = cache.get_positions(str(tmp_path / "target"), present)
    assert len(hits) > 0, "a primer that occurs in the genome must be found by scanning"
    assert not cache.missing_primers, "nothing should remain unresolved after scanning"


def test_scan_mode_separates_unknown_from_genuinely_zero(genome, seq, tmp_path):
    """This distinction is the whole point: unknown != zero."""
    present = seq[200:212]
    absent = "TTTTTTTTTTTT"
    assert absent not in seq

    cache = PositionCache(
        [str(tmp_path / "target")],
        [present, absent],
        genome_paths=[str(genome)],
        circular=True,
        on_missing="scan",
    )
    assert cache.missing_primers == [], "everything was resolved one way or the other"
    assert cache.zero_site_primers == [absent], (
        "a primer with genuinely no sites must be reported as a real zero, "
        "distinct from one that was never indexed"
    )


def test_scan_mode_requires_aligned_genome_paths(tmp_path):
    with pytest.raises(ValueError, match="genome_paths"):
        PositionCache([str(tmp_path / "a")], ["ATCGATCGATCG"], on_missing="scan")


def test_scan_mode_rejects_mismatched_lengths(genome, tmp_path):
    with pytest.raises(ValueError, match="align"):
        PositionCache(
            [str(tmp_path / "a"), str(tmp_path / "b")],
            ["ATCGATCGATCG"],
            genome_paths=[str(genome)],
            on_missing="scan",
        )


# ----------------------------------------------------------------------
# evaluate-set end to end
# ----------------------------------------------------------------------


def _run(args, cwd):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=300,
    )


def test_evaluate_set_works_with_no_params_json(tmp_path, genome, seq):
    """The headline capability: oligos + genome, nothing else."""
    primers_file = tmp_path / "primers.txt"
    primers = [seq[100:112], seq[400:412], "TTTTTTTTTTTT"]
    primers_file.write_text("\n".join(primers) + "\n")

    proc = _run(
        [
            "evaluate-set",
            "--primers-file",
            str(primers_file),
            "--genome",
            str(genome),
            "-o",
            str(tmp_path / "out"),
        ],
        cwd=str(tmp_path),
    )
    assert proc.returncode == 0, proc.stderr

    report = json.loads((tmp_path / "out" / "evaluation.json").read_text())
    assert report["num_primers"] == 3
    assert report["total_binding_sites"] > 0, (
        "primers taken from the genome must register binding sites; zero here "
        "would mean the scan fallback is not working"
    )
    assert report["fg_coverage"] > 0.0
    assert "TTTTTTTTTTTT" in report["zero_site_primers"]
    assert report["mean_binding_distance_bp"] is not None


def test_evaluate_set_reports_mean_binding_distance(tmp_path, genome, seq):
    """Binding-site density is the literature's dominant success predictor."""
    primers_file = tmp_path / "p.txt"
    primers_file.write_text(f"{seq[100:112]}\n{seq[400:412]}\n")

    proc = _run(
        [
            "evaluate-set",
            "--primers-file",
            str(primers_file),
            "--genome",
            str(genome),
            "-o",
            str(tmp_path / "o"),
        ],
        cwd=str(tmp_path),
    )
    assert proc.returncode == 0, proc.stderr
    report = json.loads((tmp_path / "o" / "evaluation.json").read_text())

    expected = round(report["genome_bp"] / report["total_binding_sites"])
    assert report["mean_binding_distance_bp"] == expected


def test_evaluate_set_needs_a_source(tmp_path):
    """Neither -j nor --genome should fail with an explanation, not a traceback."""
    proc = _run(
        ["evaluate-set", "--primers", "ATCGATCGATCG", "-o", str(tmp_path / "o")],
        cwd=str(tmp_path),
    )
    assert proc.returncode != 0
    assert "--genome" in (proc.stderr + proc.stdout)


def test_evaluate_set_is_registered_and_grouped():
    """An ungrouped command triggers a warning in create_parser."""
    from neoswga.cli_unified import COMMAND_GROUPS, create_parser

    parser = create_parser()
    assert "evaluate-set" in parser._subparsers._group_actions[0].choices
    assert "evaluate-set" in {c for _, cmds in COMMAND_GROUPS for c in cmds}


def test_expand_primers_accepts_a_candidates_file():
    """step3_df.csv must be a fallback, not a hard precondition."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    args = parser.parse_args(
        [
            "expand-primers",
            "-j",
            "p.json",
            "--fixed-primers",
            "ATCG",
            "-o",
            "out",
            "--candidates-file",
            "cands.csv",
        ]
    )
    assert args.candidates_file == "cands.csv"


def test_analyze_set_no_longer_demands_unused_arguments():
    """--fg/--fg-kmers were required and never read; --output wrote nothing."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    args = parser.parse_args(["analyze-set", "--primers", "ATCGATCGATCG"])
    assert args.fg is None and args.fg_kmers is None and args.output is None


def test_coverage_commands_use_realistic_reach_not_processivity():
    """contract-set and rescore-set both scored coverage at 70 kb."""
    import inspect

    from neoswga.cli import iterate

    src = inspect.getsource(iterate)
    assert "70000" not in src and "70_000" not in src, (
        "a hardcoded 70 kb reach is back; coverage must use "
        "polymerase_extension_reach(..., coverage_metric='realistic')"
    )
    assert "polymerase_extension_reach" in src
