"""Execution tests for `neoswga evaluate-set`.

This command exists to fix a specific silent failure: a primer set designed
elsewhere is absent from the HDF5 index, so every coverage number computed for
it is zero -- not because the primers are bad, but because nothing looked them
up. Before, that was indistinguishable from a genuinely useless set.

So the tests that matter here are the ones that check the *distinction* is
preserved all the way into the JSON, not just in a log line. A caller reading
`evaluation.json` has to be able to tell "this primer binds nowhere" from "we
never found out".
"""

import json
import os

import pytest


def _evaluate(run_cli, genome_fasta, out, primers, extra=()):
    argv = ["evaluate-set", "--genome", genome_fasta, "-o", str(out), "--primers", *primers]
    return run_cli(argv + list(extra))


# ----------------------------------------------------------------------
# It runs at all, without a prior count-kmers
# ----------------------------------------------------------------------


def test_runs_on_a_bare_fasta_with_no_prior_pipeline(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """The whole point of --genome: no init, no count-kmers, no params.json."""
    out = tmp_path / "eval"
    result = _evaluate(run_cli, genome_fasta, out, planted_primers)

    assert os.path.exists(out / "evaluation.json")
    assert result["num_primers"] == len(planted_primers)
    assert result["genome_bp"] == 20_000


def test_finds_the_planted_binding_sites(run_cli, genome_fasta, planted_primers, tmp_path):
    """Three primers planted at three positions each -- nine sites."""
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)

    per_primer = {p["primer"]: p for p in result["primers"]}
    for primer in planted_primers:
        assert per_primer[primer]["binding_sites"] >= 3, primer
    assert result["total_binding_sites"] >= 9


def test_written_json_matches_the_returned_result(run_cli, genome_fasta, planted_primers, tmp_path):
    out = tmp_path / "eval"
    result = _evaluate(run_cli, genome_fasta, out, planted_primers)
    on_disk = json.loads((out / "evaluation.json").read_text())
    assert on_disk == json.loads(json.dumps(result))


def test_output_is_valid_json_not_just_python_readable(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """Same RFC 8259 concern as the optimizer summary: no Infinity/NaN tokens."""
    out = tmp_path / "eval"
    _evaluate(run_cli, genome_fasta, out, planted_primers)

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    json.loads((out / "evaluation.json").read_text(), parse_constant=reject)


# ----------------------------------------------------------------------
# The distinction the command exists to make
# ----------------------------------------------------------------------


def test_a_primer_that_binds_nowhere_is_reported_as_such(
    run_cli, genome_fasta, planted_primers, absent_primer, tmp_path
):
    """Zero sites after a real scan is a real answer, and must be labelled."""
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers + [absent_primer])

    per_primer = {p["primer"]: p for p in result["primers"]}
    assert per_primer[absent_primer]["binding_sites"] == 0
    assert absent_primer in result["zero_site_primers"]

    # It was looked up successfully -- the answer is "binds nowhere", not
    # "unknown". Conflating the two is the bug this command was built to fix.
    assert per_primer[absent_primer]["status"] == "no_sites"


def test_primers_that_do_bind_are_not_flagged(
    run_cli, genome_fasta, planted_primers, absent_primer, tmp_path
):
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers + [absent_primer])
    per_primer = {p["primer"]: p for p in result["primers"]}
    for primer in planted_primers:
        assert per_primer[primer]["status"] == "ok", primer
        assert primer not in result["zero_site_primers"]


def test_unresolved_primers_are_distinguished_from_absent_ones(
    tmp_path, planted_primers, absent_primer
):
    """Without a FASTA to scan, zero coverage is meaningless -- say so in the JSON.

    This is the failure mode the command was written for, and it was still
    reachable through the command's own output: with no genome to scan,
    `zero_site_primers` is empty and every primer was reported as fine, while
    the coverage numbers were computed from nothing.
    """
    from neoswga.core.position_cache import PositionCache

    prefix = str(tmp_path / "nonexistent_index")
    cache = PositionCache(
        [prefix], planted_primers + [absent_primer], genome_paths=None, on_missing="warn"
    )

    # Nothing could be resolved, so nothing is a confirmed zero...
    assert cache.zero_site_primers == []
    # ...but they are all recorded as unresolved.
    unresolved = {p for _, p in cache.missing_primers}
    assert unresolved == set(planted_primers) | {absent_primer}


# ----------------------------------------------------------------------
# Numbers that feed the published-benchmark comparison
# ----------------------------------------------------------------------


def test_gap_statistics_are_finite_and_ordered(run_cli, genome_fasta, planted_primers, tmp_path):
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)

    assert result["mean_gap_bp"] > 0
    assert result["max_gap_bp"] >= result["mean_gap_bp"]
    assert 0.0 <= result["gap_gini"] <= 1.0


def test_gap_interpretation_is_present_and_measured(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """With real sites planted, the verdicts must not read 'unknown'."""
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)

    rows = result["gap_interpretation"]
    assert len(rows) == 3
    assert [r["verdict"] for r in rows] != ["unknown"] * 3


def test_coverage_is_bounded_and_uses_the_realistic_reach(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """phi29's realistic per-primer reach is ~3 kb, not the 70 kb processivity."""
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)

    assert 0.0 < result["fg_coverage"] <= 1.0
    assert result["extension_reach_bp"] == 3000


def test_linear_flag_drops_the_wraparound_gap(run_cli, genome_fasta, planted_primers, tmp_path):
    """On a linear target the origin-spanning gap is not a real gap."""
    circular = _evaluate(run_cli, genome_fasta, tmp_path / "circ", planted_primers)
    linear = _evaluate(run_cli, genome_fasta, tmp_path / "lin", planted_primers, extra=["--linear"])
    assert linear["max_gap_bp"] <= circular["max_gap_bp"]


# ----------------------------------------------------------------------
# Argument handling
# ----------------------------------------------------------------------


def test_primers_file_is_equivalent_to_primers(run_cli, genome_fasta, planted_primers, tmp_path):
    listing = tmp_path / "primers.txt"
    listing.write_text("\n".join(planted_primers) + "\n")

    from_args = _evaluate(run_cli, genome_fasta, tmp_path / "a", planted_primers)
    from_file = run_cli(
        [
            "evaluate-set",
            "--genome",
            genome_fasta,
            "-o",
            str(tmp_path / "b"),
            "--primers-file",
            str(listing),
        ]
    )
    assert from_file["total_binding_sites"] == from_args["total_binding_sites"]
    assert from_file["fg_coverage"] == from_args["fg_coverage"]


def test_no_primers_exits_with_a_message(run_cli, genome_fasta, tmp_path, caplog):
    """`collect_primers_from_args` exits rather than raising.

    Worth pinning explicitly: it means a programmatic caller of
    `run_evaluate_set` gets SystemExit, not an exception it can catch by type.
    That is the shared helper's convention across every command, so the command
    is consistent with the rest of the CLI; the guard here is that the exit is
    accompanied by an explanation rather than being silent.
    """
    import logging

    # Attach our own handler rather than relying on caplog: logging config is
    # process-global and other tests reconfigure propagation, which made this
    # pass alone and fail in the full suite.
    records = []

    class _Capture(logging.Handler):
        def emit(self, record):
            records.append(record.getMessage())

    target = logging.getLogger("neoswga.cli._common")
    handler = _Capture()
    target.addHandler(handler)
    try:
        with pytest.raises(SystemExit) as exc:
            run_cli(["evaluate-set", "--genome", genome_fasta, "-o", str(tmp_path / "e")])
    finally:
        target.removeHandler(handler)

    assert exc.value.code == 1
    assert any("No primer" in m for m in records), f"exit was not explained; captured: {records}"


def test_no_genome_and_no_params_is_a_clear_error(run_cli, tmp_path, planted_primers):
    with pytest.raises(ValueError, match="params.json|--genome"):
        run_cli(["evaluate-set", "-o", str(tmp_path / "e"), "--primers", *planted_primers])
