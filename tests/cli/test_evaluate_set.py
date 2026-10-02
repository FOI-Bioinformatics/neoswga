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


# ----------------------------------------------------------------------
# Hosts given by path, and the per-reference blocks
# ----------------------------------------------------------------------
#
# Before 2026-10-02 no command could score a set against a host that was not in
# the design's params.json: the background came only from `parameter.bg_prefixes`
# / `bg_genomes`. These pin the new route and, as importantly, that the fields
# that were already there did not move.


#: Every field `evaluation.json` carried before the per-reference blocks were
#: added. A caller reads these by name, so the list is checked rather than
#: trusted: an addition must not rename or repurpose one of them.
LEGACY_FIELDS = (
    "primers",
    "num_primers",
    "fg_coverage",
    "per_target_coverage",
    "extension_reach_bp",
    "total_binding_sites",
    "total_background_sites",
    "selectivity_ratio",
    "occupancy_selectivity_ratio",
    "mean_binding_distance_bp",
    "genome_bp",
    "mean_gap_bp",
    "max_gap_bp",
    "gap_gini",
    "gap_interpretation",
    "zero_site_primers",
    "unresolved_primers",
    "_regime_note",
    "conditions",
)


@pytest.fixture
def host_fasta(tmp_path, planted_primers, genome_seq):
    """A host carrying one of the planted primers, and 8 kb of its own sequence."""
    import random

    rng = random.Random(4242)
    seq = list("".join(rng.choice("ACGT") for _ in range(8_000)))
    host_primer = planted_primers[0]
    for offset in (500, 3_000, 6_500):
        seq[offset : offset + len(host_primer)] = list(host_primer)
    sequence = "".join(seq)
    path = tmp_path / "host.fasta"
    path.write_text(
        ">host\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    return str(path)


def test_every_legacy_field_keeps_its_name(run_cli, genome_fasta, planted_primers, tmp_path):
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)
    missing = [name for name in LEGACY_FIELDS if name not in result]
    assert not missing, f"fields disappeared from evaluation.json: {missing}"


def test_the_new_blocks_do_not_change_a_run_without_them(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """With no host anywhere, every legacy value is what it was.

    The per-reference blocks are additive by construction: the pooled fields are
    computed by the code that always computed them, and this is the check that
    the construction holds. A background that was never configured still reads
    None rather than 0 -- an absent measurement, not perfect specificity.
    """
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)

    assert result["total_background_sites"] is None
    assert result["selectivity_ratio"] is None
    assert result["occupancy_selectivity_ratio"] is None
    assert result["per_host"] == {}
    assert result["target_host_pairs"] == []
    # One target, so the worst over the targets is that target's coverage.
    assert result["worst_target_coverage"]["value"] == pytest.approx(
        result["per_target_coverage"][next(iter(result["per_target_coverage"]))], abs=1e-4
    )
    # No host, so the worst host density is unknown rather than unbounded.
    assert result["worst_host_selectivity_density"]["value"] is None
    assert result["worst_host_selectivity_density"]["unavailable"]


def test_a_host_by_path_is_measured_with_no_params_json(
    run_cli, genome_fasta, host_fasta, planted_primers, tmp_path
):
    """--background: a host nobody configured, scored anyway."""
    result = _evaluate(
        run_cli,
        genome_fasta,
        tmp_path / "eval",
        planted_primers,
        extra=["--background", host_fasta],
    )

    hosts = result["per_host"]
    assert len(hosts) == 1
    host = next(iter(hosts.values()))
    assert host["status"] == "measured"
    assert host["sites"]["value"] == 3  # one primer, three planted sites
    assert host["site_source"] == "counts"
    # Counted, not located: coverage says so instead of reading zero.
    assert host["coverage"]["value"] is None
    assert "counted rather than located" in host["coverage"]["unavailable"]

    # The pooled fields now have a measurement where they had none, and the
    # per-host rows sit beside them.
    assert result["total_background_sites"] == 3
    assert result["selectivity_ratio"] is not None


def test_scanning_a_host_by_path_gives_the_same_sites_and_adds_coverage(
    run_cli, genome_fasta, host_fasta, planted_primers, tmp_path
):
    counted = _evaluate(
        run_cli,
        genome_fasta,
        tmp_path / "counted",
        planted_primers,
        extra=["--background", host_fasta],
    )
    scanned = _evaluate(
        run_cli,
        genome_fasta,
        tmp_path / "scanned",
        planted_primers,
        extra=["--background", host_fasta, "--scan-background"],
    )

    counted_host = next(iter(counted["per_host"].values()))
    scanned_host = next(iter(scanned["per_host"].values()))
    assert scanned_host["sites"]["value"] == counted_host["sites"]["value"]
    assert scanned_host["site_source"] == "positions"
    assert scanned_host["coverage"]["value"] is not None
    assert counted["total_background_sites"] == scanned["total_background_sites"]


def test_the_pairwise_table_reports_both_figures_and_names_the_comparable_one(
    run_cli, genome_fasta, host_fasta, planted_primers, tmp_path
):
    result = _evaluate(
        run_cli,
        genome_fasta,
        tmp_path / "eval",
        planted_primers,
        extra=["--background", host_fasta],
    )
    pairs = result["target_host_pairs"]
    assert len(pairs) == 1
    pair = pairs[0]
    assert pair["selectivity_density"]["value"] is not None
    assert pair["selectivity_ratio"]["value"] is not None
    # Known Issue 6: only one of the two survives a change of host size.
    assert pair["comparable"] == "selectivity_density"
    assert result["worst_host_selectivity_density"]["value"] == pytest.approx(
        pair["selectivity_density"]["value"]
    )


def test_an_unreadable_host_leaves_the_reduction_unknown(
    run_cli, genome_fasta, host_fasta, planted_primers, tmp_path
):
    """Two hosts, one of them empty of sequence: the worst case is unknown.

    The measurable host has a density, so a reduction that skipped the other
    would return a number that reads as the panel's worst case against its hosts.
    """
    broken = tmp_path / "broken.fasta"
    broken.write_text(">broken\n\n")

    result = _evaluate(
        run_cli,
        genome_fasta,
        tmp_path / "eval",
        planted_primers,
        extra=["--background", host_fasta, str(broken)],
    )

    densities = [p["selectivity_density"]["value"] for p in result["target_host_pairs"]]
    assert None in densities, "the premise needs a host that could not be measured"
    assert [d for d in densities if d is not None], "and one that could"
    assert result["worst_host_selectivity_density"]["value"] is None
    assert result["worst_host_selectivity_density"]["unavailable"]
    # The pooled count is withheld for the same reason: a sum over the hosts that
    # answered is not this panel's host binding.
    assert result["total_background_sites"] is None


# ----------------------------------------------------------------------
# Reading the delivered set
# ----------------------------------------------------------------------


def _write_results(directory, sets):
    """A minimal step4_improved_df.csv holding several alternative sets."""
    import csv

    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "step4_improved_df.csv"
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["primer", "set_index"])
        for index, primers in enumerate(sets):
            for primer in primers:
                writer.writerow([primer, index])
    return directory


def test_from_results_reads_one_set_not_the_union(
    run_cli, genome_fasta, planted_primers, absent_primer, tmp_path
):
    """Alternatives are separate answers; pooling them describes no real tube."""
    results = _write_results(tmp_path / "run", [planted_primers, [absent_primer]])

    result = run_cli(
        [
            "evaluate-set",
            "--genome",
            genome_fasta,
            "-o",
            str(tmp_path / "eval"),
            "--from-results",
            str(results),
        ]
    )
    assert [p["primer"] for p in result["primers"]] == planted_primers
    assert absent_primer not in [p["primer"] for p in result["primers"]]


def test_from_results_honours_the_set_index(
    run_cli, genome_fasta, planted_primers, absent_primer, tmp_path
):
    results = _write_results(tmp_path / "run", [planted_primers, [absent_primer]])

    result = run_cli(
        [
            "evaluate-set",
            "--genome",
            genome_fasta,
            "-o",
            str(tmp_path / "eval"),
            "--from-results",
            str(results),
            "--set",
            "1",
        ]
    )
    assert [p["primer"] for p in result["primers"]] == [absent_primer]


def test_a_set_index_the_file_does_not_hold_is_refused(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """An empty panel and an absent one read identically downstream."""
    from neoswga.core.exceptions import ReferenceDataError

    results = _write_results(tmp_path / "run", [planted_primers])
    with pytest.raises(ReferenceDataError):
        run_cli(
            [
                "evaluate-set",
                "--genome",
                genome_fasta,
                "-o",
                str(tmp_path / "eval"),
                "--from-results",
                str(results),
                "--set",
                "7",
            ]
        )


def test_from_results_and_primers_together_are_refused(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    results = _write_results(tmp_path / "run", [planted_primers])
    with pytest.raises(ValueError, match="not both"):
        run_cli(
            [
                "evaluate-set",
                "--genome",
                genome_fasta,
                "-o",
                str(tmp_path / "eval"),
                "--from-results",
                str(results),
                "--primers",
                *planted_primers,
            ]
        )
