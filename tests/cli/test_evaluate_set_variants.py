"""`evaluate-set --variants FILE`: which sites survive in which strain.

The flag answers a question no command could ask before: given a panel and a
table of variants against the reference it was designed on, how much of that
panel still binds in each strain. It decides nothing -- no limit, no schema
key, no selection stage reads it -- so what these tests check is that the
figures reach the JSON, that they are labelled with what the route cannot see,
and that a table the reader refuses stops the command without leaving a
failure record beside somebody's finished design.
"""

import json
import os
import pathlib
import subprocess
import sys

import pytest

from neoswga.core.exceptions import ReferenceDataError
from tests.variant_invariants import assert_site_counts_add_up

ALT = {"A": "C", "C": "G", "G": "T", "T": "A"}


def _evaluate(run_cli, genome_fasta, out, primers, extra=()):
    argv = ["evaluate-set", "--genome", genome_fasta, "-o", str(out), "--primers", *primers]
    return _checked(run_cli(argv + list(extra)))


def _checked(result):
    """The site-count invariants, asserted on every result that carries strains."""
    if "per_strain" in result:
        assert_site_counts_add_up({**result["variant_route"], "per_strain": result["per_strain"]})
    return result


def _sequence(genome_fasta):
    with open(genome_fasta) as handle:
        lines = [line.strip() for line in handle if not line.startswith(">")]
    return "".join(lines)


def _tsv(path, sequence, positions, strains=("strain_a",), carried=None):
    """A chrom/pos/ref/alt table, 1-based, with one 0/1 column per strain."""
    carried = carried if carried is not None else {name: positions for name in strains}
    header = ["chrom", "pos", "ref", "alt", *strains]
    rows = []
    for position in sorted(positions):
        base = sequence[position - 1]
        row = ["target", str(position), base, ALT[base]]
        row += ["1" if position in carried.get(name, ()) else "0" for name in strains]
        rows.append("\t".join(row))
    path.write_text("\n".join(["\t".join(header), *rows]) + "\n")
    return str(path)


def _planted(tmp_path, genome_fasta, planted_primers):
    """A table that destroys the first planted site of the first primer.

    One known site, so the counts are checkable by hand rather than only
    self-consistent.
    """
    sequence = _sequence(genome_fasta)
    first = sequence.find(planted_primers[0])
    assert first >= 0
    hit = first + 1  # 1-based, inside the site
    return _tsv(tmp_path / "variants.tsv", sequence, [hit]), sequence, hit


# ----------------------------------------------------------------------
# The block reaches the JSON
# ----------------------------------------------------------------------


def test_the_per_strain_block_is_written(run_cli, genome_fasta, planted_primers, tmp_path):
    variants, _text, _hit = _planted(tmp_path, genome_fasta, planted_primers)
    out = tmp_path / "eval"
    result = _evaluate(run_cli, genome_fasta, out, planted_primers, ["--variants", variants])

    assert "per_strain" in result
    record = result["per_strain"]["strain_a"]
    assert record["status"] == "measured"
    assert record["affected_sites"]["value"] == 1.0
    assert record["intact_sites"]["value"] == result["variant_route"]["reference_sites"] - 1
    on_disk = json.loads((out / "evaluation.json").read_text())
    assert on_disk["per_strain"] == json.loads(json.dumps(result["per_strain"]))


def test_the_route_states_its_three_limits_and_its_window(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    variants, _text, _hit = _planted(tmp_path, genome_fasta, planted_primers)
    result = _evaluate(
        run_cli, genome_fasta, tmp_path / "eval", planted_primers, ["--variants", variants]
    )

    route = result["variant_route"]
    assert len(route["limits"]) == 3
    assert route["three_prime_window_nt"] == 5
    assert "reporting split" in route["three_prime_window_basis"]
    joined = " ".join(route["limits"]).lower()
    assert "gained" in joined and "indel" in joined and "absent from the reference" in joined


def test_each_primer_carries_its_intact_fraction(run_cli, genome_fasta, planted_primers, tmp_path):
    variants, _text, _hit = _planted(tmp_path, genome_fasta, planted_primers)
    result = _evaluate(
        run_cli, genome_fasta, tmp_path / "eval", planted_primers, ["--variants", variants]
    )

    per_primer = {record["primer"]: record for record in result["primers"]}
    hit = per_primer[planted_primers[0]]
    untouched = per_primer[planted_primers[1]]
    assert 0.0 < hit["intact_fraction_every_strain"] < 1.0
    assert untouched["intact_fraction_every_strain"] == 1.0
    assert set(hit["intact_fraction_by_strain"]) == {"strain_a"}


def test_without_the_flag_nothing_changes(run_cli, genome_fasta, planted_primers, tmp_path):
    """The flag is the only thing that adds these keys; an absent flag leaves
    the output it had."""
    result = _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers)
    assert "per_strain" not in result
    assert "variant_route" not in result
    assert all("intact_fraction_every_strain" not in record for record in result["primers"])


def test_a_strain_carrying_nothing_keeps_every_site(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    sequence = _sequence(genome_fasta)
    first = sequence.find(planted_primers[0]) + 1
    variants = _tsv(
        tmp_path / "two_strains.tsv",
        sequence,
        [first],
        strains=("carrier", "clean"),
        carried={"carrier": [first], "clean": []},
    )
    result = _evaluate(
        run_cli, genome_fasta, tmp_path / "eval", planted_primers, ["--variants", variants]
    )

    strains = result["per_strain"]
    assert strains["clean"]["intact_site_fraction"]["value"] == 1.0
    assert strains["carrier"]["intact_site_fraction"]["value"] < 1.0
    # The worst case is the carrier's, and it is a measurement because both
    # strains were measured.
    worst = result["variant_route"]["worst_strain_intact_fraction"]
    assert worst["value"] == strains["carrier"]["intact_site_fraction"]["value"]


# ----------------------------------------------------------------------
# Refusals, and what they must not leave behind
# ----------------------------------------------------------------------


def test_a_table_against_another_assembly_stops_the_command(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    sequence = _sequence(genome_fasta)
    path = tmp_path / "wrong.tsv"
    base = sequence[99]
    path.write_text(f"chrom\tpos\tref\talt\ntarget\t100\t{ALT[base]}\t{base}\n")

    with pytest.raises(ReferenceDataError):
        _evaluate(
            run_cli, genome_fasta, tmp_path / "eval", planted_primers, ["--variants", str(path)]
        )


def test_a_refusal_writes_no_failure_record(run_cli, genome_fasta, planted_primers, tmp_path):
    """`evaluate-set` describes an existing set and designs nothing, so a
    refusal of ITS input says nothing about the design whose directory was
    named. A `design_failure.json` left there would make `export` refuse a
    panel no run had failed on."""
    from neoswga.cli._failure import REPORT_ONLY_COMMANDS, write_failure_artifact

    assert "evaluate-set" in REPORT_ONLY_COMMANDS

    out = tmp_path / "eval"
    out.mkdir()

    class _Args:
        command = "evaluate-set"
        json_file = None
        output = str(out)

    written = write_failure_artifact(_Args(), ReferenceDataError("t", "r"))
    assert written is None
    assert not os.path.exists(out / "design_failure.json")


def test_several_foreground_genomes_are_refused_rather_than_guessed(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """A variant table is stated against one assembly. Picking one of several
    targets would report a strain's figures against the wrong coordinates."""
    second = tmp_path / "second.fasta"
    second.write_text(pathlib.Path(genome_fasta).read_text().replace(">target", ">target2"))
    sequence = _sequence(genome_fasta)
    variants = _tsv(tmp_path / "v.tsv", sequence, [sequence.find(planted_primers[0]) + 1])

    argv = [
        "evaluate-set",
        "--genome",
        genome_fasta,
        str(second),
        "-o",
        str(tmp_path / "eval"),
        "--primers",
        *planted_primers,
        "--variants",
        variants,
    ]
    with pytest.raises(ReferenceDataError, match="one assembly"):
        run_cli(argv)


def test_several_genomes_on_the_command_line_work_once_the_reference_is_named(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    """Follows the refusal above: the advice names `--variants-reference`, and
    adding it to the same command succeeds against that genome."""
    second = tmp_path / "second.fasta"
    second.write_text(pathlib.Path(genome_fasta).read_text().replace(">target", ">target2"))
    sequence = _sequence(genome_fasta)
    variants = _tsv(tmp_path / "v.tsv", sequence, [sequence.find(planted_primers[0]) + 1])
    argv = [
        "evaluate-set",
        "--genome",
        genome_fasta,
        str(second),
        "-o",
        str(tmp_path / "eval"),
        "--primers",
        *planted_primers,
        "--variants",
        variants,
    ]
    with pytest.raises(ReferenceDataError) as refused:
        run_cli(argv)
    assert "--variants-reference" in refused.value.remediation

    result = _checked(run_cli(argv + ["--variants-reference", genome_fasta]))
    assert result["variant_route"]["reference_fasta"] == genome_fasta
    assert result["per_strain"]["strain_a"]["affected_sites"]["value"] == 1.0


def test_a_variants_reference_without_a_table_is_refused(
    run_cli, genome_fasta, planted_primers, tmp_path
):
    with pytest.raises(ValueError, match="no --variants"):
        _evaluate(
            run_cli,
            genome_fasta,
            tmp_path / "eval",
            planted_primers,
            ["--variants-reference", genome_fasta],
        )


def test_the_printed_report_names_the_reference(
    run_cli, genome_fasta, planted_primers, tmp_path, capsys
):
    variants, _text, _hit = _planted(tmp_path, genome_fasta, planted_primers)
    _evaluate(run_cli, genome_fasta, tmp_path / "eval", planted_primers, ["--variants", variants])
    printed = capsys.readouterr().out
    assert f"Variants placed against {genome_fasta}" in printed


# ----------------------------------------------------------------------
# With -j: the refusal's advice, followed, has to work
# ----------------------------------------------------------------------

ROOT = pathlib.Path(__file__).resolve().parents[2]


def _run(args, cwd):
    """The CLI in a child process, against this checkout.

    A child rather than `run_cli`: `-j` loads params.json into the `parameter`
    module globals, which would outlive this test in the session.
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(filter(None, [str(ROOT), env.get("PYTHONPATH", "")]))
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=300,
        env=env,
    )


@pytest.fixture
def two_target_params(tmp_path, genome_fasta, planted_primers):
    sequence = _sequence(genome_fasta)
    second = tmp_path / "second.fasta"
    second.write_text(pathlib.Path(genome_fasta).read_text().replace(">target", ">target2"))
    data = tmp_path / "data"
    data.mkdir()
    params = {
        "fg_genomes": [genome_fasta, str(second)],
        "fg_prefixes": [str(data / "target"), str(data / "second")],
        "fg_seq_lengths": [len(sequence), len(sequence)],
        "data_dir": str(data),
        "polymerase": "phi29",
    }
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))
    variants = _tsv(tmp_path / "v.tsv", sequence, [sequence.find(planted_primers[0]) + 1])
    return {"params": str(path), "data": data, "second": str(second), "variants": variants}


def test_two_configured_genomes_are_refused_and_the_advice_works(
    two_target_params, genome_fasta, planted_primers, tmp_path
):
    base = [
        "evaluate-set",
        "-j",
        two_target_params["params"],
        "--primers",
        *planted_primers,
        "--linear",
        "--variants",
        two_target_params["variants"],
    ]
    refused = _run(base + ["-o", "refused"], cwd=str(tmp_path))
    assert refused.returncode != 0
    assert "--variants-reference" in refused.stderr
    assert genome_fasta in refused.stderr
    # A refusal of the evaluation leaves nothing beside the design.
    assert not (two_target_params["data"] / "design_failure.json").exists()

    followed = _run(
        base + ["--variants-reference", genome_fasta, "-o", "followed"], cwd=str(tmp_path)
    )
    assert followed.returncode == 0, followed.stderr
    result = json.loads((tmp_path / "followed" / "evaluation.json").read_text())
    _checked(result)
    assert result["variant_route"]["reference_fasta"] == genome_fasta
    assert f"Variants placed against {genome_fasta}" in followed.stdout


def test_one_readable_genome_of_two_is_not_chosen_silently(
    two_target_params, genome_fasta, planted_primers, tmp_path
):
    os.remove(two_target_params["second"])
    refused = _run(
        [
            "evaluate-set",
            "-j",
            two_target_params["params"],
            "--primers",
            *planted_primers,
            "--linear",
            "--variants",
            two_target_params["variants"],
            "-o",
            "out",
        ],
        cwd=str(tmp_path),
    )
    # Through the CLI, an absent configured genome is refused before the variant
    # route is reached; the unit test below pins the route's own rule.
    assert refused.returncode != 0
    assert not (tmp_path / "out" / "evaluation.json").exists()


def test_the_reference_rule_counts_configured_genomes_not_readable_ones(tmp_path, genome_fasta):
    """The rule itself, without the earlier check that happens to stand in
    front of it on the CLI: two configured, one readable, nothing named, is a
    refusal -- not the readable one."""
    from neoswga.cli.evaluate import _variant_reference

    missing = str(tmp_path / "moved.fasta")
    with pytest.raises(ReferenceDataError, match="2 foreground genomes"):
        _variant_reference(["a", "b"], [genome_fasta, missing], [100, 100])
    assert _variant_reference(["a", "b"], [genome_fasta, missing], [100, 100], genome_fasta) == (
        "a",
        genome_fasta,
        100,
    )
    with pytest.raises(ReferenceDataError, match="cannot be read"):
        _variant_reference(["a", "b"], [genome_fasta, missing], [100, 100], missing)


def test_a_genome_flag_beside_params_gets_advice_that_does_not_loop(
    two_target_params, genome_fasta, planted_primers, tmp_path
):
    """`-j` with two targets plus `--genome A` pairs two prefixes with one
    FASTA. The refusal must not tell the user to pass `--genome`, which they
    did; it tells them to leave it off, and doing so with the reference named
    succeeds."""
    common = [
        "evaluate-set",
        "-j",
        two_target_params["params"],
        "--primers",
        *planted_primers,
        "--linear",
        "--variants",
        two_target_params["variants"],
    ]
    refused = _run(common + ["--genome", genome_fasta, "-o", "mixed"], cwd=str(tmp_path))
    assert refused.returncode != 0
    assert "leave --genome off" in refused.stderr

    followed = _run(
        common + ["--variants-reference", genome_fasta, "-o", "fixed"], cwd=str(tmp_path)
    )
    assert followed.returncode == 0, followed.stderr
