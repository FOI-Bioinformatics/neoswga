"""`improve-set` end to end, on an oligo set this pipeline never saw.

Phase 2b of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.

The core rules are pinned in `tests/test_set_improvement.py`. What only a real
invocation can show is pinned here:

- a mixed-length, externally designed set is evaluated and improved with no
  `count-kmers` run, no k-mer table and no position index anywhere;
- the figures predicted for the top proposal are the figures `evaluate-set`
  reports for the edited set, compared across two separate processes;
- `fixed_oligos` in params.json reaches the code that decides;
- the command reports and writes no primer set.
"""

from __future__ import annotations

import json
import os
import random
import subprocess
import sys
from pathlib import Path

import pytest

from neoswga.core.dimer import is_dimer_fast

ROOT = Path(__file__).resolve().parent.parent
LENGTH = 60_000


def _sequence(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _compatible_oligos(lengths, seed=0):
    """Oligos of the given lengths, no pair of which dimerises at the default 3 bp.

    Chosen rather than assumed, so the dimer screen is not what decides the
    tests below.
    """
    rng = random.Random(seed)
    chosen = []
    for length in lengths:
        while True:
            oligo = "".join(rng.choice("ACGT") for _ in range(length))
            if all(
                not is_dimer_fast(oligo, other, max_dimer_bp=3)
                and not is_dimer_fast(other, oligo, max_dimer_bp=3)
                for other in chosen
            ):
                chosen.append(oligo)
                break
    return chosen


def _plant(sequence, plants):
    seq = list(sequence)
    for primer, offsets in plants.items():
        for offset in offsets:
            seq[offset : offset + len(primer)] = list(primer)
    return "".join(seq)


def _write_fasta(path, sequence, name):
    path.write_text(
        f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    return str(path)


def _run(args, cwd):
    """Run the CLI from `cwd`, against the checkout these tests were read from.

    The child runs in a temporary directory, where `neoswga` would otherwise
    resolve to whichever checkout is installed. In a git worktree that is a
    different tree from the one under test.
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
def workdir(tmp_path):
    """Two 60 kb targets and a host, as bare FASTA files. Nothing is indexed.

    The set is of mixed length: a 12-mer and a 10-mer that bind, and a 14-mer
    that binds nothing. Target B is the worse covered. `fills_b` binds only B,
    where the set does not reach; `helps_a` binds only A and would raise the
    pooled coverage more.
    """
    twelve, ten, fourteen, fills_b, helps_a = _compatible_oligos([12, 10, 14, 12, 12], seed=7)
    a = _plant(
        _sequence(LENGTH, seed=1),
        {
            twelve: range(1_000, 30_000, 4_000),
            ten: range(3_000, 30_000, 4_000),
            helps_a: range(31_000, 58_000, 2_000),
        },
    )
    b = _plant(
        _sequence(LENGTH, seed=2),
        {
            twelve: range(1_000, 12_000, 4_000),
            ten: range(3_000, 12_000, 4_000),
            fills_b: range(20_000, 32_000, 2_000),
        },
    )
    host = _plant(_sequence(90_000, seed=3), {twelve: (5_000, 40_000), fills_b: (70_000,)})
    for sequence in (a, b, host):
        assert fourteen not in sequence and _reverse_complement(fourteen) not in sequence
    assert fills_b not in a and _reverse_complement(fills_b) not in a
    assert helps_a not in b and _reverse_complement(helps_a) not in b

    data = tmp_path / "data"
    data.mkdir()
    (data / "step3_df.csv").write_text(f"primer\n{helps_a}\n{fills_b}\n")
    genomes = [
        _write_fasta(tmp_path / "target_a.fasta", a, "a"),
        _write_fasta(tmp_path / "target_b.fasta", b, "b"),
    ]
    params = {
        "fg_genomes": genomes,
        "fg_prefixes": [str(data / "target_a"), str(data / "target_b")],
        "fg_seq_lengths": [LENGTH, LENGTH],
        "data_dir": str(data),
        "polymerase": "phi29",
        # The design's own range. The 14-mer in the set is outside it, as an
        # oligo designed elsewhere may be; it is evaluated all the same.
        "min_k": 10,
        "max_k": 12,
        "min_tm": 0,
        "max_tm": 90,
    }
    (tmp_path / "params.json").write_text(json.dumps(params))
    return {
        "dir": tmp_path,
        "params": params,
        "genomes": genomes,
        "host": _write_fasta(tmp_path / "host.fasta", host, "host"),
        "set": [twelve, ten, fourteen],
        "nowhere": fourteen,
        "fills_b": fills_b,
        "helps_a": helps_a,
    }


def _report(directory):
    return json.loads((directory / "improvement_report.json").read_text())


def _shown(report, kind):
    return report["sections"][kind]["proposals"]


def _every_proposal(report):
    return [p for section in report["sections"].values() for p in section["proposals"]]


def _no_primer_set_was_written(root):
    return not list(root.rglob("step4_improved_df.csv"))


def test_an_outside_mixed_length_set_is_diagnosed_from_fasta_alone(workdir):
    """No params.json, no count-kmers: the diagnosis and the drop still come."""
    proc = _run(
        [
            "improve-set",
            "--primers",
            *workdir["set"],
            "--genome",
            *workdir["genomes"],
            "--background",
            workdir["host"],
            "--linear",
            "-o",
            str(workdir["dir"] / "out"),
        ],
        cwd=str(workdir["dir"]),
    )
    assert proc.returncode == 0, proc.stderr

    report = _report(workdir["dir"] / "out")
    assert sorted(len(p) for p in report["primers"]) == [10, 12, 14]
    assert report["settings"]["source"] == "defaults (no params.json)"
    assert len(report["current"]["per_target"]) == 2
    (host,) = report["current"]["per_host"].values()
    assert host["sites"]["value"] == 2
    assert host["coverage"]["value"] is None, "a counted host has no coverage, not zero"

    assert [p["drop"] for p in _shown(report, "drop")] == [[workdir["nowhere"]]]
    assert report["candidate_pool"]["unavailable"]
    assert all(not p["add"] for p in _every_proposal(report))

    assert "SET IMPROVEMENT" in proc.stdout
    assert "Nothing was applied" in proc.stdout
    assert _no_primer_set_was_written(workdir["dir"])


def test_the_predicted_figures_are_what_evaluate_set_then_reports(workdir):
    """Improved with no count-kmers run, and the prediction holds in a new process."""
    common = ["-j", "params.json", "--background", workdir["host"], "--linear"]
    proc = _run(
        ["improve-set", *common, "--primers", *workdir["set"], "-o", "shared"],
        cwd=str(workdir["dir"]),
    )
    assert proc.returncode == 0, proc.stderr
    report = _report(workdir["dir"] / "shared")

    assert report["settings"]["source"] == "params.json"
    top = _shown(report, "add")[0]
    assert top["rank"] == 1 and top["kind"] == "add"
    assert top["add"] == [workdir["fills_b"]], "the hole on the worst target, not the pooled gain"
    assert all(workdir["helps_a"] not in p["add"] for p in _every_proposal(report))
    assert "not an independent check" in report["prediction_note"]
    # With a pool, the oligo that binds nothing is still reported.
    assert [p["drop"] for p in _shown(report, "drop")] == [[workdir["nowhere"]]]

    # The same output directory, because a host given by path is named after
    # it, and the comparison below is of whole records including their names.
    again = _run(
        ["evaluate-set", *common, "--primers", *top["resulting_set"], "-o", "shared"],
        cwd=str(workdir["dir"]),
    )
    assert again.returncode == 0, again.stderr
    measured = json.loads((workdir["dir"] / "shared" / "evaluation.json").read_text())

    predicted = top["predicted"]
    for block in (
        "per_target",
        "per_host",
        "target_host_pairs",
        "worst_target_coverage",
        "worst_host_selectivity_density",
    ):
        assert measured[block] == predicted[block], block
    assert measured["extension_reach_bp"] == predicted["extension_reach_bp"]
    assert _no_primer_set_was_written(workdir["dir"])


def test_fixed_oligos_and_panel_limits_reach_the_command(workdir):
    """`fixed_oligos` protects an oligo; a configured limit is named, not enforced."""
    params = dict(workdir["params"], fixed_oligos=[workdir["nowhere"]], max_worst_hole=1000)
    (workdir["dir"] / "fixed.json").write_text(json.dumps(params))

    proc = _run(
        ["improve-set", "-j", "fixed.json", "--primers", *workdir["set"], "--linear", "-o", "out"],
        cwd=str(workdir["dir"]),
    )
    assert proc.returncode == 0, proc.stderr
    report = _report(workdir["dir"] / "out")

    assert report["settings"]["fixed_oligos"] == [workdir["nowhere"]]
    assert all(workdir["nowhere"] not in p["drop"] for p in _every_proposal(report))
    row = next(r for r in report["attribution"] if r["oligo"] == workdir["nowhere"])
    assert row["fixed"] is True

    top = _shown(report, "add")[0]
    assert top["add"] == [workdir["fills_b"]]
    limits = top["panel_limits"]
    assert limits["evaluated"] is True, limits
    assert limits["violations"] == ["worst hole above maximum"]
    assert limits["values"][0]["limit"] == "max_worst_hole"
    current = report["current_panel_limits"]
    assert current["violations"] == ["worst hole above maximum"]

    # The limits are judged on the evaluation's geometry, which is the flag's
    # and not params.json's `fg_circular` (absent here in both runs). Without
    # --linear the targets are circular, so the worst hole runs across the
    # origin and is a different figure.
    circular = _run(
        ["improve-set", "-j", "fixed.json", "--primers", *workdir["set"], "-o", "circular"],
        cwd=str(workdir["dir"]),
    )
    assert circular.returncode == 0, circular.stderr
    wrapped = _report(workdir["dir"] / "circular")["current_panel_limits"]
    assert wrapped["evaluated"] is True, wrapped
    assert wrapped["values"][0]["panel"] != current["values"][0]["panel"]


def test_the_limit_check_refuses_a_geometry_it_cannot_represent():
    """Circular targets, a linear configured host, and a limit on host coverage.

    The limit evaluator has one circular flag for every reference. Rather than
    judge the host on the targets' geometry, the limits are reported as not
    evaluated, with that as the reason.
    """
    from types import SimpleNamespace

    from neoswga.cli.iterate import _panel_limit_check
    from neoswga.core.pool_objective import PoolConstraints

    parameter = SimpleNamespace(
        fg_prefixes=["fg"],
        fg_seq_lengths=[1000],
        bg_prefixes=["bg"],
        bg_seq_lengths=[1000],
        fg_genomes=[],
        bg_genomes=[],
        bg_circular=False,
    )
    context = SimpleNamespace(constraints=PoolConstraints(max_host_coverage=0.5))

    check = _panel_limit_check(
        SimpleNamespace(linear=False, scan_background=False), parameter, context
    )
    outcome = check(["ACGTACGTACGT"])

    assert outcome["evaluated"] is False
    assert outcome["violations"] is None
    assert "one geometry" in outcome["unavailable"]


def test_a_refused_params_file_leaves_no_failure_record_beside_a_design(workdir):
    """A report-only command must not mark a finished design as failed.

    `export` refuses while `design_failure.json` is present, and only a design
    step clears it. So `improve-set` writes none when it refuses a params file,
    and it does not remove one a real failed design left either.
    """
    data = workdir["dir"] / "data"
    twelve, ten, _fourteen = workdir["set"]
    (data / "step4_improved_df.csv").write_text(
        f"primer,set_index,score\n{twelve},0,1.0\n{ten},0,1.0\n"
    )
    # phi29 has no recorded model for 14-mers, so the design request refuses.
    refused = dict(workdir["params"], max_k=14)
    (workdir["dir"] / "refused.json").write_text(json.dumps(refused))
    arguments = ["--from-results", str(data), "--linear", "-o", "out"]

    proc = _run(["improve-set", "-j", "refused.json", *arguments], cwd=str(workdir["dir"]))

    assert proc.returncode != 0
    assert "UnsupportedModelError" in proc.stderr, "the refusal is still printed"
    assert not list(workdir["dir"].rglob("design_failure.json"))
    assert not (workdir["dir"] / "out" / "improvement_report.json").exists()

    # A record left by a design that really failed is not this command's to touch.
    record = data / "design_failure.json"
    record.write_text(json.dumps({"stage": "optimize", "run_state": "failed"}))
    for params in ("refused.json", "params.json"):
        proc = _run(["improve-set", "-j", params, *arguments], cwd=str(workdir["dir"]))
        assert (proc.returncode == 0) == (params == "params.json"), proc.stderr
        assert json.loads(record.read_text()) == {"stage": "optimize", "run_state": "failed"}
    assert len(list(workdir["dir"].rglob("design_failure.json"))) == 1


def test_both_commands_measure_at_the_configured_coverage_reach(workdir):
    """`coverage_reach` in params.json reaches `evaluate-set` as it reaches this.

    The key was inert on `evaluate-set`: with 800 configured it reported 3000,
    so a prediction made at the design's reach was not what `evaluate-set` then
    measured.
    """
    (workdir["dir"] / "reach.json").write_text(
        json.dumps(dict(workdir["params"], coverage_reach=800))
    )
    common = ["-j", "reach.json", "--background", workdir["host"], "--linear"]

    proc = _run(
        ["improve-set", *common, "--primers", *workdir["set"], "-o", "shared"],
        cwd=str(workdir["dir"]),
    )
    assert proc.returncode == 0, proc.stderr
    report = _report(workdir["dir"] / "shared")
    assert report["settings"]["extension_reach_bp"] == 800
    assert report["current"]["extension_reach_bp"] == 800
    top = _shown(report, "add")[0]

    again = _run(
        ["evaluate-set", *common, "--primers", *top["resulting_set"], "-o", "shared"],
        cwd=str(workdir["dir"]),
    )
    assert again.returncode == 0, again.stderr
    measured = json.loads((workdir["dir"] / "shared" / "evaluation.json").read_text())

    assert measured["extension_reach_bp"] == 800
    assert measured["extension_reach_source"] == "coverage_reach in params.json"
    assert "800 bp reach" in again.stdout
    for block in ("per_target", "per_host", "target_host_pairs", "worst_target_coverage"):
        assert measured[block] == top["predicted"][block], block

    # Without the key nothing moved: the polymerase's reach, and it says so.
    plain = _run(
        ["evaluate-set", "-j", "params.json", "--primers", *workdir["set"], "--linear", "-o", "p"],
        cwd=str(workdir["dir"]),
    )
    assert plain.returncode == 0, plain.stderr
    unset = json.loads((workdir["dir"] / "p" / "evaluation.json").read_text())
    assert unset["extension_reach_bp"] == 3000
    assert unset["extension_reach_source"] == "polymerase default"
    assert unset["per_target"] != measured["per_target"], "the reach changes what is measured"


def test_max_edits_bounds_what_is_reported_and_says_what_was_considered(workdir):
    proc = _run(
        [
            "improve-set",
            "-j",
            "params.json",
            "--primers",
            *workdir["set"],
            "--linear",
            "--max-edits",
            "1",
            "-o",
            "out",
        ],
        cwd=str(workdir["dir"]),
    )
    assert proc.returncode == 0, proc.stderr
    report = _report(workdir["dir"] / "out")

    # Per section: one add and one drop, each saying what it was cut from.
    assert report["settings"]["max_edits_per_section"] == 1
    assert all(section["shown"] <= 1 for section in report["sections"].values())
    assert all(
        section["shown"] == min(1, section["considered"]) for section in report["sections"].values()
    )
    assert report["sections"]["add"]["shown"] == 1
    assert report["sections"]["drop"]["shown"] == 1
    assert "shown 1 of 1 considered" in proc.stdout


def test_a_delivered_set_is_read_through_from_results(workdir):
    """`--from-results DIR --set N` reads one set, as evaluate-set does."""
    results = workdir["dir"] / "results"
    results.mkdir()
    twelve, ten, fourteen = workdir["set"]
    (results / "step4_improved_df.csv").write_text(
        "primer,set_index,score\n"
        f"{twelve},0,1.0\n{fourteen},0,1.0\n{ten},1,0.5\n{twelve},1,0.5\n"
    )
    before = (results / "step4_improved_df.csv").read_text()

    proc = _run(
        [
            "improve-set",
            "--from-results",
            str(results),
            "--set",
            "1",
            "--genome",
            *workdir["genomes"],
            "--linear",
            "-o",
            "out",
        ],
        cwd=str(workdir["dir"]),
    )
    assert proc.returncode == 0, proc.stderr
    report = _report(workdir["dir"] / "out")

    assert sorted(report["primers"]) == sorted([ten, twelve])
    assert (results / "step4_improved_df.csv").read_text() == before
    assert not (workdir["dir"] / "out" / "step4_improved_df.csv").exists()


def test_neither_a_params_file_nor_a_genome_is_refused_by_name(tmp_path):
    """An explanation that names this command, not a traceback or a sibling's."""
    proc = _run(["improve-set", "--primers", "ATCGATCGATCG", "-o", "out"], cwd=str(tmp_path))

    assert proc.returncode != 0
    assert "improve-set needs either -j params.json or --genome" in proc.stderr
    assert "Traceback" not in proc.stderr
    assert not (tmp_path / "out" / "improvement_report.json").exists()


def test_the_command_is_registered_and_grouped():
    from neoswga.cli_unified import COMMAND_GROUPS, create_parser

    parser = create_parser()
    assert "improve-set" in parser._subparsers._group_actions[0].choices
    assert "improve-set" in {c for _, cmds in COMMAND_GROUPS for c in cmds}
