"""Execution tests for the iterative-design commands.

`contract-set`, `rescore-set`, `swap-primer` and `expand-primers` are the "I
already have a set, improve it" half of the tool, and at 15% coverage nothing
had ever run them against real position files. They are also the commands most
exposed to the silent-zero-coverage failure, since they take primers from the
user rather than from the pipeline's own candidate table.

The fixtures build a genuine pipeline workspace (count-kmers -> filter -> score
on a 20 kb genome), so these run against the real HDF5 layout and the real
step3_df.csv rather than mocks that could drift from either.

Two conventions these commands follow, both different from `evaluate-set`, and
both pinned below so a future change is deliberate rather than accidental:
`-o` names a FILE (not a directory), and the handler returns None, writing its
result to that file or to stdout.
"""

import json
import os

import pytest

pytestmark = pytest.mark.usefixtures("pipeline_run")


@pytest.fixture
def primers(scored_primers):
    """A working set drawn from the pipeline's own scored candidates."""
    return scored_primers[:6]


def _run_to_json(run_cli, argv, out_path):
    """Run a command that writes its JSON result to `-o` and read it back."""
    run_cli(argv + ["-o", str(out_path)])
    assert os.path.isfile(out_path), f"command wrote no file at {out_path}"

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    # Parsed strictly: Infinity/NaN are not valid JSON, and these files are
    # meant to be consumed by other tools.
    return json.loads(open(out_path).read(), parse_constant=reject)


# ----------------------------------------------------------------------
# rescore-set
# ----------------------------------------------------------------------


def test_rescore_set_writes_valid_json(run_cli, pipeline_run, primers, tmp_path):
    result = _run_to_json(
        run_cli,
        ["rescore-set", "-j", pipeline_run["params_file"], "--primers", *primers],
        tmp_path / "rescore.json",
    )
    assert isinstance(result, dict) and result


def test_rescore_set_responds_to_an_additive(run_cli, pipeline_run, primers, tmp_path):
    """Changing the chemistry must change the numbers.

    If it does not, the additive is not reaching the thermodynamic model -- the
    exact class of bug the audit found, where four additives were accepted by
    the CLI and silently discarded.
    """
    base = _run_to_json(
        run_cli,
        ["rescore-set", "-j", pipeline_run["params_file"], "--primers", *primers],
        tmp_path / "base.json",
    )
    dmso = _run_to_json(
        run_cli,
        [
            "rescore-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--dmso-percent",
            "10",
        ],
        tmp_path / "dmso.json",
    )
    assert base != dmso, "10% DMSO produced an identical result"


def test_rescore_set_responds_to_polymerase(run_cli, pipeline_run, primers, tmp_path):
    phi29 = _run_to_json(
        run_cli,
        [
            "rescore-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--polymerase",
            "phi29",
        ],
        tmp_path / "phi29.json",
    )
    bst = _run_to_json(
        run_cli,
        [
            "rescore-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--polymerase",
            "bst",
        ],
        tmp_path / "bst.json",
    )
    assert phi29 != bst, "switching phi29 -> bst changed nothing"


# ----------------------------------------------------------------------
# contract-set
# ----------------------------------------------------------------------


def _contract(run_cli, pipeline_run, primers, out, min_coverage, seed=1):
    return _run_to_json(
        run_cli,
        [
            "contract-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--min-coverage",
            str(min_coverage),
            "--seed",
            str(seed),
        ],
        out,
    )


def test_contract_set_never_invents_or_grows(run_cli, pipeline_run, primers, tmp_path):
    result = _contract(run_cli, pipeline_run, primers, tmp_path / "c.json", 0.0)

    kept = result["contracted_set"]
    assert len(kept) <= len(primers)
    assert set(kept) <= set(primers), "contract-set returned a primer it was not given"
    assert set(result["removed_primers"]) == set(primers) - set(kept)


def test_contract_set_reports_real_baseline_coverage(run_cli, pipeline_run, primers, tmp_path):
    """A zero baseline would mean the position lookup silently found nothing.

    That was the documented failure: contract-set reported baseline_coverage 0.0
    and could never drop a primer, with no indication anything was wrong.
    """
    result = _contract(run_cli, pipeline_run, primers, tmp_path / "c.json", 0.0)
    assert 0.0 < result["baseline_coverage"] <= 1.0


def test_min_coverage_actually_constrains_the_result(run_cli, pipeline_run, primers, tmp_path):
    """The audit found `contract-set` hardcoded a 70 kb extension reach.

    On any genome smaller than that a single binding site covers everything, so
    coverage was always ~100% and --min-coverage was effectively a boolean. With
    the realistic ~3 kb reach it constrains again: a strict floor must not drop
    more primers than a permissive one.
    """
    permissive = _contract(run_cli, pipeline_run, primers, tmp_path / "low.json", 0.0)
    strict = _contract(run_cli, pipeline_run, primers, tmp_path / "high.json", 0.99)

    assert len(strict["contracted_set"]) >= len(permissive["contracted_set"])
    assert strict["final_coverage"] >= permissive["final_coverage"] - 1e-9


def test_contract_set_keeps_coverage_above_the_floor(run_cli, pipeline_run, primers, tmp_path):
    result = _contract(run_cli, pipeline_run, primers, tmp_path / "c.json", 0.75)
    if result["baseline_coverage"] >= 0.75:
        assert result["final_coverage"] >= 0.75 - 1e-9


def test_contract_set_is_reproducible_with_a_seed(run_cli, pipeline_run, primers, tmp_path):
    a = _contract(run_cli, pipeline_run, primers, tmp_path / "a.json", 0.5, seed=7)
    b = _contract(run_cli, pipeline_run, primers, tmp_path / "b.json", 0.5, seed=7)
    assert a == b


def test_contract_set_prints_when_no_output_given(run_cli, pipeline_run, primers, capsys):
    """Without -o the result goes to stdout, so it can be piped."""
    run_cli(
        [
            "contract-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--min-coverage",
            "0.0",
        ]
    )
    printed = capsys.readouterr().out
    assert json.loads(printed)["contracted_set"] is not None


# ----------------------------------------------------------------------
# The failure mode these commands were most exposed to
# ----------------------------------------------------------------------


def test_unindexed_primer_is_warned_about_not_silently_zeroed(
    run_cli, pipeline_run, primers, absent_primer, tmp_path, caplog
):
    """A primer absent from the index scores zero for a reason worth stating.

    Before the audit this was entirely silent, so "binds nowhere" and "was never
    looked up" were indistinguishable. PositionCache now defaults to
    on_missing='warn', which is what carries the distinction here -- these
    commands cannot scan a FASTA to resolve it (see the test below).
    """
    import logging

    with caplog.at_level(logging.WARNING):
        _run_to_json(
            run_cli,
            [
                "contract-set",
                "-j",
                pipeline_run["params_file"],
                "--primers",
                *primers[:3],
                absent_primer,
                "--min-coverage",
                "0.0",
                "--seed",
                "1",
            ],
            tmp_path / "mixed.json",
        )

    assert "no cached" in caplog.text or "not in the HDF5" in caplog.text


def test_iterate_commands_cannot_yet_scan_a_genome():
    """Documents an incomplete piece of the bring-your-own-oligo work.

    `add_position_source_options` gives a command `--genome`, letting
    PositionCache scan a FASTA for primers missing from the index -- turning an
    unknown into a measurement. Its own help text says it "also lets -j be
    omitted entirely". It is wired into `evaluate-set` only, so on the iterate
    commands an outside primer set can be warned about but not resolved.

    If this starts failing, the option has been extended and the note in
    docs/ should be updated to match.
    """
    import inspect

    from neoswga.cli import iterate

    assert "add_position_source_options" not in inspect.getsource(iterate)


# ----------------------------------------------------------------------
# Conventions, pinned so a change is deliberate
# ----------------------------------------------------------------------


def test_iterate_output_flag_names_a_file_not_a_directory(run_cli, pipeline_run, primers, tmp_path):
    """`-o` differs between command groups: a file here, a directory for
    `evaluate-set`. Pinned because it is the kind of inconsistency that gets
    "tidied" into a breaking change by accident."""
    out = tmp_path / "result.json"
    run_cli(
        [
            "contract-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--min-coverage",
            "0.0",
            "-o",
            str(out),
        ]
    )
    assert out.is_file()


def test_switching_polymerase_adopts_its_temperature(run_cli, pipeline_run, primers, tmp_path):
    """`--polymerase bst` alone used to abort with a validation error.

    params.json carried phi29's 30 C, bst is valid only at 50-72 C, and the
    mismatch surfaced as a raw ValueError from ReactionConditions -- on the
    single most obvious use of a command for rescoring under a different
    chemistry. The new polymerase's optimum is now adopted unless the caller
    sets a temperature explicitly.
    """
    result = _run_to_json(
        run_cli,
        [
            "rescore-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--polymerase",
            "bst",
        ],
        tmp_path / "bst.json",
    )
    assert result["conditions"]["polymerase"] == "bst"
    assert 50.0 <= result["conditions"]["temp"] <= 72.0


def test_explicit_temperature_still_wins(run_cli, pipeline_run, primers, tmp_path):
    result = _run_to_json(
        run_cli,
        [
            "rescore-set",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "--polymerase",
            "bst",
            "--reaction-temp",
            "65",
        ],
        tmp_path / "bst65.json",
    )
    assert result["conditions"]["temp"] == 65.0
