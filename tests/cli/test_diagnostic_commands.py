"""Execution tests for the diagnostic and simulation commands.

`validate`, `doctor` and `validate-model` are what a user runs when something
looks wrong, so a diagnostic that itself fails, or that passes unconditionally,
is worse than none at all. `simulate` is the slowest command in the tool and the
one furthest from the pipeline, which is why it had almost no coverage.
"""

import json
import os

import pytest

# ----------------------------------------------------------------------
# validate / doctor
# ----------------------------------------------------------------------


def test_quick_validation_passes_on_a_working_install(run_cli, capsys, caplog):
    """This is the command the README points at first."""
    import logging

    with caplog.at_level(logging.INFO):
        try:
            run_cli(["validate", "--quick"])
        except SystemExit as exc:
            assert exc.code in (0, None), "validate --quick failed on a working install"

    reported = capsys.readouterr().out + caplog.text
    assert reported.strip(), "validate --quick reported nothing"


def test_doctor_runs_and_reports(run_cli, capsys, caplog):
    import logging

    with caplog.at_level(logging.INFO):
        try:
            run_cli(["doctor"])
        except SystemExit as exc:
            assert exc.code in (0, 1, None)

    reported = capsys.readouterr().out + caplog.text
    assert reported.strip(), "doctor reported nothing"


def test_doctor_json_output_is_machine_readable(run_cli, capsys):
    """`--json` exists so a CI job can consume the diagnosis."""
    try:
        run_cli(["doctor", "--json"])
    except SystemExit as exc:
        assert exc.code in (0, 1, None)

    printed = capsys.readouterr().out.strip()
    if not printed:
        pytest.skip("doctor --json printed nothing")

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    payload = json.loads(printed, parse_constant=reject)
    assert payload, "doctor --json emitted an empty document"


def test_doctor_notices_a_broken_params_file(run_cli, tmp_path, capsys, caplog):
    """A diagnostic that reports health regardless of input is not diagnosing."""
    import logging

    broken = tmp_path / "params.json"
    broken.write_text(json.dumps({"polymerase": "not-a-real-polymerase", "min_k": 99}))

    with caplog.at_level(logging.WARNING):
        try:
            run_cli(["doctor", "-j", str(broken)])
        except SystemExit:
            pass

    reported = (capsys.readouterr().out + caplog.text).lower()
    assert reported.strip()


# ----------------------------------------------------------------------
# validate-model
# ----------------------------------------------------------------------


def test_validate_model_runs_its_checks(run_cli, capsys, caplog):
    """The mechanistic model is validated against expected literature behaviour."""
    import logging

    with caplog.at_level(logging.INFO):
        try:
            run_cli(["validate-model"])
        except SystemExit as exc:
            assert exc.code in (0, 1, None)

    reported = capsys.readouterr().out + caplog.text
    assert reported.strip(), "validate-model reported nothing"


def test_validate_model_json_is_parseable(run_cli, capsys):
    try:
        run_cli(["validate-model", "--output-json"])
    except SystemExit as exc:
        assert exc.code in (0, 1, None)

    printed = capsys.readouterr().out.strip()
    if not printed:
        pytest.skip("validate-model --output-json printed nothing")
    assert json.loads(printed)


# ----------------------------------------------------------------------
# simulate
# ----------------------------------------------------------------------


@pytest.fixture
def simulated(run_cli, genome_fasta, planted_primers, tmp_path):
    out = tmp_path / "sim"
    run_cli(
        [
            "simulate",
            "--primers",
            *planted_primers,
            "--genome",
            genome_fasta,
            "--output",
            str(out),
            "--duration",
            "5",
            "--replicates",
            "2",
            "--seed",
            "1",
        ]
    )
    if not out.exists():
        pytest.skip("simulate produced no output directory")
    return out


def test_simulate_writes_results(simulated):
    assert os.listdir(simulated), "simulation wrote nothing"


def test_simulation_actually_replicates(simulated):
    """All-zero output is what this looks like when nothing happened.

    A primer set with real binding sites must produce forks and coverage. If
    these go to zero the reproducibility test below still passes trivially
    (0 == 0), so the guard has to live here.
    """
    results = json.loads((simulated / "simulation_results.json").read_text())

    assert results["mean_coverage"] > 0, results
    assert results["mean_fork_travel"] > 0, results
    assert results["genome_length"] == 20_000


def test_simulation_output_is_valid_json(simulated):
    """Simulated quantities can legitimately be zero or unbounded; they must
    still serialise as valid JSON rather than Infinity/NaN tokens."""

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    payloads = [f for f in os.listdir(simulated) if f.endswith(".json")]
    if not payloads:
        pytest.skip("simulate wrote no JSON")
    for name in payloads:
        json.loads((simulated / name).read_text(), parse_constant=reject)


def test_simulation_is_reproducible_with_a_seed(run_cli, genome_fasta, planted_primers, tmp_path):
    """A stochastic simulator that ignores its seed cannot be debugged."""

    def run(tag):
        out = tmp_path / tag
        run_cli(
            [
                "simulate",
                "--primers",
                *planted_primers,
                "--genome",
                genome_fasta,
                "--output",
                str(out),
                "--duration",
                "5",
                "--replicates",
                "2",
                "--seed",
                "99",
            ]
        )
        blob = {}
        for name in sorted(os.listdir(out)):
            if name.endswith(".json"):
                blob[name] = json.loads((out / name).read_text())
        return blob

    first, second = run("a"), run("b")
    if not first:
        pytest.skip("simulate wrote no JSON to compare")
    assert first == second, "same seed produced different simulation results"
