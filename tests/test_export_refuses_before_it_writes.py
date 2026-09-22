"""An unqualified pool must not reach disk as an oligo order.

`neoswga export` had the right check in the wrong place. `run_export` wrote
every order file first and consulted `_blocking_validator_findings` afterwards,
so a pool breaking the user's configured `max_dimer_bp` produced a complete
FASTA and then printed "NOT ready for ordering" beside it. The exit code was
zero either way. A person who read the message still had the file, and any
script that checked the status code had no idea.

Found 22 September 2026 while acting on an external comparative review. The
review did not contain it; it was found by reading `run_export` top to bottom
after the review pointed at the export path for a different reason.

Two other things made it worse than a misplaced call:

- `export_is_blocked`, the gate that IS correctly positioned before the writes,
  read only `design_failure.json`. That file is written only on FAILURE, always
  with `qualified: False`, and is deleted when step 4 succeeds. So the formal
  gate `design_result.recommendation_allowed(state, qualified)` never once saw a
  finished run with a real qualification value. It could not block anything.
- The set of blocking codes was a bare literal in two places, `cli/report.py`
  and `results_interpreter.py`, so a new blocking code had to be added twice or
  the two commands would disagree about the same pool.

The fix is not a new mechanism. The validator findings move INTO
`export_is_blocked`, which already runs before anything is written, and the
code set moves to one place both readers import.
"""

import json
import subprocess
import sys

import pytest

RESULTS = "step4_improved_df.csv"
VALIDATION = "step4_improved_df_validation.json"

CLEAN_PANEL = ["ACCACAGATAGC", "GTTGTAGATGGA", "ATCAGCAGACCA"]


def write_results(directory, primers=CLEAN_PANEL):
    rows = ["primer,set_index,score,coverage"]
    rows += [f"{p},0,1.0,0.5" for p in primers]
    (directory / RESULTS).write_text("\n".join(rows) + "\n")


def write_validation(directory, codes, ok=True):
    payload = {
        "ok": ok,
        "issues": [
            {"code": code, "level": "warning", "detail": f"detail for {code}"} for code in codes
        ],
    }
    (directory / VALIDATION).write_text(json.dumps(payload))


def run_export(directory, *extra):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "neoswga.cli_unified",
            "export",
            "-d",
            str(directory),
            "--format",
            "fasta",
            "-o",
            str(directory / "out"),
            *extra,
        ],
        capture_output=True,
        text=True,
        timeout=300,
    )


def exported_files(directory):
    out = directory / "out"
    return sorted(p.name for p in out.iterdir()) if out.exists() else []


# ---------------------------------------------------------------------------
# The defect
# ---------------------------------------------------------------------------


def test_a_blocking_finding_writes_no_file(tmp_path):
    """The whole point. Nothing may reach disk before the pool has passed."""
    write_results(tmp_path)
    write_validation(tmp_path, ["delivered_pool_exceeds_max_dimer_bp"])

    result = run_export(tmp_path)

    assert exported_files(tmp_path) == [], (
        "an order file was written for a pool the optimizer recorded as "
        "breaking its configured dimer limit"
    )


def test_a_blocking_finding_exits_nonzero(tmp_path):
    """A script that checks the status code must be able to tell."""
    write_results(tmp_path)
    write_validation(tmp_path, ["delivered_pool_exceeds_max_dimer_bp"])

    assert run_export(tmp_path).returncode != 0


def test_the_refusal_names_the_finding_and_the_remedy(tmp_path):
    write_results(tmp_path)
    write_validation(tmp_path, ["delivered_pool_exceeds_max_dimer_bp"])

    result = run_export(tmp_path)
    output = result.stdout + result.stderr

    assert "delivered_pool_exceeds_max_dimer_bp" in output or "dimer" in output.lower()
    assert "interpret" in output, "the user needs to be told where to see the full assessment"


# ---------------------------------------------------------------------------
# What must still work
# ---------------------------------------------------------------------------


def test_a_clean_pool_still_exports(tmp_path):
    """Guard the guard: a gate that refuses everything also passes the tests above."""
    write_results(tmp_path)
    write_validation(tmp_path, [])

    result = run_export(tmp_path)

    assert result.returncode == 0, result.stdout + result.stderr
    assert exported_files(tmp_path), "a qualifying pool must still produce an order file"


def test_a_missing_validation_file_does_not_block(tmp_path):
    """An export must not fail closed because a validator never ran.

    Directories written before the validator existed are explicitly supported,
    and refusing them would be the opposite failure: a refusal on absence of
    evidence rather than on evidence.
    """
    write_results(tmp_path)

    result = run_export(tmp_path)

    assert result.returncode == 0, result.stdout + result.stderr
    assert exported_files(tmp_path)


def test_a_saturation_warning_does_not_block(tmp_path):
    """It says a metric is untrustworthy on a small target, which is inherent
    to plasmid-scale design rather than a defect in the pool."""
    write_results(tmp_path)
    write_validation(tmp_path, ["coverage_saturated_on_small_genome"])

    result = run_export(tmp_path)

    assert result.returncode == 0, result.stdout + result.stderr
    assert exported_files(tmp_path)


# ---------------------------------------------------------------------------
# The escape hatch
# ---------------------------------------------------------------------------


def test_an_explicit_override_exports_and_warns(tmp_path):
    """A refusal with no way forward would send people to the CSV by hand.

    `--allow-unqualified` is the same shape as `--allow-dimer-relaxation`: the
    constraint is traded deliberately, by name, and the run says so loudly.
    """
    write_results(tmp_path)
    write_validation(tmp_path, ["delivered_pool_exceeds_max_dimer_bp"])

    result = run_export(tmp_path, "--allow-unqualified")
    output = result.stdout + result.stderr

    assert result.returncode == 0, output
    assert exported_files(tmp_path), "the override must actually export"
    assert "delivered_pool_exceeds_max_dimer_bp" in output or "dimer" in output.lower()


# ---------------------------------------------------------------------------
# One rule, one place
# ---------------------------------------------------------------------------


def test_the_blocking_codes_are_defined_once():
    """They were a bare literal in two modules, so a new code had to be added
    twice or `export` and `interpret` would disagree about the same pool."""
    from neoswga.core.design_result import BLOCKING_VALIDATOR_CODES

    assert "delivered_pool_exceeds_max_dimer_bp" in BLOCKING_VALIDATOR_CODES

    import inspect

    from neoswga.cli import report
    from neoswga.core import results_interpreter

    for module in (report, results_interpreter):
        source = inspect.getsource(module)
        assert "delivered_pool_exceeds_max_dimer_bp" not in source, (
            f"{module.__name__} still hardcodes a blocking code; import "
            "BLOCKING_VALIDATOR_CODES instead"
        )
