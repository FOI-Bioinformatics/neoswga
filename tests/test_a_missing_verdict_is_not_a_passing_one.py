"""An unreadable verdict must not read as a clean one.

`design_result.blocking_validator_findings` states the rule in its own
docstring -- absence of evidence is not evidence, and "A corrupt *failure*
record is the opposite case and blocks, because there the evidence exists and
cannot be read" -- and then applies it to the failure record only. Its own
`except (OSError, ValueError): return []` treats a corrupt validation record
exactly like an absent one. `results_interpreter._blocks_synthesis` does the
same.

So a `step4_improved_df_validation.json` truncated by a full disk, a killed
process or a half-finished copy made `neoswga export` print "Primers ready for
ordering!" for a pool nobody could vouch for.

Missing stays permissive, and that is deliberate rather than an oversight:
directories written before the validator existed are explicitly supported, and
refusing them would be refusing on the absence of evidence. Corrupt is the
other case. The evidence exists, it cannot be read, and the reason a reader
cannot tell whether it named a blocking finding is exactly why it must not be
assumed it did not.

The write is the third leg. `_write_validation_report` swallowed `OSError` and
continued, so a run that produced a panel and could not record its verdict left
a directory that every consumer then read as clean -- the same end state by a
different route, and the one `unified_optimizer` already converts a validator
EXCEPTION into `ModelEvaluationError` to prevent.

The unset-`data_dir` skip is left alone. That function's docstring argues it and
the argument holds: an unset `data_dir` is not a location, the previous
`os.getcwd()` fallback wrote the file into the repository root, and there is no
directory for `export` to misread afterwards.
"""

import json

import pytest

from neoswga.core.design_result import VALIDATION_FILENAME

CLEAN = {
    "optimizer": "test",
    "num_primers": 2,
    "ok": True,
    "issues": [],
}

BLOCKING = {
    "optimizer": "test",
    "num_primers": 2,
    "ok": False,
    "issues": [
        {
            "level": "error",
            "code": "delivered_pool_exceeds_max_dimer_bp",
            "detail": "worst heterodimer is 9 bp against a configured max_dimer_bp of 3",
        }
    ],
}


def _write(directory, payload):
    path = directory / VALIDATION_FILENAME
    if isinstance(payload, str):
        path.write_text(payload)
    else:
        path.write_text(json.dumps(payload))
    return path


# ---------------------------------------------------------------------------
# The gate: design_result.blocking_validator_findings
# ---------------------------------------------------------------------------


def test_an_unreadable_verdict_is_a_blocking_finding(tmp_path):
    """The defect. Truncated JSON read as a clean pool."""
    from neoswga.core.design_result import blocking_validator_findings

    _write(tmp_path, '{"optimizer": "test", "issues": [{"code":')

    findings = blocking_validator_findings(tmp_path)

    assert findings, "a validation record that cannot be read was treated as clean"
    assert any("read" in f.lower() or "unreadable" in f.lower() for f in findings), findings


def test_a_missing_verdict_is_still_not_a_finding(tmp_path):
    """Deliberate, and the reason is in this module's docstring. Directories
    written before the validator existed are supported."""
    from neoswga.core.design_result import blocking_validator_findings

    assert blocking_validator_findings(tmp_path) == []


def test_a_clean_verdict_is_still_not_a_finding(tmp_path):
    """Guard the guard: a gate that blocks everything is no better."""
    from neoswga.core.design_result import blocking_validator_findings

    _write(tmp_path, CLEAN)

    assert blocking_validator_findings(tmp_path) == []


def test_a_recorded_blocking_finding_still_blocks(tmp_path):
    from neoswga.core.design_result import blocking_validator_findings

    _write(tmp_path, BLOCKING)

    findings = blocking_validator_findings(tmp_path)
    assert len(findings) == 1
    assert "heterodimer" in findings[0]


def test_the_export_refuses_an_unreadable_verdict(tmp_path):
    """Through the command's own gate, not only the helper."""
    from neoswga.core.export import export_is_blocked

    _write(tmp_path, "not json at all")

    assert export_is_blocked(str(tmp_path))


def test_the_export_allows_a_missing_verdict(tmp_path):
    from neoswga.core.export import export_is_blocked

    assert export_is_blocked(str(tmp_path)) is None


# ---------------------------------------------------------------------------
# The second gate: results_interpreter
# ---------------------------------------------------------------------------


def test_interpret_treats_an_unreadable_verdict_as_blocking(tmp_path):
    """Two commands tell a user a pool is ready. They must not disagree."""
    from neoswga.core.results_interpreter import ResultsInterpreter

    _write(tmp_path, '{"issues": [')

    interpreter = ResultsInterpreter.__new__(ResultsInterpreter)
    interpreter.validation_file = tmp_path / VALIDATION_FILENAME

    assert interpreter._blocks_synthesis() is True


def test_interpret_treats_a_missing_verdict_as_not_blocking(tmp_path):
    from neoswga.core.results_interpreter import ResultsInterpreter

    interpreter = ResultsInterpreter.__new__(ResultsInterpreter)
    interpreter.validation_file = tmp_path / VALIDATION_FILENAME

    assert interpreter._blocks_synthesis() is False


# ---------------------------------------------------------------------------
# The write
# ---------------------------------------------------------------------------


def test_a_write_that_fails_is_not_swallowed(tmp_path, monkeypatch):
    """A run that produced a panel and could not record its verdict leaves a
    directory every consumer reads as clean. That is the same end state the
    caller already raises `ModelEvaluationError` to prevent when the validator
    itself throws."""
    from neoswga.core import parameter, step4_output

    monkeypatch.setattr(parameter, "data_dir", str(tmp_path / "no" / "such" / "dir"), raising=False)

    with pytest.raises(OSError):
        step4_output._write_validation_report(CLEAN)


def test_no_data_dir_still_skips_quietly(tmp_path, monkeypatch):
    """Left alone on purpose; the function's own docstring argues it. An unset
    `data_dir` is not a location, and there is no directory for `export` to
    misread afterwards."""
    from neoswga.core import parameter, step4_output

    monkeypatch.setattr(parameter, "data_dir", None, raising=False)

    step4_output._write_validation_report(CLEAN)  # must not raise
