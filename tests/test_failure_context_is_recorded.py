"""The next unexplained suite failure should carry its own evidence.

One roadmap finding has no cause. On 2026-09-10 a full run gave 6 failures and
4124 passes; two immediate re-runs on a byte-identical tree gave none. The only
unusual condition was 25 concurrent agent processes. The entry says a single
test name was not enough to diagnose it, which is exactly why the entry cannot
say more.

Three known sources of load-dependent failure have since been removed: the
wall-clock assertions, an over-broad caplog assertion, and subprocess tests that
skipped when a step failed. What remains, if anything does, is this one.

A plan cannot fix what it cannot reproduce, so this records what a diagnosis
would need instead. On a green run it writes nothing.
"""

import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / ".pytest_failure_context.json"


def test_the_hook_is_installed():
    """The recording happens in a pytest hook, so nothing exercises it on a
    green run. This at least pins that the hook exists and is wired."""
    import tests.conftest as conftest

    assert hasattr(conftest, "pytest_runtest_makereport")
    assert hasattr(conftest, "pytest_sessionfinish")


def test_the_report_is_gitignored():
    """It is a diagnostic, not an artifact of the run, and the suite must still
    leave the tree unchanged -- see
    tests/test_the_suite_leaves_no_files_behind.py."""
    ignored = (ROOT / ".gitignore").read_text()
    assert ".pytest_failure_context.json" in ignored


def test_a_recorded_failure_carries_what_a_diagnosis_needs(tmp_path, monkeypatch):
    """Drive the writer directly rather than wait for a flake.

    The roadmap entry says a test name alone was not enough. What was missing:
    which worker ran it, how loaded the machine was, and how many other tests
    were running beside it.
    """
    import tests.conftest as conftest

    monkeypatch.setattr(conftest, "_FAILURE_REPORT", tmp_path / "ctx.json")
    monkeypatch.setattr(
        conftest, "_FAILURES", [{"test": "tests/test_x.py::test_y", "worker": "gw3"}]
    )

    conftest._write_failure_context()

    recorded = json.loads((tmp_path / "ctx.json").read_text())
    assert recorded["failures"][0]["test"] == "tests/test_x.py::test_y"
    assert recorded["failures"][0]["worker"] == "gw3"
    assert "load_average" in recorded
    assert "workers" in recorded


def test_a_green_run_writes_nothing(tmp_path, monkeypatch):
    import tests.conftest as conftest

    monkeypatch.setattr(conftest, "_FAILURE_REPORT", tmp_path / "ctx.json")
    monkeypatch.setattr(conftest, "_FAILURES", [])

    conftest._write_failure_context()

    assert not (tmp_path / "ctx.json").exists()
