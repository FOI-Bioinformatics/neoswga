"""Nineteen test files depend on artifacts a fixture may not have built.

The roadmap recorded this as "an upper bound of 82 tests vanish silently
without jellyfish". Both halves were wrong, and the correction matters because
it changes what the fix has to do.

`examples/plasmid_example/` is COMMITTED: git tracks its README, both FASTAs,
params.json and step2_df.csv.original. Everything else there -- the k-mer
tables, the HDF5 position files, step2_df.csv, step3_df.csv -- is a build
artifact produced by `tests/conftest.py::_prime_plasmid_example`, which runs
only when jellyfish is on PATH.

So the nineteen `if not EXAMPLE_DIR.is_dir(): pytest.skip(...)` guards never
fire: the directory is always there. Without jellyfish those tests do not
vanish, they FAIL, with an error about a missing k-mer file rather than about a
missing tool. Only three files guard on `check_jellyfish_available()`.

The real figure is 19 files and 44 test functions, not 23 and 82.

This file pins the two things that make the dependency legible: a priming
failure is reported rather than swallowed, and a guard checks for what the
tests actually need.
"""

import os
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "examples" / "plasmid_example"


def test_the_example_directory_is_committed_so_its_presence_proves_nothing():
    """The premise of the correction, asserted rather than assumed.

    If this ever stops being true, `is_dir()` becomes a meaningful guard again
    and the reasoning in this file needs revisiting.
    """
    import subprocess

    tracked = subprocess.run(
        ["git", "ls-files", "examples/plasmid_example"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=60,
    ).stdout.split()

    assert tracked, "examples/plasmid_example is not tracked at all"
    assert any(name.endswith("params.json") for name in tracked)
    assert not any(
        name.endswith(("step3_df.csv", "_positions.h5")) for name in tracked
    ), "a generated artifact is committed; the guard question changes"


def test_a_priming_failure_is_reported(monkeypatch, caplog):
    """`except Exception: pass` hid a priming failure WITH jellyfish present.

    A swallowed failure sends nineteen files into a confusing error about a
    missing k-mer file, pointing at the wrong thing. Driven here rather than
    read out of the source: the priming body is a function, so a step can be
    made to raise and the handler's behaviour observed.
    """
    import logging

    import tests.conftest as conftest

    def _boom():
        raise RuntimeError("jellyfish produced no output")

    monkeypatch.setattr(conftest, "_run_priming", _boom)
    monkeypatch.setattr(conftest, "plasmid_example_ready", lambda: False)
    monkeypatch.setattr("neoswga.core.kmer_counter.check_jellyfish_available", lambda: True)

    with caplog.at_level(logging.WARNING):
        conftest._prime_plasmid_example.__wrapped__()

    assert "Priming examples/plasmid_example failed" in caplog.text
    assert (
        "jellyfish produced no output" in caplog.text
    ), "the cause is not in the report, so a reader still has to guess"


def test_the_shared_guard_checks_for_the_artifacts_not_the_directory():
    """The guard has to test what the dependent tests actually need."""
    from tests.conftest import plasmid_example_ready

    assert callable(plasmid_example_ready)


@pytest.mark.skipif(
    not (EXAMPLE / "step3_df.csv").exists(),
    reason="plasmid example not primed (jellyfish absent?)",
)
def test_when_primed_the_guard_agrees_the_example_is_usable():
    from tests.conftest import plasmid_example_ready

    assert plasmid_example_ready() is True
