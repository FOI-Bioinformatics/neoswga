"""The MILP backend is chosen, not defaulted.

python-mip's default is CBC, and constructing a CBC model terminates the
interpreter with SIGKILL on Python 3.13 (measured on macOS arm64, mip 2.0.0,
cbcbox 2.935): no exception, no traceback, no stderr. HiGHS solves the same
models on 3.13 and on 3.11, so it is preferred wherever its runtime is
importable.

The fallback is deliberately unprobed. A probe would construct a CBC model,
which is the operation that kills the process, so nothing in process can tell
a working CBC from a fatal one.
"""

import logging

import pytest

pytest.importorskip("mip", reason="python-mip ships in the 'improved' extra")

import mip

from neoswga.core import ilp_solver


def test_highs_is_preferred_when_its_runtime_is_present(monkeypatch):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: True)
    assert ilp_solver.select_solver_name() == mip.HIGHS


def test_cbc_is_the_fallback(monkeypatch):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    assert ilp_solver.select_solver_name() == mip.CBC


def test_falling_back_says_why_and_how_to_fix_it(monkeypatch, caplog):
    """Unconditional, because every supported interpreter is one where CBC has
    been measured to die. A 3.13 version guard here was dead code once
    `requires-python` became ">=3.13", and ruff said so."""
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    with caplog.at_level(logging.WARNING):
        ilp_solver.select_solver_name()
    assert "highsbox" in caplog.text, "the warning must name the package that fixes it"
    assert "SIGKILL" in caplog.text, "and what it is warning about"


def test_the_preferred_path_is_quiet(monkeypatch, caplog):
    """No warning when HiGHS is used, or every ILP run would emit one."""
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: True)
    with caplog.at_level(logging.WARNING):
        ilp_solver.select_solver_name()
    assert caplog.text == ""


def test_availability_is_decided_by_import_not_by_construction(monkeypatch):
    """Constructing a model is what dies, so availability cannot be tested
    that way. mip raises a plain ImportError naming its runtime when HiGHS is
    unusable, which makes a spec lookup sufficient and safe."""
    monkeypatch.setattr(ilp_solver, "_HIGHS_RUNTIME", "a_module_that_does_not_exist")
    assert ilp_solver.highs_is_available() is False


def test_the_optimizer_asks_for_a_solver_rather_than_defaulting():
    """Pins the PATH. A selector nothing calls is the defect this repository
    names Known Issue 8, and the benchmark script needed wiring too: fixing the
    library alone left 16 of 17 tests passing and the 17th killing the run."""
    import inspect

    from neoswga.core import dominating_set_optimizer

    source = inspect.getsource(dominating_set_optimizer)
    assert "select_solver_name()" in source
    assert "Model(sense=MAXIMIZE)" not in source, "a bare Model() takes CBC by default"
