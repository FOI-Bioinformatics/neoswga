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


def test_cbc_is_no_longer_the_silent_fallback(monkeypatch):
    """Changed deliberately on 2026-09-27; this test asserted the old policy.

    CBC was substituted with a warning. A warning that precedes a SIGKILL is
    never read in context: the process dies with no traceback, so the user sees
    a killed command and no reason to connect it to a log line. Refusing costs
    one `pip install highsbox` and says so.

    The behaviour it used to pin now lives behind an explicit request; see
    `tests/test_the_solver_is_not_substituted_silently.py`.
    """
    from neoswga.core.exceptions import UnsupportedModelError

    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.delenv(ilp_solver._ALLOW_CBC_ENV, raising=False)
    with pytest.raises(UnsupportedModelError):
        ilp_solver.select_solver_name()


def test_an_explicit_request_for_cbc_still_says_why_it_is_risky(monkeypatch, caplog):
    """The fatality measurement is one platform, so the opt-in exists -- and it
    still states what was measured, because asking for CBC is not the same as
    knowing it is safe here."""
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.setenv(ilp_solver._ALLOW_CBC_ENV, "1")
    with caplog.at_level(logging.WARNING):
        assert ilp_solver.select_solver_name() == mip.CBC
    assert "SIGKILL" in caplog.text, "the warning must say what it is warning about"


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
