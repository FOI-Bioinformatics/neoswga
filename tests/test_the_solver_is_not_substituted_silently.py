"""A MILP backend measured to kill the process is not substituted silently.

python-mip's default is CBC, and constructing a CBC model terminates the
interpreter on Python 3.13 -- SIGKILL, no traceback, nothing catchable, measured
on macOS arm64 with mip 2.0.0 and cbcbox 2.935. HiGHS solves the same models.

This used to prefer HiGHS and fall back to CBC with a warning. A warning that
precedes a SIGKILL is never read in context: the process dies with no traceback,
so the user sees a killed command and no reason to connect it to a log line.
Refusing costs one `pip install highsbox` and says so.

The fatality measurement is ONE platform, so an installation whose CBC works can
ask for it. That is an explicit request rather than a substitution nobody asked
for, which is the distinction Task 6 of the valid-design plan draws.

Nothing here constructs a model. A probe would have to, which is the operation
that kills the process, so no in-process test can tell a working CBC from a fatal
one -- and a test that tried would take the suite with it.
"""

import pytest

from neoswga.core import ilp_solver
from neoswga.core.exceptions import UnsupportedModelError

mip = pytest.importorskip("mip", reason="the MILP backend choice only matters with python-mip")


def test_highs_is_chosen_when_its_runtime_is_present(monkeypatch):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: True)
    assert ilp_solver.select_solver_name() == mip.HIGHS


def test_without_highs_it_refuses_rather_than_falling_back(monkeypatch):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.delenv(ilp_solver._ALLOW_CBC_ENV, raising=False)

    with pytest.raises(UnsupportedModelError) as caught:
        ilp_solver.select_solver_name()

    message = str(caught.value)
    assert "highsbox" in message, "the refusal must name what to install"
    assert ilp_solver._ALLOW_CBC_ENV in message, "and how to override it"
    assert "SIGKILL" in message, "and why it refuses"


def test_the_refusal_is_in_the_design_error_family(monkeypatch):
    """So the command boundary reports it and writes a failure record."""
    from neoswga.core.exceptions import DesignError

    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.delenv(ilp_solver._ALLOW_CBC_ENV, raising=False)

    with pytest.raises(DesignError):
        ilp_solver.select_solver_name()


@pytest.mark.parametrize("value", ["1", "true", "TRUE", "yes"])
def test_an_explicit_request_for_cbc_is_honoured(monkeypatch, value):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.setenv(ilp_solver._ALLOW_CBC_ENV, value)

    assert ilp_solver.select_solver_name() == mip.CBC


@pytest.mark.parametrize("value", ["", "0", "no", "maybe"])
def test_anything_other_than_an_opt_in_still_refuses(monkeypatch, value):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.setenv(ilp_solver._ALLOW_CBC_ENV, value)

    with pytest.raises(UnsupportedModelError):
        ilp_solver.select_solver_name()


def test_the_opt_in_warns_that_cbc_has_been_measured_to_die(monkeypatch, caplog):
    monkeypatch.setattr(ilp_solver, "highs_is_available", lambda: False)
    monkeypatch.setenv(ilp_solver._ALLOW_CBC_ENV, "1")

    with caplog.at_level("WARNING"):
        ilp_solver.select_solver_name()

    assert "SIGKILL" in caplog.text
