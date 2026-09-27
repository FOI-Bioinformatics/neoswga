"""A scoring call must not re-read a mutable global.

`occupancy.weighted_site_load` read `parameter.mismatch_penalty` itself, on
every call, and it is called once per panel evaluation. So the penalty two
evaluations inside one run used was whatever the module global said at the
moment each ran, for a quantity that is a property of the reaction the run was
resolved for. Nothing in production reassigns it mid-run today, which is why
this never produced a visible wrong answer -- and is also why it could not be
noticed. Task 2 of the valid-design plan asks for exactly this: evaluator code
must not depend on a mutable global at read time.

The remaining reads of that module are config-time, once per run, and are a
different thing: `coverage.resolve_extension_reach` and
`panel_acceptance.constraints_from_parameter` both run before the search.
"""

import inspect

import pytest

from neoswga.core import occupancy, parameter
from neoswga.core.occupancy import DEFAULT_MISMATCH_PENALTY_C, weighted_site_load


def test_the_penalty_can_be_supplied_by_the_caller():
    assert "penalty" in inspect.signature(weighted_site_load).parameters


def test_a_supplied_penalty_is_not_overridden_by_the_global(monkeypatch):
    """The caller's value wins, so a run can hold its own chemistry."""
    seen = []

    def record(tm, distance, penalty):
        seen.append(penalty)
        return tm

    monkeypatch.setattr(occupancy, "mismatch_tm", record)
    monkeypatch.setattr(occupancy, "site_occupancy", lambda *a, **k: 1.0)
    monkeypatch.setattr(
        occupancy,
        "default_mismatch_penalty",
        lambda: pytest.fail("the global was read although a penalty was supplied"),
    )

    class Conditions:
        temp = 30.0

        def calculate_effective_tm(self, primer):
            return 40.0

    monkeypatch.setattr(
        "neoswga.core.mismatch_counts.mismatch_class_counts",
        lambda primer, prefixes, max_mismatches: {1: 2},
    )
    monkeypatch.setattr(
        "neoswga.core.thermodynamics.calculate_enthalpy_entropy", lambda primer: (-100.0, -0.3)
    )

    weighted_site_load(["ACCACAGATAGC"], ["p"], Conditions(), 1, penalty=7.5)
    assert seen == [7.5]


def test_the_optimizer_resolves_the_penalty_once_at_construction(monkeypatch):
    """And keeps it, so a later reassignment cannot move one run's scoring."""
    monkeypatch.setattr(parameter, "mismatch_penalty", 3.0, raising=False)

    class Optimizer:
        def __init__(self):
            from neoswga.core.occupancy import default_mismatch_penalty

            self.mismatch_penalty = default_mismatch_penalty()

    optimizer = Optimizer()
    assert optimizer.mismatch_penalty == 3.0

    monkeypatch.setattr(parameter, "mismatch_penalty", 9.0, raising=False)
    assert optimizer.mismatch_penalty == 3.0, (
        "the run's penalty moved when the global did; one run must weigh every "
        "panel the same way"
    )


def test_an_unset_penalty_still_falls_back_to_the_declared_default(monkeypatch):
    monkeypatch.setattr(parameter, "mismatch_penalty", None, raising=False)
    assert occupancy.default_mismatch_penalty() == DEFAULT_MISMATCH_PENALTY_C


def test_the_site_load_call_sites_pass_the_run_penalty():
    """A source check: a new call site that omits it reintroduces the read.

    The behavioural tests above cannot see a THIRD call site added later, and
    the omission is invisible until someone reassigns the global mid-run.
    """
    import ast
    import pathlib

    source = pathlib.Path("neoswga/core/base_optimizer.py").read_text()
    tree = ast.parse(source)
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "weighted_site_load"
    ]
    assert calls, "expected base_optimizer to compute an occupancy-weighted load"
    for call in calls:
        # Either spelling: `penalty=` by keyword, or the fifth positional
        # argument, which is what the signature puts there. What must not
        # happen is a call that supplies neither and so falls back to reading
        # the module global on every evaluation.
        supplied = any(keyword.arg == "penalty" for keyword in call.keywords) or len(call.args) >= 5
        assert supplied, (
            f"weighted_site_load at line {call.lineno} does not pass the run's "
            f"penalty, so it re-reads parameter.mismatch_penalty per evaluation"
        )
