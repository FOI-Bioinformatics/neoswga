"""Running out of allowance mid-stage must not discard what the stage achieved.

`reduce_result` removes one primer per iteration, and each iteration's smaller
panel is a valid incumbent: it met the coverage target and violated no
constraint, which is why it was accepted. The SHARED allowance can run out
inside that loop, because `objective.coverage` is the budgeted call, and the
exception used to propagate out of the function. `run_panel_search`'s `execute`
then returned the PRE-STAGE incumbent, since it has no access to the stage's
partial state.

Measured before the fix, on an eight-primer panel with a twenty-evaluation
allowance: two removals were accepted and both were thrown away, so the run
delivered eight primers having found six. That is the plan's "return valid
incumbents at budget boundaries, including improvements already accepted within
a stage".

The stop reason distinguishes the two exhaustions on purpose. `shared_search_allowance`
says the RUN is out of compute; `evaluation_budget` says this loop reached its
own per-stage bound while the run had more to give. A reader deciding whether to
raise a limit needs to know which.
"""

from dataclasses import dataclass

import pytest

from neoswga.core import optimization_service as svc
from neoswga.core.search_control import SearchBudget, SearchBudgetExhausted

PANEL = tuple(f"P{i}" for i in range(8))


class Objective:
    """Every single removal qualifies, so the loop always has a proposal."""

    _inner = None

    def __init__(self):
        self.evaluation_budget = None
        self.constraints = type("C", (), {"coverage_metric": "effective"})()
        self.calls = 0

    def coverage(self, panel):
        self.calls += 1
        if self.evaluation_budget is not None:
            self.evaluation_budget.consume()
        return 1.0 - 0.01 * (len(PANEL) - len(panel))

    def violations(self, panel):
        return []


@dataclass
class Result:
    primers: tuple
    metrics: object = None
    score: float = 0.0


class Optimizer:
    config = type(
        "Cfg",
        (),
        {
            "swap_max_evaluations": 10_000,
            "swap_max_seconds": 30.0,
            "max_dimer_bp": None,
            "max_self_dimer_bp": None,
            "max_dimer_dg": None,
        },
    )()
    conditions = None

    def compute_metrics(self, panel):
        return type("M", (), {"normalized_score": lambda self: 0.5})()


@pytest.fixture
def harness(monkeypatch):
    objective = Objective()
    monkeypatch.setattr(svc, "objective_for_optimizer", lambda o, constraints=None: objective)
    monkeypatch.setattr(svc, "panel_violations", lambda o, primers: [])
    monkeypatch.setattr(svc, "panel_rank", lambda obj, panel: -len(panel))
    return objective


def test_removals_already_accepted_survive_a_spent_allowance(harness):
    budget = SearchBudget(max_evaluations=20)
    harness.evaluation_budget = budget
    diagnostics = {}

    result = svc.reduce_result(
        Result(primers=PANEL), Optimizer(), target_coverage=0.5, diagnostics=diagnostics
    )

    assert len(result.primers) < len(PANEL), (
        "the stage accepted removals and then hit the allowance; returning the "
        "original panel discards work that was already valid"
    )
    assert set(result.primers).issubset(PANEL)
    assert diagnostics["stop_reason"] == "shared_search_allowance"


def test_the_exception_does_not_escape_the_stage(harness):
    """`execute` upstream returns the pre-stage incumbent on this exception, so
    letting it escape is what lost the work."""
    harness.evaluation_budget = SearchBudget(max_evaluations=3)
    try:
        svc.reduce_result(Result(primers=PANEL), Optimizer(), target_coverage=0.5)
    except SearchBudgetExhausted:  # pragma: no cover - the defect
        pytest.fail("reduce_result let the allowance exception escape")


def test_an_allowance_that_is_never_spent_reduces_fully(harness):
    """The fix must not stop the loop early when there is compute to spend."""
    harness.evaluation_budget = SearchBudget(max_evaluations=10_000)
    diagnostics = {}

    result = svc.reduce_result(
        Result(primers=PANEL), Optimizer(), target_coverage=0.5, diagnostics=diagnostics
    )

    assert len(result.primers) == 1, "every removal qualified, so it should reduce to one"
    assert diagnostics["stop_reason"] != "shared_search_allowance"


def test_the_per_stage_bound_is_named_apart_from_the_shared_one(harness):
    """Two exhaustions, two names. They tell a reader to do different things."""
    harness.evaluation_budget = None

    class Bounded(Optimizer):
        config = type(
            "Cfg",
            (),
            {
                "swap_max_evaluations": 5,
                "swap_max_seconds": 30.0,
                "max_dimer_bp": None,
                "max_self_dimer_bp": None,
                "max_dimer_dg": None,
            },
        )()

    diagnostics = {}
    svc.reduce_result(
        Result(primers=PANEL), Bounded(), target_coverage=0.5, diagnostics=diagnostics
    )
    assert diagnostics["stop_reason"] == "evaluation_budget"


def test_fixed_primers_are_kept_even_when_the_allowance_runs_out(harness):
    """A partial reduction must still honour what the request pinned."""
    harness.evaluation_budget = SearchBudget(max_evaluations=20)

    result = svc.reduce_result(
        Result(primers=PANEL),
        Optimizer(),
        target_coverage=0.5,
        fixed_primers=("P0", "P1"),
    )

    assert {"P0", "P1"}.issubset(set(result.primers))
