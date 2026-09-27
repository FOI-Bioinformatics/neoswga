"""Alternative primer sets spend the run's allowance, not a fresh one.

`budgeted_objective` binds the shared `SearchBudget` to the evaluator and
RESTORES the previous binding when it exits, which is correct for a context
manager and was the whole problem: `run_panel_search` returns, the binding is
undone, and `collect_alternative_sets` then ran each alternative with no budget
bound at all. A run that declared `total_search_evaluations` could spend far
past it, once per alternative, and its own `except SearchBudgetExhausted` clause
could never fire -- a handler for something that could not happen.

That clause is the tell. It was written believing the ledger was in force, which
is the shape this repository calls Known Issue 8: both ends exist and nothing
walks between them.
"""

import pytest

from neoswga.core.search_control import (
    SearchBudget,
    SearchBudgetExhausted,
    budgeted_objective,
    collect_alternative_sets,
)


class Objective:
    _inner = None

    def __init__(self):
        self.evaluation_budget = None


class Result:
    def __init__(self, primers):
        self.primers = tuple(primers)


class Optimizer:
    """Each `optimize` call spends one evaluation, as a real search would."""

    name = "fake"

    def __init__(self, objective):
        self.pool_objective = objective
        self.calls = 0

    def optimize(self, candidates, target_size):
        self.calls += 1
        budget = self.pool_objective.evaluation_budget
        if budget is not None:
            budget.consume()
        return Result(candidates[:target_size])


def test_the_binding_does_not_survive_the_context():
    """The mechanism behind the defect, pinned so the fix is not misread."""
    objective = Objective()
    with budgeted_objective(objective, SearchBudget(max_evaluations=5)):
        assert objective.evaluation_budget is not None
    assert objective.evaluation_budget is None


def test_alternatives_spend_the_shared_allowance():
    objective = Objective()
    optimizer = Optimizer(objective)
    budget = SearchBudget(max_evaluations=10)
    pool = [f"P{i}" for i in range(40)]

    collect_alternative_sets(
        primary=Result(pool[:4]),
        optimizer=optimizer,
        candidates=pool,
        target_size=4,
        max_sets=5,
        max_iterations=5,
        budget=budget,
    )

    assert budget.evaluations > 0, (
        "the alternative search spent nothing from the run's ledger, so a "
        "declared total does not bound it"
    )
    assert budget.evaluations == optimizer.calls


def test_an_exhausted_allowance_stops_the_alternatives_and_keeps_the_primary():
    """A later stage must not get a fresh allowance."""
    objective = Objective()
    optimizer = Optimizer(objective)
    budget = SearchBudget(max_evaluations=1)
    pool = [f"P{i}" for i in range(40)]

    sets = collect_alternative_sets(
        primary=Result(pool[:4]),
        optimizer=optimizer,
        candidates=pool,
        target_size=4,
        max_sets=5,
        max_iterations=5,
        budget=budget,
    )

    assert sets[0] == tuple(pool[:4]), "the primary set must always stand"
    assert optimizer.calls <= 2, "the search continued past its exhausted allowance"


def test_the_budget_clause_can_now_actually_fire():
    """Before the fix nothing raised here, so the handler was unreachable."""
    objective = Objective()
    optimizer = Optimizer(objective)
    spent = SearchBudget(max_evaluations=1)
    spent.consume()

    sets = collect_alternative_sets(
        primary=Result(["P0"]),
        optimizer=optimizer,
        candidates=[f"P{i}" for i in range(20)],
        target_size=2,
        max_sets=3,
        budget=spent,
    )

    assert sets == [("P0",)], "only the primary should survive an already-spent allowance"
    with pytest.raises(SearchBudgetExhausted):
        spent.check()


def test_no_budget_leaves_the_search_unconstrained():
    """An unconstrained run is the default and must be untouched."""
    objective = Objective()
    optimizer = Optimizer(objective)
    pool = [f"P{i}" for i in range(40)]

    sets = collect_alternative_sets(
        primary=Result(pool[:4]),
        optimizer=optimizer,
        candidates=pool,
        target_size=4,
        max_sets=3,
        max_iterations=3,
    )

    assert len(sets) > 1, "alternatives should still be found with no allowance set"
    assert objective.evaluation_budget is None


def test_the_call_site_passes_the_run_allowance():
    """A source check: the parameter is useless if nothing supplies it."""
    import pathlib

    source = pathlib.Path("neoswga/core/unified_optimizer.py").read_text()
    call = source[source.index("collect_alternative_sets(") :]
    call = call[: call.index(")\n")]
    assert (
        "budget=search_budget" in call
    ), "run_optimization must hand alternatives the run's ledger"
