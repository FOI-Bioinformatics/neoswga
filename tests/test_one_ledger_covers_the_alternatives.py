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
        # The bare route, named explicitly. This pins the exception POLICY of
        # `collect_alternative_sets` with a fake too small to survive the whole
        # service, and that policy is the same whichever route the search takes.
        # `test_the_default_route_is_the_shared_contract` pins what production uses.
        through_contract=False,
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
        # The bare route, named explicitly. This pins the exception POLICY of
        # `collect_alternative_sets` with a fake too small to survive the whole
        # service, and that policy is the same whichever route the search takes.
        # `test_the_default_route_is_the_shared_contract` pins what production uses.
        through_contract=False,
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
        # The bare route, named explicitly. This pins the exception POLICY of
        # `collect_alternative_sets` with a fake too small to survive the whole
        # service, and that policy is the same whichever route the search takes.
        # `test_the_default_route_is_the_shared_contract` pins what production uses.
        through_contract=False,
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
        # The bare route, named explicitly. This pins the exception POLICY of
        # `collect_alternative_sets` with a fake too small to survive the whole
        # service, and that policy is the same whichever route the search takes.
        # `test_the_default_route_is_the_shared_contract` pins what production uses.
        through_contract=False,
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


def test_the_default_route_is_the_shared_contract():
    """Production routes alternatives through `run_panel_search`, not a bare
    `optimize`.

    That is Task 6's "route ... alternatives through the same contract", and it
    is what makes the ledger binding real: `optimize` alone consults no
    objective, so a budget bound around it bounds nothing. Measured on the
    Wolbachia design, the same five alternatives cost 0 ledger evaluations
    through the bare route and 2,988 through the contract
    (docs/validation/alternatives_through_the_contract_2026-09-28.md).

    The tests above pass `through_contract=False` deliberately, so this is what
    stops the default drifting away from what production runs.
    """
    import inspect

    from neoswga.core.search_control import collect_alternative_sets

    default = inspect.signature(collect_alternative_sets).parameters["through_contract"].default
    assert default is True, "alternatives must take the shared contract by default"


def test_the_caller_hands_alternatives_the_run_constraints():
    """An alternative assessed against nothing is not assessed."""
    import pathlib

    source = pathlib.Path("neoswga/core/unified_optimizer.py").read_text()
    call = source[source.index("collect_alternative_sets(") :]
    call = call[: call.index(")\n")]
    assert "constraints=constraints_from_parameter(parameter)" in call
