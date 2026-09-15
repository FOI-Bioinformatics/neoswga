"""Search and acceptance must mean the same thing by "coverage".

Task 5 of the condition-aware pool design plan.

The optimizers selected on one quantity and the report accepted on another. A
greedy maximising raw window coverage can prefer a panel that the
occupancy-weighted metric scores lower, because occupancy weights each site by
how much of the time it is actually bound. When the number a design is accepted
on is not the number it was chosen on, improving the search does not reliably
improve the result.

`PoolObjective` is the single contract. It wraps whatever evaluator computes a
panel's metrics and answers three questions with one set of definitions: what
the metrics are, which coverage figure counts, and which constraints a panel
violates.

The dimer guarantee is deliberately NOT folded in. It is a hard constraint on
the delivered panel enforced separately, and burying it among scoring terms is
how it came to be tradeable in the first place.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.pool_objective import PoolConstraints, PoolObjective

A, C, G = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"


def _metrics(raw, effective, density=100.0, bg_sites=0):
    return SimpleNamespace(
        fg_coverage=raw,
        effective_fg_coverage=effective,
        selectivity_density=density,
        total_bg_sites=bg_sites,
    )


def _objective(table, **limits):
    return PoolObjective(
        evaluate=lambda primers: table[tuple(sorted(primers))],
        constraints=PoolConstraints(**limits),
    )


def test_greater_raw_reach_can_score_lower_once_occupancy_is_weighted():
    """The disagreement the shared contract exists to prevent."""
    table = {
        (A, C): _metrics(raw=0.90, effective=0.40),
        (A, G): _metrics(raw=0.70, effective=0.65),
    }
    objective = _objective(table)

    assert objective.coverage([A, C]) < objective.coverage([A, G])
    assert table[(A, C)].fg_coverage > table[(A, G)].fg_coverage


def test_the_raw_metric_is_available_when_it_is_what_was_asked_for():
    table = {(A, C): _metrics(raw=0.90, effective=0.40)}
    objective = _objective(table, coverage_metric="raw")

    assert objective.coverage([A, C]) == 0.90


def test_equal_coverage_is_separated_by_background_load():
    table = {
        (A, C): _metrics(raw=0.8, effective=0.8, bg_sites=5),
        (A, G): _metrics(raw=0.8, effective=0.8, bg_sites=500),
    }
    objective = _objective(table, max_background_sites=100)

    assert objective.violations([A, C]) == ()
    assert "background sites above maximum" in objective.violations([A, G])


def test_a_selectivity_floor_is_enforced():
    table = {(A, C): _metrics(raw=0.8, effective=0.8, density=2.0)}
    objective = _objective(table, min_selectivity_density=10.0)

    assert "selectivity below minimum" in objective.violations([A, C])


def test_an_unavailable_coverage_is_a_violation_not_a_zero():
    """`None` means not measured, and must not be read as a bad panel."""
    table = {(A, C): _metrics(raw=0.8, effective=None)}
    objective = _objective(table)

    assert objective.coverage([A, C]) is None
    assert "coverage unavailable" in objective.violations([A, C])


def test_the_objective_and_the_report_read_the_same_metrics_object():
    """Acceptance must not recompute what search already decided on."""
    table = {(A, C): _metrics(raw=0.8, effective=0.7)}
    objective = _objective(table)

    metrics = objective.metrics([A, C])
    assert objective.coverage([A, C]) == metrics.effective_fg_coverage


def test_the_evaluator_is_called_once_per_distinct_panel():
    """Memoised, because search asks for the same panel repeatedly."""
    calls = []

    def evaluate(primers):
        calls.append(tuple(primers))
        return _metrics(raw=0.5, effective=0.5)

    objective = PoolObjective(evaluate=evaluate, constraints=PoolConstraints())
    objective.metrics([A, C])
    objective.coverage([A, C])
    objective.violations([A, C])

    assert len(calls) == 1


def test_panel_identity_ignores_the_order_it_was_built_in():
    calls = []

    def evaluate(primers):
        calls.append(tuple(primers))
        return _metrics(raw=0.5, effective=0.5)

    objective = PoolObjective(evaluate=evaluate, constraints=PoolConstraints())
    objective.metrics([A, C])
    objective.metrics([C, A])

    assert len(calls) == 1, "the same panel in a different order is the same panel"


def test_specificity_limits_require_a_background():
    """A limit on a quantity nothing measured cannot be satisfied or refused."""
    with pytest.raises(ValueError, match="background"):
        PoolConstraints(min_selectivity_density=10.0).require_background(available=False)

    PoolConstraints(min_selectivity_density=10.0).require_background(available=True)
    PoolConstraints().require_background(available=False)


def test_constraints_are_immutable():
    """A design must not have its acceptance criteria changed under it."""
    constraints = PoolConstraints(max_background_sites=10)
    with pytest.raises(Exception):
        constraints.max_background_sites = 99


def test_plan_pool_accepts_through_the_shared_objective():
    """Acceptance must not restate the rules in its own words.

    `plan_pool` had its own copy of the coverage choice and both specificity
    limits. Two copies of a rule drift, and the drift is invisible because each
    is self-consistent.
    """
    import ast
    import pathlib

    source = pathlib.Path("neoswga/core/pool_planner.py").read_text()
    tree = ast.parse(source)

    uses_objective = any(
        isinstance(node, ast.Name) and node.id in {"PoolObjective", "PoolConstraints"}
        for node in ast.walk(tree)
    )
    assert uses_objective, "plan_pool does not build a PoolObjective"

    plan = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "plan_pool"
    )
    # Reading a metric to RECORD it is fine. What must not be duplicated is the
    # RULE: comparing a panel's measurement against a configured limit, which is
    # the judgement the objective now owns.
    restated = [
        node
        for node in ast.walk(plan)
        if isinstance(node, ast.Compare)
        and any(
            isinstance(side, ast.Name)
            and side.id in {"min_selectivity_density", "max_background_sites"}
            for side in [node.left, *node.comparators]
        )
        and any(
            isinstance(side, ast.Attribute)
            and side.attr in {"selectivity_density", "total_bg_sites"}
            for side in [node.left, *node.comparators]
        )
    ]
    assert not restated, (
        "plan_pool still compares a panel metric against a configured limit "
        "itself instead of asking the objective which limits it violates"
    )


def test_the_dimer_guard_stays_outside_the_objective():
    """It is a hard constraint on the delivered panel, not a scoring term.

    Folding it in among the others is how it became tradeable, and the
    relaxation that followed produced an 11 bp heterodimer against a
    configured 3.
    """
    import pathlib

    source = pathlib.Path("neoswga/core/pool_objective.py").read_text()
    assert "incompatible_pairs" not in source
    assert "has_self_dimer" not in source

    planner = pathlib.Path("neoswga/core/pool_planner.py").read_text()
    assert "incompatible_pairs" in planner, "the panel dimer guard must still run"
