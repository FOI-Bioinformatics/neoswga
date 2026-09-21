"""One allowance for the whole search, and it cannot be reset by a later stage.

Task 6 of the 2026-09-21 valid-design plan. A search runs through several
stages -- greedy selection, repair, swaps, deletion, frontier refills,
ensemble members, alternatives -- and each one that carries its own counter can
spend a full allowance after the previous one has already spent its own. The
user asked for a bound on the run, not on each stage of it.

Two counters exist and they answer different questions, which is fine as long
as they are not confused:

- `SearchBudget` is the SHARED ledger. It counts uncached evaluations of the
  shared objective, across every stage, and raises when it is spent.
- `swap_max_evaluations` is a PER-STAGE allowance inside the deletion and swap
  loops. It starts at zero every time one of them is entered.

The per-stage one is not a bug; bounding one loop is a reasonable thing to do.
What would be a bug is believing it bounds the run. `total_search_evaluations`
is the only thing that does, and it is None by default, so by default there is
no total bound at all. This file pins both halves so the distinction stays
visible.
"""

import pytest

from neoswga.core.search_control import (
    SEARCH_CONTROL_DEFAULTS,
    SearchBudget,
    SearchBudgetExhausted,
    budgeted_objective,
)

# ---------------------------------------------------------------------------
# The ledger itself
# ---------------------------------------------------------------------------


def test_one_allowance_cannot_be_reset_by_the_next_stage():
    budget = SearchBudget(max_evaluations=2)
    budget.consume()
    budget.consume()
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()
    assert budget.evaluations == 2


def test_the_count_does_not_advance_past_the_cap_on_a_refused_call():
    """A refused evaluation is not an evaluation.

    If `consume` incremented before raising, a caller catching the exception
    would see a count above the cap and the ledger would misreport what was
    actually spent.
    """
    budget = SearchBudget(max_evaluations=1)
    budget.consume()
    for _ in range(3):
        with pytest.raises(SearchBudgetExhausted):
            budget.consume()
    assert budget.evaluations == 1


def test_a_spent_budget_records_why_it_stopped():
    budget = SearchBudget(max_evaluations=1)
    budget.consume()
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()

    assert budget.stop_reason == "total_evaluation_budget"
    assert budget.describe()["stop_reason"] == "total_evaluation_budget"


def test_a_zero_allowance_permits_nothing():
    """Zero is a real request, not an absent one. `None` is the absent one."""
    budget = SearchBudget(max_evaluations=0)
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()
    assert budget.evaluations == 0


def test_no_allowance_is_unbounded_rather_than_zero():
    budget = SearchBudget(max_evaluations=None)
    for _ in range(50):
        budget.consume()
    assert budget.evaluations == 50
    assert budget.stop_reason is None


@pytest.mark.parametrize("bad", [-1, 1.5, True])
def test_a_nonsense_allowance_is_refused_at_construction(bad):
    """Caught where the budget is defined, not after a long preparation step."""
    with pytest.raises(ValueError, match="total_search_evaluations"):
        SearchBudget(max_evaluations=bad)


def test_the_time_limit_is_declared_cooperative():
    """An opaque proposal generator can overrun it, and the record says so.

    Claiming a hard wall-clock limit for a call this code cannot interrupt
    would be a promise the implementation does not keep.
    """
    described = SearchBudget(max_seconds=5.0).describe()

    assert described["time_enforcement"] == "cooperative"
    assert described["evaluation_scope"] == "uncached_shared_objective"


# ---------------------------------------------------------------------------
# What the ledger counts
# ---------------------------------------------------------------------------


class Metrics:
    fg_coverage = 0.5
    effective_fg_coverage = 0.5
    total_bg_sites = 10
    selectivity_density = 30.0
    max_gap = 100.0
    mean_gap = 50.0
    gap_gini = 0.2
    bg_coverage = 0.01

    def normalized_score(self, application="balanced"):
        return 0.5


def make_objective(counter):
    from neoswga.core.panel_acceptance import PoolConstraints
    from neoswga.core.pool_objective import PoolObjective

    def evaluate(primers):
        counter.append(tuple(primers))
        return Metrics()

    return PoolObjective(evaluate, PoolConstraints(min_selectivity_density=1.0))


def test_a_cache_hit_is_not_charged():
    """Asking about the same panel twice is one evaluation, not two.

    Search revisits a panel constantly and in different orders. Charging for a
    lookup would make the allowance a function of how the search happened to
    traverse rather than of how much work it did.
    """
    evaluated = []
    objective = make_objective(evaluated)
    budget = SearchBudget(max_evaluations=10)

    with budgeted_objective(objective, budget):
        for _ in range(5):
            objective.metrics(["AAAACCCCGGGG", "TTTTGGGGCCCC"])
            objective.metrics(["TTTTGGGGCCCC", "AAAACCCCGGGG"])

    assert len(evaluated) == 1, "the same panel in a different order is the same panel"
    assert budget.evaluations == 1


def test_distinct_panels_are_each_charged_once():
    evaluated = []
    objective = make_objective(evaluated)
    budget = SearchBudget(max_evaluations=10)

    with budgeted_objective(objective, budget):
        for index in range(4):
            objective.metrics([f"AAAACCCCGG{index:02d}"])

    assert len(evaluated) == 4
    assert budget.evaluations == 4


def test_the_ledger_stops_the_search_rather_than_the_stage():
    """The property the whole file is about, on a real objective.

    Two stages, one ledger. The second must not get a fresh allowance.
    """
    evaluated = []
    objective = make_objective(evaluated)
    budget = SearchBudget(max_evaluations=3)

    with budgeted_objective(objective, budget):
        for index in range(3):
            objective.metrics([f"AAAACCCCGG{index:02d}"])

        # A second stage, entered afterwards, with its own loop.
        with pytest.raises(SearchBudgetExhausted):
            for index in range(10, 20):
                objective.metrics([f"AAAACCCCGG{index:02d}"])

    assert budget.evaluations == 3
    assert len(evaluated) == 3


def test_the_binding_is_removed_when_the_scope_ends():
    """A budget left attached would charge a later, unrelated run."""
    objective = make_objective([])
    budget = SearchBudget(max_evaluations=5)

    with budgeted_objective(objective, budget):
        assert getattr(objective, "evaluation_budget", None) is budget
    assert getattr(objective, "evaluation_budget", None) is None


def test_the_binding_reaches_a_composed_objective():
    """A deficit objective wraps the shared one; the ledger must reach through.

    `budgeted_objective` walks `_inner` for exactly this reason. Binding the
    wrapper would leave the evaluations the inner one performs uncounted, which
    is the shape of the defect recorded in `attach_search_config`: both ends
    existed and the path did not.
    """
    evaluated = []
    inner = make_objective(evaluated)

    class Wrapper:
        def __init__(self, inner):
            self._inner = inner

    wrapper = Wrapper(inner)
    budget = SearchBudget(max_evaluations=5)

    with budgeted_objective(wrapper, budget):
        assert getattr(inner, "evaluation_budget", None) is budget
        inner.metrics(["AAAACCCCGGGG"])

    assert budget.evaluations == 1


# ---------------------------------------------------------------------------
# The per-stage allowance is a different thing, and says so
# ---------------------------------------------------------------------------


def test_there_is_no_total_bound_by_default():
    """Recorded rather than implied.

    `swap_max_evaluations` defaults to 10,000 and looks like a bound on the
    run. It is a bound on one loop, reset every time that loop is entered, so
    a run with several stages can spend several multiples of it. The only
    setting that bounds the run is `total_search_evaluations`, and it is None
    unless someone sets it.
    """
    assert SEARCH_CONTROL_DEFAULTS["total_search_evaluations"] is None
    assert SEARCH_CONTROL_DEFAULTS["total_search_seconds"] is None


def test_the_per_stage_allowance_starts_fresh_each_time():
    """Pinned so the distinction from the shared ledger stays visible.

    The deletion loop reads `swap_max_evaluations` off the optimizer config and
    starts its own counter at zero. Two calls therefore spend two allowances.
    That is correct for a per-loop bound and wrong for a run bound, and the
    difference is only safe while both are written down.
    """
    import ast
    import inspect

    from neoswga.core import optimization_service

    source = inspect.getsource(optimization_service)
    tree = ast.parse(source)

    resets = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Name) and target.id == "evaluations" for target in node.targets
        )
        and isinstance(node.value, ast.Constant)
        and node.value.value == 0
    ]

    assert resets, (
        "the per-stage counter no longer starts at zero; if it became shared "
        "state this file's description of the two counters is out of date"
    )


def test_a_stage_that_stops_on_its_own_allowance_records_which_one():
    """`evaluation_budget` and `total_evaluation_budget` are different stops.

    A reader told only that the search "hit its budget" cannot tell whether to
    raise the per-loop allowance or the run allowance.
    """
    from neoswga.core import optimization_service

    source = __import__("inspect").getsource(optimization_service)

    assert '"evaluation_budget"' in source
    assert SearchBudget(max_evaluations=0).describe()["stop_reason"] is None
    spent = SearchBudget(max_evaluations=0)
    with pytest.raises(SearchBudgetExhausted):
        spent.consume()
    assert spent.stop_reason == "total_evaluation_budget"


# ---------------------------------------------------------------------------
# What the ledger does not bound, named rather than inferred
# ---------------------------------------------------------------------------


#: Functions that evaluate candidate panels in a LOOP through
#: `compute_metrics` rather than through the shared objective, so the ledger
#: never sees them. Each needs a reason and its own bound. The list can only
#: shrink: a new entry means a new stage of the search became invisible to the
#: allowance the user set.
UNCOUNTED_SEARCH_LOOPS = {
    "core/clique_optimizer.py::optimize": (
        "scores the top `max_scored_sets` dimer-free cliques to pick a winner. "
        "Bounded by that setting rather than by the shared ledger; `clique` is "
        "not in the default ensemble and runs on pools of about 200 candidates."
    ),
}


def _compute_metrics_calls_inside_loops():
    """Every `compute_metrics` call that a `for` or `while` encloses.

    A call made once per stage, after the panel is decided, is final
    assessment and deliberately uncharged. One made inside a loop is search
    work, and search work the ledger cannot see is the thing this checks for.
    """
    import ast
    import pathlib

    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    found = {}
    for path in sorted(package.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:  # pragma: no cover - would fail elsewhere
            continue
        for function in ast.walk(tree):
            if not isinstance(function, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            enclosed = {
                node.lineno
                for loop in ast.walk(function)
                if isinstance(loop, (ast.For, ast.While))
                for node in ast.walk(loop)
                if hasattr(node, "lineno")
            }
            for call in ast.walk(function):
                if not isinstance(call, ast.Call):
                    continue
                if getattr(call.func, "attr", None) != "compute_metrics":
                    continue
                if call.lineno in enclosed:
                    key = f"{path.relative_to(package).as_posix()}::{function.name}"
                    found[key] = call.lineno
    return found


def test_no_new_search_loop_escapes_the_ledger():
    """The ratchet. A stage that evaluates panels in a loop must be accounted for.

    `SearchBudget` counts uncached evaluations of the SHARED objective. A loop
    calling `compute_metrics` directly does the same work by another route and
    is invisible to the allowance, so a user who set
    `total_search_evaluations` would find the run spending past it.
    """
    unaccounted = sorted(set(_compute_metrics_calls_inside_loops()) - set(UNCOUNTED_SEARCH_LOOPS))

    assert not unaccounted, (
        "these evaluate candidate panels in a loop without going through the "
        f"shared objective, so the search allowance cannot see them: {unaccounted}. "
        "Route them through the objective, or add an entry naming the reason and "
        "the bound that does apply."
    )


def test_the_allowlist_has_no_stale_entries():
    """A fix is not finished until it also removes its excuse."""
    in_loops = set(_compute_metrics_calls_inside_loops())
    stale = sorted(set(UNCOUNTED_SEARCH_LOOPS) - in_loops)

    assert not stale, f"no longer evaluate in a loop, so the excuse can go: {stale}"


def test_the_ledger_says_what_it_does_not_bound():
    """A count lower than a reader expects needs its scope beside it."""
    described = SearchBudget(max_evaluations=10).describe()

    assert "proposal_generation" in described["uncounted_scopes"]
    assert "final_assessment" in described["uncounted_scopes"]
