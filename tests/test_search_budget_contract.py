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

import ast
import pathlib

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


#: Calls that measure a whole candidate panel. A loop over any of these is
#: search work by another route, whatever the caller names the variable.
#:
#: `optimize` is deliberately NOT here. The ensemble loop and
#: `collect_alternative_sets` both call it in a loop, and both run the inner
#: search against the SHARED ledger -- which is why `collect_alternative_sets`
#: has a `SearchBudgetExhausted` clause to catch. Listing it would report two
#: accounted loops as unaccounted, and an allowlist entry for work the ledger
#: already sees teaches the wrong thing about what this file is for.
PANEL_EVALUATING_CALLS = frozenset({"compute_metrics", "compute_pool_metrics", "evaluate_panel"})

#: Every construct that runs its body more than once. A comprehension is a loop
#: and was missing: only `ast.For` and `ast.While` were walked, so
#: `[compute_metrics(p) for p in panels]` -- the natural way to write a scan --
#: was invisible to this check. Nothing in the package does that today, which
#: is the point: the gap was in the detector, not in the code it guards.
LOOP_NODES = (ast.For, ast.While, ast.ListComp, ast.SetComp, ast.DictComp, ast.GeneratorExp)


def _called_name(call):
    """The name a call site uses, whether attribute or bare.

    `getattr(call.func, "attr", None)` alone saw `self.compute_metrics(...)`
    and missed `compute_metrics(...)`, which is the same work reached through
    an import rather than through an object.
    """
    return getattr(call.func, "attr", None) or getattr(call.func, "id", None)


def panel_evaluations_inside_loops(tree, label):
    """Every panel evaluation in `tree` that a loop encloses, by function."""
    found = {}
    for function in ast.walk(tree):
        if not isinstance(function, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        enclosed = {
            node.lineno
            for loop in ast.walk(function)
            if isinstance(loop, LOOP_NODES)
            for node in ast.walk(loop)
            if hasattr(node, "lineno")
        }
        for call in ast.walk(function):
            if not isinstance(call, ast.Call):
                continue
            if _called_name(call) not in PANEL_EVALUATING_CALLS:
                continue
            if call.lineno in enclosed:
                found[f"{label}::{function.name}"] = call.lineno
    return found


def _compute_metrics_calls_inside_loops():
    """Every panel evaluation the package makes inside a loop.

    A call made once per stage, after the panel is decided, is final
    assessment and deliberately uncharged. One made inside a loop is search
    work, and search work the ledger cannot see is the thing this checks for.
    """
    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    found = {}
    for path in sorted(package.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:  # pragma: no cover - would fail elsewhere
            continue
        found.update(panel_evaluations_inside_loops(tree, path.relative_to(package).as_posix()))
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


# ---------------------------------------------------------------------------
# The detector itself, against source it is handed
# ---------------------------------------------------------------------------
#
# Driving the detector on synthetic source rather than mutating the package.
# A ratchet that has never been shown to fire is a ratchet nobody can trust,
# and the two gaps closed here were both in the detector rather than in the
# code it guards -- so measuring the package would have shown nothing either
# way, before or after.


def _detect(source):
    return sorted(panel_evaluations_inside_loops(ast.parse(source), "probe"))


def test_a_comprehension_counts_as_a_loop():
    """The gap. `ast.For` and `ast.While` were the only constructs walked, so
    the natural way to write a scan was invisible."""
    assert _detect(
        "def scan(panels):\n" "    return [self.compute_metrics(p) for p in panels]\n"
    ) == ["probe::scan"]


def test_a_generator_expression_counts_too():
    assert _detect(
        "def scan(panels):\n" "    return max(self.compute_metrics(p) for p in panels)\n"
    ) == ["probe::scan"]


def test_a_bare_call_counts_as_well_as_an_attribute():
    """Reached through an import rather than through an object; same work."""
    assert _detect(
        "def scan(panels):\n" "    for p in panels:\n" "        compute_metrics(p)\n"
    ) == ["probe::scan"]


def test_the_other_panel_evaluators_count():
    for name in sorted(PANEL_EVALUATING_CALLS):
        assert _detect(f"def scan(ps):\n    return [{name}(p) for p in ps]\n") == [
            "probe::scan"
        ], f"{name} is listed as a panel evaluation but is not detected"


def test_a_single_call_outside_a_loop_is_not_flagged():
    """Guard the guard. Final assessment runs once per stage and is
    deliberately uncharged; a detector flagging it would be unusable."""
    assert _detect("def assess(panel):\n    return self.compute_metrics(panel)\n") == []


def test_an_unrelated_call_in_a_loop_is_not_flagged():
    assert _detect("def scan(ps):\n    return [len(p) for p in ps]\n") == []


def test_optimize_in_a_loop_is_not_flagged():
    """Deliberate, and the reason is in `PANEL_EVALUATING_CALLS`.

    The ensemble loop and `collect_alternative_sets` both call `optimize` in a
    loop, and both charge the inner search to the shared ledger -- which is why
    the latter has a `SearchBudgetExhausted` clause to catch. Flagging them
    would put two accounted loops on an allowlist for work the ledger sees.
    """
    assert _detect("def run(os_):\n    return [o.optimize(c, n) for o in os_]\n") == []
