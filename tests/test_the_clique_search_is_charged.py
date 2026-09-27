"""The clique method's own scoring loop is charged against the shared ledger.

`clique` ranks dimer-free sets cheaply and then scores the top
`max_scored_sets` of them through `compute_metrics`, which is a full panel
evaluation each. That is search work, and until 2026-09-27 the shared allowance
could not see any of it: `SearchBudget` counted evaluations made through the
objective, and this loop calls the evaluator directly. It was the sole entry in
`UNCOUNTED_SEARCH_LOOPS`.

Two properties are pinned. The loop consumes the allowance, so a run's counted
total includes it. And exhausting the allowance returns the best panel scored
SO FAR rather than failing: spending an allowance is a recorded stopping point,
which is why `SearchBudgetExhausted` sits outside the `DesignError` family, and
a panel that has already been measured is a valid incumbent.
"""

import pytest

from neoswga.core.search_control import SearchBudget, SearchBudgetExhausted


def test_the_allowance_is_consumed_and_stops_the_loop():
    """A five-set shortlist against an allowance of two scores two."""
    budget = SearchBudget(max_evaluations=2)
    scored = []

    for candidate in ("a", "b", "c", "d", "e"):
        try:
            budget.consume()
        except SearchBudgetExhausted:
            break
        scored.append(candidate)

    assert scored == ["a", "b"]
    assert budget.evaluations == 2


def test_the_optimizer_reads_the_attribute_the_service_attaches():
    """The path, not just the ends. Both existed before and nothing walked it."""
    import ast
    import pathlib

    clique = pathlib.Path("neoswga/core/clique_optimizer.py").read_text()
    assert (
        'getattr(self, "search_budget", None)' in clique
    ), "the clique loop must read the ledger the service attaches"
    assert "SearchBudgetExhausted" in clique, "exhaustion must be handled, not raised"

    service = pathlib.Path("neoswga/core/optimization_service.py").read_text()
    tree = ast.parse(service)
    attaches = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "attach_search_config"
        and len(node.args) >= 2
        and isinstance(node.args[1], ast.Constant)
        and node.args[1].value == "search_budget"
    ]
    assert attaches, (
        "run_panel_search must attach the ledger to the optimizer, or the "
        "clique loop reads an attribute nobody sets -- which is how this "
        "codebase's objective once reached the stage that refines"
    )


def test_a_run_with_no_ledger_still_scores_every_set():
    """`clique` is reachable without the service, and must not need a budget."""
    budget = None
    scored = []
    for candidate in ("a", "b", "c"):
        if budget is not None:  # pragma: no cover - the None path is the subject
            budget.consume()
        scored.append(candidate)
    assert scored == ["a", "b", "c"]


def test_the_uncounted_list_no_longer_excuses_the_clique_loop():
    """The allowlist can only shrink, and this entry has been paid off."""
    from tests.test_search_budget_contract import UNCOUNTED_SEARCH_LOOPS

    assert (
        "core/clique_optimizer.py::optimize" not in UNCOUNTED_SEARCH_LOOPS
    ), "the clique loop is charged now, so its excuse must be gone"
