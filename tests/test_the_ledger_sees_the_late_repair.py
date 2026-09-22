"""Work done after the search still spends the search's allowance.

`SearchBudget` is described as the shared ledger: it counts uncached objective
evaluations across every stage so that a later stage cannot get a fresh
allowance. It reaches an objective through exactly one seam,
`budgeted_objective`, used at exactly one site inside `_run_panel_stages`.

`panel_acceptance.apply_configured_limits` runs after that returns. It calls
`repair_result`, which can spend up to `swap_max_evaluations` (10,000 by
default) refining the delivered panel, and none of it was counted. A run that
declared `total_search_evaluations` got that allowance for the search and then
a second, undeclared one afterwards.

The work was never unbounded -- `swap_max_evaluations` and `swap_max_seconds`
bound it -- so this is about the ledger telling the truth rather than about
runaway cost. `SearchBudget.describe()` names the two scopes it knowingly does
not see, `proposal_generation` and `final_assessment`, and this was not among
them: it was not a declared exemption, it was a gap.

**What this changes for a default run: nothing.** `total_search_evaluations`
defaults to None, so there is no allowance to exhaust and the repair runs
exactly as before. It changes a run that set one, which is the run that asked
to be bounded.
"""

import pytest

from neoswga.core.pool_objective import PoolConstraints, PoolObjective
from neoswga.core.search_control import SearchBudget, budgeted_objective


class Metrics:
    fg_coverage = 0.5
    effective_fg_coverage = 0.45
    bg_coverage = 0.01
    total_fg_sites = 100
    total_bg_sites = 5
    selectivity_density = 30.0


def objective_with(calls):
    def evaluate(primers):
        calls.append(tuple(primers))
        return Metrics()

    return PoolObjective(evaluate, PoolConstraints())


# ---------------------------------------------------------------------------
# The ledger's own contract
# ---------------------------------------------------------------------------


def test_an_evaluation_inside_the_binding_is_counted():
    """Guard the guard: the mechanism this file relies on."""
    budget = SearchBudget(max_evaluations=10)
    objective = objective_with([])

    with budgeted_objective(objective, budget):
        objective.metrics(["ACGTACGTACGT"])

    assert budget.evaluations == 1


def test_an_evaluation_outside_the_binding_is_not_counted():
    """The shape of the defect, in one line."""
    budget = SearchBudget(max_evaluations=10)
    objective = objective_with([])

    objective.metrics(["ACGTACGTACGT"])

    assert budget.evaluations == 0


def test_the_binding_is_removed_when_the_scope_ends():
    budget = SearchBudget(max_evaluations=10)
    objective = objective_with([])

    with budgeted_objective(objective, budget):
        pass
    objective.metrics(["ACGTACGTACGT"])

    assert budget.evaluations == 0


# ---------------------------------------------------------------------------
# The late repair
# ---------------------------------------------------------------------------


def test_the_late_repair_accepts_a_budget():
    """`apply_configured_limits` had no way to be told about the ledger.

    Asserted on the signature because the alternative is standing up a whole
    optimizer, and what went wrong was that the parameter did not exist at all.
    """
    import inspect

    from neoswga.core.panel_acceptance import apply_configured_limits

    assert "budget" in inspect.signature(apply_configured_limits).parameters, (
        "the late repair cannot be told about the shared ledger, so up to "
        "swap_max_evaluations of objective work after the search is invisible"
    )


def test_the_holder_passes_the_budget_through():
    """Both ends existing is not the same as the path existing.

    `attach_search_config` records what that costs here: two tests asserted the
    two ends of a connection that did not exist.
    """
    import ast
    import inspect

    from neoswga.core import unified_optimizer

    source = inspect.getsource(unified_optimizer._hold_to_configured_limits)
    tree = ast.parse(source.strip())

    passes_budget = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "apply_configured_limits"
        and any(kw.arg == "budget" for kw in node.keywords)
    ]

    assert passes_budget, "_hold_to_configured_limits does not forward the budget"


def test_the_ledger_no_longer_lists_the_repair_as_unseen():
    """`describe()` names the scopes the ledger knowingly does not see. The
    late repair was never among them: it was a gap, not a declared exemption,
    and it must not become one now."""
    budget = SearchBudget(max_evaluations=10)

    scopes = budget.describe().get("uncounted_scopes", ())

    assert "late_repair" not in scopes
    assert "configured_limits" not in scopes
