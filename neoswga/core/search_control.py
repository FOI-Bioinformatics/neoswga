"""Shared search limits and candidate-frontier traversal."""

import math
import time
from contextlib import contextmanager
from dataclasses import dataclass, field

from .exceptions import InvalidDesignRequest

SEARCH_CONTROL_DEFAULTS = {
    "total_search_evaluations": None,
    "total_search_seconds": None,
    "max_frontier_refills": 4,
}


class SearchBudgetExhausted(RuntimeError):
    """The shared search allowance has been spent."""


@dataclass
class SearchBudget:
    """Count uncached shared-objective evaluations across stages and refills.

    Time limits are cooperative: checks occur before stages and evaluations.
    An opaque proposal generator or an in-flight evaluation may overrun them.
    Final assessment/reporting is deliberately outside the search allowance.
    """

    max_evaluations: int | None = None
    max_seconds: float | None = None
    evaluations: int = 0
    started: float = field(default_factory=time.monotonic)
    stop_reason: str | None = None

    def __post_init__(self):
        if self.max_evaluations is not None and (
            isinstance(self.max_evaluations, bool)
            or not isinstance(self.max_evaluations, int)
            or self.max_evaluations < 0
        ):
            raise ValueError("total_search_evaluations must be a non-negative integer or null")
        if self.max_seconds is not None and (
            not math.isfinite(self.max_seconds) or self.max_seconds < 0
        ):
            raise ValueError("total_search_seconds must be finite and non-negative or null")

    @classmethod
    def from_config(cls, config):
        return cls(
            getattr(config, "total_search_evaluations", None),
            getattr(config, "total_search_seconds", None),
        )

    def check(self):
        if self.max_evaluations is not None and self.evaluations >= self.max_evaluations:
            self.stop_reason = "total_evaluation_budget"
        elif self.max_seconds is not None and time.monotonic() - self.started >= self.max_seconds:
            self.stop_reason = "total_time_budget"
        if self.stop_reason:
            raise SearchBudgetExhausted(self.stop_reason)

    def consume(self):
        self.check()
        self.evaluations += 1

    def describe(self):
        return dict(
            evaluations=self.evaluations,
            max_evaluations=self.max_evaluations,
            max_seconds=self.max_seconds,
            stop_reason=self.stop_reason,
            evaluation_scope="uncached_shared_objective",
            time_enforcement="cooperative",
        )


@contextmanager
def budgeted_objective(objective, budget):
    """Bind one budget to the evaluator, including a composed deficit objective."""
    base = objective
    while getattr(base, "_inner", None) is not None:
        base = base._inner
    if base is None:
        yield
        return
    previous = getattr(base, "evaluation_budget", None)
    base.evaluation_budget = budget
    try:
        yield
    finally:
        base.evaluation_budget = previous


def search_frontiers(source, pool, attempt, assess, prepare, max_refills, budget=None):
    """Widen on unmet requests; preserve the best incumbent across attempts.

    `assess` returns (qualified, rank, primers). Front ends own only preparation
    and result shape; stopping, incumbent retention and accounting live here.
    """
    best = None
    best_rank = None
    records = []
    refills = 0
    reason = "qualified"
    while True:
        if budget is not None:
            try:
                budget.check()
            except SearchBudgetExhausted:
                reason = budget.stop_reason
                break
        started = time.monotonic()
        outcome = attempt(pool)
        qualified, rank, primers = assess(outcome)
        records.append(
            dict(
                frontier_size=len(pool),
                qualified=qualified,
                primers=list(primers),
                seconds=time.monotonic() - started,
            )
        )
        if best is None or rank > best_rank:
            best, best_rank = outcome, rank
        if qualified:
            break
        if budget is not None:
            try:
                budget.check()
            except SearchBudgetExhausted:
                reason = budget.stop_reason
                break
        if refills >= max_refills:
            reason = "refill_budget"
            break
        keep = assess(best)[2]
        if source is None or not source.advance(keep=keep):
            reason = "inventory_exhausted"
            break
        pool = prepare(source.frontier())
        refills += 1
    return (
        best,
        pool,
        dict(
            frontier_refills=refills,
            frontier_attempts=records,
            stop_reason=reason,
            candidate_source=source.describe() if source else None,
            search_budget=budget.describe() if budget else None,
        ),
    )


#: Search settings resolved here rather than in `parameter.py`, with the
#: default that applies when params.json is silent. Every one of them bounds or
#: steers the panel search, so they belong beside the budget they spend.
SEARCH_SETTING_DEFAULTS = {
    "total_search_evaluations": None,
    "total_search_seconds": None,
    "max_frontier_refills": 4,
    "allow_dimer_relaxation": False,
    "refinement_method": "network",
    "swap_max_evaluations": 10000,
    "swap_max_seconds": 10.0,
    "stage1_objective_width": None,
    "coverage_reach": None,
}


def resolve_search_settings(data):
    """Resolve every search setting once, from one table, with validation.

    Extracted from `parameter._apply_params_only_keys` on 2026-09-21. Two
    reasons, one of them a defect.

    The defect: three of these keys were assigned through `globals()[key] =`
    inside a loop. That binds the global at run time, but
    `tests/test_no_schema_key_is_inert.py` reads `global` declarations from the
    source, so a key bound that way is invisible to the ratchet that exists to
    catch inert keys -- and Known Issue 8 is a list of keys that were inert
    while looking wired. A binding the ratchet cannot see is the shape of the
    next entry on that list, whether or not it works today.

    The second reason is that an out-of-range budget should be refused where
    the budget is defined. `SearchBudget.__post_init__` already validates the
    two evaluation limits, but it runs when a search starts, which is after a
    long preparation step; the same check here fails the run at load time.

    Returns a new dict. The caller assigns the globals explicitly, so the
    names stay greppable.
    """
    resolved = {key: data.get(key, default) for key, default in SEARCH_SETTING_DEFAULTS.items()}

    refills = resolved["max_frontier_refills"]
    if isinstance(refills, bool) or not isinstance(refills, int) or refills < 0:
        raise InvalidDesignRequest(
            "max_frontier_refills", "must be a non-negative integer", refills
        )

    # Raises `ValueError` for a negative or non-finite allowance. Constructing
    # it is the validation; the instance is discarded because the search builds
    # its own from the same values.
    try:
        SearchBudget(resolved["total_search_evaluations"], resolved["total_search_seconds"])
    except ValueError as exc:
        raise InvalidDesignRequest(
            "total_search_evaluations/total_search_seconds", str(exc)
        ) from exc

    return resolved
