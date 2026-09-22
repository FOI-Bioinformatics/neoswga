"""Shared search limits and candidate-frontier traversal."""

import logging
import math
import time
from contextlib import contextmanager
from dataclasses import dataclass, field

from .exceptions import DesignError, InvalidDesignRequest

logger = logging.getLogger(__name__)

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
        """What was spent, and what this ledger does not claim to have bounded.

        `evaluation_scope` is narrow on purpose. Two kinds of work sit outside
        it and a reader comparing this count against a run's wall clock needs
        to know which:

        - **Proposal effort.** An optimizer that scores candidate panels
          through `compute_metrics` directly, rather than through the shared
          objective, is doing real search work this counter never sees. The
          `clique` method does exactly that, bounded by its own
          `max_scored_sets`.
        - **Final assessment.** One `compute_metrics` per stage once the panel
          is decided. Deliberately uncharged, so that reporting cannot consume
          the allowance a search needs.

        `uncounted_scopes` names both rather than leaving a reader to infer
        them from a count that is lower than they expected.
        """
        return dict(
            evaluations=self.evaluations,
            max_evaluations=self.max_evaluations,
            max_seconds=self.max_seconds,
            stop_reason=self.stop_reason,
            evaluation_scope="uncached_shared_objective",
            uncounted_scopes=("proposal_generation", "final_assessment"),
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


# ---------------------------------------------------------------------------
# Alternative sets
# ---------------------------------------------------------------------------
#
# Moved here from `unified_optimizer` on 2026-09-22, for two reasons that point
# the same way. That module sits against its size budget, and this loop runs a
# whole optimizer per iteration OUTSIDE the ledger below -- so it belongs beside
# the accounting it currently evades rather than in the module that dispatches
# optimizers. `unified_optimizer` re-exports it, so every existing importer is
# unaffected.


def collect_alternative_sets(
    primary, optimizer, candidates, target_size, max_sets=1, max_iterations=8
):
    """Up to `max_sets` distinct primer sets, best first.

    Alternatives are found by removing the primers already chosen from the
    candidate pool and running selection again, so each one is a genuinely
    different set rather than a reordering of the same oligos. `max_iterations`
    caps how many such attempts are made, which matters because a pool can run
    out of usable candidates long before `max_sets` is reached.

    Fewer than `max_sets` is a normal outcome, not a failure: a small pool
    simply cannot yield many disjoint sets. The primary result is always first.

    Both parameters were documented and inert -- `max_sets` reached only
    `search_context.BFSConfig`, which has no callers, and no optimizer reads
    `config.max_iterations`. The output format already anticipated this: the
    `set_index` column of step4_improved_df.csv was hardcoded to 0.

    Note `max_iterations` bounds the search for ALTERNATIVES only. Bounding the
    primary selection with it would cap how many primers a run can choose, so
    `iterations: 8` would quietly truncate a 96-oligo panel.
    """
    sets = [tuple(primary.primers)]
    if not primary.primers or max_sets <= 1:
        return sets

    seen = {frozenset(primary.primers)}
    remaining = [p for p in candidates if p not in set(primary.primers)]
    attempts = 0

    while len(sets) < max_sets and attempts < max(1, int(max_iterations)):
        attempts += 1
        if len(remaining) < target_size:
            break
        try:
            result = optimizer.optimize(remaining, target_size)
        except DesignError:
            # Identical policy, and for the identical reason, to the ensemble
            # path below: a failed calculation, a missing reference answer or
            # an unsupported model is not a property of THIS attempt. The
            # primary set was computed from the same data and the same models,
            # so reporting "the alternative search stopped" names the symptom
            # and hides the cause. The two sat in one module with opposite
            # policy on the same exception class until 2026-09-22.
            raise
        except SearchBudgetExhausted as exc:
            # Deliberately outside the DesignError family: spending an
            # allowance is a recorded stopping point, not a failure. The
            # primary set stands and the run continues.
            logger.info(f"Alternative set search stopped on its budget: {exc}")
            break
        except Exception as exc:
            # Anything else is unexpected here. It was logged at DEBUG, which
            # is invisible at default verbosity, so a user got fewer sets than
            # `max_sets` with no reason given anywhere they would look.
            logger.warning(f"Alternative set search stopped: {exc}")
            break

        chosen = tuple(result.primers)
        if not chosen or frozenset(chosen) in seen:
            break

        sets.append(chosen)
        seen.add(frozenset(chosen))
        remaining = [p for p in remaining if p not in set(chosen)]

    return sets
