"""The alternative-set search must not hide a failed calculation.

`collect_alternative_sets` wrapped its optimizer call in `except Exception` and
logged at DEBUG, so every failure became "Alternative set search stopped" at a
level nobody sees by default. The user got fewer primer sets and no reason.

Twenty lines below, the ensemble path does the opposite and explains itself:

    except DesignError:
        # A failed calculation, a missing reference answer or an unsupported
        # model is not a property of this ensemble member: the next member
        # would compute the same quantity from the same data.
        raise

That reasoning applies here verbatim. If computing a second primer set raises
because a reference is missing or a model is unsupported, the first set was
computed from the same data and the same models. Reporting it as "the
alternative search stopped" describes the symptom and hides the cause.

Two code paths in one module, opposite policy, on the same exception class.
This file makes them agree.

`SearchBudgetExhausted` is the case that genuinely should stop the loop rather
than fail the run. It is deliberately outside the `DesignError` family, because
spending an allowance is a recorded stopping point and not a failure, and it
now reads as one instead of as an unnamed exception at debug level.
"""

import logging

import pytest

from neoswga.core.exceptions import ReferenceDataError

# Imported from its new home; `unified_optimizer` re-exports it, and
# `test_the_ensemble_path_and_this_one_agree` checks the re-export holds.
from neoswga.core.search_control import SearchBudgetExhausted, collect_alternative_sets

PRIMARY = ("ACCACAGATAGC", "GTTGTAGATGGA")
POOL = ["ACCACAGATAGC", "GTTGTAGATGGA", "ATCAGCAGACCA", "TTGACCATGCAT", "CCATGGATCCAA"]


class Result:
    def __init__(self, primers):
        self.primers = list(primers)


class Optimizer:
    """An optimizer whose second call raises whatever it was given."""

    def __init__(self, error=None):
        self.error = error
        self.calls = 0

    def optimize(self, candidates, target_size):
        self.calls += 1
        if self.error is not None:
            raise self.error
        return Result(candidates[:target_size])


def collect(optimizer, max_sets=3):
    return collect_alternative_sets(
        Result(PRIMARY), optimizer, POOL, target_size=2, max_sets=max_sets, max_iterations=4
    )


# ---------------------------------------------------------------------------
# A design failure is not an alternative-search failure
# ---------------------------------------------------------------------------


def test_a_design_failure_propagates():
    """The headline. It was swallowed and logged at debug."""
    optimizer = Optimizer(
        ReferenceDataError("position index", "not found", "Re-run `neoswga filter`.")
    )

    with pytest.raises(ReferenceDataError):
        collect(optimizer)


def test_the_ensemble_path_and_this_one_agree():
    """Asserted against the source because the two sit in one module and
    drifted apart precisely by nobody comparing them."""
    import inspect

    from neoswga.core import unified_optimizer

    # Still reachable by its old name: the function moved to `search_control`
    # on 2026-09-22, beside the ledger its loop evades, because
    # `unified_optimizer` was at its size budget.
    assert unified_optimizer.collect_alternative_sets is collect_alternative_sets

    source = inspect.getsource(collect_alternative_sets)

    assert "DesignError" in source, (
        "the alternative-set loop still catches every exception alike, while "
        "the ensemble path in the same module deliberately re-raises DesignError"
    )


# ---------------------------------------------------------------------------
# A spent allowance is a stopping point, not a failure
# ---------------------------------------------------------------------------


def test_a_spent_budget_stops_the_search_without_failing_the_run():
    """`SearchBudgetExhausted` is deliberately outside the DesignError family:
    spending an allowance is a recorded stopping point."""
    sets = collect(Optimizer(SearchBudgetExhausted("evaluations spent")))

    assert sets == [tuple(PRIMARY)], "the primary set must survive a spent budget"


def test_a_spent_budget_is_reported_as_one(caplog):
    with caplog.at_level(logging.INFO, logger="neoswga.core.search_control"):
        collect(Optimizer(SearchBudgetExhausted("evaluations spent")))

    assert "budget" in caplog.text.lower() or "allowance" in caplog.text.lower()


# ---------------------------------------------------------------------------
# Anything else is visible
# ---------------------------------------------------------------------------


def test_an_unexpected_failure_is_not_logged_at_debug(caplog):
    """Debug is invisible at default verbosity, so the user saw fewer sets and
    no reason at all."""
    with caplog.at_level(logging.DEBUG, logger="neoswga.core.search_control"):
        collect(Optimizer(RuntimeError("something specific went wrong")))

    records = [r for r in caplog.records if "something specific" in r.getMessage()]
    assert records, "the failure was not reported at any level"
    assert all(
        r.levelno >= logging.WARNING for r in records
    ), "reported below WARNING, so a default run shows the user nothing"


# ---------------------------------------------------------------------------
# The normal path is unchanged
# ---------------------------------------------------------------------------


def test_alternatives_are_still_collected():
    """Guard the guard: a loop that always raises also passes the tests above."""
    sets = collect(Optimizer())

    assert len(sets) > 1
    assert sets[0] == tuple(PRIMARY)


def test_one_set_requested_means_no_search():
    optimizer = Optimizer()

    assert collect(optimizer, max_sets=1) == [tuple(PRIMARY)]
    assert optimizer.calls == 0
