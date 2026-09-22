"""A search may stop at the first feasible panel, or keep looking.

`search_frontiers` broke on the first qualifying attempt, unconditionally:

    if qualified:
        break

So refills existed only to rescue a search that had NOT qualified, and a run
qualifying on its opening frontier never looked wider. On the Wolbachia design
that is about 2,000 candidates examined of some 492,000 that cleared hard QC.

That is a defensible default for a fast feasible mode. What it is not is a
claim that the best small pool has been sought, and the code said nothing
either way.

`stop_on_first_qualified` names the choice. True is exactly the old behaviour
and remains the default: shipping the other would move every delivered panel,
and this project does not change a default on reasoning.

The property that makes the option safe is that `best` is already retained
across attempts by rank, so continuing can only return the same panel or a
better-ranked one. A test below pins that, because it is the whole argument.
"""

import pytest

from neoswga.core.search_control import SearchBudget, search_frontiers


class Source:
    """A frontier that can widen a fixed number of times."""

    def __init__(self, widenings):
        self.remaining = widenings
        self.width = 1

    def advance(self, keep=()):
        if self.remaining <= 0:
            return False
        self.remaining -= 1
        self.width += 1
        return True

    def frontier(self):
        return list(range(self.width))

    def describe(self):
        return {"universe": 99, "examined": self.width}


def run(ranks, *, stop_on_first_qualified=True, qualifies_from=0, max_refills=5, budget=None):
    """Drive the loop over a sequence of attempt ranks.

    `ranks[i]` is the rank of attempt i, so a later attempt can be better or
    worse than an earlier one and the incumbent logic is exercised either way.
    """
    seen = {"n": 0}

    def attempt(pool):
        index = min(seen["n"], len(ranks) - 1)
        seen["n"] += 1
        return index

    def assess(outcome):
        index = outcome if isinstance(outcome, int) else 0
        return (index >= qualifies_from, ranks[min(index, len(ranks) - 1)], [index])

    best, _pool, record = search_frontiers(
        Source(max_refills),
        [0],
        attempt,
        assess,
        lambda frontier: list(frontier),
        max_refills,
        budget=budget,
        stop_on_first_qualified=stop_on_first_qualified,
    )
    return best, record, seen["n"]


# ---------------------------------------------------------------------------
# The default is unchanged
# ---------------------------------------------------------------------------


def test_the_default_stops_at_the_first_qualifying_attempt():
    """Exactly the old behaviour, and what every shipped design still does."""
    best, record, attempts = run([1.0, 9.0, 9.0], qualifies_from=0)

    assert attempts == 1, "the default looked past the first qualifying attempt"
    assert best == 0
    assert record["stop_reason"] == "qualified"


def test_the_default_still_refills_when_nothing_qualifies():
    """Refills exist to rescue a search that has not qualified. Unchanged."""
    _best, record, attempts = run([1.0, 2.0, 3.0], qualifies_from=99, max_refills=2)

    assert attempts > 1
    assert record["frontier_refills"] == 2


# ---------------------------------------------------------------------------
# The quality policy
# ---------------------------------------------------------------------------


def test_the_quality_policy_looks_past_the_first_qualifying_attempt():
    _best, _record, attempts = run([1.0, 5.0, 9.0], stop_on_first_qualified=False, max_refills=2)

    assert attempts > 1


def test_the_quality_policy_returns_the_better_panel_when_one_exists():
    """The reason to offer it at all."""
    best, _record, _attempts = run([1.0, 5.0, 9.0], stop_on_first_qualified=False, max_refills=2)

    assert best == 2, "kept the first qualifying panel despite a better one later"


def test_the_quality_policy_never_returns_a_worse_panel():
    """The whole argument for the option being safe.

    `best` is retained across attempts by rank, so continuing can only return
    the same panel or a better-ranked one. Ranks that get WORSE after the first
    qualifying attempt are the case that would break it.
    """
    first, _r, _a = run([9.0, 5.0, 1.0], qualifies_from=0)
    kept, _r2, _a2 = run([9.0, 5.0, 1.0], qualifies_from=0, stop_on_first_qualified=False)

    assert kept == first, "looking further lost the better panel it started with"


def test_the_stop_reason_says_it_qualified_before_it_ran_out():
    """A reader must not read "refill_budget" as "never found one"."""
    _best, record, _attempts = run(
        [1.0, 5.0, 9.0], qualifies_from=0, stop_on_first_qualified=False, max_refills=2
    )

    assert "qualif" in record["stop_reason"], record["stop_reason"]


def test_the_quality_policy_respects_the_allowance():
    """It keeps looking while the budget permits, not indefinitely."""
    budget = SearchBudget(max_evaluations=1)
    budget.consume()

    _best, record, attempts = run(
        [1.0, 5.0, 9.0], stop_on_first_qualified=False, budget=budget, max_refills=5
    )

    assert attempts <= 2, "the allowance did not bound the quality policy"
