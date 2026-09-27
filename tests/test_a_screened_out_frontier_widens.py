"""A screen emptying one frontier must not refuse a design over unseen candidates.

`prepare` raises `NoCandidatesError` when a screen rejects every candidate in a
non-empty frontier. That is right when there is nothing else to reach: it names
the screen and its threshold instead of letting an optimizer complain about an
empty list. It is wrong while candidates remain unexamined. On the Wolbachia
design a frontier is about 2,000 candidates of some 492,000, so one emptied
frontier says nothing about the next.

The refusal is kept for the case it was written for, and is raised from the
WIDEST frontier reached rather than the first, so its message names the screen
that emptied everything available.
"""

import pytest

from neoswga.core.exceptions import NoCandidatesError
from neoswga.core.search_control import search_frontiers


class Source:
    """A candidate source whose frontier widens through a fixed sequence."""

    def __init__(self, frontiers):
        self._frontiers = list(frontiers)
        self.index = 0
        self.advances = 0

    def frontier(self):
        return list(self._frontiers[self.index])

    def advance(self, keep=()):
        if self.index + 1 >= len(self._frontiers):
            return False
        self.index += 1
        self.advances += 1
        return True

    def describe(self):
        return {"frontier": self.index}


def _run(source, screened, *, max_refills=4, qualify_on=None):
    """Run the loop with a `prepare` that empties the named frontiers."""

    def prepare(sequences):
        pool = list(sequences)
        if pool and tuple(pool) in screened:
            raise NoCandidatesError(len(pool), "the self-dimer screen at max_self_dimer_bp=4")
        return pool

    def attempt(pool):
        return {"pool": list(pool)}

    def assess(outcome):
        pool = outcome["pool"] if outcome else []
        qualified = bool(pool) and (qualify_on is None or tuple(pool) == qualify_on)
        return qualified, (len(pool),), pool

    return search_frontiers(
        source, prepare(source.frontier()), attempt, assess, prepare, max_refills
    )


def test_an_emptied_frontier_widens_instead_of_refusing():
    """The first frontier is screened out; the second is usable."""
    source = Source([["A"], ["B", "C"], ["D", "E", "F"]])
    # The opening pool is prepared outside the loop, so make the SECOND
    # frontier the emptied one: that is the refill path this fix is about.
    _result, pool, history = _run(source, screened={("B", "C")}, qualify_on=("D", "E", "F"))

    assert pool == ["D", "E", "F"], "the loop stopped at a screened-out frontier"
    emptied = [r for r in history["frontier_attempts"] if r.get("screened_out")]
    assert len(emptied) == 1, "the widening past an emptied frontier must be recorded"
    assert "self-dimer screen" in emptied[0]["screened_out"]
    assert history["stop_reason"] == "qualified"


def test_the_refusal_survives_when_nothing_else_can_be_reached():
    """The case the refusal was written for, unchanged."""
    source = Source([["A"], ["B", "C"]])

    # Nothing qualifies, so the loop refills into the screened frontier and
    # then has nowhere left to widen.
    with pytest.raises(NoCandidatesError, match="self-dimer screen"):
        _run(source, screened={("B", "C")}, qualify_on=("unreachable",))


def test_several_emptied_frontiers_are_each_recorded_then_refused():
    source = Source([["A"], ["B"], ["C"], ["D"]])

    with pytest.raises(NoCandidatesError):
        _run(source, screened={("B",), ("C",), ("D",)}, qualify_on=("unreachable",))

    assert source.index == 3, "every reachable frontier should have been tried"


def test_a_frontier_that_prepares_cleanly_is_untouched():
    """No behaviour change where no screen empties anything."""
    source = Source([["A"], ["B", "C"]])
    _result, pool, history = _run(source, screened=set(), qualify_on=("B", "C"))

    assert pool == ["B", "C"]
    assert not [r for r in history["frontier_attempts"] if r.get("screened_out")]


def test_widening_past_an_empty_frontier_does_not_consume_extra_refills():
    """A refill is a widening the POLICY allowed; rescuing a screened-out
    frontier inside it must not be charged as another one, or a screen would
    silently shrink how far a run may look."""
    source = Source([["A"], ["B"], ["C", "D"]])
    _result, pool, history = _run(source, screened={("B",)}, max_refills=1, qualify_on=("C", "D"))

    assert pool == ["C", "D"], "one refill should still reach past an emptied frontier"
    assert history["frontier_refills"] == 1
