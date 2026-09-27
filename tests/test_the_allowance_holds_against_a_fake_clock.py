"""The time limit and the evaluation limit, pinned deterministically.

`total_search_seconds` is enforced by reading `time.monotonic`, so a test that
waited for real time would be slow and flaky -- which is why the plan asks for a
fake clock. A test that sleeps is measuring the machine.

The limit is COOPERATIVE and that word carries weight: checks happen before
stages and before evaluations, so an opaque proposal generator or an evaluation
already in flight overruns it. `describe()` says so rather than letting a reader
assume a wall-clock guarantee, and one test below pins exactly that, because the
honest limit and a hard one behave differently on the case that matters.
"""

import pytest

from neoswga.core import search_control
from neoswga.core.search_control import SearchBudget, SearchBudgetExhausted


class Clock:
    """A monotonic clock that only moves when told to."""

    def __init__(self, start=1000.0):
        self.now = start

    def monotonic(self):
        return self.now

    def advance(self, seconds):
        self.now += seconds


@pytest.fixture
def clock(monkeypatch):
    fake = Clock()
    # The module reference, not the global `time` module, so nothing else in the
    # process sees a frozen clock.
    monkeypatch.setattr(search_control, "time", fake)
    return fake


# ---------------------------------------------------------------------------
# The time limit
# ---------------------------------------------------------------------------


def test_the_clock_starts_when_the_budget_is_built(clock):
    budget = SearchBudget(max_seconds=10.0)
    clock.advance(9.999)
    budget.check()  # inside the limit, so no exception
    assert budget.stop_reason is None


def test_reaching_the_limit_exactly_stops_the_search(clock):
    budget = SearchBudget(max_seconds=10.0)
    clock.advance(10.0)

    with pytest.raises(SearchBudgetExhausted):
        budget.check()
    assert budget.stop_reason == "total_time_budget"


def test_a_spent_clock_cannot_be_reset_by_a_later_stage(clock):
    """The whole point of one shared ledger."""
    budget = SearchBudget(max_seconds=5.0)
    clock.advance(6.0)
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()

    # A later stage asks again; the answer must not improve.
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()
    assert budget.evaluations == 0


def test_the_evaluation_limit_is_reported_when_both_are_exceeded(clock):
    """Two reasons are true; the reported one must be stable, not incidental."""
    budget = SearchBudget(max_evaluations=1, max_seconds=5.0)
    budget.consume()
    clock.advance(99.0)

    with pytest.raises(SearchBudgetExhausted):
        budget.check()
    assert budget.stop_reason == "total_evaluation_budget"


def test_no_time_limit_means_the_clock_is_never_consulted(clock):
    budget = SearchBudget(max_evaluations=3)
    clock.advance(10_000.0)
    budget.consume()
    assert budget.stop_reason is None


def test_an_evaluation_already_in_flight_is_not_interrupted(clock):
    """Cooperative, measured rather than assumed.

    The check happens BEFORE the work. A call that starts inside the limit and
    runs long finishes, and the budget reports the overrun only when next asked.
    A reader told "time limit" would expect otherwise, which is why `describe()`
    names the enforcement.
    """
    budget = SearchBudget(max_seconds=10.0)
    budget.consume()  # starts inside the limit
    clock.advance(50.0)  # the "evaluation" ran long

    assert budget.evaluations == 1, "the completed evaluation still counts"
    with pytest.raises(SearchBudgetExhausted):
        budget.check()


def test_the_enforcement_is_named_in_the_record(clock):
    described = SearchBudget(max_seconds=10.0).describe()
    assert (
        described["time_enforcement"] == "cooperative"
    ), "a reader must not read this as a hard wall-clock guarantee"
    assert described["max_seconds"] == 10.0


# ---------------------------------------------------------------------------
# The counted allowance is never exceeded
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("allowance", [0, 1, 2, 7])
def test_the_count_reaches_the_allowance_and_never_passes_it(allowance):
    """`consume` checks before incrementing, so the count tops out exactly."""
    budget = SearchBudget(max_evaluations=allowance)

    for _ in range(allowance):
        budget.consume()
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()

    assert budget.evaluations == allowance
    assert budget.describe()["evaluations"] == allowance


def test_every_termination_carries_a_reason():
    for budget, expected in (
        (SearchBudget(max_evaluations=0), "total_evaluation_budget"),
        (SearchBudget(max_seconds=0.0), "total_time_budget"),
    ):
        with pytest.raises(SearchBudgetExhausted):
            budget.check()
        assert budget.describe()["stop_reason"] == expected


def test_a_negative_or_non_finite_limit_is_refused():
    """A limit nobody can satisfy is a configuration error, not a stop reason."""
    with pytest.raises(ValueError, match="non-negative integer"):
        SearchBudget(max_evaluations=-1)
    with pytest.raises(ValueError, match="finite and non-negative"):
        SearchBudget(max_seconds=float("nan"))
    with pytest.raises(ValueError, match="finite and non-negative"):
        SearchBudget(max_seconds=-1.0)
    with pytest.raises(ValueError, match="non-negative integer"):
        SearchBudget(max_evaluations=True)


def test_the_start_instant_is_read_through_the_module(clock):
    """Guard the indirection every test above depends on.

    `default_factory=time.monotonic` captures the function object when the
    dataclass is created, so the start instant would come from the real clock
    while `check` reads the substituted one. The two time bases then differ by
    however long the process has been alive, the subtraction goes negative, and
    the limit never fires -- so every test above would pass vacuously while the
    feature was broken.

    A source check as well as a behavioural one, because the behavioural form is
    exactly what stops working if the indirection is removed.
    """
    import ast
    import pathlib

    budget = SearchBudget(max_seconds=1.0)
    assert budget.started == clock.now, "the start instant ignored the substituted clock"

    source = pathlib.Path("neoswga/core/search_control.py").read_text()
    tree = ast.parse(source)
    assert (
        "default_factory=time.monotonic)" not in source
    ), "capturing the function object makes the time limit untestable"
    assert "default_factory=lambda: time.monotonic()" in source
    assert isinstance(tree, ast.Module)
