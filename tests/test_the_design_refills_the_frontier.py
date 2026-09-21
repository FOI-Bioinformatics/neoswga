"""A row that cannot qualify asks for more candidates before giving up.

Phase 4 increment 5 of the plan for `docs/validation/pipeline_audit_2026-09-16/`,
the half that makes the retained inventory reachable.

`advance()` now widens the frontier and `exhaustion()` says whether there is
anything behind it, but a capability nothing calls is the defect this audit is
named after. This pins the caller: `plan_pool` refills when a size row fails
its constraints and the source still holds candidates, and stops refilling when
the inventory is genuinely spent.

Stub optimizer and stub source on purpose. What is under test is the decision to
refill, how many times, and what gets reported, not what a real search finds
with more candidates. That is a measurement and belongs on the real pool.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.candidate_source import FRONTIER_EXHAUSTED, INVENTORY_EXHAUSTED
from neoswga.core.pool_planner import plan_pool

# Distinct 12-mers with no long complementary run against each other.
POOL = [
    "AAAAAAAAAAAA",
    "CCCCCCCCCCCC",
    "ACACACACACAC",
    "AGAGAGAGAGAG",
    "ATATATATATAT",
    "CACACACACACA",
]


class _GrowingSource:
    """A frontier that starts narrow and can be widened on request."""

    kind = "inventory"

    def __init__(self, sequences, frontier):
        self._sequences = list(sequences)
        self._frontier = frontier
        self._examined = []
        self.advances = 0
        self.vouched = []

    def initial(self, limit=None):
        self._examined = self._sequences[: self._frontier]
        return list(self._examined)

    def frontier(self):
        return list(self._examined)

    def advance(self, keep=()):
        if len(self._examined) >= len(self._sequences):
            return False
        self.advances += 1
        self._frontier = min(len(self._sequences), max(self._frontier * 2, 1))
        self._examined = self._sequences[: self._frontier]
        return True

    def universe_size(self):
        return len(self._sequences)

    def examined(self):
        return len(self._examined)

    def exhausted(self):
        return self.examined() >= self.universe_size()

    def exhaustion(self):
        return INVENTORY_EXHAUSTED if self.exhausted() else FRONTIER_EXHAUSTED

    def attach_positions(self, cache):
        self.position_cache = cache

    def ensure_positions(self, sequences):
        self.vouched.append(tuple(sequences))

    def describe(self):
        return {
            "kind": self.kind,
            "universe": self.universe_size(),
            "frontier": self._frontier,
            "examined": self.examined(),
            "unexamined": self.universe_size() - self.examined(),
            "exhausted": self.exhausted(),
            "rank_key": "search_rank",
        }


class _Optimizer:
    """Qualifies only once the frontier reaches the candidate it needs.

    The panel it returns is the last primer of whatever pool it is handed, and
    only `CACACACACACA` clears the density floor. So the row fails until the
    frontier has grown to include it, which is exactly the situation a refill
    exists for.
    """

    name = "test"
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    conditions = object()
    # A cache must be present for the position check to run at all, which is
    # deliberate: see `_vet_frontier`. Its contents are never read here.
    cache = object()

    def __init__(self, config=None):
        self.config = config or OptimizerConfig(objective_scan_width=None)
        self.pools_seen = []

    def optimize(self, candidates, target_size):
        self.pools_seen.append(list(candidates))
        return SimpleNamespace(
            primers=[candidates[-1]], status=OptimizationStatus.PARTIAL, message=""
        )

    def compute_metrics(self, primers):
        good = primers == ["CACACACACACA"]
        return SimpleNamespace(
            effective_fg_coverage=0.9 if good else 0.6,
            fg_coverage=0.9 if good else 0.6,
            selectivity_density=50.0 if good else 1.0,
            total_bg_sites=5,
            max_gap=100,
        )


def _plan(source, optimizer=None, **kwargs):
    optimizer = optimizer or _Optimizer()
    return optimizer, plan_pool(
        optimizer,
        source,
        [1],
        [0.5],
        primer_length=12,
        min_selectivity_density=10.0,
        repair=False,
        **kwargs,
    )


# -- the refill happens ----------------------------------------------------


def test_a_failing_row_widens_the_frontier():
    """The wiring. Without it the retained candidates stay unreachable."""
    source = _GrowingSource(POOL, frontier=2)

    _plan(source)

    assert source.advances > 0


def test_the_refill_reaches_the_candidate_that_qualifies():
    """And the row then reports itself as eligible."""
    source = _GrowingSource(POOL, frontier=2)

    _, plan = _plan(source)

    row = plan["rows"][0]
    assert row["primers"] == ["CACACACACACA"]
    assert not row["failed_constraints"]


def test_a_qualifying_row_does_not_refill():
    """Refilling a design that already met its constraints is waste."""
    source = _GrowingSource(POOL, frontier=len(POOL))

    _plan(source)

    assert source.advances == 0


def test_every_widened_frontier_is_vouched_for():
    """Increment 3's position check has to run on the new candidates too.

    A refill admits candidates the cache was not built over, which is the case
    `ensure_positions` and `PositionCache.load` were written for.
    """
    source = _GrowingSource(POOL, frontier=2)

    _plan(source)

    assert len(source.vouched) == source.advances + 1, (
        f"{len(source.vouched)} frontiers vouched for across "
        f"{source.advances} refills; each widened frontier needs checking"
    )


def test_the_optimizer_sees_a_bigger_pool_each_time():
    source = _GrowingSource(POOL, frontier=2)

    optimizer, _ = _plan(source)

    sizes = [len(pool) for pool in optimizer.pools_seen]
    assert sizes == sorted(sizes), sizes
    assert sizes[-1] > sizes[0]


# -- and it stops ----------------------------------------------------------


def test_refilling_stops_when_the_inventory_is_spent():
    """No infinite loop on a design nothing in the universe can satisfy."""
    source = _GrowingSource(POOL[:3], frontier=1)

    _, plan = _plan(source)

    assert source.exhausted() is True
    assert plan["rows"][0]["failed_constraints"], "the fixture should not qualify"


def test_the_refill_count_is_reported_per_row():
    """A reader cannot otherwise tell a first answer from a fourth."""
    source = _GrowingSource(POOL, frontier=2)

    _, plan = _plan(source)

    assert plan["rows"][0]["frontier_refills"] == source.advances


def test_the_row_says_which_exhaustion_ended_it():
    """`inventory_exhausted` is a different fact from `frontier_exhausted`."""
    source = _GrowingSource(POOL[:3], frontier=1)

    _, plan = _plan(source)

    assert plan["rows"][0]["candidates_exhausted"] == INVENTORY_EXHAUSTED


def test_a_budget_of_zero_refills_preserves_the_old_behaviour():
    """So the change can be turned off, and the row says it was."""
    source = _GrowingSource(POOL, frontier=2)
    optimizer = _Optimizer(OptimizerConfig(objective_scan_width=None, max_frontier_refills=0))

    _, plan = _plan(source, optimizer=optimizer)

    assert source.advances == 0
    assert plan["rows"][0]["frontier_refills"] == 0
    assert plan["rows"][0]["candidates_exhausted"] == FRONTIER_EXHAUSTED


def _many_12mers(n):
    """Distinct 12-mers, so the planner's length filter keeps them."""
    letters = "ACGT"
    out = []
    for i in range(n):
        value, digits = i, []
        for _ in range(12):
            digits.append(letters[value % 4])
            value //= 4
        out.append("".join(digits))
    return out


def test_the_refill_budget_is_honoured():
    source = _GrowingSource(_many_12mers(4096), frontier=1)
    optimizer = _Optimizer(OptimizerConfig(objective_scan_width=None, max_frontier_refills=3))

    _plan(source, optimizer=optimizer)

    assert source.advances == 3


@pytest.mark.parametrize("budget", [-1, 2.5, "3"])
def test_a_refill_budget_that_is_not_a_non_negative_integer_is_refused(budget):
    with pytest.raises(ValueError):
        OptimizerConfig(max_frontier_refills=budget).validate()


def test_feasible_pool_below_coverage_target_refills():
    class CoverageOptimizer(_Optimizer):
        def compute_metrics(self, primers):
            metrics = super().compute_metrics(primers)
            metrics.selectivity_density = 20
            metrics.effective_fg_coverage = 0.95 if POOL[2] in primers else 0.8
            return metrics

    source = _GrowingSource(POOL[:3], frontier=1)
    plan = plan_pool(
        CoverageOptimizer(), source, [1], [0.9], min_selectivity_density=10, repair=False
    )
    assert plan["rows"][0]["coverage"] == 0.95
    assert source.advances == 2


def test_refill_keeps_better_feasible_incumbent():
    class CoverageOptimizer(_Optimizer):
        def compute_metrics(self, primers):
            metrics = super().compute_metrics(primers)
            metrics.selectivity_density = 20
            metrics.effective_fg_coverage = 0.95 if POOL[0] in primers else 0.8
            return metrics

    source = _GrowingSource(POOL[:3], frontier=1)
    plan = plan_pool(
        CoverageOptimizer(), source, [1], [0.99], min_selectivity_density=10, repair=False
    )
    row = plan["rows"][0]
    assert row["primers"] == [POOL[0]]
    assert row["coverage"] == 0.95
    assert len(row["frontier_attempts"]) == 3
