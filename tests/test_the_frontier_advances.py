"""The frontier moves, and says which kind of exhaustion stopped it.

Phase 4 increment 5 of the plan for `docs/validation/pipeline_audit_2026-09-16/`.

`advance()` returned False from the day it was written. The docstring said so,
and called that honest, which it was as far as it went: the counts reported how
much of the universe went unexamined. But a source that never advances makes the
universe decorative. On the Wolbachia design 20,670 candidates are eligible,
2,000 were searched, and 18,670 could not affect any panel.

The other half is the vocabulary. A search that ran out of the candidates it was
handed reported `inventory_exhausted`, which is true of the batch and false of
the inventory. The two want different responses: refill and continue, against
stop because there is genuinely nothing left. So `frontier_exhausted` is the
answer when the source still holds candidates.
"""

import pytest

from neoswga.core.candidate_source import (
    FRONTIER_EXHAUSTED,
    INVENTORY_EXHAUSTED,
    InventoryCandidateSource,
    ListCandidateSource,
)


class _StubInventory:
    """An ordered eligible set. Ordering itself is tested elsewhere."""

    def __init__(self, sequences):
        self._sequences = list(sequences)
        self.queries = 0

    def current_policy(self, condition_id):
        return None

    def iter_eligible(self, condition_id, lengths, *rest):
        self.queries += 1
        return list(self._sequences)


def _source(universe=40, frontier=8):
    from neoswga.core.candidate_provider import CandidateProvider

    sequences = [f"P{i:03d}" for i in range(universe)]
    provider = CandidateProvider(_StubInventory(sequences), "cond", [12])
    return InventoryCandidateSource(provider, frontier=frontier), sequences


# -- advancing -------------------------------------------------------------


def test_the_frontier_starts_at_the_configured_size():
    source, sequences = _source(universe=40, frontier=8)

    assert source.initial() == sequences[:8]
    assert source.examined() == 8
    assert source.universe_size() == 40


def test_advance_widens_the_frontier_by_a_batch():
    """The behaviour the method promised and did not have."""
    source, sequences = _source(universe=40, frontier=8)
    source.initial()

    assert source.advance() is True
    assert source.frontier() == sequences[:16]


def test_advance_keeps_the_recorded_order():
    """The optimizers are order-sensitive, so a refill extends, not reshuffles.

    Asserted as a prefix rather than as a size, because the growth rule is a
    separate decision and is pinned by its own test below.
    """
    source, sequences = _source(universe=40, frontier=8)
    source.initial()
    source.advance()
    source.advance()

    frontier = source.frontier()
    assert frontier == sequences[: len(frontier)], "the refill reordered the window"
    assert len(frontier) > 8


def test_a_refill_doubles_rather_than_adding_a_fixed_batch():
    """So reaching a distant candidate costs log refills, not linear ones.

    On the Wolbachia design the frontier opens at 2,000 of 20,670 eligible
    candidates. Doubling reaches the whole universe in four refills; adding
    2,000 at a time would take nine.
    """
    source, _ = _source(universe=40, frontier=5)
    source.initial()

    sizes = [source.examined()]
    while source.advance():
        sizes.append(source.examined())

    assert sizes == [5, 10, 20, 40]


def test_advance_stops_at_the_universe_and_says_so():
    source, sequences = _source(universe=20, frontier=8)
    source.initial()

    assert source.advance() is True  # 16
    assert source.advance() is True  # 20, clamped
    assert source.frontier() == sequences
    assert source.advance() is False
    assert source.exhausted() is True


def test_an_unbounded_frontier_is_already_exhausted():
    """Nothing to refill when the first draw took everything."""
    source, sequences = _source(universe=20, frontier=None)
    source.initial()

    assert source.frontier() == sequences
    assert source.exhausted() is True
    assert source.advance() is False


def test_a_list_source_never_advances():
    """A `--candidates` file names the pool; there is nothing beyond it."""
    source = ListCandidateSource(["AAA", "CCC"])
    source.initial()

    assert source.advance() is False
    assert source.exhausted() is True


def test_advance_reports_the_examined_count_it_reached():
    source, _ = _source(universe=40, frontier=8)
    source.initial()
    source.advance()

    described = source.describe()
    assert described["examined"] == 16
    assert described["unexamined"] == 24
    assert described["exhausted"] is False


def test_the_kept_primers_stay_in_the_frontier():
    """A refill must not drop what the search has already chosen.

    `keep` exists so a caller can name the incumbent panel. Those primers are
    in the frontier already when they came from it, and this pins that an
    explicit keep cannot lose them either.
    """
    source, sequences = _source(universe=40, frontier=8)
    source.initial()

    source.advance(keep=[sequences[3], sequences[30]])

    assert sequences[3] in source.frontier()
    assert sequences[30] in source.frontier()


def test_the_inventory_is_queried_once_however_often_it_advances():
    """The ordered set is read once; advancing is a window over it."""
    from neoswga.core.candidate_provider import CandidateProvider

    inventory = _StubInventory([f"P{i:03d}" for i in range(40)])
    source = InventoryCandidateSource(
        CandidateProvider(inventory, "cond", [12]), frontier=8
    )
    source.initial()
    for _ in range(4):
        source.advance()

    assert inventory.queries == 1


# -- which exhaustion ------------------------------------------------------


def test_a_frontier_with_more_behind_it_is_not_an_exhausted_inventory():
    source, _ = _source(universe=40, frontier=8)
    source.initial()

    assert source.exhaustion() == FRONTIER_EXHAUSTED


def test_a_frontier_covering_the_universe_is_an_exhausted_inventory():
    source, _ = _source(universe=8, frontier=8)
    source.initial()

    assert source.exhaustion() == INVENTORY_EXHAUSTED


def test_a_list_source_is_always_an_exhausted_inventory():
    """Its universe and its frontier are the same set by definition."""
    source = ListCandidateSource(["AAA", "CCC"])
    source.initial()

    assert source.exhaustion() == INVENTORY_EXHAUSTED


@pytest.mark.parametrize("name", [FRONTIER_EXHAUSTED, INVENTORY_EXHAUSTED])
def test_the_two_reasons_are_distinguishable_strings(name):
    """They are reported to users and read by callers, so they must differ."""
    assert isinstance(name, str) and name
    assert FRONTIER_EXHAUSTED != INVENTORY_EXHAUSTED
