"""The fg/bg prefilter stops deleting candidates and orders the scan instead.

Phase 4 increment 6 of the plan for `docs/validation/pipeline_audit_2026-09-16/`,
which states the decision: "The background prefilter becomes an ordering
heuristic. It stops deleting candidates and instead orders the scan."

Two reasons, and the measurement is the first.

`--min-fg-bg-ratio` was inert. `_prefilter_by_background` kept every candidate
at or above the ratio, and then, if that removed more than
`max_removal_fraction` of them, threw the threshold away and kept the top
80% by ratio instead. On the real Wolbachia shortlist the threshold would remove
64.8% at its default of 1.0, so the clause fired and exactly 400 of 2,000 were
removed. It fired identically at 2.0 and 5.0. So the flag a user sets changed
nothing above about 1.0, and the rule in force was "always drop the worst 20%".

A bound on the fraction of a batch is not a rule about candidates, either. Which
candidates survive depended on how many others happened to be in the same batch.

Ordering has neither problem, and after increment 5 a candidate placed at the
back of the scan is still reachable: the frontier refills. Deleting one is
final, which is the shape of Known Issue 9, where an evenness gate removed
primers the optimizer had already selected.
"""

import pytest

from neoswga.core.unified_optimizer import order_candidates_by_background


class _Cache:
    """Positions per (prefix, primer), so fg/bg ratios are controllable."""

    def __init__(self, mapping):
        self._mapping = mapping

    def get_positions(self, prefix, primer, strand="both"):
        return list(self._mapping.get((prefix, primer), ()))


FG, BG = ["fg"], ["bg"]


def _cache(spec):
    """`spec` maps primer to (foreground sites, background sites)."""
    mapping = {}
    for primer, (fg_sites, bg_sites) in spec.items():
        mapping[("fg", primer)] = list(range(fg_sites))
        mapping[("bg", primer)] = list(range(bg_sites))
    return _Cache(mapping)


SPEC = {
    "AAAA": (10, 0),  # ratio 10.0, well above
    "CCCC": (10, 4),  # ratio 2.0, above
    "GGGG": (2, 9),  # ratio 0.2, below
    "TTTT": (1, 19),  # ratio 0.05, below
}
POOL = ["AAAA", "CCCC", "GGGG", "TTTT"]


# -- nothing is deleted ----------------------------------------------------


def test_every_candidate_survives():
    """The property the whole change is for."""
    ordered, rejected = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=1.0)

    assert sorted(ordered) == sorted(POOL)


def test_the_ones_below_the_threshold_go_to_the_back():
    ordered, rejected = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=1.0)

    assert ordered[:2] == ["AAAA", "CCCC"]
    assert sorted(ordered[2:]) == ["GGGG", "TTTT"]


def test_the_ones_below_the_threshold_are_reported_with_a_reason():
    """Deprioritised is not the same as deleted, and a reader needs to know."""
    _, rejected = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=1.0)

    assert set(rejected) == {"GGGG", "TTTT"}
    for primer, reason in rejected.items():
        assert "ratio" in reason, reason
        assert "1.0" in reason or "1" in reason


def test_the_order_within_each_group_is_the_order_it_arrived_in():
    """A stable partition, not a re-sort.

    Increment 5 established the inventory's `search_rank` as the traversal
    order, and the optimizers are order-sensitive. Re-sorting the whole pool by
    ratio would discard that ranking; partitioning keeps it inside each group.
    """
    pool = ["CCCC", "AAAA", "TTTT", "GGGG"]

    ordered, _ = order_candidates_by_background(_cache(SPEC), pool, FG, BG, min_ratio=1.0)

    assert ordered == ["CCCC", "AAAA", "TTTT", "GGGG"]


def test_a_threshold_nothing_meets_still_returns_everything():
    """And in the order it arrived, since the partition is then a no-op."""
    ordered, rejected = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=1000.0)

    assert ordered == POOL
    assert len(rejected) == len(POOL)


def test_a_threshold_everything_meets_rejects_nothing():
    ordered, rejected = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=0.0)

    assert ordered == POOL
    assert rejected == {}


def test_an_empty_pool_is_not_an_error():
    ordered, rejected = order_candidates_by_background(_cache(SPEC), [], FG, BG, min_ratio=1.0)

    assert ordered == []
    assert rejected == {}


def test_no_background_prefixes_means_no_ordering():
    """There is no ratio to order on, and inventing one would be worse."""
    ordered, rejected = order_candidates_by_background(_cache(SPEC), POOL, FG, [], min_ratio=1.0)

    assert ordered == POOL
    assert rejected == {}


# -- the fraction clause is gone -------------------------------------------


def test_the_outcome_does_not_depend_on_how_many_others_share_the_batch():
    """What the removed clause got wrong.

    Under `max_removal_fraction` the survivors depended on the size of the
    batch: the same candidate could be kept in one pool and dropped in another
    purely because of how many others were below the threshold. Ordering is a
    per-candidate decision.
    """
    small = ["AAAA", "GGGG"]
    large = POOL + ["GGGG2", "TTTT2"]
    spec = dict(SPEC, GGGG2=(2, 9), TTTT2=(1, 19))

    _, rejected_small = order_candidates_by_background(_cache(spec), small, FG, BG, min_ratio=1.0)
    _, rejected_large = order_candidates_by_background(_cache(spec), large, FG, BG, min_ratio=1.0)

    assert ("GGGG" in rejected_small) == ("GGGG" in rejected_large)


def test_the_threshold_actually_changes_the_partition():
    """The flag was inert above about 1.0 on a real pool; it must not be now."""
    at_one, _ = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=1.0)
    at_five, _ = order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=5.0)

    assert at_one[:2] == ["AAAA", "CCCC"]
    assert at_five[0] == "AAAA", "a stricter threshold left the partition unchanged"
    assert "CCCC" in at_five[1:], "CCCC has ratio 2.0 and must fall behind at 5.0"


@pytest.mark.parametrize("bad", [-1.0, float("nan"), float("inf")])
def test_a_threshold_that_is_not_a_usable_number_is_refused(bad):
    with pytest.raises(ValueError):
        order_candidates_by_background(_cache(SPEC), POOL, FG, BG, min_ratio=bad)


def test_the_old_deleting_prefilter_is_gone():
    """So a caller cannot quietly keep the behaviour the measurement retired."""
    import neoswga.core.unified_optimizer as module

    assert not hasattr(module, "_prefilter_by_background"), (
        "the deleting prefilter is back. It made --min-fg-bg-ratio inert above "
        "about 1.0 and let batch size decide which candidates survived."
    )
