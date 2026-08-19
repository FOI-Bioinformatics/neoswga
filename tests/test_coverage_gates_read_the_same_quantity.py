"""Three 0.95 gates, three different quantities.

"Coverage" has six meanings in this tree. The authoritative one is the
fraction of genome *bases* within reach of a selected binding site
(`base_optimizer._union_coverage`), and that is what `fg_coverage` reports and
what a reader compares against a target. Three decisions are made against a
0.95 threshold, and none of them reads it:

- `dominating_set_adapter` marks a run PARTIAL below 0.95 of
  `result["coverage"]`, which is `len(covered_regions) / len(graph.regions)`.
  `graph.regions` only gains a bin when some *candidate* reaches it, so the
  denominator is the part of the genome the pool can touch, not the genome. A
  pool that reaches a tenth of the target and covers all of it scores 1.000 and
  reports SUCCESS. That same number is also returned as `score`.

- `hybrid_optimizer._prune_background` holds its floor against a binned
  estimate its own docstring describes as "an APPROXIMATION used for progress
  reporting".

- `background_aware_optimizer` holds the same floor against a third definition
  (directional interval merge over concatenated length, reach hardcoded 3000),
  and its module header advertises ">95% target coverage" on the strength of it.

The failure is quiet in the way that matters: every one of these is a number
between 0 and 1 that rises as the design improves, so nothing looks wrong until
the figure is compared with the one in the report.
"""

import numpy as np
import pytest


class _Cache:
    """Position cache stand-in. Sites are clustered in the first tenth of the
    genome, so the candidate pool cannot reach the rest of it."""

    def __init__(self, positions_by_primer):
        self._by_primer = positions_by_primer

    def get_positions(self, prefix, primer, strand="both"):
        pos = self._by_primer.get(primer, [])
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.asarray(pos, dtype=np.int64)


GENOME = 1_000_000
REACH = 3000
# Four primers, all binding inside the first 100 kb.
CLUSTERED = {
    "AAAACCCCGGGG": [5_000, 25_000],
    "TTTTGGGGCCCC": [45_000, 65_000],
    "ACGTACGTACGT": [15_000, 35_000],
    "TGCATGCATGCA": [55_000, 85_000],
}


@pytest.fixture
def optimizer():
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    return DominatingSetOptimizer(
        cache=_Cache(CLUSTERED),
        fg_prefixes=["fg"],
        fg_seq_lengths=[GENOME],
        bin_size=REACH // 4,
        extension_reach=REACH,
    )


def _true_coverage(primers):
    """Union of extension windows over the whole genome -- the authoritative
    definition, computed here independently of the code under test."""
    occupied = np.zeros(GENOME, dtype=bool)
    for primer in primers:
        for pos in CLUSTERED[primer]:
            occupied[max(0, pos - REACH) : min(GENOME, pos + REACH)] = True
    return occupied.sum() / GENOME


def test_reported_coverage_is_a_fraction_of_the_genome(optimizer):
    """The regression. Every site sits in the first 100 kb of a 1 Mb genome, so
    no design over this pool can exceed about 10% -- but the reachable-bin
    denominator made it 1.000."""
    result = optimizer.optimize_greedy(candidates=list(CLUSTERED), max_primers=4, verbose=False)

    assert _true_coverage(result["primers"]) < 0.15, "fixture no longer isolates the fault"
    assert result["coverage"] == pytest.approx(_true_coverage(result["primers"]), abs=0.02)


def test_the_reachable_fraction_is_still_available_and_is_higher(optimizer):
    """It is a useful diagnostic -- "the pool cannot reach the rest" is exactly
    what a user needs told -- so it is reported alongside rather than dropped."""
    result = optimizer.optimize_greedy(candidates=list(CLUSTERED), max_primers=4, verbose=False)

    assert result["reachable_coverage"] > result["coverage"]
    assert result["reachable_coverage"] == pytest.approx(1.0, abs=0.05)


def test_a_pool_that_cannot_reach_the_genome_is_not_reported_as_success():
    """The consequence a user sees: SUCCESS on a design covering a tenth of the
    target."""
    from neoswga.core.base_optimizer import OptimizationStatus
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter, DominatingSetConfig

    adapter = DominatingSetAdapter(
        position_cache=_Cache(CLUSTERED),
        fg_prefixes=["fg"],
        fg_seq_lengths=[GENOME],
        config=DominatingSetConfig(
            target_set_size=4,
            extension_reach=REACH,
            bin_size=REACH // 4,
            verbose=False,
        ),
    )
    result = adapter.optimize(list(CLUSTERED))

    assert result.status is OptimizationStatus.PARTIAL
    assert result.score < 0.15


def test_the_stop_condition_targets_the_genome_not_the_pool(optimizer):
    """`min_coverage=0.95` means 95% of the target. A pool that cannot reach
    that must exhaust its budget, not stop early declaring victory."""
    result = optimizer.optimize_greedy(
        candidates=list(CLUSTERED), max_primers=4, min_coverage=0.95, verbose=False
    )

    assert len(result["primers"]) == 4


def test_the_stop_condition_still_fires_when_it_is_genuinely_met():
    """Guard against the opposite failure: an achievable target must still stop
    the greedy early rather than spending the whole budget."""
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    # Two primers tile a 60 kb genome at a 30 kb forward reach; the other two
    # are redundant, so a working stop condition never reaches them.
    spread = {
        "AAAACCCCGGGG": [0],
        "TTTTGGGGCCCC": [30_000],
        "ACGTACGTACGT": [1_000],
        "TGCATGCATGCA": [31_000],
    }
    optimizer = DominatingSetOptimizer(
        cache=_Cache(spread),
        fg_prefixes=["fg"],
        fg_seq_lengths=[60_000],
        bin_size=REACH // 4,
        extension_reach=30_000,
    )
    result = optimizer.optimize_greedy(
        candidates=list(spread), max_primers=4, min_coverage=0.90, verbose=False
    )

    assert len(result["primers"]) < 4
    assert result["coverage"] >= 0.90
