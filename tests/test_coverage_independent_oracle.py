"""Coverage checked against a base-by-base oracle written in this file.

Task 5 of the 2026-09-21 valid-design plan. Coverage is the headline number of
a design, and this codebase has produced disagreeing versions of it before: one
audit found three different semantics under one name.

The oracle here allocates a boolean per base and sets it. It does not call
`merged_window_intervals`, `_mark_window`, `_union_coverage` or anything they
call, so agreement is evidence rather than a restatement. It is deliberately
the slow implementation the production code replaced, which is what makes it
worth having: the fast one accumulates over window edges in log space, and an
edge-accounting error is invisible from inside that formulation.

Window convention, established by measurement rather than assumed: a site at
`pos` with reach `r` covers the half-open span `[pos - r, pos + r)`, so a
window is `2r` bases wide and `r = 20` covers 40 of 100.
"""

import math

import pytest

from neoswga.core.coverage import _mark_window, merged_window_intervals
from neoswga.core.occupancy_coverage import occupancy_weighted_coverage

LENGTH = 100


# ---------------------------------------------------------------------------
# The oracle
# ---------------------------------------------------------------------------


def covered_bases(sites, reach, length, circular, record_starts=None):
    """Which bases a primer's sites cover. One boolean per base, set directly.

    `record_starts` confines a window to the record holding its site, which is
    what `_mark_window` does when it is given them. Passing None is the
    unconfined convention, which is what the interval form implements.
    """
    marked = [False] * length
    for pos in sites:
        low, high = pos - reach, pos + reach
        if record_starts:
            index = _record_index(record_starts, pos)
            floor = record_starts[index]
            ceiling = record_starts[index + 1] if index + 1 < len(record_starts) else length
            low, high = max(low, floor), min(high, ceiling)
            for base in range(low, high):
                marked[base] = True
            continue
        for base in range(low, high):
            if 0 <= base < length:
                marked[base] = True
            elif circular:
                marked[base % length] = True
    return marked


def _record_index(record_starts, pos):
    """Which record `pos` falls in. A linear scan, so no bisect is borrowed."""
    index = 0
    for candidate, start in enumerate(record_starts):
        if start <= pos:
            index = candidate
    return index


def oracle_coverage(sites, reach, length, circular, record_starts=None):
    marked = covered_bases(sites, reach, length, circular, record_starts)
    return sum(marked) / length


def oracle_occupancy_coverage(positions_by_primer, thetas, reach, length, circular):
    """Probability a base is covered by at least one primer, base by base.

        P(covered at x) = 1 - PRODUCT over PRIMERS p reaching x of (1 - theta_p)

    The product is over primers, not sites: one primer's overlapping windows
    are unioned first and its occupancy applied once, because a primer does not
    stack with itself. That grouping is the model the production code states,
    and computing it here the slow way is how it gets checked.
    """
    total = 0.0
    reach_sets = {
        primer: covered_bases(sites, reach, length, circular)
        for primer, sites in positions_by_primer.items()
    }
    for base in range(length):
        uncovered = 1.0
        for primer, marked in reach_sets.items():
            if marked[base]:
                uncovered *= 1.0 - thetas[primer]
        total += 1.0 - uncovered
    return total / length


def intervals_to_bases(spans, length):
    marked = [False] * length
    for start, end in spans:
        for base in range(start, end):
            marked[base] = True
    return marked


# ---------------------------------------------------------------------------
# Geometry: the interval form against the oracle
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "sites,reach,circular",
    [
        ([50], 5, False),
        ([50], 5, True),
        ([0], 5, False),
        ([0], 5, True),
        ([99], 5, False),
        ([99], 5, True),
        ([2], 5, True),
        ([98], 5, True),
        ([50, 52], 5, False),
        ([50, 52, 54, 56], 3, False),
        ([10, 90], 20, True),
        ([], 5, False),
        ([50], 0, False),
        ([50], 500, False),
        ([50], 500, True),
        (list(range(0, 100, 7)), 4, True),
    ],
)
def test_the_interval_form_marks_the_bases_the_oracle_marks(sites, reach, circular):
    """Linear ends, circular wrap, no sites, overlaps and a reach past the genome."""
    spans = merged_window_intervals(sites, reach, LENGTH, circular)

    assert intervals_to_bases(spans, LENGTH) == covered_bases(sites, reach, LENGTH, circular)


def test_the_spans_are_disjoint_and_sorted():
    """Overlapping windows must merge, or a weighted accumulation double-counts."""
    spans = merged_window_intervals([10, 12, 14, 60], 5, LENGTH, False)

    assert spans == sorted(spans)
    for (_first_start, first_end), (second_start, _second_end) in zip(spans, spans[1:]):
        assert first_end < second_start, spans


def test_a_reach_covering_the_genome_on_a_circle_covers_everything():
    assert oracle_coverage([50], 500, LENGTH, True) == 1.0
    assert (
        intervals_to_bases(merged_window_intervals([50], 500, LENGTH, True), LENGTH)
        == [True] * LENGTH
    )


def test_no_sites_covers_nothing_however_long_the_reach():
    assert merged_window_intervals([], 10_000, LENGTH, True) == []
    assert oracle_coverage([], 10_000, LENGTH, True) == 0.0


# ---------------------------------------------------------------------------
# Record boundaries
# ---------------------------------------------------------------------------


RECORD_STARTS = [0, 40, 70]


@pytest.mark.parametrize("pos", [0, 1, 38, 39, 40, 41, 55, 68, 69, 70, 71, 98, 99])
def test_a_confined_window_stops_at_its_own_record_edge(pos):
    """Every boundary, and the base either side of each one.

    A window that crosses a join credits a panel with covering bases on a
    contig its site is not on, which on a fragmented assembly inflates coverage
    throughout.
    """
    import numpy as np

    marks = np.zeros(LENGTH, dtype=bool)
    _mark_window(marks, pos, 10, LENGTH, False, record_starts=RECORD_STARTS)

    assert list(marks) == covered_bases([pos], 10, LENGTH, False, RECORD_STARTS)


def test_the_interval_form_deliberately_does_not_confine():
    """Documented, and the reason the two production paths can disagree.

    `merged_window_intervals` takes no `record_starts`, with a stated reason:
    neither `_union_coverage` nor `_compute_effective_coverage` passes them to
    `_mark_window` today, so accepting them would let a caller believe a
    confinement that is not applied. The consequence is measured in the next
    test rather than left implicit.
    """
    import inspect

    assert "record_starts" not in inspect.signature(merged_window_intervals).parameters


def test_the_two_production_paths_disagree_across_a_record_join():
    """Named and measured, not fixed here.

    `compute_per_prefix_coverage` marks windows through `_mark_window` WITH
    record starts, so a window stops at a contig edge. `_union_coverage` and
    the occupancy path go through `merged_window_intervals`, which cannot
    confine. On a multi-record reference the same panel therefore has two
    coverage figures, and which one a reader sees depends on which code path
    produced it.

    This test exists so the size of that gap is on record. Changing it is
    Phase 6's subject; asserting it here stops it growing unnoticed.
    """
    import numpy as np

    site = 38  # 2 bases before the join at 40
    reach = 10

    confined = np.zeros(LENGTH, dtype=bool)
    _mark_window(confined, site, reach, LENGTH, False, record_starts=RECORD_STARTS)

    unconfined = intervals_to_bases(merged_window_intervals([site], reach, LENGTH, False), LENGTH)

    assert sum(confined) == 12, "28..40, stopping at the join"
    assert sum(unconfined) == 20, "28..48, running 8 bases into the next record"
    assert sum(unconfined) > sum(confined)


# ---------------------------------------------------------------------------
# Occupancy weighting
# ---------------------------------------------------------------------------


class Conditions:
    """Every primer melts exactly at the reaction temperature, so theta = 0.5.

    Supplied rather than derived, so the expected value below is arithmetic a
    reader can check rather than the output of the thing under test.
    """

    temp = 30.0

    def calculate_effective_tm(self, sequence):
        return self.temp


def test_no_sites_is_zero_even_with_a_long_reach():
    assert (
        occupancy_weighted_coverage(
            {}, 100, extension_reach=1000, circular=True, conditions=object()
        )
        == 0.0
    )


def test_no_conditions_is_unavailable_rather_than_zero():
    """None and 0.0 mean opposite things and must not share a representation.

    Without a temperature there is no occupancy to evaluate, and a fabricated
    zero would read exactly like a panel that binds nothing.
    """
    assert (
        occupancy_weighted_coverage(
            {"ACGTACGTACGT": [50]}, 100, extension_reach=10, circular=False, conditions=None
        )
        is None
    )


def assert_two_independent_primers_cover_half_a_window(actual):
    # Two distinct primers, each theta=0.5, both reaching 40 of 100 bases.
    assert actual == pytest.approx(0.4 * (1 - (1 - 0.5) ** 2))


def test_two_distinct_primers_on_the_same_window_compose_independently():
    """The grouping the docstring states: a product over primers."""
    positions = {"AAAACCCCGGGG": [50], "TTTTGGGGCCCC": [50]}

    actual = occupancy_weighted_coverage(
        positions, LENGTH, extension_reach=20, circular=False, conditions=Conditions()
    )

    assert_two_independent_primers_cover_half_a_window(actual)


def test_one_primers_overlapping_sites_do_not_stack_with_themselves():
    """The correction recorded as audit F5.

    Per-site independence would multiply (1 - theta) in once per overlapping
    window of the SAME primer and report a larger number. One primer's windows
    are unioned first and its occupancy applied once.
    """
    one_site = occupancy_weighted_coverage(
        {"AAAACCCCGGGG": [50]}, LENGTH, extension_reach=20, circular=False, conditions=Conditions()
    )
    overlapping = occupancy_weighted_coverage(
        {"AAAACCCCGGGG": [50, 52]},
        LENGTH,
        extension_reach=20,
        circular=False,
        conditions=Conditions(),
    )

    # 40 bases at theta 0.5 gives 0.20; adding an overlapping site widens the
    # union to 42 bases and must scale linearly with it, not compound.
    assert one_site == pytest.approx(0.40 * 0.5)
    assert overlapping == pytest.approx(0.42 * 0.5)


@pytest.mark.parametrize(
    "positions,reach,circular",
    [
        ({"AAAACCCCGGGG": [50]}, 20, False),
        ({"AAAACCCCGGGG": [50], "TTTTGGGGCCCC": [50]}, 20, False),
        ({"AAAACCCCGGGG": [10], "TTTTGGGGCCCC": [90]}, 20, True),
        ({"AAAACCCCGGGG": [0], "TTTTGGGGCCCC": [99]}, 5, True),
        ({"AAAACCCCGGGG": [10, 12, 14], "TTTTGGGGCCCC": [11, 13]}, 6, False),
        ({"AAAACCCCGGGG": [], "TTTTGGGGCCCC": [50]}, 10, False),
    ],
)
def test_the_edge_accumulation_agrees_with_the_base_by_base_oracle(positions, reach, circular):
    """The check the fast formulation most needs.

    Production accumulates log(1 - theta) at window edges and walks them. An
    error in that bookkeeping produces a plausible number, and nothing inside
    that formulation can see it. The oracle sets one boolean per base.
    """
    thetas = {primer: 0.5 for primer in positions}

    actual = occupancy_weighted_coverage(
        positions, LENGTH, extension_reach=reach, circular=circular, conditions=Conditions()
    )
    expected = oracle_occupancy_coverage(positions, thetas, reach, LENGTH, circular)

    assert actual == pytest.approx(expected)


def test_occupancy_coverage_never_exceeds_geometric_coverage():
    """A probability-weighted union cannot exceed the union it weights.

    Both are unions of the same windows; one multiplies each by a number in
    [0, 1]. A result above the geometric figure would mean the weighting had
    added coverage, which is not a thing it can do.
    """
    positions = {"AAAACCCCGGGG": [10, 40], "TTTTGGGGCCCC": [41, 80]}
    weighted = occupancy_weighted_coverage(
        positions, LENGTH, extension_reach=8, circular=False, conditions=Conditions()
    )
    geometric = oracle_coverage(
        sorted(site for sites in positions.values() for site in sites), 8, LENGTH, False
    )

    assert 0.0 <= weighted <= geometric


def test_a_saturated_primer_is_counted_rather_than_summed():
    """theta == 1 makes log(1 - theta) negative infinity.

    Adding it at one edge and subtracting it at another gives NaN, so a
    fully-bound primer is counted separately. No real reaction reaches this,
    and it is handled because the alternative failure is a silent NaN.
    """

    class Saturating(Conditions):
        def calculate_effective_tm(self, sequence):
            return 10_000.0

    value = occupancy_weighted_coverage(
        {"AAAACCCCGGGG": [50]},
        LENGTH,
        extension_reach=20,
        circular=False,
        conditions=Saturating(),
    )

    assert not math.isnan(value)
    assert value == pytest.approx(0.40)
