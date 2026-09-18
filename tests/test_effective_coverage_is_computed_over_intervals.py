"""Occupancy-weighted coverage, accumulated over window edges not over bases.

The rewrite the measurements in
`docs/validation/parallelism_opportunities_2026-09-17.md` and
`docs/validation/objective_evaluation_cost_2026-09-17.md` argue for.

`_compute_effective_coverage` was the dominant term in one objective
evaluation: 95% of it, against 16% for the 144 Mb background everyone assumed
was responsible. It made two passes over the entire target per primer, one to
reset a boolean window and one masked multiply, so its cost was linear in the
target length and linear in the panel size and independent of how many sites
there actually were.

The same quantity is a sum over segments between window endpoints. This pins
the equality, because a faster wrong answer is worse than a slow right one, and
because the failure mode of interval arithmetic is an off-by-one at a boundary
that produces a plausible number.

The oracle here is deliberately NOT the old loop copied across. It computes the
formula the docstring states, over an explicit float64 mask, using the shared
`_mark_window` for geometry:

    P(covered at x) = 1 - PRODUCT over primers p reaching x of (1 - theta_p)

with the union taken per primer first, since a primer does not stack with
itself.
"""

import math
import random
from types import SimpleNamespace

import numpy as np
import pytest

from neoswga.core.base_optimizer import BaseOptimizer
from neoswga.core.coverage import _mark_window
from neoswga.core.occupancy import site_occupancy
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import calculate_enthalpy_entropy

# Real 12-mers, self-dimer-free, spanning a range of occupancies at 30 C.
PRIMERS = [
    "GCTAAAGACAAT",
    "TACATAACATAC",
    "ACGTCAGCACGA",
    "CAGTCAGGATCA",
    "TTGACAGTCAAG",
]

_compute = BaseOptimizer._compute_effective_coverage


def _optimizer(reach, circular, conditions=None):
    """The smallest thing `_compute_effective_coverage` needs.

    It reads `self.conditions` and two fields of `self.config` and nothing
    else, so a real optimizer with an HDF5 index behind it would make these
    tests slower without making them stricter.
    """
    return SimpleNamespace(
        conditions=ReactionConditions(temp=30.0) if conditions is None else conditions,
        config=SimpleNamespace(extension_reach=reach, fg_circular=circular),
    )


def _oracle(optimizer, positions_by_primer, total_length):
    """The formula in the docstring, over a float64 mask. Slow and direct."""
    if optimizer.conditions is None:
        return None
    if not positions_by_primer or total_length <= 0:
        return 0.0
    reach = optimizer.config.extension_reach
    circular = optimizer.config.fg_circular
    not_covered = np.ones(total_length, dtype=np.float64)
    for primer, positions in positions_by_primer.items():
        if not positions:
            continue
        tm = optimizer.conditions.calculate_effective_tm(primer)
        dh, _ = calculate_enthalpy_entropy(primer)
        theta = site_occupancy(dh, tm, optimizer.conditions.temp)
        if theta <= 0.0:
            continue
        window = np.zeros(total_length, dtype=bool)
        for pos in positions:
            _mark_window(window, int(pos), reach, total_length, circular)
        not_covered[window] *= 1.0 - theta
    return float((1.0 - not_covered).sum()) / total_length


def _agree(optimizer, pbp, length, tol=1e-12):
    got = _compute(optimizer, pbp, length)
    want = _oracle(optimizer, pbp, length)
    assert got == pytest.approx(want, abs=tol), f"got {got!r}, oracle {want!r}"
    return got


# -- the degenerate answers ------------------------------------------------


def test_no_reaction_conditions_gives_no_number():
    """Occupancy needs a temperature; inventing one would look like a result."""
    optimizer = _optimizer(300, False, conditions=False)
    optimizer.conditions = None

    assert _compute(optimizer, {PRIMERS[0]: [10]}, 1000) is None


@pytest.mark.parametrize(
    "pbp,length",
    [
        ({}, 1000),
        ({PRIMERS[0]: []}, 1000),
        ({PRIMERS[0]: [10]}, 0),
        ({PRIMERS[0]: [10]}, -5),
    ],
)
def test_nothing_to_cover_is_zero(pbp, length):
    assert _compute(_optimizer(300, False), pbp, length) == 0.0


# -- linear geometry ------------------------------------------------------


def test_one_site_in_the_middle_of_a_linear_target():
    _agree(_optimizer(300, False), {PRIMERS[0]: [5000]}, 20000)


@pytest.mark.parametrize("pos", [0, 1, 299, 300, 301, 19699, 19700, 19999])
def test_a_site_near_either_end_of_a_linear_target_is_clipped(pos):
    """Where an off-by-one in the interval arithmetic would show up."""
    _agree(_optimizer(300, False), {PRIMERS[0]: [pos]}, 20000)


def test_one_primer_does_not_stack_with_itself():
    """Two overlapping windows of one primer apply its occupancy once.

    This is the grouping the docstring specifies and the property a naive
    per-site accumulation gets wrong, by multiplying (1 - theta) in for every
    overlapping window and understating coverage where a primer binds densely.
    """
    optimizer = _optimizer(300, False)
    close = _agree(optimizer, {PRIMERS[0]: [5000, 5010]}, 20000)

    theta_once = close * 20000 / 610.0
    assert theta_once == pytest.approx(
        _agree(optimizer, {PRIMERS[0]: [5000]}, 20000) * 20000 / 600.0, rel=1e-9
    ), "overlapping windows of one primer changed its per-base occupancy"


def test_two_primers_do_stack():
    """Distinct primers compound, which is what makes the product a product."""
    optimizer = _optimizer(300, False)
    one = _agree(optimizer, {PRIMERS[0]: [5000]}, 20000)
    two = _agree(optimizer, {PRIMERS[0]: [5000], PRIMERS[1]: [5000]}, 20000)

    assert two > one


# -- circular geometry ----------------------------------------------------


@pytest.mark.parametrize("pos", [0, 1, 10, 299, 19700, 19990, 19999])
def test_a_site_near_either_end_of_a_circular_target_wraps(pos):
    _agree(_optimizer(300, True), {PRIMERS[0]: [pos]}, 20000)


def test_wrapping_covers_more_than_clipping():
    """Guard the guard: a fixture where wrap is redundant proves nothing.

    A first attempt at this used sites within one reach of BOTH ends, where
    each site's window already covers the other end, so circular and linear
    agreed exactly and the test passed without exercising the wrap.
    """
    pbp = {PRIMERS[0]: [10]}
    wrapped = _agree(_optimizer(300, True), pbp, 20000)
    clipped = _agree(_optimizer(300, False), pbp, 20000)

    assert wrapped > clipped * 1.5


def test_a_window_spanning_the_whole_circle_covers_all_of_it():
    """2 * reach >= length on a circle. `_mark_window` marks everything."""
    optimizer = _optimizer(600, True)
    covered = _agree(optimizer, {PRIMERS[0]: [500]}, 1000)

    tm = optimizer.conditions.calculate_effective_tm(PRIMERS[0])
    dh, _ = calculate_enthalpy_entropy(PRIMERS[0])
    assert covered == pytest.approx(site_occupancy(dh, tm, 30.0), rel=1e-9)


def test_sites_at_both_ends_of_a_circle_do_not_double_count():
    """Their wrapped windows overlap; the primer still applies once."""
    _agree(_optimizer(300, True), {PRIMERS[0]: [5, 19995]}, 20000)


# -- the shape of a real panel --------------------------------------------


def test_a_panel_with_many_sites_each():
    optimizer = _optimizer(3000, True)
    rng = random.Random(20260917)
    pbp = {primer: sorted(rng.randrange(0, 1_267_782) for _ in range(20)) for primer in PRIMERS}

    _agree(optimizer, pbp, 1_267_782, tol=1e-9)


def test_dense_sites_whose_windows_all_overlap():
    """Merging has to collapse them, and the answer must not move."""
    optimizer = _optimizer(3000, False)
    pbp = {PRIMERS[0]: list(range(10_000, 20_000, 100))}

    _agree(optimizer, pbp, 100_000, tol=1e-9)


@pytest.mark.parametrize("seed", range(12))
def test_randomised_layouts_agree(seed):
    """Positions, panel size, reach and geometry all varied together."""
    rng = random.Random(seed)
    length = rng.choice([500, 2_000, 50_000])
    reach = rng.choice([1, 7, 250, 3_000])
    circular = rng.random() < 0.5
    panel = PRIMERS[: rng.randint(1, len(PRIMERS))]
    pbp = {primer: [rng.randrange(0, length) for _ in range(rng.randint(0, 6))] for primer in panel}

    _agree(_optimizer(reach, circular), pbp, length, tol=1e-10)


# -- occupancy edge cases -------------------------------------------------


class _FixedTm:
    """Conditions that report one melting temperature, whatever the primer.

    `ReactionConditions` refuses a temperature outside the polymerase's band,
    which is correct and means the two occupancy extremes cannot be reached by
    moving the temperature. They are reachable by moving the Tm: within the
    valid phi29 band the most extreme real 12-mer gives 0.999999999999 and
    5.5e-05, neither of which is exactly 1 or exactly 0.
    """

    def __init__(self, tm, temp=30.0):
        self.temp = temp
        self._tm = tm

    def calculate_effective_tm(self, primer):
        return self._tm


def test_a_primer_that_is_never_bound_contributes_nothing():
    """theta <= 0 is skipped, so it must not open an interval either."""
    optimizer = _optimizer(300, False)
    optimizer.conditions = _FixedTm(-1e6)

    pbp = {PRIMERS[0]: [5000]}
    assert _compute(optimizer, pbp, 20000) == 0.0
    assert _oracle(optimizer, pbp, 20000) == 0.0


def test_a_fully_occupied_primer_covers_its_window_completely():
    """theta == 1.0 makes the log-space weight infinite; it must still work.

    Accumulating log(1 - theta) is what makes the sweep possible, and at
    theta == 1 that is negative infinity. Adding it at one endpoint and
    subtracting it at another gives NaN, which is the worst kind of wrong
    answer, so saturation is counted rather than summed.

    No valid reaction reaches this: it needs a Tm near 200 C, and a 12-mer
    tops out around 70. It is handled because it is representable, not because
    it is expected.
    """
    optimizer = _optimizer(300, False)
    optimizer.conditions = _FixedTm(1e6)

    pbp = {PRIMERS[0]: [5000], PRIMERS[1]: [9000]}
    covered = _compute(optimizer, pbp, 20000)

    assert covered is not None and not math.isnan(covered)
    assert covered == pytest.approx(_oracle(optimizer, pbp, 20000), abs=1e-12)
    # Two disjoint 600 bp windows, each occupied all of the time.
    assert covered == pytest.approx(1200 / 20000, rel=1e-12)


def test_a_saturated_primer_mixed_with_an_ordinary_one():
    """The two accumulators have to combine, not shadow each other."""
    optimizer = _optimizer(300, False)

    class _PerPrimer:
        temp = 30.0

        def calculate_effective_tm(self, primer):
            return 1e6 if primer == PRIMERS[0] else 40.0

    optimizer.conditions = _PerPrimer()
    # Overlapping windows: one primer saturates the overlap, the other adds
    # its own occupancy outside it.
    _agree(optimizer, {PRIMERS[0]: [5000], PRIMERS[1]: [5100]}, 20000)


# -- the reason for the rewrite -------------------------------------------


def test_the_cost_does_not_grow_with_the_target_length():
    """The property the rewrite exists for, asserted as work not as wall clock.

    A timing assertion here would flake under load, and two files in this
    repository already do that. What is actually claimed is that the number of
    arithmetic steps depends on the sites rather than on the genome, so this
    counts the segments the sweep evaluates by watching `math.exp`, which it
    calls once per segment.
    """
    optimizer = _optimizer(300, False)
    pbp = {PRIMERS[0]: [5000], PRIMERS[1]: [9000]}

    calls = []
    real_exp = math.exp

    import neoswga.core.base_optimizer as module

    monkey = getattr(module, "math", math)
    original = monkey.exp
    monkey.exp = lambda x, _c=calls, _f=real_exp: (_c.append(x), _f(x))[1]
    try:
        small = _compute(optimizer, pbp, 20_000)
        at_small = len(calls)
        calls.clear()
        large = _compute(optimizer, pbp, 2_000_000)
        at_large = len(calls)
    finally:
        monkey.exp = original

    assert at_small == at_large, (
        f"segment count moved with the target length ({at_small} then {at_large}); "
        "the sweep is still doing work proportional to the genome"
    )
    assert at_small > 0, "no segments were evaluated; the probe is not measuring anything"
    assert large < small, "a longer genome with the same sites must have lower coverage"


# -- the geometry contract ------------------------------------------------


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("length,reach", [(50, 1), (50, 7), (50, 25), (50, 40), (1000, 300)])
def test_the_intervals_cover_exactly_what_mark_window_marks(circular, length, reach):
    """Two descriptions of one window, compared base by base.

    `merged_window_intervals` exists so coverage can be accumulated over edges
    instead of over bases, which is only sound if it describes the same window.
    A disagreement here is an off-by-one that would show up downstream as a
    coverage figure that is wrong without looking wrong.
    """
    from neoswga.core.coverage import merged_window_intervals

    rng = random.Random((length, reach, circular).__hash__() & 0xFFFF)
    for _ in range(40):
        positions = [rng.randrange(0, length) for _ in range(rng.randint(1, 6))]

        marked = np.zeros(length, dtype=bool)
        for pos in positions:
            _mark_window(marked, pos, reach, length, circular)

        swept = np.zeros(length, dtype=bool)
        for start, end in merged_window_intervals(positions, reach, length, circular):
            swept[start:end] = True

        assert np.array_equal(marked, swept), (
            f"positions {positions} at reach {reach} on length {length} "
            f"(circular={circular}): mask has {marked.sum()} bases, "
            f"intervals have {swept.sum()}"
        )


@pytest.mark.parametrize("circular", [False, True])
def test_the_intervals_are_disjoint_and_ordered(circular):
    """What the sweep relies on, and what merging is for."""
    from neoswga.core.coverage import merged_window_intervals

    spans = merged_window_intervals([10, 12, 500, 505, 990], 50, 1000, circular)

    assert spans, "no spans returned for sites that plainly have windows"
    for (_, end), (start, _) in zip(spans, spans[1:]):
        assert start > end, f"spans touch or overlap: {spans}"
    for start, end in spans:
        assert 0 <= start < end <= 1000


def test_no_positions_and_no_length_give_no_intervals():
    from neoswga.core.coverage import merged_window_intervals

    assert merged_window_intervals([], 300, 1000, False) == []
    assert merged_window_intervals([10], 300, 0, False) == []
    assert merged_window_intervals([10], -1, 1000, False) == []
