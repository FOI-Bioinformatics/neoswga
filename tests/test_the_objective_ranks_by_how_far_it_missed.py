"""Among panels that fail, the one closer to passing is the better one.

Found by increment 5's measurement
(`docs/validation/frontier_refill_2026-09-17.md`). Widening the candidate
frontier moved a panel AWAY from the selectivity floor it was chasing: density
fell from 20.9 to 14.6 while coverage rose from 0.740 to 0.765, on a run whose
only failed constraint was that floor.

The cause is that `violations()` returns names and the searches ranked on
`len(...)` of them. Two panels failing the same single constraint therefore tie,
coverage breaks the tie, and the deciding metric is free to drift. More
candidates give the coverage term more to work with, so more search made it
worse.

`refine_by_swaps` already states the intended rule: "before a panel is feasible
the useful direction is out of violation rather than up the coverage curve."
Counting implements that only across different numbers of violated constraints,
never within one. `PoolObjective.shortfall` measures how far a panel is from
feasibility so the rule holds inside one constraint too.

The contract that matters most is the boundary: shortfall is zero exactly when
nothing is violated. If those two ever disagree, a feasible panel could be
ranked behind an infeasible one, which is the failure this ordering exists to
prevent.
"""

import math
from types import SimpleNamespace

import pytest

from neoswga.core.pool_objective import PoolConstraints, PoolObjective


def _objective(table, **limits):
    """A real objective over a table of panel metrics."""

    def evaluate(primers):
        coverage, density, sites = table[tuple(sorted(primers))]
        return SimpleNamespace(
            effective_fg_coverage=coverage,
            fg_coverage=coverage,
            selectivity_density=density,
            total_bg_sites=sites,
            max_gap=100,
        )

    return PoolObjective(evaluate, PoolConstraints(**limits))


A, B, C = ["AAA"], ["BBB"], ["CCC"]
KEY = lambda panel: tuple(sorted(panel))  # noqa: E731


# -- the boundary ----------------------------------------------------------


def test_a_feasible_panel_has_no_shortfall():
    objective = _objective({KEY(A): (0.9, 50.0, 10)}, min_selectivity_density=40.0)

    assert objective.shortfall(A) == 0.0
    assert objective.violations(A) == ()


@pytest.mark.parametrize(
    "density,sites,floor,ceiling",
    [
        (50.0, 10, 40.0, None),
        (40.0, 10, 40.0, None),  # exactly at the floor is not a violation
        (39.9, 10, 40.0, None),
        (50.0, 10, None, 10),  # exactly at the ceiling is not a violation
        (50.0, 11, None, 10),
        (10.0, 99, 40.0, 10),
        (50.0, 1, 40.0, 10),
    ],
)
def test_shortfall_is_zero_exactly_when_nothing_is_violated(density, sites, floor, ceiling):
    """The property the whole ordering rests on.

    If a feasible panel could carry a positive shortfall, or an infeasible one a
    zero, the search could rank a panel that passes behind one that does not.
    """
    limits = {}
    if floor is not None:
        limits["min_selectivity_density"] = floor
    if ceiling is not None:
        limits["max_background_sites"] = ceiling
    objective = _objective({KEY(A): (0.9, density, sites)}, **limits)

    assert (objective.shortfall(A) == 0.0) is (objective.violations(A) == ())


def test_no_constraints_means_no_shortfall():
    objective = _objective({KEY(A): (0.9, 0.001, 10_000)})

    assert objective.shortfall(A) == 0.0


# -- the magnitude ---------------------------------------------------------


def test_a_panel_closer_to_the_density_floor_has_less_shortfall():
    """The regression from the real pool, in miniature."""
    objective = _objective(
        {KEY(A): (0.740, 20.852, 151), KEY(B): (0.765, 14.589, 211)},
        min_selectivity_density=100.0,
    )

    assert objective.shortfall(A) < objective.shortfall(B)


def test_a_panel_closer_to_the_background_ceiling_has_less_shortfall():
    objective = _objective(
        {KEY(A): (0.9, 50.0, 12), KEY(B): (0.9, 50.0, 40)},
        max_background_sites=10,
    )

    assert objective.shortfall(A) < objective.shortfall(B)


def test_two_violated_constraints_weigh_more_than_one():
    objective = _objective(
        {KEY(A): (0.9, 20.0, 5), KEY(B): (0.9, 20.0, 500)},
        min_selectivity_density=40.0,
        max_background_sites=10,
    )

    assert objective.violations(A) == ("selectivity below minimum",)
    assert len(objective.violations(B)) == 2
    assert objective.shortfall(B) > objective.shortfall(A)


def test_the_shortfall_is_relative_so_the_two_limits_are_comparable():
    """Sites and densities are different units; the sum has to mean something.

    Missing a floor of 100 by half is the same distance as exceeding a ceiling
    of 10 by half, so neither constraint dominates merely by being measured on
    a larger scale.
    """
    floor_miss = _objective({KEY(A): (0.9, 50.0, 0)}, min_selectivity_density=100.0)
    ceiling_miss = _objective({KEY(A): (0.9, 1e9, 15)}, max_background_sites=10)

    assert floor_miss.shortfall(A) == pytest.approx(0.5)
    assert ceiling_miss.shortfall(A) == pytest.approx(0.5)


def test_an_unmeasurable_coverage_dominates_every_other_shortfall():
    """It is not a distance. A panel that cannot be scored is not nearly-feasible."""
    objective = _objective({KEY(A): (None, 1.0, 10_000)}, min_selectivity_density=100.0)

    assert objective.shortfall(A) == math.inf
    assert "coverage unavailable" in objective.violations(A)


# -- what the searches do with it -----------------------------------------


class _Dimers:
    def dimerises(self, incoming, retained):
        return False


def test_the_swap_loop_moves_towards_feasibility_not_towards_coverage():
    """The regression. Coverage used to decide between two failing panels."""
    from neoswga.core.swap_refinement import refine_by_swaps

    start, closer, further = ["P0"], ["P1"], ["P2"]
    objective = _objective(
        {
            KEY(start): (0.70, 10.0, 0),
            KEY(closer): (0.74, 20.852, 0),  # nearer the floor, less coverage
            KEY(further): (0.765, 14.589, 0),  # more coverage, further away
        },
        min_selectivity_density=100.0,
    )
    bins = {"P0": {0}, "P1": {1}, "P2": {2}}
    weights = {0: 10, 1: 10, 2: 10}

    result = refine_by_swaps(
        start, ["P0", "P1", "P2"], bins, weights, _Dimers(), objective=objective
    )

    assert (
        list(result.primers) == closer
    ), "the swap took the higher-coverage panel that sits further from the floor"


def test_a_feasible_panel_is_still_never_traded_for_an_infeasible_one():
    """The property the count ordering got right, which must survive."""
    from neoswga.core.swap_refinement import refine_by_swaps

    start, tempting = ["P0"], ["P1"]
    objective = _objective(
        {KEY(start): (0.60, 150.0, 0), KEY(tempting): (0.99, 10.0, 0)},
        min_selectivity_density=100.0,
    )

    result = refine_by_swaps(
        start, ["P0", "P1"], {"P0": {0}, "P1": {1}}, {0: 10, 1: 10}, _Dimers(), objective=objective
    )

    assert list(result.primers) == start
    assert objective.violations(list(result.primers)) == ()


def test_the_beam_ranks_the_closer_panel_first():
    """The beam had the same count ordering, and needs the same fix."""
    from neoswga.core.panel_beam import _rank

    closer, further = ["P1"], ["P2"]
    objective = _objective(
        {KEY(closer): (0.74, 20.852, 0), KEY(further): (0.765, 14.589, 0)},
        min_selectivity_density=100.0,
    )

    closer_key = _rank(objective, tuple(closer))
    further_key = _rank(objective, tuple(further))

    assert closer_key < further_key, "the beam preferred the panel further from the floor"


def test_the_beam_still_puts_any_feasible_panel_ahead_of_any_infeasible_one():
    from neoswga.core.panel_beam import _rank

    feasible, infeasible = ["P1"], ["P2"]
    objective = _objective(
        {KEY(feasible): (0.10, 150.0, 0), KEY(infeasible): (0.99, 99.0, 0)},
        min_selectivity_density=100.0,
    )

    feasible_key = _rank(objective, tuple(feasible))
    infeasible_key = _rank(objective, tuple(infeasible))

    assert feasible_key < infeasible_key


# -- what makes the ordering safe -----------------------------------------


def test_a_failed_repair_returns_the_panel_it_was_given():
    """Otherwise chasing feasibility can cost coverage for nothing.

    Ranking on distance from feasibility is right while feasibility is
    reachable. On a constraint no panel can meet, it will trade real coverage
    for a step toward a floor it never reaches, and the row is rejected either
    way. Keeping the original panel is what makes the ordering safe: progress
    is pursued, and failing to get there costs nothing.

    Two existing guards in `tests/test_pool_plan_repair.py` depend on this and
    pass unchanged because of it.
    """
    from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
    from neoswga.core.pool_planner import plan_pool

    start = ["AAAAAAAAAAAA", "CCCCCCCCCCCC"]
    nearer = ["AAAAAAAAAAAA", "ACACACACACAC"]
    table = {
        KEY(start): (0.90, 2.0, 40),  # best coverage, furthest from the floor
        KEY(nearer): (0.40, 3.0, 45),  # nearer the floor, much less coverage
        KEY(["CCCCCCCCCCCC", "ACACACACACAC"]): (0.30, 1.0, 50),
        KEY(["AAAAAAAAAAAA"]): (0.30, 2.0, 20),
        KEY(["CCCCCCCCCCCC"]): (0.20, 1.0, 25),
        KEY(["ACACACACACAC"]): (0.10, 1.5, 30),
    }

    class _Opt:
        name = "test"
        bg_prefixes = ["bg"]
        bg_seq_lengths = [1000]
        conditions = object()
        cache = None
        config = OptimizerConfig(objective_scan_width=None, max_frontier_refills=0)

        def optimize(self, candidates, target_size):
            return SimpleNamespace(primers=start, status=OptimizationStatus.PARTIAL, message="")

        def compute_metrics(self, primers):
            coverage, density, sites = table[KEY(primers)]
            return SimpleNamespace(
                effective_fg_coverage=coverage,
                fg_coverage=coverage,
                selectivity_density=density,
                total_bg_sites=sites,
                max_gap=100,
            )

    row = plan_pool(
        _Opt(),
        [*start, "ACACACACACAC"],
        [2],
        [0.8],
        primer_length=12,
        min_selectivity_density=10.0,
    )["rows"][0]

    assert row["repair"]["succeeded"] is False
    assert set(row["primers"]) == set(start), (
        "a failed repair delivered the panel the search wandered to, which here "
        "costs 0.50 of coverage for 1.0 of density against a floor of 10"
    )
