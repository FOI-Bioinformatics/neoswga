"""A beam keeps more than one partial panel alive.

Task 6 of the condition-aware pool design plan.

The greedy commits to its first choice. Where a constraint is not monotonic
under additions -- selectivity density is the one that matters here -- the
locally best first primer can be the one that makes the qualifying panel
unreachable, and the greedy has no way back.

A bounded beam keeps several partial panels and lets a worse-looking prefix
survive long enough to be completed. Bounded is the operative word: the cost is
the beam width times the pool size times the panel size in objective
evaluations, so every test here also states what it spends.

The fixtures are small enough to enumerate exhaustively, so "the beam found the
qualifying panel" is checked against a known answer rather than against the
beam's own opinion.
"""

import itertools
from types import SimpleNamespace

import pytest

from neoswga.core.panel_beam import beam_search
from neoswga.core.pool_objective import PoolConstraints, PoolObjective

A, B, C, D = "AAAACCCCGGGG", "CCCCGGGGTTTT", "GGGGTTTTAAAA", "TTTTAAAACCCC"


def _objective(table, constraints=None, default_density=100.0):
    """Panel -> metrics, from an explicit table keyed on the sorted panel."""

    def evaluate(primers):
        key = tuple(sorted(primers))
        row = table.get(key, {})
        coverage = row.get("coverage", 0.0)
        return SimpleNamespace(
            fg_coverage=coverage,
            effective_fg_coverage=coverage,
            selectivity_density=row.get("density", default_density),
            total_bg_sites=row.get("bg", 0),
        )

    return PoolObjective(evaluate, constraints or PoolConstraints())


def _every_panel(pool, size):
    """All panels up to `size`, so a fixture cannot be silently incomplete."""
    return [panel for n in range(1, size + 1) for panel in itertools.combinations(sorted(pool), n)]


def _oracle(pool, size, objective):
    """Every panel of exactly `size`, best feasible first."""
    feasible = [
        panel
        for panel in itertools.combinations(sorted(pool), size)
        if not objective.violations(list(panel))
    ]
    return sorted(feasible, key=lambda p: (-objective.coverage(list(p)), p))


def test_a_width_one_beam_is_the_greedy():
    """The degenerate case, so a difference later is the width and not the rule."""
    table = {
        (A,): {"coverage": 0.30},
        (B,): {"coverage": 0.10},
        (A, B): {"coverage": 0.35},
    }
    result = beam_search([A, B], _objective(table), size=2, beam_width=1)

    assert set(result.primers) == {A, B}
    assert result.status == "target_met"


def test_a_wider_beam_finds_the_panel_the_greedy_cannot_reach():
    """The case the beam exists for.

    A is the best single primer, and every panel containing it stays below the
    density floor. B and C are individually worse and jointly qualify. A greedy
    takes A first and never recovers.
    """
    table = {
        (A,): {"coverage": 0.40, "density": 1.0},
        (B,): {"coverage": 0.20, "density": 1.0},
        (C,): {"coverage": 0.15, "density": 1.0},
        (A, B): {"coverage": 0.50, "density": 1.0},
        (A, C): {"coverage": 0.45, "density": 1.0},
        (B, C): {"coverage": 0.30, "density": 9.0},
    }
    constraints = PoolConstraints(min_selectivity_density=5.0)
    objective = _objective(table, constraints)

    narrow = beam_search([A, B, C], objective, size=2, beam_width=1)
    wide = beam_search([A, B, C], objective, size=2, beam_width=3)

    assert narrow.violations, "width 1 should be stuck on the locally best prefix"
    assert set(wide.primers) == {B, C}
    assert wide.violations == ()
    assert set(wide.primers) == set(_oracle([A, B, C], 2, objective)[0])


def test_a_density_violation_does_not_prune_the_partial_panel():
    """Density is a ratio, so an addition can raise it back over the floor."""
    table = {
        (B,): {"coverage": 0.20, "density": 1.0},
        (C,): {"coverage": 0.15, "density": 1.0},
        (B, C): {"coverage": 0.30, "density": 9.0},
    }
    objective = _objective(table, PoolConstraints(min_selectivity_density=5.0))

    result = beam_search([B, C], objective, size=2, beam_width=2)

    assert result.violations == ()
    assert result.pruned == 0, "a recoverable violation was pruned"


def test_a_background_cap_violation_is_pruned():
    """Background sites only accumulate, so that panel cannot come back."""
    table = {
        (A,): {"coverage": 0.40, "bg": 500},
        (B,): {"coverage": 0.20, "bg": 1},
        (C,): {"coverage": 0.15, "bg": 1},
        (B, C): {"coverage": 0.30, "bg": 2},
    }
    objective = _objective(table, PoolConstraints(max_background_sites=100))

    result = beam_search([A, B, C], objective, size=2, beam_width=3)

    assert set(result.primers) == {B, C}
    assert result.pruned >= 1, "the over-cap prefix was carried anyway"


def test_pruning_costs_fewer_evaluations_than_not_pruning():
    """The point of pruning, stated as a number rather than assumed.

    Pruning only saves work when it leaves the beam holding fewer panels than
    its width. Three of the four candidates are over the cap here, so the beam
    carries one prefix instead of four and the next round is a quarter of the
    size. At a width the surviving prefixes still fill, a beam costs the same
    whether or not anything was pruned, which is what makes its cost bounded.
    """
    loads = {A: 500, B: 500, C: 500, D: 1}
    table = {
        panel: {"coverage": 0.1 * len(panel), "bg": sum(loads[p] for p in panel)}
        for panel in _every_panel([A, B, C, D], 2)
    }
    capped = _objective(table, PoolConstraints(max_background_sites=100))
    uncapped = _objective(table, PoolConstraints())

    with_prune = beam_search([A, B, C, D], capped, size=2, beam_width=4)
    without = beam_search([A, B, C, D], uncapped, size=2, beam_width=4)

    # Three over-cap singletons in the first round, then all three pairs that
    # extend the one survivor, because the load only accumulates.
    assert with_prune.pruned == 6
    assert with_prune.evaluations < without.evaluations


def test_a_dimerising_pair_never_reaches_a_returned_panel():
    """Excluded throughout, not scored down and outbid."""
    table = {
        (A,): {"coverage": 0.40},
        (B,): {"coverage": 0.39},
        (C,): {"coverage": 0.05},
        (A, B): {"coverage": 0.90},
        (A, C): {"coverage": 0.42},
        (B, C): {"coverage": 0.41},
    }

    def dimerises(candidate, selected):
        return {candidate, *selected} >= {A, B}

    result = beam_search([A, B, C], _objective(table), size=2, beam_width=3, dimerises=dimerises)

    assert {A, B} - set(result.primers), "the dimerising pair was delivered"
    assert len(result.primers) == 2


def test_the_evaluation_budget_is_reported_rather_than_exceeded():
    table = {(p,): {"coverage": 0.1} for p in (A, B, C, D)}
    result = beam_search([A, B, C, D], _objective(table), size=3, beam_width=4, max_evaluations=3)

    assert result.evaluations <= 3
    assert result.status == "budget_exhausted"


def test_a_pool_smaller_than_the_request_says_so():
    table = {(A,): {"coverage": 0.4}, (B,): {"coverage": 0.2}, (A, B): {"coverage": 0.5}}
    result = beam_search([A, B], _objective(table), size=5, beam_width=2)

    assert result.status == "inventory_exhausted"
    assert len(result.primers) == 2


def test_an_empty_pool_is_not_an_empty_panel():
    result = beam_search([], _objective({}), size=2, beam_width=2)

    assert result.status == "no_candidates"
    assert result.primers == ()


def test_a_smaller_feasible_panel_beats_a_larger_infeasible_one():
    """A constraint is not a scoring term the final size can outbid."""
    table = {
        (A,): {"coverage": 0.40, "bg": 1},
        (B,): {"coverage": 0.20, "bg": 1},
        (A, B): {"coverage": 0.60, "bg": 500},
    }
    objective = _objective(table, PoolConstraints(max_background_sites=100))

    result = beam_search([A, B], objective, size=2, beam_width=2)

    assert result.primers == (A,)
    assert result.violations == ()
    assert result.status == "inventory_exhausted", (
        "no qualifying panel of the requested size was reached, and saying "
        "otherwise would read as a proof that none exists"
    )


def test_the_result_does_not_depend_on_the_candidate_order():
    """Ties are broken on the sequence, so a reordered pool gives one answer."""
    table = {
        (A,): {"coverage": 0.20},
        (B,): {"coverage": 0.20},
        (C,): {"coverage": 0.20},
        (A, B): {"coverage": 0.40},
        (A, C): {"coverage": 0.40},
        (B, C): {"coverage": 0.40},
    }
    orders = [[A, B, C], [C, B, A], [B, A, C]]
    panels = {
        tuple(beam_search(o, _objective(table), size=2, beam_width=2).primers) for o in orders
    }

    assert len(panels) == 1, f"candidate order changed the panel: {panels}"


def test_the_beam_agrees_with_an_exhaustive_oracle_on_a_small_pool():
    """Four candidates, all fifteen panels enumerable."""
    coverages = {
        (A,): 0.10,
        (B,): 0.22,
        (C,): 0.18,
        (D,): 0.05,
        (A, B): 0.30,
        (A, C): 0.24,
        (A, D): 0.14,
        (B, C): 0.45,
        (B, D): 0.26,
        (C, D): 0.21,
    }
    table = {k: {"coverage": v} for k, v in coverages.items()}
    objective = _objective(table)

    result = beam_search([A, B, C, D], objective, size=2, beam_width=4)

    assert set(result.primers) == set(_oracle([A, B, C, D], 2, objective)[0])


def test_a_negative_budget_is_refused():
    with pytest.raises(ValueError):
        beam_search([A], _objective({}), size=1, max_evaluations=-1)


def test_a_beam_width_below_one_is_refused():
    with pytest.raises(ValueError):
        beam_search([A], _objective({}), size=1, beam_width=0)
