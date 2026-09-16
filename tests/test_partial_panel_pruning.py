"""Which violating partial panels can be abandoned, and which cannot.

Task 6 of the condition-aware pool design plan.

A search over partial panels has to decide whether a panel that already violates
a constraint is worth extending. Getting that wrong in one direction wastes
work; in the other it discards panels that would have qualified.

The two constraints behave differently under additions, and the difference is
not a detail:

background site COUNT is monotonic
    Every primer added contributes its own sites, so the count never falls. A
    partial panel already above the cap cannot be rescued, and pruning it is
    safe.

selectivity DENSITY is not
    It is a ratio of target to background site density. Adding a primer that
    binds the target well and the background little RAISES it. A partial panel
    below the floor may be above it two primers later, so pruning on density is
    discarding solutions.

Dimer violations are excluded throughout rather than pruned: a panel containing
an incompatible pair is not a partial solution on its way anywhere, because no
addition removes a pair that is already in it.
"""

import pytest

from neoswga.core.partial_panel import can_prune

A, C = "AAAAAAAAAAAA", "CCCCCCCCCCCC"


def test_a_panel_over_the_background_cap_is_pruned():
    """Monotonic: additions only add sites."""
    assert (
        can_prune(
            violations=("background sites above maximum",),
            max_background_sites=100,
            min_selectivity_density=None,
        )
        is True
    )


def test_a_panel_below_the_density_floor_is_not_pruned():
    """Not monotonic: a later primer can raise the ratio."""
    assert (
        can_prune(
            violations=("selectivity below minimum",),
            max_background_sites=None,
            min_selectivity_density=10.0,
        )
        is False
    )


def test_a_panel_failing_both_is_pruned_on_the_monotonic_one():
    assert (
        can_prune(
            violations=("selectivity below minimum", "background sites above maximum"),
            max_background_sites=100,
            min_selectivity_density=10.0,
        )
        is True
    )


def test_a_dimer_violation_is_never_a_partial_solution():
    """No addition removes a pair already in the panel."""
    assert (
        can_prune(
            violations=("dimer constraint",),
            max_background_sites=None,
            min_selectivity_density=None,
        )
        is True
    )


def test_a_feasible_panel_is_not_pruned():
    assert can_prune(violations=(), max_background_sites=100, min_selectivity_density=10.0) is False


def test_an_unavailable_coverage_is_not_grounds_for_pruning():
    """It means not measured, which is a reason to look rather than to stop."""
    assert (
        can_prune(
            violations=("coverage unavailable",),
            max_background_sites=None,
            min_selectivity_density=None,
        )
        is False
    )


def test_the_cap_must_be_configured_for_the_count_rule_to_apply():
    """Without a cap there is no monotonic bound to prune against."""
    assert (
        can_prune(
            violations=("background sites above maximum",),
            max_background_sites=None,
            min_selectivity_density=None,
        )
        is False
    )


@pytest.mark.parametrize("unknown", ["something new", "a future constraint"])
def test_an_unrecognised_violation_is_never_pruned(unknown):
    """Conservative by default: a constraint nobody has reasoned about here
    might well be non-monotonic, and discarding solutions is the worse error."""
    assert (
        can_prune(violations=(unknown,), max_background_sites=100, min_selectivity_density=10.0)
        is False
    )


def test_the_greedy_stops_once_a_monotonic_limit_is_passed():
    """The rule has to be applied, or it is a module nothing calls.

    A greedy grows ONE panel by addition, which is exactly the situation the
    monotonic argument covers: once the panel is over the background cap, every
    further addition leaves it over. Continuing spends the whole size budget on
    panels that cannot qualify.

    This is deliberately NOT applied to `plan_pool`'s size sweep, where each
    size is optimised independently and a larger request is not a superset of a
    smaller one. The argument is about additions, not about requests.
    """
    from types import SimpleNamespace

    from neoswga.core.dominating_set_optimizer import _should_stop_extending
    from neoswga.core.pool_objective import PoolConstraints, PoolObjective

    def evaluate(primers):
        return SimpleNamespace(
            fg_coverage=0.5,
            effective_fg_coverage=0.5,
            selectivity_density=100.0,
            total_bg_sites=1000 if len(primers) > 1 else 5,
        )

    objective = PoolObjective(evaluate, PoolConstraints(max_background_sites=100))

    assert _should_stop_extending(objective, [A]) is False
    assert _should_stop_extending(objective, [A, C]) is True


def test_a_density_floor_never_stops_the_greedy():
    from types import SimpleNamespace

    from neoswga.core.dominating_set_optimizer import _should_stop_extending
    from neoswga.core.pool_objective import PoolConstraints, PoolObjective

    def evaluate(primers):
        return SimpleNamespace(
            fg_coverage=0.5,
            effective_fg_coverage=0.5,
            selectivity_density=1.0,
            total_bg_sites=0,
        )

    objective = PoolObjective(evaluate, PoolConstraints(min_selectivity_density=10.0))

    assert _should_stop_extending(objective, [A, C]) is False


def test_without_an_objective_nothing_stops_early():
    from neoswga.core.dominating_set_optimizer import _should_stop_extending

    assert _should_stop_extending(None, [A, C]) is False


def test_the_greedy_actually_stops_when_the_cap_is_passed(tmp_path):
    """End to end through `optimize_greedy`, not just the helper.

    A rule applied nowhere is the failure this project has met before, so this
    drives the real loop and checks it delivered fewer primers than the budget.
    """
    import numpy as np

    h5py = pytest.importorskip("h5py")
    from types import SimpleNamespace

    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
    from neoswga.core.pool_objective import PoolConstraints, PoolObjective
    from neoswga.core.position_cache import PositionCache

    primers = ["AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"]
    prefix = str(tmp_path / "t")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w") as db:
        for i, primer in enumerate(primers):
            db.create_dataset(primer, data=np.array([i * 9000 + 100], dtype=np.int64))

    optimizer = DominatingSetOptimizer(
        PositionCache([prefix], primers),
        fg_prefixes=[prefix],
        fg_seq_lengths=[40_000],
        extension_reach=3_000,
    )

    # Over the cap as soon as a second primer is added.
    def evaluate(panel):
        return SimpleNamespace(
            fg_coverage=0.1 * len(panel),
            effective_fg_coverage=0.1 * len(panel),
            selectivity_density=100.0,
            total_bg_sites=5 if len(panel) < 2 else 5000,
        )

    objective = PoolObjective(evaluate, PoolConstraints(max_background_sites=100))
    result = optimizer.optimize_greedy(primers, max_primers=4, verbose=False, objective=objective)

    assert len(result["primers"]) < 4, "the greedy spent its whole budget past the cap"
