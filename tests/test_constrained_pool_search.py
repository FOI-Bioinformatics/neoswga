"""The greedy should choose on the quantity the design is judged on.

Task 6 of the condition-aware pool design plan.

The greedy counted NEW BINS: how many coverage regions a candidate adds that
nothing selected already covers. The design is accepted on occupancy-weighted
coverage, which weights each site by how much of the time it is actually bound.
Those disagree, and they disagree in a specific direction: a primer with many
sites and a melting temperature far below the reaction temperature touches many
bins and contributes little amplification.

So the search could prefer a panel the acceptance metric scores lower, and
improving the search did not reliably improve the result. Task 5 gave both ends
one definition; this uses it in the inner loop.

The fixtures here are small enough to enumerate exhaustively, so "the heuristic
found the qualifying pool" is checked against a known answer rather than against
the heuristic's own opinion.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.pool_objective import PoolConstraints, PoolObjective

WIDE = "AAAAAAAAAAAA"  # many bins, poorly bound
NARROW = "CCCCCCCCCCCC"  # fewer bins, well bound
OTHER = "GGGGGGGGGGGG"


class _Graph:
    """Just the mapping the scan reads."""

    def __init__(self, mapping):
        self.primer_to_regions = mapping


class _Dimers:
    """The matrix interface `_would_dimerise` uses."""

    def __init__(self, conflicting):
        self.conflicting = set(conflicting)

    def dimerises(self, candidate, selected):  # noqa: ARG002
        return candidate in self.conflicting


def _scan(objective=None, selected=None, covered=None, order=None, mapping=None):
    optimizer = DominatingSetOptimizer.__new__(DominatingSetOptimizer)
    optimizer.max_dimer_bp = None
    return DominatingSetOptimizer._select_next_primer(
        optimizer,
        order or [WIDE, NARROW],
        selected or [],
        covered or set(),
        _Graph(mapping or {WIDE: {1, 2, 3, 4}, NARROW: {5, 6}}),
        None,
        objective=objective,
    )


def _objective(coverages, background=None):
    """Panel -> coverage, with an optional background load per panel."""
    background = background or {}

    def evaluate(primers):
        key = tuple(sorted(primers))
        return SimpleNamespace(
            fg_coverage=coverages[key],
            effective_fg_coverage=coverages[key],
            selectivity_density=100.0,
            total_bg_sites=background.get(key, 0),
        )

    return PoolObjective(evaluate, PoolConstraints())


def test_without_an_objective_the_scan_still_counts_bins():
    """Legacy callers must be untouched."""
    best, gain, _ = _scan()

    assert best == WIDE
    assert gain == 4


def test_with_an_objective_the_scan_prefers_the_better_bound_primer():
    """The case the bin count gets wrong.

    WIDE touches twice as many bins as NARROW and contributes half the
    occupancy-weighted coverage, because its sites are mostly not bound at the
    reaction temperature.
    """
    objective = _objective({(): 0.0, (WIDE,): 0.10, (NARROW,): 0.30, (NARROW, WIDE): 0.35})

    best, gain, _ = _scan(objective=objective)

    assert best == NARROW, "the scan chose on bins, not on the accepted metric"
    assert gain == pytest.approx(0.30)


def test_an_exhaustive_oracle_agrees_with_the_scan():
    """Two candidates, both enumerable; the better single pick is known."""
    coverages = {(): 0.0, (WIDE,): 0.10, (NARROW,): 0.30, (NARROW, WIDE): 0.35}
    objective = _objective(coverages)

    oracle = max([WIDE, NARROW], key=lambda p: coverages[(p,)])
    best, _, _ = _scan(objective=objective)

    assert best == oracle


def test_ties_break_on_lower_background_load():
    """Equal coverage is separated by how much host the panel picks up."""
    objective = _objective(
        {(): 0.0, (WIDE,): 0.20, (NARROW,): 0.20},
        background={(WIDE,): 50, (NARROW,): 5},
    )

    best, _, _ = _scan(objective=objective)

    assert best == NARROW


def test_a_tie_on_both_falls_back_to_a_stable_order():
    """Reproducibility: the same pool must give the same panel."""
    objective = _objective({(): 0.0, (WIDE,): 0.20, (NARROW,): 0.20})

    first, _, _ = _scan(objective=objective, order=[WIDE, NARROW])
    second, _, _ = _scan(objective=objective, order=[NARROW, WIDE])

    assert first == second == min(WIDE, NARROW)


def test_a_candidate_adding_nothing_is_not_chosen():
    objective = _objective({(): 0.0, (WIDE,): 0.0, (NARROW,): 0.0})

    best, gain, _ = _scan(objective=objective)

    assert best is None or gain == pytest.approx(0.0)


def test_the_dimer_guard_still_runs_before_scoring():
    """Order matters: a primer the guard rejects must not return via scoring."""
    # `selected` is [OTHER], so every panel the scan asks about includes it.
    objective = _objective(
        {
            (OTHER,): 0.05,
            (OTHER, WIDE): 0.90,
            (NARROW, OTHER): 0.10,
        }
    )
    optimizer = DominatingSetOptimizer.__new__(DominatingSetOptimizer)
    optimizer.max_dimer_bp = 3

    best, _, skipped = DominatingSetOptimizer._select_next_primer(
        optimizer,
        [WIDE, NARROW],
        [OTHER],
        set(),
        _Graph({WIDE: {1}, NARROW: {2}}),
        _Dimers({WIDE}),  # WIDE dimerises with what is selected
        objective=objective,
    )

    assert best == NARROW
    assert skipped is True


def test_the_stop_reason_names_why_the_search_ended():
    """Exhausting a heuristic is not a proof that no panel exists.

    The greedy stopped for several different reasons and reported all of them
    the same way: a panel shorter than requested. A caller could not tell "the
    target was met", "the budget ran out", "the eligible inventory ran out" and
    "there were no candidates" apart, and the remedy differs for each.
    """
    from neoswga.core.dominating_set_optimizer import STOP_REASONS

    assert STOP_REASONS == (
        "target_met",
        "budget_exhausted",
        "inventory_exhausted",
        "no_qc_candidates",
    )


def test_a_full_coverage_stop_is_reported_as_the_target_being_met(tmp_path):
    import numpy as np

    h5py = pytest.importorskip("h5py")
    from neoswga.core.position_cache import PositionCache

    prefix = str(tmp_path / "t")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w") as db:
        db.create_dataset(WIDE, data=np.array([0, 5000], dtype=np.int64))
    string = pytest.importorskip("neoswga.core.string_search")

    optimizer = DominatingSetOptimizer(
        PositionCache([prefix], [WIDE]),
        fg_prefixes=[prefix],
        fg_seq_lengths=[10_000],
        extension_reach=10_000,
    )
    result = optimizer.optimize_greedy([WIDE], max_primers=5, verbose=False)

    assert result["stop_reason"] == "target_met"
    assert result["examined_candidates"] == 1
    assert result["eligible_candidates"] == 1
    assert result["search_budget"] == 5


def test_an_empty_pool_is_reported_as_having_no_candidates(tmp_path):
    h5py = pytest.importorskip("h5py")
    from neoswga.core.position_cache import PositionCache

    prefix = str(tmp_path / "t")
    with h5py.File(f"{prefix}_12mer_positions.h5", "w"):
        pass

    optimizer = DominatingSetOptimizer(
        PositionCache([prefix], []),
        fg_prefixes=[prefix],
        fg_seq_lengths=[10_000],
        extension_reach=3_000,
    )
    result = optimizer.optimize_greedy([], max_primers=5, verbose=False)

    assert result["stop_reason"] == "no_qc_candidates"
    assert result["eligible_candidates"] == 0


def _swap_objective(coverages, violations=None, background=None):
    violations = violations or {}
    background = background or {}

    class _Obj:
        def coverage(self, primers):
            return coverages[tuple(sorted(primers))]

        def violations(self, primers):
            return tuple(violations.get(tuple(sorted(primers)), ()))

        def shortfall(self, primers):
            """The searches rank on this; see `PoolObjective.shortfall`.

            This stub's violations are named in a table and carry no magnitude,
            so the count is the only distance available. Zero exactly when
            nothing is violated, which is the property the ordering needs.
            """
            return float(len(self.violations(primers)))

        def metrics(self, primers):
            key = tuple(sorted(primers))
            return SimpleNamespace(total_bg_sites=background.get(key, 0))

    return _Obj()


def test_a_feasible_panel_is_never_swapped_for_an_infeasible_one():
    """The incumbent has no violations; the alternative covers more and fails a limit.

    Coverage alone would take the swap. A constraint is not a scoring term to be
    outbid: once a panel satisfies the limits, leaving them is not an
    improvement however much coverage it buys.
    """
    from neoswga.core.swap_refinement import refine_by_swaps

    class _NoDimers:
        def dimerises(self, candidate, retained):  # noqa: ARG002
            return False

    objective = _swap_objective(
        coverages={(NARROW,): 0.50, (WIDE,): 0.95},
        violations={(WIDE,): ("background sites above maximum",)},
    )

    result = refine_by_swaps(
        [NARROW],
        [NARROW, WIDE],
        {NARROW: {1}, WIDE: {1, 2}},
        {1: 100, 2: 100},
        _NoDimers(),
        objective=objective,
    )

    assert result.primers == (NARROW,)


def test_an_infeasible_panel_moves_towards_feasibility_first():
    """Fewer violations wins even when coverage falls.

    Before a panel is feasible the useful direction is out of violation, not up
    the coverage curve.
    """
    from neoswga.core.swap_refinement import refine_by_swaps

    class _NoDimers:
        def dimerises(self, candidate, retained):  # noqa: ARG002
            return False

    objective = _swap_objective(
        coverages={(WIDE,): 0.95, (NARROW,): 0.50},
        violations={(WIDE,): ("selectivity below minimum",)},
    )

    result = refine_by_swaps(
        [WIDE],
        [WIDE, NARROW],
        {WIDE: {1, 2}, NARROW: {1}},
        {1: 100, 2: 100},
        _NoDimers(),
        objective=objective,
    )

    assert result.primers == (NARROW,)


def test_without_an_objective_the_swap_rule_is_unchanged():
    """Legacy callers keep the raw lexicographic rule."""
    from neoswga.core.swap_refinement import refine_by_swaps

    class _NoDimers:
        def dimerises(self, candidate, retained):  # noqa: ARG002
            return False

    result = refine_by_swaps(
        [NARROW],
        [NARROW, WIDE],
        {NARROW: {1}, WIDE: {1, 2}},
        {1: 100, 2: 100},
        _NoDimers(),
    )

    assert result.primers == (WIDE,), "raw covered bases should still take the swap"
