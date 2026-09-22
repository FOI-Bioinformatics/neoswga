"""The exact solver answers the fixed-budget coverage question.

`DominatingSetOptimizer.optimize_ilp` used to minimise primer *count* subject to
covering every reachable bin. That is the opposite of what every greedy
optimizer in the tree does, so it was infeasible below the minimum cover and
could not bound the greedy result at any budget a user actually asks for.

These tests pin the max-coverage formulation instead:

    maximise    sum_b w_b * y_b
    subject to  y_b <= sum_{i covers b} x_i
                sum_i x_i <= S

with `w_b` the bases the bin spans.
"""

from unittest.mock import Mock

import numpy as np
import pytest

from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

pytest.importorskip("mip", reason="python-mip ships in the 'improved' extra")


def _optimizer(bin_size=10000, extension_reach=0, genome_length=100000):
    """Four primers over a 100 kb genome binned at 10 kb.

    PRIMER1 -> bins 0,1,2   PRIMER2 -> bins 3,4
    PRIMER3 -> bins 1,3     PRIMER4 -> bin 5
    """
    positions_map = {
        "PRIMER1": np.array([5000, 15000, 25000]),
        "PRIMER2": np.array([35000, 45000]),
        "PRIMER3": np.array([15000, 35000]),
        "PRIMER4": np.array([55000]),
    }
    cache = Mock()
    cache.get_positions = lambda prefix, primer, strand: positions_map.get(primer, np.array([]))
    return DominatingSetOptimizer(
        cache=cache,
        fg_prefixes=["genome1"],
        fg_seq_lengths=[genome_length],
        bin_size=bin_size,
        extension_reach=extension_reach,
    )


CANDIDATES = ["PRIMER1", "PRIMER2", "PRIMER3", "PRIMER4"]


class TestFixedBudgetFormulation:
    """The budget is a constraint, not a source of infeasibility."""

    def test_budget_below_minimum_cover_returns_a_solution(self):
        """One primer cannot cover the reachable genome; that is not an error.

        Covering every reachable bin needs three primers. The old formulation
        called a budget of one infeasible and returned None, which is why no
        caller could use it to bound a 6- or 12-primer design.
        """
        result = _optimizer().optimize_ilp(CANDIDATES, max_primers=1, verbose=False)

        assert result is not None
        assert result["n_primers"] == 1

    def test_picks_the_primer_covering_the_most_bases(self):
        """At a budget of one, PRIMER1 (3 bins) beats PRIMER2 (2 bins)."""
        result = _optimizer().optimize_ilp(CANDIDATES, max_primers=1, verbose=False)

        assert result["primers"] == ["PRIMER1"]

    def test_never_exceeds_the_budget(self):
        result = _optimizer().optimize_ilp(CANDIDATES, max_primers=2, verbose=False)

        assert len(result["primers"]) <= 2

    def test_reports_the_solve_as_proven_optimal(self):
        result = _optimizer().optimize_ilp(CANDIDATES, max_primers=2, verbose=False)

        assert result["proven_optimal"] is True


class TestItBoundsTheGreedyResult:
    """The point of the exact solve: an upper bound greedy cannot beat."""

    @pytest.mark.parametrize("budget", [1, 2, 3])
    def test_optimum_is_at_least_the_greedy_coverage(self, budget):
        optimizer = _optimizer()

        greedy = optimizer.optimize_greedy(CANDIDATES, max_primers=budget, verbose=False)
        exact = optimizer.optimize_ilp(CANDIDATES, max_primers=budget, verbose=False)

        assert exact["coverage"] >= greedy["coverage"] - 1e-9

    def test_lp_relaxation_bounds_the_integer_optimum(self):
        optimizer = _optimizer()

        exact = optimizer.optimize_ilp(CANDIDATES, max_primers=2, verbose=False)
        bound = optimizer.coverage_upper_bound(CANDIDATES, max_primers=2)

        assert bound >= exact["coverage"] - 1e-9


class TestCoverageUnits:
    """Coverage is base-weighted, and it is a fraction of the genome."""

    def test_short_final_bin_is_not_counted_as_a_full_one(self):
        """A 15 kb genome at 10 kb bins has a full bin and a 5 kb one.

        Covering only the short bin is 1 of 2 bins but 5000 of 15000 bases.
        Counting bins would read 0.5 -- the same units confusion measured at
        5.3% on the shipped plasmid example.
        """
        cache = Mock()
        cache.get_positions = lambda prefix, primer, strand: np.array([12000])
        optimizer = DominatingSetOptimizer(
            cache=cache,
            fg_prefixes=["genome1"],
            fg_seq_lengths=[15000],
            bin_size=10000,
        )

        result = optimizer.optimize_ilp(["TAILBIN"], max_primers=1, verbose=False)

        assert result["coverage"] == pytest.approx(5000 / 15000)

    def test_coverage_is_a_fraction_of_the_genome_not_of_the_reachable_pool(self):
        """Bins 6-9 are unreachable by any candidate and still count against us."""
        result = _optimizer().optimize_ilp(CANDIDATES, max_primers=4, verbose=False)

        # Every candidate selected reaches bins 0-5 only: 60 kb of a 100 kb genome.
        assert result["coverage"] == pytest.approx(0.6)


class TestExtensionReachIsHonoured:
    """The bound must be computed over the bins the greedy optimizer sees."""

    def test_reach_increases_the_achievable_coverage(self):
        """Without the reach the ILP bounds a different, smaller problem.

        `optimize_greedy` passes `extension_reach` into the bipartite graph;
        the old `optimize_ilp` did not, so its bins were binding sites alone.
        With the reduced bin size that reach implies, ignoring it made the
        reported coverage go *down* when reach was switched on.
        """
        without = _optimizer(extension_reach=0)
        with_reach = _optimizer(extension_reach=20000)

        flat = without.optimize_ilp(CANDIDATES, max_primers=3, verbose=False)
        extended = with_reach.optimize_ilp(CANDIDATES, max_primers=3, verbose=False)

        assert extended["coverage"] > flat["coverage"]


def _conflicting_optimizer():
    pair = ["AAGGTGCGAATA", "TATTCGCACCTT"]
    cache = Mock()
    cache.get_positions = lambda prefix, primer, strand: np.array(
        [5000 if primer == pair[0] else 15000]
    )
    return DominatingSetOptimizer(cache, ["fg"], [20000], max_dimer_bp=3), pair


def test_exact_solver_enforces_dimer_constraints():
    opt, pair = _conflicting_optimizer()
    strict = opt.optimize_ilp(pair, max_primers=2, verbose=False)
    coverage_only = opt.optimize_ilp(pair, max_primers=2, verbose=False, enforce_dimers=False)
    assert strict["coverage"] == pytest.approx(0.5)
    assert strict["n_primers"] == 1
    assert coverage_only["coverage"] == pytest.approx(1.0)
    assert opt.coverage_upper_bound(pair, max_primers=2) == pytest.approx(0.5)


def test_fixed_primers_are_included_in_total_budget():
    opt, pair = _conflicting_optimizer()
    result = opt.optimize_ilp([pair[0]], max_primers=1, fixed_primers=[pair[1]], verbose=False)
    assert result["primers"] == [pair[1]]


def test_incompatible_fixed_panel_reports_infeasibility_without_fake_zero():
    opt, pair = _conflicting_optimizer()
    result = opt.optimize_ilp(pair, max_primers=2, fixed_primers=pair, verbose=False)
    assert not result["feasible"]
    assert result["coverage"] is None
    assert result["coverage_upper_bound"] is None
    assert result["status"] in {"INFEASIBLE", "INT_INFEASIBLE"}


def test_upper_bound_uses_solver_bound_not_incumbent(monkeypatch):
    opt = _optimizer()
    monkeypatch.setattr(
        opt,
        "_solve_max_coverage",
        lambda *a, **k: {
            "coverage": 0.3,
            "coverage_upper_bound": 0.8,
            "proven_optimal": False,
        },
    )
    assert opt.coverage_upper_bound(CANDIDATES) == 0.8


def test_zero_coverage_fixed_primer_still_consumes_budget():
    opt, pair = _conflicting_optimizer()
    opt.cache.get_positions = lambda prefix, primer, strand: np.array([])
    result = opt.optimize_ilp(pair, max_primers=0, fixed_primers=[pair[0]], verbose=False)
    assert not result["feasible"]


def test_independent_benchmark_matches_constrained_library():
    import importlib.util
    import sys
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "scripts/benchmarking/max_coverage_bound.py"
    spec = importlib.util.spec_from_file_location("coverage_benchmark_test", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    opt, pair = _conflicting_optimizer()
    benchmark = module.coverage_bounds(opt, pair, budget=2, max_seconds=10)
    library = opt.optimize_ilp(pair, max_primers=2, verbose=False)
    assert benchmark["ilp"].coverage == pytest.approx(library["coverage"])
    assert benchmark["lp"].coverage_upper_bound == pytest.approx(0.5)
