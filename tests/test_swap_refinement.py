"""Swap search improves a greedy local choice without relaxing compatibility."""

import pytest

from neoswga.core.swap_refinement import refine_by_swaps


class Compatibility:
    def __init__(self, conflicts=()):
        self.conflicts = {frozenset(pair) for pair in conflicts}

    def dimerises(self, primer, others):
        return any(frozenset((primer, p)) in self.conflicts for p in others)


def run(selected=("A", "B"), **kwargs):
    bins = {"A": {0, 1}, "B": {1, 2}, "C": {2, 3}, "D": {4}}
    return refine_by_swaps(
        selected,
        list(bins),
        bins,
        {b: 10 for b in range(5)},
        kwargs.pop("dimers", Compatibility()),
        **kwargs,
    )


def test_swap_recovers_coverage_from_an_excluded_candidate():
    result = run()
    assert result.covered_bases == 40
    assert result.swaps == 1
    assert result.primers == ("A", "C")


def test_fixed_primers_remain_selected():
    result = run(fixed_primers=("A", "B"))
    assert result.primers == ("A", "B")
    assert result.evaluations == 0


def test_incompatible_replacements_are_rejected():
    result = run(dimers=Compatibility([("A", "C"), ("B", "C")]))
    assert "C" not in result.primers


def test_background_breaks_equal_coverage_ties():
    bins = {"A": {0}, "B": {0}}
    result = refine_by_swaps(
        ["A"], ["A", "B"], bins, {0: 5}, Compatibility(), background_sites={"A": 10, "B": 1}
    )
    assert result.primers == ("B",)
    assert result.background_sites == 1


@pytest.mark.parametrize("budget", [0, 1, 3])
def test_evaluation_budget_is_respected(budget):
    result = run(max_evaluations=budget)
    assert result.evaluations <= budget
    assert result.covered_bases >= 30
    assert result.stop_reason == "evaluation_limit"


def test_zero_time_returns_original_panel():
    result = run(max_seconds=0)
    assert result.primers == ("A", "B")
    assert result.stop_reason == "time_limit"


def test_short_bins_use_base_weights():
    bins = {"A": {0}, "B": {1, 2}}
    result = refine_by_swaps(["A"], ["A", "B"], bins, {0: 10, 1: 2, 2: 2}, Compatibility())
    assert result.primers == ("A",)


def test_missing_fixed_primer_is_rejected():
    with pytest.raises(ValueError, match="Fixed primers"):
        run(fixed_primers=["C"])
