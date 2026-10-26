"""Swap search improves a greedy local choice without relaxing compatibility."""

from types import SimpleNamespace

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


# ----------------------------------------------------------------------
# The host term as a weighted axis (Known Issue 14)
# ----------------------------------------------------------------------
#
# Without `background_weight` the host load is the lexicographic tie-break the
# score has always carried: it only decides between panels of EQUAL foreground
# coverage. With `background_weight > 0` the host term instead becomes a
# weighted axis, so a swap to a less-covering but much-less-host-binding primer
# is taken even though it loses foreground coverage. Both coverages are 0..1
# fractions, so the weight is a genuine trade. The feature ships off; these pin
# that off is unchanged and on does what it says. The score runs only when an
# `objective` is supplied, so these use a minimal stub.


class _TradeObjective:
    """Feasible everywhere, with a foreground and a host coverage per primer so
    a swap trades one against the other. ``total_bg_sites`` is the lexicographic
    fallback's tie-breaker, kept distinct from ``bg_coverage`` so a test can
    tell which axis decided."""

    #: primer -> (fg_coverage, bg_coverage, total_bg_sites)
    _TABLE = {"A": (0.60, 0.90, 9), "B": (0.50, 0.10, 1)}

    def shortfall(self, panel):
        return 0.0

    def coverage(self, panel):
        return max(self._TABLE[p][0] for p in panel)

    def metrics(self, panel):
        return SimpleNamespace(
            bg_coverage=max(self._TABLE[p][1] for p in panel),
            total_bg_sites=sum(self._TABLE[p][2] for p in panel),
        )


def _trade_run(background_weight):
    return refine_by_swaps(
        ["A"],
        ["A", "B"],
        None,
        None,
        Compatibility(),
        objective=_TradeObjective(),
        background_weight=background_weight,
    )


def test_off_by_default_keeps_the_higher_coverage_primer():
    # B covers less (0.50 < 0.60), so on pure coverage A is never dropped --
    # the host difference is invisible to the lexicographic score.
    assert _trade_run(0.0).primers == ("A",)


def test_a_weighted_host_term_trades_coverage_for_a_lower_host_panel():
    # 0.60 - 0.8*0.90 = -0.12 for A; 0.50 - 0.8*0.10 = 0.42 for B. The weighted
    # axis now prefers B even though it covers less foreground.
    assert _trade_run(0.8).primers == ("B",)


class _UnmeasuredHostObjective(_TradeObjective):
    """Host coverage cannot be measured (no background genome)."""

    def metrics(self, panel):
        return SimpleNamespace(
            bg_coverage=None,
            total_bg_sites=sum(self._TABLE[p][2] for p in panel),
        )


def test_unmeasurable_host_coverage_falls_back_without_inventing_a_zero():
    # background_weight > 0 but bg_coverage is None: the weighted form is NOT
    # applied (that would be `coverage - w*None`, a TypeError, or a fabricated
    # zero). It falls back to the lexicographic tuple, which on unequal coverage
    # keeps the higher-coverage A -- and, crucially, does not raise.
    result = refine_by_swaps(
        ["A"],
        ["A", "B"],
        None,
        None,
        Compatibility(),
        objective=_UnmeasuredHostObjective(),
        background_weight=0.8,
    )
    assert result.primers == ("A",)
