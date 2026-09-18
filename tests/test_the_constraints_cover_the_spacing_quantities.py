"""A user should be able to constrain the properties nobody had a line for.

`docs/validation/getting_ahead_on_spacing_2026-09-18.md` rules out NeoSWGA
picking a spacing threshold: no reach-derived limit separates the 18 published
sets with wet-lab outcomes, and a fitted weight is wrong for one of the two
benchmarks either way. Letting a USER set one is a different claim. They know
their target, their regime and what they will trade, and item 1's report now
tells them which properties had no reference at all.

So these limits exist and every one of them is unset by default. A new limit
that changed a delivered panel without being asked for would be the scoring
change both benchmarks refuse.

Two quantities are deliberately NOT constrainable, and
`test_the_strand_metrics_are_not_constrainable` pins that.
"""

import math

import pytest

from neoswga.core.pool_objective import PoolConstraints, PoolObjective


class _Metrics:
    def __init__(
        self,
        *,
        effective_fg_coverage=0.85,
        fg_coverage=0.90,
        selectivity_density=45.0,
        total_bg_sites=120,
        max_gap=31_000.0,
        mean_gap=4_800.0,
        gap_gini=0.48,
        bg_coverage=0.02,
    ):
        self.effective_fg_coverage = effective_fg_coverage
        self.fg_coverage = fg_coverage
        self.selectivity_density = selectivity_density
        self.total_bg_sites = total_bg_sites
        self.max_gap = max_gap
        self.mean_gap = mean_gap
        self.gap_gini = gap_gini
        self.bg_coverage = bg_coverage


def _objective(metrics=None, **limits):
    constraints = PoolConstraints(**limits)
    return PoolObjective(lambda _primers: metrics or _Metrics(), constraints)


PRIMERS = ["ACGTACGTAC", "TTGCATGCAT"]


class TestNothingIsLimitedByDefault:
    """The whole set of new limits has to be inert unless asked for."""

    def test_an_unconstrained_panel_violates_nothing(self):
        assert _objective().violations(PRIMERS) == ()

    def test_an_unconstrained_panel_has_no_shortfall(self):
        assert _objective().shortfall(PRIMERS) == 0.0

    @pytest.mark.parametrize(
        "field",
        ["max_worst_hole", "max_mean_gap", "max_evenness", "max_host_coverage"],
    )
    def test_each_new_limit_defaults_to_unset(self, field):
        assert getattr(PoolConstraints(), field) is None


class TestEachSpacingLimitBinds:
    @pytest.mark.parametrize(
        "limits,expected",
        [
            ({"max_worst_hole": 20_000.0}, "worst hole above maximum"),
            ({"max_mean_gap": 3_000.0}, "mean gap above maximum"),
            ({"max_evenness": 0.30}, "evenness above maximum"),
            ({"max_host_coverage": 0.01}, "host coverage above maximum"),
        ],
    )
    def test_a_panel_over_the_limit_is_named_as_violating(self, limits, expected):
        assert expected in _objective(**limits).violations(PRIMERS)

    @pytest.mark.parametrize(
        "limits",
        [
            {"max_worst_hole": 40_000.0},
            {"max_mean_gap": 6_000.0},
            {"max_evenness": 0.60},
            {"max_host_coverage": 0.05},
        ],
    )
    def test_a_panel_inside_the_limit_violates_nothing(self, limits):
        assert _objective(**limits).violations(PRIMERS) == ()


class TestShortfallStaysConsistentWithViolations:
    """The two must agree at the boundary or a feasible panel can rank behind an
    infeasible one, which is the ordering `shortfall` exists to prevent."""

    @pytest.mark.parametrize(
        "limits",
        [
            {"max_worst_hole": 20_000.0},
            {"max_mean_gap": 3_000.0},
            {"max_evenness": 0.30},
            {"max_host_coverage": 0.01},
            {"max_worst_hole": 40_000.0},
            {"max_evenness": 0.60},
        ],
    )
    def test_shortfall_is_zero_exactly_when_nothing_is_violated(self, limits):
        objective = _objective(**limits)

        assert (objective.shortfall(PRIMERS) == 0.0) == (objective.violations(PRIMERS) == ())

    def test_each_term_is_relative_to_its_own_limit(self):
        """A hole over its limit by half must score like evenness over by half,
        or the criterion measured in base pairs dominates every comparison."""
        hole = _objective(_Metrics(max_gap=30_000.0), max_worst_hole=20_000.0)
        evenness = _objective(_Metrics(gap_gini=0.45), max_evenness=0.30)

        assert hole.shortfall(PRIMERS) == pytest.approx(evenness.shortfall(PRIMERS))

    def test_two_failed_limits_are_worse_than_one(self):
        one = _objective(max_worst_hole=20_000.0)
        two = _objective(max_worst_hole=20_000.0, max_evenness=0.30)

        assert two.shortfall(PRIMERS) > one.shortfall(PRIMERS)

    def test_an_unmeasurable_gap_is_infinitely_short_not_merely_large(self):
        """`PrimerSetMetrics.empty` carries `inf` gaps. A panel that could not be
        measured is not nearly acceptable, and a finite value would let coverage
        trade against it."""
        objective = _objective(_Metrics(max_gap=math.inf), max_worst_hole=20_000.0)

        assert objective.shortfall(PRIMERS) == math.inf
        assert "worst hole above maximum" in objective.violations(PRIMERS)


class TestABackgroundLimitNeedsABackground:
    """A limit on a quantity nothing measured would pass every panel, which
    reads as specificity rather than as an absent measurement. Known Issues 5,
    6 and 13."""

    def test_a_host_coverage_limit_requires_one(self):
        constraints = PoolConstraints(max_host_coverage=0.01)

        assert constraints.needs_background is True
        with pytest.raises(ValueError, match="specificity limit"):
            constraints.require_background(available=False)

    def test_a_spacing_limit_on_the_target_does_not(self):
        constraints = PoolConstraints(max_worst_hole=20_000.0, max_evenness=0.5)

        assert constraints.needs_background is False
        constraints.require_background(available=False)


class TestWhatIsDeliberatelyNotConstrainable:
    def test_the_strand_metrics_are_not_constrainable(self):
        """`strand_coverage_ratio` and `strand_alternation_score` are 0.0 both
        when measured zero and when the position cache could not supply them,
        and nothing distinguishes the two. A limit on such a quantity would
        reject a panel for a MISSING MEASUREMENT while reporting a violated
        constraint, which is the failure Known Issues 5, 6 and 13 record.

        Remove this test only alongside a metrics change that makes an
        uncomputed strand statistic distinguishable from a measured zero.
        """
        fields = set(PoolConstraints().__dataclass_fields__)

        assert not {f for f in fields if "strand" in f}

    def test_a_dimer_limit_is_still_outside_the_objective(self):
        """It is a hard constraint on the delivered panel, not a scoring term.
        Folding it in is how it became tradeable and produced an 11 bp
        heterodimer against a configured 3."""
        fields = set(PoolConstraints().__dataclass_fields__)

        assert not {f for f in fields if "dimer" in f}


class TestTheConstraintsStayFrozen:
    def test_a_limit_cannot_be_changed_under_a_running_design(self):
        constraints = PoolConstraints(max_worst_hole=20_000.0)

        with pytest.raises(Exception):
            constraints.max_worst_hole = 10_000.0
