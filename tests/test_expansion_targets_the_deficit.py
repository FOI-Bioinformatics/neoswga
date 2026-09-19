"""Expansion should chase the depth a region is missing, not gap membership.

Finding F7 of the 2026-09-16 pipeline audit. Three defects, one cause: the
objective asked whether a binding SITE fell inside a gap rather than how much
of the missing depth a candidate's WINDOW would recover.

- **The reach was never consulted.** A candidate binding just outside a 50 kb
  gap, whose modelled window would blanket it, was discarded. One binding at
  the gap's last base and extending away was kept.
- **The filter was all-or-nothing.** When fewer than `target_new` candidates
  survived it, the whole filter was abandoned for the unfiltered pool, so the
  run silently stopped targeting gaps at all.
- **`gap_coverage` counted gaps, not bases.** `1 - len(after)/len(before)` goes
  NEGATIVE when one long gap splits into two, which is progress reported as
  regression.

A weight per base replaces all three: `max(0, 1 - depth/desired)`, zero where
the depth was never evaluable, because absent evidence is not a deficit. That
last point is what connects this to the record binding: a record the BAM does
not cover has unknown depth, and chasing it would design primers for a region
nothing measured.
"""

import numpy as np
import pytest

from neoswga.core.deficit_objective import (
    covered_mask,
    deficit_gain,
    deficit_weights,
    recovered_deficit,
)


class TestTheWeights:
    def test_full_depth_is_no_deficit(self):
        weights = deficit_weights(np.full(10, 30), desired_depth=30)

        assert weights.sum() == 0.0

    def test_more_than_enough_depth_is_not_a_negative_deficit(self):
        weights = deficit_weights(np.full(10, 100), desired_depth=30)

        assert (weights >= 0).all()
        assert weights.sum() == 0.0

    def test_no_depth_is_a_full_deficit(self):
        weights = deficit_weights(np.zeros(10), desired_depth=30)

        assert weights.tolist() == [1.0] * 10

    def test_half_the_depth_is_half_the_deficit(self):
        weights = deficit_weights(np.full(4, 15), desired_depth=30)

        assert weights.tolist() == pytest.approx([0.5] * 4)

    def test_an_unevaluable_base_carries_no_deficit(self):
        """The link to the record binding: a record the BAM does not cover has
        UNKNOWN depth. Treating unknown as zero would design primers for a
        region nothing measured, which is the defect, not the fix."""
        depth = np.zeros(6)
        evaluable = np.array([True, True, True, False, False, False])

        weights = deficit_weights(depth, desired_depth=30, evaluable=evaluable)

        assert weights.tolist() == [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]

    def test_a_desired_depth_of_zero_is_refused(self):
        with pytest.raises(ValueError, match="desired_depth"):
            deficit_weights(np.zeros(4), desired_depth=0)


class TestRecoveryIsMeasuredInBases:
    def _weights(self, length=1000, deficit_span=(400, 600)):
        weights = np.zeros(length)
        weights[deficit_span[0] : deficit_span[1]] = 1.0
        return weights

    def test_a_window_over_the_deficit_recovers_it(self):
        weights = self._weights()

        recovered = recovered_deficit(
            positions=[500], weights=weights, extension=50, length=1000, circular=False
        )

        assert recovered == pytest.approx(100.0), "450..550 is 100 deficit bases"

    def test_covering_nothing_recovers_nothing(self):
        weights = self._weights()

        recovered = recovered_deficit(
            positions=[50], weights=weights, extension=10, length=1000, circular=False
        )

        assert recovered == 0.0

    def test_two_primers_do_not_double_count_an_overlap(self):
        """A union, not a sum. Two windows over the same deficit recover it
        once, which is what makes the figure a number of bases."""
        weights = self._weights()

        one = recovered_deficit([500], weights, extension=50, length=1000, circular=False)
        both = recovered_deficit(
            [495, 505], weights, extension=50, length=1000, circular=False
        )

        assert both < 2 * one

    def test_recovery_can_never_be_negative(self):
        """`1 - len(after)/len(before)` went negative when one gap split into
        two. A count of recovered bases cannot."""
        weights = self._weights()

        for positions in ([], [0], [999], [400, 600]):
            assert (
                recovered_deficit(positions, weights, extension=30, length=1000, circular=False)
                >= 0.0
            )


class TestTheReachIsConsulted:
    """The F7 headline: a site outside a gap whose WINDOW blankets it."""

    def _weights(self):
        weights = np.zeros(100_000)
        weights[20_000:70_000] = 1.0  # a 50 kb deficit
        return weights

    def test_a_site_outside_the_gap_still_recovers_it(self):
        """Binding at 19,000 with a 3 kb reach reaches 22,000, so it recovers
        2,000 deficit bases. The old screen discarded this candidate for
        binding outside the gap."""
        recovered = recovered_deficit(
            [19_000], self._weights(), extension=3_000, length=100_000, circular=False
        )

        assert recovered == pytest.approx(2_000.0)

    def test_a_site_at_the_edge_extending_away_recovers_little(self):
        """Binding at the gap's last base and extending outward was KEPT by the
        old screen, though it recovers almost nothing."""
        recovered = recovered_deficit(
            [69_999], self._weights(), extension=3_000, length=100_000, circular=False
        )

        assert recovered < 3_100

    def test_the_two_are_ranked_by_what_they_recover(self):
        """Not by whether the site happens to sit inside the interval."""
        weights = self._weights()
        outside = recovered_deficit([19_000], weights, 3_000, 100_000, False)
        deep_inside = recovered_deficit([45_000], weights, 3_000, 100_000, False)

        assert deep_inside > outside


class TestWindowsRespectRecordBoundaries:
    def test_a_window_does_not_reach_into_the_next_record(self):
        """Two molecules. A polymerase cannot travel between them, so a primer
        near the end of one cannot recover deficit in the other."""
        weights = np.zeros(2000)
        weights[1000:1200] = 1.0  # deficit lives entirely in the second record

        recovered = recovered_deficit(
            [980],
            weights,
            extension=500,
            length=2000,
            circular=False,
            record_starts=[0, 1000],
        )

        assert recovered == 0.0

    def test_without_record_starts_the_window_spans_the_join(self):
        """The behaviour that made this necessary, pinned so the difference is
        visible rather than assumed."""
        weights = np.zeros(2000)
        weights[1000:1200] = 1.0

        recovered = recovered_deficit(
            [980], weights, extension=500, length=2000, circular=False
        )

        assert recovered > 0.0


class TestMarginalGain:
    def _weights(self):
        weights = np.zeros(1000)
        weights[400:600] = 1.0
        return weights

    def test_a_candidate_adding_nothing_new_gains_nothing(self):
        weights = self._weights()
        pool = covered_mask([500], extension=100, length=1000, circular=False)

        assert deficit_gain([500], pool, weights, 100, 1000, False) == 0.0

    def test_a_candidate_reaching_untouched_deficit_gains_it(self):
        weights = self._weights()
        pool = covered_mask([450], extension=25, length=1000, circular=False)

        gain = deficit_gain([560], pool, weights, 25, 1000, False)

        assert gain == pytest.approx(50.0)

    def test_gain_is_never_negative(self):
        weights = self._weights()
        pool = covered_mask([500], extension=300, length=1000, circular=False)

        assert deficit_gain([10], pool, weights, 5, 1000, False) >= 0.0


class TestTheFilterIsAPrescreenNotAGate:
    def test_the_prescreen_dilates_by_the_reach(self):
        """A candidate binding within one reach of a gap can blanket part of
        it, so the prescreen must keep it. The old screen tested membership."""
        from neoswga.core.deficit_objective import dilate_intervals

        assert dilate_intervals([(20_000, 70_000)], reach=3_000, length=100_000) == [
            (17_000, 73_000)
        ]

    def test_dilation_is_clamped_to_the_sequence(self):
        from neoswga.core.deficit_objective import dilate_intervals

        assert dilate_intervals([(10, 50)], reach=100, length=200) == [(0, 150)]

    def test_overlapping_dilated_intervals_merge(self):
        from neoswga.core.deficit_objective import dilate_intervals

        assert dilate_intervals([(0, 10), (15, 20)], reach=5, length=100) == [(0, 25)]

    def test_no_intervals_dilate_to_nothing(self):
        from neoswga.core.deficit_objective import dilate_intervals

        assert dilate_intervals([], reach=5, length=100) == []
