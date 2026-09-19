"""A pool where every candidate is always bound cannot tell a mismatch.

Known Issue 17. The entry claimed two things and measurement on 2026-09-19
supports only one of them.

**Not supported: that the Tm floor pads the pool with primers that are almost
never bound.** phi29's default `primer_tm_range` floor does sit 10 C below the
reaction temperature, but at k = 12 there is almost nothing down there to
admit: of 40,000 random 12-mers, 8 fall in the Tm 20-25 band, which is 0.02%.
Moving the floor changes essentially nothing.

**Supported, and severe: that nothing excludes saturated primers.** At phi29
30 C, 85.6% of random 12-mers sit above 0.99 occupancy and mean discrimination
-- matched over single-mismatch occupancy -- is 1.065. A primer that binds its
true site 1.07 times better than a near-miss is not selecting anything.

**But a gate is the wrong remedy, and that was measured too.** Capping
candidate occupancy on the real Wolbachia pool made the delivered panel worse
on both axes: at n=12, coverage 0.7334 -> 0.4397 and selectivity density
25.62 -> 6.60, with host sites RISING from 149 to 237. Discrimination is
concentrated in a tail too small to build a panel from, and forcing selection
into it costs more than it returns.

What does move the pool is chemistry, which is what this module reports:

| reaction | occupancy > 0.99 | mean discrimination |
|---|---|---|
| phi29 30 C | 85.6% | 1.065 |
| phi29 30 C, DMSO 10% + betaine 1.5 M | 53.5% | 1.374 |
| equiphi29 42 C | 24.2% | 2.190 |

So the honest output is a measurement the user can act on -- change the
polymerase or add an additive -- not a filter that silently shrinks the pool.
"""

import pytest

from neoswga.core.occupancy import DISCRIMINATION_FLOOR, discrimination_profile
from neoswga.core.reaction_conditions import ReactionConditions


def _random_12mers(n=300, seed=11):
    """A representative pool, because a hand-picked one is not.

    The first draft of this file used six chosen 12-mers including
    `TTTTTTTTTTTT` and `ATATATATATAT`. Those melt far below the reaction
    temperature and so discriminate beautifully: the fixture measured 2.92
    where 40,000 random 12-mers measure 1.065, and every assertion about the
    saturated regime failed on a pool that was not in it.
    """
    import random

    rng = random.Random(seed)
    out, seen = [], set()
    while len(out) < n:
        candidate = "".join(rng.choice("ACGT") for _ in range(12))
        if candidate in seen:
            continue
        seen.add(candidate)
        out.append(candidate)
    return out


PRIMERS = _random_12mers()


class TestItMeasuresTheRegime:
    def test_it_reports_a_mean_discrimination(self):
        profile = discrimination_profile(PRIMERS, ReactionConditions(temp=30.0))

        assert profile.mean_discrimination > 0
        assert profile.n == len(PRIMERS)

    def test_it_reports_the_saturated_fraction(self):
        profile = discrimination_profile(PRIMERS, ReactionConditions(temp=30.0))

        assert 0.0 <= profile.saturated_fraction <= 1.0

    def test_a_warm_reaction_discriminates_better_than_a_cold_one(self):
        """The finding the whole entry turns on: the same primers discriminate
        far better at 42 C than at 30 C, because occupancy is no longer pinned
        at 1. Nothing about the sequences changed.

        `ReactionConditions` validates the temperature against the polymerase,
        so the warm case needs the polymerase that actually runs there -- 55 C
        under phi29 raises rather than measuring anything.
        """
        cold = discrimination_profile(PRIMERS, ReactionConditions(temp=30.0))
        warm = discrimination_profile(
            PRIMERS, ReactionConditions(temp=42.0, polymerase="equiphi29")
        )

        assert warm.mean_discrimination > cold.mean_discrimination
        assert warm.saturated_fraction < cold.saturated_fraction

    def test_an_empty_pool_is_not_an_exception(self):
        profile = discrimination_profile([], ReactionConditions(temp=30.0))

        assert profile.n == 0
        assert profile.mean_discrimination == 0.0
        assert profile.saturated is False

    def test_without_conditions_it_declines_to_guess(self):
        """No reaction temperature means no occupancy, and a fabricated number
        here would read exactly like a measured one."""
        profile = discrimination_profile(PRIMERS, None)

        assert profile.n == 0
        assert profile.saturated is False


class TestItSaysWhenThePoolCannotDiscriminate:
    def test_a_cold_phi29_pool_is_flagged(self):
        """85.6% of random 12-mers above 0.99 occupancy, mean discrimination
        1.065. This is the shipped default for phi29."""
        profile = discrimination_profile(PRIMERS, ReactionConditions(temp=30.0))

        assert profile.saturated is True
        assert profile.mean_discrimination < DISCRIMINATION_FLOOR

    def test_a_warm_pool_is_not(self):
        profile = discrimination_profile(
            PRIMERS, ReactionConditions(temp=42.0, polymerase="equiphi29")
        )

        assert profile.saturated is False

    def test_the_advice_names_the_lever_that_was_measured(self):
        """Chemistry, because that is what the measurement supports. It must
        not suggest tightening the Tm window, which was measured and makes the
        delivered panel worse."""
        profile = discrimination_profile(PRIMERS, ReactionConditions(temp=30.0))

        advice = profile.advice.lower()
        assert "polymerase" in advice
        assert "additive" in advice
        # It may MENTION the Tm window -- it warns against narrowing it -- but
        # it must not offer it as the remedy.
        assert "narrow" in advice and "worse panel" in advice

    def test_a_pool_that_discriminates_offers_no_advice(self):
        profile = discrimination_profile(
            PRIMERS, ReactionConditions(temp=42.0, polymerase="equiphi29")
        )

        assert profile.advice == ""


class TestTheFloorIsWhereMeasurementPutIt:
    def test_it_sits_between_the_two_measured_regimes(self):
        """phi29 at 30 C measures 1.065 and equiphi29 at 42 C measures 2.190,
        so a floor between them separates the case worth reporting from the
        case that is fine. Pinned so a later edit has to argue with the
        numbers."""
        assert 1.065 < DISCRIMINATION_FLOOR < 2.190
