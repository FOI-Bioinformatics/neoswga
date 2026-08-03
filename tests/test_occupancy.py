"""Per-site occupancy: the term that lets conditions affect specificity.

Selectivity in this project is a count of exact k-mer matches
(`base_optimizer.py:777`), so every site contributes 1.0 no matter how stable
that duplex is, and no reaction condition can change the answer. Additives move
Tm, and Tm is not in the formula.

This module supplies the missing term. For a duplex melting at `tm` under the
configured conditions, at reaction temperature `temp`, the two-state fraction
bound is

    theta = 1 / (1 + exp( (dH/R) * (1/T - 1/Tm) ))

with `dH` from the nearest-neighbour model and both temperatures in Kelvin. At
`T == Tm` it is 0.5 by construction, which is what Tm means.

Nothing here is new thermodynamics. `dH` comes from
`thermodynamics.calculate_enthalpy_entropy` and the additive- and salt-corrected
Tm from `ReactionConditions.calculate_effective_tm`; occupancy is layered onto
corrections that already exist and were only ever used per-primer.

The property these tests exist for is the last one: **the gap between a perfect
match and a mismatched one widens as stringency rises.** That is the entire
mechanism by which an additive can buy specificity — mismatched background sites
fall off the sigmoid faster than perfect foreground ones. If it does not hold,
nothing downstream can work, so it is asserted directly rather than inferred
from a primer set.
"""

import pytest

from neoswga.core.occupancy import mismatch_tm, site_occupancy
from neoswga.core.thermodynamics import calculate_enthalpy_entropy

PRIMER = "ATCGATCGAT"
DH, _DS = calculate_enthalpy_entropy(PRIMER)


# ----------------------------------------------------------------------
# Shape
# ----------------------------------------------------------------------


def test_half_bound_at_the_melting_temperature():
    """The definition of Tm, and the anchor for everything else."""
    assert site_occupancy(DH, tm=45.0, temp=45.0) == pytest.approx(0.5)


def test_bounded_to_a_probability():
    for temp in (-20.0, 0.0, 30.0, 45.0, 60.0, 100.0):
        theta = site_occupancy(DH, tm=45.0, temp=temp)
        assert 0.0 <= theta <= 1.0


def test_falls_with_temperature():
    """Heat melts duplexes. A model where it did not would be unusable."""
    temps = [20.0, 30.0, 40.0, 50.0, 60.0]
    thetas = [site_occupancy(DH, tm=45.0, temp=t) for t in temps]

    for hotter, cooler in zip(thetas, thetas[1:]):
        assert cooler < hotter


def test_rises_with_melting_temperature():
    """A more stable duplex is more occupied at fixed temperature."""
    tms = [30.0, 40.0, 50.0, 60.0]
    thetas = [site_occupancy(DH, tm=tm, temp=42.0) for tm in tms]

    for lower, higher in zip(thetas, thetas[1:]):
        assert higher > lower


def test_a_stronger_duplex_is_sharper():
    """Larger |dH| gives a steeper transition.

    This is why the discrimination below works at all: a long, high-enthalpy
    duplex holds on while a weak one has already let go.
    """
    weak_dh, _ = calculate_enthalpy_entropy("ATCGAT")
    strong_dh, _ = calculate_enthalpy_entropy("GGGGCCCCGGGGCCCC")

    above = 5.0
    weak = site_occupancy(weak_dh, tm=45.0, temp=45.0 + above)
    strong = site_occupancy(strong_dh, tm=45.0, temp=45.0 + above)

    assert strong < weak, "the higher-enthalpy duplex should fall off faster past Tm"


def test_zero_enthalpy_is_handled_rather_than_dividing_by_nothing():
    """A degenerate input must not raise; an unresolvable duplex is simply
    uninformative, and 0.5 is the honest answer for no information."""
    assert site_occupancy(0.0, tm=45.0, temp=30.0) == pytest.approx(0.5)


# ----------------------------------------------------------------------
# Mismatches
# ----------------------------------------------------------------------


def test_a_mismatch_lowers_the_melting_temperature():
    assert mismatch_tm(50.0, n_mismatches=1, penalty=4.0) == 46.0
    assert mismatch_tm(50.0, n_mismatches=2, penalty=4.0) == 42.0


def test_a_perfect_match_is_unchanged():
    assert mismatch_tm(50.0, n_mismatches=0, penalty=4.0) == 50.0


def test_the_penalty_comes_from_configuration():
    """`mismatch_penalty` has been a params.json key with no reader at all --
    three occurrences in the package, none of them a use. This is its first
    consumer and its intended meaning."""
    from neoswga.core.occupancy import default_mismatch_penalty

    assert default_mismatch_penalty() > 0


# ----------------------------------------------------------------------
# The thesis
# ----------------------------------------------------------------------


def test_stringency_widens_the_gap_between_perfect_and_mismatched():
    """The mechanism the whole plan rests on.

    A perfect foreground site and a mismatched background site both lose
    occupancy as temperature rises, but the mismatched one loses it faster.
    The ratio between them is what an additive or a hotter reaction buys, and
    it has to grow monotonically with stringency.
    """
    perfect_tm = 50.0
    mismatched_tm = mismatch_tm(perfect_tm, n_mismatches=1, penalty=4.0)

    ratios = []
    for temp in (30.0, 35.0, 40.0, 45.0, 50.0):
        perfect = site_occupancy(DH, tm=perfect_tm, temp=temp)
        mismatched = site_occupancy(DH, tm=mismatched_tm, temp=temp)
        ratios.append(perfect / mismatched)

    for cooler, hotter in zip(ratios, ratios[1:]):
        assert hotter > cooler, f"discrimination did not improve with stringency: {ratios}"


def test_lowering_tm_with_an_additive_has_the_same_effect_as_heating():
    """Additives act by lowering effective Tm, so they must buy discrimination
    the same way raising the temperature does. This is what makes
    `--dmso-percent` a specificity control rather than a Tm decoration.
    """
    temp = 30.0
    perfect_tm = 50.0

    def discrimination(tm_shift):
        perfect = site_occupancy(DH, tm=perfect_tm - tm_shift, temp=temp)
        mismatched = site_occupancy(DH, tm=mismatch_tm(perfect_tm - tm_shift, 1, 4.0), temp=temp)
        return perfect / mismatched

    # 0 C shift is no additive; 8 C is roughly 15% DMSO or 3 M betaine.
    assert discrimination(8.0) > discrimination(0.0)


def test_discrimination_is_worthless_when_everything_is_saturated():
    """The other side of the trade, and the reason an optimum is interior.

    Far below Tm both perfect and mismatched sites are fully bound, so there is
    no discrimination to be had -- which is why a 30 C phi29 reaction with
    12-mers cannot be made selective by chemistry alone.
    """
    perfect_tm = 70.0
    perfect = site_occupancy(DH, tm=perfect_tm, temp=20.0)
    mismatched = site_occupancy(DH, tm=mismatch_tm(perfect_tm, 1, 4.0), temp=20.0)

    assert perfect / mismatched < 1.05
