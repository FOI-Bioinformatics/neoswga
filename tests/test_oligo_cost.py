"""What a primer set costs, and why the modification dominates for SWGA.

Set size is the main lever on coverage, so "how many primers" is really "how
many oligos will you buy" -- and the frontier has a missing axis until that is
answerable. The previous estimate was `len(primers) * 5.0`, which ignores both
length and modifications.

Modifications are not a rounding error here. SWGA primers carry 3'
phosphorothioate bonds because phi29's 3'->5' exonuclease degrades unprotected
primers, and a PTO coupling costs roughly ten times a standard base. For the
short primers this tool designs, the fixed and modification charges dominate
and the sequence itself is the smallest component -- the opposite of the
intuition carried over from PCR primers, where length drives price.

That has a design consequence worth stating: buying coverage by adding primers
is close to linear in primer count and barely sensitive to primer length, so
there is little cost argument for preferring 10-mers over 12-mers.

The absolute level is not defensible to better than a factor of a few -- oligo
pricing varies severalfold with vendor, scale and contract -- so every figure
is labelled an estimate. The comparisons are what the model is for, and those
are robust to the level being wrong.
"""

import pytest

from neoswga.core.oligo_cost import (
    OligoCostModel,
    cost_from_modifications,
    cost_per_coverage_point,
    format_cost,
    set_cost,
)

TWELVE = ["ACGTACGTACGT"] * 8


# ----------------------------------------------------------------------
# The arithmetic
# ----------------------------------------------------------------------


def test_cost_is_setup_plus_bases_plus_modifications():
    model = OligoCostModel(
        setup_cost=1.0, per_base_cost=0.1, per_pto_cost=2.0, five_prime_block_cost=10.0
    )

    # 10-mer, 2 PTO bonds: 1.0 + 10*0.1 + 2*2.0 = 6.0
    assert model.oligo_cost("A" * 10, pto_bonds=2) == pytest.approx(6.0)
    assert model.oligo_cost("A" * 10, pto_bonds=2, five_prime_block="c18") == pytest.approx(16.0)


def test_the_breakdown_sums_to_the_total():
    cost = set_cost(TWELVE, pto_bonds=2)

    assert cost.setup + cost.bases + cost.modifications == pytest.approx(cost.total)
    assert cost.per_primer_mean == pytest.approx(cost.total / cost.n_primers)


def test_cost_scales_with_set_size():
    """The lever the frontier needs. Adding primers is the way to buy coverage,
    so its price has to be visible."""
    eight = set_cost(["ACGTACGTACGT"] * 8)
    sixteen = set_cost(["ACGTACGTACGT"] * 16)

    assert sixteen.total == pytest.approx(2 * eight.total)


def test_an_empty_set_costs_nothing_rather_than_raising():
    cost = set_cost([])

    assert cost.total == 0.0
    assert cost.n_primers == 0
    assert "nothing to cost" in format_cost(cost)


# ----------------------------------------------------------------------
# Modifications, which is the point
# ----------------------------------------------------------------------


def test_modifications_dominate_for_short_primers():
    """The fact that makes a length-blind estimate wrong for SWGA.

    At default prices a 10-mer's two PTO bonds cost more than all ten of its
    bases, so an estimate that ignores modifications understates the price and
    misattributes what drives it.
    """
    cost = set_cost(["A" * 10] * 6, pto_bonds=2)

    assert cost.modifications > cost.bases
    assert cost.modification_fraction > 0.4


def test_unmodified_primers_cost_less():
    """`--modification-profile none` is a real option and must be priced as one."""
    modified = set_cost(TWELVE, pto_bonds=2)
    plain = set_cost(TWELVE, pto_bonds=0)

    assert plain.total < modified.total
    assert plain.modifications == 0.0


def test_a_five_prime_block_is_charged():
    """The low-input profile adds a C18 spacer, which is a separate charge."""
    standard = set_cost(TWELVE, pto_bonds=2)
    low_input = set_cost(TWELVE, pto_bonds=2, five_prime_block="c18")

    assert low_input.total > standard.total


def test_pto_bonds_cannot_exceed_the_linkages_a_primer_has():
    """An n-mer has n-1 linkages. Charging for more would silently inflate the
    estimate, and at the 6-8 bp end of the SWGA range asking for more bonds
    than the primer can carry is a reachable mistake rather than a hypothetical
    one."""
    model = OligoCostModel(setup_cost=0.0, per_base_cost=0.0, per_pto_cost=1.0)

    assert model.oligo_cost("ACGTAA", pto_bonds=99) == pytest.approx(5.0)
    assert model.oligo_cost("A", pto_bonds=2) == pytest.approx(0.0)


def test_it_reads_an_existing_modification_profile():
    """The ordering profile and the price of that profile must not drift apart:
    a set exported as low-input has to be costed with the spacer it will be
    ordered with."""
    from neoswga.core.export import PrimerModifications

    plain = cost_from_modifications(TWELVE, PrimerModifications.from_profile("none"))
    standard = cost_from_modifications(TWELVE, PrimerModifications.from_profile("standard"))
    low_input = cost_from_modifications(TWELVE, PrimerModifications.from_profile("low-input"))

    assert plain.total < standard.total < low_input.total
    assert plain.pto_bonds == 0
    assert low_input.five_prime_block == "c18"


def test_primer_length_barely_moves_the_price():
    """The design consequence. Because fixed and modification charges dominate,
    a 12-mer set is only slightly dearer than a 10-mer set -- so cost is not an
    argument for shorter primers, which is where the specificity is worse."""
    ten = set_cost(["A" * 10] * 8)
    twelve = set_cost(["A" * 12] * 8)

    assert twelve.total > ten.total
    assert twelve.total / ten.total < 1.15


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------


def test_the_report_says_it_is_an_estimate():
    """Oligo pricing varies severalfold with vendor and contract. A bare dollar
    figure would be quoted back as though it came from one."""
    text = format_cost(set_cost(TWELVE))

    assert "not a quote" in text
    assert "estimate" in text.lower()
    assert set_cost(TWELVE).to_dict()["estimate"] is True


def test_the_report_shows_the_breakdown_not_just_a_total():
    text = format_cost(set_cost(TWELVE, pto_bonds=2))

    for part in ("setup", "bases", "modifications", "PTO"):
        assert part in text


def test_cost_per_coverage_point_makes_the_frontier_actionable():
    """What the next percentage point costs is what decides whether to buy it."""
    cost = set_cost(["ACGTACGTACGT"] * 10)

    assert cost_per_coverage_point(cost, 0.50) == pytest.approx(cost.total / 50.0)


def test_zero_coverage_sorts_last_rather_than_crashing():
    """A failed design must not take a comparison down with it."""
    assert cost_per_coverage_point(set_cost(TWELVE), 0.0) == float("inf")


def test_prices_can_be_overridden_from_a_real_quote():
    """The defaults are order-of-magnitude. A user with a contract price should
    be able to use it."""
    model = OligoCostModel.from_overrides(per_pto_cost=0.0, setup_cost=None)

    assert model.per_pto_cost == 0.0
    assert model.setup_cost == OligoCostModel().setup_cost
