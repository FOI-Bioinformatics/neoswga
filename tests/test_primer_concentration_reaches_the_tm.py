"""The configured primer concentration must reach the melting temperature.

Audit finding F1. `ReactionConditions.calculate_effective_tm` took
`primer_conc` as an argument defaulting to 0.5 uM, the constructor discarded the
configured value, and filtering, occupancy and effective coverage all called the
method without passing one. So `primer_conc` in params.json was accepted,
validated against a schema range, assigned a module global -- and never reached
the calculation it names.

The spread is not marginal. Across the configurable range the same oligo moves
about 10 C, and melting temperature decides which candidates survive the Tm gate
in `filter`, so a design run at a non-default concentration was filtered against
the wrong window.

**The figure is per oligo, not per pool.** Nearest-neighbour Tm depends on the
concentration of the annealing strand, so a 96-oligo panel at 0.5 uM each is the
supported experiment. A fixed total concentration shared across a growing panel
is a different experiment and is NOT modelled here: raising panel size would
lower per-oligo concentration and shift every Tm, and nothing in the pipeline
does that.
"""

import pytest

from neoswga.core.reaction_conditions import ReactionConditions

PRIMER = "ACGTACGTACGT"


def test_the_constructor_retains_the_configured_concentration():
    assert ReactionConditions(temp=30.0, polymerase="phi29", primer_conc=5e-6).primer_conc == 5e-6


def test_configuring_a_concentration_changes_the_effective_tm():
    """Probe from the audit: 0.05 uM against 5 uM returned an identical Tm."""
    low = ReactionConditions(temp=30.0, polymerase="phi29", primer_conc=0.05e-6)
    high = ReactionConditions(temp=30.0, polymerase="phi29", primer_conc=5e-6)

    assert low.calculate_effective_tm(PRIMER) != pytest.approx(high.calculate_effective_tm(PRIMER))
    # Higher concentration favours the duplex, so Tm rises.
    assert high.calculate_effective_tm(PRIMER) > low.calculate_effective_tm(PRIMER)


def test_an_explicit_argument_still_wins_over_the_configured_value():
    conditions = ReactionConditions(temp=30.0, polymerase="phi29", primer_conc=0.05e-6)

    assert conditions.calculate_effective_tm(PRIMER, primer_conc=5e-6) == pytest.approx(
        ReactionConditions(temp=30.0, polymerase="phi29", primer_conc=5e-6).calculate_effective_tm(
            PRIMER
        )
    )


def test_the_default_is_unchanged_when_nothing_is_configured():
    """Existing designs must not move because this gap was closed."""
    assert ReactionConditions(temp=30.0, polymerase="phi29").calculate_effective_tm(
        PRIMER
    ) == pytest.approx(
        ReactionConditions(
            temp=30.0, polymerase="phi29", primer_conc=0.5e-6
        ).calculate_effective_tm(PRIMER)
    )


def test_it_survives_the_builder_the_pipeline_actually_uses():
    """Tested on the real path, not only the low-level function.

    `filter`, occupancy and effective coverage all call
    `conditions.calculate_effective_tm(primer)` with no concentration argument,
    so the configured value has to arrive on the conditions object or it does
    not arrive at all. `build_reaction_conditions` is what assembles that object
    from configuration, and it is where the value was being dropped.
    """
    from neoswga.core.reaction_conditions import build_reaction_conditions

    def build(primer_conc):
        return build_reaction_conditions(
            None, polymerase="phi29", temp=30.0, primer_conc=primer_conc
        )

    # Compared builder against builder: the builder also applies the
    # polymerase's buffer defaults, so a bare `ReactionConditions` is not the
    # same reaction and would differ for reasons unrelated to concentration.
    low, high = build(0.05e-6), build(5e-6)

    assert (low.primer_conc, high.primer_conc) == (0.05e-6, 5e-6)
    assert high.calculate_effective_tm(PRIMER) > low.calculate_effective_tm(PRIMER)
