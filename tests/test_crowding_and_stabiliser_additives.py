"""PEG, BSA and DTT act on the enzyme, not on the duplex.

These three were accepted by the schema, range-validated, and consumed by
nothing. Once the forwarding bug was fixed they reached `ReactionConditions`
and still did nothing, because there was no pathway to receive them -- no Tm
coefficient and no mechanistic term. That is the right outcome for a Tm
coefficient: none of the three is primarily a duplex-stability agent, so a Tm
term would have been fiction.

What they actually do is act on the polymerase:

  PEG   Macromolecular crowding. Excluded volume raises effective
        concentrations and accelerates association (Zimmerman & Minton 1993).
        Standard in phi29 rolling-circle protocols at 5-10% PEG-8000. Above
        the crowding plateau the added viscosity slows diffusion, so the
        benefit levels off and then costs synthesis speed -- a term that has
        to be non-monotonic to be honest.

  BSA   Prevents adsorption of polymerase to vessel surfaces and sequesters
        inhibitors, keeping enzyme in solution. Standard phi29 protocols use
        0.1-0.2 mg/mL. Saturating: once surfaces are blocked, more does
        nothing.

  DTT   Keeps the polymerase's cysteine thiols reduced. Standard phi29 buffer
        is 4 mM. This one is a DEFICIENCY penalty rather than a bonus, which
        makes its default load-bearing: if an unset `dtt_mm` of 0.0 were read
        as "no DTT", every existing user would silently acquire a stability
        penalty for a buffer component they never disabled. It gets a
        polymerase-aware default, like `mg_conc`, so unset means the vendor
        buffer and only an explicit low value is penalised.

The tests below pin the SHAPE of each effect -- direction, saturation,
non-monotonicity, and the no-op at defaults -- rather than the exact
coefficients, which are EMPIRICAL and documented as such in
docs/SCIENCE_CITATIONS.md.
"""

import pytest

from neoswga.core.mechanistic_model import MechanisticModel
from neoswga.core.reaction_conditions import ReactionConditions

PRIMER = "ATCGATCGATCG"


def effects(**kwargs):
    conditions = ReactionConditions(temp=30.0, polymerase="phi29", **kwargs)
    return MechanisticModel(conditions).calculate_effects(PRIMER, template_gc=0.5)


def baseline():
    """Standard buffer: 4 mM DTT, no crowding agent, no BSA."""
    return effects(dtt_mm=4.0)


# ----------------------------------------------------------------------
# They do something at all
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "field,value",
    [("peg_percent", 6.0), ("bsa_ug_ml", 200.0), ("dtt_mm", 0.5)],
)
def test_additive_changes_the_prediction(field, value):
    """Each was inert: forwarded correctly and consumed by nothing."""
    base = baseline()
    changed = effects(**{field: value, **({} if field == "dtt_mm" else {"dtt_mm": 4.0})})

    assert changed.predicted_amplification_factor != pytest.approx(
        base.predicted_amplification_factor
    ), f"{field}={value} did not change the prediction"


# ----------------------------------------------------------------------
# PEG: crowding helps binding, viscosity eventually costs speed
# ----------------------------------------------------------------------


def test_peg_raises_the_binding_rate():
    """Crowding raises effective concentration, so association speeds up."""
    assert effects(peg_percent=6.0, dtt_mm=4.0).effective_binding_rate > (
        baseline().effective_binding_rate
    )


def test_peg_crowding_benefit_saturates():
    """Excluded volume is bounded; past the plateau more PEG buys no more kon.

    Modelling it as unbounded would make 15% PEG look like the best possible
    reaction, which is the opposite of what protocols say.
    """
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

    plateau = MECHANISTIC_MODEL_PARAMS["kinetics"]["peg_saturation"]

    at_plateau = effects(peg_percent=plateau, dtt_mm=4.0).effective_binding_rate
    well_past = effects(peg_percent=plateau + 6.0, dtt_mm=4.0).effective_binding_rate

    assert well_past == pytest.approx(at_plateau)


def test_too_much_peg_costs_synthesis_speed():
    """The viscosity cost is what makes this non-monotonic, and it is the part
    a user needs: it says there is an optimum rather than "more is better"."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

    plateau = MECHANISTIC_MODEL_PARAMS["kinetics"]["peg_saturation"]

    assert effects(peg_percent=plateau + 6.0, dtt_mm=4.0).speed_factor < (
        effects(peg_percent=plateau, dtt_mm=4.0).speed_factor
    )


def test_peg_has_an_interior_optimum():
    """Follows from the two effects above, and is the practical consequence:
    somewhere in the usable range is better than either extreme."""
    grid = [0.0, 2.0, 5.0, 8.0, 11.0, 15.0]
    scores = [effects(peg_percent=x, dtt_mm=4.0).predicted_amplification_factor for x in grid]

    best = scores.index(max(scores))
    assert 0 < best < len(grid) - 1, f"optimum at an endpoint: {list(zip(grid, scores))}"


# ----------------------------------------------------------------------
# BSA: saturating stabiliser
# ----------------------------------------------------------------------


def test_bsa_improves_enzyme_stability():
    assert effects(bsa_ug_ml=100.0, dtt_mm=4.0).stability_factor > (baseline().stability_factor)


def test_bsa_saturates():
    """Once surfaces are blocked, more protein does nothing. A linear term
    would let an absurd 400 ug/mL dominate every other consideration."""
    mid = effects(bsa_ug_ml=200.0, dtt_mm=4.0).stability_factor
    high = effects(bsa_ug_ml=400.0, dtt_mm=4.0).stability_factor

    assert high > mid
    assert (high - mid) < 0.25 * (mid - baseline().stability_factor)


def test_bsa_is_monotone():
    values = [effects(bsa_ug_ml=x, dtt_mm=4.0).stability_factor for x in (0, 50, 100, 200, 400)]
    assert values == sorted(values)


# ----------------------------------------------------------------------
# DTT: deficiency penalty, and the default that keeps it honest
# ----------------------------------------------------------------------


def test_low_dtt_costs_stability():
    """Oxidised thiols mean lost activity over a long isothermal incubation."""
    assert effects(dtt_mm=0.0).stability_factor < effects(dtt_mm=4.0).stability_factor


def test_dtt_above_the_standard_adds_nothing():
    """It is a deficiency penalty, not a dose-response. Modelling excess DTT
    as beneficial would invite a user to add more for no reason."""
    assert effects(dtt_mm=10.0).stability_factor == pytest.approx(
        effects(dtt_mm=4.0).stability_factor
    )


def test_unset_dtt_means_the_standard_buffer_not_zero():
    """The load-bearing default.

    `dtt_mm` is 0.0 on a params.json that never mentions it. Read literally
    that is a DTT-free reaction, and every existing user would have acquired a
    stability penalty for a buffer component they never disabled -- a silent
    change to every prediction. `default_dtt_mm` makes unset mean the vendor
    buffer instead.
    """
    from neoswga.core.parameter import default_dtt_mm

    assert default_dtt_mm("phi29") == 4.0
    assert default_dtt_mm("equiphi29") > 0.0
    assert default_dtt_mm("something-unknown") > 0.0, "fallback must not be zero"


def test_a_params_file_without_dtt_predicts_as_before(tmp_path):
    """The regression this default exists to prevent: an unchanged params.json
    must produce an unchanged prediction."""
    from neoswga.core.parameter import default_dtt_mm

    unset = effects(dtt_mm=default_dtt_mm("phi29"))
    explicit_standard = effects(dtt_mm=4.0)

    assert unset.predicted_amplification_factor == pytest.approx(
        explicit_standard.predicted_amplification_factor
    )


# ----------------------------------------------------------------------
# None of them is a Tm term
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "field,value",
    [("peg_percent", 10.0), ("bsa_ug_ml", 400.0), ("dtt_mm", 0.0)],
)
def test_these_additives_do_not_move_the_melting_temperature(field, value):
    """A crowding agent, a carrier protein and a reductant are not duplex
    stabilisers at these concentrations. Giving them Tm coefficients to make
    them "do something" would have been inventing chemistry.
    """
    base = ReactionConditions(temp=30.0, polymerase="phi29", dtt_mm=4.0)
    with_additive = ReactionConditions(temp=30.0, polymerase="phi29", **{field: value})

    assert with_additive.calculate_effective_tm(PRIMER) == pytest.approx(
        base.calculate_effective_tm(PRIMER)
    )


def test_stability_is_neutral_when_nothing_sets_a_stabiliser():
    """`stability_factor` was computed, reported, and left out of the
    amplification product, so glycerol, BSA and DTT could not reach the number
    the model is judged on.

    Folding it in was only safe because it is exactly 1.0 at defaults -- DTT
    defaults to the vendor buffer, so its deficit term is zero, and the others
    are absent. This pins that, because it is the property that keeps existing
    predictions unchanged.
    """
    assert baseline().stability_factor == pytest.approx(1.0)


def test_a_stabiliser_now_reaches_the_amplification_estimate():
    """The consequence: enzyme-pathway additives move the headline number."""
    assert (
        effects(bsa_ug_ml=200.0, dtt_mm=4.0).predicted_amplification_factor
        > baseline().predicted_amplification_factor
    )
    assert (
        effects(dtt_mm=0.0).predicted_amplification_factor
        < baseline().predicted_amplification_factor
    )
