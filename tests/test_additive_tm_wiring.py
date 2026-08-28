"""Guard: an additive with a Tm model must actually reach both Tm code paths.

`AdditiveConcentrations.calculate_tm_correction` has two implementations: a
legacy fixed-coefficient sum, and an Arrhenius temperature-dependent path used
whenever a reaction temperature is supplied (which is the normal case, since
`ReactionConditions` always passes one).

The Arrhenius path built its additive dictionary from a hardcoded if-chain. That
list had already fallen behind, so it was possible for an additive to be defined
in `ADDITIVE_TM_PARAMS`, accepted by the dataclass, applied on the legacy path,
and still contribute exactly nothing on the path actually taken. It is now
derived from `TM_PARAM_FIELDS`; these tests keep that mapping honest.
"""

import pytest

from neoswga.core.additives import TM_PARAM_FIELDS, AdditiveConcentrations
from neoswga.core.mechanistic_params import ADDITIVE_TM_PARAMS
from neoswga.core.reaction_conditions import ReactionConditions


def test_every_modelled_additive_has_a_field_mapping():
    missing = set(ADDITIVE_TM_PARAMS) - set(TM_PARAM_FIELDS)
    assert not missing, (
        f"{sorted(missing)} have Tm parameters but no field mapping, so they "
        f"would contribute nothing on the Arrhenius path"
    )


def test_field_mapping_has_no_phantoms():
    from dataclasses import fields

    declared = {f.name for f in fields(AdditiveConcentrations)}
    phantom = {f for f in TM_PARAM_FIELDS.values() if f not in declared}
    assert not phantom, f"TM_PARAM_FIELDS names non-existent field(s): {sorted(phantom)}"

    unknown = set(TM_PARAM_FIELDS) - set(ADDITIVE_TM_PARAMS)
    assert not unknown, f"TM_PARAM_FIELDS maps additives with no Tm model: {sorted(unknown)}"


@pytest.mark.parametrize("key,field", sorted(TM_PARAM_FIELDS.items()))
def test_each_additive_shifts_tm_on_both_paths(key, field):
    """A modelled additive must move Tm whether or not a temperature is given."""
    conc = min(ADDITIVE_TM_PARAMS[key]["max_concentration"], 1.0)
    additives = AdditiveConcentrations(**{field: conc})

    legacy = additives.calculate_tm_correction(gc_content=0.5, primer_length=12)
    arrhenius = additives.calculate_tm_correction(
        gc_content=0.5, primer_length=12, reaction_temp_celsius=37.0
    )

    assert legacy != 0.0, f"{key} has no effect on the legacy path"
    assert arrhenius != 0.0, f"{key} has no effect on the Arrhenius path"
    # At the 37 C reference the two should broadly agree in sign and scale.
    assert legacy * arrhenius > 0, f"{key} disagrees in sign between the two paths"


def test_propanediol_matches_the_published_measurement():
    """Horakova et al. (2011) BMC Biotechnol 11:41: 4.9-5.9 C at 1 M, 45.5% GC.

    Checked at the 37 C reference temperature, which is where the coefficient
    was derived; the Arrhenius correction moves it at other temperatures.
    """
    conditions = ReactionConditions(temp=37.0, polymerase="klenow", propanediol_m=1.0)
    shift = conditions.calculate_tm_correction(gc_content=0.455, primer_length=12)
    assert -5.9 <= shift <= -4.9, f"expected the published range, got {shift:.2f} C"


def test_propanediol_is_dose_dependent():
    shifts = [
        ReactionConditions(temp=37.0, polymerase="klenow", propanediol_m=m).calculate_tm_correction(
            gc_content=0.455, primer_length=12
        )
        for m in (0.0, 0.5, 1.0, 1.5)
    ]
    assert shifts[0] == 0.0
    assert shifts[1] > shifts[2] > shifts[3], "Tm depression must grow with concentration"


def test_unsourced_gc_rich_additives_are_absent():
    """Ethylene glycol and sulfolane are deliberately NOT modelled.

    Both are real, validated GC-rich enhancers and both are recommended
    alongside 1,2-propanediol. Neither has a Tm coefficient we could source, and
    a guessed coefficient is worse than none: it acquires false authority the
    moment someone reads it out of the table. They are documented as
    real-but-unmodelled in the swga-chemistry skill instead.

    If a sourced value turns up, add it and delete this test.
    """
    assert "ethylene_glycol" not in ADDITIVE_TM_PARAMS
    assert "sulfolane" not in ADDITIVE_TM_PARAMS
