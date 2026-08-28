"""Guard: every shipped reaction preset must describe a reaction that can run.

Two defects, both found while writing user-facing guidance rather than by any
test:

1. **Five of ten presets shipped with `mg_conc = 0.0`** - `high_gc_genome`,
   `q_solution`, `crude_sample`, `low_temp` and `extreme_gc`. They never set the
   field, and `ReactionConditions.__init__` defaulted it to zero. No polymerase
   reaction proceeds without magnesium, so those presets described nothing
   runnable. Two more used 2.5 mM, a PCR figure rather than the ~10 mM of the
   phi29-family buffer.

2. **`show-presets` printed values that `--preset` did not apply.** The listing
   read a hand-maintained dict in the CLI layer while `--preset` resolved through
   the `reaction_conditions` factories. They disagreed for three of four named
   presets - `long_primers_15mer` advertised 45 C / 7% DMSO / 1.5 M betaine and
   applied 42 C / 5% / 1.0 M.

The root fix for (1) was making `mg_conc=None` mean "this polymerase's buffer
value" rather than defaulting to zero; for (2), deriving the display dict from
the same factories `--preset` uses.
"""

import pytest

from neoswga.cli._common import PRESETS, load_preset_conditions
from neoswga.core import reaction_conditions as rc
from neoswga.core.registry import POLYMERASES

ALL_PRESETS = {
    "standard_phi29": rc.get_standard_conditions,
    "enhanced_equiphi29": rc.get_enhanced_conditions,
    "high_gc_genome": rc.get_high_gc_conditions,
    "q_solution": rc.get_q_solution_equivalent,
    "gc_melt": rc.get_gc_melt_conditions,
    "crude_sample": rc.get_crude_sample_conditions,
    "low_temp": rc.get_low_temp_conditions,
    "bst": rc.get_bst_conditions,
    "klenow": rc.get_klenow_conditions,
    "extreme_gc": rc.get_extreme_gc_conditions,
}


@pytest.mark.parametrize("name", sorted(ALL_PRESETS))
def test_preset_has_usable_magnesium(name):
    """A preset with no magnesium is not a reaction."""
    conditions = ALL_PRESETS[name]()
    assert conditions.mg_conc > 1.0, (
        f"preset {name!r} has mg_conc={conditions.mg_conc} mM - a polymerase "
        f"reaction does not proceed at that concentration"
    )


@pytest.mark.parametrize("name", sorted(ALL_PRESETS))
def test_preset_magnesium_is_in_the_models_usable_band(name):
    """The enzyme-activity model must not penalise a shipped preset."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

    enzyme = MECHANISTIC_MODEL_PARAMS["enzyme"]
    conditions = ALL_PRESETS[name]()
    assert enzyme["mg_low_threshold"] <= conditions.mg_conc <= enzyme["mg_high_threshold"], (
        f"preset {name!r} at {conditions.mg_conc} mM Mg2+ falls outside the "
        f"model's usable band - the tool would penalise its own preset"
    )


@pytest.mark.parametrize("polymerase", sorted(POLYMERASES))
def test_reaction_conditions_defaults_to_a_real_buffer(polymerase):
    """Constructing without mg_conc must give the polymerase's buffer value."""
    spec = POLYMERASES[polymerase]
    conditions = rc.ReactionConditions(temp=spec.optimal_temp, polymerase=polymerase)
    assert conditions.mg_conc == spec.mg_default_mm
    assert conditions.mg_conc > 0.0


def test_explicit_zero_magnesium_is_still_honoured():
    """The default is a default, not a floor - callers can still force 0."""
    conditions = rc.ReactionConditions(temp=30.0, polymerase="phi29", mg_conc=0.0)
    assert conditions.mg_conc == 0.0


@pytest.mark.parametrize("name", sorted(PRESETS))
def test_displayed_preset_matches_applied_preset(name):
    """`show-presets` must not advertise conditions `--preset` will not apply."""
    shown = PRESETS[name]
    applied = load_preset_conditions(name)

    for shown_key, applied_key in [
        ("temperature", "reaction_temp"),
        ("polymerase", "polymerase"),
        ("dmso_percent", "dmso_percent"),
        ("betaine_m", "betaine_m"),
        ("na_conc", "na_conc"),
        ("mg_conc", "mg_conc"),
    ]:
        assert shown[shown_key] == applied[applied_key], (
            f"preset {name!r}: show-presets reports {shown_key}="
            f"{shown[shown_key]} but --preset applies {applied[applied_key]}"
        )


def test_preset_display_dict_is_derived_not_handwritten():
    """Guard against reintroducing a parallel hand-maintained table."""
    import inspect

    from neoswga.cli import _common

    src = inspect.getsource(_common)
    assert "_build_presets" in src, (
        "PRESETS should be derived from load_preset_conditions; a literal dict "
        "will drift from what --preset applies, as it did before"
    )
