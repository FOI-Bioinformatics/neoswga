"""The two commands that recommend magnesium recommended a PCR figure.

`neoswga suggest` and `neoswga optimize-conditions` are where a user goes to be
told what buffer to run. Both carried magnesium grids from PCR practice:

    condition_suggester.MG_RECOMMENDATIONS   0.5 - 2.5 mM
    AdditiveOptimizer.SEARCH_RANGES["mg_conc"]   1.5 - 4.0 mM

Strand-displacing isothermal amplification is not PCR. The standard phi29 buffer
is 10 mM MgCl2 (Thermo MAN0030290, NEB phi29 buffer), which is what the registry
already supplies through `mg_default_mm` and what every shipped preset applies.
The mechanistic model was recalibrated to match -- `mg_optimal = 10.0`,
`mg_low_threshold = 4.0` -- when the PCR-era 2.5 mM optimum was removed along
with its mis-attributed citation (registry/INCONSISTENCIES.md #6).

The two grids were missed. Every value in both sat at or below the model's own
low threshold, so `suggest` recommended magnesium the rest of the tool scores as
deficient, and `optimize-conditions` could not search the band containing the
optimum however long it ran. A user following either would set up a reaction the
tool's own enzyme model then penalises.

The GC-dependent modulation is kept -- more magnesium for AT-rich template is a
reasonable heuristic and is documented as one in SCIENCE_CITATIONS.md -- but it
now modulates around the polymerase's buffer value rather than around a PCR
midpoint.
"""

import pytest

from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

MG_LOW = MECHANISTIC_MODEL_PARAMS["enzyme"]["mg_low_threshold"]
MG_HIGH = MECHANISTIC_MODEL_PARAMS["enzyme"]["mg_high_threshold"]
MG_OPTIMAL = MECHANISTIC_MODEL_PARAMS["enzyme"]["mg_optimal"]


# ----------------------------------------------------------------------
# suggest
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "gc_content",
    [0.19, 0.30, 0.41, 0.50, 0.65, 0.72],
    ids=["p_falciparum", "at_rich", "prevotella", "e_coli", "gc_rich", "extreme_gc"],
)
def test_suggested_magnesium_is_inside_the_models_usable_band(gc_content):
    """The failure at its simplest: whatever `suggest` prints, the enzyme model
    must not treat it as a deficiency."""
    from neoswga.core.condition_suggester import ConditionSuggester

    rec = ConditionSuggester(genome_gc=gc_content).suggest()

    assert MG_LOW <= rec.mg_conc <= MG_HIGH, (
        f"suggest recommends {rec.mg_conc} mM at GC {gc_content}, outside the "
        f"model's usable band [{MG_LOW}, {MG_HIGH}]"
    )


def test_the_gc_heuristic_still_runs_in_the_documented_direction():
    """More magnesium for AT-rich template, less for GC-rich. Anchoring the
    scale correctly should not flatten the heuristic into a constant."""
    from neoswga.core.condition_suggester import ConditionSuggester

    at_rich = ConditionSuggester(genome_gc=0.22).suggest().mg_conc
    balanced = ConditionSuggester(genome_gc=0.50).suggest().mg_conc
    gc_rich = ConditionSuggester(genome_gc=0.72).suggest().mg_conc

    assert at_rich > balanced > gc_rich


def test_a_balanced_genome_gets_the_buffer_value():
    """With no compositional reason to move, the recommendation should be the
    vendor buffer rather than an interpolation of its own."""
    from neoswga.core.condition_suggester import ConditionSuggester
    from neoswga.core.parameter import default_mg_conc

    rec = ConditionSuggester(genome_gc=0.50).suggest()

    assert rec.mg_conc == default_mg_conc(rec.polymerase)


# ----------------------------------------------------------------------
# optimize-conditions
# ----------------------------------------------------------------------


@pytest.mark.parametrize("grid", ["SEARCH_RANGES", "QUICK_RANGES"])
def test_the_search_grid_spans_the_usable_band(grid):
    """A grid entirely below `mg_low_threshold` cannot find the optimum, and
    reports whichever of its deficient values scored least badly."""
    from neoswga.core.additive_optimizer import AdditiveOptimizer

    values = getattr(AdditiveOptimizer, grid)["mg_conc"]

    assert min(values) >= MG_LOW, f"{grid} searches below the usable band: {values}"
    assert max(values) <= MG_HIGH, f"{grid} searches above the usable band: {values}"
    assert MG_OPTIMAL in values, f"{grid} cannot reach the optimum {MG_OPTIMAL}: {values}"


def test_the_optimizer_can_return_the_buffer_value():
    """End to end: the recommendation a user would actually act on."""
    from neoswga.core.additive_optimizer import AdditiveOptimizer

    rec = AdditiveOptimizer(polymerase="phi29").optimize(
        primer_length=12, template_gc=0.5, quick=True
    )

    assert MG_LOW <= rec.mg_conc <= MG_HIGH
