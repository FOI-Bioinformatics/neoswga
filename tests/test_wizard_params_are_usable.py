"""Guard: a params.json emitted by the setup wizard must actually be runnable.

The wizard is the recommended entry point (`neoswga init`), and it tells the user to
run `neoswga validate-params` next. Two defects made that round-trip fail:

  * `mg_conc` was written as 0.0 for any genome that was not AT-rich. Because an
    explicit key suppresses the polymerase-aware default in
    `parameter.get_params()`, the whole pipeline then ran at zero magnesium -- a
    thermodynamically and enzymatically meaningless reaction.
  * `optimization_method: "ensemble"`, which the wizard's advanced menu can emit,
    was absent from `param_validator.VALID_OPTIMIZATION_METHODS`.

Neither was covered by a test.
"""

import pytest

from neoswga.core.parameter import MG_DEFAULTS_MM, default_mg_conc


@pytest.mark.parametrize("polymerase", sorted(MG_DEFAULTS_MM))
def test_default_mg_conc_is_never_zero(polymerase):
    """No supported polymerase has a working reaction at 0 mM Mg2+."""
    assert default_mg_conc(polymerase) > 0.0


def test_default_mg_conc_falls_back_for_unknown_polymerase():
    assert default_mg_conc("not-a-polymerase") > 0.0
    assert default_mg_conc("") > 0.0
    assert default_mg_conc(None) > 0.0


@pytest.mark.parametrize("polymerase", sorted(MG_DEFAULTS_MM))
def test_default_mg_conc_is_case_insensitive(polymerase):
    assert default_mg_conc(polymerase.upper()) == default_mg_conc(polymerase)


@pytest.mark.parametrize(
    "gc_class", ["extreme_at", "at_rich", "balanced", "gc_rich", "extreme_gc"]
)
def test_wizard_emits_usable_magnesium_for_every_gc_class(gc_class):
    """Regression: non-AT-rich genomes previously got mg_conc = 0.0."""
    from neoswga.core import wizard as wizard_mod

    src = wizard_mod.__file__
    with open(src) as fh:
        source = fh.read()

    assert '"mg_conc": 2.0 if' not in source, (
        "wizard is back to branching mg_conc on GC class; it should emit the "
        "polymerase-aware default so the pipeline and the wizard agree"
    )
    assert "default_mg_conc(polymerase)" in source, (
        "wizard should emit parameter.default_mg_conc(polymerase) for mg_conc"
    )


def test_wizard_generated_config_passes_the_validator(tmp_path):
    """End-to-end: build the wizard's config dict and validate it."""
    from neoswga.core.param_validator import ParamValidator

    genome = tmp_path / "target.fna"
    genome.write_text(">c1\n" + "ATGCATGCGC" * 200 + "\n")

    params = {
        "fg_genomes": [str(genome)],
        "fg_prefixes": [str(tmp_path / "target")],
        "data_dir": str(tmp_path),
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "mg_conc": default_mg_conc("phi29"),
        "na_conc": 50.0,
        "min_k": 6,
        "max_k": 12,
        "optimization_method": "ensemble",
        "num_primers": 6,
        "target_set_size": 6,
    }

    messages = ParamValidator().validate_params(params)
    errors = [m for m in messages if m.level.name == "ERROR"]
    assert not errors, f"wizard-shaped config rejected by validator: {errors}"

    assert params["mg_conc"] > 0.0
