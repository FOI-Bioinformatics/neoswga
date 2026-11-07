"""Guard: `--preset <name>` must actually apply the preset it names.

`load_preset_conditions` flattens a `ReactionConditions` preset into a dict that the
CLI merges onto the global parameter module. It hand-listed the fields to copy and
had gone stale, omitting `urea_m`, `tmac_m`, `ethanol_percent` and
`formamide_percent`. The visible consequence: `--preset extreme_gc` dropped
`urea_m=0.5` and `tmac_m=0.05`, i.e. the two additives that distinguish that preset
from `high_gc_genome`, so the flag quietly did something other than advertised.

Also covers the `filter --preset` NameError: `run_step2` called
`load_preset_conditions` without importing it into `neoswga/cli/pipeline.py`, so the
flag raised `NameError` for every preset. No test passed `--preset` at the time.
"""

import inspect
import json

import pytest

from neoswga.cli._common import PRESET_ADDITIVE_FIELDS, load_preset_conditions

PRESET_NAMES = [
    "standard_phi29",
    "enhanced_equiphi29",
    "high_gc_genome",
    "long_primers_15mer",
    "q_solution",
    "gc_melt",
    "crude_sample",
    "low_temp",
    "bst",
    "klenow",
    "extreme_gc",
]


def _preset_funcs():
    """Preset name to the builder `load_preset_conditions` uses for it."""
    from neoswga.core import reaction_conditions as rc

    return {
        "standard_phi29": rc.get_standard_conditions,
        "enhanced_equiphi29": rc.get_enhanced_conditions,
        "high_gc_genome": rc.get_high_gc_conditions,
        "long_primers_15mer": rc.get_enhanced_conditions,
        "q_solution": rc.get_q_solution_equivalent,
        "gc_melt": rc.get_gc_melt_conditions,
        "crude_sample": rc.get_crude_sample_conditions,
        "low_temp": rc.get_low_temp_conditions,
        "bst": rc.get_bst_conditions,
        "klenow": rc.get_klenow_conditions,
        "extreme_gc": rc.get_extreme_gc_conditions,
    }


_PRESET_FUNCS = _preset_funcs()


def _reaction_conditions_additive_params():
    """Concentration knobs accepted by ReactionConditions, excluding salts/flags."""
    from neoswga.core.reaction_conditions import ReactionConditions

    sig = inspect.signature(ReactionConditions.__init__)
    skip = {"self", "temp", "polymerase", "na_conc", "mg_conc", "ssb", "primer_conc"}
    return {
        name
        for name, p in sig.parameters.items()
        if name not in skip
        and p.kind is not inspect.Parameter.VAR_KEYWORD
        and isinstance(p.default, (int, float))
        and not isinstance(p.default, bool)
    }


def test_preset_additive_fields_cover_reaction_conditions():
    """Every concentration knob must be carried through a preset."""
    expected = _reaction_conditions_additive_params()
    missing = expected - set(PRESET_ADDITIVE_FIELDS)
    assert not missing, (
        f"PRESET_ADDITIVE_FIELDS is missing {sorted(missing)}; a preset setting any "
        f"of these would silently lose it on the way to the pipeline"
    )


def test_preset_additive_fields_have_no_phantoms():
    from neoswga.core.reaction_conditions import ReactionConditions

    sig = inspect.signature(ReactionConditions.__init__)
    phantom = set(PRESET_ADDITIVE_FIELDS) - set(sig.parameters)
    assert not phantom, f"PRESET_ADDITIVE_FIELDS names non-existent field(s): {sorted(phantom)}"


def _is_dropped(name, field):
    """Whether this preset's value for `field` is one a design refuses.

    Such a field is deliberately left OUT of the flattened dict since
    2026-11-01, so the two assertions below exempt exactly it and nothing
    else -- a preset that dropped `urea_m` still fails.
    """
    from neoswga.core.model_evidence import refuses_additive

    conditions = _PRESET_FUNCS[name]()
    return bool(getattr(conditions, field, 0.0)) and refuses_additive(field)


@pytest.mark.parametrize("name", PRESET_NAMES)
def test_every_preset_resolves_and_is_complete(name):
    got = load_preset_conditions(name)
    assert got, f"preset {name!r} did not resolve"
    for field in PRESET_ADDITIVE_FIELDS:
        if _is_dropped(name, field):
            continue
        assert field in got, f"preset {name!r} dropped {field}"
    for key in ("reaction_temp", "polymerase", "na_conc", "mg_conc"):
        assert key in got


@pytest.mark.parametrize("name", PRESET_NAMES)
def test_an_additive_a_design_refuses_is_absent_rather_than_zeroed(name):
    """Left out, not carried as 0.0.

    `resolve_design_request` reads the params FILE and a preset reaches the
    `parameter` globals, so a preset's additives never passed that refusal:
    `crude_sample` and `q_solution` set `glycerol_percent`, whose duplex model
    the registry records as `absent`, and the value was accepted,
    range-checked, written into the run manifest's `effective_conditions` and
    applied to no Tm. Carrying it as 0.0 would be a different claim -- "this
    reaction has no glycerol" rather than "this preset does not speak about
    glycerol" -- and would leave an explicit value in params.json unable to
    override it.
    """
    got = load_preset_conditions(name)

    for field in PRESET_ADDITIVE_FIELDS:
        if _is_dropped(name, field):
            assert field not in got, f"{name!r} should not carry {field} at all"


def test_at_least_one_preset_exercises_the_drop():
    """Otherwise the two tests above are exempting nothing and prove nothing."""
    dropped = [
        (name, field)
        for name in PRESET_NAMES
        for field in PRESET_ADDITIVE_FIELDS
        if _is_dropped(name, field)
    ]

    assert dropped, (
        "no preset sets an additive a design refuses, so the drop path is "
        "untested; if the registry or the presets changed, re-point this"
    )


def test_extreme_gc_preserves_its_defining_additives():
    """Regression: urea_m and tmac_m were dropped, making the preset a no-op."""
    from neoswga.core.reaction_conditions import get_extreme_gc_conditions

    conditions = get_extreme_gc_conditions()
    got = load_preset_conditions("extreme_gc")

    assert conditions.urea_m > 0, "fixture assumption: extreme_gc sets urea"
    assert conditions.tmac_m > 0, "fixture assumption: extreme_gc sets TMAC"
    assert got["urea_m"] == pytest.approx(conditions.urea_m)
    assert got["tmac_m"] == pytest.approx(conditions.tmac_m)


@pytest.mark.parametrize("name", PRESET_NAMES)
def test_preset_values_match_the_conditions_object(name):
    """The flattened dict must not distort the values it copies."""
    from neoswga.cli._common import load_preset_conditions as load

    got = load(name)
    conditions = _PRESET_FUNCS[name]()
    for field in PRESET_ADDITIVE_FIELDS:
        if _is_dropped(name, field):
            continue
        assert got[field] == pytest.approx(getattr(conditions, field)), field
    assert got["reaction_temp"] == pytest.approx(conditions.temp)
    assert got["polymerase"] == conditions.polymerase


def test_unknown_preset_returns_empty_rather_than_raising():
    assert load_preset_conditions("no-such-preset") == {}


def test_filter_step_can_resolve_load_preset_conditions():
    """Regression: `neoswga filter --preset X` raised NameError."""
    from neoswga.cli import pipeline as pipeline_cli

    assert "load_preset_conditions(args.preset)" in inspect.getsource(pipeline_cli.run_step2)
    assert hasattr(pipeline_cli, "load_preset_conditions"), (
        "run_step2 calls load_preset_conditions but it is not importable from "
        "neoswga.cli.pipeline - `filter --preset` would raise NameError"
    )


def test_cli_preset_choices_all_resolve():
    """Every --preset value argparse offers must map to a real preset."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    subparsers = parser._subparsers._group_actions[0].choices
    for command in ("filter", "prepare-candidates", "optimize"):
        if command not in subparsers:
            continue
        action = next((a for a in subparsers[command]._actions if a.dest == "preset"), None)
        if action is None or not action.choices:
            continue
        for choice in action.choices:
            assert load_preset_conditions(choice), (
                f"`{command} --preset {choice}` is offered by the CLI but does not "
                f"resolve to a preset"
            )


# ----------------------------------------------------------------------
# A preset is part of the configuration, not a late overwrite (2026-11-01)
# ----------------------------------------------------------------------
#
# A preset used to reach the `parameter` globals by `setattr` and nothing
# else, so every guard that asks "did the user say this" by testing key
# presence in `_json_data` was blind to it. Two consequences, both measured in
# `docs/validation/2026-11-01-the-preset-path-bypassed-every-guard.md`.


class _Args:
    """Every CLI option unset, so only the file and the overlay decide.

    `__getattr__` rather than a listed namespace: `get_params` reads a few
    dozen attributes and a list of them would rot into a confusing
    `AttributeError` the first time it read one more.
    """

    def __init__(self, json_file):
        self.json_file = str(json_file)

    def __getattr__(self, name):
        return None


def _args(json_file):
    return _Args(json_file)


@pytest.fixture
def params_file(tmp_path):
    """A params.json that names no additive, so GC-adaptive is free to."""
    seq = "G" * 700 + "A" * 300
    fasta = tmp_path / "fg.fna"
    fasta.write_text(">seq\n" + seq + "\n")
    path = tmp_path / "params.json"
    path.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "src_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "fg")],
                "fg_seq_lengths": [len(seq)],
                "cpus": 1,
            }
        )
    )
    return path


def test_a_preset_lands_in_the_raw_configuration(params_file):
    """`_json_data` is where the key-presence guards look."""
    from neoswga.core import parameter

    parameter.reset_to_defaults()
    parameter._preset_overlay = load_preset_conditions("low_temp")
    try:
        parameter.get_params(_args(params_file))
        assert parameter._json_data.get("betaine_m") == pytest.approx(1.0)
        assert parameter._json_data.get("dmso_percent") == pytest.approx(3.0)
    finally:
        parameter._preset_overlay = {}
        parameter.reset_to_defaults()


def test_without_a_preset_the_raw_configuration_is_only_the_file(params_file):
    """The overlay must not invent keys the user never supplied, or
    GC-adaptive would read its own silence as a user choice."""
    from neoswga.core import parameter

    parameter.reset_to_defaults()
    parameter._preset_overlay = {}
    try:
        parameter.get_params(_args(params_file))
        assert "betaine_m" not in parameter._json_data
        assert "dmso_percent" not in parameter._json_data
    finally:
        parameter.reset_to_defaults()


def test_a_presets_explicit_zero_is_a_choice_not_silence(params_file):
    """The defect: `standard_phi29` writes betaine 0.0, GC-adaptive could not
    tell that from an absent key, logged "Setting betaine to 1.0M", and the
    preset then overwrote it. The log said one thing and the run used
    another."""
    from neoswga.core import parameter

    parameter.reset_to_defaults()
    parameter._preset_overlay = load_preset_conditions("standard_phi29")
    try:
        parameter.get_params(_args(params_file))
        assert "betaine_m" in parameter._json_data
        assert parameter._json_data["betaine_m"] == pytest.approx(0.0)
    finally:
        parameter._preset_overlay = {}
        parameter.reset_to_defaults()


def test_one_runs_preset_cannot_leak_into_the_next(params_file):
    """`_preset_overlay` is module state, so the CLI assigns it every run."""
    import inspect as _inspect

    from neoswga.cli import pipeline as pipeline_cli

    source = _inspect.getsource(pipeline_cli.run_step2)
    assert "parameter._preset_overlay = load_preset_conditions(_preset) if _preset else {}" in (
        source
    ), "run_step2 must assign the overlay unconditionally, not only when a preset is given"


# ----------------------------------------------------------------------
# The effective configuration meets the guards (2026-11-01)
# ----------------------------------------------------------------------
#
# `resolve_design_request` validates the params FILE, so none of the four
# guards in `model_evidence.require_model_support` ever saw a preset's
# chemistry. The live one is oligo length: a preset sets the POLYMERASE while
# `min_k`/`max_k` come from the file, and `--preset bst` computed Tm for
# 7-12mers at 63 C from a parameter set covering 15-25 nt. That is
# `_require_length_support`'s own stated defect.

#: preset -> the enzyme's supported oligo range, for presets whose enzyme does
#: not cover the documented 6-12 window.
OUTSIDE_THE_DEFAULT_WINDOW = {
    "enhanced_equiphi29": "10-18",
    "high_gc_genome": "10-18",
    "long_primers_15mer": "10-18",
    "q_solution": "10-18",
    "gc_melt": "10-18",
    "extreme_gc": "10-18",
    "bst": "15-25",
    "klenow": "8-15",
}

#: The rest, whose enzyme is phi29 and so covers 6-12 exactly.
INSIDE_THE_DEFAULT_WINDOW = ("standard_phi29", "crude_sample", "low_temp")


def _effective(preset, min_k, max_k):
    """A `parameter`-shaped stand-in carrying the resolved configuration."""
    overlay = load_preset_conditions(preset)
    param = _Args(None)
    param._json_data = {
        "fg_prefixes": ["t"],
        "fg_seq_lengths": [1_270_000],
        "min_k": min_k,
        "max_k": max_k,
    }
    param.min_k, param.max_k = min_k, max_k
    for key, value in overlay.items():
        setattr(param, key, value)
    return param


def test_every_preset_is_classified():
    assert set(OUTSIDE_THE_DEFAULT_WINDOW) | set(INSIDE_THE_DEFAULT_WINDOW) == set(PRESET_NAMES)


@pytest.mark.parametrize("preset,domain", sorted(OUTSIDE_THE_DEFAULT_WINDOW.items()))
def test_a_preset_whose_enzyme_excludes_the_window_is_refused(preset, domain):
    """Accepted silently before: the guard reads the file, and the preset
    never reached it."""
    from neoswga.cli._common import validate_effective_chemistry
    from neoswga.core.exceptions import DesignError

    with pytest.raises(DesignError) as caught:
        validate_effective_chemistry(_effective(preset, 6, 12), preset)

    message = str(caught.value)
    assert "oligo length" in message
    assert domain in message, f"the message should name the enzyme's {domain} nt range"


@pytest.mark.parametrize("preset", INSIDE_THE_DEFAULT_WINDOW)
def test_a_phi29_preset_still_passes_the_documented_window(preset):
    """phi29 covers 6-12, so these must not start failing."""
    from neoswga.cli._common import validate_effective_chemistry

    validate_effective_chemistry(_effective(preset, 6, 12), preset)


@pytest.mark.parametrize(
    "preset,min_k,max_k", (("bst", 15, 25), ("enhanced_equiphi29", 10, 18), ("klenow", 8, 15))
)
def test_a_matching_window_is_accepted(preset, min_k, max_k):
    """The remedy the refusal names has to work, or the message is a dead end."""
    from neoswga.cli._common import validate_effective_chemistry

    validate_effective_chemistry(_effective(preset, min_k, max_k), preset)


def test_the_validation_runs_only_when_a_preset_was_applied():
    """Without a preset the effective configuration IS the file, which
    `optimize` already resolves, so no run that works today starts failing."""
    source = inspect.getsource(__import__("neoswga.cli.pipeline", fromlist=["run_step2"]).run_step2)

    assert "if _preset:\n            validate_effective_chemistry(parameter, _preset)" in source


def test_the_glycerol_presets_are_refused_for_length_not_for_glycerol():
    """`q_solution` sets glycerol AND an enzyme outside the window. The
    glycerol half was fixed by dropping it from the preset, so what remains
    must be the length refusal -- otherwise this test would pass on the old
    bug."""
    from neoswga.cli._common import validate_effective_chemistry
    from neoswga.core.exceptions import DesignError

    with pytest.raises(DesignError) as caught:
        validate_effective_chemistry(_effective("q_solution", 6, 12), "q_solution")

    assert "glycerol" not in str(caught.value)
    assert "oligo length" in str(caught.value)
