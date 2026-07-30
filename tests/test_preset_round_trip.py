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


@pytest.mark.parametrize("name", PRESET_NAMES)
def test_every_preset_resolves_and_is_complete(name):
    got = load_preset_conditions(name)
    assert got, f"preset {name!r} did not resolve"
    for field in PRESET_ADDITIVE_FIELDS:
        assert field in got, f"preset {name!r} dropped {field}"
    for key in ("reaction_temp", "polymerase", "na_conc", "mg_conc"):
        assert key in got


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
    from neoswga.core import reaction_conditions as rc

    preset_map = {
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
    conditions = preset_map[name]()
    for field in PRESET_ADDITIVE_FIELDS:
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
    for command in ("filter", "score", "optimize"):
        if command not in subparsers:
            continue
        action = next(
            (a for a in subparsers[command]._actions if a.dest == "preset"), None
        )
        if action is None or not action.choices:
            continue
        for choice in action.choices:
            assert load_preset_conditions(choice), (
                f"`{command} --preset {choice}` is offered by the CLI but does not "
                f"resolve to a preset"
            )
