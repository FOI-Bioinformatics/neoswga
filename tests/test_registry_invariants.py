"""Invariants over the polymerase registry, plus the ledger of known violations.

Two jobs:

1. **Prove the seeding is faithful.** The registry replaced ~7 hardcoded tables. The
   `test_registry_reproduces_*` tests assert it reproduces each one exactly, which is
   what makes the refactor safe to land with the rest of the suite untouched.

2. **Track the contradictions it inherited.** The old tables disagreed with each other
   and, in two cases, with physics. Those values were seeded verbatim rather than
   quietly fixed, because changing them re-ranks primer sets. The violated invariants
   are marked `xfail(strict=True)`: fixing one turns the test green and *fails the
   suite* until this file and `INCONSISTENCIES.md` are updated to match.

See `neoswga/core/registry/INCONSISTENCIES.md` for the full write-up of each.
"""

import pytest

from neoswga.core import registry as R
from neoswga.core.registry import POLYMERASES


def _contains(outer, inner):
    """True if the closed interval `inner` lies within `outer`."""
    return outer[0] <= inner[0] and inner[1] <= outer[1]


# ----------------------------------------------------------------------
# Seeding fidelity: the registry must reproduce what it replaced
# ----------------------------------------------------------------------


def test_registry_reproduces_polymerase_characteristics():
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS

    assert R.as_characteristics() == POLYMERASE_CHARACTERISTICS


def test_registry_reproduces_validator_tables():
    from neoswga.core.param_validator import (
        POLYMERASE_PRIMER_LENGTHS,
        POLYMERASE_TEMP_RANGES,
    )

    assert R.warn_temp_ranges() == POLYMERASE_TEMP_RANGES
    assert R.primer_length_ranges() == POLYMERASE_PRIMER_LENGTHS


def test_registry_reproduces_mechanistic_enzyme_params():
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

    existing = {
        k: v for k, v in MECHANISTIC_MODEL_PARAMS["enzyme"].items() if isinstance(v, dict)
    }
    assert R.mechanistic_enzyme_params() == existing


def test_registry_reproduces_mg_defaults():
    from neoswga.core.parameter import MG_DEFAULTS_MM

    assert R.mg_defaults() == MG_DEFAULTS_MM


def test_registry_reproduces_additive_optimizer_constraints():
    from neoswga.core.additive_optimizer import AdditiveOptimizer

    assert R.additive_optimizer_constraints() == AdditiveOptimizer.POLYMERASE_CONSTRAINTS


def test_registry_reproduces_hybrid_presets():
    from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS

    for key, spec in POLYMERASES.items():
        cfg = POLYMERASE_PRESETS[key]
        assert cfg.max_extension == spec.legacy_hybrid_max_extension, key
        assert cfg.thermo_filter == spec.thermo_filter, key
        assert cfg.primer_multiplier == spec.primer_multiplier, key
        assert cfg.reaction_temp == spec.preset_reaction_temp, key
        assert (cfg.min_primer_tm, cfg.max_primer_tm) == spec.primer_tm_range, key
        assert (cfg.min_gc, cfg.max_gc) == spec.gc_range, key


# ----------------------------------------------------------------------
# Lookup behaviour
# ----------------------------------------------------------------------


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_get_polymerase_round_trips(key):
    assert R.get_polymerase(key).key == key
    assert R.get_polymerase(key.upper()).key == key


def test_get_polymerase_raises_on_unknown_by_default():
    with pytest.raises(ValueError, match="Unknown polymerase"):
        R.get_polymerase("nonesuch")


def test_get_polymerase_honours_explicit_default():
    assert R.get_polymerase("nonesuch", default="phi29").key == "phi29"


def test_polymerase_names_are_stable_and_ordered():
    assert R.polymerase_names() == list(POLYMERASES)
    assert {"phi29", "equiphi29", "bst", "klenow"} <= set(R.polymerase_names())


# ----------------------------------------------------------------------
# Invariants that hold today
# ----------------------------------------------------------------------


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_warn_band_contains_optimal_band(key):
    spec = POLYMERASES[key]
    assert _contains(spec.temp_warn_range, spec.temp_optimal_range) or _contains(
        spec.temp_hard_range, spec.temp_optimal_range
    ), f"{key}: vendor-optimal band sits outside both the warn and hard bands"


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_optimal_temp_sits_inside_the_hard_band(key):
    spec = POLYMERASES[key]
    lo, hi = spec.temp_hard_range
    assert lo <= spec.optimal_temp <= hi, f"{key}: optimal_temp outside the reject band"


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_coverage_reach_does_not_exceed_processivity(key):
    """Per-primer reach must fit inside the single-molecule maximum.

    Exempt for distributive enzymes: they dissociate after a short run and
    re-bind, so effective reach over a long incubation legitimately exceeds
    single-event processivity. Klenow (40 nt processivity, ~500 bp effective
    reach) is the case this exists for.
    """
    spec = POLYMERASES[key]
    if spec.distributive:
        pytest.skip(f"{key} is distributive; reach is not bounded by processivity")
    assert spec.typical_amplicon_bp <= spec.processivity_bp, key


def test_distributive_enzymes_are_flagged():
    """A very low processivity must be accompanied by the distributive flag."""
    for key, spec in POLYMERASES.items():
        if spec.processivity_bp < 1000:
            assert spec.distributive, (
                f"{key} has processivity {spec.processivity_bp} bp but is not "
                f"flagged distributive; reach/coverage reasoning will be wrong"
            )


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_default_primer_length_sits_in_the_design_range(key):
    spec = POLYMERASES[key]
    lo, hi = spec.primer_length_range
    assert lo <= spec.default_primer_length <= hi, key


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_every_polymerase_has_a_usable_magnesium_default(key):
    assert POLYMERASES[key].mg_default_mm > 0.0, key


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_dmso_response_is_complete(key):
    assert set(POLYMERASES[key].dmso_response) == {
        "threshold",
        "mild_coef",
        "steep_coef",
    }, key


@pytest.mark.parametrize("key", sorted(POLYMERASES))
def test_additive_limits_are_complete(key):
    assert set(POLYMERASES[key].additive_limits) >= {"dmso_percent", "betaine_m"}, key


# ----------------------------------------------------------------------
# Ledger: invariants VIOLATED by the seeded data
#
# strict=True on purpose. If one starts passing, the underlying value was fixed --
# which means INCONSISTENCIES.md and this marker must be updated in the same change.
# ----------------------------------------------------------------------


def _ledger(key, entry, detail):
    """Mark exactly the polymerase that violates a ledgered invariant."""
    return pytest.param(
        key,
        marks=pytest.mark.xfail(
            strict=True, reason=f"INCONSISTENCIES.md #{entry}: {detail}"
        ),
    )


@pytest.mark.parametrize(
    "key",
    [
        _ledger(
            "bst",
            1,
            "bst legacy_hybrid_max_extension is 10000 but bst processivity is 2000 - "
            "the set-cover stage credits bst primers with 5x the reach the enzyme "
            "has. Fixing re-ranks every bst hybrid-optimizer result.",
        ),
        _ledger(
            "klenow",
            1,
            "klenow legacy_hybrid_max_extension is 5000 against a corrected 40 nt "
            "processivity and ~500 bp effective reach - a 10x over-credit. This "
            "violation was HIDDEN by the old 10000 bp processivity figure and was "
            "exposed by correcting it.",
        ),
        "phi29",
        "equiphi29",
    ],
)
def test_set_cover_reach_does_not_exceed_processivity(key):
    spec = POLYMERASES[key]
    assert spec.legacy_hybrid_max_extension <= spec.processivity_bp, key


@pytest.mark.parametrize(
    "key",
    [
        _ledger(
            "klenow",
            2,
            "klenow temp_warn_range (20,42) is wider than temp_hard_range (25,40), "
            "so a klenow run at 41 C passes validate-params and then raises inside "
            "ReactionConditions.__init__.",
        ),
        "phi29",
        "equiphi29",
        "bst",
    ],
)
def test_warn_band_is_contained_in_hard_band(key):
    spec = POLYMERASES[key]
    assert _contains(spec.temp_hard_range, spec.temp_warn_range), (
        f"{key}: validate-params would accept a temperature that "
        f"ReactionConditions rejects"
    )


@pytest.mark.parametrize(
    "key",
    [
        _ledger(
            "bst",
            3,
            "bst preset_reaction_temp is 60.0 while its optimal_temp is 63.0 - "
            "hybrid_optimizer and additive_optimizer run bst 3 C below the optimum "
            "the rest of the codebase advertises.",
        ),
        "phi29",
        "equiphi29",
        "klenow",
    ],
)
def test_preset_reaction_temp_matches_optimal_temp(key):
    spec = POLYMERASES[key]
    assert spec.preset_reaction_temp == spec.optimal_temp, key


def test_klenow_processivity_matches_literature():
    """Fixed: was 10000 bp citing Bambara 1978, which does not support it."""
    assert POLYMERASES["klenow"].processivity_bp <= 100


def test_phi29_extension_rate_matches_literature():
    """Fixed: was 150 nt/s against Blanco 1989's ~53 nt/s."""
    assert POLYMERASES["phi29"].extension_rate_nt_s == pytest.approx(53, abs=10)


def test_magnesium_model_covers_the_default_concentration():
    """Fixed: the model was PCR-calibrated (2.5 mM) and penalised its own default."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS

    enzyme = MECHANISTIC_MODEL_PARAMS["enzyme"]
    assert POLYMERASES["phi29"].mg_default_mm <= enzyme["mg_high_threshold"]


# ----------------------------------------------------------------------
# Preset isolation
# ----------------------------------------------------------------------


def test_polymerase_config_is_not_shared_between_callers():
    """Regression: `_get_polymerase_config` handed out the shared preset object.

    `PolymeraseConfig` is a mutable dataclass and `HybridOptimizer._adjust_for_gc`
    writes `min_gc`/`max_gc` on it, so one genome's GC adaptation permanently
    mutated the module-level preset for every optimizer built later in the same
    process -- across ensemble methods, multi-genome targets and long sessions.
    """
    from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS, _get_polymerase_config

    original = (
        POLYMERASE_PRESETS["equiphi29"].min_gc,
        POLYMERASE_PRESETS["equiphi29"].max_gc,
    )

    cfg = _get_polymerase_config("equiphi29")
    assert cfg is not POLYMERASE_PRESETS["equiphi29"], "config must be a copy"

    cfg.min_gc, cfg.max_gc = 0.20, 0.80

    assert (
        POLYMERASE_PRESETS["equiphi29"].min_gc,
        POLYMERASE_PRESETS["equiphi29"].max_gc,
    ) == original, "mutating a returned config leaked into the shared preset"


def test_two_optimizers_do_not_share_gc_adaptation():
    """An AT-rich run must not widen the GC window for a later GC-rich run."""
    from neoswga.core.hybrid_optimizer import _get_polymerase_config

    a = _get_polymerase_config("phi29")
    b = _get_polymerase_config("phi29")
    assert a is not b

    a.min_gc = 0.20
    assert b.min_gc != 0.20, "two configs for the same polymerase are aliased"


def test_unknown_polymerase_falls_back_to_phi29_without_aliasing():
    from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS, _get_polymerase_config

    cfg = _get_polymerase_config("nonesuch")
    assert cfg.max_extension == POLYMERASE_PRESETS["phi29"].max_extension
    assert cfg is not POLYMERASE_PRESETS["phi29"]


# ----------------------------------------------------------------------
# Downstream surfaces must agree with the registry
# ----------------------------------------------------------------------


def test_schema_polymerase_enum_matches_registry():
    """params.schema.json must offer exactly the registry's polymerases."""
    import json
    from pathlib import Path

    schema_path = (
        Path(__file__).resolve().parents[1]
        / "neoswga"
        / "core"
        / "schema"
        / "params.schema.json"
    )
    schema = json.loads(schema_path.read_text())
    assert set(schema["properties"]["polymerase"]["enum"]) == set(POLYMERASES), (
        "params.schema.json polymerase enum has drifted from the registry"
    )


def test_cli_polymerase_choices_match_registry():
    """Every --polymerase flag must offer exactly the registry's polymerases."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    subparsers = parser._subparsers._group_actions[0].choices

    checked = 0
    for command, sub in subparsers.items():
        action = next(
            (a for a in sub._actions if a.dest == "polymerase" and a.choices), None
        )
        if action is None:
            continue
        checked += 1
        assert set(action.choices) == set(POLYMERASES), (
            f"`{command} --polymerase` choices have drifted from the registry"
        )
    assert checked > 0, "no --polymerase flags found; did the CLI change?"


def test_no_module_still_hardcodes_the_polymerase_list():
    """Catch a fifth copy of the polymerase list being reintroduced."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "neoswga"
    literal = '"phi29", "equiphi29", "bst", "klenow"'
    offenders = [
        str(p.relative_to(root))
        for p in root.rglob("*.py")
        # the registry itself is allowed to name them
        if "registry/" not in str(p.relative_to(root))
        and literal in p.read_text()
    ]
    assert not offenders, (
        f"polymerase list hardcoded again in: {offenders}. Use "
        f"neoswga.core.registry.polymerase_names()."
    )
