"""The Tm window has to match the enzyme it is filtering for.

A primer is usable when it binds at the reaction temperature. Those
temperatures differ by more than thirty degrees across the supported enzymes --
phi29 runs at 30 C, bst at 63 C -- so a Tm window that does not move with the
polymerase is filtering for the wrong chemistry.

`filter.py` resolved the window as:

    tm_min = getattr(parameter, "min_tm", None) or 15
    tm_max = getattr(parameter, "max_tm", None) or 55

Two faults in two lines. The fallback is a fixed 15-55 whatever the enzyme, so
a params.json that does not mention Tm filtered a bst design at 63 C through a
window built for phi29 -- keeping primers that cannot prime at that
temperature and rejecting ones that can. And `or` treats a configured 0.0 as
absent, the same sentinel-versus-value confusion this audit has found
repeatedly.

The registry already carries `primer_tm_range` per enzyme, which is where the
optimizer's own presets come from. Defaulting from it makes the window
polymerase-aware and removes a fourth independent copy of the same concept:
params.json, this fallback, and the two optimizer presets were all separately
answering "what Tm is acceptable".
"""

import json

import pytest

from neoswga.core import parameter
from neoswga.core.registry import POLYMERASES

ENZYMES = ["phi29", "equiphi29", "bst", "bst3.0", "bsu", "klenow"]


# ----------------------------------------------------------------------
# The default tracks the enzyme
# ----------------------------------------------------------------------


@pytest.mark.parametrize("polymerase", ENZYMES)
def test_default_window_matches_the_registry(polymerase):
    """One source of truth, not a hand-copied table."""
    low, high = parameter.default_tm_range(polymerase)
    assert (low, high) == POLYMERASES[polymerase].primer_tm_range


def test_unknown_polymerase_still_yields_a_usable_window():
    low, high = parameter.default_tm_range("something-unsupported")
    assert low < high
    assert high > 0


def test_the_window_actually_differs_between_enzymes():
    """Guards the point of the change. If every enzyme resolved to the same
    window the parametrised tests above would pass while filtering stayed
    polymerase-blind.
    """
    phi29 = parameter.default_tm_range("phi29")
    bst = parameter.default_tm_range("bst")

    assert bst[0] > phi29[1] - 15, f"bst window {bst} barely differs from phi29 {phi29}"
    assert bst != phi29


@pytest.mark.parametrize("polymerase", ENZYMES)
def test_the_window_brackets_the_reaction_temperature(polymerase):
    """A primer has to be stable at the temperature the reaction runs at.

    The upper bound may sit above it -- that is the selectivity trade -- but a
    window entirely below the reaction temperature would reject every primer
    that can work.
    """
    low, high = parameter.default_tm_range(polymerase)
    reaction_temp = parameter.default_reaction_temp(polymerase)

    assert high > reaction_temp, (
        f"{polymerase}: whole Tm window {(low, high)} sits below its "
        f"{reaction_temp} C reaction temperature"
    )


# ----------------------------------------------------------------------
# params.json still wins where it speaks
# ----------------------------------------------------------------------


def _load(tmp_path, **overrides):
    import neoswga.core.pipeline as pipeline_mod

    genome = tmp_path / "g.fasta"
    genome.write_text(">g\n" + "ACGT" * 400 + "\n")
    params = {
        "schema_version": 2,
        "data_dir": str(tmp_path),
        "fg_genomes": [str(genome)],
        "fg_prefixes": [str(tmp_path / "g")],
        "fg_seq_lengths": [1600],
    }
    params.update(overrides)
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))

    parameter.json_file = str(path)
    pipeline_mod._initialized = False
    pipeline_mod._initialize()


@pytest.fixture(autouse=True)
def restore(tmp_path):
    """Put the `parameter` module back exactly as it was.

    `get_params` writes dozens of module globals -- polymerase, min_tm,
    mg_conc, genome paths -- so restoring only `json_file` leaves a bst
    configuration in place for whatever test runs next. That is how this file
    broke the plasmid golden snapshot on its first run: the snapshot passed
    alone and failed in the suite.
    """
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import filter as filter_mod

    snapshot = dict(vars(parameter))
    yield

    current = vars(parameter)
    for key, value in snapshot.items():
        current[key] = value
    for key in [k for k in current if k not in snapshot]:
        del current[key]

    # `pipeline` caches per-run state on the module object and `filter` caches
    # a ReactionConditions singleton, so restoring `parameter` alone is not
    # enough -- the next test would read this file's genome config and its
    # bst-derived conditions. Mirrors tests/integration/conftest.py's
    # reset_pipeline_state.
    pipeline_mod._initialized = False
    for attr in (
        "fg_prefixes",
        "bg_prefixes",
        "fg_genomes",
        "bg_genomes",
        "fg_seq_lengths",
        "bg_seq_lengths",
        "fg_circular",
        "bg_circular",
    ):
        setattr(pipeline_mod, attr, None)
    filter_mod.reset_reaction_conditions()


def test_an_unset_window_takes_the_polymerase_default(tmp_path):
    _load(tmp_path, polymerase="bst")

    assert (parameter.min_tm, parameter.max_tm) == parameter.default_tm_range("bst")


def test_an_explicit_window_is_respected(tmp_path):
    """A user narrowing the window on purpose must not be overridden by the
    enzyme default."""
    _load(tmp_path, polymerase="bst", min_tm=55.0, max_tm=70.0)

    assert (parameter.min_tm, parameter.max_tm) == (55.0, 70.0)


def test_a_zero_lower_bound_is_honoured_not_treated_as_absent(tmp_path):
    """`or` read a configured 0.0 as unset and substituted 15.

    0 C is a legitimate lower bound for someone who wants no lower filtering,
    and silently replacing it is the same class of fault as every other
    sentinel-versus-value confusion in this codebase.
    """
    _load(tmp_path, polymerase="phi29", min_tm=0.0, max_tm=40.0)

    assert parameter.min_tm == 0.0


# ----------------------------------------------------------------------
# The filter uses the resolved window, not its own copy
# ----------------------------------------------------------------------


def test_filter_does_not_hardcode_its_own_fallback():
    """The fallback was a fourth independent answer to "what Tm is
    acceptable", alongside params.json and the two optimizer presets.

    Checked against the parsed source rather than the text, so the docstring
    that quotes the old lines for context does not trip it.
    """
    import ast
    from pathlib import Path

    tree = ast.parse((Path(__file__).parent.parent / "neoswga" / "core" / "filter.py").read_text())

    offenders = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.BoolOp):
            continue
        if not isinstance(node.value.op, ast.Or):
            continue
        dumped = ast.dump(node.value)
        if "min_tm" in dumped or "max_tm" in dumped:
            offenders.append(ast.unparse(node))

    assert not offenders, f"filter.py still resolves the Tm window with `or`: {offenders}"


def test_filter_window_follows_the_polymerase(tmp_path):
    """End to end: the value the filter actually applies moves with the enzyme."""
    from neoswga.core import filter as filter_mod

    _load(tmp_path, polymerase="bst")
    filter_mod.reset_reaction_conditions()
    bst_window = filter_mod._resolve_tm_window()

    _load(tmp_path, polymerase="phi29")
    filter_mod.reset_reaction_conditions()
    phi29_window = filter_mod._resolve_tm_window()

    assert bst_window != phi29_window
    assert bst_window[0] > phi29_window[0]


# ----------------------------------------------------------------------
# The optimizer's own thermodynamic screen
# ----------------------------------------------------------------------


def _hybrid(**kwargs):
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    return HybridOptimizer(position_cache=None, fg_prefixes=["x"], fg_seq_lengths=[1000], **kwargs)


def test_optimizer_screen_respects_an_explicit_window():
    """Stage-0 thermodynamic screening built its criteria from the polymerase
    preset alone, so a user who widened `max_tm` in params.json had it silently
    narrowed back to the preset. The filter step honoured the request and the
    optimizer then undid it.
    """
    optimizer = _hybrid(polymerase="equiphi29", min_tm=25.0, max_tm=70.0)
    criteria = optimizer._thermo_criteria()

    assert (criteria.min_tm, criteria.max_tm) == (25.0, 70.0)


def test_optimizer_screen_falls_back_to_the_preset():
    """With nothing configured, the enzyme's own range is the right answer."""
    from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS

    optimizer = _hybrid(polymerase="equiphi29")
    criteria = optimizer._thermo_criteria()
    preset = POLYMERASE_PRESETS["equiphi29"]

    assert (criteria.min_tm, criteria.max_tm) == (
        preset.min_primer_tm,
        preset.max_primer_tm,
    )


def test_optimizer_screen_does_not_run_at_zero_magnesium():
    """`mg_conc=0.0` was hardcoded into the screening criteria.

    Every hairpin, homodimer and heterodimer energy this screen rejects
    primers on was therefore computed for a reaction containing no magnesium --
    the same zero-magnesium fault schema v2 corrected in the presets. phi29
    buffer is 10 mM, and Mg2+ is the dominant term in the salt correction.
    """
    optimizer = _hybrid(polymerase="phi29")
    criteria = optimizer._thermo_criteria()

    assert criteria.mg_conc > 0.0, "thermodynamic screening still assumes no magnesium"


def test_optimizer_screen_uses_the_configured_buffer():
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=30.0, polymerase="phi29", mg_conc=6.0, na_conc=80.0)
    criteria = _hybrid(polymerase="phi29", conditions=conditions)._thermo_criteria()

    assert criteria.mg_conc == 6.0
    assert criteria.na_conc == 80.0
