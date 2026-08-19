"""Changing the enzyme after `get_params` has to move what the enzyme decided.

`get_params` derives a group of values from the polymerase: the Tm window, the
primer-length range, magnesium, DTT and the reaction temperature. Anything the
user set explicitly wins; the rest are filled in from the registry.

Two paths then change the polymerase *after* that has happened, and neither
re-derived anything:

- `filter --preset bst` / `filter --polymerase bst`. `pipeline._initialize()`
  runs `get_params` first, so the windows are phi29's; the merge that follows
  swaps only `parameter.polymerase`. A bst design was filtered through phi29's
  20-50 C Tm window and 6-12 bp length range -- a window built for a 30 C
  reaction applied to one running at 63 C.
- The GC-adaptive path in `core/pipeline.py`, which can select equiphi29 on a
  GC-rich target. It overwrote `min_k`/`max_k` with its own table but left
  `min_tm`/`max_tm` at phi29's 20-50, where equiphi29's window is 37-62. Its
  temperature branch was dead as well: it applied 42 C only when
  `reaction_temp is None`, and `get_params` always assigns a float. So the run
  was equiphi29, at 30 C, on phi29's Tm window.

The distinction that makes this fixable is between a value the user chose and
one the polymerase supplied. `get_params` records the latter in
`polymerase_derived_fields`, and `retune_for_polymerase` re-derives exactly
those. A user who pinned `min_tm` keeps it across an enzyme change; a user who
did not gets the window that matches the enzyme actually being used.
"""

import json

import pytest

from neoswga.core import parameter
from neoswga.core.registry import POLYMERASES

DERIVED = ("min_tm", "max_tm", "min_k", "max_k", "mg_conc", "reaction_temp", "dtt_mm")


@pytest.fixture(autouse=True)
def restore():
    """`get_params` writes dozens of module globals; put them all back.

    Same reasoning as `tests/test_tm_window_follows_the_polymerase.py`: leaving
    a bst configuration behind breaks whatever runs next, and does so only when
    the suite runs in full.
    """
    snapshot = dict(vars(parameter))
    yield
    current = vars(parameter)
    for key, value in snapshot.items():
        current[key] = value
    for key in [k for k in current if k not in snapshot]:
        del current[key]


@pytest.fixture
def loaded(tmp_path):
    """Load a minimal params.json through the real lazy-init path.

    The polymerase is pinned unless a test says otherwise, so the GC-adaptive
    strategy defers and these tests measure retuning alone. Its own interaction
    with retuning is covered by `tests/test_gc_adaptive_swap_retunes.py`.
    """
    import neoswga.core.pipeline as pipeline_mod

    def load(**overrides):
        genome = tmp_path / "g.fasta"
        genome.write_text(">g\n" + "ACGT" * 400 + "\n")
        params = {
            "schema_version": 2,
            "data_dir": str(tmp_path),
            "fg_genomes": [str(genome)],
            "fg_prefixes": [str(tmp_path / "g")],
            "fg_seq_lengths": [1600],
            "polymerase": "phi29",
        }
        params.update(overrides)
        path = tmp_path / "params.json"
        path.write_text(json.dumps(params))

        parameter.json_file = str(path)
        pipeline_mod._initialized = False
        pipeline_mod._initialize()
        return parameter

    return load


# ----------------------------------------------------------------------
# What the polymerase supplied is re-derived
# ----------------------------------------------------------------------


@pytest.mark.parametrize("target", sorted(POLYMERASES))
def test_retuning_moves_the_tm_window_to_the_new_enzyme(loaded, target):
    """The regression: a bst design filtered on phi29's window."""
    param = loaded()  # defaults to phi29
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("phi29")

    parameter.retune_for_polymerase(param, target)

    assert (param.min_tm, param.max_tm) == parameter.default_tm_range(target)
    assert param.polymerase == target


class Stand_in:
    """A parameter-module stand-in with everything still at its default.

    The module-level fixture cannot serve here: `_apply_gc_adaptive_defaults`
    always chooses a k-mer range, so after a real init `min_k`/`max_k` are a
    deliberate choice rather than a registry default and correctly stop moving.
    This isolates the mechanics from that.
    """

    def __init__(self, polymerase="phi29", derived=None):
        self.polymerase = polymerase
        self.polymerase_derived_fields = (
            set(parameter.POLYMERASE_DERIVED_FIELDS) if derived is None else set(derived)
        )
        for name, value in parameter.polymerase_defaults(polymerase).items():
            setattr(self, name, value)


@pytest.mark.parametrize("target", sorted(POLYMERASES))
def test_retuning_moves_the_primer_length_range(target):
    from neoswga.core.registry import views

    param = Stand_in()
    parameter.retune_for_polymerase(param, target)

    assert (param.min_k, param.max_k) == views.primer_length_ranges()[target]


def test_a_deliberate_primer_length_range_is_left_alone(loaded):
    """Once the GC-adaptive strategy or the user has chosen a range, an enzyme
    change must not silently widen it back to the registry default."""
    param = loaded()
    chosen = (param.min_k, param.max_k)
    assert "min_k" not in param.polymerase_derived_fields

    parameter.retune_for_polymerase(param, "equiphi29")

    assert (param.min_k, param.max_k) == chosen


@pytest.mark.parametrize("target", sorted(POLYMERASES))
def test_retuning_moves_the_buffer_and_the_temperature(loaded, target):
    param = loaded()
    parameter.retune_for_polymerase(param, target)

    assert param.mg_conc == parameter.default_mg_conc(target)
    assert param.dtt_mm == parameter.default_dtt_mm(target)
    assert param.reaction_temp == parameter.default_reaction_temp(target)


def test_the_retuned_temperature_is_inside_the_new_enzymes_hard_range(loaded):
    """`bst` at phi29's 30 C is rejected by `ReactionConditions`, so leaving the
    temperature behind turns an enzyme change into a crash or a silent fallback
    to baseline effects."""
    from neoswga.core.reaction_conditions import ReactionConditions

    param = loaded()
    for target in sorted(POLYMERASES):
        parameter.retune_for_polymerase(param, target)
        # Constructs without raising, which is the whole assertion.
        ReactionConditions(temp=param.reaction_temp, polymerase=target)


# ----------------------------------------------------------------------
# What the user chose is left alone
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "field,value",
    [
        ("min_tm", 33.0),
        ("max_tm", 71.0),
        ("min_k", 9),
        ("max_k", 14),
        ("mg_conc", 6.5),
        ("reaction_temp", 39.0),
        ("dtt_mm", 2.0),
    ],
)
def test_an_explicit_value_survives_an_enzyme_change(loaded, field, value):
    """Re-deriving everything would silently discard a pinned value, which is
    the opposite failure and just as quiet."""
    param = loaded(**{field: value})
    assert field not in param.polymerase_derived_fields

    parameter.retune_for_polymerase(param, "bst")

    assert getattr(param, field) == value


def test_the_derived_set_is_recorded_when_the_user_sets_nothing(loaded):
    """Everything the user left to the registry is recorded as movable.

    `min_k`/`max_k` are the exception on any real params.json, because
    `_apply_gc_adaptive_defaults` always picks a range and that is a choice, not
    a default. The chemistry fields have no such second author.
    """
    param = loaded()

    assert {"min_tm", "max_tm", "mg_conc", "reaction_temp", "dtt_mm"} <= (
        param.polymerase_derived_fields
    )


def test_retuning_to_the_same_enzyme_changes_nothing(loaded):
    param = loaded(polymerase="equiphi29")
    before = {name: getattr(param, name) for name in DERIVED}

    parameter.retune_for_polymerase(param, "equiphi29")

    assert {name: getattr(param, name) for name in DERIVED} == before


def test_an_unknown_enzyme_is_refused_rather_than_silently_applied(loaded):
    # Explicit, because the GC-adaptive path selects equiphi29 for anything at
    # or above ~40% GC and this test is about the refusal, not about that.
    param = loaded(polymerase="phi29")
    before = {name: getattr(param, name) for name in DERIVED}

    parameter.retune_for_polymerase(param, "taq")

    assert param.polymerase == "phi29"
    assert {name: getattr(param, name) for name in DERIVED} == before


def test_an_alias_resolves_before_retuning(loaded):
    param = loaded()

    parameter.retune_for_polymerase(param, "bst2.0")

    assert param.polymerase == "bst"
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("bst")
