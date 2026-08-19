"""Selecting a different enzyme has to bring that enzyme's window with it.

`_apply_gc_adaptive_defaults` picks equiphi29 for any target at or above ~41%
GC, which is most bacterial genomes -- so this is the common path, not an edge
case. It set `parameter.polymerase` and overwrote `min_k`/`max_k` from its own
table, and left everything else at phi29's values:

- `min_tm`/`max_tm` stayed at phi29's 20-50 C. equiphi29's window is 37-62.
  A design meant for a 42 C reaction was filtered for primers that melt as low
  as 20 C.
- `reaction_temp` stayed at 30 C. The branch meant to set 42 applied only when
  `reaction_temp is None`, and `get_params` always assigns a float, so it never
  fired. The run was equiphi29 at phi29's temperature.
- `mg_conc` and `dtt_mm` likewise.

The k-mer range is the one thing it did move, and that exposed the second
problem: `gc_adaptive_strategy` carries its own polymerase table, and its ranges
disagree with the registry's (phi29 8-11 against 6-12, equiphi29 11-15 against
10-18). Two tables, no test comparing their values.

The strategy's narrowing is worth keeping -- it is GC-specific where the
registry is only enzyme-specific -- so the resolution is to retune from the
registry first and then let the strategy narrow *within* the enzyme's range,
rather than have the two overwrite each other in whichever order they happen to
run.
"""

import json

import pytest

from neoswga.core import parameter
from neoswga.core.registry import views


@pytest.fixture(autouse=True)
def restore():
    snapshot = dict(vars(parameter))
    yield
    current = vars(parameter)
    for key, value in snapshot.items():
        current[key] = value
    for key in [k for k in current if k not in snapshot]:
        del current[key]


@pytest.fixture
def load(tmp_path):
    """Run the real lazy-init path on a params.json of a given GC content."""
    import neoswga.core.pipeline as pipeline_mod

    def _load(gc_fraction, **overrides):
        # A sequence with the requested GC, so `genome_gc` is computed not faked.
        gc_bases = int(round(gc_fraction * 100))
        unit = "GC" * (gc_bases // 2) + "AT" * ((100 - gc_bases) // 2)
        genome = tmp_path / "g.fasta"
        genome.write_text(">g\n" + unit * 16 + "\n")

        params = {
            "schema_version": 2,
            "data_dir": str(tmp_path),
            "fg_genomes": [str(genome)],
            "fg_prefixes": [str(tmp_path / "g")],
            "fg_seq_lengths": [len(unit) * 16],
        }
        params.update(overrides)
        (tmp_path / "params.json").write_text(json.dumps(params))

        parameter.json_file = str(tmp_path / "params.json")
        pipeline_mod._initialized = False
        pipeline_mod._initialize()
        return parameter

    return _load


def test_a_gc_adaptive_swap_brings_the_tm_window_with_it(load):
    """The regression: equiphi29 selected, phi29's Tm window left in place."""
    param = load(0.62)

    assert param.polymerase == "equiphi29"
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("equiphi29")


def test_a_gc_adaptive_swap_brings_the_temperature_with_it(load):
    """The dead branch: 42 C was applied only when `reaction_temp is None`,
    which `get_params` guarantees it is not."""
    param = load(0.62)

    assert param.reaction_temp == parameter.default_reaction_temp("equiphi29")


def test_a_gc_adaptive_swap_brings_the_buffer_with_it(load):
    param = load(0.62)

    assert param.mg_conc == parameter.default_mg_conc("equiphi29")
    assert param.dtt_mm == parameter.default_dtt_mm("equiphi29")


def test_the_resulting_conditions_are_constructible(load):
    """The practical consequence of leaving the temperature behind: for an
    enzyme whose hard range excludes 30 C this is a crash, and for equiphi29 it
    is a silent 12 C error."""
    from neoswga.core.reaction_conditions import ReactionConditions

    param = load(0.62)
    conditions = ReactionConditions(temp=param.reaction_temp, polymerase=param.polymerase)

    assert conditions.temp == parameter.default_reaction_temp("equiphi29")


def test_the_adaptive_kmer_range_stays_inside_the_enzymes_range(load):
    """The two tables disagree; neither should be able to contradict the other.
    The strategy narrows within the registry's range rather than replacing it.
    """
    param = load(0.62)

    reg_min, reg_max = views.primer_length_ranges()["equiphi29"]
    assert reg_min <= param.min_k <= param.max_k <= reg_max


def test_the_adaptive_narrowing_survives_a_later_retune(load):
    """Once the strategy has chosen a range deliberately, it is no longer a
    registry default and must not be reset by a subsequent enzyme change."""
    param = load(0.62)
    narrowed = (param.min_k, param.max_k)

    parameter.retune_for_polymerase(param, "equiphi29")

    assert (param.min_k, param.max_k) == narrowed


def test_an_explicit_polymerase_is_not_overridden(load):
    """The strategy defers to the user; that behaviour predates this change and
    has to survive it."""
    param = load(0.62, polymerase="phi29")

    assert param.polymerase == "phi29"
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("phi29")


def test_an_explicit_window_is_not_overridden(load):
    param = load(0.62, min_tm=33.0, max_tm=48.0)

    assert param.polymerase == "equiphi29"
    assert (param.min_tm, param.max_tm) == (33.0, 48.0)


def test_an_at_rich_target_stays_on_phi29(load):
    param = load(0.24)

    assert param.polymerase == "phi29"
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("phi29")
