"""`filter --polymerase` and `filter --preset` have to move the windows too.

`run_step2` calls `pipeline._initialize()` first, which runs `get_params` and
derives the Tm window, primer-length range and buffer from whatever polymerase
params.json named. The preset application and the `--polymerase` merge happen
*after* that, and both set `parameter.polymerase` on its own.

So `neoswga filter -j params.json --polymerase bst` filtered a bst design
through phi29's 20-50 C Tm window and 6-12 bp length range: a window built for a
30 C reaction applied to one running at 63 C, on primer lengths bst cannot use.
Nothing warned, and the resulting candidate pool looks perfectly ordinary.

Same mechanism as the `optimize --polymerase` and GC-adaptive faults; the shared
resolution is `parameter.retune_for_polymerase`, which moves what the enzyme
decides and leaves what the user pinned.
"""

import json

import pytest

from neoswga.core import parameter
from neoswga.core.registry import POLYMERASES


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
def run_filter(tmp_path, monkeypatch):
    """Drive `run_step2` up to the point the filtering itself would start."""
    import neoswga.cli.pipeline as cli_pipeline
    import neoswga.core.pipeline as pipeline_mod

    def _run(*, params=None, **cli_args):
        genome = tmp_path / "g.fasta"
        genome.write_text(">g\n" + "ACGT" * 400 + "\n")
        payload = {
            "schema_version": 2,
            "data_dir": str(tmp_path),
            "fg_genomes": [str(genome)],
            "fg_prefixes": [str(tmp_path / "g")],
            "fg_seq_lengths": [1600],
            "polymerase": "phi29",
        }
        payload.update(params or {})
        path = tmp_path / "params.json"
        path.write_text(json.dumps(payload))

        # The filtering itself needs k-mer count files this fixture does not
        # build, and everything under test has happened before it runs. A no-op
        # rather than a raise, because `run_step2` wraps its body in broad
        # handlers that would turn an exception into SystemExit and hide it.
        monkeypatch.setattr(pipeline_mod, "step2", lambda *a, **k: None, raising=False)
        pipeline_mod._initialized = False

        args = _Args(json_file=str(path), **cli_args)
        try:
            cli_pipeline.run_step2(args)
        except SystemExit:  # pragma: no cover - only if a later stage objects
            pass
        return parameter

    return _run


class _Args:
    """argparse-style namespace: every unknown attribute reads as None."""

    def __init__(self, **kwargs):
        self.quiet = True
        self.__dict__.update(kwargs)

    def __getattr__(self, name):
        return None


@pytest.mark.parametrize("target", sorted(set(POLYMERASES) - {"phi29"}))
def test_the_polymerase_flag_moves_the_tm_window(run_filter, target):
    """The regression."""
    param = run_filter(polymerase=target)

    assert param.polymerase == target
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range(target)


@pytest.mark.parametrize("target", sorted(set(POLYMERASES) - {"phi29"}))
def test_the_polymerase_flag_moves_the_buffer_and_temperature(run_filter, target):
    param = run_filter(polymerase=target)

    assert param.mg_conc == parameter.default_mg_conc(target)
    assert param.reaction_temp == parameter.default_reaction_temp(target)


def test_a_preset_moves_the_window_with_its_enzyme(run_filter):
    """`--preset bst` names an enzyme just as surely as `--polymerase bst`."""
    param = run_filter(preset="bst")

    assert param.polymerase == "bst"
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("bst")


def test_an_explicit_window_survives_the_flag(run_filter):
    """A user narrowing the window on purpose is not overridden by the enzyme
    default -- the same guarantee `get_params` already gives."""
    param = run_filter(params={"min_tm": 52.0, "max_tm": 68.0}, polymerase="bst")

    assert param.polymerase == "bst"
    assert (param.min_tm, param.max_tm) == (52.0, 68.0)


def test_a_flag_supplied_window_survives_the_flag(run_filter):
    param = run_filter(polymerase="bst", min_tm=51.0, max_tm=69.0)

    assert (param.min_tm, param.max_tm) == (51.0, 69.0)


def test_no_polymerase_flag_leaves_the_configured_enzyme(run_filter):
    param = run_filter(params={"polymerase": "equiphi29"})

    assert param.polymerase == "equiphi29"
    assert (param.min_tm, param.max_tm) == parameter.default_tm_range("equiphi29")
