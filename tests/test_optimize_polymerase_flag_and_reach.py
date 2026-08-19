"""Two more places the polymerase choice did not arrive intact.

**`optimize --polymerase` was accepted and discarded.** The subparser offers all
six enzymes, but `run_step4` is decorated `@params_command(merge=None)` so the
flag was never merged onto `parameter`. It printed no warning, so the only way
to notice was to compare two runs.

**`max_extension` used a sentinel-by-value.** `HybridOptimizer` treated the
literal 70000 as "caller did not choose", so a caller who *did* choose 70000 --
the phi29 processivity, a perfectly reasonable explicit value -- had it silently
replaced by the enzyme preset. For bst that is 2000, a 35-fold change to the
amplification-network reach, applied on the strength of a coincidence between a
default and a real value. `None` says "not chosen" without colliding with any
number a caller might mean.
"""

import json

import pytest

from neoswga.core.registry import POLYMERASES


class _Args:
    """argparse-style namespace: unknown attributes read as None."""

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __getattr__(self, name):
        return None


# ----------------------------------------------------------------------
# optimize --polymerase
# ----------------------------------------------------------------------


@pytest.fixture
def params_file(tmp_path):
    def write(**overrides):
        payload = {"schema_version": 2, "data_dir": str(tmp_path)}
        payload.update(overrides)
        path = tmp_path / "params.json"
        path.write_text(json.dumps(payload))
        return str(path)

    return write


@pytest.mark.parametrize("target", sorted(set(POLYMERASES) - {"phi29"}))
def test_the_optimize_flag_beats_the_params_file(params_file, target):
    """A flag the parser offers has to reach the optimizer, and has to win over
    params.json the way every other flag does."""
    from neoswga.cli._params_preread import polymerase_from_params

    parameter = _FakeParameter(params_file(polymerase="phi29"))
    args = _Args(polymerase=target)

    assert polymerase_from_params(parameter, args=args) == target


def test_no_flag_falls_back_to_the_params_file(params_file):
    from neoswga.cli._params_preread import polymerase_from_params

    parameter = _FakeParameter(params_file(polymerase="bst"))

    assert polymerase_from_params(parameter, args=_Args()) == "bst"


def test_an_unknown_flag_value_does_not_win(params_file, caplog):
    """`argparse` choices already reject this, but the resolver is called from
    more than one place and must not silently apply a name it cannot resolve."""
    import logging

    from neoswga.cli._params_preread import polymerase_from_params

    parameter = _FakeParameter(params_file(polymerase="bst"))
    with caplog.at_level(logging.WARNING):
        assert polymerase_from_params(parameter, args=_Args(polymerase="taq")) == "bst"

    assert "taq" in caplog.text


def test_the_optimize_subparser_offers_exactly_the_registry():
    """If the flag works, its choices must not be a fourth enzyme list."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    action = next(
        a
        for a in parser._subparsers._group_actions[0].choices["optimize"]._actions
        if a.dest == "polymerase"
    )

    assert set(action.choices) == set(POLYMERASES)


class _FakeParameter:
    def __init__(self, json_file):
        self.json_file = json_file
        self._json_data = {}
        self.polymerase = "phi29"


# ----------------------------------------------------------------------
# max_extension
# ----------------------------------------------------------------------


def _optimizer(**kwargs):
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    return HybridOptimizer(
        position_cache=None, fg_prefixes=["fg"], fg_seq_lengths=[100_000], **kwargs
    )


@pytest.mark.parametrize("polymerase", sorted(POLYMERASES))
def test_an_unset_max_extension_takes_the_enzyme_preset(polymerase):
    optimizer = _optimizer(polymerase=polymerase)

    assert optimizer.max_extension == POLYMERASES[polymerase].processivity_bp


@pytest.mark.parametrize("polymerase", sorted(POLYMERASES))
def test_an_explicit_max_extension_is_honoured(polymerase):
    """Including 70000, which used to be indistinguishable from "unset"."""
    for value in (70000, 12345):
        assert _optimizer(polymerase=polymerase, max_extension=value).max_extension == value
