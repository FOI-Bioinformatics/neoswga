"""`optimization_method` in params.json must select the optimizer.

Known Issue 8, and the last surviving instance of this repository's recurring
defect: a key declared in `params.schema.json`, documented in CLAUDE.md,
accepted by the validator, and read by nothing.

`get_params` put it in the returned dict and assigned no module global, so
`hasattr(parameter, "optimization_method")` was False. `run_step4` passed
`args.optimization_method` straight through, and that argument's argparse
default was the literal `"hybrid"`. The flag's default therefore beat the
config every time:

    params.json "optimization_method": "dominating-set"
      -> parameter global : <UNSET>
      -> optimizer run    : hybrid

It cost more than provenance. `hybrid` returns a set identical to
`dominating-set` (Jaccard 1.000) at 7.8x the cost at 32 primers and 260x at 128
-- 2,239 s against 8.6 s -- so every params.json user was pinned to the slowest
method for the same answer.

The sentinel matters. `"hybrid"` is both the argparse default AND a legitimate
explicit choice, so "did the user type it?" cannot be answered by comparing to
the default. The flag now defaults to None and the resolution is explicit:
flag, then params.json, then `"hybrid"`.

`design` compounded it from a second direction: that subparser offers no
`--optimization-method` at all and `run_design` filled the missing attribute
with a hardcoded `"hybrid"`, so the one-shot path was pinned by two independent
routes.
"""

import json

import pytest

from neoswga.core import parameter


class _Args:
    """The attributes the resolver reads."""

    optimization_method = None


_MINIMAL = {
    "schema_version": 2,
    "fg_genomes": ["fg.fna"],
    "fg_prefixes": ["fg"],
    "bg_genomes": [],
    "bg_prefixes": [],
    "data_dir": "./",
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "min_fg_freq": 1e-5,
    "max_bg_freq": 5e-6,
    "min_tm": 15,
    "max_tm": 45,
    "max_gini": 0.6,
    "max_primer": 500,
    "cpus": 1,
    "fg_seq_lengths": [1_000_000],
}


def _load(tmp_path, monkeypatch, **overrides):
    params = tmp_path / "params.json"
    params.write_text(json.dumps({**_MINIMAL, **overrides}))
    monkeypatch.chdir(tmp_path)

    class _A:
        """Every flag unset; params.json is the only source. Same idiom as
        `test_params_json_routes_optional_keys.py`."""

        json_file = str(params)

        def __getattr__(self, name):
            return None

    return parameter.get_params(_A())


# ---------------------------------------------------------------------------
# The routing
# ---------------------------------------------------------------------------


def test_the_key_reaches_the_parameter_module(tmp_path, monkeypatch):
    """`get_params` assigned no global, so every read site defaulted."""
    _load(tmp_path, monkeypatch, optimization_method="dominating-set")

    assert getattr(parameter, "optimization_method", None) == "dominating-set"


def test_absent_leaves_it_unset_rather_than_guessing(tmp_path, monkeypatch):
    """Loaded after a run that set it, so a leaked global would show."""
    _load(tmp_path, monkeypatch, optimization_method="clique")
    _load(tmp_path, monkeypatch)

    assert getattr(parameter, "optimization_method", None) is None


# ---------------------------------------------------------------------------
# The resolution, flag over config over default
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "flag,configured,expected",
    [
        (None, None, "hybrid"),  # nothing set anywhere
        (None, "dominating-set", "dominating-set"),  # params.json alone
        ("network", None, "network"),  # flag alone
        ("network", "dominating-set", "network"),  # flag beats params.json
        ("hybrid", "dominating-set", "hybrid"),  # the sentinel case
    ],
)
def test_the_flag_beats_the_config_beats_the_default(monkeypatch, flag, configured, expected):
    """The last row is why a sentinel is needed rather than a default string.

    `"hybrid"` typed explicitly must win over a configured `dominating-set`.
    With `"hybrid"` as the argparse default there is no way to tell that apart
    from the user saying nothing.
    """
    from neoswga.cli.pipeline import resolve_optimization_method

    monkeypatch.setattr(parameter, "optimization_method", configured, raising=False)
    args = _Args()
    args.optimization_method = flag

    assert resolve_optimization_method(args) == expected


def test_the_flag_default_is_a_sentinel_not_a_method():
    """If argparse defaults this to a real method name again, the resolution
    above becomes undecidable and the config silently loses."""
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    subparsers = [a for a in parser._subparsers._group_actions][0]
    optimize = subparsers.choices["optimize"]

    default = {a.dest: a.default for a in optimize._actions}["optimization_method"]
    assert default is None, (
        f"--optimization-method defaults to {default!r}; it must default to None "
        "so an explicit choice can be told from an absent one"
    )


def test_design_does_not_hardcode_the_method(monkeypatch):
    """`design` offers no flag, so it must fall through to params.json rather
    than pinning `hybrid` a second time."""
    from neoswga.cli.pipeline import resolve_optimization_method

    monkeypatch.setattr(parameter, "optimization_method", "clique", raising=False)

    class _DesignArgs:
        pass  # the subparser defines no such attribute at all

    assert resolve_optimization_method(_DesignArgs()) == "clique"
