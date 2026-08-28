"""params.json's polymerase has to survive to the optimizer.

The same defect that was fixed for the set size, twenty lines further up the
same function, was left in place for the enzyme. `run_step4` resolved it with

    polymerase = getattr(parameter, "polymerase", "phi29")

before anything had read params.json. That module global carries a default of
"phi29" and is only overwritten when `get_params()` runs, which happens lazily
inside `optimize_step4` -- after the value has been resolved and passed down as
an explicit kwarg. `unified_optimizer.run_optimization` then prefers the kwarg
over the module global, because the stale "phi29" is truthy:

    polymerase = kwargs.get("polymerase") or getattr(parameter, "polymerase", "phi29")

So a params.json asking for bst was optimized as phi29: bst's presets, Tm window
and 2 kb network reach were all replaced by phi29's, and the coverage reach used
for set-cover selection became 3000 bp rather than bst's 1000. That changes
which primers are selected, silently, with no warning anywhere.

The comment on the old line read "use adapted value, not raw JSON, so
GC-adaptive strategy recommendations are respected". That intent is real -- the
GC-adaptive path in `core/pipeline.py` can swap the enzyme -- but at this point
in a fresh `optimize` process neither params.json nor the adaptive path has run,
so there was no adapted value to respect. The resolution order below keeps that
intent where it applies: once `get_params()` has run, `parameter.polymerase` is
authoritative because it carries any adaptation; before that, params.json is.
"""

import json

import pytest

from neoswga.cli.pipeline import _polymerase_from_params


class FakeParameter:
    """Stands in for the parameter module, which is a module of globals."""

    def __init__(self, json_file=None, json_data=None, polymerase="phi29"):
        self.json_file = json_file
        self._json_data = json_data if json_data is not None else {}
        # The default that used to win over params.json.
        self.polymerase = polymerase


@pytest.fixture
def params_file(tmp_path):
    def write(payload):
        path = tmp_path / "params.json"
        path.write_text(json.dumps(payload))
        return str(path)

    return write


# ----------------------------------------------------------------------
# The regression
# ----------------------------------------------------------------------


@pytest.mark.parametrize("name", ["bst", "bst3.0", "bsu", "klenow", "equiphi29"])
def test_polymerase_from_params_json_is_used(params_file, name):
    """The bug at its smallest: `_json_data` is empty at this point in the
    command, which is exactly when the old code consulted `parameter`."""
    parameter = FakeParameter(params_file({"polymerase": name}))

    assert _polymerase_from_params(parameter) == name


def test_an_alias_resolves_to_its_canonical_key(params_file):
    """`bst2.0` is an accepted alias in the registry. Passing the alias down
    unresolved would miss every `POLYMERASE_PRESETS[...]` lookup keyed on the
    canonical name."""
    parameter = FakeParameter(params_file({"polymerase": "bst2.0"}))

    assert _polymerase_from_params(parameter) == "bst"


def test_an_already_loaded_json_defers_to_the_adapted_value(params_file):
    """Once `get_params()` has run, `parameter.polymerase` carries whatever the
    GC-adaptive path decided on top of params.json, and that is what the rest of
    the run sees. Re-reading the file here would discard the adaptation -- which
    is the behaviour the original comment was reaching for."""
    parameter = FakeParameter(
        params_file({"polymerase": "phi29"}),
        json_data={"polymerase": "phi29"},
        polymerase="equiphi29",
    )

    assert _polymerase_from_params(parameter) == "equiphi29"


# ----------------------------------------------------------------------
# It falls back rather than failing
# ----------------------------------------------------------------------


def test_a_missing_params_file_falls_back_to_the_default():
    """`optimize` can run without -j."""
    parameter = FakeParameter(json_file=None)

    assert _polymerase_from_params(parameter) == "phi29"


def test_an_absent_key_falls_back_to_the_default(params_file):
    parameter = FakeParameter(params_file({"min_k": 8}))

    assert _polymerase_from_params(parameter) == "phi29"


def test_an_unreadable_params_file_does_not_raise_here(tmp_path):
    """`get_params()` would raise on this file in a moment with a far better
    message than this helper could produce."""
    path = tmp_path / "params.json"
    path.write_text("{not json")
    parameter = FakeParameter(str(path))

    assert _polymerase_from_params(parameter) == "phi29"


@pytest.mark.parametrize("bad", ["taq", "", None, 42, True])
def test_an_unknown_enzyme_falls_back_and_says_so(params_file, caplog, bad):
    """A hand-edited params.json should not silently select a different enzyme.
    `param_validator` rejects an unknown name with a much better message once
    `get_params()` runs; this helper only has to avoid pre-empting it with a
    wrong answer."""
    import logging

    parameter = FakeParameter(params_file({"polymerase": bad}))
    with caplog.at_level(logging.WARNING):
        assert _polymerase_from_params(parameter) == "phi29"

    if bad:
        assert "polymerase" in caplog.text
