"""Two settings behaved differently in params.json than on the command line.

Both are the same shape: a fix reached one route and not the other, and nothing
warned, because `params.schema.json` sets `additionalProperties: true`.

**`coverage_reach`.** `unified_optimizer.run_optimization` reads it with
`getattr(parameter, "coverage_reach", None)`, but `get_params` never assigned
that global -- the value arrived in the returned dict and in `_json_data`, and
stopped there. So `"coverage_reach": 20000` in params.json produced output
byte-identical to omitting the key, while `--coverage-reach 20000` worked. On
Prevotella vs human chr21 that is the difference between a reported
`fg_coverage` of 0.3816 and 0.9513 for the same primers, since coverage is
roughly linear in reach over this range. It is the single largest lever on any
coverage figure the tool prints, so a route that silently ignores it is worse
than one that rejects it.

**`num_primers` / `target_set_size`.** `PARAM_RANGES` capped both at 50 while
`params.schema.json` allows 200. The validator gate runs first, so a params.json
asking for 64 exited with "Value 64 outside valid range [1, 50]" while
`--num-primers 64` ran fine. The ceiling was raised to 200 in the 2026-08-19
audit precisely because a >95% design on a 3.2 Mb target needs 31-96 primers,
and that change did not reach the validator -- so the capability existed and was
unreachable from the documented configuration file. Measured on the same
dataset: 64 primers reach 0.661 coverage against 0.538 at 32, and 128 reach
0.810.

See docs/AUDIT_2026-08_alternatives_and_scaling.md (F1, F2).
"""

import json

import pytest

from neoswga.core import parameter
from neoswga.core.param_validator import PARAM_RANGES, ParamValidator, ValidationLevel

# ---------------------------------------------------------------------------
# coverage_reach
# ---------------------------------------------------------------------------


class _Args:
    """Stand-in for the argparse namespace `get_params` consumes.

    Every flag is unset, which is the case that matters here: with no
    `--coverage-reach` on the command line, params.json is the only source.
    """

    def __getattr__(self, name):
        return None


_MINIMAL_PARAMS = {
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
    "iterations": 8,
    "max_sets": 5,
    "fg_seq_lengths": [1_000_000],
}


def _load_params(tmp_path, monkeypatch, **overrides):
    """Run `get_params` against a params.json holding `overrides`."""
    params_file = tmp_path / "params.json"
    params_file.write_text(json.dumps({**_MINIMAL_PARAMS, **overrides}))

    monkeypatch.chdir(tmp_path)
    args = _Args()
    args.json_file = str(params_file)
    return parameter.get_params(args)


def test_coverage_reach_in_params_json_reaches_the_parameter_module(tmp_path, monkeypatch):
    """The value must land where `run_optimization` looks for it.

    `unified_optimizer.py` reads `getattr(parameter, "coverage_reach", None)`,
    so putting the key only in the returned dict leaves the optimizer on the
    polymerase default with no warning.
    """
    _load_params(tmp_path, monkeypatch, coverage_reach=20_000)

    assert getattr(parameter, "coverage_reach", None) == 20_000


def test_coverage_reach_absent_from_params_json_stays_unset(tmp_path, monkeypatch):
    """Absent means "use the polymerase default", not zero.

    `resolve_coverage_reach` treats `None` as "no override". A leaked 0 from a
    previous run would read as a reach of zero and collapse coverage.
    """
    _load_params(tmp_path, monkeypatch, coverage_reach=20_000)
    _load_params(tmp_path, monkeypatch)

    assert getattr(parameter, "coverage_reach", None) is None


# ---------------------------------------------------------------------------
# num_primers / target_set_size
# ---------------------------------------------------------------------------

_REQUIRED = {
    "fg_genomes": ["fg.fna"],
    "fg_prefixes": ["fg"],
    "data_dir": "results",
}


def _range_errors(params, name):
    messages = ParamValidator().validate_params({**_REQUIRED, **params})
    return [
        m
        for m in messages
        if m.level == ValidationLevel.ERROR and name in m.parameter and "range" in m.message
    ]


@pytest.mark.parametrize("name", ["num_primers", "target_set_size"])
@pytest.mark.parametrize("size", [51, 64, 96, 128, 200])
def test_set_sizes_the_schema_allows_are_accepted(name, size):
    """The schema allows up to 200; the validator must not be stricter."""
    assert _range_errors({name: size}, name) == []


@pytest.mark.parametrize("name", ["num_primers", "target_set_size"])
def test_set_sizes_past_the_schema_ceiling_are_still_rejected(name):
    """Widening the validator must not remove the ceiling altogether."""
    assert _range_errors({name: 201}, name) != []


@pytest.mark.parametrize("name", ["num_primers", "target_set_size"])
def test_validator_ceiling_matches_the_schema(name):
    """The two must not drift apart again -- that drift is the whole defect."""
    schema = json.loads(
        (
            __import__("pathlib").Path(parameter.__file__).parent / "schema" / "params.schema.json"
        ).read_text()
    )
    assert PARAM_RANGES[name][1] == schema["properties"][name]["maximum"]
