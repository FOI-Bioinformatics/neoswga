"""Four more params.json keys read from a global `get_params` never assigned.

Same shape as `coverage_reach` (see
`test_params_json_routes_reach_the_pipeline.py`): the read site uses
`getattr(parameter, name, default)`, `get_params` puts the value in the returned
dict and nowhere else, and `additionalProperties: true` means nothing warns. The
key is documented, accepted, and inert.

- `occupancy_ranking` and `occupancy_shortlist` (`pipeline.py`) are declared in
  `params.schema.json` with defaults, and reachable only by monkeypatch.
- `max_mismatches` (`pipeline.py`) is not in the schema at all.
- `sampled_index_path` (`filter.py`) has a CLI route that works, and
  `filter.py`'s own "index not found" message tells the user to set it in
  params.json -- a route that did not exist.

`sampled_index_path` raises a precedence question, and the answer is ordering:
`run_step2` calls `pipeline._initialize()` -- and so `get_params` -- before it
merges any flag onto the parameter module, so a flag set afterwards wins without
`get_params` needing to know about it. That is why these keys can be assigned
unconditionally, the way `coverage_reach` is, rather than needing a
"the CLI already set this" registry.

See docs/AUDIT_2026-08_alternatives_and_scaling.md (F1).
"""

import json
import pathlib

import pytest

from neoswga.core import parameter
from neoswga.core.pipeline import DEFAULT_OCCUPANCY_SHORTLIST, _configured_positive_int


def _schema() -> dict:
    path = pathlib.Path(parameter.__file__).parent / "schema" / "params.schema.json"
    return json.loads(path.read_text())


def _configured_positive_int_accepts(value) -> bool:
    """Whether `_configured_positive_int` keeps `value` or falls back."""
    sentinel = -12345
    try:
        parameter.probe_key = value
        return _configured_positive_int("probe_key", sentinel) != sentinel
    finally:
        del parameter.probe_key


class _Args:
    """Stand-in for the argparse namespace `get_params` consumes.

    Every flag unset -- the case that matters here, where params.json is the
    only source.
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


class TestOccupancyRanking:
    """`pipeline.py` reads `getattr(parameter, "occupancy_ranking", True)`."""

    def test_disabling_it_in_params_json_reaches_the_parameter_module(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, occupancy_ranking=False)

        assert getattr(parameter, "occupancy_ranking", None) is False

    def test_absent_leaves_ranking_enabled(self, tmp_path, monkeypatch):
        """The schema's default is true, and the read site's fallback is True.

        Loaded after a run that disabled it, so a leaked global would show.
        """
        _load_params(tmp_path, monkeypatch, occupancy_ranking=False)
        _load_params(tmp_path, monkeypatch)

        assert getattr(parameter, "occupancy_ranking", True) is not False


class TestOccupancyShortlist:
    """Read through `_configured_positive_int`, so assert at the read site."""

    def test_shortlist_size_in_params_json_reaches_the_read_site(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, occupancy_shortlist=1000)

        assert _configured_positive_int("occupancy_shortlist", DEFAULT_OCCUPANCY_SHORTLIST) == 1000

    def test_absent_falls_back_to_the_documented_default(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, occupancy_shortlist=1000)
        _load_params(tmp_path, monkeypatch)

        assert (
            _configured_positive_int("occupancy_shortlist", DEFAULT_OCCUPANCY_SHORTLIST)
            == DEFAULT_OCCUPANCY_SHORTLIST
        )


class TestMaxMismatches:
    """Sets how many mismatch classes the occupancy ranking counts."""

    def test_max_mismatches_in_params_json_reaches_the_read_site(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, max_mismatches=2)

        assert _configured_positive_int("max_mismatches", 1) == 2

    def test_absent_falls_back_to_one_mismatch_class(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, max_mismatches=2)
        _load_params(tmp_path, monkeypatch)

        assert _configured_positive_int("max_mismatches", 1) == 1

    def test_declared_in_the_schema(self):
        """Undeclared keys are accepted silently and cannot be validated."""
        assert "max_mismatches" in _schema()["properties"]

    @pytest.mark.parametrize("key", ["max_mismatches", "occupancy_shortlist"])
    def test_schema_floor_matches_what_the_read_site_accepts(self, key):
        """A schema that admits a value the read site drops is the same defect.

        `_configured_positive_int` substitutes its default for anything that is
        not a positive int, silently. Advertising `minimum: 0` would make
        `max_mismatches: 0` a documented setting that does nothing.
        """
        floor = _schema()["properties"][key]["minimum"]

        assert _configured_positive_int_accepts(floor)
        assert not _configured_positive_int_accepts(floor - 1)


class TestSampledIndexPath:
    """`filter.py` tells the user to set this in params.json. It must work."""

    def test_path_in_params_json_reaches_the_parameter_module(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, sampled_index_path="bg_sampled.pkl")

        assert getattr(parameter, "sampled_index_path", None) == "bg_sampled.pkl"

    def test_absent_stays_unset(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, sampled_index_path="bg_sampled.pkl")
        _load_params(tmp_path, monkeypatch)

        assert getattr(parameter, "sampled_index_path", None) is None

    def test_declared_in_the_schema(self):
        assert "sampled_index_path" in _schema()["properties"]


class TestCommandLineStillWins:
    """A flag beats params.json, and ordering is what delivers that.

    `run_step2` calls `pipeline._initialize()` (hence `get_params`) before it
    merges flags onto the parameter module, so the flag lands last. This pins
    the ordering rather than the merge: assigning these keys unconditionally in
    `get_params` is only safe while `get_params` runs first.
    """

    @pytest.mark.parametrize(
        "name, cli_value",
        [
            ("sampled_index_path", "from_cli.pkl"),
            ("bloom_filter_path", "from_cli_bloom.pkl"),
            ("use_bloom_filter", True),
        ],
    )
    def test_a_flag_applied_after_get_params_wins(self, name, cli_value, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, sampled_index_path="from_params.pkl")

        # What cli/pipeline.py does, at the point it does it.
        monkeypatch.setattr(parameter, name, cli_value, raising=False)

        assert getattr(parameter, name, None) == cli_value
