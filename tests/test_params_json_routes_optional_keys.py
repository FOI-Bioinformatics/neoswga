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

See docs/reports/AUDIT_2026-08_alternatives_and_scaling.md (F1).
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


class TestMismatchModel:
    """Which mismatch model weights a site. Added 2026-10-02 (Phase 4b).

    params.json is the only route: there is no CLI flag, deliberately, because
    the model is an extrapolation with no demonstrated benefit and the key is
    the documented way in.
    """

    def test_the_model_in_params_json_reaches_the_parameter_module(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, mismatch_model="position-dependent")

        assert getattr(parameter, "mismatch_model", None) == "position-dependent"

    def test_it_reaches_the_resolver_every_consumer_asks(self, tmp_path, monkeypatch):
        from neoswga.core.mismatch_model import POSITION_DEPENDENT, resolve_mismatch_model

        _load_params(tmp_path, monkeypatch, mismatch_model="position-dependent")

        assert resolve_mismatch_model() == POSITION_DEPENDENT

    def test_absent_leaves_the_uniform_model(self, tmp_path, monkeypatch):
        """Loaded after a run that set it, so a leaked global would show."""
        from neoswga.core.mismatch_model import UNIFORM, resolve_mismatch_model

        _load_params(tmp_path, monkeypatch, mismatch_model="position-dependent")
        _load_params(tmp_path, monkeypatch)

        assert getattr(parameter, "mismatch_model", None) is None
        assert resolve_mismatch_model() == UNIFORM

    def test_declared_in_the_schema_with_the_values_the_code_accepts(self):
        from neoswga.core.mismatch_model import MISMATCH_MODELS

        declared = _schema()["properties"]["mismatch_model"]

        assert tuple(declared["enum"]) == MISMATCH_MODELS

    def test_the_optimizer_resolves_it_once_at_construction(self):
        """Like the penalty beside it, and for the same reason: re-reading the
        global per panel evaluation let two evaluations in one run weigh their
        sites differently. Source-pinned, because the behaviour it guards
        against is invisible unless something reassigns the global mid-run."""
        import ast
        import inspect

        from neoswga.core.base_optimizer import BaseOptimizer

        source = inspect.getsource(BaseOptimizer.__init__)
        tree = ast.parse(source.strip())
        assigned = {
            target.attr
            for node in ast.walk(tree)
            if isinstance(node, ast.Assign)
            for target in node.targets
            if isinstance(target, ast.Attribute)
        }

        assert {"mismatch_penalty", "mismatch_model"} <= assigned


class TestBackgroundAggregate:
    """What "the background" is when several hosts are listed. Phase 7a and 7b.

    params.json is the only route, as for the key above and for the same
    reason: `worst-case` is unmeasured as a default, and the key is the
    documented way in. No flag was added on `filter` when Phase 7b wired the
    gate there, deliberately, so there is no precedence between a flag and the
    file to preserve.
    """

    def test_the_value_in_params_json_reaches_the_parameter_module(self, tmp_path, monkeypatch):
        _load_params(tmp_path, monkeypatch, background_aggregate="worst-case")

        assert getattr(parameter, "background_aggregate", None) == "worst-case"

    def test_it_reaches_the_resolver_every_consumer_asks(self, tmp_path, monkeypatch):
        from neoswga.core.background_aggregate import WORST_CASE, resolve_background_aggregate

        _load_params(tmp_path, monkeypatch, background_aggregate="worst-case")

        assert resolve_background_aggregate() == WORST_CASE

    def test_absent_leaves_the_pooled_background(self, tmp_path, monkeypatch):
        """Loaded after a run that set it, so a leaked global would show."""
        from neoswga.core.background_aggregate import SUM, resolve_background_aggregate

        _load_params(tmp_path, monkeypatch, background_aggregate="worst-case")
        _load_params(tmp_path, monkeypatch)

        assert getattr(parameter, "background_aggregate", None) is None
        assert resolve_background_aggregate() == SUM

    def test_declared_in_the_schema_with_the_values_the_code_accepts(self):
        from neoswga.core.background_aggregate import BACKGROUND_AGGREGATES

        declared = _schema()["properties"]["background_aggregate"]

        assert tuple(declared["enum"]) == BACKGROUND_AGGREGATES
        assert declared["default"] == "sum"

    def test_the_value_from_params_json_reaches_the_step2_gate(self, tmp_path, monkeypatch):
        """The route that matters as of Phase 7b: `filter.get_all_rates` with
        no explicit argument takes the value the file set."""
        from neoswga.core import filter as filter_module
        from neoswga.core.background_aggregate import WORST_CASE

        _load_params(tmp_path, monkeypatch, background_aggregate="worst-case")
        monkeypatch.setattr(parameter, "use_bloom_filter", False, raising=False)
        monkeypatch.setattr(parameter, "bloom_filter_path", None, raising=False)
        monkeypatch.setattr(parameter, "min_fg_freq", 1e-9, raising=False)
        monkeypatch.setattr(parameter, "max_bg_freq", 5e-6, raising=False)
        monkeypatch.setattr(
            filter_module,
            "_get_rate_for_one_file",
            lambda task: {primer: 1 for primer in task[0]},
        )

        frame = filter_module.get_all_rates(
            ["AAAACCCCGG"],
            fg_prefixes=["fg"],
            bg_prefixes=["small", "large"],
            fg_total_length=100_000,
            bg_total_length=2_040_000,
            bg_seq_lengths=[40_000, 2_000_000],
        )

        assert frame.attrs["background_aggregate"] == WORST_CASE
        # 2 sites over 2,040,000 bp passes the pooled gate; 1 site over the
        # small host's own 40,000 bp does not, so the verdict is the key's.
        assert int(frame.loc[0, "bg_count"]) == 2
        assert bool(frame.loc[0, "bg_bool"]) is False


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
