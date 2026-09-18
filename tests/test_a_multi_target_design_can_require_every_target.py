"""A pan-target design should be able to require coverage on every target.

Item 5 of `docs/validation/getting_ahead_on_spacing_2026-09-18.md`.

`per_target_coverage` is computed in `run_optimization` after the search and
consumed only by `OptimizationResult.validate`, which emits a warning. No
selection stage sees it and no limit enforces it, so on aggregate coverage a
panel covering two targets 0.9 and 0.1 beats one covering them 0.5 and 0.5, and
the run says so only in passing.

Multi-target designs are a capability none of the three published SWGA tools
has, so this extends a lead rather than closing a gap.

**It is checked and reported, and deliberately NOT repaired.** The repair works
through `PoolObjective`, whose evaluator is `compute_metrics`, and
`compute_metrics` does not populate `per_target_coverage` -- it is filled in by
the caller so that all methods get it uniformly. Chasing this floor through the
repair would therefore score every candidate panel against an empty dict. That
is the same reason `pool_planner.repair_panel` refuses to repair a dimer
violation, and it is stated rather than hidden.
"""

import json
import pathlib

import pytest

from neoswga.core.panel_acceptance import check_per_target_coverage

ROOT = pathlib.Path(__file__).resolve().parents[1]


class _Metrics:
    def __init__(self, per_target_coverage=None):
        self.per_target_coverage = {} if per_target_coverage is None else dict(per_target_coverage)


class TestItIsOffByDefault:
    def test_the_schema_declares_it(self):
        schema = json.loads(
            (ROOT / "neoswga" / "core" / "schema" / "params.schema.json").read_text()
        )

        assert "min_per_target_coverage" in schema["properties"]

    def test_the_parameter_global_is_unset(self):
        from neoswga.core import parameter

        assert parameter.min_per_target_coverage is None

    def test_the_flag_default_is_the_none_sentinel(self):
        """A real default of 0.0 would beat a configured value on every run,
        which is Known Issue 8's shape. The flag carried 0.0 before this.
        """
        import argparse

        from neoswga.cli._optimize_parser import _add_optimize_selection_groups

        parser = argparse.ArgumentParser()
        _add_optimize_selection_groups(parser)
        actions = {a.dest: a for a in parser._actions}

        assert actions["min_per_target_coverage"].default is None

    def test_no_floor_yields_no_check(self):
        assert check_per_target_coverage(_Metrics({"a": 0.1}), None) is None

    def test_a_zero_floor_yields_no_check(self):
        """0.0 is how the flag has always spelled "disabled"."""
        assert check_per_target_coverage(_Metrics({"a": 0.0}), 0.0) is None

    def test_a_single_genome_run_yields_no_check(self):
        """`per_target_coverage` is empty in single-genome mode, and a floor on
        an absent measurement must not read as satisfied."""
        assert check_per_target_coverage(_Metrics({}), 0.5) is None


class TestTheFloorBinds:
    def test_a_target_below_the_floor_is_named(self):
        report = check_per_target_coverage(_Metrics({"a": 0.9, "b": 0.1}), 0.5)

        assert report.below == ("b",)
        assert report.met is False

    def test_every_target_above_the_floor_passes(self):
        report = check_per_target_coverage(_Metrics({"a": 0.9, "b": 0.7}), 0.5)

        assert report.below == ()
        assert report.met is True

    def test_the_worst_target_is_identified(self):
        report = check_per_target_coverage(_Metrics({"a": 0.9, "b": 0.4, "c": 0.2}), 0.5)

        assert report.worst_target == "c"
        assert report.worst_coverage == pytest.approx(0.2)

    def test_a_target_exactly_on_the_floor_passes(self):
        report = check_per_target_coverage(_Metrics({"a": 0.5}), 0.5)

        assert report.met is True

    def test_the_aggregate_cannot_hide_a_starved_target(self):
        """The defect this exists for: two panels with the same mean, one of
        which abandons a target."""
        balanced = check_per_target_coverage(_Metrics({"a": 0.5, "b": 0.5}), 0.4)
        lopsided = check_per_target_coverage(_Metrics({"a": 0.9, "b": 0.1}), 0.4)

        assert balanced.met is True
        assert lopsided.met is False


class TestTheReport:
    def test_it_names_every_target_and_its_coverage(self):
        report = check_per_target_coverage(_Metrics({"alpha": 0.9, "beta": 0.1}), 0.5)
        lines = "\n".join(report.lines())

        assert "alpha" in lines
        assert "beta" in lines
        assert "0.1" in lines

    def test_it_says_the_floor_is_not_repaired_and_why(self):
        """Silence here would read as "nothing could be done", when the truth
        is that the repair cannot see this quantity at all."""
        report = check_per_target_coverage(_Metrics({"a": 0.1}), 0.5)
        lines = "\n".join(report.lines()).lower()

        assert "not repaired" in lines

    def test_a_passing_report_still_shows_the_numbers(self):
        report = check_per_target_coverage(_Metrics({"a": 0.9, "b": 0.7}), 0.5)

        assert report.lines()


class TestItReachesTheOptimizePath:
    def test_the_acceptance_module_is_called_from_the_optimize_path(self):
        import ast
        import inspect

        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "check_per_target_coverage" in names

    def test_the_floor_resolves_from_params_json_when_no_flag_is_given(self):
        """With the flag at the None sentinel, the configured value must win."""
        from neoswga.core.panel_acceptance import per_target_floor

        class _Args:
            min_per_target_coverage = None

        class _Params:
            min_per_target_coverage = 0.6

        assert per_target_floor(_Args(), _Params()) == 0.6

    def test_an_explicit_flag_beats_the_configured_value(self):
        from neoswga.core.panel_acceptance import per_target_floor

        class _Args:
            min_per_target_coverage = 0.2

        class _Params:
            min_per_target_coverage = 0.6

        assert per_target_floor(_Args(), _Params()) == 0.2

    def test_neither_set_is_off(self):
        from neoswga.core.panel_acceptance import per_target_floor

        class _Empty:
            pass

        assert per_target_floor(_Empty(), _Empty()) is None
