"""A pan-target design should be able to require coverage on every target.

Item 5 of `docs/validation/getting_ahead_on_spacing_2026-09-18.md`, finished by
Phase 6c of the 2026-10-02 diversity plan.

On aggregate coverage a panel covering two targets 0.9 and 0.1 beats one
covering them 0.5 and 0.5. Multi-target designs are a capability none of the
three published SWGA tools has, so this extends a lead rather than closing a
gap.

**It was checked and reported, and deliberately NOT repaired, and that has
changed.** The stated reason was mechanical: the repair works through
`PoolObjective`, whose evaluator is `compute_metrics`, and `compute_metrics` did
not populate `per_target_coverage` -- it was filled in afterwards by
`unified_optimizer` so that all methods got it uniformly. Chasing the floor
through the repair would have scored every candidate panel against an empty
dict.

`compute_metrics` now produces the per-target dict and a scalar
`worst_target_coverage` from the same loop that produces `fg_coverage`
(`core/target_coverage.py`), which removes that reason. So
`min_per_target_coverage` is a row in `pool_objective._LIMITS` read against
`worst_target_coverage`, a field on `PoolConstraints`, and a member of
`panel_acceptance.LIMIT_KEYS`: one bounded repair is attempted and its outcome
is reported. The measured loss it acts on is in
`docs/validation/2026-10-02-diversity-baseline.md` (Q2, Q3).

The key is still unset by default and 0.0 still means disabled, so no delivered
panel moves unless someone sets it. Two guards keep it from passing on an
absence: an empty per-target dict with a floor requested is refused, and a
`worst_target_coverage` of `None` -- which is what an unmeasured target
produces -- misses the limit rather than clearing it.
"""

import json
import pathlib

import pytest

from neoswga.core.panel_acceptance import check_per_target_coverage

ROOT = pathlib.Path(__file__).resolve().parents[1]


class _Metrics:
    def __init__(self, per_target_coverage=None, worst_target_coverage="derive"):
        self.per_target_coverage = {} if per_target_coverage is None else dict(per_target_coverage)
        if worst_target_coverage == "derive":
            worst_target_coverage = (
                min(self.per_target_coverage.values()) if self.per_target_coverage else None
            )
        self.worst_target_coverage = worst_target_coverage


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

    def test_it_says_what_the_floor_is_read_against(self):
        """A figure reported without saying how it was produced reads as a
        property of the panel. This one is a minimum over several."""
        report = check_per_target_coverage(_Metrics({"a": 0.9, "b": 0.1}), 0.5)
        lines = "\n".join(report.lines()).lower()

        assert "worst_target_coverage" in lines
        assert "minimum" in lines
        assert "could not be measured" in lines

    def test_it_no_longer_claims_the_floor_is_not_repaired(self):
        """The claim was true until `compute_metrics` populated the figure, and
        a stale "nothing can be done" sends the reader to the wrong lever."""
        report = check_per_target_coverage(_Metrics({"a": 0.1}), 0.5)
        lines = "\n".join(report.lines()).lower()

        assert "not repaired" not in lines
        assert "bounded repair" in lines

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


# ----------------------------------------------------------------------
# Phase 6c: the figure comes from the evaluator, and the floor is enforced
# ----------------------------------------------------------------------

TARGET_LENGTH = 60_000
REACH = 3_000

#: Dense enough to cover one target end to end at `REACH`.
DENSE = sorted(set(range(REACH, TARGET_LENGTH, 2 * REACH - 500)) | {TARGET_LENGTH - REACH})
#: One window only, so a target reached by this alone sits well under 1.0.
SPARSE = [TARGET_LENGTH // 2]

#: Picked so no pair shares a complementary run the default screen rejects,
#: which would stop the repair using them before the floor could bind.
A_PRIMERS = ["AAACAAAGAAAT", "AAAGAAATAAAC"]
B_PRIMER = "GGGTGGGAGGGA"


class _TwoTargetCache:
    """`A_PRIMERS` bind densely in target A only; `B_PRIMER` sparsely in B."""

    def __init__(self, binds):
        self._binds = binds

    def get_positions(self, prefix, primer, strand="both"):
        import numpy as np

        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.array(self._binds.get((prefix, primer), ()), dtype=np.int64)


def _probe_class():
    from neoswga.core.base_optimizer import BaseOptimizer

    class _Probe(BaseOptimizer):
        @property
        def name(self):
            return "probe"

        @property
        def description(self):
            return "test double"

        def optimize(self, candidates, target_size=None, **kwargs):  # pragma: no cover
            raise NotImplementedError

    return _Probe


def _two_target_optimizer(seq_lengths=(TARGET_LENGTH, TARGET_LENGTH)):
    from neoswga.core.base_optimizer import OptimizerConfig

    binds = {("A", primer): DENSE for primer in A_PRIMERS}
    binds[("B", B_PRIMER)] = SPARSE

    return _probe_class()(
        position_cache=_TwoTargetCache(binds),
        fg_prefixes=["A", "B"],
        fg_seq_lengths=list(seq_lengths),
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(
            target_set_size=2,
            extension_reach=REACH,
            # Explicit, because without reaction conditions the default
            # effective coverage is None and every panel would carry a
            # "coverage unavailable" violation as well.
            coverage_metric="raw",
            max_dimer_bp=3,
            max_self_dimer_bp=4,
        ),
    )


class TestTheEvaluatorProducesTheFigure:
    """The change that makes the rest of this file possible: the quantity the
    floor is read against comes from the code that measures coverage."""

    def test_compute_metrics_populates_every_target(self):
        metrics = _two_target_optimizer().compute_metrics(A_PRIMERS)

        assert set(metrics.per_target_coverage) == {"A", "B"}
        assert metrics.per_target_coverage["A"] == pytest.approx(1.0, abs=0.01)
        assert metrics.per_target_coverage["B"] == pytest.approx(0.0)

    def test_the_worst_target_is_the_minimum_over_them(self):
        metrics = _two_target_optimizer().compute_metrics([*A_PRIMERS, B_PRIMER])
        per_target = metrics.per_target_coverage

        assert metrics.worst_target_coverage == pytest.approx(min(per_target.values()))
        assert metrics.worst_target_coverage == pytest.approx(per_target["B"])

    def test_the_aggregate_is_the_length_weighted_mean_of_them(self):
        """`fg_coverage` is not redefined: the two come from one loop, and a
        discrepancy would mean the floor is read against a different
        measurement from the one the headline reports."""
        metrics = _two_target_optimizer().compute_metrics([*A_PRIMERS, B_PRIMER])
        expected = sum(metrics.per_target_coverage.values()) / 2

        assert metrics.fg_coverage == pytest.approx(expected)

    def test_the_focused_evaluator_agrees_on_the_worst_target(self):
        """`objective_for_optimizer` hands the search `compute_pool_metrics`
        when no gap or host limit is set, which is the case here."""
        optimizer = _two_target_optimizer()
        panel = [*A_PRIMERS, B_PRIMER]

        full = optimizer.compute_metrics(panel)
        light = optimizer.compute_pool_metrics(panel)

        assert light.worst_target_coverage == pytest.approx(full.worst_target_coverage)

    def test_a_single_target_run_reports_its_one_target(self):
        from neoswga.core.base_optimizer import OptimizerConfig

        optimizer = _probe_class()(
            position_cache=_TwoTargetCache({("A", primer): DENSE for primer in A_PRIMERS}),
            fg_prefixes=["A"],
            fg_seq_lengths=[TARGET_LENGTH],
            bg_prefixes=[],
            bg_seq_lengths=[],
            config=OptimizerConfig(target_set_size=2, extension_reach=REACH),
        )
        metrics = optimizer.compute_metrics(A_PRIMERS)

        assert metrics.worst_target_coverage == pytest.approx(metrics.fg_coverage)

    def test_the_summary_json_carries_it(self):
        """`step4_improved_df_summary.json` is what the report reads."""
        metrics = _two_target_optimizer().compute_metrics([*A_PRIMERS, B_PRIMER])
        as_dict = metrics.to_dict()

        assert as_dict["worst_target_coverage"] == pytest.approx(
            min(as_dict["per_target_coverage"].values())
        )
        json.dumps(as_dict["worst_target_coverage"])


class TestAnUnmeasuredTargetMakesTheWorstUnknown:
    """A worst case over the members that answered is not a worst case."""

    def test_a_target_with_no_usable_length_is_absent_rather_than_zero(self):
        metrics = _two_target_optimizer(seq_lengths=(TARGET_LENGTH, 0)).compute_metrics(A_PRIMERS)

        assert set(metrics.per_target_coverage) == {"A"}

    def test_the_worst_is_none_and_not_the_lowest_of_the_rest(self):
        metrics = _two_target_optimizer(seq_lengths=(TARGET_LENGTH, 0)).compute_metrics(A_PRIMERS)

        assert metrics.worst_target_coverage is None

    def test_a_floor_on_an_unknown_worst_is_missed_and_not_cleared(self):
        """`_unmet` is where the refusal happens. A `None` that passed would
        report compliance for a measurement nobody made."""
        from neoswga.core.pool_objective import _AT_LEAST, _unmet

        assert _unmet(None, 0.5, _AT_LEAST) is True

    def test_the_objective_names_it_as_a_violation(self):
        from neoswga.core.panel_refinement import objective_for_optimizer
        from neoswga.core.pool_objective import PoolConstraints

        optimizer = _two_target_optimizer(seq_lengths=(TARGET_LENGTH, 0))
        objective = objective_for_optimizer(
            optimizer, PoolConstraints(min_per_target_coverage=0.5, coverage_metric="raw")
        )

        assert "per-target coverage below minimum" in objective.violations(A_PRIMERS)

    def test_the_limits_table_shows_it_as_not_measured(self):
        from neoswga.core.panel_acceptance import enforce_constraints
        from neoswga.core.pool_objective import PoolConstraints

        optimizer = _two_target_optimizer(seq_lengths=(TARGET_LENGTH, 0))
        report = enforce_constraints(
            A_PRIMERS,
            optimizer,
            candidates=[*A_PRIMERS, B_PRIMER],
            constraints=PoolConstraints(min_per_target_coverage=0.5, coverage_metric="raw"),
        )
        lines = "\n".join(report.lines())

        assert "not measured" in lines
        assert "min_per_target_coverage" in lines

    def test_the_limits_table_says_the_figure_is_a_minimum(self):
        from neoswga.core.panel_acceptance import enforce_constraints
        from neoswga.core.pool_objective import PoolConstraints

        optimizer = _two_target_optimizer()
        report = enforce_constraints(
            A_PRIMERS,
            optimizer,
            candidates=list(A_PRIMERS),
            constraints=PoolConstraints(min_per_target_coverage=0.5, coverage_metric="raw"),
        )
        lines = "\n".join(report.lines())

        assert "worst_target_coverage" in lines
        assert "MINIMUM" in lines


class TestAnEmptyPerTargetDictCannotPassVacuously:
    """`docs/development/design-contracts.md` records this defect:
    `check_per_target_coverage` and `result_validation` both gate on a
    non-empty dict, so an empty one made a REQUESTED floor pass without a
    word."""

    def test_a_requested_floor_with_nothing_measured_is_refused(self):
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.unified_optimizer import _require_per_target

        with pytest.raises(ModelEvaluationError):
            _require_per_target({}, None, ["A", "B"], 0.5)

    def test_an_unmeasured_target_with_a_requested_floor_is_refused(self):
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.unified_optimizer import _require_per_target

        with pytest.raises(ModelEvaluationError):
            _require_per_target({"A": 0.9}, None, ["A", "B"], 0.5)

    def test_the_refusal_names_the_unmeasured_target(self):
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.unified_optimizer import _require_per_target

        with pytest.raises(ModelEvaluationError) as excinfo:
            _require_per_target({"A": 0.9}, None, ["A", "B"], 0.5)

        assert "B" in str(excinfo.value)

    def test_no_floor_requested_refuses_nothing(self):
        """What keeps this step from moving a delivered panel: with the key
        unset, an absent per-target measurement is not an error."""
        from neoswga.core.unified_optimizer import _require_per_target

        assert _require_per_target({}, None, ["A", "B"], None) is None
        assert _require_per_target({}, None, ["A", "B"], 0.0) is None

    def test_a_measured_worst_passes(self):
        from neoswga.core.unified_optimizer import _require_per_target

        assert _require_per_target({"A": 0.9, "B": 0.1}, 0.1, ["A", "B"], 0.5) is None

    def test_the_reporter_still_says_nothing_on_an_empty_dict(self):
        """It is a reporter, and "nothing to say" and "all clear" print the
        same. The refusal above is what makes that safe."""
        assert check_per_target_coverage(_Metrics({}), 0.5) is None


class TestTheThreeNameSetsAgree:
    """One tuple, one dataclass and one limit table. A key declared in one and
    spelled differently in another enforces the wrong limit silently."""

    def test_every_limit_key_has_a_row_in_the_limit_table(self):
        from neoswga.core.panel_acceptance import LIMIT_KEYS
        from neoswga.core.pool_objective import _LIMITS

        assert {field for field, *_rest in _LIMITS} == set(LIMIT_KEYS)

    def test_every_limit_key_is_a_field_of_the_constraints(self):
        from neoswga.core.panel_acceptance import LIMIT_KEYS
        from neoswga.core.pool_objective import PoolConstraints

        assert set(LIMIT_KEYS) <= set(PoolConstraints().__dataclass_fields__)

    def test_the_configured_keys_are_the_limits_plus_the_dimer_floor(self):
        from neoswga.core.panel_acceptance import CONFIGURED_LIMIT_KEYS, LIMIT_KEYS

        assert set(CONFIGURED_LIMIT_KEYS) - set(LIMIT_KEYS) == {"max_dimer_dg"}

    def test_the_per_target_floor_sits_where_a_repaired_limit_belongs(self):
        """It was in `CONFIGURED_LIMIT_KEYS` only, which is the configuration
        plumbing. Enforcement reads `LIMIT_KEYS`."""
        from neoswga.core.panel_acceptance import LIMIT_KEYS

        assert "min_per_target_coverage" in LIMIT_KEYS

    def test_it_is_read_against_the_worst_target(self):
        from neoswga.core.pool_objective import _AT_LEAST, _LIMITS

        row = next(r for r in _LIMITS if r[0] == "min_per_target_coverage")

        assert row[1] == "worst_target_coverage"
        assert row[2] == _AT_LEAST
        assert row[4] is False, "it is not measured against the background"

    def test_the_first_two_rows_keep_their_position(self):
        """`pool_planner._repair` records `reasons[0]`, so the new row was
        appended rather than inserted."""
        from neoswga.core.pool_objective import _LIMITS

        assert [row[0] for row in _LIMITS][:2] == [
            "min_selectivity_density",
            "max_background_sites",
        ]


class TestNothingMovesWithTheKeyUnset:
    def test_an_unset_key_builds_no_constraints_at_all(self):
        """Building an objective is what would change a delivered panel."""
        from neoswga.core.panel_acceptance import constraints_from_parameter

        class _Nothing:
            pass

        assert constraints_from_parameter(_Nothing()) is None
        assert constraints_from_parameter(_Nothing(), {}) is None

    def test_a_zero_floor_still_builds_no_constraints(self):
        """0.0 is how this key has always spelled "disabled"; reading it as a
        floor of zero would give a run that asked for nothing an objective."""
        from neoswga.core.panel_acceptance import constraints_from_parameter

        class _Zero:
            min_per_target_coverage = 0.0

        assert constraints_from_parameter(_Zero()) is None
        assert constraints_from_parameter(_Zero(), {"min_per_target_coverage": 0.0}) is None

    def test_a_set_floor_does_build_them(self):
        from neoswga.core.panel_acceptance import constraints_from_parameter

        class _Set:
            min_per_target_coverage = 0.4

        constraints = constraints_from_parameter(_Set())

        assert constraints is not None
        assert constraints.min_per_target_coverage == 0.4


class TestTheFlagReachesTheEnforcingCode:
    """Known Issue 8's shape. The flag printed a per-target table; if it did
    not also reach `constraints_from_parameter` it would not be repaired, while
    the same value in params.json would be."""

    def test_the_flag_builds_the_constraints(self):
        from neoswga.core.panel_acceptance import constraints_from_parameter

        class _Nothing:
            pass

        constraints = constraints_from_parameter(_Nothing(), {"min_per_target_coverage": 0.7})

        assert constraints is not None
        assert constraints.min_per_target_coverage == 0.7

    def test_the_flag_beats_the_file(self):
        from neoswga.core.panel_acceptance import constraints_from_parameter

        class _File:
            min_per_target_coverage = 0.3

        constraints = constraints_from_parameter(_File(), {"min_per_target_coverage": 0.7})

        assert constraints.min_per_target_coverage == 0.7

    def test_an_absent_flag_does_not_beat_the_file(self):
        from neoswga.core.panel_acceptance import constraints_from_parameter

        class _File:
            min_per_target_coverage = 0.3

        constraints = constraints_from_parameter(_File(), {"min_per_target_coverage": None})

        assert constraints.min_per_target_coverage == 0.3

    def test_the_run_hands_its_flags_to_the_acceptance_step(self):
        import inspect

        from neoswga.core import unified_optimizer

        source = inspect.getsource(unified_optimizer.run_optimization)

        assert "flags=kwargs" in source


class TestTheFloorTriggersOneBoundedRepair:
    """The plan's verification for 6c: a floor above the delivered worst target
    triggers exactly one repair and the outcome is reported."""

    FLOOR = 0.05

    def _apply(self, monkeypatch, panel, candidates, floor=None):
        from neoswga.core import pool_planner
        from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus
        from neoswga.core.panel_acceptance import apply_configured_limits
        from neoswga.core.pool_objective import PoolConstraints

        optimizer = _two_target_optimizer()
        calls = []
        original = pool_planner.repair_panel

        def counted(*args, **kwargs):
            calls.append(list(args[3] if len(args) > 3 else kwargs.get("reasons") or ()))
            return original(*args, **kwargs)

        monkeypatch.setattr(pool_planner, "repair_panel", counted)

        result = OptimizationResult(
            primers=list(panel),
            score=1.0,
            status=OptimizationStatus.SUCCESS,
            metrics=optimizer.compute_metrics(list(panel)),
            iterations=1,
            optimizer_name="probe",
        )
        replacement, report = apply_configured_limits(
            result,
            optimizer,
            candidates=list(candidates),
            config=optimizer.config,
            constraints=PoolConstraints(
                min_per_target_coverage=self.FLOOR if floor is None else floor,
                coverage_metric="raw",
            ),
        )
        return replacement, report, calls

    def test_the_delivered_panel_starves_one_target(self):
        """Guard the guard: without this the repair below proves nothing."""
        metrics = _two_target_optimizer().compute_metrics(A_PRIMERS)

        assert metrics.worst_target_coverage == pytest.approx(0.0)
        assert metrics.worst_target_coverage < self.FLOOR

    def test_a_floor_above_the_worst_target_triggers_exactly_one_repair(self, monkeypatch):
        _replacement, _report, calls = self._apply(monkeypatch, A_PRIMERS, [*A_PRIMERS, B_PRIMER])

        assert len(calls) == 1, f"expected one bounded repair, got {len(calls)}"

    def test_the_repair_is_asked_about_the_per_target_limit(self, monkeypatch):
        _replacement, _report, calls = self._apply(monkeypatch, A_PRIMERS, [*A_PRIMERS, B_PRIMER])

        assert "per-target coverage below minimum" in calls[0]

    def test_the_repair_reaches_the_starved_target(self, monkeypatch):
        replacement, report, _calls = self._apply(monkeypatch, A_PRIMERS, [*A_PRIMERS, B_PRIMER])

        assert report.repaired is True
        assert B_PRIMER in report.primers
        assert report.violations == ()
        assert replacement is not None, "a changed panel must replace the result"

    def test_the_outcome_is_reported_when_the_repair_cannot_succeed(self, monkeypatch):
        """A pool with nothing that reaches target B cannot meet the floor, and
        the report has to say so rather than return a bare None."""
        _replacement, report, calls = self._apply(monkeypatch, A_PRIMERS, list(A_PRIMERS))

        assert len(calls) == 1
        assert report.repaired is False
        assert "per-target coverage below minimum" in report.violations
        assert any(name == "min_per_target_coverage" for name, _v, _l in report.values)

    def test_a_panel_already_above_the_floor_is_not_repaired(self, monkeypatch):
        _replacement, report, _calls = self._apply(
            monkeypatch, [*A_PRIMERS, B_PRIMER], [*A_PRIMERS, B_PRIMER], floor=0.001
        )

        assert report.violations == ()
        assert report.repaired is False
