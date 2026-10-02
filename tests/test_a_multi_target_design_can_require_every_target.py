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


def _two_target_optimizer(
    seq_lengths=(TARGET_LENGTH, TARGET_LENGTH), floor=0.05, geometry="symmetric"
):
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
            # A floor in force, which is what makes `compute_metrics` produce
            # the per-target figure at all: with none configured nothing reads
            # it during the search and it is not computed, which
            # `TestTheFigureIsComputedOnlyWhenSomethingReadsIt` pins.
            min_per_target_coverage=floor,
            coverage_geometry=geometry,
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

    def test_the_aggregate_matches_them_with_no_record_joins(self):
        """With no record starts the two geometries coincide, so the per-target
        figures and `fg_coverage` agree -- which covers every single-record
        reference in this repository. Where they diverge is the record join, and
        `TestTheFloorReadsTheConfinedFigure` is where that is pinned."""
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
            config=OptimizerConfig(
                target_set_size=2, extension_reach=REACH, min_per_target_coverage=0.05
            ),
        )
        metrics = optimizer.compute_metrics(A_PRIMERS)

        assert metrics.worst_target_coverage == pytest.approx(metrics.fg_coverage)

    def test_the_html_report_shows_it_and_says_it_is_a_minimum(self):
        """Review finding F6: the report rendered the per-target dict and not
        the figure the limit is read against, leaving a reader to take the
        minimum themselves. A reduction over several measurements has to say
        that it is one."""
        from neoswga.core.report.technical_report import _render_per_target_coverage

        html = _render_per_target_coverage({"a": 0.9, "b": 0.42}, 0.42)

        assert "worst_target_coverage" in html
        assert "42.0%" in html
        assert "MINIMUM" in html
        assert "confined" in html

    def test_the_html_report_says_not_measured_rather_than_a_number(self):
        """0.0 there would read as a starved target rather than as an absent
        measurement, which is the defect this whole step is built around."""
        from neoswga.core.report.technical_report import _render_per_target_coverage

        html = _render_per_target_coverage({"a": 0.9, "b": 0.42}, None)
        worst_row = html.split("worst_target_coverage</strong>")[1].split("</tr>")[0]

        assert "not measured" in worst_row
        assert "%" not in worst_row

    def test_the_report_metrics_do_not_default_the_figure_to_zero(self):
        """A summary written before the field existed must not read as a panel
        that reaches one target nowhere."""
        from neoswga.core.report.metrics import CoverageMetrics

        assert CoverageMetrics().worst_target_coverage is None

    def test_the_module_is_watched_for_unreachability(self):
        """Review finding F8. Two of its four public names are reached from one
        call site each, behind a limit that is off by default, so a default run
        would notice nothing if either became orphaned."""
        from tests.test_no_capability_is_unreachable import WATCHED

        assert "core/target_coverage.py" in WATCHED

    def test_the_summary_json_carries_it(self):
        """`step4_improved_df_summary.json` is what the report reads."""
        metrics = _two_target_optimizer().compute_metrics([*A_PRIMERS, B_PRIMER])
        as_dict = metrics.to_dict()

        assert as_dict["worst_target_coverage"] == pytest.approx(
            min(as_dict["per_target_coverage"].values())
        )
        json.dumps(as_dict["worst_target_coverage"])


#: A draft assembly: 500 contigs of 4 kb, one site 100 bp inside each contig
#: end, so every window would cross a join if nothing stopped it.
DRAFT_LENGTH = 2_000_000
DRAFT_CONTIG = DRAFT_LENGTH // 500
DRAFT_STARTS = [i * DRAFT_CONTIG for i in range(500)]
DRAFT_SITES = [start + DRAFT_CONTIG - 100 for start in DRAFT_STARTS]


class _DraftCache:
    """One prefix of 500 contigs, with its record starts, as a current
    `neoswga filter` index carries them."""

    def __init__(self, record_starts=DRAFT_STARTS, sites=DRAFT_SITES):
        self._starts = record_starts
        self._sites = sites

    def get_positions(self, prefix, primer, strand="both"):
        import numpy as np

        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.array(self._sites, dtype=np.int64)

    def get_record_starts(self, prefix):
        return None if self._starts is None else list(self._starts)


def _draft_optimizer(cache=None, floor=0.9, geometry="symmetric"):
    from neoswga.core.base_optimizer import OptimizerConfig

    return _probe_class()(
        position_cache=cache or _DraftCache(),
        fg_prefixes=["draft"],
        fg_seq_lengths=[DRAFT_LENGTH],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(
            target_set_size=1,
            extension_reach=REACH,
            coverage_metric="raw",
            min_per_target_coverage=floor,
            coverage_geometry=geometry,
        ),
    )


class TestTheFloorReadsTheConfinedFigure:
    """Review finding F2, and the reason the first version of Phase 6c was
    wrong. It took the per-target figure from the loop that produces
    `fg_coverage`, which is `coverage.union_coverage` and does not confine a
    window to the record holding its site. Audit finding F3 of
    `docs/validation/additives_coverage_audit_2026-09-14/` records that as a P1
    error for target coverage on a fragmented reference, fixed and pinned by
    `tests/test_coverage_windows_stop_at_record_edges.py`.

    On this 500-contig draft the two read 0.775 and 0.9996, so a floor of 0.9
    passed a panel that does not meet it -- a false pass on 449,100 bp credited
    across joins the polymerase cannot cross.
    """

    def test_the_unconfined_figure_would_be_a_false_pass(self):
        """Guard the guard: without this spread the rest proves nothing."""
        metrics = _draft_optimizer().compute_metrics([A_PRIMERS[0]])

        assert metrics.fg_coverage == pytest.approx(0.9995, abs=1e-4)
        assert metrics.fg_coverage > 0.9, "the unconfined figure must clear the floor"

    def test_the_per_target_figure_is_the_confined_one(self):
        metrics = _draft_optimizer().compute_metrics([A_PRIMERS[0]])

        assert metrics.per_target_coverage["draft"] == pytest.approx(0.775)
        assert metrics.worst_target_coverage == pytest.approx(0.775)

    def test_a_floor_of_zero_point_nine_fails_on_this_draft(self):
        """The false pass, as a verdict."""
        from neoswga.core.panel_refinement import objective_for_optimizer
        from neoswga.core.pool_objective import PoolConstraints

        optimizer = _draft_optimizer()
        objective = objective_for_optimizer(
            optimizer, PoolConstraints(min_per_target_coverage=0.9, coverage_metric="raw")
        )

        assert "per-target coverage below minimum" in objective.violations([A_PRIMERS[0]])

    def test_the_focused_evaluator_fails_it_too(self):
        """The search reads `compute_pool_metrics`; a confined full evaluator
        and an unconfined light one would accept and refuse the same panel."""
        light = _draft_optimizer().compute_pool_metrics([A_PRIMERS[0]])

        assert light.worst_target_coverage == pytest.approx(0.775)

    def test_it_is_the_figure_evaluate_set_reports(self):
        """One quantity under one key name across the two commands. They used
        to disagree: `optimize` could accept a panel at 0.9996 and
        `evaluate-set` report 0.775 for the same panel and reference."""
        from neoswga.core.coverage import compute_per_prefix_coverage

        metrics = _draft_optimizer().compute_metrics([A_PRIMERS[0]])
        _aggregate, confined = compute_per_prefix_coverage(
            cache=_DraftCache(),
            primers=[A_PRIMERS[0]],
            prefixes=["draft"],
            seq_lengths=[DRAFT_LENGTH],
            extension=REACH,
        )

        assert metrics.per_target_coverage["draft"] == pytest.approx(confined["draft"])

    def test_an_index_without_record_starts_still_agrees_with_the_aggregate(self):
        """The prepared indexes in this repository store none, which is why no
        recorded figure moved either way."""
        optimizer = _draft_optimizer(_DraftCache(record_starts=None))
        metrics = optimizer.compute_metrics([A_PRIMERS[0]])

        assert metrics.per_target_coverage["draft"] == pytest.approx(metrics.fg_coverage)
        assert metrics.per_target_coverage["draft"] == pytest.approx(0.9995, abs=1e-4)

    def test_one_marking_implementation_serves_both_commands(self):
        """Three implementations of one coverage number is how this codebase
        produced disagreeing figures before; `union_coverage` says so."""
        import inspect

        from neoswga.core import coverage

        assert "confined_covered_bases" in inspect.getsource(coverage.compute_per_prefix_coverage)


class TestTheFigureIsComputedOnlyWhenSomethingReadsIt:
    """A bitmap per target is linear in the sequence -- measured at +40 ms, about
    18 percent of one objective evaluation, on a 145 Mb foreground -- and a
    search performs very many evaluations. Paying that so an unused figure could
    be reported is the wrong trade, so `compute_metrics` produces it only when a
    per-target floor is configured, and `attach_per_target` fills it in once on
    the delivered panel for every other run.
    """

    def test_with_no_floor_the_search_does_not_compute_it(self):
        metrics = _two_target_optimizer(floor=None).compute_metrics([*A_PRIMERS, B_PRIMER])

        assert metrics.per_target_coverage == {}
        assert metrics.worst_target_coverage is None

    def test_a_zero_floor_counts_as_no_floor(self):
        """0.0 has meant disabled for this key since before it was a limit, and
        the budget decision must read it the same way the limit does."""
        metrics = _two_target_optimizer(floor=0.0).compute_metrics([*A_PRIMERS, B_PRIMER])

        assert metrics.per_target_coverage == {}

    def test_with_a_floor_the_search_computes_it(self):
        metrics = _two_target_optimizer(floor=0.05).compute_metrics([*A_PRIMERS, B_PRIMER])

        assert set(metrics.per_target_coverage) == {"A", "B"}
        assert metrics.worst_target_coverage is not None

    def test_not_computed_is_absent_and_never_zero(self):
        """The whole point: an uncomputed figure must not look like a starved
        target, or a floor would read a fabricated number out of it."""
        from neoswga.core.target_coverage import not_computed

        absent = not_computed()

        assert absent.per_target == {}
        assert absent.worst is None

    def test_a_requested_floor_is_still_refused_on_an_uncomputed_figure(self):
        """Belt and braces: if the gate and the floor ever disagreed, the floor
        must refuse rather than pass."""
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.target_coverage import require_per_target

        with pytest.raises(ModelEvaluationError):
            require_per_target({}, None, ["A", "B"], 0.5)

    def test_the_light_evaluator_makes_the_same_decision(self):
        """The search reads `compute_pool_metrics`. One evaluator computing the
        figure and the other not is how two definitions of one number start."""
        optimizer = _two_target_optimizer(floor=None)

        assert optimizer.compute_metrics([*A_PRIMERS, B_PRIMER]).per_target_coverage == {}
        assert optimizer.compute_pool_metrics([*A_PRIMERS, B_PRIMER]).worst_target_coverage is None

    def test_one_reader_decides_whether_a_floor_is_in_force(self):
        from neoswga.core.base_optimizer import OptimizerConfig
        from neoswga.core.target_coverage import per_target_floor_in_force

        assert per_target_floor_in_force(OptimizerConfig()) is None
        assert per_target_floor_in_force(OptimizerConfig(min_per_target_coverage=0.4)) == 0.4

    def test_the_config_field_is_resolved_from_the_run(self):
        """`_build_optimizer_config` must set it, or the gate reads None on
        every run and the figure is never computed during a search."""
        import inspect

        from neoswga.core import unified_optimizer

        source = inspect.getsource(unified_optimizer._build_optimizer_config)

        assert 'min_per_target_coverage=pick("min_per_target_coverage", None)' in source


class TestTheReportStillCarriesTheFigureWithNoFloor:
    """A reported figure must not disappear because an unrelated limit is
    unset. `attach_per_target` fills it in once, on the delivered panel."""

    def _result(self, primers, metrics):
        from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus

        return OptimizationResult(
            primers=list(primers),
            score=1.0,
            status=OptimizationStatus.SUCCESS,
            metrics=metrics,
            iterations=1,
            optimizer_name="probe",
        )

    def _attach(self, optimizer, primers, geometry="symmetric"):
        from neoswga.core.target_coverage import attach_per_target

        result = self._result(primers, optimizer.compute_metrics(list(primers)))
        return attach_per_target(
            result,
            optimizer.cache,
            optimizer.fg_prefixes,
            optimizer.fg_seq_lengths,
            extension=REACH,
            circular=False,
            geometry=geometry,
        )

    def test_it_is_filled_in_once_when_the_search_did_not(self):
        optimizer = _two_target_optimizer(floor=None)
        attached = self._attach(optimizer, [*A_PRIMERS, B_PRIMER])

        assert set(attached.metrics.per_target_coverage) == {"A", "B"}
        assert attached.metrics.worst_target_coverage is not None

    def test_the_filled_in_figure_equals_the_searched_one(self):
        """The two paths must agree, or the number a floor was enforced against
        and the number reported would be different quantities."""
        searched = _two_target_optimizer(floor=0.05).compute_metrics([*A_PRIMERS, B_PRIMER])
        attached = self._attach(_two_target_optimizer(floor=None), [*A_PRIMERS, B_PRIMER])

        assert attached.metrics.per_target_coverage == pytest.approx(searched.per_target_coverage)
        assert attached.metrics.worst_target_coverage == pytest.approx(
            searched.worst_target_coverage
        )

    def test_it_never_overwrites_what_the_search_measured(self):
        """Asserted by OBJECT IDENTITY, not by comparing values: both paths
        compute the same confined quantity, so a recomputation would give equal
        numbers and an equality check would pass while the result the floor was
        enforced against had been replaced."""
        from neoswga.core.target_coverage import attach_per_target

        optimizer = _two_target_optimizer(floor=0.05)
        searched = self._result(
            [*A_PRIMERS, B_PRIMER], optimizer.compute_metrics([*A_PRIMERS, B_PRIMER])
        )
        attached = attach_per_target(
            searched,
            optimizer.cache,
            optimizer.fg_prefixes,
            optimizer.fg_seq_lengths,
            extension=REACH,
            circular=False,
            geometry="symmetric",
        )

        assert attached is searched, "a populated result must be returned untouched"

    def test_the_run_assigns_the_attached_result_back(self):
        """`attach_per_target(result, ...)` called and discarded would leave
        every no-floor run reporting nothing, and an AST check that the name is
        called cannot tell the two apart."""
        import inspect

        from neoswga.core import unified_optimizer

        source = inspect.getsource(unified_optimizer.run_optimization)

        assert "result = attach_per_target(" in source

    def test_an_unmeasurable_target_stays_absent_rather_than_zero(self):
        """`compute_per_prefix_coverage` writes 0.0 for a zero-length prefix,
        which is the convention this module exists to refuse."""
        optimizer = _two_target_optimizer(seq_lengths=(TARGET_LENGTH, 0), floor=None)
        attached = self._attach(optimizer, A_PRIMERS)

        assert "B" not in attached.metrics.per_target_coverage
        assert attached.metrics.worst_target_coverage is None

    def test_the_optimize_path_calls_it(self):
        import ast
        import inspect

        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        called = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "attach_per_target" in called


class TestTheDirectionalGeometryIsRefusedNotApproximated:
    """`coverage.confined_covered_bases` marks `[pos-reach, pos+reach)`, which is
    the symmetric geometry and what `evaluate-set` reports. Under
    `coverage_geometry: directional` a site reaches one way only and half the
    width, so reporting the symmetric number there is the same class of error as
    reporting an unconfined one -- the error this step was reopened to fix.
    """

    def test_the_figure_is_absent_under_directional(self):
        metrics = _two_target_optimizer(geometry="directional").compute_metrics(
            [*A_PRIMERS, B_PRIMER]
        )

        assert metrics.per_target_coverage == {}
        assert metrics.worst_target_coverage is None

    def test_the_reason_is_carried_rather_than_a_silent_absence(self):
        from neoswga.core.target_coverage import confined_coverage_by_target

        optimizer = _two_target_optimizer(geometry="directional")
        coverage = confined_coverage_by_target(optimizer, {"A": DENSE, "B": SPARSE})

        assert coverage.unavailable is not None
        assert "directional" in coverage.unavailable
        assert coverage.worst is None

    def test_a_reason_suppresses_the_worst_even_beside_measured_targets(self):
        """`not_computed()` always leaves `per_target` empty, so the empty-dict
        branch would hide a missing `unavailable` check. Constructed directly:
        a reason must win over figures that happen to be present, or a future
        caller reporting a partial measurement alongside a reason would publish
        a minimum over part of the panel."""
        from neoswga.core.target_coverage import TargetCoverage

        partial = TargetCoverage(0.5, {"a": 0.5, "b": 0.9}, (), "a reason")

        assert partial.worst is None
        assert partial.worst_target is None
        assert TargetCoverage(0.5, {"a": 0.5, "b": 0.9}, ()).worst == pytest.approx(0.5)

    def test_a_floor_is_refused_under_directional(self):
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.target_coverage import DIRECTIONAL_UNAVAILABLE, require_per_target

        with pytest.raises(ModelEvaluationError) as excinfo:
            require_per_target({}, None, ["A", "B"], 0.5, DIRECTIONAL_UNAVAILABLE)

        assert "directional" in str(excinfo.value)

    def test_the_objective_treats_it_as_unmet_rather_than_satisfied(self):
        from neoswga.core.panel_refinement import objective_for_optimizer
        from neoswga.core.pool_objective import PoolConstraints

        optimizer = _two_target_optimizer(geometry="directional")
        objective = objective_for_optimizer(
            optimizer, PoolConstraints(min_per_target_coverage=0.5, coverage_metric="raw")
        )

        assert "per-target coverage below minimum" in objective.violations(A_PRIMERS)

    def test_the_reporting_pass_does_not_substitute_a_symmetric_number(self):
        from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus
        from neoswga.core.target_coverage import attach_per_target

        optimizer = _two_target_optimizer(geometry="directional")
        result = OptimizationResult(
            primers=[*A_PRIMERS, B_PRIMER],
            score=1.0,
            status=OptimizationStatus.SUCCESS,
            metrics=optimizer.compute_metrics([*A_PRIMERS, B_PRIMER]),
            iterations=1,
            optimizer_name="probe",
        )
        attached = attach_per_target(
            result,
            optimizer.cache,
            optimizer.fg_prefixes,
            optimizer.fg_seq_lengths,
            extension=REACH,
            circular=False,
            geometry="directional",
        )

        assert attached.metrics.per_target_coverage == {}

    def test_the_optimize_path_passes_the_reason_through(self):
        import inspect

        from neoswga.core import unified_optimizer

        source = inspect.getsource(unified_optimizer.run_optimization)

        assert "DIRECTIONAL_UNAVAILABLE" in source

    def test_the_schema_says_so(self):
        import json

        schema = json.loads(
            (ROOT / "neoswga" / "core" / "schema" / "params.schema.json").read_text()
        )
        description = schema["properties"]["min_per_target_coverage"]["description"]

        assert "directional" in description
        assert "REFUSED" in description


class TestTheRecordedMetricIsTheOneThatApplied:
    """`constraints_from_parameter` never set `coverage_metric`, so a recorded
    request whose run selected on RAW coverage named `effective` anyway -- and
    `constraints` is hashed, so the CLI overlay made the hash claim to cover the
    constraints faithfully while one field was a placeholder."""

    def _params(self, **extra):
        return dict(
            schema_version=2,
            fg_prefixes=["/tmp/t"],
            fg_seq_lengths=[100_000],
            fg_genomes=["/tmp/t.fna"],
            polymerase="phi29",
            reaction_temp=30.0,
            num_primers=6,
            target_set_size=6,
            min_per_target_coverage=0.5,
            **extra,
        )

    def test_a_stated_metric_reaches_the_recorded_constraints(self):
        from neoswga.core.design_request import resolve_design_request

        request = resolve_design_request(self._params(coverage_metric="raw"))

        assert request.constraints.coverage_metric == "raw"

    def test_a_stated_metric_moves_the_request_hash(self):
        from neoswga.core.design_request import resolve_design_request

        raw = resolve_design_request(self._params(coverage_metric="raw"))
        effective = resolve_design_request(self._params(coverage_metric="effective"))

        assert raw.request_hash != effective.request_hash

    def test_an_unstated_metric_is_recorded_as_resolved_later(self):
        """Not silently named: the record says the choice was not made here."""
        from neoswga.core.design_request import resolve_design_request

        request = resolve_design_request(self._params())

        assert "resolved at selection" in request.default_sources["coverage_metric"]

    def test_a_stated_metric_is_recorded_as_coming_from_the_request(self):
        from neoswga.core.design_request import resolve_design_request

        request = resolve_design_request(self._params(coverage_metric="raw"))

        assert request.default_sources["coverage_metric"] == "request"

    def test_no_limit_names_no_metric_at_all(self):
        """With no limit there is no constraints object, so nothing can be
        misdescribed."""
        from neoswga.core.design_request import resolve_design_request

        params = self._params()
        del params["min_per_target_coverage"]

        assert resolve_design_request(params).constraints is None


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
        from neoswga.core.target_coverage import require_per_target

        with pytest.raises(ModelEvaluationError):
            require_per_target({}, None, ["A", "B"], 0.5)

    def test_the_empty_dict_refusal_says_that_is_what_happened(self):
        """The empty-dict guard, ISOLATED. Removing only its `raise` leaves the
        `worst is None` guard to catch the same input, so `pytest.raises` alone
        cannot tell the two apart -- the review's mutation M2 showed that, and
        corrected an overstated row in the first report. Asserting the distinct
        message is what separates them: with the first raise gone this reads
        "coverage is unavailable for 2 of 2 targets" instead.
        """
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.target_coverage import require_per_target

        with pytest.raises(ModelEvaluationError) as excinfo:
            require_per_target({}, None, ["A", "B"], 0.5)

        assert "no target's coverage could be measured" in str(excinfo.value)

    def test_an_unmeasured_target_with_a_requested_floor_is_refused(self):
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.target_coverage import require_per_target

        with pytest.raises(ModelEvaluationError):
            require_per_target({"A": 0.9}, None, ["A", "B"], 0.5)

    def test_the_refusal_names_the_unmeasured_target(self):
        from neoswga.core.exceptions import ModelEvaluationError
        from neoswga.core.target_coverage import require_per_target

        with pytest.raises(ModelEvaluationError) as excinfo:
            require_per_target({"A": 0.9}, None, ["A", "B"], 0.5)

        assert "B" in str(excinfo.value)

    def test_no_floor_requested_refuses_nothing(self):
        """What keeps this step from moving a delivered panel: with the key
        unset, an absent per-target measurement is not an error."""
        from neoswga.core.target_coverage import require_per_target

        assert require_per_target({}, None, ["A", "B"], None) is None
        assert require_per_target({}, None, ["A", "B"], 0.0) is None

    def test_a_measured_worst_passes(self):
        from neoswga.core.target_coverage import require_per_target

        assert require_per_target({"A": 0.9, "B": 0.1}, 0.1, ["A", "B"], 0.5) is None

    def test_the_guard_is_called_from_the_optimize_path(self):
        """A guard nothing calls is the same as no guard."""
        import ast
        import inspect

        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        called = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "require_per_target" in called

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
