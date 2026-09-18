"""A delivered panel should say which criterion limits it, and which have none.

Two wet-lab benchmarks in `docs/validation/published_primer_sets.md` disagree
about which spacing property predicts success, and the conclusion drawn there is
that what predicts is whichever property is currently LIMITING. Nothing reported
which one that was.

This is the diagnostic form of that conclusion, and the tests below exist mostly
to stop it becoming something else. Two failure modes are specifically guarded:

- **Manufacturing a reference.** `max_gap` has no non-fitted line to compare
  against: measured on the 18 published sets with wet-lab outcomes, every one
  exceeds twice the calibrated reach, winners included
  (`docs/validation/getting_ahead_on_spacing_2026-09-18.md`). A criterion with
  no reference must be REPORTED and must never be ranked, or the diagnostic
  would name the worst hole as limiting on every real design and say nothing.
- **Compressing the criteria into one number.** That is what the same document
  rules out. Every criterion stays separately addressable.
"""

import pytest

from neoswga.core.panel_regime import PanelRegime, assess_panel, format_regime


class _Metrics:
    """The fields `assess_panel` reads, and nothing else.

    A stub rather than a real `PrimerSetMetrics` so a test states exactly which
    quantities the diagnostic depends on. A field this stub lacks is a field the
    production code must not read.
    """

    def __init__(
        self,
        *,
        effective_fg_coverage=0.85,
        fg_coverage=0.90,
        selectivity_ratio=8.0,
        selectivity_density=45.0,
        total_bg_sites=120,
        max_gap=31_000.0,
        mean_gap=4_800.0,
        gap_gini=0.48,
        bg_coverage=0.02,
        strand_coverage_ratio=0.9,
        strand_alternation_score=0.7,
    ):
        self.effective_fg_coverage = effective_fg_coverage
        self.fg_coverage = fg_coverage
        self.selectivity_ratio = selectivity_ratio
        self.selectivity_density = selectivity_density
        self.total_bg_sites = total_bg_sites
        self.max_gap = max_gap
        self.mean_gap = mean_gap
        self.gap_gini = gap_gini
        self.bg_coverage = bg_coverage
        self.strand_coverage_ratio = strand_coverage_ratio
        self.strand_alternation_score = strand_alternation_score


def _assess(metrics=None, **overrides):
    """One assessment with the defaults a plain `optimize` run would supply."""
    kwargs = dict(
        coverage_target=0.80,
        min_fg_bg_ratio=5.0,
        requested_size=12,
        delivered_size=12,
        coverage_reach=3_000,
        genome_length=3_200_000,
    )
    kwargs.update(overrides)
    return assess_panel(metrics or _Metrics(), **kwargs)


def _named(regime, name):
    for criterion in regime.criteria:
        if criterion.name == name:
            return criterion
    raise AssertionError(f"no criterion named {name!r} in {[c.name for c in regime.criteria]}")


class TestWhatCarriesAReference:
    """Only a line the user or physics drew counts as a reference."""

    def test_a_panel_clearing_every_reference_fails_nothing(self):
        regime = _assess()

        assert regime.failing == ()

    def test_the_tightest_criterion_is_named_even_when_all_of_them_pass(self):
        """`limiting` answers "what is closest to binding", not "what broke".

        A design that clears everything still has one criterion nearest its
        limit, and that is the one worth knowing about before changing anything.
        """
        regime = _assess()

        assert regime.limiting is not None
        assert _named(regime, regime.limiting).slack is not None

    def test_a_panel_under_the_coverage_target_reports_coverage_failing(self):
        regime = _assess(_Metrics(effective_fg_coverage=0.40))

        assert "coverage" in regime.failing
        assert regime.limiting == "coverage"

    def test_a_short_panel_reports_its_size_failing(self):
        regime = _assess(delivered_size=7)

        assert "panel_size" in regime.failing

    def test_slack_is_relative_so_two_criteria_are_comparable(self):
        """Missing a coverage target by half must score like missing a ratio
        floor by half.

        Without this the criterion measured on the larger scale always looks
        worse, which is the defect `PoolObjective.shortfall` was given relative
        terms to avoid. Here it would make `limiting` a function of units.
        """
        half_coverage = _assess(_Metrics(effective_fg_coverage=0.40))
        half_ratio = _assess(_Metrics(selectivity_ratio=2.5))

        assert _named(half_coverage, "coverage").slack == pytest.approx(
            _named(half_ratio, "fg_bg_ratio").slack
        )

    def test_an_unset_optional_limit_creates_no_criterion(self):
        """A limit nobody set is not a criterion this panel passed.

        Reporting it as satisfied would read as specificity where there is only
        an absent constraint, which is the shape of Known Issues 5, 6 and 13.
        """
        regime = _assess()

        assert "selectivity_density" not in [c.name for c in regime.criteria]

    def test_a_configured_density_floor_becomes_a_criterion(self):
        regime = _assess(min_selectivity_density=60.0)

        criterion = _named(regime, "selectivity_density")
        assert criterion.reference == 60.0
        assert criterion.slack is not None
        assert "selectivity_density" in regime.failing


class TestWhatCarriesNoReference:
    """The spacing properties, which is the point of the exercise."""

    @pytest.mark.parametrize("name", ["worst_hole", "mean_gap", "evenness", "host_coverage"])
    def test_it_is_reported_without_a_reference(self, name):
        regime = _assess()

        criterion = _named(regime, name)
        assert criterion.reference is None
        assert criterion.slack is None
        assert name in regime.unreferenced

    def test_limiting_never_names_an_unreferenced_criterion(self):
        """The guard that keeps this a diagnostic rather than a slogan.

        Every one of the 18 published sets, winners included, has a worst hole
        beyond twice the reach. Ranking an unreferenced criterion would name
        `worst_hole` as limiting on essentially every design.
        """
        regime = _assess(_Metrics(max_gap=2_000_000.0, mean_gap=900_000.0))

        assert regime.limiting not in regime.unreferenced

    def test_the_worst_hole_is_reported_in_multiples_of_the_reach(self):
        """A bare base count cannot be read without knowing the chemistry."""
        regime = _assess(_Metrics(max_gap=30_000.0), coverage_reach=3_000)

        assert "10" in _named(regime, "worst_hole").note

    def test_host_coverage_is_reported_at_all(self):
        """`bg_coverage` had no reader anywhere; see Known Issue 18.

        It is the only computed quantity that sees background site POSITION, so
        losing it again would take the diagnostic's one positional term with it.
        """
        regime = _assess(_Metrics(bg_coverage=0.037))

        assert _named(regime, "host_coverage").value == pytest.approx(0.037)


class TestCoverageProvenance:
    """Which of the two coverage figures was judged, said rather than inferred."""

    def test_the_occupancy_weighted_figure_is_preferred(self):
        regime = _assess(_Metrics(effective_fg_coverage=0.55, fg_coverage=0.95))

        assert _named(regime, "coverage").value == pytest.approx(0.55)

    def test_it_falls_back_to_raw_coverage_and_says_so(self):
        """`None` means the evaluator could not compute it, not zero."""
        regime = _assess(_Metrics(effective_fg_coverage=None, fg_coverage=0.95))

        criterion = _named(regime, "coverage")
        assert criterion.value == pytest.approx(0.95)
        assert "raw" in criterion.note


class TestTheReport:
    def test_every_criterion_reaches_the_formatted_output(self):
        """No criterion may be computed and then left out of what is printed."""
        regime = _assess(min_selectivity_density=60.0, max_background_sites=50)
        lines = "\n".join(format_regime(regime))

        for criterion in regime.criteria:
            assert criterion.name in lines

    def test_the_output_names_the_limiting_criterion(self):
        regime = _assess(_Metrics(effective_fg_coverage=0.40))
        lines = "\n".join(format_regime(regime))

        assert "coverage" in lines
        assert regime.limiting in lines

    def test_the_output_says_which_criteria_nobody_constrained(self):
        regime = _assess()
        lines = "\n".join(format_regime(regime)).lower()

        assert "no reference" in lines

    def test_a_regime_with_nothing_to_report_formats_to_nothing(self):
        empty = PanelRegime(criteria=(), limiting=None, failing=(), unreferenced=())

        assert format_regime(empty) == []


class TestItReachesTheRun:
    """The diagnostic has to be wired, or it is the defect it was written after.

    Known Issue 16 is a parameter that exists, is documented and is passed by
    nothing. These assert the path rather than the pieces.
    """

    def test_a_failed_assessment_logs_nothing_and_does_not_raise(self):
        """`assess_from_parameter` returns None when it cannot resolve a
        reference, and a report of nothing must stay silent rather than print a
        limit nothing established."""
        from neoswga.core.panel_regime import log_regime

        log_regime(None)

    def test_an_unknown_application_yields_no_assessment(self):
        from neoswga.core.panel_regime import assess_from_parameter

        class _Params:
            polymerase = "phi29"
            coverage_reach = None
            fg_seq_lengths = [3_200_000]
            num_primers = 12

        assert (
            assess_from_parameter(_Metrics(), _Params(), delivered_size=12, application="nonsense")
            is None
        )

    def test_the_summary_payload_round_trips_as_json(self):
        """It is written into `step4_improved_df_summary.json`, so every value
        has to be JSON-serialisable without a custom encoder."""
        import json

        from neoswga.core.panel_regime import criteria_for_summary

        payload = criteria_for_summary(_assess(min_selectivity_density=60.0))

        assert json.loads(json.dumps(payload))["limiting"] == payload["limiting"]

    def test_the_optimize_path_passes_a_regime_to_the_summary_writer(self):
        """Asserted on the source of the call, because the two ends existing is
        what Known Issue 16 is about."""
        import ast
        import inspect

        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        passes = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and getattr(node.func, "id", None) == "save_results"
            and any(kw.arg == "regime" for kw in node.keywords)
        ]
        assert passes, "save_results is called without a regime, so the summary has none"

    def test_the_summary_writer_accepts_and_records_it(self, tmp_path):
        import inspect

        from neoswga.core.step4_output import save_results

        assert "regime" in inspect.signature(save_results).parameters


class TestTheWordingIsTrueInBothDirections:
    """`limiting` is the minimum slack, which means two different things.

    On a panel that fails something it is the criterion missed by the widest
    relative margin. On a panel that clears everything it is the one with the
    least headroom. Calling both "closest to binding" is wrong in the first
    case, and that is the case a user sees when something needs fixing.
    """

    def test_a_failing_panel_says_the_criterion_is_furthest_from_its_reference(self):
        regime = _assess(_Metrics(effective_fg_coverage=0.10))
        lines = "\n".join(format_regime(regime)).lower()

        assert "furthest from" in lines
        assert "closest to binding" not in lines

    def test_a_passing_panel_says_the_criterion_is_closest_to_binding(self):
        regime = _assess()
        lines = "\n".join(format_regime(regime)).lower()

        assert regime.failing == ()
        assert "closest to binding" in lines
        assert "furthest from" not in lines


class TestAnUncomputedZeroIsNotAMeasurement:
    """Three quantities default to 0.0 when they could not be computed.

    `PrimerSetMetrics.empty` sets `bg_coverage`, `strand_coverage_ratio` and
    `strand_alternation_score` to 0.0, and `_compute_metrics` leaves the two
    strand figures at 0.0 when the cache cannot supply them. A zero there is
    indistinguishable from a measured zero, which is the shape of Known Issues
    5, 6 and 13. The diagnostic must not launder it as a measurement.
    """

    def test_a_host_coverage_zero_carries_the_caveat(self):
        """`bg_coverage` is still 0.0 both when measured zero and when nothing
        measured it, so its zero still needs the caveat."""
        regime = _assess(_Metrics(bg_coverage=0.0))

        assert "may mean" in _named(regime, "host_coverage").note

    def test_a_non_zero_host_coverage_does_not(self):
        regime = _assess(_Metrics(bg_coverage=0.25))

        assert "may mean" not in _named(regime, "host_coverage").note

    @pytest.mark.parametrize(
        "name,field",
        [
            ("strand_balance", "strand_coverage_ratio"),
            ("strand_alternation", "strand_alternation_score"),
        ],
    )
    def test_the_strand_zeros_no_longer_need_a_caveat(self, name, field):
        """`core/strand_metrics.py` made these `None` when uncomputed, so a 0.0
        is now a measurement and captioning it would be wrong."""
        regime = _assess(_Metrics(**{field: 0.0}))

        assert "may mean" not in _named(regime, name).note

    @pytest.mark.parametrize(
        "name,field",
        [
            ("strand_balance", "strand_coverage_ratio"),
            ("strand_alternation", "strand_alternation_score"),
        ],
    )
    def test_an_uncomputed_strand_value_reads_as_not_measured(self, name, field):
        regime = _assess(_Metrics(**{field: None}))

        assert _named(regime, name).value is None


class TestThePrintedTableDoesNotLaunderAnAmbiguousZero:
    """The note reaches the summary JSON. The table shows only the value.

    A reader looking at `host_coverage 0` in the terminal would read a
    measurement, which is the thing the note exists to prevent. The caveat has
    to survive into the printed output too.
    """

    def test_the_printed_output_names_the_ambiguous_zero(self):
        regime = _assess(_Metrics(bg_coverage=0.0))
        lines = "\n".join(format_regime(regime))

        assert "host_coverage" in lines
        assert "may not have been computed" in lines

    def test_it_says_nothing_when_no_value_is_ambiguous(self):
        regime = _assess(
            _Metrics(bg_coverage=0.02, strand_coverage_ratio=0.0, strand_alternation_score=0.0)
        )
        lines = "\n".join(format_regime(regime))

        assert "may not have been computed" not in lines


class TestTheHostConvergentGapReachesTheReport:
    """`strand_alternation_gap_max` on the HOST is the closest quantity here to
    the off-target amplification mechanism, and swga 2.0 fits its proxy against
    measured sequencing breadth. Item 3 started computing it; it has to arrive
    somewhere a user sees."""

    def test_it_is_reported_when_the_host_was_measured(self):
        metrics = _Metrics()
        metrics.strand_stats = {
            "target": {"strand_alternation_gap_max": 41_000.0},
            "host": {"strand_alternation_gap_max": 900.0},
        }

        regime = _assess(metrics, bg_prefixes=["host"])

        criterion = _named(regime, "host_convergent_gap")
        assert criterion.value == pytest.approx(900.0)
        assert criterion.reference is None

    def test_it_is_absent_when_no_host_was_measured(self):
        """Absent rather than zero: an unmeasured host must not read as one
        whose sites are all adjacent."""
        metrics = _Metrics()
        metrics.strand_stats = {"target": {"strand_alternation_gap_max": 41_000.0}}

        regime = _assess(metrics, bg_prefixes=["host"])

        assert "host_convergent_gap" not in [c.name for c in regime.criteria]

    def test_the_target_gap_is_reported_too(self):
        metrics = _Metrics()
        metrics.strand_stats = {"target": {"strand_alternation_gap_max": 41_000.0}}

        regime = _assess(metrics, fg_prefixes=["target"])

        assert _named(regime, "convergent_gap").value == pytest.approx(41_000.0)
