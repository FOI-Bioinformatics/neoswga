"""swga 2.0's second-largest fitted weight, computed here for the first time.

Phase 7.2 of `docs/reports/PLAN_2026-10-28_pool_design_trustworthiness.md`;
measured in `docs/validation/2026-10-29-off-gap-gini.md`.

`off_gap_gini` is a Gini index of the distances between background binding
sites, computed per strand, averaged across strands and then across genomes.
`docs/validation/pool_selection_audit_2026-09-18.md` records it as carrying
the second-largest weight in the only set-level model fitted against measured
sequencing breadth, and `docs/validation/tool_comparison.md` recorded three
times that nothing here computed it.

It is REPORTED and never scored. The sample behind the weight is 46 sets with
one response variable, read from a table nobody here could open, and COATswga
computes no gap statistic at all on the stated ground that a per-primer Gini
cannot speak for a set. `test_it_does_not_enter_normalized_score` pins that
position, because the cheapest way to get this wrong is to add a term.

**The zero is the unfavourable end of this scale**, which makes "unknown is not
zero" sharper than usual here. A positive fitted weight means an UNEVEN
background is better, since clustered host sites leave most of the host out of
a polymerase's reach. So a Gini of 0.0 reads as "perfectly even host, bad
news", and a strand nobody could measure must not land there. A strand with
fewer than two sites has no distance at all: it is left out of the average, not
averaged in as zero.
"""

import pytest

from neoswga.core.background_evenness import GINI_KEY, gini, off_gap_gini, per_strand_gap_gini


class TestTheGiniItself:
    def test_no_distances_is_not_a_measurement(self):
        assert gini([]) is None

    def test_equal_distances_are_perfectly_even(self):
        assert gini([100.0, 100.0, 100.0]) == 0.0

    def test_one_distance_is_even_and_measured(self):
        """A single distance is a measurement, unlike none at all."""
        assert gini([500.0]) == 0.0

    def test_clustered_distances_score_higher_than_dispersed(self):
        clustered = gini([1.0, 1.0, 1.0, 997.0])
        dispersed = gini([250.0, 250.0, 250.0, 250.0])

        assert clustered > dispersed == 0.0

    def test_all_zero_distances_is_zero_not_a_division_error(self):
        assert gini([0.0, 0.0]) == 0.0


class TestThePerStrandAverage:
    def test_both_strands_are_averaged(self):
        value = per_strand_gap_gini(
            {"forward": [0, 1, 2, 1000], "reverse": [0, 500, 1000]},
        )

        expected = (gini([1.0, 1.0, 998.0]) + 0.0) / 2
        assert value == pytest.approx(expected)

    def test_a_strand_with_one_site_is_left_out_not_averaged_as_zero(self):
        """Zero is the end of this scale that reads as bad news, so averaging
        an unmeasurable strand in at 0.0 would halve a real measurement."""
        measurable = {"forward": [0, 1, 2, 1000]}
        with_a_lone_reverse_site = {"forward": [0, 1, 2, 1000], "reverse": [7]}

        assert per_strand_gap_gini(with_a_lone_reverse_site) == pytest.approx(
            per_strand_gap_gini(measurable)
        )

    def test_neither_strand_measurable_is_absent(self):
        assert per_strand_gap_gini({"forward": [7], "reverse": []}) is None

    def test_an_empty_mapping_is_absent(self):
        assert per_strand_gap_gini({}) is None

    def test_duplicate_positions_do_not_invent_a_zero_distance(self):
        """Two primers binding the same site is one site, not a gap of zero."""
        assert per_strand_gap_gini({"forward": [10, 10, 10]}) is None


class TestTheBackgroundAverage:
    def test_a_design_with_no_host_has_no_host_evenness(self):
        assert off_gap_gini({"target": {GINI_KEY: 0.4}}, []) is None

    def test_an_unmeasured_host_is_absent_rather_than_zero(self):
        assert off_gap_gini({"host": {GINI_KEY: None}}, ["host"]) is None

    def test_a_host_missing_from_the_stats_is_absent(self):
        assert off_gap_gini({}, ["host"]) is None

    def test_two_hosts_are_averaged(self):
        value = off_gap_gini({"h1": {GINI_KEY: 0.2}, "h2": {GINI_KEY: 0.6}}, ["h1", "h2"])

        assert value == pytest.approx(0.4)

    def test_only_the_background_prefixes_count(self):
        """The target's own figure is swga 2.0's `on_gap_gini`, a different
        term with a different sign, and must not leak into this one."""
        stats = {"target": {GINI_KEY: 1.0}, "host": {GINI_KEY: 0.3}}

        assert off_gap_gini(stats, ["host"]) == pytest.approx(0.3)

    def test_one_measured_host_among_two_is_not_diluted_by_the_other(self):
        stats = {"h1": {GINI_KEY: 0.5}, "h2": {GINI_KEY: None}}

        assert off_gap_gini(stats, ["h1", "h2"]) == pytest.approx(0.5)


class TestItReachesTheArtifacts:
    def test_the_metrics_carry_it(self):
        import dataclasses

        from neoswga.core.base_optimizer import PrimerSetMetrics

        metrics = dataclasses.replace(PrimerSetMetrics.empty(), off_gap_gini=0.37)

        assert metrics.to_dict()["off_gap_gini"] == 0.37

    def test_an_unmeasured_panel_serialises_it_as_absent(self):
        from neoswga.core.base_optimizer import PrimerSetMetrics

        assert PrimerSetMetrics.empty().to_dict()["off_gap_gini"] is None

    def test_it_is_reported_as_a_panel_criterion_with_no_reference(self):
        from neoswga.core.panel_regime import _spacing_criteria

        class _M:
            max_gap = 4000.0
            mean_gap = 500.0
            gap_gini = 0.4
            bg_coverage = 0.01
            strand_coverage_ratio = 0.9
            strand_alternation_score = 0.6
            off_gap_gini = 0.37

        named = {c.name: c for c in _spacing_criteria(_M(), 3000, 1_270_000)}

        assert named["host_gap_evenness"].value == pytest.approx(0.37)
        assert "No reference" in named["host_gap_evenness"].note

    def test_an_unmeasured_host_is_reported_as_absent_not_as_even(self):
        from neoswga.core.panel_regime import _spacing_criteria

        class _M:
            max_gap = 4000.0
            mean_gap = 500.0
            gap_gini = 0.4
            bg_coverage = None
            strand_coverage_ratio = None
            strand_alternation_score = None
            off_gap_gini = None

        named = {c.name: c for c in _spacing_criteria(_M(), 3000, 1_270_000)}

        assert named["host_gap_evenness"].value is None

    def test_the_strand_section_renders_the_per_genome_figure(self):
        from neoswga.core.report.strand_section import render_strand_balance

        class _U:
            from_optimizer = True
            strand_alternation_score = 0.6
            strand_coverage_ratio = 0.9
            strand_stats = {"host": {GINI_KEY: 0.37}}
            strand_roles = {"host": "host"}

        html = render_strand_balance(_U(), 3000)

        assert "Gap evenness per strand" in html
        assert "0.37" in html


class TestItIsNotScored:
    def test_it_does_not_enter_normalized_score(self):
        """The sample does not support a weight. See the module docstring of
        `core/background_evenness.py` and Known Issue 18."""
        import dataclasses

        from neoswga.core.base_optimizer import PrimerSetMetrics

        def _scored(value):
            return dataclasses.replace(
                PrimerSetMetrics.empty(),
                fg_coverage=0.5,
                selectivity_ratio=10.0,
                dimer_risk_score=0.1,
                gap_gini=0.4,
                off_gap_gini=value,
            ).normalized_score()

        assert _scored(0.0) == _scored(0.9) == _scored(None)

    def test_the_search_objective_does_not_carry_it(self):
        """`PoolMetrics` has five fields on purpose; a sixth is an
        `AttributeError` rather than a default. See
        `docs/validation/2026-10-29-stage2-background-weight-does-not-run.md`."""
        from neoswga.core.pool_metrics import PoolMetrics

        assert "off_gap_gini" not in getattr(PoolMetrics, "__annotations__", {})


class TestTheSourceComputesIt:
    """Against a real position index, which is where the positions come from."""

    @staticmethod
    def _cache(tmp_path, sites_by_primer, k=10):
        h5py = pytest.importorskip("h5py")
        import numpy as np

        from neoswga.core.position_cache import PositionCache

        tmp_path.mkdir(parents=True, exist_ok=True)
        prefix = str(tmp_path / "g")
        primers = list(sites_by_primer)
        with h5py.File(f"{prefix}_{k}mer_positions.h5", "w") as db:
            for primer, positions in sites_by_primer.items():
                db.create_dataset(primer, data=np.array(positions, dtype=np.int64))
            db.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))
        return PositionCache([prefix], primers), prefix

    def test_the_key_is_returned(self, tmp_path):
        cache, prefix = self._cache(tmp_path, {"ACGTACGTAC": [100, 2000, 40_000]})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], 50_000)

        assert GINI_KEY in stats
        assert stats[GINI_KEY] is not None

    def test_no_sites_at_all_is_absent(self, tmp_path):
        cache, prefix = self._cache(tmp_path, {"ACGTACGTAC": []})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], 50_000)

        assert stats[GINI_KEY] is None

    def test_a_clustered_genome_scores_higher_than_a_dispersed_one(self, tmp_path):
        even, _ = self._cache(tmp_path / "even", {"ACGTACGTAC": [0, 10_000, 20_000, 30_000]})
        clumped, _ = self._cache(tmp_path / "clumped", {"ACGTACGTAC": [0, 10, 20, 30_000]})
        even_prefix = str(tmp_path / "even" / "g")
        clumped_prefix = str(tmp_path / "clumped" / "g")

        a = even.compute_strand_alternation_stats(even_prefix, ["ACGTACGTAC"], 50_000)
        b = clumped.compute_strand_alternation_stats(clumped_prefix, ["ACGTACGTAC"], 50_000)

        assert b[GINI_KEY] > a[GINI_KEY]


class TestItReachesTheMetricsFromTheSearch:
    """The PATH, not the two ends.

    `off_gap_gini` existed on `PrimerSetMetrics`, was serialised, and was
    reported, with nothing asserting that `_compute_metrics` actually fills it
    from the background's strand stats. Replacing the call with a literal
    `None` passed every other test in this file, which is the Known Issue 8
    shape exactly, so this drives the real aggregation.
    """

    @staticmethod
    def _optimizer(tmp_path):
        h5py = pytest.importorskip("h5py")
        import numpy as np

        from neoswga.core.base_optimizer import BaseOptimizer, OptimizationResult
        from neoswga.core.position_cache import PositionCache

        primer = "ACGTACGTAC"
        rc = "GTACGTACGT"
        # The host's sites are clustered on one strand and spread on the
        # other, so its per-strand figure is neither 0 nor the target's.
        planted = {
            "target": {primer: [100, 9_000, 18_000], rc: [4_000, 8_000, 12_000]},
            "host": {primer: [10, 20, 30, 40_000], rc: [0, 15_000, 30_000]},
        }
        prefixes = {}
        for name, entries in planted.items():
            prefix = str(tmp_path / name)
            with h5py.File(f"{prefix}_10mer_positions.h5", "w") as db:
                for key, positions in entries.items():
                    db.create_dataset(key, data=np.array(positions, dtype=np.int64))
                db.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))
            prefixes[name] = prefix

        class _Minimal(BaseOptimizer):
            def optimize(self, candidates, target_size=10, **kwargs):
                return OptimizationResult.failure("minimal", "not used")

            def _score_set(self, primers):
                return 0.0

            @property
            def name(self) -> str:
                return "minimal"

        cache = PositionCache([prefixes["target"], prefixes["host"]], [primer])
        optimizer = _Minimal(
            cache,
            [prefixes["target"]],
            [50_000],
            bg_prefixes=[prefixes["host"]],
            bg_seq_lengths=[60_000],
        )
        return optimizer, [primer], prefixes

    def test_the_scalar_is_filled_from_the_background_stats(self, tmp_path):
        optimizer, primers, prefixes = self._optimizer(tmp_path)

        metrics = optimizer.compute_metrics(primers)

        expected = off_gap_gini(metrics.strand_stats, [prefixes["host"]])
        assert expected is not None
        assert metrics.off_gap_gini == pytest.approx(expected)

    def test_it_is_the_hosts_figure_and_not_the_targets(self, tmp_path):
        """swga 2.0's `on_gap_gini` is a different term with a different sign."""
        optimizer, primers, prefixes = self._optimizer(tmp_path)

        metrics = optimizer.compute_metrics(primers)
        target = metrics.strand_stats[prefixes["target"]][GINI_KEY]

        assert metrics.off_gap_gini is not None
        assert target is not None
        assert metrics.off_gap_gini != pytest.approx(target)

    def test_a_design_with_no_background_reports_no_host_evenness(self, tmp_path):
        optimizer, primers, _prefixes = self._optimizer(tmp_path)
        optimizer.bg_prefixes = []
        optimizer.bg_seq_lengths = []

        assert optimizer.compute_metrics(primers).off_gap_gini is None
