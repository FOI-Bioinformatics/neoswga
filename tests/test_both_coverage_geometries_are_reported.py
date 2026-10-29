"""The headline coverage figure says which window produced it.

Phase 5.2 of `docs/reports/PLAN_2026-10-28_pool_design_trustworthiness.md`;
decision and figures in
`docs/validation/2026-10-29-the-coverage-geometry-decision.md`.

`fg_coverage` is computed under one of two geometries and reported as a bare
number. Measured on 2026-10-29, the shipped `symmetric` window overstates the
directional figure by 1.27x to 1.85x on three real panels at the shipped 3 kb
reach (`docs/validation/2026-10-29-directional-coverage-measured.md`). A reader
of a summary could not tell which convention produced the number, nor by how
much the other would differ.

Two things follow, and neither changes the default:

**Both figures are reported for the delivered panel, labelled.** Once, after
the search: `_compute_metrics` runs per candidate panel and fetches the second
orientation only when the configured geometry asks for it, because
`get_positions` is the hottest method in the package.

**A directional run must name its own reach.** The polymerase defaults were
fitted under symmetric windows -- 2,866-4,551 bp enters and leaves the
published breadth band, against 4,354-6,710 bp directional
(`docs/validation/2026-09-23-reach-refit-directional.md`) -- so a directional
run at the shipped 3,000 reads 0.48 on the set whose wet-lab outcome recorded
60-75%. That is not a more correct model, it is a model mis-calibrated by
construction, and it was reachable from params.json with one key. It is refused
rather than given a substitute: the directional band was fitted for phi29
against one outcome and none exists for `bst`, `equiphi29` or `klenow`, so a
default would invent three numbers from one measurement.
"""

import pytest

from neoswga.core.coverage_geometries import compare_geometries
from neoswga.core.design_request import resolve_design_request
from neoswga.core.exceptions import InvalidDesignRequest

GENOME = 60_000
PRIMER = "ACGTACGTAC"


class _Cache:
    """Answers with planted positions, per prefix and strand."""

    def __init__(self, planted):
        self._planted = planted

    def get_positions(self, prefix, primer, strand="both"):
        return list(self._planted.get((prefix, strand), ()))


def _cache(forward=(1_000, 20_000, 40_000), reverse=(10_000, 30_000)):
    return _Cache({("t", "forward"): forward, ("t", "reverse"): reverse})


class TestTheComparison:
    def test_the_wider_window_covers_more(self):
        out = compare_geometries(_cache(), ["t"], [GENOME], [PRIMER], reach=2_000)

        assert out["symmetric"] > out["directional"]
        assert out["inflation"] == pytest.approx(out["symmetric"] / out["directional"])

    def test_an_isolated_site_is_exactly_twice_as_wide(self):
        """The analytic case: one site marks 2r symmetrically and r
        directionally, so with no overlap the ratio is exactly 2."""
        cache = _Cache({("t", "forward"): (30_000,), ("t", "reverse"): ()})

        out = compare_geometries(cache, ["t"], [GENOME], [PRIMER], reach=1_000)

        assert out["symmetric"] == pytest.approx(2_000 / GENOME)
        assert out["directional"] == pytest.approx(1_000 / GENOME)
        assert out["inflation"] == pytest.approx(2.0)

    def test_the_label_names_the_geometry_the_run_scored_on(self):
        out = compare_geometries(
            _cache(), ["t"], [GENOME], [PRIMER], reach=2_000, reported="directional"
        )

        assert out["reported_geometry"] == "directional"
        assert "selected and scored under 'directional'" in out["note"]

    def test_the_reach_is_recorded_beside_the_figures(self):
        """Coverage is meaningless without the reach it was computed at."""
        out = compare_geometries(_cache(), ["t"], [GENOME], [PRIMER], reach=2_345)

        assert out["reach"] == 2_345

    def test_an_empty_panel_has_no_comparison(self):
        assert compare_geometries(_cache(), ["t"], [GENOME], [], reach=2_000) is None

    def test_no_cache_has_no_comparison(self):
        assert compare_geometries(None, ["t"], [GENOME], [PRIMER], reach=2_000) is None

    def test_no_target_length_is_absent_rather_than_zero_coverage(self):
        out = compare_geometries(_cache(), ["t"], [0], [PRIMER], reach=2_000)

        assert out is None

    def test_a_panel_that_binds_nowhere_covers_nothing_and_has_no_ratio(self):
        """A ratio over zero is not a large number, it is no number."""
        cache = _Cache({})

        out = compare_geometries(cache, ["t"], [GENOME], [PRIMER], reach=2_000)

        assert out["symmetric"] == 0.0
        assert out["directional"] == 0.0
        assert out["inflation"] is None

    def test_two_targets_are_weighted_by_length(self):
        cache = _Cache(
            {
                ("a", "forward"): (5_000,),
                ("b", "forward"): (),
            }
        )

        out = compare_geometries(cache, ["a", "b"], [10_000, 30_000], [PRIMER], reach=1_000)

        # 1,000 bp covered in `a` directionally, nothing in `b`, over 40,000.
        assert out["directional"] == pytest.approx(1_000 / 40_000)


class TestItAgreesWithTheFigureTheRunReports:
    """The PATH: the symmetric column must BE `fg_coverage`, not merely
    resemble it. Three implementations of one quantity is how this codebase has
    produced disagreeing coverage numbers before."""

    def test_the_symmetric_column_equals_the_optimizers_own_coverage(self, tmp_path):
        h5py = pytest.importorskip("h5py")
        import numpy as np

        from neoswga.core.base_optimizer import (
            BaseOptimizer,
            OptimizationResult,
            OptimizerConfig,
        )
        from neoswga.core.position_cache import PositionCache
        from neoswga.core.thermodynamics import reverse_complement

        prefix = str(tmp_path / "t")
        rc = reverse_complement(PRIMER)
        with h5py.File(f"{prefix}_10mer_positions.h5", "w") as db:
            db.create_dataset(PRIMER, data=np.array([1_000, 20_000, 40_000], dtype=np.int64))
            db.create_dataset(rc, data=np.array([10_000, 30_000], dtype=np.int64))
            db.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))

        class _Minimal(BaseOptimizer):
            def optimize(self, candidates, target_size=10, **kwargs):
                return OptimizationResult.failure("minimal", "not used")

            def _score_set(self, primers):
                return 0.0

            @property
            def name(self) -> str:
                return "minimal"

        cache = PositionCache([prefix], [PRIMER])
        optimizer = _Minimal(
            cache,
            [prefix],
            [GENOME],
            config=OptimizerConfig(target_set_size=1, extension_reach=2_000),
        )
        metrics = optimizer.compute_metrics([PRIMER])

        out = compare_geometries(cache, [prefix], [GENOME], [PRIMER], reach=2_000)

        assert out["symmetric"] == pytest.approx(metrics.fg_coverage)
        assert out["directional"] < metrics.fg_coverage


def _params(**overrides):
    params = {
        "fg_prefixes": ["target"],
        "fg_seq_lengths": [1_270_000],
        "polymerase": "phi29",
    }
    params.update(overrides)
    return params


class TestADirectionalRunNamesItsOwnReach:
    def test_directional_without_a_reach_is_refused(self):
        with pytest.raises(InvalidDesignRequest) as caught:
            resolve_design_request(_params(coverage_geometry="directional"))

        message = str(caught.value)
        assert "coverage_reach" in message
        assert "4,354-6,710" in message
        assert "2026-09-23-reach-refit-directional" in message

    def test_directional_with_a_reach_is_accepted(self):
        request = resolve_design_request(
            _params(coverage_geometry="directional", coverage_reach=4400)
        )

        assert request.coverage_geometry == "directional"
        assert request.coverage_reach == 4400
        assert request.default_sources["coverage_reach"] == "request"

    def test_the_symmetric_default_still_takes_the_polymerase_reach(self):
        request = resolve_design_request(_params())

        assert request.coverage_geometry == "symmetric"
        assert request.coverage_reach == 3000
        assert request.default_sources["coverage_reach"] == "polymerase:phi29"
        assert request.default_sources["coverage_geometry"] == "default:symmetric"

    def test_symmetric_named_explicitly_is_recorded_as_a_request(self):
        request = resolve_design_request(_params(coverage_geometry="symmetric"))

        assert request.default_sources["coverage_geometry"] == "request"
        assert request.coverage_reach == 3000

    def test_an_unknown_geometry_is_refused_before_any_index_is_opened(self):
        """`coverage.site_spans` raises on it, but only mid-search."""
        with pytest.raises(InvalidDesignRequest, match="coverage_geometry"):
            resolve_design_request(_params(coverage_geometry="radial"))

    def test_the_geometry_moves_the_request_hash(self):
        baseline = resolve_design_request(_params()).request_hash
        other = resolve_design_request(
            _params(coverage_geometry="directional", coverage_reach=4400)
        ).request_hash

        assert baseline != other


def _result():
    from neoswga.core.base_optimizer import (
        OptimizationResult,
        OptimizationStatus,
        PrimerSetMetrics,
    )

    return OptimizationResult(
        primers=(PRIMER,),
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="minimal",
    )


class TestTheSummaryCarriesIt:
    def test_nothing_is_written_when_there_is_no_comparison(self, tmp_path):
        """One figure labelled as both would be worse than no figure."""
        import json

        from neoswga.core.step4_output import save_results

        result = _result()
        out = tmp_path / "step4_improved_df.csv"
        save_results(result, str(out), coverage_geometries=None)

        summary = json.loads((tmp_path / "step4_improved_df_summary.json").read_text())
        assert "coverage_geometries" not in summary

    def test_the_comparison_is_written_when_there_is_one(self, tmp_path):
        import json

        from neoswga.core.step4_output import save_results

        result = _result()
        comparison = compare_geometries(_cache(), ["t"], [GENOME], [PRIMER], reach=2_000)
        out = tmp_path / "step4_improved_df.csv"
        save_results(result, str(out), coverage_geometries=comparison)

        summary = json.loads((tmp_path / "step4_improved_df_summary.json").read_text())
        assert summary["coverage_geometries"]["inflation"] == pytest.approx(comparison["inflation"])
        assert summary["coverage_geometries"]["reported_geometry"] == "symmetric"
