"""End-to-end test: extension_reach propagates from optimizer to report.

Confirms the chain works without running the full pipeline:
    PrimerSetMetrics  ->  OptimizationResult.to_dict()  ->  summary JSON
        ->  report._load_optimizer_summary  ->  CoverageMetrics.extension_reach.

A user reading the rendered report should see the same per-primer reach value
the optimizer used to compute fg_coverage. Without this plumbing, "95%
coverage" is ambiguous (bin coverage at 10 kb? interval coverage at 3 kb?
theoretical reach at 70 kb?).
"""

import json
from pathlib import Path

from neoswga.core.base_optimizer import (
    OptimizationResult,
    OptimizationStatus,
    PrimerSetMetrics,
)


def _make_metrics(extension_reach=4000, fg_coverage=0.95):
    return PrimerSetMetrics(
        fg_coverage=fg_coverage,
        bg_coverage=0.0,
        coverage_uniformity=0.3,
        total_fg_sites=200,
        total_bg_sites=0,
        selectivity_ratio=999.0,
        mean_tm=50.0,
        tm_range=(48.0, 52.0),
        dimer_risk_score=0.05,
        mean_gap=4000.0,
        max_gap=12000.0,
        gap_gini=0.3,
        gap_entropy=8.0,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
        extension_reach=extension_reach,
    )


def test_extension_reach_round_trips_through_summary_json(tmp_path):
    metrics = _make_metrics(extension_reach=4000)
    result = OptimizationResult(
        primers=("ATCGATCGATCG", "GCTAGCTAGCTA"),
        score=0.95,
        status=OptimizationStatus.SUCCESS,
        metrics=metrics,
        iterations=10,
        optimizer_name="test",
    )

    summary_path = tmp_path / "step4_improved_df_summary.json"
    summary_path.write_text(json.dumps(result.to_dict(), default=str))

    from neoswga.core.report import metrics as report_metrics
    loaded = report_metrics._load_optimizer_summary(tmp_path)

    assert loaded["metrics"]["extension_reach"] == 4000


def test_coverage_metrics_carries_reach_after_full_load(tmp_path):
    """Construct a results dir with a step4 CSV + summary JSON and verify
    collect_pipeline_metrics surfaces the extension_reach on CoverageMetrics."""
    metrics = _make_metrics(extension_reach=3000, fg_coverage=0.88)
    result = OptimizationResult(
        primers=("ATCGATCGATCG",),
        score=0.88,
        status=OptimizationStatus.SUCCESS,
        metrics=metrics,
        iterations=5,
        optimizer_name="test",
    )
    (tmp_path / "step4_improved_df_summary.json").write_text(
        json.dumps(result.to_dict(), default=str)
    )

    # Minimal step4 CSV so collect_pipeline_metrics has primer rows.
    (tmp_path / "step4_improved_df.csv").write_text(
        "sequence,score,fg_freq,bg_freq,tm,gini,gc,fg_count,bg_count\n"
        "ATCGATCGATCG,0.88,0.001,0.0,50.0,0.3,0.5,200,0\n"
    )
    # Minimal params so coverage.total_bases is set.
    (tmp_path / "params.json").write_text(json.dumps({
        "fg_size": 5_000_000,
        "polymerase": "phi29",
    }))

    from neoswga.core.report.metrics import collect_pipeline_metrics
    pm = collect_pipeline_metrics(str(tmp_path))

    assert pm.coverage is not None
    assert pm.coverage.from_optimizer is True
    assert pm.coverage.extension_reach == 3000
    assert pm.coverage.overall_coverage == 0.88
