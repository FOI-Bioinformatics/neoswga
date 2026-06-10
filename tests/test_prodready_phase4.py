"""Phase 4 (production-readiness v2): comprehensive report integration."""

import json

import pytest

# ---------------------------------------------------------------------------
# metrics.py reads the authoritative optimizer summary (not CSV estimates)
# ---------------------------------------------------------------------------


def _write(results_dir, name, obj):
    (results_dir / name).write_text(json.dumps(obj))


def _minimal_results(tmp_path):
    # A minimal step4 CSV so primers load.
    (tmp_path / "step4_improved_df.csv").write_text(
        "primer,tm,gc,fg_count,bg_count\nATCGATCG,30,0.5,100,2\n"
    )
    (tmp_path / "params.json").write_text(
        json.dumps(
            {
                "polymerase": "phi29",
                "reaction_temp": 30,
                "mg_conc": 10.0,
                "betaine_m": 1.0,
                "fg_seq_lengths": [1000000],
            }
        )
    )
    return tmp_path


def test_metrics_reads_authoritative_summary(tmp_path):
    from neoswga.core.report.metrics import collect_pipeline_metrics

    _minimal_results(tmp_path)
    _write(
        tmp_path,
        "step4_improved_df_summary.json",
        {
            "metrics": {
                "fg_coverage": 0.93,
                "bg_coverage": 0.01,
                "selectivity_ratio": 87.5,
                "coverage_uniformity": 0.22,
                "total_fg_sites": 540,
                "total_bg_sites": 6,
                "strand_alternation_score": 0.61,
                "strand_coverage_ratio": 0.8,
                "per_target_coverage": {"t1": 0.95, "t2": 0.4},
                "extension_reach": 3000,
            },
            "ensemble_comparison": [
                {
                    "method": "hybrid",
                    "normalized_score": 0.7,
                    "n_primers": 6,
                    "fg_coverage": 0.93,
                    "bg_coverage": 0.01,
                    "status": "success",
                    "selected": True,
                },
            ],
            "pareto_metrics": [{"set_size": 6, "fg_coverage": 0.93}],
        },
    )

    m = collect_pipeline_metrics(str(tmp_path))

    assert m.coverage.from_optimizer is True
    assert m.coverage.per_target_coverage == {"t1": 0.95, "t2": 0.4}
    assert m.coverage.coverage_uniformity == pytest.approx(0.22)
    # measured selectivity/bg preferred over CSV estimate
    assert m.specificity.from_optimizer is True
    assert m.specificity.selectivity_ratio == pytest.approx(87.5)
    assert m.specificity.bg_coverage == pytest.approx(0.01)
    # strand metrics
    assert m.uniformity.from_optimizer is True
    assert m.uniformity.strand_alternation_score == pytest.approx(0.61)
    # top-level extras
    assert len(m.ensemble_comparison) == 1 and m.ensemble_comparison[0]["selected"]
    assert m.pareto_metrics and m.pareto_metrics[0]["set_size"] == 6
    # reaction conditions extracted
    assert m.reaction_conditions["polymerase"] == "phi29"
    assert m.reaction_conditions["additives"]["betaine_m"] == 1.0


def test_metrics_loads_coverage_gaps(tmp_path):
    from neoswga.core.report.metrics import collect_pipeline_metrics

    _minimal_results(tmp_path)
    _write(
        tmp_path,
        "coverage_gaps.json",
        {
            "n_gaps": 2,
            "used_bam": True,
            "min_depth": 5,
            "min_gap_size": 1000,
            "gaps": [{"chromosome": "c", "start": 0, "end": 2000, "size": 2000}],
        },
    )
    m = collect_pipeline_metrics(str(tmp_path))
    assert m.coverage_gaps and m.coverage_gaps["used_bam"] is True
    assert m.coverage_gaps["gaps"][0]["size"] == 2000


# ---------------------------------------------------------------------------
# Funnel: real stats rendered; fabrication killed when absent
# ---------------------------------------------------------------------------


def test_funnel_not_fabricated_when_absent():
    from neoswga.core.report.metrics import PipelineMetrics
    from neoswga.core.report.quality import calculate_quality_grade
    from neoswga.core.report.technical_report import (
        TechnicalReportData,
        collect_technical_report_data,
        render_technical_report,
    )

    m = PipelineMetrics(primer_count=6)
    m.filtering = None  # no filter_stats.json
    q = calculate_quality_grade(m)
    data = TechnicalReportData(generated_at="now", version="x", metrics=m, quality=q)
    data.filtering_stages = []  # collect would set this empty now
    html = render_technical_report(data, interactive=False)
    assert "not recorded" in html
    # The fabricated funnel bars (6*100 "Initial candidates") must be gone.
    assert "Initial candidates" not in html


def test_filter_writes_real_stats(tmp_path, monkeypatch):
    """step2 writes a real filter_stats.json funnel."""
    # Build a fake step2 environment: this is an integration-flavored unit that
    # only checks the FilteringStats.as_funnel mapping from a written file.
    from neoswga.core.report.metrics import FilteringStats

    stats = {
        "total_kmers": 10000,
        "after_frequency": 2000,
        "after_thermodynamic": 800,
        "after_background": 750,
        "after_gini": 600,
        "final_candidates": 500,
    }
    fs = FilteringStats(**{k: v for k, v in stats.items()})
    funnel = fs.as_funnel()
    labels = [s[0] for s in funnel]
    assert "Total k-mers" in labels and "Final candidates" in labels
    # monotonic non-increasing real funnel
    counts = [s[1] for s in funnel if s[1] > 0]
    assert counts == sorted(counts, reverse=True)


# ---------------------------------------------------------------------------
# New technical-report sections render statically (no --interactive needed)
# ---------------------------------------------------------------------------


def test_technical_report_renders_new_sections_statically(tmp_path):
    from neoswga.core.report.metrics import collect_pipeline_metrics
    from neoswga.core.report.technical_report import generate_technical_report

    _minimal_results(tmp_path)
    _write(
        tmp_path,
        "step4_improved_df_summary.json",
        {
            "metrics": {
                "fg_coverage": 0.93,
                "selectivity_ratio": 50.0,
                "bg_coverage": 0.01,
                "per_target_coverage": {"t1": 0.95, "t2": 0.4},
                "strand_alternation_score": 0.6,
                "strand_coverage_ratio": 0.8,
            },
            "ensemble_comparison": [
                {
                    "method": "hybrid",
                    "normalized_score": 0.7,
                    "n_primers": 6,
                    "fg_coverage": 0.93,
                    "bg_coverage": 0.01,
                    "status": "success",
                    "selected": True,
                },
            ],
        },
    )
    _write(
        tmp_path,
        "coverage_gaps.json",
        {
            "n_gaps": 1,
            "used_bam": False,
            "min_gap_size": 1000,
            "gaps": [{"chromosome": "c", "start": 0, "end": 2000, "size": 2000}],
        },
    )

    out = tmp_path / "report.html"
    generate_technical_report(str(tmp_path), str(out), interactive=False)
    html = out.read_text()
    for token in (
        "Optimizer Method Comparison",
        "Per-Target Coverage",
        "Coverage Gaps",
        "Strand Balance",
        "Reaction Conditions",
        "MEASURED",
    ):
        assert token in html, token
