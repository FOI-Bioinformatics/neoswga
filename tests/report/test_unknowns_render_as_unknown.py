"""The report shows unmeasured figures and sentinels as what they are.

Three defects (2026-11-23 sweep, tier 2):
- A missing or NaN per-primer Gini became 0.0, the best evenness, so the
  primer's uniformity score was 1.0. Known Issue 9: Gini is NaN below
  `min_gini_sites`.
- The technical report printed the `MAX_SELECTIVITY` sentinel (1e6, "no
  background binding") as "1000000x" enrichment, where the executive summary
  and the grade already say "no background binding".
- The strand section showed the genome-length placeholder as a measured gap for
  a host with no sites, where the strand ratio is None rather than 0.0.
"""

from neoswga.core.report.metrics import PipelineMetrics, PrimerMetrics, SpecificityMetrics
from neoswga.core.report.quality import calculate_quality_grade
from neoswga.core.report.strand_section import _notes
from neoswga.core.report.technical_report import _contribution_score, _render_summary_metrics
from neoswga.core.selectivity import MAX_SELECTIVITY


def test_a_nan_gini_is_unknown():
    primer = PrimerMetrics.from_row({"primer": "ACGTACGTAC", "gini": "nan"})
    assert primer.gini is None
    assert PrimerMetrics.from_row({"primer": "ACGTACGTAC"}).gini is None
    assert PrimerMetrics.from_row({"primer": "ACGTACGTAC", "gini": "0.0"}).gini == 0.0


def test_an_unknown_gini_does_not_score_as_perfectly_even():
    measured_even = PrimerMetrics.from_row({"primer": "ACGTACGTAC", "gini": "0.0"})
    unknown = PrimerMetrics.from_row({"primer": "ACGTACGTAC", "gini": "nan"})
    assert _contribution_score(unknown, amp_available=False) != _contribution_score(
        measured_even, amp_available=False
    )


def test_the_sentinel_is_not_printed_as_an_enrichment():
    metrics = PipelineMetrics(primer_count=6)
    metrics.specificity = SpecificityMetrics(enrichment_ratio=MAX_SELECTIVITY)
    html = _render_summary_metrics(calculate_quality_grade(metrics))
    assert "1000000x" not in html
    assert "no background binding" in html


def test_a_host_with_no_sites_gets_no_gap_note():
    stats = {"pLTR": {"strand_coverage_ratio": None, "strand_alternation_gap_max": 6258}}
    assert _notes(stats, {"pLTR": "host"}, 3000) == ""


def test_a_measured_gap_still_gets_its_note():
    stats = {"pcDNA": {"strand_coverage_ratio": 0.5, "strand_alternation_gap_max": 9000}}
    assert "widest stretch" in _notes(stats, {"pcDNA": "target"}, 3000)


def _unmeasured_uniformity_metrics():
    from neoswga.core.report.metrics import _calculate_uniformity_metrics

    primers = [
        PrimerMetrics.from_row({"primer": "ACGTACGTAC", "gini": "nan"}),
        PrimerMetrics.from_row({"primer": "TTGACCATGA", "gini": "nan"}),
    ]
    metrics = PipelineMetrics(primer_count=2)
    metrics.primers = primers
    metrics.uniformity = _calculate_uniformity_metrics(primers)
    return metrics


def test_the_set_gini_is_unknown_when_no_primer_measured_one():
    """Deferred tier-2 minor: the set's mean and max Gini fell to 0.0, the
    perfectly even value, when no primer had enough sites for one."""
    uniformity = _unmeasured_uniformity_metrics().uniformity
    assert uniformity.mean_gini is None and uniformity.max_gini is None


def test_the_grade_leaves_an_unmeasured_uniformity_out():
    quality = calculate_quality_grade(_unmeasured_uniformity_metrics())
    uniformity = next(c for c in quality.components if c.name == "Uniformity")
    assert uniformity.rating == "Not measured"
    measured = [c for c in quality.components if c.name != "Uniformity"]
    expected = sum(c.weight * c.normalized_score for c in measured) / sum(
        c.weight for c in measured
    )
    assert abs(quality.composite_score - expected) < 1e-9


def test_both_reports_render_an_unmeasured_uniformity():
    from neoswga.core.report.executive_summary import ExecutiveSummary, render_executive_summary
    from neoswga.core.report.technical_report import TechnicalReportData, render_technical_report

    metrics = _unmeasured_uniformity_metrics()
    quality = calculate_quality_grade(metrics)
    summary = render_executive_summary(ExecutiveSummary(metrics, quality, "now"))
    technical = render_technical_report(
        TechnicalReportData(generated_at="now", metrics=metrics, quality=quality)
    )
    for html in (summary, technical):
        assert "nan" not in html.lower().replace("nanomolar", "")
        assert "not measured" in html.lower()


def test_the_json_report_writes_an_unmeasured_figure_as_null():
    """Final review: the JSON report wrote NaN, which `jq` and JSON.parse
    reject, beside a normalized score of 0.0 that read as a measured zero."""
    import json

    from neoswga.cli.report import _component_record

    quality = calculate_quality_grade(_unmeasured_uniformity_metrics())
    records = [_component_record(c) for c in quality.components]
    json.dumps(records, allow_nan=False)
    uniformity = next(r for r in records if r["name"] == "Uniformity")
    assert uniformity["raw_value"] is None and uniformity["measured"] is False
    assert uniformity["normalized_score"] is None
