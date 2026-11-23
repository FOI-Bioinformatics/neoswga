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
