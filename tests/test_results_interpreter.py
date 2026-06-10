"""Tests for neoswga.core.results_interpreter module."""

import csv
import json
import pytest

from neoswga.core.results_interpreter import (
    QualityRating,
    MetricAssessment,
    ResultsReport,
    ResultsInterpreter,
    rate_metric,
    interpret_results,
    COVERAGE_THRESHOLDS,
    ENRICHMENT_THRESHOLDS,
    UNIFORMITY_THRESHOLDS,
    DIMER_SCORE_THRESHOLDS,
)


# ---------------------------------------------------------------------------
# QualityRating enum
# ---------------------------------------------------------------------------

class TestQualityRating:
    def test_enum_values(self):
        assert QualityRating.EXCELLENT.value == "EXCELLENT"
        assert QualityRating.GOOD.value == "GOOD"
        assert QualityRating.ACCEPTABLE.value == "ACCEPTABLE"
        assert QualityRating.POOR.value == "POOR"
        assert QualityRating.CRITICAL.value == "CRITICAL"

    def test_enum_member_count(self):
        assert len(QualityRating) == 5


# ---------------------------------------------------------------------------
# Dataclass construction
# ---------------------------------------------------------------------------

class TestDataclasses:
    def test_metric_assessment_creation(self):
        ma = MetricAssessment(
            name="Coverage",
            value=0.9,
            rating=QualityRating.GOOD,
            unit="%",
            context="test context",
            threshold_info="test thresholds",
        )
        assert ma.name == "Coverage"
        assert ma.value == 0.9
        assert ma.rating == QualityRating.GOOD

    def test_results_report_creation(self):
        report = ResultsReport(
            primer_count=8,
            assessments=[],
            overall_rating=QualityRating.GOOD,
            recommendation="Go ahead",
            next_steps=["Order primers"],
            warnings=[],
        )
        assert report.primer_count == 8
        assert report.enrichment_estimate is None

    def test_results_report_with_enrichment_estimate(self):
        report = ResultsReport(
            primer_count=6,
            assessments=[],
            overall_rating=QualityRating.EXCELLENT,
            recommendation="Ready",
            next_steps=[],
            warnings=[],
            enrichment_estimate={"estimated_fold_enrichment": 500},
        )
        assert report.enrichment_estimate["estimated_fold_enrichment"] == 500


# ---------------------------------------------------------------------------
# rate_metric (higher-is-better and lower-is-better)
# ---------------------------------------------------------------------------

class TestRateMetric:
    # Coverage thresholds: EXCELLENT >= 0.95, GOOD >= 0.85, ACCEPTABLE >= 0.70, POOR >= 0.50
    def test_higher_is_better_excellent(self):
        assert rate_metric(0.96, COVERAGE_THRESHOLDS) == QualityRating.EXCELLENT

    def test_higher_is_better_good(self):
        assert rate_metric(0.90, COVERAGE_THRESHOLDS) == QualityRating.GOOD

    def test_higher_is_better_acceptable(self):
        assert rate_metric(0.75, COVERAGE_THRESHOLDS) == QualityRating.ACCEPTABLE

    def test_higher_is_better_poor(self):
        assert rate_metric(0.55, COVERAGE_THRESHOLDS) == QualityRating.POOR

    def test_higher_is_better_critical(self):
        assert rate_metric(0.30, COVERAGE_THRESHOLDS) == QualityRating.CRITICAL

    def test_higher_is_better_boundary_excellent(self):
        assert rate_metric(0.95, COVERAGE_THRESHOLDS) == QualityRating.EXCELLENT

    # Uniformity thresholds (lower is better): EXCELLENT <= 0.3, GOOD <= 0.45, etc.
    def test_lower_is_better_excellent(self):
        assert rate_metric(0.2, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.EXCELLENT

    def test_lower_is_better_good(self):
        assert rate_metric(0.4, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.GOOD

    def test_lower_is_better_acceptable(self):
        assert rate_metric(0.55, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.ACCEPTABLE

    def test_lower_is_better_poor(self):
        assert rate_metric(0.7, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.POOR

    def test_lower_is_better_critical(self):
        assert rate_metric(0.9, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.CRITICAL

    def test_lower_is_better_boundary_excellent(self):
        assert rate_metric(0.3, UNIFORMITY_THRESHOLDS, lower_is_better=True) == QualityRating.EXCELLENT

    # Enrichment thresholds
    def test_enrichment_excellent(self):
        assert rate_metric(250, ENRICHMENT_THRESHOLDS) == QualityRating.EXCELLENT

    def test_enrichment_critical(self):
        assert rate_metric(5, ENRICHMENT_THRESHOLDS) == QualityRating.CRITICAL

    # Dimer score thresholds (lower is better)
    def test_dimer_excellent(self):
        assert rate_metric(0.05, DIMER_SCORE_THRESHOLDS, lower_is_better=True) == QualityRating.EXCELLENT

    def test_dimer_poor(self):
        assert rate_metric(0.45, DIMER_SCORE_THRESHOLDS, lower_is_better=True) == QualityRating.POOR


# ---------------------------------------------------------------------------
# Helper: create a step4 CSV file with controllable columns
# ---------------------------------------------------------------------------

def _write_step4_csv(path, rows):
    """Write a step4_improved_df.csv with given rows (list of dicts)."""
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _make_primer_rows(
    n=6,
    coverage=0.92,
    fg_freq=1e-3,
    bg_freq=1e-5,
    gini=0.35,
    dimer_score=0.1,
    primer_prefix="ATCGATCG",
):
    """Generate n primer rows with sensible defaults."""
    rows = []
    for i in range(n):
        rows.append({
            "primer": f"{primer_prefix}{i:02d}"[-len(primer_prefix):],
            "fg_freq": str(fg_freq),
            "bg_freq": str(bg_freq),
            "gini": str(gini),
            "dimer_score": str(dimer_score),
            "coverage": str(coverage),
        })
    return rows


# ---------------------------------------------------------------------------
# ResultsInterpreter with synthetic data
# ---------------------------------------------------------------------------

class TestResultsInterpreter:
    def test_analyze_with_step4(self, tmp_path):
        rows = _make_primer_rows(n=6, coverage=0.92, gini=0.35, dimer_score=0.1)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        interpreter = ResultsInterpreter(str(tmp_path))
        report = interpreter.analyze()

        assert report.primer_count == 6
        assert isinstance(report.overall_rating, QualityRating)
        assert len(report.assessments) > 0

    def test_analyze_with_step3_fallback(self, tmp_path):
        rows = _make_primer_rows(n=4, coverage=0.80)
        # Write as step3_df.csv (no step4 present)
        _write_step4_csv(tmp_path / "step3_df.csv", rows)

        interpreter = ResultsInterpreter(str(tmp_path))
        report = interpreter.analyze()

        assert report.primer_count == 4

    def test_analyze_no_files_raises(self, tmp_path):
        interpreter = ResultsInterpreter(str(tmp_path))
        with pytest.raises(FileNotFoundError, match="No results found"):
            interpreter.analyze()

    def test_coverage_assessment(self, tmp_path):
        rows = _make_primer_rows(coverage=0.98)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        cov = [a for a in report.assessments if a.name == "Genome Coverage"]
        assert len(cov) == 1
        assert cov[0].rating == QualityRating.EXCELLENT
        assert cov[0].value == pytest.approx(0.98)

    def test_enrichment_from_fg_bg_freq(self, tmp_path):
        rows = _make_primer_rows(fg_freq=1e-3, bg_freq=1e-5)
        # Remove explicit coverage/enrichment so enrichment is calculated from freq
        for r in rows:
            r.pop("coverage", None)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        enr = [a for a in report.assessments if a.name == "Enrichment Ratio"]
        assert len(enr) == 1
        assert enr[0].value == pytest.approx(100.0)  # 1e-3 / 1e-5

    def test_uniformity_assessment(self, tmp_path):
        rows = _make_primer_rows(gini=0.55)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        uni = [a for a in report.assessments if a.name == "Binding Uniformity"]
        assert len(uni) == 1
        assert uni[0].rating == QualityRating.ACCEPTABLE

    def test_dimer_assessment(self, tmp_path):
        rows = _make_primer_rows(dimer_score=0.4)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        dim = [a for a in report.assessments if a.name == "Dimer Risk"]
        assert len(dim) == 1
        assert dim[0].rating == QualityRating.POOR

    def test_overall_excellent(self, tmp_path):
        rows = _make_primer_rows(coverage=0.98, gini=0.2, dimer_score=0.05,
                                 fg_freq=1e-2, bg_freq=1e-5)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.overall_rating == QualityRating.EXCELLENT

    def test_overall_poor(self, tmp_path):
        rows = _make_primer_rows(coverage=0.40, gini=0.8, dimer_score=0.6,
                                 fg_freq=1e-5, bg_freq=1e-5)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.overall_rating in (QualityRating.POOR, QualityRating.CRITICAL)

    def test_warnings_generated_for_poor_metrics(self, tmp_path):
        rows = _make_primer_rows(coverage=0.40, gini=0.8, dimer_score=0.6)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert len(report.warnings) > 0
        warning_text = " ".join(report.warnings)
        assert "coverage" in warning_text.lower() or "gap" in warning_text.lower()

    def test_no_assessments_yield_poor_overall(self, tmp_path):
        """If CSV has no recognized metric columns, overall should be POOR."""
        rows = [{"primer": "ATCGATCG", "other_col": "123"}]
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.overall_rating == QualityRating.POOR
        assert report.primer_count == 1


# ---------------------------------------------------------------------------
# Recommendation generation
# ---------------------------------------------------------------------------

class TestRecommendations:
    def test_excellent_recommendation(self, tmp_path):
        rows = _make_primer_rows(coverage=0.98, gini=0.2, dimer_score=0.05,
                                 fg_freq=1e-2, bg_freq=1e-5)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert "synthesis" in report.recommendation.lower() or "ready" in report.recommendation.lower()

    def test_critical_recommendation(self, tmp_path):
        rows = _make_primer_rows(coverage=0.10, gini=0.9, dimer_score=0.8,
                                 fg_freq=1e-6, bg_freq=1e-5)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert "not proceed" in report.recommendation.lower() or "issues" in report.recommendation.lower()

    def test_next_steps_not_empty(self, tmp_path):
        rows = _make_primer_rows()
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert len(report.next_steps) > 0


# ---------------------------------------------------------------------------
# print_report (capsys)
# ---------------------------------------------------------------------------

class TestPrintReport:
    def test_print_report_does_not_crash(self, tmp_path, capsys):
        rows = _make_primer_rows(coverage=0.92, gini=0.4, dimer_score=0.15)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        interpreter = ResultsInterpreter(str(tmp_path))
        report = interpreter.analyze()
        interpreter.print_report(report)

        captured = capsys.readouterr()
        assert "PRIMER SET QUALITY ASSESSMENT" in captured.out
        assert "Overall rating:" in captured.out
        assert "Recommendation" in captured.out

    def test_print_report_with_warnings(self, tmp_path, capsys):
        rows = _make_primer_rows(coverage=0.30, gini=0.85, dimer_score=0.7)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        interpreter = ResultsInterpreter(str(tmp_path))
        report = interpreter.analyze()
        interpreter.print_report(report)

        captured = capsys.readouterr()
        assert "Warnings" in captured.out

    def test_print_report_with_enrichment_estimate(self, tmp_path, capsys):
        """Print a report that has an enrichment estimate dict."""
        report = ResultsReport(
            primer_count=6,
            assessments=[],
            overall_rating=QualityRating.GOOD,
            recommendation="Ready",
            next_steps=["Order"],
            warnings=[],
            enrichment_estimate={
                "estimated_fold_enrichment": 350.0,
                "mean_amplification_factor": 0.55,
                "polymerase": "phi29",
                "reaction_temp": 30.0,
                "confidence": "moderate",
            },
        )
        interpreter = ResultsInterpreter(str(tmp_path))
        interpreter.print_report(report)

        captured = capsys.readouterr()
        assert "Predicted Enrichment" in captured.out
        assert "350" in captured.out


# ---------------------------------------------------------------------------
# interpret_results convenience function
# ---------------------------------------------------------------------------

class TestInterpretResults:
    def test_interpret_results_verbose(self, tmp_path, capsys):
        rows = _make_primer_rows()
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = interpret_results(str(tmp_path), verbose=True)
        assert isinstance(report, ResultsReport)
        captured = capsys.readouterr()
        assert "PRIMER SET QUALITY ASSESSMENT" in captured.out

    def test_interpret_results_quiet(self, tmp_path, capsys):
        rows = _make_primer_rows()
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = interpret_results(str(tmp_path), verbose=False)
        assert isinstance(report, ResultsReport)
        captured = capsys.readouterr()
        assert captured.out == ""


# ---------------------------------------------------------------------------
# Graceful handling of edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_empty_csv_file(self, tmp_path):
        """Empty CSV (header only, no rows) should yield primer_count 0."""
        csv_path = tmp_path / "step4_improved_df.csv"
        csv_path.write_text("primer,fg_freq,bg_freq,gini\n")

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.primer_count == 0
        assert report.overall_rating == QualityRating.POOR

    def test_missing_columns(self, tmp_path):
        """CSV with no recognized metric columns still produces a report."""
        rows = [{"primer": "ATCG", "name": "test"}]
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.primer_count == 1
        assert len(report.assessments) == 0

    def test_zero_bg_freq_enrichment(self, tmp_path):
        """Zero background frequency should not cause division-by-zero crash."""
        rows = [{"primer": "ATCG", "fg_freq": "0.001", "bg_freq": "0"}]
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        enr = [a for a in report.assessments if a.name == "Enrichment Ratio"]
        # Should handle gracefully (either inf or skip)
        if enr:
            assert enr[0].value == float("inf")

    def test_enrichment_estimate_without_params_json(self, tmp_path):
        """No params.json => enrichment_estimate is None."""
        rows = _make_primer_rows()
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.enrichment_estimate is None

    def test_enrichment_estimate_with_params_json(self, tmp_path):
        """With a valid params.json, enrichment estimate may be produced."""
        rows = _make_primer_rows()
        _write_step4_csv(tmp_path / "step4_improved_df.csv", rows)

        params = {
            "polymerase": "phi29",
            "reaction_temp": 30.0,
            "na_conc": 50.0,
            "mg_conc": 2.5,
            "fg_gc": 0.5,
        }
        with open(tmp_path / "params.json", "w") as f:
            json.dump(params, f)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        # May or may not succeed depending on MechanisticModel availability;
        # at minimum it should not crash.
        if report.enrichment_estimate is not None:
            assert "estimated_fold_enrichment" in report.enrichment_estimate
            assert report.enrichment_estimate["polymerase"] == "phi29"

    def test_step4_preferred_over_step3(self, tmp_path):
        """When both step3 and step4 exist, step4 is used."""
        step4_rows = _make_primer_rows(n=5, coverage=0.95)
        step3_rows = _make_primer_rows(n=10, coverage=0.60)
        _write_step4_csv(tmp_path / "step4_improved_df.csv", step4_rows)
        _write_step4_csv(tmp_path / "step3_df.csv", step3_rows)

        report = ResultsInterpreter(str(tmp_path)).analyze()
        assert report.primer_count == 5  # step4 used, not step3
