"""
Unit tests for neoswga.core.report.quality module.

Tests quality grading, rating functions, and recommendations.
"""

import pytest

from neoswga.core.report.quality import (
    QualityGrade,
    GradeComponent,
    QualityAssessment,
    _safe_divide,
    _rate_value,
    _score_to_grade,
    _generate_recommendation,
    calculate_quality_grade,
    format_grade_display,
    COVERAGE_THRESHOLDS,
    ENRICHMENT_THRESHOLDS,
    UNIFORMITY_THRESHOLDS,
    TM_RANGE_THRESHOLDS,
    DIMER_THRESHOLDS,
)
from neoswga.core.report.metrics import (
    PipelineMetrics,
    CoverageMetrics,
    SpecificityMetrics,
    ThermodynamicMetrics,
    UniformityMetrics,
    PrimerMetrics,
)


class TestSafeDivide:
    """Tests for _safe_divide() function."""

    def test_normal_division(self):
        """Normal division works."""
        assert _safe_divide(10.0, 2.0) == 5.0
        assert _safe_divide(1.0, 4.0) == 0.25

    def test_zero_denominator_returns_default(self):
        """Zero denominator returns default."""
        assert _safe_divide(10.0, 0.0) == 0.0
        assert _safe_divide(10.0, 0.0, default=1.0) == 1.0

    def test_very_small_denominator_returns_default(self):
        """Very small denominator treated as zero."""
        assert _safe_divide(10.0, 1e-15) == 0.0
        assert _safe_divide(10.0, 1e-11) == 0.0

    def test_custom_default(self):
        """Custom default value works."""
        assert _safe_divide(10.0, 0.0, default=-1.0) == -1.0
        assert _safe_divide(0.0, 0.0, default=999.0) == 999.0


class TestRateValue:
    """Tests for _rate_value() function."""

    def test_excellent_rating_higher_is_better(self):
        """Excellent rating for high values (higher is better)."""
        rating, score = _rate_value(1.0, COVERAGE_THRESHOLDS)
        assert rating == "Excellent"
        assert score == 1.0

        rating, score = _rate_value(0.96, COVERAGE_THRESHOLDS)
        assert rating == "Excellent"
        assert score == 1.0

    def test_good_rating_higher_is_better(self):
        """Good rating for intermediate values."""
        rating, score = _rate_value(0.90, COVERAGE_THRESHOLDS)
        assert rating == "Good"
        assert 0.75 <= score < 1.0

    def test_acceptable_rating_higher_is_better(self):
        """Acceptable rating for lower values."""
        rating, score = _rate_value(0.75, COVERAGE_THRESHOLDS)
        assert rating == "Acceptable"
        assert 0.50 <= score < 0.75

    def test_poor_rating_higher_is_better(self):
        """Poor rating for low values."""
        rating, score = _rate_value(0.55, COVERAGE_THRESHOLDS)
        assert rating == "Poor"
        assert 0.25 <= score < 0.50

    def test_critical_rating_higher_is_better(self):
        """Critical rating for very low values."""
        rating, score = _rate_value(0.3, COVERAGE_THRESHOLDS)
        assert rating == "Critical"
        assert 0.0 <= score < 0.25

    def test_excellent_rating_lower_is_better(self):
        """Excellent rating for low values (lower is better)."""
        rating, score = _rate_value(2.0, TM_RANGE_THRESHOLDS, lower_is_better=True)
        assert rating == "Excellent"
        assert score == 1.0

    def test_good_rating_lower_is_better(self):
        """Good rating for intermediate values (lower is better)."""
        rating, score = _rate_value(4.0, TM_RANGE_THRESHOLDS, lower_is_better=True)
        assert rating == "Good"
        assert 0.75 <= score < 1.0

    def test_critical_rating_lower_is_better(self):
        """Critical rating for high values (lower is better)."""
        rating, score = _rate_value(15.0, TM_RANGE_THRESHOLDS, lower_is_better=True)
        assert rating == "Critical"
        assert 0.0 <= score < 0.25

    def test_score_interpolation(self):
        """Score is interpolated within ranges."""
        # Midpoint of good range (0.85-0.95 for coverage)
        rating, score = _rate_value(0.90, COVERAGE_THRESHOLDS)
        assert rating == "Good"
        # Score should be around 0.875 (midpoint of 0.75-1.0)
        assert 0.80 <= score <= 0.95


class TestScoreToGrade:
    """Tests for _score_to_grade() function."""

    def test_grade_a(self):
        """Score >= 0.90 gets grade A."""
        assert _score_to_grade(1.0) == QualityGrade.A
        assert _score_to_grade(0.95) == QualityGrade.A
        assert _score_to_grade(0.90) == QualityGrade.A

    def test_grade_b(self):
        """Score 0.75-0.89 gets grade B."""
        assert _score_to_grade(0.89) == QualityGrade.B
        assert _score_to_grade(0.80) == QualityGrade.B
        assert _score_to_grade(0.75) == QualityGrade.B

    def test_grade_c(self):
        """Score 0.60-0.74 gets grade C."""
        assert _score_to_grade(0.74) == QualityGrade.C
        assert _score_to_grade(0.65) == QualityGrade.C
        assert _score_to_grade(0.60) == QualityGrade.C

    def test_grade_d(self):
        """Score 0.40-0.59 gets grade D."""
        assert _score_to_grade(0.59) == QualityGrade.D
        assert _score_to_grade(0.50) == QualityGrade.D
        assert _score_to_grade(0.40) == QualityGrade.D

    def test_grade_f(self):
        """Score < 0.40 gets grade F."""
        assert _score_to_grade(0.39) == QualityGrade.F
        assert _score_to_grade(0.20) == QualityGrade.F
        assert _score_to_grade(0.0) == QualityGrade.F


class TestGenerateRecommendation:
    """Tests for _generate_recommendation() function."""

    def test_grade_a_recommendation(self):
        """Grade A gets proceed to synthesis recommendation."""
        components = [
            GradeComponent("Coverage", 0.35, 0.95, 1.0, "Excellent", ""),
        ]
        short, detailed, considerations = _generate_recommendation(
            QualityGrade.A, components
        )

        assert "PROCEED TO SYNTHESIS" in short
        assert "excellent" in detailed.lower()
        assert len(considerations) == 0  # No weak components

    def test_grade_b_recommendation(self):
        """Grade B gets proceed with monitoring recommendation."""
        components = []
        short, detailed, considerations = _generate_recommendation(
            QualityGrade.B, components
        )

        assert "MONITORING" in short
        assert "good" in detailed.lower()

    def test_grade_c_recommendation(self):
        """Grade C gets optimization recommendation."""
        components = []
        short, detailed, considerations = _generate_recommendation(
            QualityGrade.C, components
        )

        assert "OPTIMIZATION" in short
        assert "acceptable" in detailed.lower()

    def test_grade_d_recommendation(self):
        """Grade D gets strong optimization recommendation."""
        components = []
        short, detailed, considerations = _generate_recommendation(
            QualityGrade.D, components
        )

        assert "RECOMMENDED" in short
        assert "significant" in detailed.lower()

    def test_grade_f_recommendation(self):
        """Grade F gets do not proceed recommendation."""
        components = []
        short, detailed, considerations = _generate_recommendation(
            QualityGrade.F, components
        )

        assert "DO NOT PROCEED" in short
        assert "critical" in detailed.lower()

    def test_weak_component_considerations(self):
        """Weak components generate considerations."""
        components = [
            GradeComponent("Coverage", 0.35, 0.50, 0.40, "Poor", ""),
            GradeComponent("Specificity", 0.30, 30.0, 0.45, "Poor", ""),
            GradeComponent("Uniformity", 0.20, 0.45, 0.55, "Acceptable", ""),
            GradeComponent("Thermodynamics", 0.10, 10.0, 0.35, "Poor", ""),
            GradeComponent("Dimer Risk", 0.05, 0.6, 0.30, "Poor", ""),
        ]

        _, _, considerations = _generate_recommendation(QualityGrade.D, components)

        # All components with score < 0.6 should generate considerations
        assert len(considerations) == 5
        assert any("Coverage" in c for c in considerations)
        assert any("Enrichment" in c for c in considerations)
        assert any("Gini" in c for c in considerations)
        assert any("Tm range" in c for c in considerations)
        assert any("dimer" in c.lower() for c in considerations)


class TestCalculateQualityGrade:
    """Tests for calculate_quality_grade() function."""

    def test_excellent_grade(self, sample_pipeline_metrics):
        """Excellent metrics get grade A."""
        # Modify to ensure excellent scores
        sample_pipeline_metrics.coverage.overall_coverage = 0.98
        sample_pipeline_metrics.specificity.enrichment_ratio = 600
        sample_pipeline_metrics.uniformity.max_gini = 0.10
        sample_pipeline_metrics.thermodynamics.tm_range = 1.0

        assessment = calculate_quality_grade(sample_pipeline_metrics)

        assert assessment.grade == QualityGrade.A
        assert assessment.composite_score >= 0.90

    def test_poor_grade(self, sample_pipeline_metrics):
        """Poor metrics get low grade."""
        sample_pipeline_metrics.coverage.overall_coverage = 0.30
        sample_pipeline_metrics.specificity.enrichment_ratio = 10
        sample_pipeline_metrics.uniformity.max_gini = 0.70
        sample_pipeline_metrics.thermodynamics.tm_range = 15.0

        assessment = calculate_quality_grade(sample_pipeline_metrics)

        assert assessment.grade in [QualityGrade.D, QualityGrade.F]
        assert assessment.composite_score < 0.60

    def test_component_weights_sum_to_one(self, sample_pipeline_metrics):
        """Component weights sum to 1.0."""
        assessment = calculate_quality_grade(sample_pipeline_metrics)

        total_weight = sum(c.weight for c in assessment.components)
        assert total_weight == pytest.approx(1.0, abs=0.001)

    def test_all_components_present(self, sample_pipeline_metrics):
        """All five components are present."""
        assessment = calculate_quality_grade(sample_pipeline_metrics)

        component_names = {c.name for c in assessment.components}
        assert component_names == {
            "Coverage",
            "Specificity",
            "Uniformity",
            "Thermodynamics",
            "Dimer Risk",
        }

    def test_empty_metrics(self, minimal_pipeline_metrics):
        """Handle empty metrics gracefully."""
        assessment = calculate_quality_grade(minimal_pipeline_metrics)

        # Should not raise, should produce valid assessment
        assert assessment.grade is not None
        assert len(assessment.components) == 5

    def test_recommendation_generated(self, sample_pipeline_metrics):
        """Recommendation is generated."""
        assessment = calculate_quality_grade(sample_pipeline_metrics)

        assert assessment.recommendation != ""
        assert assessment.recommendation_details != ""

    def test_enrichment_capped_for_scoring(self, sample_pipeline_metrics):
        """Very high enrichment is capped for scoring."""
        sample_pipeline_metrics.specificity.enrichment_ratio = 10000.0

        assessment = calculate_quality_grade(sample_pipeline_metrics)

        # Should still be scored (not overflow)
        spec_component = next(
            c for c in assessment.components if c.name == "Specificity"
        )
        assert spec_component.normalized_score == 1.0
        # But raw value should be preserved
        assert spec_component.raw_value == 10000.0

    def test_dimer_risk_from_primers(self, sample_pipeline_metrics):
        """Dimer risk is calculated from primer dimer scores."""
        # Modify primers to have significant dimer scores
        for primer in sample_pipeline_metrics.primers:
            primer.dimer_score = -8.0

        assessment = calculate_quality_grade(sample_pipeline_metrics)

        dimer_component = next(
            c for c in assessment.components if c.name == "Dimer Risk"
        )
        # -8.0 / 10.0 = 0.8 risk
        assert dimer_component.raw_value == pytest.approx(0.8, abs=0.01)


class TestFormatGradeDisplay:
    """Tests for format_grade_display() function."""

    def test_all_grades_formatted(self):
        """All grades have formatted display."""
        assert format_grade_display(QualityGrade.A) == "A (Excellent)"
        assert format_grade_display(QualityGrade.B) == "B (Good)"
        assert format_grade_display(QualityGrade.C) == "C (Acceptable)"
        assert format_grade_display(QualityGrade.D) == "D (Poor)"
        assert format_grade_display(QualityGrade.F) == "F (Critical)"


class TestThresholds:
    """Tests for threshold configurations."""

    def test_coverage_thresholds_ordered(self):
        """Coverage thresholds are in order."""
        assert COVERAGE_THRESHOLDS["excellent"] > COVERAGE_THRESHOLDS["good"]
        assert COVERAGE_THRESHOLDS["good"] > COVERAGE_THRESHOLDS["acceptable"]
        assert COVERAGE_THRESHOLDS["acceptable"] > COVERAGE_THRESHOLDS["poor"]

    def test_enrichment_thresholds_ordered(self):
        """Enrichment thresholds are in order."""
        assert ENRICHMENT_THRESHOLDS["excellent"] > ENRICHMENT_THRESHOLDS["good"]
        assert ENRICHMENT_THRESHOLDS["good"] > ENRICHMENT_THRESHOLDS["acceptable"]

    def test_tm_range_thresholds_ordered(self):
        """Tm range thresholds are in order (lower is better)."""
        assert TM_RANGE_THRESHOLDS["excellent"] < TM_RANGE_THRESHOLDS["good"]
        assert TM_RANGE_THRESHOLDS["good"] < TM_RANGE_THRESHOLDS["acceptable"]
        assert TM_RANGE_THRESHOLDS["acceptable"] < TM_RANGE_THRESHOLDS["poor"]


class TestEdgeCases:
    """Edge case tests for quality grading."""

    def test_zero_coverage(self, minimal_pipeline_metrics):
        """Handle zero coverage."""
        minimal_pipeline_metrics.coverage = CoverageMetrics(overall_coverage=0.0)

        assessment = calculate_quality_grade(minimal_pipeline_metrics)

        coverage_comp = next(
            c for c in assessment.components if c.name == "Coverage"
        )
        assert coverage_comp.rating == "Critical"

    def test_perfect_scores(self, sample_pipeline_metrics):
        """Handle perfect scores."""
        sample_pipeline_metrics.coverage.overall_coverage = 1.0
        sample_pipeline_metrics.specificity.enrichment_ratio = 1000.0
        sample_pipeline_metrics.uniformity.max_gini = 0.0  # Perfect uniformity
        sample_pipeline_metrics.thermodynamics.tm_range = 0.0
        sample_pipeline_metrics.primers = []  # No dimers

        assessment = calculate_quality_grade(sample_pipeline_metrics)

        assert assessment.grade == QualityGrade.A
        assert assessment.composite_score >= 0.95

    def test_missing_thermodynamics(self, minimal_pipeline_metrics):
        """Handle missing thermodynamics."""
        minimal_pipeline_metrics.thermodynamics = None

        assessment = calculate_quality_grade(minimal_pipeline_metrics)

        # Should use default tm_range=10.0
        thermo_comp = next(
            c for c in assessment.components if c.name == "Thermodynamics"
        )
        assert thermo_comp.raw_value == 10.0

    def test_missing_uniformity(self, minimal_pipeline_metrics):
        """Handle missing uniformity."""
        minimal_pipeline_metrics.uniformity = None

        assessment = calculate_quality_grade(minimal_pipeline_metrics)

        # Should use default gini=0.5
        uniform_comp = next(
            c for c in assessment.components if c.name == "Uniformity"
        )
        assert uniform_comp.raw_value == 0.5  # 1 - 0.5 gini

    def test_negative_dimer_scores(self, sample_pipeline_metrics):
        """Handle negative dimer scores correctly."""
        # Set very negative dimer scores
        for primer in sample_pipeline_metrics.primers:
            primer.dimer_score = -12.0

        assessment = calculate_quality_grade(sample_pipeline_metrics)

        dimer_comp = next(
            c for c in assessment.components if c.name == "Dimer Risk"
        )
        # abs(-12.0) / 10.0 = 1.2, capped to 1.0
        assert dimer_comp.raw_value == 1.0
        assert dimer_comp.rating == "Critical"
