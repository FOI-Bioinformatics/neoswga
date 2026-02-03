"""
Quality grading system for SWGA primer sets.

Provides letter grades (A-F) based on weighted metrics:
- Coverage (35%)
- Specificity (30%)
- Uniformity (20%)
- Thermodynamics (10%)
- Dimer Risk (5%)
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Tuple, Optional
import logging

from neoswga.core.report.metrics import PipelineMetrics

logger = logging.getLogger(__name__)


class QualityGrade(Enum):
    """Quality grade levels A-F."""
    A = "A"  # Excellent - ready for synthesis
    B = "B"  # Good - minor concerns
    C = "C"  # Acceptable - consider optimization
    D = "D"  # Poor - optimization recommended
    F = "F"  # Critical - do not proceed


@dataclass
class GradeComponent:
    """A single component of the quality grade."""
    name: str
    weight: float
    raw_value: float
    normalized_score: float  # 0-1 scale
    rating: str  # "Excellent", "Good", etc.
    description: str


@dataclass
class QualityAssessment:
    """Complete quality assessment with grade and components."""
    grade: QualityGrade
    composite_score: float  # 0-1 scale
    components: List[GradeComponent]
    recommendation: str
    recommendation_details: str
    considerations: List[str]


# Thresholds for each metric (values that achieve each rating level)
COVERAGE_THRESHOLDS = {
    "excellent": 0.95,
    "good": 0.85,
    "acceptable": 0.70,
    "poor": 0.50,
}

ENRICHMENT_THRESHOLDS = {
    "excellent": 500.0,
    "good": 100.0,
    "acceptable": 50.0,
    "poor": 20.0,
}

# Uniformity is measured by Gini (lower is better)
# We report 1-Gini as "uniformity score"
UNIFORMITY_THRESHOLDS = {
    "excellent": 0.85,  # Gini < 0.15
    "good": 0.70,       # Gini < 0.30
    "acceptable": 0.55, # Gini < 0.45
    "poor": 0.40,       # Gini < 0.60
}

TM_RANGE_THRESHOLDS = {
    "excellent": 3.0,
    "good": 5.0,
    "acceptable": 8.0,
    "poor": 12.0,
}

# Dimer risk thresholds (lower is better)
DIMER_THRESHOLDS = {
    "excellent": 0.1,
    "good": 0.2,
    "acceptable": 0.35,
    "poor": 0.5,
}


def _safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """Safe division that returns default on zero denominator."""
    if abs(denominator) < 1e-10:
        return default
    return numerator / denominator


def _rate_value(value: float, thresholds: Dict[str, float],
                lower_is_better: bool = False) -> Tuple[str, float]:
    """
    Rate a value against thresholds.

    Returns:
        Tuple of (rating string, normalized score 0-1)
    """
    if lower_is_better:
        if value <= thresholds["excellent"]:
            return "Excellent", 1.0
        elif value <= thresholds["good"]:
            # Interpolate between excellent and good
            denom = thresholds["good"] - thresholds["excellent"]
            score = 0.75 + 0.25 * _safe_divide(thresholds["good"] - value, denom, 0.0)
            return "Good", score
        elif value <= thresholds["acceptable"]:
            denom = thresholds["acceptable"] - thresholds["good"]
            score = 0.50 + 0.25 * _safe_divide(thresholds["acceptable"] - value, denom, 0.0)
            return "Acceptable", score
        elif value <= thresholds["poor"]:
            denom = thresholds["poor"] - thresholds["acceptable"]
            score = 0.25 + 0.25 * _safe_divide(thresholds["poor"] - value, denom, 0.0)
            return "Poor", score
        else:
            # Critical range: score decreases as value exceeds poor threshold
            denom = thresholds["poor"] if thresholds["poor"] != 0 else 1.0
            excess_ratio = _safe_divide(value - thresholds["poor"], denom, 0.0)
            return "Critical", max(0.0, 0.25 * (1.0 - excess_ratio))
    else:
        if value >= thresholds["excellent"]:
            return "Excellent", 1.0
        elif value >= thresholds["good"]:
            denom = thresholds["excellent"] - thresholds["good"]
            score = 0.75 + 0.25 * _safe_divide(value - thresholds["good"], denom, 0.0)
            return "Good", score
        elif value >= thresholds["acceptable"]:
            denom = thresholds["good"] - thresholds["acceptable"]
            score = 0.50 + 0.25 * _safe_divide(value - thresholds["acceptable"], denom, 0.0)
            return "Acceptable", score
        elif value >= thresholds["poor"]:
            denom = thresholds["acceptable"] - thresholds["poor"]
            score = 0.25 + 0.25 * _safe_divide(value - thresholds["poor"], denom, 0.0)
            return "Poor", score
        else:
            # Critical range: score based on how far below poor threshold
            denom = thresholds["poor"] if thresholds["poor"] != 0 else 1.0
            return "Critical", max(0.0, 0.25 * _safe_divide(value, denom, 0.0))


def _score_to_grade(score: float) -> QualityGrade:
    """Convert composite score to letter grade."""
    if score >= 0.90:
        return QualityGrade.A
    elif score >= 0.75:
        return QualityGrade.B
    elif score >= 0.60:
        return QualityGrade.C
    elif score >= 0.40:
        return QualityGrade.D
    else:
        return QualityGrade.F


def _generate_recommendation(
    grade: QualityGrade,
    components: List[GradeComponent],
) -> Tuple[str, str, List[str]]:
    """
    Generate recommendation based on grade and components.

    Returns:
        Tuple of (short recommendation, detailed recommendation, considerations)
    """
    considerations = []

    # Find weakest components
    weak_components = [c for c in components if c.normalized_score < 0.6]

    for comp in weak_components:
        if comp.name == "Coverage":
            considerations.append(
                f"Coverage ({comp.raw_value:.1%}) is below target. "
                "Consider adding more primers or adjusting parameters."
            )
        elif comp.name == "Specificity":
            considerations.append(
                f"Enrichment ({comp.raw_value:.0f}x) could be improved. "
                "Consider background-aware optimization."
            )
        elif comp.name == "Uniformity":
            considerations.append(
                f"Binding uniformity (Gini={1-comp.raw_value:.2f}) shows clustering. "
                "Some regions may be under-amplified."
            )
        elif comp.name == "Thermodynamics":
            considerations.append(
                f"Tm range ({comp.raw_value:.1f}C) is wide. "
                "Consider tighter Tm filtering."
            )
        elif comp.name == "Dimer Risk":
            considerations.append(
                "Primer-dimer interactions detected. "
                "Monitor amplification efficiency."
            )

    if grade == QualityGrade.A:
        short_rec = "PROCEED TO SYNTHESIS"
        detailed = (
            "This primer set demonstrates excellent performance across all metrics. "
            "The predicted coverage and specificity suggest efficient selective "
            "amplification of target DNA."
        )
    elif grade == QualityGrade.B:
        short_rec = "PROCEED WITH MONITORING"
        detailed = (
            "This primer set shows good overall quality with minor areas for "
            "potential improvement. Suitable for experimental validation."
        )
    elif grade == QualityGrade.C:
        short_rec = "CONSIDER OPTIMIZATION"
        detailed = (
            "This primer set has acceptable quality but shows moderate concerns. "
            "Consider re-running optimization or adjusting parameters before synthesis."
        )
    elif grade == QualityGrade.D:
        short_rec = "OPTIMIZATION RECOMMENDED"
        detailed = (
            "This primer set has significant quality issues. Optimization is "
            "strongly recommended before proceeding with synthesis."
        )
    else:  # F
        short_rec = "DO NOT PROCEED"
        detailed = (
            "This primer set has critical quality issues and is not suitable "
            "for SWGA. Review parameters and genome suitability."
        )

    return short_rec, detailed, considerations


def calculate_quality_grade(metrics: PipelineMetrics) -> QualityAssessment:
    """
    Calculate quality grade for a primer set.

    Uses weighted scoring:
    - Coverage: 35%
    - Specificity: 30%
    - Uniformity: 20%
    - Thermodynamics: 10%
    - Dimer Risk: 5%

    Args:
        metrics: PipelineMetrics from collect_pipeline_metrics()

    Returns:
        QualityAssessment with grade, score, and recommendations
    """
    components = []

    # Coverage (35%)
    coverage = metrics.coverage.overall_coverage if metrics.coverage else 0
    rating, score = _rate_value(coverage, COVERAGE_THRESHOLDS)
    components.append(GradeComponent(
        name="Coverage",
        weight=0.35,
        raw_value=coverage,
        normalized_score=score,
        rating=rating,
        description=f"{coverage:.1%} of target genome covered",
    ))

    # Specificity (30%)
    enrichment = metrics.specificity.enrichment_ratio if metrics.specificity else 0
    # Cap enrichment for scoring (very high values all get perfect score)
    capped_enrichment = min(enrichment, 1000)
    rating, score = _rate_value(capped_enrichment, ENRICHMENT_THRESHOLDS)
    components.append(GradeComponent(
        name="Specificity",
        weight=0.30,
        raw_value=enrichment,
        normalized_score=score,
        rating=rating,
        description=f"{enrichment:.0f}x target/background enrichment",
    ))

    # Uniformity (20%)
    # Convert Gini to uniformity score (1 - Gini)
    gini = metrics.uniformity.max_gini if metrics.uniformity else 0.5
    uniformity = 1.0 - gini
    rating, score = _rate_value(uniformity, UNIFORMITY_THRESHOLDS)
    components.append(GradeComponent(
        name="Uniformity",
        weight=0.20,
        raw_value=uniformity,
        normalized_score=score,
        rating=rating,
        description=f"Gini index {gini:.2f} (lower is more uniform)",
    ))

    # Thermodynamics (10%)
    tm_range = metrics.thermodynamics.tm_range if metrics.thermodynamics else 10.0
    rating, score = _rate_value(tm_range, TM_RANGE_THRESHOLDS, lower_is_better=True)
    components.append(GradeComponent(
        name="Thermodynamics",
        weight=0.10,
        raw_value=tm_range,
        normalized_score=score,
        rating=rating,
        description=f"Tm range {tm_range:.1f}C across primers",
    ))

    # Dimer Risk (5%)
    # Use normalized dimer score (0-1 scale)
    # Note: dimer scores are typically negative (kcal/mol), more negative = worse
    dimer_risk = 0.0
    if metrics.primers:
        # Use min() to get the most negative (worst) dimer score
        worst_dimer = min((p.dimer_score for p in metrics.primers), default=0)
        # Normalize: assume -10 kcal/mol is significant risk
        dimer_risk = min(abs(worst_dimer) / 10.0, 1.0)
    rating, score = _rate_value(dimer_risk, DIMER_THRESHOLDS, lower_is_better=True)
    components.append(GradeComponent(
        name="Dimer Risk",
        weight=0.05,
        raw_value=dimer_risk,
        normalized_score=score,
        rating=rating,
        description="Primer-primer interaction risk",
    ))

    # Calculate composite score
    composite = sum(c.weight * c.normalized_score for c in components)

    # Get grade
    grade = _score_to_grade(composite)

    # Generate recommendation
    short_rec, detailed_rec, considerations = _generate_recommendation(
        grade, components
    )

    return QualityAssessment(
        grade=grade,
        composite_score=composite,
        components=components,
        recommendation=short_rec,
        recommendation_details=detailed_rec,
        considerations=considerations,
    )


def format_grade_display(grade: QualityGrade) -> str:
    """Format grade for display with description."""
    descriptions = {
        QualityGrade.A: "Excellent",
        QualityGrade.B: "Good",
        QualityGrade.C: "Acceptable",
        QualityGrade.D: "Poor",
        QualityGrade.F: "Critical",
    }
    return f"{grade.value} ({descriptions[grade]})"
