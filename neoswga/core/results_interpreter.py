"""
Results interpretation for neoswga primer design output.

Analyzes pipeline output files and provides:
1. Quality ratings for coverage, enrichment, uniformity
2. Comparison to typical SWGA performance
3. Actionable recommendations
4. Go/no-go decision support
"""

import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class QualityRating(Enum):
    """Quality rating levels."""
    EXCELLENT = "EXCELLENT"
    GOOD = "GOOD"
    ACCEPTABLE = "ACCEPTABLE"
    POOR = "POOR"
    CRITICAL = "CRITICAL"


@dataclass
class MetricAssessment:
    """Assessment of a single metric."""
    name: str
    value: float
    rating: QualityRating
    unit: str
    context: str
    threshold_info: str


@dataclass
class ResultsReport:
    """Complete results interpretation report."""
    primer_count: int
    assessments: List[MetricAssessment]
    overall_rating: QualityRating
    recommendation: str
    next_steps: List[str]
    warnings: List[str]


# Thresholds for quality ratings
COVERAGE_THRESHOLDS = {
    QualityRating.EXCELLENT: 0.95,
    QualityRating.GOOD: 0.85,
    QualityRating.ACCEPTABLE: 0.70,
    QualityRating.POOR: 0.50,
}

ENRICHMENT_THRESHOLDS = {
    QualityRating.EXCELLENT: 200.0,
    QualityRating.GOOD: 100.0,
    QualityRating.ACCEPTABLE: 50.0,
    QualityRating.POOR: 20.0,
}

UNIFORMITY_THRESHOLDS = {  # Gini index (lower is better)
    QualityRating.EXCELLENT: 0.3,
    QualityRating.GOOD: 0.45,
    QualityRating.ACCEPTABLE: 0.6,
    QualityRating.POOR: 0.75,
}

DIMER_SCORE_THRESHOLDS = {  # Lower is better
    QualityRating.EXCELLENT: 0.1,
    QualityRating.GOOD: 0.2,
    QualityRating.ACCEPTABLE: 0.35,
    QualityRating.POOR: 0.5,
}


def rate_metric(value: float, thresholds: Dict[QualityRating, float],
                lower_is_better: bool = False) -> QualityRating:
    """
    Rate a metric value against thresholds.

    Args:
        value: The metric value
        thresholds: Dictionary mapping ratings to threshold values
        lower_is_better: If True, lower values get better ratings

    Returns:
        QualityRating for the value
    """
    if lower_is_better:
        if value <= thresholds[QualityRating.EXCELLENT]:
            return QualityRating.EXCELLENT
        elif value <= thresholds[QualityRating.GOOD]:
            return QualityRating.GOOD
        elif value <= thresholds[QualityRating.ACCEPTABLE]:
            return QualityRating.ACCEPTABLE
        elif value <= thresholds[QualityRating.POOR]:
            return QualityRating.POOR
        else:
            return QualityRating.CRITICAL
    else:
        if value >= thresholds[QualityRating.EXCELLENT]:
            return QualityRating.EXCELLENT
        elif value >= thresholds[QualityRating.GOOD]:
            return QualityRating.GOOD
        elif value >= thresholds[QualityRating.ACCEPTABLE]:
            return QualityRating.ACCEPTABLE
        elif value >= thresholds[QualityRating.POOR]:
            return QualityRating.POOR
        else:
            return QualityRating.CRITICAL


class ResultsInterpreter:
    """
    Interprets neoswga pipeline results.

    Usage:
        interpreter = ResultsInterpreter('results/')
        report = interpreter.analyze()
        interpreter.print_report(report)
    """

    def __init__(self, results_dir: str):
        """
        Initialize interpreter with results directory.

        Args:
            results_dir: Path to directory containing pipeline output
        """
        self.results_dir = Path(results_dir)
        self.step4_file = self.results_dir / 'step4_improved_df.csv'
        self.step3_file = self.results_dir / 'step3_df.csv'

    def analyze(self) -> ResultsReport:
        """
        Analyze pipeline results and generate report.

        Returns:
            ResultsReport with quality assessments
        """
        # Load results
        if self.step4_file.exists():
            primers = self._load_step4_results()
        elif self.step3_file.exists():
            primers = self._load_step3_results()
            logger.warning("Using step3 results (optimize step not run)")
        else:
            raise FileNotFoundError(
                f"No results found in {self.results_dir}. "
                "Run the pipeline first."
            )

        assessments = []
        warnings = []

        # Assess coverage
        coverage = self._calculate_coverage(primers)
        if coverage is not None:
            rating = rate_metric(coverage, COVERAGE_THRESHOLDS)
            assessments.append(MetricAssessment(
                name="Genome Coverage",
                value=coverage,
                rating=rating,
                unit="%",
                context=f"Percentage of target genome covered by primer binding",
                threshold_info=f">95% excellent, >85% good, >70% acceptable"
            ))
            if rating in (QualityRating.POOR, QualityRating.CRITICAL):
                warnings.append("Low coverage may result in amplification gaps")

        # Assess enrichment
        enrichment = self._calculate_enrichment(primers)
        if enrichment is not None:
            rating = rate_metric(enrichment, ENRICHMENT_THRESHOLDS)
            assessments.append(MetricAssessment(
                name="Enrichment Ratio",
                value=enrichment,
                rating=rating,
                unit="x",
                context="Expected fold-enrichment of target vs background",
                threshold_info=">200x excellent, >100x good, >50x acceptable"
            ))
            if rating in (QualityRating.POOR, QualityRating.CRITICAL):
                warnings.append("Low enrichment may result in high background")

        # Assess uniformity
        uniformity = self._calculate_uniformity(primers)
        if uniformity is not None:
            rating = rate_metric(uniformity, UNIFORMITY_THRESHOLDS, lower_is_better=True)
            assessments.append(MetricAssessment(
                name="Binding Uniformity",
                value=uniformity,
                rating=rating,
                unit="Gini",
                context="Evenness of primer binding across genome (0=uniform, 1=clustered)",
                threshold_info="<0.3 excellent, <0.45 good, <0.6 acceptable"
            ))
            if rating in (QualityRating.POOR, QualityRating.CRITICAL):
                warnings.append("Uneven binding may cause amplification bias")

        # Assess dimer score
        dimer_score = self._calculate_dimer_score(primers)
        if dimer_score is not None:
            rating = rate_metric(dimer_score, DIMER_SCORE_THRESHOLDS, lower_is_better=True)
            assessments.append(MetricAssessment(
                name="Dimer Risk",
                value=dimer_score,
                rating=rating,
                unit="score",
                context="Risk of primer-dimer formation (lower is better)",
                threshold_info="<0.1 excellent, <0.2 good, <0.35 acceptable"
            ))
            if rating in (QualityRating.POOR, QualityRating.CRITICAL):
                warnings.append("High dimer risk may reduce amplification efficiency")

        # Calculate overall rating
        if not assessments:
            overall = QualityRating.POOR
        else:
            rating_scores = {
                QualityRating.EXCELLENT: 5,
                QualityRating.GOOD: 4,
                QualityRating.ACCEPTABLE: 3,
                QualityRating.POOR: 2,
                QualityRating.CRITICAL: 1
            }
            avg_score = sum(rating_scores[a.rating] for a in assessments) / len(assessments)
            if avg_score >= 4.5:
                overall = QualityRating.EXCELLENT
            elif avg_score >= 3.5:
                overall = QualityRating.GOOD
            elif avg_score >= 2.5:
                overall = QualityRating.ACCEPTABLE
            elif avg_score >= 1.5:
                overall = QualityRating.POOR
            else:
                overall = QualityRating.CRITICAL

        # Generate recommendation and next steps
        recommendation, next_steps = self._generate_recommendation(overall, warnings)

        return ResultsReport(
            primer_count=len(primers),
            assessments=assessments,
            overall_rating=overall,
            recommendation=recommendation,
            next_steps=next_steps,
            warnings=warnings
        )

    def _load_step4_results(self) -> List[Dict]:
        """Load results from step4 output file."""
        primers = []
        with open(self.step4_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                primers.append(row)
        return primers

    def _load_step3_results(self) -> List[Dict]:
        """Load results from step3 output file."""
        primers = []
        with open(self.step3_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                primers.append(row)
        return primers

    def _calculate_coverage(self, primers: List[Dict]) -> Optional[float]:
        """Calculate genome coverage from primer data."""
        # Look for coverage field
        for key in ['coverage', 'genome_coverage', 'fg_coverage']:
            values = [float(p.get(key, 0)) for p in primers if key in p]
            if values:
                return max(values)  # Return best coverage

        # Calculate from set_coverage if available
        for key in ['set_coverage']:
            values = [float(p.get(key, 0)) for p in primers if key in p]
            if values:
                return values[0]  # First set's coverage

        return None

    def _calculate_enrichment(self, primers: List[Dict]) -> Optional[float]:
        """Calculate enrichment ratio from primer data."""
        # Look for enrichment field
        for key in ['enrichment', 'enrichment_ratio', 'fg_bg_ratio']:
            values = [float(p.get(key, 0)) for p in primers if key in p]
            if values:
                return sum(values) / len(values)  # Average enrichment

        # Calculate from fg_freq and bg_freq
        try:
            fg_freqs = [float(p.get('fg_freq', 0)) for p in primers]
            bg_freqs = [float(p.get('bg_freq', 0)) for p in primers]

            # Filter out zero/near-zero background values to avoid division by zero
            valid_pairs = [(fg, bg) for fg, bg in zip(fg_freqs, bg_freqs) if bg > 1e-10]

            if valid_pairs:
                ratios = [fg / bg for fg, bg in valid_pairs]
                return sum(ratios) / len(ratios)
            elif any(fg > 0 for fg in fg_freqs):
                # Foreground signal but no background = very high enrichment
                return float('inf')
        except (ValueError, ZeroDivisionError):
            pass

        return None

    def _calculate_uniformity(self, primers: List[Dict]) -> Optional[float]:
        """Calculate binding uniformity (Gini index) from primer data."""
        for key in ['gini', 'gini_index', 'uniformity']:
            values = [float(p.get(key, 0)) for p in primers if key in p]
            if values:
                # Use max (worst-case) Gini for quality assessment
                # Higher Gini indicates more uneven binding, so worst-case matters
                return max(values)
        return None

    def _calculate_dimer_score(self, primers: List[Dict]) -> Optional[float]:
        """Calculate dimer risk score from primer data."""
        for key in ['dimer_score', 'dimer_risk', 'heterodimer_score']:
            values = [float(p.get(key, 0)) for p in primers if key in p]
            if values:
                return max(values)  # Worst case dimer risk
        return None

    def _generate_recommendation(self, overall: QualityRating,
                                  warnings: List[str]) -> Tuple[str, List[str]]:
        """Generate recommendation and next steps based on assessment."""
        if overall == QualityRating.EXCELLENT:
            recommendation = "Ready for synthesis. Excellent primer set quality."
            next_steps = [
                "Order primers for experimental validation",
                "Consider running simulation for additional confidence",
                "Prepare SWGA reaction with recommended conditions"
            ]
        elif overall == QualityRating.GOOD:
            recommendation = "Ready for synthesis. Good primer set quality."
            next_steps = [
                "Order primers for experimental validation",
                "Monitor coverage uniformity during amplification",
                "Consider backup primer designs if results are variable"
            ]
        elif overall == QualityRating.ACCEPTABLE:
            recommendation = "Consider optimization before synthesis."
            next_steps = [
                "Re-run optimization with different parameters",
                "Try background-aware optimization for better specificity",
                "Consider increasing primer pool size",
                "If time-critical, proceed with caution"
            ]
        elif overall == QualityRating.POOR:
            recommendation = "Optimization recommended before synthesis."
            next_steps = [
                "Review parameter settings (GC bounds, Tm range)",
                "Try different optimization method (dominating-set, background-aware)",
                "Increase max_primer to expand candidate pool",
                "Check if genome is suitable for SWGA (run analyze-genome)"
            ]
        else:  # CRITICAL
            recommendation = "Do not proceed. Significant issues detected."
            next_steps = [
                "Review genome suitability with analyze-genome",
                "Check parameter configuration with validate-params",
                "Consider if target genome is appropriate for SWGA",
                "Contact support if issues persist"
            ]

        return recommendation, next_steps

    def print_report(self, report: ResultsReport) -> None:
        """Print formatted report to console."""
        print("\n" + "=" * 60)
        print("PRIMER SET QUALITY ASSESSMENT")
        print("=" * 60)

        print(f"\nPrimer count: {report.primer_count}")
        print(f"Overall rating: {report.overall_rating.value}")

        if report.assessments:
            print("\n--- Metrics ---")
            for assessment in report.assessments:
                print(f"\n{assessment.name}: {assessment.value:.2f}{assessment.unit}")
                print(f"  Rating: {assessment.rating.value}")
                print(f"  Context: {assessment.context}")
                print(f"  Thresholds: {assessment.threshold_info}")

        if report.warnings:
            print("\n--- Warnings ---")
            for warning in report.warnings:
                print(f"  - {warning}")

        print(f"\n--- Recommendation ---")
        print(f"  {report.recommendation}")

        print(f"\n--- Next Steps ---")
        for i, step in enumerate(report.next_steps, 1):
            print(f"  {i}. {step}")

        print()


def interpret_results(results_dir: str, verbose: bool = True) -> ResultsReport:
    """
    Interpret pipeline results.

    Args:
        results_dir: Path to results directory
        verbose: Print report to console

    Returns:
        ResultsReport with assessments and recommendations
    """
    interpreter = ResultsInterpreter(results_dir)
    report = interpreter.analyze()

    if verbose:
        interpreter.print_report(report)

    return report


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python results_interpreter.py <results_dir>")
        sys.exit(1)

    try:
        interpret_results(sys.argv[1])
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
