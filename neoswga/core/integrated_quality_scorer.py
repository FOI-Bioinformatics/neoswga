"""
Integrated Quality Scoring System for SWGA Primers.

CRITICAL INSIGHT:
Individual quality checks (strand bias, dimers, 3' stability, etc.) are
valuable, but primer sets need HOLISTIC evaluation combining ALL quality
dimensions into a single, actionable score.

This module provides the integration layer that:
1. Combines all quality metrics into unified scores
2. Ranks primers by overall quality
3. Identifies specific weaknesses for targeted improvement
4. Provides actionable recommendations

This is the "glue" that makes Sprint 1 & 2 features immediately usable
in production pipelines.

Expected Impact:
- Simplifies primer selection (single quality score vs. 10+ metrics)
- Enables intelligent ranking and filtering
- Identifies root causes of poor primer quality
- Guides iterative improvement

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.3 - Integration Layer
"""

import logging
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass, field
from enum import Enum
import numpy as np

from neoswga.core.strand_bias_analyzer import (
    StrandBiasAnalyzer,
    create_strand_bias_analyzer,
    StrandBindingSite
)
from neoswga.core.dimer_network_analyzer import (
    DimerNetworkAnalyzer,
    create_dimer_network_analyzer
)
from neoswga.core.three_prime_stability import (
    ThreePrimeStabilityAnalyzer,
    create_three_prime_analyzer,
    create_three_prime_analyzer_adaptive
)
from neoswga.core.reaction_conditions import ReactionConditions

logger = logging.getLogger(__name__)


class QualityDimension(Enum):
    """Quality dimensions for primer evaluation."""
    STRAND_BIAS = "strand_bias"
    DIMER_NETWORK = "dimer_network"
    THREE_PRIME_STABILITY = "three_prime_stability"
    SEQUENCE_COMPLEXITY = "sequence_complexity"
    THERMODYNAMICS = "thermodynamics"


@dataclass
class PrimerQualityScore:
    """Comprehensive quality score for a single primer."""
    primer: str

    # Individual dimension scores (0-1, higher is better)
    strand_bias_score: float
    dimer_score: float  # Inverse of dimer severity
    three_prime_score: float
    complexity_score: float
    thermo_score: float

    # Overall composite score (0-1)
    overall_score: float

    # Pass/fail flags
    passes_strand_bias: bool
    passes_dimer: bool
    passes_three_prime: bool
    passes_all: bool

    # Failure reasons (if any)
    failure_reasons: List[str] = field(default_factory=list)

    # Rank among all primers (1 = best)
    rank: Optional[int] = None

    def __str__(self):
        status = "PASS" if self.passes_all else "FAIL"
        return (f"{self.primer} ({status}): Overall={self.overall_score:.3f}, Rank={self.rank}\n"
                f"  Strand: {self.strand_bias_score:.2f}, "
                f"Dimer: {self.dimer_score:.2f}, "
                f"3': {self.three_prime_score:.2f}\n"
                f"  Failures: {', '.join(self.failure_reasons) if self.failure_reasons else 'None'}")


@dataclass
class SetQualityScore:
    """Comprehensive quality score for a primer set."""
    primers: List[str]

    # Set-level scores (0-1)
    mean_overall_score: float
    min_overall_score: float  # Weakest link
    set_dimer_score: float  # Network-level dimer quality
    set_strand_bias_score: float

    # Pass/fail
    passes: bool
    failure_reason: Optional[str] = None

    # Number of primers failing each dimension
    num_failing_strand: int = 0
    num_failing_dimer: int = 0
    num_failing_three_prime: int = 0

    def __str__(self):
        status = "PASS" if self.passes else f"FAIL ({self.failure_reason})"
        return (f"Primer Set ({len(self.primers)} primers): {status}\n"
                f"  Mean quality: {self.mean_overall_score:.3f}, "
                f"Min quality: {self.min_overall_score:.3f}\n"
                f"  Failures: Strand={self.num_failing_strand}, "
                f"Dimer={self.num_failing_dimer}, "
                f"3'={self.num_failing_three_prime}")


class IntegratedQualityScorer:
    """
    Integrated quality scoring combining all quality dimensions.

    Combines:
    1. Strand bias (from strand_bias_analyzer)
    2. Dimer network (from dimer_network_analyzer)
    3. 3' stability (from three_prime_stability)
    4. Sequence complexity
    5. Thermodynamics

    Algorithm:
    1. Score each primer on all dimensions (0-1 scale)
    2. Calculate weighted composite score
    3. Rank primers by overall quality
    4. Evaluate set-level quality
    5. Provide actionable recommendations

    Scoring weights (default):
    - Dimer network: 35% (largest impact on success)
    - 3' stability: 25% (critical for extension)
    - Strand bias: 20% (affects coverage uniformity)
    - Thermodynamics: 15% (binding strength)
    - Complexity: 5% (background binding proxy)
    """

    def __init__(self,
                 conditions: Optional[ReactionConditions] = None,
                 stringency: str = 'moderate',
                 weights: Optional[Dict[str, float]] = None,
                 genome_gc: Optional[float] = None):
        """
        Initialize integrated quality scorer.

        Args:
            conditions: Reaction conditions
            stringency: 'lenient', 'moderate', or 'strict'
            weights: Custom weights for dimensions (default: dimer=0.35, 3'=0.25,
                    strand=0.20, thermo=0.15, complexity=0.05)
            genome_gc: Target genome GC content (0-1, optional). If provided,
                      uses genome-adaptive QA thresholds for extreme GC genomes.

        Examples:
            # Standard QA (balanced genome)
            >>> scorer = IntegratedQualityScorer(stringency='moderate')

            # Genome-adaptive QA for AT-rich genome (Francisella)
            >>> scorer = IntegratedQualityScorer(
            ...     stringency='moderate',
            ...     genome_gc=0.32
            ... )
        """
        self.conditions = conditions
        self.stringency = stringency
        self.genome_gc = genome_gc

        # Initialize component analyzers
        self.strand_analyzer = create_strand_bias_analyzer(stringency)
        self.dimer_analyzer = create_dimer_network_analyzer(stringency, conditions)

        # Use genome-adaptive analyzer if genome_gc provided
        if genome_gc is not None:
            self.three_prime_analyzer = create_three_prime_analyzer_adaptive(
                genome_gc=genome_gc,
                stringency=stringency,
                conditions=conditions
            )
            logger.info(f"Using genome-adaptive QA for {genome_gc:.1%} GC genome")
        else:
            self.three_prime_analyzer = create_three_prime_analyzer(stringency, conditions)

        # Primer selection weights emphasize dimer avoidance and 3' stability,
        # which are the primary causes of SWGA failure at the primer level.
        # This is distinct from report/quality.py which weights pipeline results
        # by coverage (35%) and specificity (30%) for grading experimental outcomes.
        if weights is None:
            self.weights = {
                'dimer': 0.35,
                'three_prime': 0.25,
                'strand_bias': 0.20,
                'thermodynamics': 0.15,
                'complexity': 0.05
            }
        else:
            self.weights = weights

        # Validate weights sum to 1
        total = sum(self.weights.values())
        if not np.isclose(total, 1.0):
            logger.warning(f"Weights sum to {total}, normalizing to 1.0")
            self.weights = {k: v/total for k, v in self.weights.items()}

    def score_primer(self,
                    primer: str,
                    binding_sites: Optional[List[StrandBindingSite]] = None) -> PrimerQualityScore:
        """
        Calculate comprehensive quality score for a single primer.

        Args:
            primer: Primer sequence
            binding_sites: Strand binding sites (optional, for strand bias)

        Returns:
            PrimerQualityScore with all metrics
        """
        failure_reasons = []

        # 1. Strand bias score
        if binding_sites:
            strand_bias = self.strand_analyzer.analyze_primer(primer, binding_sites)
            # Convert bias score (0=good, 1=bad) to quality score (0=bad, 1=good)
            strand_score = 1.0 - strand_bias.bias_score
            passes_strand = strand_bias.passes
            if not passes_strand:
                failure_reasons.append(f"Strand: {strand_bias.failure_reason}")
        else:
            # No binding sites provided, assume balanced
            strand_score = 1.0
            passes_strand = True

        # 2. 3' stability score
        three_prime = self.three_prime_analyzer.analyze_primer(primer)
        three_prime_score = three_prime.stability_score
        passes_three_prime = three_prime.passes
        if not passes_three_prime:
            failure_reasons.append(f"3': {three_prime.failure_reason}")

        # 3. Sequence complexity score
        complexity_score = self._calculate_complexity_score(primer)

        # 4. Thermodynamic score
        thermo_score = self._calculate_thermo_score(primer)

        # 5. Dimer score (individual, will be updated by set analysis)
        # For now, assume no dimers (will be set by analyze_primer_set)
        dimer_score = 1.0
        passes_dimer = True

        # Calculate overall composite score
        overall_score = (
            self.weights['strand_bias'] * strand_score +
            self.weights['three_prime'] * three_prime_score +
            self.weights['complexity'] * complexity_score +
            self.weights['thermodynamics'] * thermo_score +
            self.weights['dimer'] * dimer_score
        )

        return PrimerQualityScore(
            primer=primer,
            strand_bias_score=strand_score,
            dimer_score=dimer_score,
            three_prime_score=three_prime_score,
            complexity_score=complexity_score,
            thermo_score=thermo_score,
            overall_score=overall_score,
            passes_strand_bias=passes_strand,
            passes_dimer=passes_dimer,
            passes_three_prime=passes_three_prime,
            passes_all=(passes_strand and passes_dimer and passes_three_prime),
            failure_reasons=failure_reasons
        )

    def analyze_primer_set(self,
                          primers: List[str],
                          binding_sites_dict: Optional[Dict[str, List[StrandBindingSite]]] = None,
                          verbose: bool = False) -> Tuple[List[PrimerQualityScore], SetQualityScore]:
        """
        Comprehensive quality analysis of primer set.

        Args:
            primers: List of primer sequences
            binding_sites_dict: Dict mapping primer -> binding sites (optional)
            verbose: Print detailed analysis

        Returns:
            (primer_scores, set_score)
        """
        logger.info(f"Integrated quality scoring for {len(primers)} primers...")

        # Step 1: Score individual primers (without dimer scores)
        primer_scores = {}
        for primer in primers:
            sites = binding_sites_dict.get(primer, None) if binding_sites_dict else None
            score = self.score_primer(primer, sites)
            primer_scores[primer] = score

        # Step 2: Network dimer analysis
        dimer_metrics, dimer_profiles, _ = self.dimer_analyzer.analyze_primer_set(
            primers, verbose=False
        )

        # Update dimer scores based on network analysis
        for primer in primers:
            profile = dimer_profiles[primer]

            # Dimer score: inverse of mean severity (0=bad, 1=good)
            # max_severity also considered
            dimer_score = 1.0 - min(1.0, (profile.mean_severity + profile.max_severity) / 2)

            primer_scores[primer].dimer_score = dimer_score
            primer_scores[primer].passes_dimer = not profile.is_hub

            # Recalculate overall score with updated dimer score
            score = primer_scores[primer]
            score.overall_score = (
                self.weights['strand_bias'] * score.strand_bias_score +
                self.weights['three_prime'] * score.three_prime_score +
                self.weights['complexity'] * score.complexity_score +
                self.weights['thermodynamics'] * score.thermo_score +
                self.weights['dimer'] * score.dimer_score
            )

            # Update pass/fail
            score.passes_all = (score.passes_strand_bias and
                              score.passes_dimer and
                              score.passes_three_prime)

            if not score.passes_dimer:
                score.failure_reasons.append(f"Dimer: Hub primer ({profile.num_interactions} interactions)")

        # Step 3: Rank primers by overall quality
        sorted_scores = sorted(primer_scores.values(),
                             key=lambda x: x.overall_score,
                             reverse=True)

        for rank, score in enumerate(sorted_scores, 1):
            score.rank = rank

        # Step 4: Calculate set-level metrics
        overall_scores = [s.overall_score for s in primer_scores.values()]
        mean_overall = np.mean(overall_scores)
        min_overall = min(overall_scores)

        # Set-level dimer score
        set_dimer_score = 1.0 - dimer_metrics.mean_severity

        # Set-level strand bias (if available)
        if binding_sites_dict:
            strand_biases = [self.strand_analyzer.analyze_primer(p, binding_sites_dict.get(p, []))
                           for p in primers]
            set_strand_metrics = self.strand_analyzer.analyze_primer_set(strand_biases)
            set_strand_score = 1.0 - set_strand_metrics.mean_bias_score
        else:
            set_strand_score = 1.0

        # Count failures
        num_failing_strand = sum(1 for s in primer_scores.values() if not s.passes_strand_bias)
        num_failing_dimer = sum(1 for s in primer_scores.values() if not s.passes_dimer)
        num_failing_three_prime = sum(1 for s in primer_scores.values() if not s.passes_three_prime)

        # Set-level pass/fail
        passes = dimer_metrics.passes and min_overall >= 0.5

        failure_reason = None
        if not passes:
            if not dimer_metrics.passes:
                failure_reason = f"Dimer network: {dimer_metrics.failure_reason}"
            elif min_overall < 0.5:
                failure_reason = f"Min quality {min_overall:.2f} < 0.5 (weakest link)"

        set_score = SetQualityScore(
            primers=primers,
            mean_overall_score=mean_overall,
            min_overall_score=min_overall,
            set_dimer_score=set_dimer_score,
            set_strand_bias_score=set_strand_score,
            passes=passes,
            failure_reason=failure_reason,
            num_failing_strand=num_failing_strand,
            num_failing_dimer=num_failing_dimer,
            num_failing_three_prime=num_failing_three_prime
        )

        if verbose:
            logger.info(f"\n{set_score}")
            logger.info("\nTop 5 primers by quality:")
            for score in sorted_scores[:5]:
                logger.info(f"  {score.rank}. {score.primer[:10]}... "
                          f"(score={score.overall_score:.3f})")

            if num_failing_strand + num_failing_dimer + num_failing_three_prime > 0:
                logger.info(f"\nQuality issues detected:")
                logger.info(f"  {num_failing_strand} primers with strand bias")
                logger.info(f"  {num_failing_dimer} primers with dimer issues")
                logger.info(f"  {num_failing_three_prime} primers with poor 3' stability")

        return list(primer_scores.values()), set_score

    def _calculate_complexity_score(self, primer: str) -> float:
        """
        Calculate sequence complexity score (k-mer diversity).

        Low complexity sequences (e.g., AAAAAAA) have more non-specific binding.

        Args:
            primer: Primer sequence

        Returns:
            Score 0-1 (1 = high complexity, diverse k-mers)
        """
        if len(primer) < 4:
            return 1.0

        # Count unique 4-mers
        kmers = set()
        for i in range(len(primer) - 3):
            kmer = primer[i:i+4]
            kmers.add(kmer)

        # Maximum possible unique 4-mers for this length
        max_kmers = len(primer) - 3

        # Complexity = observed / maximum
        complexity = len(kmers) / max_kmers

        return complexity

    def _calculate_thermo_score(self, primer: str) -> float:
        """
        Calculate thermodynamic score.

        Optimal primers have Tm in target range for reaction temperature.

        Args:
            primer: Primer sequence

        Returns:
            Score 0-1 (1 = optimal Tm range)
        """
        # Simple GC-based Tm estimate (Wallace's rule for short primers)
        gc_count = sum(1 for b in primer if b in 'GC')
        at_count = len(primer) - gc_count

        tm_estimate = 2 * at_count + 4 * gc_count

        # Optimal Tm range: 50-65°C for SWGA
        # But reaction temp is 30-42°C, so primers should be ~10-20°C above
        target_temp = self.conditions.temp if self.conditions else 30.0
        optimal_tm_low = target_temp + 15
        optimal_tm_high = target_temp + 30

        if optimal_tm_low <= tm_estimate <= optimal_tm_high:
            score = 1.0
        elif tm_estimate < optimal_tm_low:
            # Too low
            score = max(0.0, tm_estimate / optimal_tm_low)
        else:
            # Too high
            score = max(0.0, 1.0 - (tm_estimate - optimal_tm_high) / 20.0)

        return score

    def get_recommendations(self,
                          primer_scores: List[PrimerQualityScore],
                          set_score: SetQualityScore) -> List[str]:
        """
        Generate actionable recommendations for improving primer set.

        Args:
            primer_scores: Individual primer scores
            set_score: Set-level score

        Returns:
            List of recommendations
        """
        recommendations = []

        # Check if set passes
        if not set_score.passes:
            recommendations.append(f"SET FAILS: {set_score.failure_reason}")

        # Identify weakest dimension
        mean_scores = {
            'strand': np.mean([s.strand_bias_score for s in primer_scores]),
            'dimer': np.mean([s.dimer_score for s in primer_scores]),
            'three_prime': np.mean([s.three_prime_score for s in primer_scores]),
            'complexity': np.mean([s.complexity_score for s in primer_scores]),
            'thermo': np.mean([s.thermo_score for s in primer_scores])
        }

        weakest_dim = min(mean_scores, key=mean_scores.get)
        weakest_score = mean_scores[weakest_dim]

        if weakest_score < 0.7:
            recommendations.append(
                f"Weakest dimension: {weakest_dim} (mean={weakest_score:.2f})"
            )

        # Specific recommendations
        if set_score.num_failing_dimer > 0:
            recommendations.append(
                f"Replace {set_score.num_failing_dimer} hub primers with low-dimer alternatives"
            )

        if set_score.num_failing_three_prime > 0:
            recommendations.append(
                f"Replace {set_score.num_failing_three_prime} primers with poor 3' stability"
            )

        if set_score.num_failing_strand > 0:
            recommendations.append(
                f"Address strand bias in {set_score.num_failing_strand} primers"
            )

        # Weakest link
        worst_primer = min(primer_scores, key=lambda x: x.overall_score)
        if worst_primer.overall_score < 0.5:
            recommendations.append(
                f"Weakest link: {worst_primer.primer[:10]}... (score={worst_primer.overall_score:.2f})"
            )

        return recommendations


def create_quality_scorer(stringency: str = 'moderate',
                         conditions: Optional[ReactionConditions] = None,
                         weights: Optional[Dict[str, float]] = None,
                         genome_gc: Optional[float] = None) -> IntegratedQualityScorer:
    """
    Factory function to create IntegratedQualityScorer with preset stringency.

    Args:
        stringency: 'lenient', 'moderate' (default), or 'strict'
        conditions: Reaction conditions (optional)
        weights: Custom scoring weights (optional)
        genome_gc: Target genome GC content (0-1, optional). Enables genome-adaptive QA.

    Returns:
        Configured IntegratedQualityScorer instance

    Examples:
        # Standard QA
        >>> scorer = create_quality_scorer('strict')
        >>> quality = scorer.score_primer('ACGTACGTACGTGC')
        >>> print(f"Quality score: {quality.overall_score:.3f}")

        # Genome-adaptive QA for Francisella
        >>> scorer = create_quality_scorer('moderate', genome_gc=0.32)
        >>> quality = scorer.score_primer('ATATATATATAT')  # AT-rich primer
    """
    return IntegratedQualityScorer(
        conditions=conditions,
        stringency=stringency,
        weights=weights,
        genome_gc=genome_gc
    )


def quick_score_primers(primers: List[str],
                       binding_sites_dict: Optional[Dict[str, List[StrandBindingSite]]] = None,
                       stringency: str = 'moderate',
                       conditions: Optional[ReactionConditions] = None,
                       genome_gc: Optional[float] = None) -> Tuple[List[PrimerQualityScore], SetQualityScore]:
    """
    Quick utility for integrated quality scoring.

    Args:
        primers: Primer sequences
        binding_sites_dict: Optional binding sites
        stringency: Quality stringency
        conditions: Reaction conditions
        genome_gc: Target genome GC content (0-1, optional). Enables genome-adaptive QA.

    Returns:
        (primer_scores, set_score)
    """
    scorer = IntegratedQualityScorer(conditions, stringency, genome_gc=genome_gc)
    return scorer.analyze_primer_set(primers, binding_sites_dict, verbose=False)


if __name__ == "__main__":
    # Example usage
    print("Testing Integrated Quality Scorer...\n")

    # Create scorer
    scorer = IntegratedQualityScorer(stringency='moderate')

    # Example primer set
    test_primers = [
        "ACGTACGTACGTGC",  # Good primer
        "AAAAAAAAAAAAA",   # Low complexity
        "GCGCGCGCGCGCG",   # High GC, potential dimers
        "ACGTACGTACGTAT",  # Marginal 3' stability
        "TGCATGCATGCATG",  # Balanced
    ]

    print(f"Analyzing {len(test_primers)} primers...\n")

    # Analyze
    primer_scores, set_score = scorer.analyze_primer_set(test_primers, verbose=True)

    print("\n" + "="*60)
    print("RECOMMENDATIONS:")
    print("="*60)
    recommendations = scorer.get_recommendations(primer_scores, set_score)
    for i, rec in enumerate(recommendations, 1):
        print(f"{i}. {rec}")

    print("\n" + "="*60)
    print("Integrated Quality Scorer ready for pipeline integration!")
