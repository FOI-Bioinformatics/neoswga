"""
Pipeline integration for quality assurance modules.

This module provides integration points for QA modules into the SWGA pipeline:
- Post-Step2: Individual primer quality filtering
- Step3: Combined RF + QA scoring
- Step4: Network-aware optimization

Integration Points:
    1. apply_post_step2_qa_filter() - Filter primers after Step 2
    2. combine_rf_qa_scores() - Merge RF and QA scores in Step 3
    3. create_qa_aware_optimizer() - Network-aware optimization for Step 4
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

from neoswga.core.dimer_network_analyzer import (
    DimerNetworkAnalyzer,
    create_dimer_network_analyzer,
    filter_primer_set_by_dimer_network,
)
from neoswga.core.integrated_quality_scorer import IntegratedQualityScorer, create_quality_scorer
from neoswga.core.io_utils import primer_column
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.strand_bias_analyzer import StrandBiasAnalyzer, create_strand_bias_analyzer
from neoswga.core.three_prime_stability import (
    ThreePrimeStabilityAnalyzer,
    create_three_prime_analyzer,
)


@dataclass
class QAFilterConfig:
    """Configuration for QA filtering."""

    enable_three_prime: bool = True
    enable_strand_bias: bool = True
    enable_dimer_network: bool = True
    enable_integrated_scorer: bool = True

    three_prime_stringency: str = "moderate"
    strand_bias_stringency: str = "moderate"
    dimer_stringency: str = "moderate"

    qa_weight: float = 0.3  # Weight for QA score in combined scoring
    rf_weight: float = 0.7  # Weight for RF score in combined scoring

    max_hub_degree: int = 5  # Maximum dimer interactions for network optimization


@dataclass
class QAFilterResult:
    """Results from QA filtering."""

    primers_in: int
    primers_out: int
    primers_filtered: int
    filter_reasons: Dict[str, int]
    mean_qa_score: float
    filtered_df: pd.DataFrame


def apply_post_step2_qa_filter(
    step2_df: pd.DataFrame,
    config: Optional[QAFilterConfig] = None,
    conditions: Optional[ReactionConditions] = None,
    verbose: bool = True,
) -> QAFilterResult:
    """
    Apply QA filtering to primers after Step 2.

    This is the primary integration point - filters primers by quality
    metrics before they enter Step 3 (RF scoring).

    Args:
        step2_df: DataFrame from Step 2 with 'seq' column
        config: QA filter configuration (default: moderate stringency)
        conditions: Reaction conditions for thermodynamic calculations
        verbose: Print filtering progress

    Returns:
        QAFilterResult with filtered DataFrame and statistics

    Example:
        >>> df = pd.read_csv('step2_df.csv')
        >>> result = apply_post_step2_qa_filter(df)
        >>> result.filtered_df.to_csv('step2_qa_filtered.csv')
        >>> print(f"Filtered {result.primers_filtered} primers")
    """
    if config is None:
        config = QAFilterConfig()

    if conditions is None:
        conditions = ReactionConditions()

    primers_in = len(step2_df)
    filter_reasons = {}

    # Work on a copy
    df = step2_df.copy()
    # Resolve the primer-sequence column once ('primer' canonical, 'seq' legacy).
    seq_col = primer_column(df)
    df["qa_pass"] = True
    df["qa_score"] = 1.0
    df["filter_reason"] = ""

    if verbose:
        logger.info(f"\n{'='*80}")
        logger.info("POST-STEP2 QA FILTERING")
        logger.info(f"{'='*80}")
        logger.info(f"Input primers: {primers_in}")

    # 1. 3' Stability Filter
    if config.enable_three_prime:
        analyzer = create_three_prime_analyzer(config.three_prime_stringency, conditions)

        if verbose:
            logger.info(f"\n[1/3] 3' Stability Analysis ({config.three_prime_stringency})...")

        # Vectorized approach: analyze only primers that pass, in batch
        passing_mask = df["qa_pass"]
        if passing_mask.any():
            # Analyze all passing primers at once
            passing_primers = df.loc[passing_mask, seq_col].tolist()

            # Batch analyze (still sequential, but avoids iterrows overhead)
            stability_results = {
                primer: analyzer.analyze_primer(primer) for primer in passing_primers
            }

            # Apply results using vectorized operations
            def apply_stability(row):
                if not row["qa_pass"]:
                    return row["qa_pass"], row["qa_score"], row["filter_reason"]
                stability = stability_results.get(row[seq_col])
                if stability and not stability.passes:
                    return False, row["qa_score"], "3prime_stability"
                elif stability:
                    return True, row["qa_score"] * stability.stability_score, row["filter_reason"]
                return row["qa_pass"], row["qa_score"], row["filter_reason"]

            results = df.apply(apply_stability, axis=1, result_type="expand")
            df["qa_pass"] = results[0]
            df["qa_score"] = results[1]
            df["filter_reason"] = results[2]

            failed_3p = len(passing_primers) - df["qa_pass"].sum()
        else:
            failed_3p = 0

        filter_reasons["3prime_stability"] = failed_3p
        if verbose:
            logger.info(f"   Filtered: {failed_3p} primers with poor 3' stability")
            logger.info(f"   Remaining: {df['qa_pass'].sum()}")

    # 2. Dimer Network Filter
    if config.enable_dimer_network:
        if verbose:
            logger.info(f"\n[2/3] Dimer Network Analysis ({config.dimer_stringency})...")

        # Only analyze primers that passed previous filters
        passing_primers = df[df["qa_pass"]][seq_col].tolist()

        if len(passing_primers) > 0:
            analyzer = create_dimer_network_analyzer(config.dimer_stringency, conditions)

            # Analyze network
            metrics, profiles, matrix = analyzer.analyze_primer_set(passing_primers)

            # Filter out hub primers (many dimer partners).
            # analyze_primer_set returns profiles as Dict[str, PrimerDimerProfile],
            # so iterate .values(); iterating the dict itself yields primer strings.
            # The per-primer partner count is `num_interactions` (there is no
            # `degree` attribute).
            failed_dimer = 0
            for profile in profiles.values():
                if profile.num_interactions > config.max_hub_degree:
                    # Find this primer in df and mark as failed
                    mask = (df[seq_col] == profile.primer) & df["qa_pass"]
                    df.loc[mask, "qa_pass"] = False
                    df.loc[mask, "filter_reason"] = "dimer_hub"
                    failed_dimer += 1
                else:
                    # Update score based on dimer tendency
                    mask = (df[seq_col] == profile.primer) & df["qa_pass"]
                    if mask.any():
                        # Lower score for primers with more interactions
                        dimer_score = max(0.0, 1.0 - profile.num_interactions / 10.0)
                        df.loc[mask, "qa_score"] *= dimer_score

            filter_reasons["dimer_hub"] = failed_dimer
            if verbose:
                logger.info(f"   Network size: {metrics.num_primers} primers")
                logger.info(f"   Hub primers filtered: {failed_dimer}")
                logger.info(f"   Remaining: {df['qa_pass'].sum()}")
        else:
            if verbose:
                logger.info("   No primers to analyze")

    # 3. Integrated Quality Score
    if config.enable_integrated_scorer:
        if verbose:
            logger.info(f"\n[3/3] Integrated Quality Scoring...")

        scorer = create_quality_scorer("moderate", conditions)

        # Batch scoring: pre-compute all scores for passing primers
        passing_mask = df["qa_pass"]
        if passing_mask.any():
            passing_primers = df.loc[passing_mask, seq_col].tolist()

            # Batch score all primers
            quality_scores = {}
            for primer in passing_primers:
                try:
                    quality = scorer.score_primer(primer)
                    quality_scores[primer] = quality.overall_score
                except (KeyError, ValueError, AttributeError, TypeError):
                    quality_scores[primer] = 1.0  # Default score on error

            # Vectorized update using map (much faster than iterrows)
            df.loc[passing_mask, "qa_score"] *= df.loc[passing_mask, seq_col].map(quality_scores)

        if verbose:
            mean_score = df[df["qa_pass"]]["qa_score"].mean()
            logger.info(f"   Mean QA score: {mean_score:.3f}")

    # Final results
    filtered_df = df[df["qa_pass"]].copy()
    primers_out = len(filtered_df)
    primers_filtered = primers_in - primers_out
    mean_qa_score = filtered_df["qa_score"].mean() if primers_out > 0 else 0.0

    if verbose:
        logger.info(f"\n{'='*80}")
        logger.info("QA FILTERING COMPLETE")
        logger.info(f"{'='*80}")
        logger.info(f"Input primers: {primers_in}")
        logger.info(f"Output primers: {primers_out}")
        logger.info(f"Filtered: {primers_filtered} ({100*primers_filtered/primers_in:.1f}%)")
        logger.info(f"Mean QA score: {mean_qa_score:.3f}")
        logger.info(f"\nFilter breakdown:")
        for reason, count in filter_reasons.items():
            logger.info(f"  {reason}: {count}")

    return QAFilterResult(
        primers_in=primers_in,
        primers_out=primers_out,
        primers_filtered=primers_filtered,
        filter_reasons=filter_reasons,
        mean_qa_score=mean_qa_score,
        filtered_df=filtered_df,
    )


def combine_rf_qa_scores(
    step3_df: pd.DataFrame,
    qa_scores: Dict[str, float],
    rf_weight: float = 0.7,
    qa_weight: float = 0.3,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Combine Random Forest and QA scores for Step 3.

    Creates a composite score that balances ML predictions with
    quality metrics.

    Args:
        step3_df: DataFrame from Step 3 with 'seq' and 'score' columns
        qa_scores: Dict mapping primer sequences to QA scores (0-1)
        rf_weight: Weight for RF score (default: 0.7)
        qa_weight: Weight for QA score (default: 0.3)
        verbose: Print scoring statistics

    Returns:
        DataFrame with added 'qa_score' and 'composite_score' columns

    Example:
        >>> df = pd.read_csv('step3_df.csv')
        >>> qa_scores = {'ACGTACGT': 0.85, 'TGCATGCA': 0.92}
        >>> df = combine_rf_qa_scores(df, qa_scores)
        >>> df.sort_values('composite_score', ascending=False)
    """
    if verbose:
        logger.info(f"\n{'='*80}")
        logger.info("COMBINING RF AND QA SCORES")
        logger.info(f"{'='*80}")
        logger.info(f"RF weight: {rf_weight:.2f}")
        logger.info(f"QA weight: {qa_weight:.2f}")

    df = step3_df.copy()
    seq_col = primer_column(df)

    # Add QA scores
    df["qa_score"] = df[seq_col].map(qa_scores).fillna(0.5)

    # Normalize RF scores to 0-1 range
    rf_min = df["score"].min()
    rf_max = df["score"].max()
    if rf_max > rf_min:
        df["rf_score_norm"] = (df["score"] - rf_min) / (rf_max - rf_min)
    else:
        df["rf_score_norm"] = 1.0

    # Compute composite score
    df["composite_score"] = rf_weight * df["rf_score_norm"] + qa_weight * df["qa_score"]

    if verbose:
        logger.info(f"\nScore statistics:")
        logger.info(f"  RF scores: {df['score'].min():.2f} - {df['score'].max():.2f}")
        logger.info(f"  QA scores: {df['qa_score'].min():.3f} - {df['qa_score'].max():.3f}")
        logger.info(
            f"  Composite: {df['composite_score'].min():.3f} - {df['composite_score'].max():.3f}"
        )

        # Show correlation
        correlation = df["rf_score_norm"].corr(df["qa_score"])
        logger.info(f"\nRF-QA correlation: {correlation:.3f}")

    return df


def create_qa_aware_optimizer(
    config: Optional[QAFilterConfig] = None, conditions: Optional[ReactionConditions] = None
):
    """
    Create an optimizer wrapper that uses QA-aware primer selection.

    This wraps the standard Step 4 optimization to prefer primers
    with better quality metrics during greedy search.

    Args:
        config: QA filter configuration
        conditions: Reaction conditions

    Returns:
        QAAwareOptimizer instance

    Example:
        >>> optimizer = create_qa_aware_optimizer()
        >>> best_set = optimizer.optimize(candidate_primers, target_set_size=5)
    """
    if config is None:
        config = QAFilterConfig()

    if conditions is None:
        conditions = ReactionConditions()

    return QAAwareOptimizer(config, conditions)


class QAAwareOptimizer:
    """
    Optimizer that incorporates quality metrics into primer set selection.

    This extends the standard BFS optimization by:
    - Filtering out high-hub primers
    - Preferring primers with better QA scores
    - Maintaining network diversity
    """

    def __init__(self, config: QAFilterConfig, conditions: ReactionConditions):
        self.config = config
        self.conditions = conditions
        self.dimer_analyzer = create_dimer_network_analyzer(config.dimer_stringency, conditions)

    def prefilter_candidates(
        self, candidates: List[str], verbose: bool = True
    ) -> Tuple[List[str], Dict]:
        """
        Pre-filter candidates before optimization.

        Args:
            candidates: List of candidate primer sequences
            verbose: Print filtering statistics

        Returns:
            (filtered_primers, metadata) tuple
        """
        # FIXED: Handle empty candidate list
        if not candidates:
            if verbose:
                logger.info("\nWarning: Empty candidate list provided to prefilter")
            return [], {"original_count": 0, "filtered_count": 0, "hub_count": 0}

        if verbose:
            logger.info(f"\n{'='*80}")
            logger.info("QA-AWARE PRE-FILTERING")
            logger.info(f"{'='*80}")
            logger.info(f"Input candidates: {len(candidates)}")

        # Analyze dimer network
        metrics, profiles, matrix = self.dimer_analyzer.analyze_primer_set(candidates)

        # Filter out hub primers
        hub_threshold = self.config.max_hub_degree
        filtered = []
        hub_primers = []

        # profiles is Dict[str, PrimerDimerProfile]; the per-primer partner count
        # is `num_interactions`.
        for profile in profiles.values():
            if profile.num_interactions <= hub_threshold:
                filtered.append(profile.primer)
            else:
                hub_primers.append((profile.primer, profile.num_interactions))

        if verbose:
            # DimerNetworkMetrics has no mean_degree; derive it from the profiles.
            mean_degree = (
                sum(p.num_interactions for p in profiles.values()) / len(profiles)
                if profiles
                else 0.0
            )
            logger.info(f"\nNetwork analysis:")
            logger.info(f"  Total interactions: {metrics.total_interactions}")
            logger.info(f"  Mean degree: {mean_degree:.2f}")
            logger.info(f"  Hub primers (>{hub_threshold} interactions): {len(hub_primers)}")
            logger.info(f"  Remaining candidates: {len(filtered)}")

            if hub_primers:
                logger.info(f"\nTop hub primers removed:")
                for primer, degree in sorted(hub_primers, key=lambda x: x[1], reverse=True)[:5]:
                    logger.info(f"    {primer}: {degree} interactions")

        metadata = {
            "original_count": len(candidates),
            "filtered_count": len(filtered),
            "hub_count": len(hub_primers),
            "network_metrics": metrics,
        }

        return filtered, metadata

    def rank_by_quality(
        self, candidates: List[str], verbose: bool = True
    ) -> List[Tuple[str, float]]:
        """
        Rank candidates by quality score.

        Args:
            candidates: List of primer sequences
            verbose: Print ranking statistics

        Returns:
            List of (primer, quality_score) tuples, sorted by score
        """
        # FIXED: Handle empty candidate list
        if not candidates:
            if verbose:
                logger.info("\nWarning: Empty candidate list provided to rank_by_quality")
            return []

        scorer = create_quality_scorer("moderate", self.conditions)

        scores = []
        for primer in candidates:
            try:
                quality = scorer.score_primer(primer)
                scores.append((primer, quality.overall_score))
            except (KeyError, ValueError, AttributeError, TypeError):
                scores.append((primer, 0.5))  # Default score if scoring fails

        scores.sort(key=lambda x: x[1], reverse=True)

        if verbose and scores:  # FIXED: Check scores is not empty
            logger.info(f"\nQuality ranking:")
            logger.info(f"  Top score: {scores[0][1]:.3f}")
            logger.info(f"  Median score: {scores[len(scores)//2][1]:.3f}")
            logger.info(f"  Bottom score: {scores[-1][1]:.3f}")

        return scores


def save_qa_report(result: QAFilterResult, output_path: Path, verbose: bool = True):
    """
    Save QA filtering report to file.

    Args:
        result: QAFilterResult from apply_post_step2_qa_filter()
        output_path: Path to save report
        verbose: Print save confirmation
    """
    report_lines = [
        "NEOSWGA QUALITY ASSURANCE FILTERING REPORT",
        "=" * 80,
        "",
        f"Input primers: {result.primers_in}",
        f"Output primers: {result.primers_out}",
        f"Filtered: {result.primers_filtered} ({100*result.primers_filtered/result.primers_in:.1f}%)",
        f"Mean QA score: {result.mean_qa_score:.3f}",
        "",
        "Filter breakdown:",
    ]

    for reason, count in result.filter_reasons.items():
        pct = 100 * count / result.primers_in if result.primers_in > 0 else 0
        report_lines.append(f"  {reason}: {count} ({pct:.1f}%)")

    report_lines.extend(["", "Top 10 primers by QA score:", ""])

    top_primers = result.filtered_df.nlargest(10, "qa_score")
    seq_col = primer_column(top_primers)
    for idx, row in top_primers.iterrows():
        report_lines.append(f"  {row[seq_col]}: {row['qa_score']:.3f}")

    output_path = Path(output_path)
    output_path.write_text("\n".join(report_lines))

    if verbose:
        logger.info(f"\nQA report saved to: {output_path}")
