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

The `--enable-qa` CLI flag reaches these through the three step wrappers at the
end of this module (apply_qa_to_step2_output, apply_qa_to_step3_output,
qa_prefiltered_candidates), which own the CSV read/write so the step handlers
stay a single call each.
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


def rf_score_column(df: pd.DataFrame) -> str:
    """Return the name of the random-forest score column.

    step3_df.csv names it ``'on.target.pred'``; ``'score'`` is the name used by
    hand-built frames and by the older QA examples.

    Raises:
        KeyError: if neither column is present.
    """
    for name in ("score", "on.target.pred"):
        if name in df.columns:
            return name
    raise KeyError(
        "DataFrame has neither a 'score' nor an 'on.target.pred' column; "
        f"columns present: {list(df.columns)}"
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
        step3_df: DataFrame from Step 3 with a primer column and an RF score
            column (see `rf_score_column`)
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
    score_col = rf_score_column(df)

    # Add QA scores
    df["qa_score"] = df[seq_col].map(qa_scores).fillna(0.5)

    # Normalize RF scores to 0-1 range
    rf_min = df[score_col].min()
    rf_max = df[score_col].max()
    if rf_max > rf_min:
        df["rf_score_norm"] = (df[score_col] - rf_min) / (rf_max - rf_min)
    else:
        df["rf_score_norm"] = 1.0

    # Compute composite score
    df["composite_score"] = rf_weight * df["rf_score_norm"] + qa_weight * df["qa_score"]

    if verbose:
        logger.info(f"\nScore statistics:")
        logger.info(f"  RF scores ({score_col}): {rf_min:.2f} - {rf_max:.2f}")
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


def _resolve_conditions(conditions: Optional[ReactionConditions]) -> ReactionConditions:
    """Fall back to the reaction the run is configured for, not the defaults.

    `ReactionConditions()` would silently score primers at 30 C with no
    additives regardless of what the run set, so QA would disagree with the
    filter and optimizer steps it sits between.
    """
    if conditions is not None:
        return conditions

    from neoswga.core.reaction_conditions import build_reaction_conditions

    return build_reaction_conditions()


def _read_step_csv(path: Path) -> Tuple[pd.DataFrame, bool]:
    """Read a pipeline CSV and report whether it carries a written index.

    step2_df.csv is written with an unnamed positional index and step3_df.csv
    with the primer as a named index. Returning the flag lets the caller write
    the file back in the shape it found it.
    """
    df = pd.read_csv(path)
    has_index = str(df.columns[0]).startswith("Unnamed:")
    if has_index:
        df = pd.read_csv(path, index_col=0)
    return df, has_index


def apply_qa_to_step2_output(
    data_dir: str,
    config: Optional[QAFilterConfig] = None,
    conditions: Optional[ReactionConditions] = None,
    verbose: bool = True,
) -> QAFilterResult:
    """Filter `step2_df.csv` in place with the QA metrics, and write a report.

    This is what `neoswga filter --enable-qa` runs after the standard step-2
    filtering. The `qa_score` column is kept in the rewritten CSV so
    `score --enable-qa` can blend it with the random-forest score.

    Args:
        data_dir: Directory holding step2_df.csv.
        config: QA filter configuration (default: moderate stringency).
        conditions: Reaction conditions (default: the run's configured ones).
        verbose: Print filtering progress.

    Returns:
        The QAFilterResult for the filtering that was applied.

    Raises:
        ValueError: if QA rejects every candidate. step2_df.csv is left as it
            was, because an empty one fails step 3 with a less specific error.
    """
    step2_path = Path(data_dir) / "step2_df.csv"
    step2_df, has_index = _read_step_csv(step2_path)

    result = apply_post_step2_qa_filter(
        step2_df,
        config=config,
        conditions=_resolve_conditions(conditions),
        verbose=verbose,
    )

    if result.primers_out == 0:
        raise ValueError(
            f"QA filtering rejected all {result.primers_in} candidate primers "
            f"({result.filter_reasons}). Relax the QA stringency or drop "
            f"--enable-qa; {step2_path} is unchanged."
        )

    filtered = result.filtered_df.drop(columns=["qa_pass", "filter_reason"], errors="ignore")
    filtered.to_csv(step2_path, index=has_index)
    save_qa_report(result, Path(data_dir) / "qa_report.txt", verbose=verbose)
    _record_qa_funnel_stage(Path(data_dir) / "filter_stats.json", result.primers_out)

    logger.info(
        f"QA filtering kept {result.primers_out}/{result.primers_in} candidates "
        f"(mean QA score {result.mean_qa_score:.3f})"
    )
    return result


def _record_qa_funnel_stage(stats_path: Path, primers_out: int) -> None:
    """Point the funnel's last stage at the pool QA actually left behind.

    The filter step writes filter_stats.json before QA runs, so its
    `final_candidates` counts primers QA has since removed. `after_qa` is
    recorded alongside for provenance; the report's FilteringStats ignores
    fields it does not declare.
    """
    if not stats_path.is_file():
        return
    import json

    try:
        stats = json.loads(stats_path.read_text())
        stats["after_qa"] = primers_out
        stats["final_candidates"] = primers_out
        stats_path.write_text(json.dumps(stats, indent=2))
    except (json.JSONDecodeError, OSError, TypeError) as e:
        logger.debug(f"Could not update filter_stats.json after QA: {e}")


def apply_qa_to_step3_output(
    data_dir: str,
    config: Optional[QAFilterConfig] = None,
    conditions: Optional[ReactionConditions] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Blend QA scores into `step3_df.csv` in place.

    This is what `neoswga score --enable-qa` runs after the random-forest
    scoring. QA scores come from step2_df.csv when `filter --enable-qa`
    produced them, and are computed here otherwise, so the flag is usable on
    the score step alone.

    The rewritten pool is ordered by `composite_score` because step 4 reads its
    candidates in file order.

    Args:
        data_dir: Directory holding step2_df.csv and step3_df.csv.
        config: QA configuration supplying the RF/QA score weights.
        conditions: Reaction conditions (default: the run's configured ones).
        verbose: Print scoring statistics.

    Returns:
        The step-3 DataFrame with `qa_score` and `composite_score` added.
    """
    if config is None:
        config = QAFilterConfig()
    conditions = _resolve_conditions(conditions)

    step3_path = Path(data_dir) / "step3_df.csv"
    step3_df, has_index = _read_step_csv(step3_path)
    primers = step3_df[primer_column(step3_df)].astype(str).tolist()

    qa_scores = _step2_qa_scores(Path(data_dir) / "step2_df.csv")
    if not qa_scores:
        logger.info("No QA scores in step2_df.csv; scoring the step-3 pool directly")
        optimizer = QAAwareOptimizer(config, conditions)
        qa_scores = dict(optimizer.rank_by_quality(primers, verbose=verbose))

    combined = combine_rf_qa_scores(
        step3_df,
        qa_scores,
        rf_weight=config.rf_weight,
        qa_weight=config.qa_weight,
        verbose=verbose,
    )
    combined = combined.sort_values("composite_score", ascending=False)
    combined.to_csv(step3_path, index=has_index)
    return combined


def _step2_qa_scores(step2_path: Path) -> Dict[str, float]:
    """QA scores carried over from `filter --enable-qa`, empty if there are none."""
    if not step2_path.is_file():
        return {}
    step2_df = pd.read_csv(step2_path)
    if "qa_score" not in step2_df.columns:
        return {}
    seq_col = primer_column(step2_df)
    return dict(zip(step2_df[seq_col].astype(str), step2_df["qa_score"]))


def qa_prefiltered_candidates(
    parameter_module,
    config: Optional[QAFilterConfig] = None,
    conditions: Optional[ReactionConditions] = None,
    verbose: bool = True,
) -> Optional[List[str]]:
    """Step-4 candidate pool with dimer-hub primers removed, or None.

    None tells the optimizer to load its own pool: either `--enable-qa` was not
    given, step3_df.csv is missing, or the pre-filter left nothing to optimize.
    The pre-filter is pairwise, so it costs O(n^2) dimer calculations and runs
    only when asked for.

    Args:
        parameter_module: The pipeline parameter module (read for `enable_qa`
            and `data_dir`).
        config: QA configuration supplying `max_hub_degree`.
        conditions: Reaction conditions (default: the run's configured ones).
        verbose: Print pre-filtering statistics.

    Returns:
        The retained candidate sequences, or None.
    """
    if not getattr(parameter_module, "enable_qa", False):
        return None

    step3_path = Path(getattr(parameter_module, "data_dir", ".")) / "step3_df.csv"
    if not step3_path.is_file():
        logger.warning(f"QA pre-filter skipped: {step3_path} not found")
        return None

    step3_df = pd.read_csv(step3_path)
    primers = step3_df[primer_column(step3_df)].astype(str).tolist()

    optimizer = create_qa_aware_optimizer(config, _resolve_conditions(conditions))
    filtered, _metadata = optimizer.prefilter_candidates(primers, verbose=verbose)

    if not filtered:
        logger.warning(
            f"QA pre-filter rejected all {len(primers)} candidates; "
            "optimizing over the unfiltered pool instead"
        )
        return None

    logger.info(f"QA pre-filter: {len(filtered)}/{len(primers)} candidates retained")
    return filtered
