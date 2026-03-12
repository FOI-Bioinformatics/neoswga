"""
Unified efficiency prediction for SWGA primer sets.

Combines multiple prediction methods into a single confidence score:
- ML score: Random forest amplification prediction
- Network score: Amplification network connectivity
- Simulation score: Agent-based replication simulation (optional)

Provides clear recommendation: SYNTHESIZE / CAUTION / DO_NOT_SYNTHESIZE

Usage:
    from neoswga.core.efficiency_predictor import EfficiencyPredictor

    predictor = EfficiencyPredictor(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
    )

    result = predictor.predict(
        primers=['ATCGATCG', 'GCTAGCTA'],
        run_simulation=True,
        genome_sequence=sequence,
    )

    print(f"Confidence: {result.confidence_score:.2f} ({result.recommendation})")
    print(f"Predicted enrichment: {result.predicted_enrichment:.0f}x")
"""

import logging
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from enum import Enum
import numpy as np

logger = logging.getLogger(__name__)


class Recommendation(Enum):
    """Synthesis recommendation based on confidence score."""
    SYNTHESIZE = "SYNTHESIZE"
    CAUTION = "CAUTION"
    DO_NOT_SYNTHESIZE = "DO_NOT_SYNTHESIZE"


@dataclass
class EfficiencyPrediction:
    """
    Unified efficiency prediction for a primer set.

    Attributes:
        confidence_score: 0-1 composite score (higher is better)
        recommendation: SYNTHESIZE / CAUTION / DO_NOT_SYNTHESIZE
        ml_score: Normalized ML model prediction (0-1)
        network_score: Amplification network connectivity (0-1)
        simulation_score: Simulation-based prediction (0-1, optional)
        predicted_enrichment: Estimated enrichment factor (fold)
        predicted_coverage: Expected genome coverage fraction
        limiting_factors: Issues that limit predicted performance
        improvement_suggestions: Suggestions to improve the set
    """
    confidence_score: float
    recommendation: Recommendation
    ml_score: float
    network_score: float
    simulation_score: Optional[float]
    predicted_enrichment: float
    predicted_coverage: float
    limiting_factors: List[str]
    improvement_suggestions: List[str]

    # Confidence interval for enrichment
    enrichment_ci_low: float = 0.0
    enrichment_ci_high: float = 0.0

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""

        def _safe_float(value):
            """Convert value to JSON-safe float (handles inf, -inf, nan)."""
            if value is None:
                return None
            import math
            if math.isnan(value):
                return None
            if math.isinf(value):
                return 1e308 if value > 0 else -1e308  # Max JSON float
            return float(value)

        return {
            'confidence_score': _safe_float(self.confidence_score),
            'recommendation': self.recommendation.value,
            'ml_score': _safe_float(self.ml_score),
            'network_score': _safe_float(self.network_score),
            'simulation_score': _safe_float(self.simulation_score),
            'predicted_enrichment': _safe_float(self.predicted_enrichment),
            'enrichment_ci': [
                _safe_float(self.enrichment_ci_low),
                _safe_float(self.enrichment_ci_high)
            ],
            'predicted_coverage': _safe_float(self.predicted_coverage),
            'limiting_factors': self.limiting_factors,
            'improvement_suggestions': self.improvement_suggestions,
        }

    def __str__(self) -> str:
        lines = [
            f"Efficiency Prediction: {self.confidence_score:.2f} ({self.recommendation.value})",
            f"",
            f"Component Scores:",
            f"  ML score: {self.ml_score:.2f}",
            f"  Network score: {self.network_score:.2f}",
        ]
        if self.simulation_score is not None:
            lines.append(f"  Simulation score: {self.simulation_score:.2f}")

        lines.extend([
            f"",
            f"Predictions:",
            f"  Enrichment: {self.predicted_enrichment:.0f}x "
            f"(95% CI: {self.enrichment_ci_low:.0f}-{self.enrichment_ci_high:.0f}x)",
            f"  Coverage: {self.predicted_coverage:.1%}",
        ])

        if self.limiting_factors:
            lines.append(f"")
            lines.append(f"Limiting factors:")
            for factor in self.limiting_factors:
                lines.append(f"  - {factor}")

        if self.improvement_suggestions:
            lines.append(f"")
            lines.append(f"Suggestions:")
            for suggestion in self.improvement_suggestions:
                lines.append(f"  - {suggestion}")

        return "\n".join(lines)


class EfficiencyPredictor:
    """
    Unified efficiency predictor for SWGA primer sets.

    Combines ML, network, and simulation predictions into a
    single confidence score with clear synthesis recommendation.
    """

    # Thresholds for recommendations
    SYNTHESIZE_THRESHOLD = 0.7
    CAUTION_THRESHOLD = 0.4

    # Component weights for composite score
    WEIGHT_ML = 0.35
    WEIGHT_NETWORK = 0.35
    WEIGHT_SIMULATION = 0.30

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        max_extension: int = 70000,
    ):
        """
        Initialize efficiency predictor.

        Args:
            position_cache: PositionCache with primer binding positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers (optional)
            bg_seq_lengths: Background genome lengths (optional)
            max_extension: Maximum polymerase extension distance (bp)
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.max_extension = max_extension
        self.total_fg_length = sum(fg_seq_lengths)
        self.total_bg_length = sum(bg_seq_lengths) if bg_seq_lengths else 0

    def predict(
        self,
        primers: List[str],
        run_simulation: bool = False,
        genome_sequence: Optional[str] = None,
        simulation_replicates: int = 3,
        verbose: bool = True,
    ) -> EfficiencyPrediction:
        """
        Generate unified efficiency prediction for primer set.

        Args:
            primers: List of primer sequences
            run_simulation: Whether to run simulation (slower but more accurate)
            genome_sequence: Genome sequence (required if run_simulation=True)
            simulation_replicates: Number of simulation replicates
            verbose: Print progress

        Returns:
            EfficiencyPrediction with confidence score and details
        """
        if not primers:
            return self._empty_prediction("No primers provided")

        if verbose:
            logger.info("="*60)
            logger.info("EFFICIENCY PREDICTION")
            logger.info("="*60)
            logger.info(f"Analyzing {len(primers)} primers")

        limiting_factors = []
        suggestions = []

        # 1. ML Score - from random forest model
        if verbose:
            logger.info("\n1. ML-based prediction...")

        ml_score, ml_enrichment = self._calculate_ml_score(primers, verbose)

        if ml_score < 0.5:
            limiting_factors.append("Low ML amplification prediction")
            suggestions.append("Consider primers with higher on-target prediction scores")

        # 2. Network Score - amplification connectivity
        if verbose:
            logger.info("\n2. Network connectivity analysis...")

        network_score, network_stats = self._calculate_network_score(primers, verbose)

        if network_score < 0.5:
            limiting_factors.append("Poor amplification network connectivity")
            suggestions.append("Add primers to improve strand pairing opportunities")

        # 3. Coverage analysis
        if verbose:
            logger.info("\n3. Coverage analysis...")

        coverage, coverage_stats = self._calculate_coverage(primers, verbose)

        if coverage < 0.8:
            limiting_factors.append(f"Incomplete coverage ({coverage:.1%})")
            suggestions.append("Add primers to cover gap regions")

        # 4. Dimer risk analysis
        if verbose:
            logger.info("\n4. Dimer risk analysis...")

        dimer_risk, dimer_issues = self._calculate_dimer_risk(primers, verbose)

        if dimer_risk > 0.3:
            limiting_factors.append(f"High primer-dimer risk ({dimer_risk:.0%} of pairs)")
            suggestions.append("Consider removing primers with high dimer potential")

        # 5. Simulation (optional)
        simulation_score = None
        simulation_coverage = None

        if run_simulation:
            if genome_sequence is None:
                if verbose:
                    logger.warning("Simulation requested but no genome sequence provided")
            else:
                if verbose:
                    logger.info("\n5. Simulation validation...")

                simulation_score, simulation_coverage = self._run_simulation(
                    primers, genome_sequence, simulation_replicates, verbose
                )

                if simulation_score < 0.5:
                    limiting_factors.append("Poor simulated amplification")
                    suggestions.append("Primers may not cooperate well in practice")

        # Calculate composite confidence score
        if simulation_score is not None:
            # Use all three components
            confidence = (
                self.WEIGHT_ML * ml_score +
                self.WEIGHT_NETWORK * network_score +
                self.WEIGHT_SIMULATION * simulation_score
            )
        else:
            # Redistribute simulation weight
            ml_weight = self.WEIGHT_ML + self.WEIGHT_SIMULATION / 2
            network_weight = self.WEIGHT_NETWORK + self.WEIGHT_SIMULATION / 2
            confidence = ml_weight * ml_score + network_weight * network_score

        # Apply penalties
        confidence *= (1 - dimer_risk * 0.3)  # Dimer penalty
        confidence *= min(1.0, coverage / 0.8)  # Coverage penalty

        # Clamp to 0-1
        confidence = max(0.0, min(1.0, confidence))

        # Determine recommendation
        if confidence >= self.SYNTHESIZE_THRESHOLD:
            recommendation = Recommendation.SYNTHESIZE
        elif confidence >= self.CAUTION_THRESHOLD:
            recommendation = Recommendation.CAUTION
        else:
            recommendation = Recommendation.DO_NOT_SYNTHESIZE

        # Estimate enrichment with confidence interval
        predicted_enrichment = ml_enrichment * (1 + network_score)
        enrichment_ci_low = predicted_enrichment * 0.5
        enrichment_ci_high = predicted_enrichment * 2.0

        if simulation_coverage is not None:
            # Adjust based on simulation
            predicted_enrichment *= simulation_score

        if verbose:
            logger.info("\n" + "="*60)
            logger.info("PREDICTION SUMMARY")
            logger.info("="*60)
            logger.info(f"Confidence: {confidence:.2f} ({recommendation.value})")
            logger.info(f"Predicted enrichment: {predicted_enrichment:.0f}x")
            logger.info(f"Coverage: {coverage:.1%}")

        return EfficiencyPrediction(
            confidence_score=confidence,
            recommendation=recommendation,
            ml_score=ml_score,
            network_score=network_score,
            simulation_score=simulation_score,
            predicted_enrichment=predicted_enrichment,
            predicted_coverage=coverage,
            limiting_factors=limiting_factors,
            improvement_suggestions=suggestions,
            enrichment_ci_low=enrichment_ci_low,
            enrichment_ci_high=enrichment_ci_high,
        )

    def _empty_prediction(self, message: str) -> EfficiencyPrediction:
        """Create empty prediction for error cases."""
        return EfficiencyPrediction(
            confidence_score=0.0,
            recommendation=Recommendation.DO_NOT_SYNTHESIZE,
            ml_score=0.0,
            network_score=0.0,
            simulation_score=None,
            predicted_enrichment=1.0,
            predicted_coverage=0.0,
            limiting_factors=[message],
            improvement_suggestions=[],
        )

    def _calculate_ml_score(
        self, primers: List[str], verbose: bool
    ) -> Tuple[float, float]:
        """
        Calculate ML-based score using random forest model.

        Returns:
            Tuple of (normalized_score, predicted_enrichment)

        Note:
            Uses pickle deserialization for the pre-trained model.
            Only loads from the package's models directory for security.
        """
        try:
            import pickle
            import os
            from neoswga.core.rf_preprocessing import compute_primer_features

            # Load the random forest model from package directory only
            package_dir = os.path.dirname(os.path.abspath(__file__))
            model_path = os.path.join(package_dir, 'models', 'random_forest_filter.p')

            # Security: Verify model path is within package directory
            model_path = os.path.abspath(model_path)
            if not model_path.startswith(package_dir):
                logger.error("Model path traversal attempt blocked")
                return 0.5, 30.0

            if not os.path.exists(model_path):
                if verbose:
                    logger.warning("RF model not found, using estimate")
                return 0.6, 50.0  # Default estimate

            # Note: pickle.load is used here for sklearn model compatibility.
            # The model file is shipped with the package and should not be
            # modified by users. For user-provided models, consider using
            # a safer serialization format like ONNX.
            logger.warning(
                f"Loading pickle file from {model_path}. "
                "Only load files from trusted sources."
            )
            with open(model_path, 'rb') as f:
                model = pickle.load(f)

            # Get predictions for each primer
            scores = []
            for primer in primers:
                try:
                    features = compute_primer_features(primer)
                    pred = model.predict_proba([features])[0][1]  # Probability of good primer
                    scores.append(pred)
                except Exception:
                    scores.append(0.5)  # Default if features fail

            # Average score
            ml_score = np.mean(scores)

            # Estimate enrichment from score
            # Empirical: score 0.8 -> ~100x enrichment
            enrichment = 10 ** (ml_score * 2.5)

            if verbose:
                logger.info(f"  Mean ML score: {ml_score:.2f}")
                logger.info(f"  Score range: {min(scores):.2f} - {max(scores):.2f}")

            return ml_score, enrichment

        except Exception as e:
            if verbose:
                logger.warning(f"ML scoring failed: {e}")
            return 0.5, 30.0  # Default

    def _calculate_network_score(
        self, primers: List[str], verbose: bool
    ) -> Tuple[float, Dict]:
        """
        Calculate network connectivity score.

        Returns:
            Tuple of (normalized_score, stats_dict)
        """
        try:
            from neoswga.core.network_optimizer import AmplificationNetwork

            network = AmplificationNetwork(max_extension=self.max_extension)

            for primer in primers:
                for prefix in self.fg_prefixes:
                    positions_fwd = self.cache.get_positions(prefix, primer, 'forward')
                    positions_rev = self.cache.get_positions(prefix, primer, 'reverse')

                    if len(positions_fwd) > 0:
                        network.add_primer_sites(primer, positions_fwd, '+')
                    if len(positions_rev) > 0:
                        network.add_primer_sites(primer, positions_rev, '-')

            network.build_edges()
            stats = network.get_statistics()

            # Normalize connectivity to 0-1
            connectivity = stats.get('connectivity', 0)
            predicted_amp = stats.get('predicted_amplification', 1)

            # Score based on connectivity and predicted amplification
            network_score = min(1.0, connectivity * 0.5 + min(predicted_amp / 200, 0.5))

            if verbose:
                logger.info(f"  Connectivity: {connectivity:.2f}")
                logger.info(f"  Largest component: {stats.get('largest_component', 0)} sites")
                logger.info(f"  Predicted amplification: {predicted_amp:.0f}x")

            return network_score, stats

        except Exception as e:
            if verbose:
                logger.warning(f"Network analysis failed: {e}")
            return 0.5, {}

    def _calculate_coverage(
        self, primers: List[str], verbose: bool
    ) -> Tuple[float, Dict]:
        """
        Calculate genome coverage.

        Returns:
            Tuple of (coverage_fraction, stats_dict)
        """
        try:
            from neoswga.core.dominating_set_optimizer import BipartiteGraph

            bin_size = 10000
            graph = BipartiteGraph(bin_size=bin_size)

            for primer in primers:
                for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                    positions = self.cache.get_positions(prefix, primer, 'both')
                    if len(positions) > 0:
                        graph.add_primer_coverage(primer, positions, prefix, length)

            total_bins = sum(
                (length + bin_size - 1) // bin_size
                for length in self.fg_seq_lengths
            )

            coverage = len(graph.regions) / total_bins if total_bins > 0 else 0.0

            stats = {
                'covered_bins': len(graph.regions),
                'total_bins': total_bins,
            }

            if verbose:
                logger.info(f"  Coverage: {coverage:.1%}")
                logger.info(f"  Bins covered: {len(graph.regions)}/{total_bins}")

            return coverage, stats

        except Exception as e:
            if verbose:
                logger.warning(f"Coverage analysis failed: {e}")
            return 0.5, {}

    # Maximum primers to analyze for O(n^2) operations
    MAX_DIMER_ANALYSIS_PRIMERS = 50

    def _calculate_dimer_risk(
        self, primers: List[str], verbose: bool
    ) -> Tuple[float, List[str]]:
        """
        Calculate primer-dimer risk.

        For sets larger than MAX_DIMER_ANALYSIS_PRIMERS, samples a subset
        to bound computation time.

        Returns:
            Tuple of (risk_fraction, list_of_issues)
        """
        try:
            from neoswga.core import dimer
            import random

            n = len(primers)
            if n < 2:
                return 0.0, []

            # Bound O(n^2) by sampling for large sets
            if n > self.MAX_DIMER_ANALYSIS_PRIMERS:
                if verbose:
                    logger.info(
                        f"  Sampling {self.MAX_DIMER_ANALYSIS_PRIMERS} primers "
                        f"for dimer analysis (total: {n})"
                    )
                primers_to_check = random.sample(primers, self.MAX_DIMER_ANALYSIS_PRIMERS)
            else:
                primers_to_check = primers

            n_check = len(primers_to_check)
            dimer_pairs = 0
            total_pairs = 0
            issues = []

            for i in range(n_check):
                for j in range(i + 1, n_check):
                    total_pairs += 1
                    if dimer.is_dimer_fast(primers_to_check[i], primers_to_check[j], max_dimer_bp=4):
                        dimer_pairs += 1
                        if len(issues) < 3:  # Limit reported issues
                            issues.append(
                                f"{primers_to_check[i][:8]}... + {primers_to_check[j][:8]}..."
                            )

            risk = dimer_pairs / total_pairs if total_pairs > 0 else 0.0

            if verbose:
                logger.info(f"  Dimer risk: {risk:.0%} ({dimer_pairs}/{total_pairs} pairs)")

            return risk, issues

        except Exception as e:
            if verbose:
                logger.warning(f"Dimer analysis failed: {e}")
            return 0.1, []  # Default low risk

    def _run_simulation(
        self,
        primers: List[str],
        genome_sequence: str,
        replicates: int,
        verbose: bool,
    ) -> Tuple[float, float]:
        """
        Run simulation to validate primer set.

        Returns:
            Tuple of (normalized_score, simulated_coverage)
        """
        try:
            from neoswga.core.simulation_fitness import SimulationBasedEvaluator

            evaluator = SimulationBasedEvaluator(
                genome_sequence=genome_sequence,
                genome_length=len(genome_sequence),
                position_cache=self.cache,
                n_replicates=replicates,
            )

            fitness = evaluator.evaluate(primers, verbose=verbose)

            # Normalize fitness to 0-1
            sim_score = min(1.0, fitness.fitness_score)
            coverage = fitness.mean_coverage

            if verbose:
                logger.info(f"  Simulated coverage: {coverage:.1%}")
                logger.info(f"  Uniformity: {fitness.coverage_uniformity:.2f}")
                logger.info(f"  Fitness: {fitness.fitness_score:.3f}")

            return sim_score, coverage

        except Exception as e:
            if verbose:
                logger.warning(f"Simulation failed: {e}")
            return None, None


def predict_efficiency(
    params_path: str,
    primers: List[str],
    run_simulation: bool = False,
    output_path: Optional[str] = None,
    verbose: bool = True,
) -> EfficiencyPrediction:
    """
    Convenience function to predict efficiency using params.json.

    Args:
        params_path: Path to params.json
        primers: Primer sequences to evaluate
        run_simulation: Whether to run simulation
        output_path: Path to save JSON result (optional)
        verbose: Print progress

    Returns:
        EfficiencyPrediction with results
    """
    import os
    import json
    from neoswga.core.position_cache import PositionCache
    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline

    # Load parameters
    parameter.json_file = params_path
    core_pipeline._initialize()

    fg_prefixes = core_pipeline.fg_prefixes
    bg_prefixes = core_pipeline.bg_prefixes
    fg_seq_lengths = core_pipeline.fg_seq_lengths
    bg_seq_lengths = core_pipeline.bg_seq_lengths

    # Load position cache
    cache = PositionCache(fg_prefixes, primers)

    # Load genome sequence if simulation requested
    genome_sequence = None
    if run_simulation:
        fg_genomes = getattr(parameter, 'fg_genomes', [])
        if fg_genomes:
            from neoswga.core.genome_io import read_genome
            genome_sequence = read_genome(fg_genomes[0])

    # Create predictor
    predictor = EfficiencyPredictor(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bg_prefixes=bg_prefixes,
        bg_seq_lengths=bg_seq_lengths,
    )

    # Run prediction
    result = predictor.predict(
        primers=primers,
        run_simulation=run_simulation,
        genome_sequence=genome_sequence,
        verbose=verbose,
    )

    # Save result if output path specified
    if output_path:
        with open(output_path, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)

    return result
