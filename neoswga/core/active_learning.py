"""
Active Learning Framework for iterative primer set optimization.

Integrates experimental feedback into the design loop:
1. Propose primer sets using current model
2. Test experimentally (enrichment, uniformity)
3. Update model with feedback
4. Refine and propose better sets

Uses Bayesian optimization to efficiently explore design space
and converge to optimal solutions with minimal experiments.
"""

import numpy as np
import networkx as nx
import logging
from typing import List, Dict, Tuple, Optional, Callable
from dataclasses import dataclass, field
from datetime import datetime
import json
import os

logger = logging.getLogger(__name__)

try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, Matern
    HAS_GP = True
except ImportError:
    HAS_GP = False
    logger.warning("scikit-learn not installed. Active learning unavailable.")


@dataclass
class ExperimentalResult:
    """
    Record of experimental primer set performance.

    Stores measured outcomes from actual amplification experiments
    for model training and validation.
    """
    primer_set: List[str]
    timestamp: str

    # Primary outcomes
    enrichment_fold: float  # Target/background ratio
    uniformity_score: float  # Coefficient of variation (lower better)

    # Secondary metrics
    total_amplification: Optional[float] = None  # Total DNA yield
    off_target_fraction: Optional[float] = None  # Background contamination

    # Experimental conditions
    temperature: float = 30.0  # Amplification temperature (C)
    time_hours: float = 4.0  # Amplification time
    primer_concentration: float = 1.0  # uM

    # Quality indicators
    replicate_cv: Optional[float] = None  # Coefficient of variation across replicates
    passed_qc: bool = True  # Whether result is reliable

    # Notes
    notes: str = ""

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization"""
        return {
            'primer_set': self.primer_set,
            'timestamp': self.timestamp,
            'enrichment_fold': self.enrichment_fold,
            'uniformity_score': self.uniformity_score,
            'total_amplification': self.total_amplification,
            'off_target_fraction': self.off_target_fraction,
            'temperature': self.temperature,
            'time_hours': self.time_hours,
            'primer_concentration': self.primer_concentration,
            'replicate_cv': self.replicate_cv,
            'passed_qc': self.passed_qc,
            'notes': self.notes
        }

    @classmethod
    def from_dict(cls, data: Dict):
        """Load from dictionary"""
        return cls(**data)


@dataclass
class PrimerSetFeatures:
    """
    Feature representation of primer set for machine learning.

    Extracts features that correlate with experimental performance
    for use in Gaussian Process regression.
    """
    # Network features
    n_primers: int
    network_connectivity: float  # Average node degree
    n_components: int  # Number of connected components
    largest_component_size: int

    # Coverage features
    target_coverage: float  # Fraction of genome covered
    coverage_uniformity: float  # 1 - CV of coverage

    # Primer properties
    mean_tm: float  # Average melting temperature
    std_tm: float  # Std dev of Tm
    mean_gc: float  # Average GC content

    # Background binding
    background_sites: int  # Total off-target sites

    def to_vector(self) -> np.ndarray:
        """Convert to feature vector for GP"""
        return np.array([
            self.n_primers,
            self.network_connectivity,
            self.n_components,
            self.largest_component_size,
            self.target_coverage,
            self.coverage_uniformity,
            self.mean_tm,
            self.std_tm,
            self.mean_gc,
            np.log1p(self.background_sites)  # Log transform
        ])

    @staticmethod
    def get_feature_names() -> List[str]:
        """Get names of features"""
        return [
            'n_primers', 'network_connectivity', 'n_components',
            'largest_component_size', 'target_coverage', 'coverage_uniformity',
            'mean_tm', 'std_tm', 'mean_gc', 'log_background_sites'
        ]


class BayesianOptimizer:
    """
    Bayesian optimization for primer set selection.

    Uses Gaussian Process to model relationship between primer set
    features and experimental performance. Guides selection of next
    experiments using acquisition function (expected improvement).
    """

    def __init__(self, n_features: int = 10):
        """
        Initialize Bayesian optimizer.

        Args:
            n_features: Number of features in representation
        """
        if not HAS_GP:
            raise ImportError(
                "scikit-learn required for Bayesian optimization. "
                "Install: pip install scikit-learn"
            )

        self.n_features = n_features

        # Gaussian Process with RBF kernel
        kernel = ConstantKernel(1.0) * Matern(
            length_scale=np.ones(n_features),
            nu=2.5
        )

        self.gp = GaussianProcessRegressor(
            kernel=kernel,
            alpha=0.1,  # Noise level
            n_restarts_optimizer=10,
            normalize_y=True
        )

        # Training data
        self.X_train = []
        self.y_train = []

        # Best observed
        self.best_score = -np.inf
        self.best_features = None

    def add_observation(self, features: PrimerSetFeatures, score: float):
        """
        Add experimental observation.

        Args:
            features: Feature representation of primer set
            score: Measured performance (e.g., enrichment fold)
        """
        x = features.to_vector()
        self.X_train.append(x)
        self.y_train.append(score)

        # Update best
        if score > self.best_score:
            self.best_score = score
            self.best_features = features

        # Retrain GP
        if len(self.X_train) >= 3:  # Need minimum observations
            X = np.array(self.X_train)
            y = np.array(self.y_train)
            self.gp.fit(X, y)

    def predict(self, features: PrimerSetFeatures) -> Tuple[float, float]:
        """
        Predict performance for primer set.

        Args:
            features: Feature representation

        Returns:
            (mean, std): Predicted mean and uncertainty
        """
        if len(self.X_train) < 3:
            # Not enough data - return uninformative prior
            return 0.0, 1.0

        x = features.to_vector().reshape(1, -1)
        mean, std = self.gp.predict(x, return_std=True)

        return float(mean[0]), float(std[0])

    def acquisition_function(self, features: PrimerSetFeatures,
                           exploration_weight: float = 2.0) -> float:
        """
        Upper confidence bound acquisition function.

        UCB = mean + exploration_weight * std

        Balances exploitation (mean) and exploration (std).

        Args:
            features: Feature representation
            exploration_weight: Weight for exploration term

        Returns:
            Acquisition value (higher is better)
        """
        mean, std = self.predict(features)
        return mean + exploration_weight * std

    def select_next_experiment(self, candidates: List[PrimerSetFeatures],
                              exploration_weight: float = 2.0) -> int:
        """
        Select most promising candidate for next experiment.

        Uses acquisition function to balance exploration and exploitation.

        Args:
            candidates: List of candidate primer set features
            exploration_weight: Weight for exploration

        Returns:
            Index of selected candidate
        """
        if len(self.X_train) < 3:
            # Not enough data - random selection
            return np.random.randint(len(candidates))

        # Evaluate acquisition function
        acq_values = []
        for features in candidates:
            acq = self.acquisition_function(features, exploration_weight)
            acq_values.append(acq)

        # Select best
        return int(np.argmax(acq_values))


class ActiveLearningLoop:
    """
    Coordinates iterative primer design with experimental feedback.

    Workflow:
    1. Generate candidate primer sets using current optimizer
    2. Select most informative set using Bayesian optimization
    3. User tests experimentally
    4. Feed results back to update model
    5. Repeat until convergence
    """

    def __init__(self, cache, fg_prefixes: List[str], bg_prefixes: List[str],
                 fg_seq_lengths: List[int], bg_seq_lengths: List[int],
                 results_dir: str = "./active_learning_results"):
        """
        Initialize active learning loop.

        Args:
            cache: PositionCache with primer positions
            fg_prefixes: Target genome prefixes
            bg_prefixes: Background genome prefixes
            fg_seq_lengths: Target genome lengths
            bg_seq_lengths: Background genome lengths
            results_dir: Directory for storing results
        """
        self.cache = cache
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_seq_lengths = bg_seq_lengths
        self.results_dir = results_dir

        # Create results directory
        os.makedirs(results_dir, exist_ok=True)

        # Bayesian optimizer
        self.optimizer = BayesianOptimizer(n_features=10)

        # Experimental history
        self.experiments: List[ExperimentalResult] = []

        # Load previous results if exist
        self._load_history()

    def extract_features(self, primer_set: List[str]) -> PrimerSetFeatures:
        """
        Extract features from primer set.

        Args:
            primer_set: List of primer sequences

        Returns:
            Feature representation
        """
        from .network_optimizer import AmplificationNetwork
        import neoswga.core.primer_attributes as pa

        # Build network
        network = AmplificationNetwork(max_extension=70000)

        for primer in primer_set:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')

                for pos in pos_fwd:
                    network.add_binding_site(primer, prefix, int(pos), '+')
                for pos in pos_rev:
                    network.add_binding_site(primer, prefix, int(pos), '-')

        network.build_amplification_edges()

        # Network features
        if len(network.graph.nodes()) > 0:
            degrees = [network.graph.degree(n) for n in network.graph.nodes()]
            network_connectivity = np.mean(degrees) if degrees else 0.0
            n_components = network.num_components()

            # Largest component
            if n_components > 0:
                components = list(nx.connected_components(network.graph))
                largest_component_size = len(max(components, key=len))
            else:
                largest_component_size = 0
        else:
            network_connectivity = 0.0
            n_components = 0
            largest_component_size = 0

        # Coverage
        covered_bases = 0
        total_length = sum(self.fg_seq_lengths)

        for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
            positions = []
            for primer in primer_set:
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')
                positions.extend(pos_fwd)
                positions.extend(pos_rev)

            if positions:
                covered_bases += len(set(positions))

        target_coverage = covered_bases / total_length if total_length > 0 else 0.0

        # Primer properties
        tms = []
        gcs = []
        for primer in primer_set:
            tm = pa.calc_tm(primer)
            gc = pa.calc_gc(primer)
            tms.append(tm)
            gcs.append(gc)

        mean_tm = np.mean(tms) if tms else 0.0
        std_tm = np.std(tms) if tms else 0.0
        mean_gc = np.mean(gcs) if gcs else 0.0

        # Background binding
        background_sites = 0
        for prefix in self.bg_prefixes:
            for primer in primer_set:
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')
                background_sites += len(pos_fwd) + len(pos_rev)

        return PrimerSetFeatures(
            n_primers=len(primer_set),
            network_connectivity=network_connectivity,
            n_components=n_components,
            largest_component_size=largest_component_size,
            target_coverage=target_coverage,
            coverage_uniformity=0.8,  # Placeholder - would compute from bins
            mean_tm=mean_tm,
            std_tm=std_tm,
            mean_gc=mean_gc,
            background_sites=background_sites
        )

    def propose_candidates(self, candidates: List[str], n_sets: int = 10,
                          max_primers: int = 15) -> List[Tuple[List[str], PrimerSetFeatures]]:
        """
        Generate candidate primer sets for evaluation.

        Uses greedy network optimization with randomization to
        generate diverse candidates.

        Args:
            candidates: Pool of candidate primers
            n_sets: Number of sets to generate
            max_primers: Max primers per set

        Returns:
            List of (primer_set, features) tuples
        """
        from .network_optimizer import NetworkOptimizer

        optimizer = NetworkOptimizer(
            cache=self.cache,
            fg_prefixes=self.fg_prefixes,
            bg_prefixes=self.bg_prefixes,
            fg_seq_lengths=self.fg_seq_lengths,
            bg_seq_lengths=self.bg_seq_lengths
        )

        primer_sets = []

        for i in range(n_sets):
            # Randomize candidate order for diversity
            shuffled = candidates.copy()
            np.random.shuffle(shuffled)

            # Take subset
            subset = shuffled[:min(100, len(shuffled))]

            # Optimize
            result = optimizer.optimize(
                primers=subset,
                num_primers=max_primers,
                verbose=False
            )

            primer_set = result['selected_primers']
            features = self.extract_features(primer_set)

            primer_sets.append((primer_set, features))

        return primer_sets

    def recommend_next_experiment(self, candidates: List[str],
                                 n_candidates: int = 10,
                                 max_primers: int = 15,
                                 exploration_weight: float = 2.0) -> Dict:
        """
        Recommend next primer set to test experimentally.

        Args:
            candidates: Pool of candidate primers
            n_candidates: Number of candidates to generate
            max_primers: Max primers per set
            exploration_weight: Balance exploration vs exploitation

        Returns:
            Dictionary with recommended set and predictions
        """
        logger.info(f"Generating {n_candidates} candidate primer sets...")

        # Generate candidates
        primer_sets = self.propose_candidates(
            candidates=candidates,
            n_sets=n_candidates,
            max_primers=max_primers
        )

        # Extract features
        features_list = [f for _, f in primer_sets]

        # Select best using acquisition function
        logger.info("Selecting most informative candidate...")
        best_idx = self.optimizer.select_next_experiment(
            features_list,
            exploration_weight=exploration_weight
        )

        best_set, best_features = primer_sets[best_idx]

        # Get prediction
        pred_mean, pred_std = self.optimizer.predict(best_features)

        result = {
            'primer_set': best_set,
            'features': best_features,
            'predicted_enrichment': pred_mean,
            'prediction_uncertainty': pred_std,
            'n_experiments_so_far': len(self.experiments),
            'acquisition_value': self.optimizer.acquisition_function(
                best_features, exploration_weight
            )
        }

        logger.info(f"Recommended set: {len(best_set)} primers")
        logger.info(f"Predicted enrichment: {pred_mean:.1f} ± {pred_std:.1f}x")

        return result

    def add_experimental_result(self, result: ExperimentalResult):
        """
        Add experimental result and update model.

        Args:
            result: Experimental outcome
        """
        self.experiments.append(result)

        # Extract features
        features = self.extract_features(result.primer_set)

        # Update Bayesian optimizer
        score = result.enrichment_fold  # Primary objective
        self.optimizer.add_observation(features, score)

        # Save to disk
        self._save_history()

        logger.info(f"Added result: {result.enrichment_fold:.1f}x enrichment")
        logger.info(f"Total experiments: {len(self.experiments)}")

    def _save_history(self):
        """Save experimental history to disk"""
        history_file = os.path.join(self.results_dir, 'experiments.json')

        data = {
            'experiments': [exp.to_dict() for exp in self.experiments],
            'best_score': float(self.optimizer.best_score),
            'n_experiments': len(self.experiments)
        }

        with open(history_file, 'w') as f:
            json.dump(data, f, indent=2)

    def _load_history(self):
        """Load experimental history from disk"""
        history_file = os.path.join(self.results_dir, 'experiments.json')

        if not os.path.exists(history_file):
            return

        try:
            with open(history_file, 'r') as f:
                data = json.load(f)

            self.experiments = [
                ExperimentalResult.from_dict(exp)
                for exp in data['experiments']
            ]

            # Rebuild optimizer
            for exp in self.experiments:
                if exp.passed_qc:
                    features = self.extract_features(exp.primer_set)
                    self.optimizer.add_observation(features, exp.enrichment_fold)

            logger.info(f"Loaded {len(self.experiments)} previous experiments")

        except Exception as e:
            logger.warning(f"Could not load history: {e}")


if __name__ == "__main__":
    print("Active Learning Framework for Primer Optimization")
    print("\nWorkflow:")
    print("  1. Generate candidate primer sets")
    print("  2. Select most informative set (Bayesian optimization)")
    print("  3. Test experimentally")
    print("  4. Update model with results")
    print("  5. Repeat until convergence")
    print("\nFeatures:")
    print("  - Gaussian Process regression")
    print("  - Upper confidence bound acquisition")
    print("  - Balances exploration and exploitation")
    print("  - Persistent experimental history")
    print("\nUsage:")
    print("  from neoswga.core.active_learning import ActiveLearningLoop")
    print("  loop = ActiveLearningLoop(cache, fg_prefixes, bg_prefixes,")
    print("                           fg_seq_lengths, bg_seq_lengths)")
    print("  recommendation = loop.recommend_next_experiment(candidates)")
    print("  # Test recommendation['primer_set'] experimentally")
    print("  # Then add result:")
    print("  result = ExperimentalResult(")
    print("      primer_set=recommendation['primer_set'],")
    print("      enrichment_fold=1500.0,  # Measured")
    print("      uniformity_score=0.15")
    print("  )")
    print("  loop.add_experimental_result(result)")
