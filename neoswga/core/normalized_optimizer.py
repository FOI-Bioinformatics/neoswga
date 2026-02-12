"""
Normalized primer set optimizer with strategy presets.

DESIGN RATIONALE:
The existing optimizers (network, hybrid, dominating-set, background-aware)
each optimize for different objectives using incompatible scoring scales.
This normalized optimizer provides:

1. NORMALIZED SCORING: All components (coverage, connectivity, amplification,
   background) are normalized to 0-1 scale before combination.

2. STRATEGY PRESETS: Pre-configured weight combinations for common use cases
   (clinical, discovery, fast, balanced) with clear documentation.

3. CUSTOMIZABLE WEIGHTS: Users can adjust individual component weights to
   tune optimization for specific applications.

4. TRANSPARENT METRICS: Each primer set evaluation returns all component
   scores, making optimization behavior interpretable.

This addresses two critical gaps identified in the codebase:
- Problem A: Algorithm gaps (e.g., hybrid optimizer mixes incompatible units)
- Problem D: No guidance for optimizer selection (users don't know which to use)
"""

import logging
import math
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np

from neoswga.core.network_optimizer import AmplificationNetwork

logger = logging.getLogger(__name__)


# =============================================================================
# Strategy Presets
# =============================================================================

STRATEGY_PRESETS = {
    'clinical': {
        'description': 'High specificity for diagnostic applications',
        'weights': {
            'coverage': 0.2,
            'connectivity': 0.2,
            'amplification': 0.2,
            'background': 0.4,  # Heavily penalize background
        },
        'recommended_for': 'Diagnostics, clinical samples, low false-positive tolerance',
    },
    'discovery': {
        'description': 'Maximum coverage for pathogen discovery',
        'weights': {
            'coverage': 0.4,
            'connectivity': 0.3,
            'amplification': 0.2,
            'background': 0.1,  # Accept some background for coverage
        },
        'recommended_for': 'Pathogen discovery, unknown samples, metagenomics',
    },
    'fast': {
        'description': 'Quick optimization prioritizing coverage',
        'weights': {
            'coverage': 0.5,
            'connectivity': 0.3,
            'amplification': 0.2,
            'background': 0.0,  # Ignore background for speed
        },
        'recommended_for': 'Quick screening, initial design, large primer pools',
    },
    'balanced': {
        'description': 'Equal weight across all objectives',
        'weights': {
            'coverage': 0.25,
            'connectivity': 0.25,
            'amplification': 0.25,
            'background': 0.25,
        },
        'recommended_for': 'General purpose, uncertain requirements',
    },
    'enrichment': {
        'description': 'Optimize for sequencing enrichment',
        'weights': {
            'coverage': 0.3,
            'connectivity': 0.25,
            'amplification': 0.3,
            'background': 0.15,
        },
        'recommended_for': 'Targeted sequencing, library prep, enrichment',
    },
}


@dataclass
class NormalizedScores:
    """
    Normalized scores for a primer set, all on 0-1 scale.

    Higher is always better for each component.
    """
    coverage: float        # 0-1: fraction of genome covered
    connectivity: float    # 0-1: normalized network connectivity
    amplification: float   # 0-1: normalized predicted amplification
    background: float      # 0-1: inverted background (1 = no background binding)

    # Raw values for inspection
    raw_coverage: float = 0.0
    raw_connectivity: float = 0.0
    raw_amplification: float = 0.0
    raw_background: float = 0.0

    def composite(self, weights: Dict[str, float]) -> float:
        """
        Weighted combination of normalized scores.

        Args:
            weights: Dictionary with keys 'coverage', 'connectivity',
                    'amplification', 'background'. Values should sum to ~1.0.

        Returns:
            Composite score in range 0-1.
        """
        return (
            weights.get('coverage', 0.25) * self.coverage +
            weights.get('connectivity', 0.25) * self.connectivity +
            weights.get('amplification', 0.25) * self.amplification +
            weights.get('background', 0.25) * self.background
        )

    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary for serialization."""
        return {
            'coverage': self.coverage,
            'connectivity': self.connectivity,
            'amplification': self.amplification,
            'background': self.background,
            'raw_coverage': self.raw_coverage,
            'raw_connectivity': self.raw_connectivity,
            'raw_amplification': self.raw_amplification,
            'raw_background': self.raw_background,
        }


@dataclass
class NormalizedOptimizationResult:
    """Result from normalized optimization."""
    primers: List[str]
    scores: NormalizedScores
    composite_score: float
    strategy: str
    weights: Dict[str, float]

    # Optional detailed metrics
    stage_history: List[Dict] = field(default_factory=list)
    runtime_seconds: float = 0.0
    candidates_evaluated: int = 0

    def summary(self) -> str:
        """Human-readable summary of optimization result."""
        lines = [
            "=" * 60,
            "NORMALIZED OPTIMIZATION RESULT",
            "=" * 60,
            f"Strategy: {self.strategy}",
            f"Primers selected: {len(self.primers)}",
            f"Composite score: {self.composite_score:.3f}",
            "",
            "Component Scores (normalized 0-1):",
            f"  Coverage:      {self.scores.coverage:.3f} (raw: {self.scores.raw_coverage:.2%})",
            f"  Connectivity:  {self.scores.connectivity:.3f} (raw: {self.scores.raw_connectivity:.3f})",
            f"  Amplification: {self.scores.amplification:.3f} (raw: {self.scores.raw_amplification:.1f}x)",
            f"  Background:    {self.scores.background:.3f} (raw: {self.scores.raw_background:.1f}x)",
            "",
            "Weights used:",
        ]
        for key, val in self.weights.items():
            lines.append(f"  {key}: {val:.2f}")

        lines.extend([
            "",
            f"Runtime: {self.runtime_seconds:.2f}s",
            f"Candidates evaluated: {self.candidates_evaluated}",
            "=" * 60,
        ])

        return "\n".join(lines)


class NormalizedScorer:
    """
    Normalized scoring for primer sets.

    All scores are normalized to 0-1 scale to enable meaningful
    weighted combination without dimensional inconsistency.
    """

    # Normalization constants (based on typical SWGA scenarios)
    MAX_CONNECTIVITY = 5.0      # Algebraic connectivity rarely exceeds this
    MAX_AMPLIFICATION = 1e6     # Realistic upper bound (2^20)
    MIN_AMPLIFICATION = 1.0     # No amplification

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        max_extension: int = 70000,
        bin_size: int = 10000,
    ):
        """
        Initialize normalized scorer.

        Args:
            position_cache: Cache with primer binding positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers
            bg_seq_lengths: Background genome lengths
            max_extension: Maximum polymerase extension distance (bp)
            bin_size: Bin size for coverage calculation (bp)
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.max_extension = max_extension
        self.bin_size = bin_size

        # Calculate total genome sizes
        self.total_fg_length = sum(fg_seq_lengths)
        self.total_bg_length = sum(bg_seq_lengths) if bg_seq_lengths else 0

        # Calculate total bins for coverage normalization
        self.total_fg_bins = sum(
            (length + bin_size - 1) // bin_size
            for length in fg_seq_lengths
        )

    def score(self, primers: List[str]) -> NormalizedScores:
        """
        Calculate normalized scores for a primer set.

        Args:
            primers: List of primer sequences

        Returns:
            NormalizedScores with all components on 0-1 scale.
        """
        if not primers:
            return NormalizedScores(
                coverage=0.0, connectivity=0.0,
                amplification=0.0, background=1.0
            )

        # Build networks
        fg_network = self._build_network(primers, self.fg_prefixes)
        bg_network = self._build_network(primers, self.bg_prefixes) if self.bg_prefixes else None

        # Calculate raw scores
        raw_coverage = self._calculate_coverage(primers)
        raw_connectivity = fg_network.connectivity_score()
        raw_amplification = fg_network.predict_amplification_fold()
        raw_background = bg_network.predict_amplification_fold() if bg_network else 1.0

        # Normalize to 0-1
        norm_coverage = raw_coverage  # Already 0-1
        norm_connectivity = min(1.0, raw_connectivity / self.MAX_CONNECTIVITY)
        norm_amplification = self._normalize_amplification(raw_amplification)
        norm_background = self._normalize_background(raw_background)

        return NormalizedScores(
            coverage=norm_coverage,
            connectivity=norm_connectivity,
            amplification=norm_amplification,
            background=norm_background,
            raw_coverage=raw_coverage,
            raw_connectivity=raw_connectivity,
            raw_amplification=raw_amplification,
            raw_background=raw_background,
        )

    def _normalize_amplification(self, raw: float) -> float:
        """
        Normalize amplification to 0-1 using log scale.

        Uses log scale because amplification spans many orders of magnitude.
        Maps: 1x -> 0.0, 1e6x -> 1.0
        """
        if raw <= self.MIN_AMPLIFICATION:
            return 0.0

        log_raw = math.log10(max(raw, 1.0))
        log_max = math.log10(self.MAX_AMPLIFICATION)

        return min(1.0, log_raw / log_max)

    def _normalize_background(self, raw: float) -> float:
        """
        Normalize background to 0-1, inverted (higher = less background).

        Uses same log scale as amplification, but inverted so that
        lower background binding gives higher score.
        """
        if raw <= self.MIN_AMPLIFICATION:
            return 1.0  # No background = perfect score

        log_raw = math.log10(max(raw, 1.0))
        log_max = math.log10(self.MAX_AMPLIFICATION)

        # Invert: high background -> low score
        return max(0.0, 1.0 - log_raw / log_max)

    def _build_network(
        self, primers: List[str], prefixes: List[str]
    ) -> AmplificationNetwork:
        """Build amplification network for primer set."""
        network = AmplificationNetwork(max_extension=self.max_extension)

        for primer in primers:
            for prefix in prefixes:
                # Use 'forward'/'reverse' for PositionCache API
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')

                if len(pos_fwd) > 0:
                    network.add_primer_sites(primer, pos_fwd, '+')
                if len(pos_rev) > 0:
                    network.add_primer_sites(primer, pos_rev, '-')

        network.build_edges()
        return network

    def _calculate_coverage(self, primers: List[str]) -> float:
        """Calculate genome coverage fraction (0-1)."""
        covered_bins = set()

        for primer in primers:
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                pos_fwd = self.cache.get_positions(prefix, primer, 'forward')
                pos_rev = self.cache.get_positions(prefix, primer, 'reverse')
                positions = np.concatenate([pos_fwd, pos_rev]) if len(pos_fwd) > 0 or len(pos_rev) > 0 else np.array([])

                # Mark bins as covered (extending by max_extension)
                for pos in positions:
                    # Coverage extends from position to position + max_extension
                    start_bin = int(pos) // self.bin_size
                    end_bin = (int(pos) + self.max_extension) // self.bin_size
                    max_bin = (length + self.bin_size - 1) // self.bin_size
                    for b in range(start_bin, min(end_bin + 1, max_bin)):
                        covered_bins.add((prefix, b))

        if self.total_fg_bins == 0:
            return 0.0

        return len(covered_bins) / self.total_fg_bins


class NormalizedOptimizer:
    """
    Normalized primer set optimizer with strategy presets.

    Replaces the need to choose between different optimizers by providing
    a single interface with configurable optimization goals.

    Usage:
        optimizer = NormalizedOptimizer(cache, fg_prefixes, fg_seq_lengths)
        result = optimizer.optimize(candidates, target_size=10, strategy='clinical')
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        max_extension: int = 70000,
        bin_size: int = 10000,
    ):
        """
        Initialize normalized optimizer.

        Args:
            position_cache: Cache with primer binding positions
            fg_prefixes: Target genome identifiers
            fg_seq_lengths: Target genome lengths
            bg_prefixes: Background genome identifiers
            bg_seq_lengths: Background genome lengths
            max_extension: Maximum polymerase extension distance (bp)
            bin_size: Bin size for coverage calculation (bp)
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.max_extension = max_extension
        self.bin_size = bin_size

        # Initialize scorer
        self.scorer = NormalizedScorer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            max_extension=max_extension,
            bin_size=bin_size,
        )

    def optimize(
        self,
        candidates: List[str],
        target_size: int = 10,
        strategy: str = 'balanced',
        custom_weights: Optional[Dict[str, float]] = None,
        fixed_primers: Optional[List[str]] = None,
        verbose: bool = True,
    ) -> NormalizedOptimizationResult:
        """
        Optimize primer set using greedy selection with normalized scoring.

        Args:
            candidates: List of candidate primers
            target_size: Number of primers to select
            strategy: Strategy preset name ('clinical', 'discovery', 'fast',
                     'balanced', 'enrichment') or 'custom' if using custom_weights
            custom_weights: Optional custom weight dictionary (overrides strategy)
            fixed_primers: Optional primers that must be included
            verbose: Print progress

        Returns:
            NormalizedOptimizationResult with selected primers and scores.
        """
        start_time = time.time()

        # Get weights
        if custom_weights:
            weights = custom_weights
            strategy = 'custom'
        elif strategy in STRATEGY_PRESETS:
            weights = STRATEGY_PRESETS[strategy]['weights']
        else:
            raise ValueError(
                f"Unknown strategy '{strategy}'. "
                f"Available: {list(STRATEGY_PRESETS.keys())}"
            )

        if verbose:
            logger.info("=" * 60)
            logger.info(f"NORMALIZED OPTIMIZATION (Strategy: {strategy})")
            logger.info("=" * 60)
            logger.info(f"Candidates: {len(candidates)}")
            logger.info(f"Target size: {target_size}")
            logger.info(f"Weights: {weights}")

        # Initialize with fixed primers
        fixed_primers = [p.upper() for p in (fixed_primers or [])]
        selected = list(fixed_primers)
        remaining = [c for c in candidates if c.upper() not in fixed_primers]

        if verbose and fixed_primers:
            logger.info(f"Fixed primers: {len(fixed_primers)}")

        # Greedy selection
        stage_history = []
        candidates_evaluated = 0

        while len(selected) < target_size and remaining:
            best_primer = None
            best_score = -float('inf')
            best_scores = None

            for primer in remaining:
                test_set = selected + [primer]
                scores = self.scorer.score(test_set)
                composite = scores.composite(weights)
                candidates_evaluated += 1

                if composite > best_score:
                    best_score = composite
                    best_primer = primer
                    best_scores = scores

            if best_primer is None:
                break

            selected.append(best_primer)
            remaining.remove(best_primer)

            # Record history
            stage_history.append({
                'iteration': len(selected),
                'primer_added': best_primer,
                'composite_score': best_score,
                'scores': best_scores.to_dict(),
            })

            if verbose:
                logger.info(
                    f"  [{len(selected)}/{target_size}] Added {best_primer[:8]}... "
                    f"score={best_score:.4f}"
                )

        # Final scoring
        final_scores = self.scorer.score(selected)
        final_composite = final_scores.composite(weights)

        runtime = time.time() - start_time

        if verbose:
            logger.info("-" * 60)
            logger.info(f"Optimization complete: {len(selected)} primers")
            logger.info(f"Final composite score: {final_composite:.4f}")
            logger.info(f"Runtime: {runtime:.2f}s")
            logger.info("=" * 60)

        return NormalizedOptimizationResult(
            primers=selected,
            scores=final_scores,
            composite_score=final_composite,
            strategy=strategy,
            weights=weights,
            stage_history=stage_history,
            runtime_seconds=runtime,
            candidates_evaluated=candidates_evaluated,
        )

    @staticmethod
    def list_strategies() -> Dict[str, Dict]:
        """List available optimization strategies with descriptions."""
        return {
            name: {
                'description': preset['description'],
                'recommended_for': preset['recommended_for'],
                'weights': preset['weights'],
            }
            for name, preset in STRATEGY_PRESETS.items()
        }

    @staticmethod
    def recommend_strategy(
        application: Optional[str] = None,
        background_available: bool = True,
        time_sensitive: bool = False,
    ) -> str:
        """
        Recommend an optimization strategy based on use case.

        Args:
            application: Application type ('clinical', 'discovery', 'sequencing')
            background_available: Whether background genome is available
            time_sensitive: Whether quick results are needed

        Returns:
            Recommended strategy name.
        """
        if time_sensitive:
            return 'fast'

        if not background_available:
            return 'fast'  # Can't optimize background without data

        if application == 'clinical':
            return 'clinical'
        elif application == 'discovery':
            return 'discovery'
        elif application == 'sequencing':
            return 'enrichment'

        return 'balanced'


# =============================================================================
# Factory Registration - BaseOptimizer Interface
# =============================================================================

from neoswga.core.base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig
)
from neoswga.core.optimizer_factory import OptimizerFactory


@OptimizerFactory.register('normalized', aliases=['normalized-optimizer', 'strategy-based'])
class NormalizedBaseOptimizer(BaseOptimizer):
    """
    Normalized optimizer implementing BaseOptimizer interface.

    Provides normalized scoring with configurable strategy presets.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        **kwargs
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config
        )
        self._normalized = NormalizedOptimizer(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            max_extension=kwargs.get('max_extension', 70000),
            bin_size=kwargs.get('bin_size', 10000),
        )
        self._strategy = kwargs.get('strategy', 'balanced')
        self._custom_weights = kwargs.get('custom_weights')

    @property
    def name(self) -> str:
        """Optimizer identifier for logging and factory registration."""
        return "normalized"

    @property
    def description(self) -> str:
        """One-line summary of the optimization strategy."""
        return f"Normalized optimizer with strategy presets (strategy: {self._strategy})"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs
    ) -> OptimizationResult:
        """Run normalized optimization."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        # Allow strategy override in kwargs
        strategy = kwargs.get('strategy', self._strategy)
        custom_weights = kwargs.get('custom_weights', self._custom_weights)

        if self.config.verbose:
            logger.info(f"Running normalized optimization: {len(candidates)} candidates")
            logger.info(f"Strategy: {strategy}")

        try:
            result = self._normalized.optimize(
                candidates=candidates,
                target_size=target,
                strategy=strategy,
                custom_weights=custom_weights,
                fixed_primers=fixed_primers,
                verbose=self.config.verbose,
            )

            metrics = self.compute_metrics(result.primers)

            return OptimizationResult(
                primers=tuple(result.primers),
                score=result.composite_score,
                status=OptimizationStatus.SUCCESS if result.primers else OptimizationStatus.NO_CONVERGENCE,
                metrics=metrics,
                iterations=len(result.stage_history),
                optimizer_name=self.name,
                message=(
                    f"Strategy: {result.strategy}, "
                    f"Coverage: {result.scores.coverage:.1%}, "
                    f"Connectivity: {result.scores.connectivity:.2f}, "
                    f"Amplification: {result.scores.amplification:.2f}, "
                    f"Background: {result.scores.background:.2f}"
                ),
            )

        except Exception as e:
            logger.error(f"Normalized optimization failed: {e}")
            return OptimizationResult.failure(self.name, str(e))
