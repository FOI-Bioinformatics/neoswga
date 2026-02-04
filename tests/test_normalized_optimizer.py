"""Tests for the normalized optimizer with strategy presets."""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from neoswga.core.normalized_optimizer import (
    NormalizedOptimizer,
    NormalizedScorer,
    NormalizedScores,
    NormalizedOptimizationResult,
    STRATEGY_PRESETS,
)


class MockPositionCache:
    """Mock position cache for testing."""

    def __init__(self, positions_dict=None):
        """
        Initialize mock cache.

        Args:
            positions_dict: Dict mapping (prefix, primer, strand) to positions
        """
        self.positions_dict = positions_dict or {}

    def get_positions(self, prefix, primer, strand='both'):
        key = (prefix, primer, strand)
        if key in self.positions_dict:
            return np.array(self.positions_dict[key])
        # Default: return some positions for testing
        if strand in ('forward', '+'):
            return np.array([1000, 5000, 15000, 30000])
        elif strand in ('reverse', '-'):
            return np.array([2000, 8000, 20000, 40000])
        return np.array([])


class TestNormalizedScores:
    """Tests for NormalizedScores dataclass."""

    def test_composite_equal_weights(self):
        """Test composite score with equal weights."""
        scores = NormalizedScores(
            coverage=0.8,
            connectivity=0.6,
            amplification=0.7,
            background=0.9,
        )
        weights = {
            'coverage': 0.25,
            'connectivity': 0.25,
            'amplification': 0.25,
            'background': 0.25,
        }
        composite = scores.composite(weights)
        expected = 0.25 * (0.8 + 0.6 + 0.7 + 0.9)
        assert abs(composite - expected) < 0.001

    def test_composite_clinical_weights(self):
        """Test composite score with clinical (background-heavy) weights."""
        scores = NormalizedScores(
            coverage=0.8,
            connectivity=0.6,
            amplification=0.7,
            background=0.9,  # High background score = low background binding
        )
        weights = STRATEGY_PRESETS['clinical']['weights']
        composite = scores.composite(weights)

        # Clinical weights: background=0.4, others=0.2 each
        expected = 0.2 * 0.8 + 0.2 * 0.6 + 0.2 * 0.7 + 0.4 * 0.9
        assert abs(composite - expected) < 0.001

    def test_composite_discovery_weights(self):
        """Test composite score with discovery (coverage-heavy) weights."""
        scores = NormalizedScores(
            coverage=0.9,  # High coverage
            connectivity=0.6,
            amplification=0.5,
            background=0.3,  # Lower background score (more background binding)
        )
        weights = STRATEGY_PRESETS['discovery']['weights']
        composite = scores.composite(weights)

        # Discovery weights: coverage=0.4, connectivity=0.3, amplification=0.2, background=0.1
        expected = 0.4 * 0.9 + 0.3 * 0.6 + 0.2 * 0.5 + 0.1 * 0.3
        assert abs(composite - expected) < 0.001

    def test_to_dict(self):
        """Test conversion to dictionary."""
        scores = NormalizedScores(
            coverage=0.8,
            connectivity=0.6,
            amplification=0.7,
            background=0.9,
            raw_coverage=0.8,
            raw_connectivity=2.5,
            raw_amplification=1000.0,
            raw_background=5.0,
        )
        d = scores.to_dict()
        assert d['coverage'] == 0.8
        assert d['raw_amplification'] == 1000.0


class TestNormalizedScorer:
    """Tests for NormalizedScorer class."""

    def test_scorer_initialization(self):
        """Test scorer initializes correctly."""
        cache = MockPositionCache()
        scorer = NormalizedScorer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
            bg_prefixes=['test_bg'],
            bg_seq_lengths=[200000],
        )
        assert scorer.total_fg_length == 100000
        assert scorer.total_bg_length == 200000

    def test_score_empty_primers(self):
        """Test scoring empty primer list."""
        cache = MockPositionCache()
        scorer = NormalizedScorer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )
        scores = scorer.score([])
        assert scores.coverage == 0.0
        assert scores.connectivity == 0.0
        assert scores.amplification == 0.0
        assert scores.background == 1.0  # No background = perfect

    def test_normalize_amplification(self):
        """Test amplification normalization."""
        cache = MockPositionCache()
        scorer = NormalizedScorer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        # Test edge cases
        assert scorer._normalize_amplification(1.0) == 0.0  # Min
        assert scorer._normalize_amplification(1e6) == 1.0  # Max
        assert 0.0 < scorer._normalize_amplification(1000) < 1.0  # Middle

    def test_normalize_background(self):
        """Test background normalization (inverted)."""
        cache = MockPositionCache()
        scorer = NormalizedScorer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        # No background = perfect score
        assert scorer._normalize_background(1.0) == 1.0

        # High background = low score
        bg_score = scorer._normalize_background(1e6)
        assert bg_score == 0.0

        # Middle values
        mid_score = scorer._normalize_background(1000)
        assert 0.0 < mid_score < 1.0


class TestNormalizedOptimizer:
    """Tests for NormalizedOptimizer class."""

    def test_optimizer_initialization(self):
        """Test optimizer initializes correctly."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )
        assert optimizer.scorer is not None

    def test_list_strategies(self):
        """Test listing available strategies."""
        strategies = NormalizedOptimizer.list_strategies()
        assert 'clinical' in strategies
        assert 'discovery' in strategies
        assert 'balanced' in strategies
        assert 'fast' in strategies
        assert 'enrichment' in strategies

        # Each strategy should have description and weights
        for name, info in strategies.items():
            assert 'description' in info
            assert 'weights' in info
            assert 'recommended_for' in info

    def test_recommend_strategy_clinical(self):
        """Test strategy recommendation for clinical use."""
        strategy = NormalizedOptimizer.recommend_strategy(application='clinical')
        assert strategy == 'clinical'

    def test_recommend_strategy_discovery(self):
        """Test strategy recommendation for discovery."""
        strategy = NormalizedOptimizer.recommend_strategy(application='discovery')
        assert strategy == 'discovery'

    def test_recommend_strategy_time_sensitive(self):
        """Test time-sensitive recommendation."""
        strategy = NormalizedOptimizer.recommend_strategy(time_sensitive=True)
        assert strategy == 'fast'

    def test_recommend_strategy_no_background(self):
        """Test recommendation when no background available."""
        strategy = NormalizedOptimizer.recommend_strategy(background_available=False)
        assert strategy == 'fast'

    def test_optimize_basic(self):
        """Test basic optimization."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        candidates = ['ATCGATCG', 'GCTAGCTA', 'TTAATTAA', 'CCGGCCGG']
        result = optimizer.optimize(
            candidates=candidates,
            target_size=2,
            strategy='balanced',
            verbose=False,
        )

        assert isinstance(result, NormalizedOptimizationResult)
        assert len(result.primers) == 2
        assert 0.0 <= result.composite_score <= 1.0
        assert result.strategy == 'balanced'

    def test_optimize_with_fixed_primers(self):
        """Test optimization with fixed primers."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        candidates = ['ATCGATCG', 'GCTAGCTA', 'TTAATTAA', 'CCGGCCGG']
        result = optimizer.optimize(
            candidates=candidates,
            target_size=3,
            strategy='balanced',
            fixed_primers=['ATCGATCG'],
            verbose=False,
        )

        assert 'ATCGATCG' in result.primers
        assert len(result.primers) == 3

    def test_optimize_custom_weights(self):
        """Test optimization with custom weights."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        custom_weights = {
            'coverage': 0.5,
            'connectivity': 0.3,
            'amplification': 0.1,
            'background': 0.1,
        }

        candidates = ['ATCGATCG', 'GCTAGCTA', 'TTAATTAA']
        result = optimizer.optimize(
            candidates=candidates,
            target_size=2,
            custom_weights=custom_weights,
            verbose=False,
        )

        assert result.strategy == 'custom'
        assert result.weights == custom_weights

    def test_optimize_invalid_strategy(self):
        """Test that invalid strategy raises error."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        with pytest.raises(ValueError, match="Unknown strategy"):
            optimizer.optimize(
                candidates=['ATCGATCG'],
                target_size=1,
                strategy='invalid_strategy',
                verbose=False,
            )


class TestStrategyPresets:
    """Tests for strategy preset configurations."""

    def test_presets_have_required_keys(self):
        """Test all presets have required keys."""
        for name, preset in STRATEGY_PRESETS.items():
            assert 'description' in preset, f"{name} missing description"
            assert 'weights' in preset, f"{name} missing weights"
            assert 'recommended_for' in preset, f"{name} missing recommended_for"

    def test_weights_sum_to_one(self):
        """Test all preset weights sum to approximately 1.0."""
        for name, preset in STRATEGY_PRESETS.items():
            weights = preset['weights']
            total = sum(weights.values())
            assert abs(total - 1.0) < 0.001, f"{name} weights sum to {total}"

    def test_weights_have_all_components(self):
        """Test all presets have all four weight components."""
        required_keys = {'coverage', 'connectivity', 'amplification', 'background'}
        for name, preset in STRATEGY_PRESETS.items():
            weights = preset['weights']
            assert set(weights.keys()) == required_keys, f"{name} missing weight keys"

    def test_weights_are_valid(self):
        """Test all weights are in valid range [0, 1]."""
        for name, preset in STRATEGY_PRESETS.items():
            for key, value in preset['weights'].items():
                assert 0.0 <= value <= 1.0, f"{name}.{key} = {value} out of range"


class TestNormalizedOptimizationResult:
    """Tests for NormalizedOptimizationResult dataclass."""

    def test_summary_generation(self):
        """Test summary generates valid output."""
        scores = NormalizedScores(
            coverage=0.8,
            connectivity=0.6,
            amplification=0.7,
            background=0.9,
            raw_coverage=0.8,
            raw_connectivity=2.5,
            raw_amplification=1000.0,
            raw_background=5.0,
        )
        result = NormalizedOptimizationResult(
            primers=['ATCGATCG', 'GCTAGCTA'],
            scores=scores,
            composite_score=0.75,
            strategy='clinical',
            weights=STRATEGY_PRESETS['clinical']['weights'],
            runtime_seconds=1.5,
            candidates_evaluated=100,
        )

        summary = result.summary()
        assert 'NORMALIZED OPTIMIZATION RESULT' in summary
        assert 'clinical' in summary
        assert '2' in summary  # Number of primers
        assert '0.75' in summary  # Composite score


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_single_candidate(self):
        """Test optimization with single candidate."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        result = optimizer.optimize(
            candidates=['ATCGATCG'],
            target_size=1,
            strategy='balanced',
            verbose=False,
        )

        assert len(result.primers) == 1
        assert result.primers[0] == 'ATCGATCG'

    def test_target_exceeds_candidates(self):
        """Test when target size exceeds available candidates."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
        )

        candidates = ['ATCGATCG', 'GCTAGCTA']
        result = optimizer.optimize(
            candidates=candidates,
            target_size=5,  # More than available
            strategy='balanced',
            verbose=False,
        )

        # Should select all available
        assert len(result.primers) == 2

    def test_no_background_genome(self):
        """Test optimization without background genome."""
        cache = MockPositionCache()
        optimizer = NormalizedOptimizer(
            position_cache=cache,
            fg_prefixes=['test_fg'],
            fg_seq_lengths=[100000],
            bg_prefixes=None,  # No background
            bg_seq_lengths=None,
        )

        result = optimizer.optimize(
            candidates=['ATCGATCG', 'GCTAGCTA'],
            target_size=2,
            strategy='balanced',
            verbose=False,
        )

        # Should still work, background score should be 1.0 (perfect)
        assert len(result.primers) == 2
        assert result.scores.background == 1.0
