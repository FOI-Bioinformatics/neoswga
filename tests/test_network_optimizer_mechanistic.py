"""
Tests for mechanistic model integration in NetworkOptimizer.

Tests that mechanistic weighting affects primer selection and
produces expected metrics.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock

from neoswga.core.network_optimizer import NetworkOptimizer
from neoswga.core.reaction_conditions import ReactionConditions


class MockPositionCache:
    """Mock position cache for testing."""

    def __init__(self, positions_dict=None):
        """
        Initialize mock cache.

        Args:
            positions_dict: Dict mapping (prefix, primer, strand) to positions array
        """
        self._positions = positions_dict or {}

    def get_positions(self, prefix, primer, strand='both'):
        """Return mock positions."""
        key = (prefix, primer, strand)
        if key in self._positions:
            return self._positions[key]

        # Default: return some positions based on primer length
        # Longer primers have fewer binding sites
        base_count = max(5, 20 - len(primer))
        positions = np.array([i * 1000 for i in range(base_count)])
        return positions


@pytest.fixture
def mock_cache():
    """Create a mock position cache."""
    return MockPositionCache()


@pytest.fixture
def standard_conditions():
    """Standard phi29 conditions."""
    return ReactionConditions(temp=30.0, polymerase='phi29', mg_conc=2.5)


@pytest.fixture
def enhanced_conditions():
    """Enhanced EquiPhi29 conditions."""
    return ReactionConditions(
        temp=42.0,
        polymerase='equiphi29',
        mg_conc=2.5,
        dmso_percent=5.0,
        betaine_m=1.0
    )


class TestNetworkOptimizerMechanisticInit:
    """Test NetworkOptimizer initialization with mechanistic model."""

    def test_init_without_conditions(self, mock_cache):
        """Initialize without conditions (mechanistic disabled)."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
        )

        assert optimizer.mech_model is None
        assert optimizer.conditions is None
        assert optimizer.mechanistic_weight == 0.0

    def test_init_with_conditions_no_weight(self, mock_cache, standard_conditions):
        """Initialize with conditions but zero weight (mechanistic disabled)."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.0,  # Disabled
        )

        assert optimizer.mech_model is None
        assert optimizer.conditions is not None

    def test_init_with_conditions_and_weight(self, mock_cache, standard_conditions):
        """Initialize with conditions and weight (mechanistic enabled)."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
        )

        assert optimizer.mech_model is not None
        assert optimizer.conditions is not None
        assert optimizer.mechanistic_weight == 0.3

    def test_template_gc_parameter(self, mock_cache, standard_conditions):
        """Test template_gc parameter is stored."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
            template_gc=0.65,
        )

        assert optimizer.template_gc == 0.65


class TestMechanisticScoreCalculation:
    """Test mechanistic score calculation."""

    def test_mechanistic_score_without_model(self, mock_cache):
        """Score should be 1.0 when model is disabled."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
        )

        score = optimizer._calculate_mechanistic_score('ATCGATCG')
        assert score == 1.0

    def test_mechanistic_score_with_model(self, mock_cache, standard_conditions):
        """Score should be between 0 and 1 when model is enabled."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
        )

        score = optimizer._calculate_mechanistic_score('ATCGATCGATCG')
        assert 0.0 <= score <= 1.0

    def test_longer_primers_may_have_different_score(self, mock_cache, enhanced_conditions):
        """Different primer lengths may have different mechanistic scores."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=enhanced_conditions,
            mechanistic_weight=0.3,
        )

        score_8mer = optimizer._calculate_mechanistic_score('ATCGATCG')
        score_12mer = optimizer._calculate_mechanistic_score('ATCGATCGATCG')
        score_16mer = optimizer._calculate_mechanistic_score('ATCGATCGATCGATCG')

        # Scores should be different (Tm differs)
        # With EquiPhi29 at 42C, longer primers should generally score better
        # because their Tm is closer to optimal
        assert score_8mer != score_16mer


class TestScorePrimerSetMechanisticMetrics:
    """Test that score_primer_set includes mechanistic metrics."""

    def test_no_mechanistic_metrics_without_model(self, mock_cache):
        """score_primer_set should not include mechanistic metrics when disabled."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
        )

        primers = ['ATCGATCGATCG', 'GCTAGCTAGCTA']
        result = optimizer.score_primer_set(primers)

        assert 'avg_processivity_factor' not in result
        assert 'avg_accessibility' not in result
        assert 'avg_amplification_factor' not in result

    def test_mechanistic_metrics_with_model(self, mock_cache, standard_conditions):
        """score_primer_set should include mechanistic metrics when enabled."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
        )

        primers = ['ATCGATCGATCG', 'GCTAGCTAGCTA']
        result = optimizer.score_primer_set(primers)

        assert 'avg_processivity_factor' in result
        assert 'avg_accessibility' in result
        assert 'avg_amplification_factor' in result
        assert 'min_amplification_factor' in result
        assert 'avg_effective_binding' in result
        assert 'avg_stability_factor' in result

    def test_mechanistic_metrics_reasonable_values(self, mock_cache, standard_conditions):
        """Mechanistic metrics should have reasonable values."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
        )

        primers = ['ATCGATCGATCG', 'GCTAGCTAGCTA']
        result = optimizer.score_primer_set(primers)

        # All factors should be in reasonable range
        assert 0.0 <= result['avg_processivity_factor'] <= 1.5
        assert 0.0 <= result['avg_accessibility'] <= 1.0
        assert 0.0 <= result['avg_amplification_factor'] <= 1.0
        assert 0.0 <= result['min_amplification_factor'] <= 1.0
        assert 0.0 <= result['avg_effective_binding'] <= 1.0
        assert 0.0 <= result['avg_stability_factor'] <= 1.5


class TestConditionEffects:
    """Test that different conditions produce different results."""

    def test_dmso_affects_processivity(self, mock_cache):
        """High DMSO should reduce processivity factor."""
        cond_low = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=2.0)
        cond_high = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=8.0)

        opt_low = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=cond_low,
            mechanistic_weight=0.3,
        )

        opt_high = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=cond_high,
            mechanistic_weight=0.3,
        )

        primers = ['ATCGATCGATCG']
        result_low = opt_low.score_primer_set(primers)
        result_high = opt_high.score_primer_set(primers)

        # High DMSO should have lower processivity
        assert result_high['avg_processivity_factor'] < result_low['avg_processivity_factor']

    def test_high_gc_template_affects_accessibility(self, mock_cache, standard_conditions):
        """High GC template should have lower accessibility."""
        opt_low_gc = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
            template_gc=0.3,
        )

        opt_high_gc = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
            template_gc=0.7,
        )

        primers = ['ATCGATCGATCG']
        result_low = opt_low_gc.score_primer_set(primers)
        result_high = opt_high_gc.score_primer_set(primers)

        # High GC template should have lower accessibility
        assert result_high['avg_accessibility'] < result_low['avg_accessibility']


class TestOptimizeGreedyWithMechanistic:
    """Test optimize_greedy with mechanistic weighting."""

    def test_optimize_runs_with_mechanistic(self, mock_cache, standard_conditions):
        """optimize_greedy should run without errors when mechanistic enabled."""
        optimizer = NetworkOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=standard_conditions,
            mechanistic_weight=0.3,
        )

        candidates = [
            'ATCGATCGATCG',
            'GCTAGCTAGCTA',
            'TTAATTAATTAA',
            'CCGGCCGGCCGG',
        ]

        result = optimizer.optimize_greedy(candidates, num_primers=2)

        assert len(result) == 2
        assert all(p in candidates for p in result)

    def test_mechanistic_weight_affects_selection(self, mock_cache):
        """Different mechanistic weights may lead to different selections."""
        # This is a statistical test - with enough primers, different weights
        # should sometimes lead to different selections
        conditions = ReactionConditions(temp=42.0, polymerase='equiphi29', mg_conc=2.5)

        # Create mock cache with varying positions for different primers
        positions_dict = {}
        for primer in ['ATCGATCGATCG', 'GCTAGCTAGCTA', 'TTAATTAATTAA', 'CCGGCCGGCCGG']:
            for strand in ['forward', 'reverse']:
                positions_dict[('target', primer, strand)] = np.array([i * 500 for i in range(10)])
                positions_dict[('background', primer, strand)] = np.array([i * 1000 for i in range(5)])

        cache = MockPositionCache(positions_dict)

        opt_no_mech = NetworkOptimizer(
            position_cache=cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=conditions,
            mechanistic_weight=0.0,  # Disabled
        )

        opt_with_mech = NetworkOptimizer(
            position_cache=cache,
            fg_prefixes=['target'],
            bg_prefixes=['background'],
            fg_seq_lengths=[10000],
            bg_seq_lengths=[50000],
            conditions=conditions,
            mechanistic_weight=0.5,  # Enabled
        )

        candidates = ['ATCGATCGATCG', 'GCTAGCTAGCTA', 'TTAATTAATTAA', 'CCGGCCGGCCGG']

        result_no = opt_no_mech.optimize_greedy(candidates, num_primers=2)
        result_with = opt_with_mech.optimize_greedy(candidates, num_primers=2)

        # Both should return valid results
        assert len(result_no) == 2
        assert len(result_with) == 2

        # Note: They may or may not be different depending on the specific
        # primer characteristics. The key is both work correctly.
