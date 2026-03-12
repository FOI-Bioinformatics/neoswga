"""
Unit tests for refactored modules.

Tests the new infrastructure:
- base_optimizer.py: BaseOptimizer interface and result types
- exceptions.py: Custom exception hierarchy
- optimizer_factory.py: Optimizer registration and creation
- search_context.py: Typed search contexts
- additives.py: Additive concentrations
"""

import pytest
import numpy as np
from typing import List, Optional

# Import modules under test
from neoswga.core.base_optimizer import (
    BaseOptimizer, OptimizationResult, OptimizationStatus,
    PrimerSetMetrics, OptimizerConfig, CompositeOptimizer
)
from neoswga.core.exceptions import (
    NeoSWGAError, InvalidPrimerError, OptimizerNotFoundError,
    NoCandidatesError, ConfigurationError, InvalidParameterError
)
from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry
from neoswga.core.search_context import (
    GenomeInfo, PositionData, DimerConstraints, BFSConfig,
    BFSSearchContext, SearchState, SearchResult
)
from neoswga.core.additives import AdditiveConcentrations, estimate_combined_effect


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def sample_primers():
    """Sample primer sequences for testing."""
    return ['ATCGATCG', 'GCTAGCTA', 'AATTCCGG', 'TTAACCGG', 'GGCCAATT']


@pytest.fixture
def sample_genome_info():
    """Sample genome information."""
    return GenomeInfo.from_lists(
        prefixes=['data/target'],
        seq_lengths=[1000000],
        circular=True
    )


@pytest.fixture
def mock_position_cache():
    """Mock position cache for testing."""
    class MockCache:
        def __init__(self):
            self.positions = {}

        def get_positions(self, prefix, primer, strand='both'):
            key = (prefix, primer, strand)
            if key in self.positions:
                return self.positions[key]
            # Generate random positions
            return np.random.randint(0, 1000000, size=100)

    return MockCache()


# =============================================================================
# BaseOptimizer Tests
# =============================================================================

class TestOptimizationResult:
    """Tests for OptimizationResult dataclass."""

    def test_creation(self):
        """Test basic creation."""
        metrics = PrimerSetMetrics(
            fg_coverage=0.95,
            bg_coverage=0.05,
            coverage_uniformity=0.2,
            total_fg_sites=1000,
            total_bg_sites=50,
            selectivity_ratio=20.0,
            mean_tm=35.0,
            tm_range=(30.0, 40.0),
            dimer_risk_score=0.1,
            mean_gap=1000.0,
            max_gap=5000.0,
            gap_gini=0.3,
        )

        result = OptimizationResult(
            primers=('ATCGATCG', 'GCTAGCTA'),
            score=0.95,
            status=OptimizationStatus.SUCCESS,
            metrics=metrics,
            iterations=10,
            optimizer_name='test',
        )

        assert result.num_primers == 2
        assert result.is_success
        assert result.score == 0.95

    def test_failure_result(self):
        """Test failure result creation."""
        result = OptimizationResult.failure('test', 'Something went wrong')

        assert not result.is_success
        assert result.status == OptimizationStatus.ERROR
        assert 'Something went wrong' in result.message
        assert result.num_primers == 0

    def test_to_dict(self):
        """Test serialization to dict."""
        metrics = PrimerSetMetrics.empty()
        result = OptimizationResult(
            primers=('ATCGATCG',),
            score=0.5,
            status=OptimizationStatus.PARTIAL,
            metrics=metrics,
            iterations=5,
            optimizer_name='test',
        )

        d = result.to_dict()
        assert 'primers' in d
        assert 'score' in d
        assert 'metrics' in d
        assert d['num_primers'] == 1


class TestOptimizerConfig:
    """Tests for OptimizerConfig."""

    def test_default_values(self):
        """Test default configuration values."""
        config = OptimizerConfig()

        assert config.target_set_size == 6
        assert config.max_iterations == 100
        assert config.verbose is True

    def test_validation(self):
        """Test parameter validation."""
        with pytest.raises(ValueError):
            config = OptimizerConfig(target_set_size=0)
            config.validate()

        with pytest.raises(ValueError):
            config = OptimizerConfig(min_tm=50.0, max_tm=30.0)
            config.validate()


# =============================================================================
# Exception Tests
# =============================================================================

class TestExceptions:
    """Tests for custom exception hierarchy."""

    def test_base_exception(self):
        """Test NeoSWGAError is base of all custom exceptions."""
        assert issubclass(InvalidPrimerError, NeoSWGAError)
        assert issubclass(OptimizerNotFoundError, NeoSWGAError)
        assert issubclass(ConfigurationError, NeoSWGAError)

    def test_invalid_primer_error(self):
        """Test InvalidPrimerError contains primer info."""
        error = InvalidPrimerError('ATXG', 'Contains invalid character X')

        assert error.primer == 'ATXG'
        assert 'invalid character' in str(error)

    def test_optimizer_not_found_error(self):
        """Test OptimizerNotFoundError lists available optimizers."""
        error = OptimizerNotFoundError('unknown', ['greedy', 'genetic'])

        assert 'unknown' in str(error)
        assert 'greedy' in str(error)

    def test_no_candidates_error(self):
        """Test NoCandidatesError with filter stage."""
        error = NoCandidatesError(100, 'Tm filtering')

        assert error.original_count == 100
        assert 'Tm filtering' in str(error)


# =============================================================================
# Optimizer Factory Tests
# =============================================================================

class TestOptimizerRegistry:
    """Tests for OptimizerRegistry."""

    def setup_method(self):
        """Save registry state before each test."""
        # Save current state
        self._saved_registry = dict(OptimizerRegistry._registry)
        self._saved_descriptions = dict(OptimizerRegistry._descriptions)
        self._saved_aliases = dict(OptimizerRegistry._aliases)
        # Clear for isolated testing
        OptimizerRegistry.clear()

    def teardown_method(self):
        """Restore registry state after each test."""
        OptimizerRegistry.clear()
        OptimizerRegistry._registry.update(self._saved_registry)
        OptimizerRegistry._descriptions.update(self._saved_descriptions)
        OptimizerRegistry._aliases.update(self._saved_aliases)

    def test_register_decorator(self):
        """Test optimizer registration via decorator."""
        @OptimizerRegistry.register('test-optimizer')
        class TestOptimizer(BaseOptimizer):
            @property
            def name(self):
                return 'test'

            def optimize(self, candidates, target_size=None, **kwargs):
                pass

        assert OptimizerRegistry.is_registered('test-optimizer')

    def test_register_with_aliases(self):
        """Test registration with aliases."""
        @OptimizerRegistry.register('full-name', aliases=['short', 'alias'])
        class AliasedOptimizer(BaseOptimizer):
            @property
            def name(self):
                return 'full-name'

            def optimize(self, candidates, target_size=None, **kwargs):
                pass

        assert OptimizerRegistry.is_registered('full-name')
        assert OptimizerRegistry.is_registered('short')
        assert OptimizerRegistry.is_registered('alias')

    def test_get_unregistered_raises(self):
        """Test getting unregistered optimizer raises error."""
        with pytest.raises(OptimizerNotFoundError):
            OptimizerRegistry.get('nonexistent')

    def test_list_all(self):
        """Test listing all registered optimizers."""
        @OptimizerRegistry.register('opt1')
        class Opt1(BaseOptimizer):
            """First optimizer."""
            @property
            def name(self):
                return 'opt1'
            def optimize(self, candidates, target_size=None, **kwargs):
                pass

        @OptimizerRegistry.register('opt2')
        class Opt2(BaseOptimizer):
            """Second optimizer."""
            @property
            def name(self):
                return 'opt2'
            def optimize(self, candidates, target_size=None, **kwargs):
                pass

        all_opts = OptimizerRegistry.list_all()
        assert 'opt1' in all_opts
        assert 'opt2' in all_opts


# =============================================================================
# Search Context Tests
# =============================================================================

class TestGenomeInfo:
    """Tests for GenomeInfo dataclass."""

    def test_creation(self):
        """Test basic creation."""
        info = GenomeInfo(
            prefixes=('data/genome1', 'data/genome2'),
            seq_lengths=(1000000, 500000),
            circular=True
        )

        assert info.total_length == 1500000
        assert info.num_sequences == 2

    def test_from_lists(self):
        """Test creation from lists."""
        info = GenomeInfo.from_lists(
            prefixes=['a', 'b'],
            seq_lengths=[100, 200]
        )

        assert isinstance(info.prefixes, tuple)
        assert info.total_length == 300

    def test_validation(self):
        """Test mismatched lengths raise error."""
        with pytest.raises(ValueError):
            GenomeInfo(
                prefixes=('a', 'b', 'c'),
                seq_lengths=(100, 200)  # Mismatched!
            )


class TestPositionData:
    """Tests for PositionData class."""

    def test_empty_creation(self):
        """Test creating empty position data."""
        data = PositionData.empty(3)

        assert len(data.positions) == 3
        assert data.total_sites() == 0

    def test_add_positions(self):
        """Test adding positions."""
        data = PositionData.empty(2)

        data.add_positions(0, np.array([100, 200]), np.array([150, 250]))
        data.add_positions(1, np.array([300]), np.array([350]))

        assert data.total_sites() == 6

    def test_copy(self):
        """Test deep copy."""
        data = PositionData.empty(1)
        data.add_positions(0, np.array([100]), np.array([200]))

        copy = data.copy()
        copy.add_positions(0, np.array([300]), np.array([400]))

        # Original should be unchanged
        assert data.total_sites() == 2
        assert copy.total_sites() == 4


class TestBFSConfig:
    """Tests for BFSConfig."""

    def test_defaults(self):
        """Test default values."""
        config = BFSConfig()

        assert config.max_sets == 10
        assert config.iterations == 10
        assert config.verbose is True

    def test_validation(self):
        """Test invalid config raises error."""
        with pytest.raises(ValueError):
            BFSConfig(max_sets=0)

        with pytest.raises(ValueError):
            BFSConfig(selection_method='invalid')


# =============================================================================
# Additives Tests
# =============================================================================

class TestAdditiveConcentrations:
    """Tests for AdditiveConcentrations."""

    def test_default_is_empty(self):
        """Test default has no additives."""
        additives = AdditiveConcentrations()

        assert additives.dmso_percent == 0.0
        assert additives.betaine_m == 0.0

    def test_validation(self):
        """Test out-of-range values raise error."""
        with pytest.raises(ValueError):
            AdditiveConcentrations(dmso_percent=15.0)  # Max is 10%

        with pytest.raises(ValueError):
            AdditiveConcentrations(betaine_m=-0.5)  # Min is 0

    def test_tm_correction_dmso(self):
        """Test DMSO Tm correction."""
        additives = AdditiveConcentrations(dmso_percent=5.0)
        correction = additives.calculate_tm_correction()

        # DMSO should lower Tm (~0.6C per %)
        assert correction < 0
        assert abs(correction - (-3.0)) < 0.5

    def test_tm_correction_betaine(self):
        """Test betaine Tm correction."""
        additives = AdditiveConcentrations(betaine_m=1.0)

        # Betaine should lower Tm uniformly
        correction = additives.calculate_tm_correction()
        assert correction < 0

    def test_gc_normalization(self):
        """Test GC-dependent corrections."""
        additives = AdditiveConcentrations(betaine_m=2.0)

        # High GC should get negative correction
        gc_rich = additives.calculate_tm_correction(gc_content=0.7, primer_length=10)

        # Low GC should get positive correction
        at_rich = additives.calculate_tm_correction(gc_content=0.3, primer_length=10)

        assert gc_rich < at_rich

    def test_factory_methods(self):
        """Test factory method creation."""
        standard = AdditiveConcentrations.for_standard_phi29()
        assert standard.dmso_percent == 0.0

        enhanced = AdditiveConcentrations.for_enhanced_equiphi29()
        assert enhanced.dmso_percent == 5.0
        assert enhanced.betaine_m == 1.0

        extreme = AdditiveConcentrations.for_extreme_gc()
        assert extreme.tmac_m > 0

    def test_max_primer_length(self):
        """Test max primer length calculation."""
        standard = AdditiveConcentrations.for_standard_phi29()
        assert standard.max_supported_primer_length() == 12

        enhanced = AdditiveConcentrations.for_enhanced_equiphi29()
        assert enhanced.max_supported_primer_length() > 12

        long_primers = AdditiveConcentrations.for_long_primers()
        assert long_primers.max_supported_primer_length() >= 17

    def test_gc_range(self):
        """Test GC content range recommendations."""
        standard = AdditiveConcentrations.for_standard_phi29()
        gc_range = standard.gc_content_range()
        assert gc_range[1] - gc_range[0] < 0.3  # Narrow range

        high_gc = AdditiveConcentrations.for_high_gc()
        gc_range = high_gc.gc_content_range()
        assert gc_range[1] - gc_range[0] >= 0.3  # Wider range

    def test_to_dict(self):
        """Test serialization."""
        additives = AdditiveConcentrations(dmso_percent=5.0, betaine_m=1.0)
        d = additives.to_dict()

        assert d['dmso_percent'] == 5.0
        assert d['betaine_m'] == 1.0

    def test_from_dict(self):
        """Test deserialization."""
        d = {'dmso_percent': 3.0, 'trehalose_m': 0.5}
        additives = AdditiveConcentrations.from_dict(d)

        assert additives.dmso_percent == 3.0
        assert additives.trehalose_m == 0.5

    def test_repr(self):
        """Test string representation."""
        additives = AdditiveConcentrations(dmso_percent=5.0, betaine_m=1.0)
        s = repr(additives)

        assert 'DMSO=5.0%' in s
        assert 'betaine=1.0M' in s

        empty = AdditiveConcentrations()
        assert 'none' in repr(empty)


class TestEstimateCombinedEffect:
    """Tests for estimate_combined_effect function."""

    def test_returns_all_fields(self):
        """Test all expected fields are returned."""
        additives = AdditiveConcentrations.for_enhanced_equiphi29()
        effects = estimate_combined_effect(additives, gc_content=0.5, primer_length=12)

        assert 'tm_correction' in effects
        assert 'max_primer_length' in effects
        assert 'gc_range' in effects
        assert 'polymerase_activity_factor' in effects
        assert 'secondary_structure_reduction' in effects
        assert 'gc_normalization' in effects

    def test_activity_factor_range(self):
        """Test activity factor is in reasonable range."""
        additives = AdditiveConcentrations.for_enhanced_equiphi29()
        effects = estimate_combined_effect(additives, gc_content=0.5, primer_length=12)

        assert 0.5 <= effects['polymerase_activity_factor'] <= 1.5


# =============================================================================
# Integration Tests
# =============================================================================

class TestReactionConditionsIntegration:
    """Test integration with ReactionConditions."""

    def test_additives_property(self):
        """Test ReactionConditions.additives property."""
        from neoswga.core.reaction_conditions import ReactionConditions

        conditions = ReactionConditions(
            temp=42.0,
            dmso_percent=5.0,
            betaine_m=1.0,
            polymerase='equiphi29'
        )

        additives = conditions.additives

        assert isinstance(additives, AdditiveConcentrations)
        assert additives.dmso_percent == 5.0
        assert additives.betaine_m == 1.0

    def test_from_additives(self):
        """Test creating ReactionConditions from AdditiveConcentrations."""
        from neoswga.core.reaction_conditions import ReactionConditions

        additives = AdditiveConcentrations.for_extreme_gc()
        conditions = ReactionConditions.from_additives(
            additives,
            temp=42.0,
            polymerase='equiphi29'
        )

        assert conditions.temp == 42.0
        assert conditions.dmso_percent == additives.dmso_percent
        assert conditions.betaine_m == additives.betaine_m


# =============================================================================
# BoundedScoreCache Tests
# =============================================================================

class TestBoundedScoreCache:
    """Tests for BoundedScoreCache LRU implementation."""

    def test_basic_get_set(self):
        """Test basic get/set operations."""
        from neoswga.core.greedy_optimizer import BoundedScoreCache

        cache = BoundedScoreCache(max_size=100)
        key = frozenset(['ATCG', 'GCTA'])

        cache.set(key, 0.95)
        assert cache.get(key) == 0.95

    def test_eviction_on_max_size(self):
        """Test LRU eviction when max size is reached."""
        from neoswga.core.greedy_optimizer import BoundedScoreCache

        cache = BoundedScoreCache(max_size=3)

        cache.set(frozenset(['A']), 1.0)
        cache.set(frozenset(['B']), 2.0)
        cache.set(frozenset(['C']), 3.0)

        # Cache is full, adding D should evict A (oldest)
        cache.set(frozenset(['D']), 4.0)

        assert len(cache) == 3
        assert cache.get(frozenset(['A'])) is None  # Evicted
        assert cache.get(frozenset(['B'])) == 2.0
        assert cache.get(frozenset(['C'])) == 3.0
        assert cache.get(frozenset(['D'])) == 4.0

    def test_access_updates_lru_order(self):
        """Test that accessing an item moves it to end (most recent)."""
        from neoswga.core.greedy_optimizer import BoundedScoreCache

        cache = BoundedScoreCache(max_size=3)

        cache.set(frozenset(['A']), 1.0)
        cache.set(frozenset(['B']), 2.0)
        cache.set(frozenset(['C']), 3.0)

        # Access A - moves it to end
        cache.get(frozenset(['A']))

        # Add D - should evict B (now oldest)
        cache.set(frozenset(['D']), 4.0)

        assert cache.get(frozenset(['B'])) is None  # Evicted
        assert cache.get(frozenset(['A'])) == 1.0   # Still there

    def test_hit_rate_tracking(self):
        """Test hit/miss tracking."""
        from neoswga.core.greedy_optimizer import BoundedScoreCache

        cache = BoundedScoreCache(max_size=10)
        key = frozenset(['X'])

        cache.set(key, 1.0)

        # 2 hits
        cache.get(key)
        cache.get(key)

        # 1 miss
        cache.get(frozenset(['Y']))

        assert cache.hits == 2
        assert cache.misses == 1
        assert abs(cache.hit_rate() - 0.667) < 0.01

    def test_clear_stats(self):
        """Test clearing statistics."""
        from neoswga.core.greedy_optimizer import BoundedScoreCache

        cache = BoundedScoreCache(max_size=10)
        cache.set(frozenset(['X']), 1.0)
        cache.get(frozenset(['X']))
        cache.get(frozenset(['Y']))

        cache.clear_stats()

        assert cache.hits == 0
        assert cache.misses == 0
        assert cache.hit_rate() == 0.0


class TestMultiAgentOrchestrator:
    """Tests for MultiAgentOrchestrator thread pool reuse."""

    def test_executor_lazy_initialization(self, mock_position_cache):
        """Test that executor is not created until needed."""
        from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

        orchestrator = MultiAgentOrchestrator(
            position_cache=mock_position_cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
        )

        # Executor should not exist yet
        assert orchestrator._executor is None

        # Get executor - should create it
        executor = orchestrator._get_executor(4)
        assert executor is not None
        assert orchestrator._executor is not None
        assert orchestrator._executor._max_workers >= 4

        # Cleanup
        orchestrator.shutdown()

    def test_executor_reuse(self, mock_position_cache):
        """Test that executor is reused across calls."""
        from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

        orchestrator = MultiAgentOrchestrator(
            position_cache=mock_position_cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
        )

        # Get executor twice with same worker count
        executor1 = orchestrator._get_executor(4)
        executor2 = orchestrator._get_executor(4)

        # Should be the same instance
        assert executor1 is executor2

        # Cleanup
        orchestrator.shutdown()

    def test_executor_upgrade_on_larger_workers(self, mock_position_cache):
        """Test that executor is recreated if more workers needed."""
        from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

        orchestrator = MultiAgentOrchestrator(
            position_cache=mock_position_cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
        )

        # Get executor with 2 workers
        executor1 = orchestrator._get_executor(2)
        assert executor1._max_workers == 2

        # Get executor with 4 workers - should create new one
        executor2 = orchestrator._get_executor(4)
        assert executor2._max_workers == 4

        # Should not be the same instance (old one shutdown)
        assert executor1 is not executor2

        # Cleanup
        orchestrator.shutdown()

    def test_context_manager_support(self, mock_position_cache):
        """Test context manager properly shuts down executor."""
        from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

        with MultiAgentOrchestrator(
            position_cache=mock_position_cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
        ) as orchestrator:
            # Create executor
            executor = orchestrator._get_executor(4)
            assert executor is not None

        # After context exit, executor should be None
        assert orchestrator._executor is None

    def test_shutdown_idempotent(self, mock_position_cache):
        """Test that shutdown can be called multiple times safely."""
        from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

        orchestrator = MultiAgentOrchestrator(
            position_cache=mock_position_cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
        )

        # Create executor
        orchestrator._get_executor(4)

        # Shutdown multiple times - should not raise
        orchestrator.shutdown()
        orchestrator.shutdown()
        orchestrator.shutdown()

        assert orchestrator._executor is None

    def test_executor_with_no_workers(self, mock_position_cache):
        """Test executor creation with 0 workers falls back to default."""
        from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

        orchestrator = MultiAgentOrchestrator(
            position_cache=mock_position_cache,
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
        )

        # Getting executor with 0 workers should still work
        executor = orchestrator._get_executor(0)
        assert executor is not None

        # Cleanup
        orchestrator.shutdown()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
