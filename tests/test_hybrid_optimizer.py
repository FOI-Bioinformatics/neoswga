"""
Unit tests for hybrid_optimizer module.

Tests:
- HybridOptimizer initialization with different polymerase configs
- PolymeraseConfig presets
- Two-stage optimization workflow
- Stage 1 vs Stage 2 primer count relationship
- Edge cases (empty candidates, fewer candidates than requested)
- Background pruning option
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch, MagicMock
from typing import List

from neoswga.core.hybrid_optimizer import (
    HybridOptimizer,
    HybridResult,
    PolymeraseConfig,
    POLYMERASE_PRESETS,
    _get_polymerase_config,
)


# =============================================================================
# Helpers
# =============================================================================

def _make_mock_cache(positions_map=None):
    """Create a mock PositionCache with configurable position data.

    Args:
        positions_map: dict mapping (prefix, primer, strand) -> np.ndarray.
            If a key is missing, returns positions based on primer name
            using a deterministic default scheme.
    """
    cache = Mock()
    positions_map = positions_map or {}

    def get_positions(prefix, primer, strand='both'):
        key = (prefix, primer, strand)
        if key in positions_map:
            return positions_map[key]
        # Default: return deterministic positions based on primer hash
        if strand == 'both':
            fwd = get_positions(prefix, primer, 'forward')
            rev = get_positions(prefix, primer, 'reverse')
            return np.concatenate([fwd, rev])
        # Use hash to generate reproducible positions
        h = hash((primer, strand)) % 10000
        return np.array([h, h + 10000, h + 20000], dtype=np.int32)

    cache.get_positions = Mock(side_effect=get_positions)
    return cache


def _make_primers(n=20):
    """Generate n deterministic primer sequences."""
    bases = 'ACGT'
    primers = []
    for i in range(n):
        seq = ''
        val = i
        for _ in range(8):
            seq += bases[val % 4]
            val //= 4
        primers.append(seq)
    return primers


# =============================================================================
# PolymeraseConfig Tests
# =============================================================================

class TestPolymeraseConfig:
    """Tests for PolymeraseConfig and presets."""

    def test_phi29_preset_defaults(self):
        """Test that phi29 preset has expected defaults."""
        config = POLYMERASE_PRESETS['phi29']

        assert config.max_extension == 70000
        assert config.thermo_filter is False
        assert config.primer_multiplier == 1.0
        assert config.reaction_temp == 30.0

    def test_equiphi29_preset(self):
        """Test that equiphi29 preset has higher temperature and thermo filter."""
        config = POLYMERASE_PRESETS['equiphi29']

        assert config.max_extension == 80000
        assert config.thermo_filter is True
        assert config.primer_multiplier == 0.85
        assert config.reaction_temp == 42.0
        assert config.min_primer_tm == 37.0
        assert config.max_primer_tm == 62.0

    def test_bst_preset(self):
        """Test that bst preset has thermostable settings."""
        config = POLYMERASE_PRESETS['bst']

        assert config.max_extension == 10000
        assert config.thermo_filter is True
        assert config.reaction_temp == 60.0
        assert config.min_primer_tm == 50.0

    def test_klenow_preset(self):
        """Test that klenow preset has lower processivity."""
        config = POLYMERASE_PRESETS['klenow']

        assert config.max_extension == 5000
        assert config.thermo_filter is False
        assert config.reaction_temp == 37.0

    def test_all_presets_present(self):
        """Test that all expected polymerase presets are defined."""
        expected = {'phi29', 'equiphi29', 'bst', 'klenow'}
        assert set(POLYMERASE_PRESETS.keys()) == expected

    def test_get_polymerase_config_known(self):
        """Test that known polymerase names return correct config."""
        config = _get_polymerase_config('equiphi29')
        assert config.reaction_temp == 42.0

    def test_get_polymerase_config_unknown_falls_back(self):
        """Test that unknown polymerase falls back to phi29."""
        config = _get_polymerase_config('unknown_enzyme')
        assert config.max_extension == 70000
        assert config.reaction_temp == 30.0


# =============================================================================
# HybridOptimizer Initialization Tests
# =============================================================================

class TestHybridOptimizerInit:
    """Tests for HybridOptimizer initialization."""

    def test_init_default_phi29(self):
        """Test initialization with default phi29 configuration."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            bin_size=10000,
        )

        assert optimizer.polymerase == 'phi29'
        assert optimizer.max_extension == 70000
        assert optimizer.bin_size == 10000
        assert optimizer.fg_prefixes == ['genome1']
        assert optimizer.fg_seq_lengths == [100000]
        assert optimizer.bg_prefixes == []
        assert optimizer.bg_seq_lengths == []
        assert optimizer.background_pruning is False

    def test_init_equiphi29_with_gc_adjustment(self):
        """Test initialization with equiphi29 and GC-rich genome adjustment."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            polymerase='equiphi29',
            genome_gc_content=0.70,
        )

        assert optimizer.polymerase == 'equiphi29'
        # equiphi29 overrides default max_extension
        assert optimizer.max_extension == 80000
        # GC-rich genome widens max_gc acceptance
        assert optimizer.poly_config.max_gc == 0.80

    def test_init_equiphi29_at_rich_adjustment(self):
        """Test initialization with equiphi29 and AT-rich genome adjustment."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            polymerase='equiphi29',
            genome_gc_content=0.30,
        )

        # AT-rich genome widens min_gc acceptance
        assert optimizer.poly_config.min_gc == 0.20

    def test_init_custom_max_extension_overrides_preset(self):
        """Test that explicit max_extension overrides polymerase preset."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            polymerase='bst',
            max_extension=50000,
        )

        # Explicit value should be used, not bst default of 10000
        assert optimizer.max_extension == 50000

    def test_init_with_background_pruning(self):
        """Test initialization with background pruning enabled."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            bg_prefixes=['bg_genome'],
            bg_seq_lengths=[3000000],
            background_pruning=True,
            background_weight=3.0,
            min_coverage_threshold=0.90,
        )

        assert optimizer.background_pruning is True
        assert optimizer.background_weight == 3.0
        assert optimizer.min_coverage_threshold == 0.90
        assert optimizer.bg_prefixes == ['bg_genome']
        assert optimizer.bg_seq_lengths == [3000000]

    def test_init_no_gc_adjustment_without_thermo_filter(self):
        """Test that GC adjustment is skipped when polymerase has no thermo filter."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            polymerase='phi29',
            genome_gc_content=0.70,
        )

        # phi29 has thermo_filter=False, so GC adjustment should not apply
        assert optimizer.poly_config.max_gc == 0.75  # unchanged default


# =============================================================================
# Optimization Tests
# =============================================================================

class TestHybridOptimization:
    """Tests for the optimize method."""

    @pytest.fixture
    def optimizer(self):
        """Create a HybridOptimizer with mock cache for optimization tests."""
        cache = _make_mock_cache()
        return HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            bin_size=10000,
        )

    def test_basic_optimization_returns_result(self, optimizer):
        """Test that optimization returns a valid HybridResult."""
        candidates = _make_primers(20)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            verbose=False,
        )

        assert isinstance(result, HybridResult)
        assert len(result.primers) > 0
        assert len(result.primers) <= 6
        assert result.total_runtime > 0

    def test_stage1_has_more_primers_than_stage2(self, optimizer):
        """Test that Stage 1 selects more primers than the final Stage 2 count."""
        candidates = _make_primers(30)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            verbose=False,
        )

        # Stage 1 should produce a superset for Stage 2 to refine from
        assert len(result.stage1_primers) >= len(result.stage2_primers)

    def test_result_contains_expected_fields(self, optimizer):
        """Test that the result contains all expected metric fields."""
        candidates = _make_primers(15)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            verbose=False,
        )

        assert result.stage1_coverage >= 0.0
        assert result.final_coverage >= 0.0
        assert result.final_connectivity >= 0.0
        assert result.final_predicted_amplification >= 0.0
        assert result.runtime_stage1 >= 0.0
        assert result.runtime_stage2 >= 0.0

    def test_custom_stage1_count(self, optimizer):
        """Test optimization with explicit stage1_count."""
        candidates = _make_primers(25)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            stage1_count=15,
            verbose=False,
        )

        assert isinstance(result, HybridResult)
        assert len(result.primers) <= 6

    def test_fixed_primers_included_in_result(self, optimizer):
        """Test that fixed primers appear in the final result."""
        candidates = _make_primers(20)
        fixed = [candidates[0], candidates[1]]

        result = optimizer.optimize(
            candidates=candidates,
            final_count=8,
            fixed_primers=fixed,
            verbose=False,
        )

        for fp in fixed:
            assert fp.upper() in [p.upper() for p in result.primers]

    def test_fixed_primers_exceed_target_returns_fixed_only(self, optimizer):
        """Test that when fixed primers >= target, only fixed primers are returned."""
        candidates = _make_primers(20)
        fixed = candidates[:8]

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            fixed_primers=fixed,
            verbose=False,
        )

        # Should return just the fixed primers
        assert set(p.upper() for p in result.primers) == set(p.upper() for p in fixed)


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    @pytest.fixture
    def optimizer(self):
        """Create optimizer for edge case testing."""
        cache = _make_mock_cache()
        return HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            bin_size=10000,
        )

    def test_fewer_candidates_than_requested(self, optimizer):
        """Test with fewer candidates than the requested final count."""
        candidates = _make_primers(3)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=10,
            verbose=False,
        )

        # Should return all available candidates (cannot exceed input)
        assert isinstance(result, HybridResult)
        assert len(result.primers) <= 3

    def test_empty_candidate_list(self, optimizer):
        """Test with an empty candidate list."""
        result = optimizer.optimize(
            candidates=[],
            final_count=6,
            verbose=False,
        )

        assert isinstance(result, HybridResult)
        assert len(result.primers) == 0

    def test_single_candidate(self, optimizer):
        """Test with a single candidate primer."""
        candidates = _make_primers(1)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            verbose=False,
        )

        assert isinstance(result, HybridResult)
        assert len(result.primers) <= 1

    def test_candidates_equal_to_final_count(self, optimizer):
        """Test when candidates exactly equal the final count."""
        candidates = _make_primers(6)

        result = optimizer.optimize(
            candidates=candidates,
            final_count=6,
            verbose=False,
        )

        assert isinstance(result, HybridResult)
        assert len(result.primers) <= 6


# =============================================================================
# Background Pruning Tests
# =============================================================================

class TestBackgroundPruning:
    """Tests for the background pruning stage."""

    def test_background_pruning_enabled(self):
        """Test that background pruning runs when enabled with bg data."""
        # Create cache with both fg and bg positions
        positions_map = {}
        primers = _make_primers(20)

        for primer in primers:
            # Foreground: spread across genome
            h = hash(primer) % 10000
            positions_map[('genome1', primer, 'forward')] = np.array(
                [h, h + 10000, h + 20000], dtype=np.int32
            )
            positions_map[('genome1', primer, 'reverse')] = np.array(
                [h + 5000, h + 15000], dtype=np.int32
            )
            # Background: some primers have many bg sites
            bg_count = (hash(primer) % 50)
            positions_map[('bg_genome', primer, 'forward')] = np.arange(
                0, bg_count * 1000, 1000, dtype=np.int32
            )
            positions_map[('bg_genome', primer, 'reverse')] = np.array(
                [], dtype=np.int32
            )

        cache = _make_mock_cache(positions_map)

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            bg_prefixes=['bg_genome'],
            bg_seq_lengths=[3000000],
            background_pruning=True,
            background_weight=2.0,
            min_coverage_threshold=0.5,
            bin_size=10000,
        )

        result = optimizer.optimize(
            candidates=primers,
            final_count=6,
            verbose=False,
        )

        assert isinstance(result, HybridResult)
        assert len(result.primers) <= 6

    def test_background_pruning_disabled_by_default(self):
        """Test that background pruning is not applied by default."""
        cache = _make_mock_cache()

        optimizer = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
        )

        assert optimizer.background_pruning is False


# =============================================================================
# HybridResult Tests
# =============================================================================

class TestHybridResult:
    """Tests for HybridResult dataclass."""

    def test_str_representation(self):
        """Test string representation of HybridResult."""
        result = HybridResult(
            primers=['ATCGATCG', 'GCTAGCTA'],
            stage1_primers=['ATCGATCG', 'GCTAGCTA', 'AAATTTCCC'],
            stage1_coverage=0.85,
            stage1_regions_covered=8,
            stage2_primers=['ATCGATCG', 'GCTAGCTA'],
            stage2_connectivity=0.45,
            stage2_predicted_amplification=12.5,
            stage2_largest_component=5,
            final_coverage=0.80,
            final_connectivity=0.45,
            final_predicted_amplification=12.5,
            total_runtime=1.23,
        )

        s = str(result)
        assert 'Final primers: 2' in s
        assert '80.0%' in s
        assert '12.5' in s

    def test_str_with_simulation(self):
        """Test string representation includes simulation data when present."""
        sim_fitness = Mock()
        sim_fitness.mean_coverage = 0.75
        sim_fitness.std_coverage = 0.05
        sim_fitness.coverage_uniformity = 0.90
        sim_fitness.fitness_score = 0.82

        result = HybridResult(
            primers=['ATCGATCG'],
            stage1_primers=['ATCGATCG'],
            stage1_coverage=0.85,
            stage1_regions_covered=8,
            stage2_primers=['ATCGATCG'],
            stage2_connectivity=0.45,
            stage2_predicted_amplification=12.5,
            stage2_largest_component=5,
            final_coverage=0.80,
            final_connectivity=0.45,
            final_predicted_amplification=12.5,
            simulation_fitness=sim_fitness,
            total_runtime=2.0,
        )

        s = str(result)
        assert 'Simulation' in s
        assert '0.820' in s


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
