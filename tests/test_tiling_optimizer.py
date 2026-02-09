"""
Tests for the interval-based tiling optimizer.

Tests cover:
- IntervalCoverage: interval merging, coverage fraction, marginal gain
- TilingOptimizer: greedy selection, strand modes, multi-seed, specificity
"""

import pytest
import numpy as np
from unittest.mock import Mock

from neoswga.core.tiling_optimizer import (
    IntervalCoverage,
    TilingConfig,
    TilingOptimizer,
    _merge_intervals,
    _uncovered_length,
)
from neoswga.core.base_optimizer import OptimizationStatus
from neoswga.core.optimizer_factory import OptimizerFactory

# Valid DNA sequences used as primer identifiers in tests
_SEQ = [
    'ATCGATCG',     # 0
    'GCTAGCTA',     # 1
    'AATTCCGG',     # 2
    'TTGGCCAA',     # 3
    'CCGGAATT',     # 4
    'GGCCTTAA',     # 5
    'ATATATCG',     # 6
    'CGCGATAT',     # 7
    'TGCATGCA',     # 8
    'ACGTACGT',     # 9
]


# =============================================================================
# IntervalCoverage Tests
# =============================================================================


class TestIntervalCoverage:
    """Tests for IntervalCoverage class."""

    def test_empty_coverage(self):
        cov = IntervalCoverage(100000)
        assert cov.coverage_fraction == 0.0
        assert cov.uncovered_bases == 100000

    def test_single_interval(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([10000]), fragment_length=5000)
        assert cov.coverage_fraction == pytest.approx(5000 / 100000)
        assert cov.uncovered_bases == 95000

    def test_overlapping_intervals_merge(self):
        cov = IntervalCoverage(100000)
        # Two positions 2000 bp apart with 5000 bp fragments overlap
        cov.add_intervals(np.array([10000, 12000]), fragment_length=5000)
        # Merged: [10000, 17000) = 7000 bases
        assert cov.coverage_fraction == pytest.approx(7000 / 100000)

    def test_adjacent_intervals_merge(self):
        cov = IntervalCoverage(100000)
        # Adjacent: [10000, 15000) and [15000, 20000)
        cov.add_intervals(np.array([10000, 15000]), fragment_length=5000)
        assert cov.coverage_fraction == pytest.approx(10000 / 100000)

    def test_nonoverlapping_intervals_sum(self):
        cov = IntervalCoverage(100000)
        # Two non-overlapping intervals
        cov.add_intervals(np.array([10000, 50000]), fragment_length=5000)
        # [10000, 15000) + [50000, 55000) = 10000 bases
        assert cov.coverage_fraction == pytest.approx(10000 / 100000)

    def test_full_genome_coverage(self):
        cov = IntervalCoverage(10000)
        # One position at 0 with fragment covering entire genome
        cov.add_intervals(np.array([0]), fragment_length=10000)
        assert cov.coverage_fraction == pytest.approx(1.0)
        assert cov.uncovered_bases == 0

    def test_clamped_to_genome_length(self):
        cov = IntervalCoverage(10000)
        # Position near end, fragment extends past genome
        cov.add_intervals(np.array([9000]), fragment_length=5000)
        # Clamped to [9000, 10000) = 1000 bases
        assert cov.coverage_fraction == pytest.approx(1000 / 10000)

    def test_marginal_gain_returns_new_bases(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([10000]), fragment_length=5000)
        # Already covered: [10000, 15000)
        # Query: position 12000, fragment 5000 -> [12000, 17000)
        # New bases: [15000, 17000) = 2000
        gain = cov.marginal_gain(np.array([12000]), fragment_length=5000)
        assert gain == 2000

    def test_marginal_gain_no_overlap(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([10000]), fragment_length=5000)
        # Query in completely different region
        gain = cov.marginal_gain(np.array([50000]), fragment_length=5000)
        assert gain == 5000

    def test_marginal_gain_full_overlap(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([10000]), fragment_length=5000)
        # Query same region
        gain = cov.marginal_gain(np.array([10000]), fragment_length=5000)
        assert gain == 0

    def test_marginal_gain_does_not_modify_state(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([10000]), fragment_length=5000)
        original_covered = cov._covered_bases
        original_intervals = list(cov._intervals)

        cov.marginal_gain(np.array([50000]), fragment_length=5000)

        assert cov._covered_bases == original_covered
        assert cov._intervals == original_intervals

    def test_incremental_additions(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([0]), fragment_length=10000)
        assert cov.coverage_fraction == pytest.approx(0.1)

        cov.add_intervals(np.array([20000]), fragment_length=10000)
        assert cov.coverage_fraction == pytest.approx(0.2)

        cov.add_intervals(np.array([10000]), fragment_length=10000)
        # Now [0, 30000) is covered = 0.3
        assert cov.coverage_fraction == pytest.approx(0.3)

    def test_empty_positions(self):
        cov = IntervalCoverage(100000)
        cov.add_intervals(np.array([]), fragment_length=5000)
        assert cov.coverage_fraction == 0.0

    def test_zero_length_genome(self):
        cov = IntervalCoverage(0)
        assert cov.coverage_fraction == 0.0
        cov.add_intervals(np.array([0]), fragment_length=5000)
        assert cov.coverage_fraction == 0.0


# =============================================================================
# Helper Function Tests
# =============================================================================


class TestMergeIntervals:
    """Tests for _merge_intervals helper."""

    def test_empty(self):
        assert _merge_intervals([]) == []

    def test_single(self):
        assert _merge_intervals([(0, 10)]) == [(0, 10)]

    def test_nonoverlapping(self):
        result = _merge_intervals([(0, 10), (20, 30)])
        assert result == [(0, 10), (20, 30)]

    def test_overlapping(self):
        result = _merge_intervals([(0, 15), (10, 25)])
        assert result == [(0, 25)]

    def test_contained(self):
        result = _merge_intervals([(0, 30), (10, 20)])
        assert result == [(0, 30)]

    def test_unsorted_input(self):
        result = _merge_intervals([(20, 30), (0, 10), (5, 25)])
        assert result == [(0, 30)]


class TestUncoveredLength:
    """Tests for _uncovered_length helper."""

    def test_no_existing(self):
        assert _uncovered_length([], 0, 100) == 100

    def test_fully_covered(self):
        assert _uncovered_length([(0, 200)], 50, 150) == 0

    def test_partial_overlap_left(self):
        # Existing [0, 50), query [25, 75): uncovered is [50, 75) = 25
        assert _uncovered_length([(0, 50)], 25, 75) == 25

    def test_partial_overlap_right(self):
        # Existing [50, 100), query [25, 75): uncovered is [25, 50) = 25
        assert _uncovered_length([(50, 100)], 25, 75) == 25

    def test_gap_in_middle(self):
        # Existing [0, 30) and [70, 100), query [0, 100)
        # Uncovered: [30, 70) = 40
        assert _uncovered_length([(0, 30), (70, 100)], 0, 100) == 40


# =============================================================================
# TilingOptimizer Tests
# =============================================================================


def _make_mock_cache(positions_map):
    """Create a mock PositionCache from a primer->positions mapping.

    positions_map: dict of primer -> positions (used for all prefixes/strands)
    """
    cache = Mock()

    def get_positions(prefix, primer, strand):
        return positions_map.get(primer, np.array([], dtype=np.int64))

    cache.get_positions = get_positions
    return cache


def _make_strand_mock_cache(positions_map):
    """Mock cache where positions depend on strand.

    positions_map: dict of (primer, strand) -> positions
    Falls back to primer-only key, then empty array.
    """
    cache = Mock()

    def get_positions(prefix, primer, strand):
        key = (primer, strand)
        if key in positions_map:
            return positions_map[key]
        if primer in positions_map:
            return positions_map[primer]
        return np.array([], dtype=np.int64)

    cache.get_positions = get_positions
    return cache


class TestTilingOptimizer:
    """Tests for TilingOptimizer class."""

    def test_returns_optimization_result(self):
        cache = _make_mock_cache({
            _SEQ[0]: np.array([0, 10000, 20000, 30000]),
            _SEQ[1]: np.array([40000, 50000, 60000]),
        })
        config = TilingConfig(
            fragment_length=10000,
            n_seeds=1,
            strand_separate=False,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize([_SEQ[0], _SEQ[1]], target_size=6)
        assert hasattr(result, 'primers')
        assert hasattr(result, 'score')
        assert hasattr(result, 'status')
        assert hasattr(result, 'metrics')

    def test_respects_target_size(self):
        primers = {_SEQ[i]: np.array([i * 5000]) for i in range(10)}
        cache = _make_mock_cache(primers)
        config = TilingConfig(
            fragment_length=5000,
            n_seeds=1,
            strand_separate=False,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize(list(primers.keys()), target_size=3)
        assert len(result.primers) <= 3

    def test_achieves_coverage(self):
        cache = _make_mock_cache({
            _SEQ[0]: np.array([0]),
            _SEQ[1]: np.array([5000]),
        })
        config = TilingConfig(
            fragment_length=5000,
            n_seeds=1,
            strand_separate=False,
            min_coverage=0.99,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[10000],
            config=config,
        )
        result = opt.optimize([_SEQ[0], _SEQ[1]], target_size=10)
        assert result.score == pytest.approx(1.0)
        assert result.status == OptimizationStatus.SUCCESS

    def test_strand_separate_mode(self):
        positions = {
            (_SEQ[0], 'forward'): np.array([0, 10000]),
            (_SEQ[0], 'reverse'): np.array([20000]),
            (_SEQ[0], 'both'): np.array([0, 10000, 20000]),
            (_SEQ[1], 'forward'): np.array([30000]),
            (_SEQ[1], 'reverse'): np.array([40000, 50000]),
            (_SEQ[1], 'both'): np.array([30000, 40000, 50000]),
        }
        cache = _make_strand_mock_cache(positions)
        config = TilingConfig(
            fragment_length=10000,
            n_seeds=1,
            strand_separate=True,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[60000],
            config=config,
        )
        result = opt.optimize([_SEQ[0], _SEQ[1]], target_size=10)
        assert len(result.primers) > 0

    def test_strand_combined_mode(self):
        cache = _make_mock_cache({
            _SEQ[0]: np.array([0, 20000, 40000]),
            _SEQ[1]: np.array([60000, 80000]),
        })
        config = TilingConfig(
            fragment_length=20000,
            n_seeds=1,
            strand_separate=False,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize([_SEQ[0], _SEQ[1]], target_size=10)
        assert len(result.primers) > 0
        assert result.score > 0.0

    def test_multi_seed_runs(self):
        primers = {_SEQ[i]: np.array([i * 10000]) for i in range(10)}
        cache = _make_mock_cache(primers)
        config = TilingConfig(
            fragment_length=10000,
            n_seeds=3,
            strand_separate=False,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize(list(primers.keys()), target_size=5)
        assert len(result.primers) > 0
        assert result.iterations == 3

    def test_specificity_penalizes_high_background(self):
        fg_map = {
            _SEQ[0]: np.array([0, 10000, 20000]),
            _SEQ[1]: np.array([0, 10000, 20000]),
        }
        bg_map = {
            _SEQ[0]: np.array([]),
            _SEQ[1]: np.array(list(range(0, 100000, 100))),  # 1000 bg sites
        }
        cache = Mock()

        def get_positions(prefix, primer, strand):
            if prefix == 'fg':
                return fg_map.get(primer, np.array([], dtype=np.int64))
            elif prefix == 'bg':
                return bg_map.get(primer, np.array([], dtype=np.int64))
            return np.array([], dtype=np.int64)

        cache.get_positions = get_positions

        config = TilingConfig(
            fragment_length=10000,
            n_seeds=1,
            strand_separate=False,
            specificity_weight=5.0,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['fg'],
            fg_seq_lengths=[50000],
            bg_prefixes=['bg'],
            bg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize([_SEQ[0], _SEQ[1]], target_size=1)
        # _SEQ[0] should be preferred because _SEQ[1] has high background
        assert _SEQ[0] in result.primers

    def test_adaptive_threshold_relaxation(self):
        # _SEQ[0] covers a lot, _SEQ[1] covers very little
        cache = _make_mock_cache({
            _SEQ[0]: np.array([0, 10000, 20000, 30000, 40000]),
            _SEQ[1]: np.array([90000]),
        })
        config = TilingConfig(
            fragment_length=10000,
            n_seeds=1,
            strand_separate=False,
            adaptive_thresholds=(0.50, 0.25, 0.10, 0.0),
            min_coverage=0.99,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize([_SEQ[0], _SEQ[1]], target_size=10)
        # Both should be selected (_SEQ[1] via threshold relaxation)
        assert _SEQ[0] in result.primers
        assert _SEQ[1] in result.primers

    def test_empty_candidates_raises(self):
        cache = _make_mock_cache({})
        config = TilingConfig(strand_separate=False, n_seeds=1)
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        with pytest.raises(ValueError):
            opt.optimize([], target_size=5)

    def test_fragment_length_affects_coverage(self):
        cache = _make_mock_cache({
            _SEQ[0]: np.array([50000]),
        })
        config_short = TilingConfig(
            fragment_length=1000, n_seeds=1, strand_separate=False,
        )
        opt_short = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config_short,
        )
        result_short = opt_short.optimize([_SEQ[0]], target_size=5)

        config_long = TilingConfig(
            fragment_length=100000, n_seeds=1, strand_separate=False,
        )
        opt_long = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config_long,
        )
        result_long = opt_long.optimize([_SEQ[0]], target_size=5)

        assert result_long.score > result_short.score

    def test_fixed_primers_included(self):
        cache = _make_mock_cache({
            _SEQ[0]: np.array([0]),
            _SEQ[1]: np.array([50000]),
            _SEQ[2]: np.array([80000]),
        })
        config = TilingConfig(
            fragment_length=10000, n_seeds=1, strand_separate=False,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize(
            [_SEQ[0], _SEQ[1], _SEQ[2]],
            target_size=5,
            fixed_primers=[_SEQ[0]],
        )
        assert _SEQ[0] in result.primers

    def test_partial_status_below_min_coverage(self):
        cache = _make_mock_cache({
            _SEQ[0]: np.array([0]),
        })
        config = TilingConfig(
            fragment_length=1000,
            n_seeds=1,
            strand_separate=False,
            min_coverage=0.95,
        )
        opt = TilingOptimizer(
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
            config=config,
        )
        result = opt.optimize([_SEQ[0]], target_size=1)
        assert result.status == OptimizationStatus.PARTIAL


# =============================================================================
# Registration Tests
# =============================================================================


class TestTilingRegistration:
    """Tests for optimizer factory registration."""

    def test_registered_as_tiling(self):
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()
        assert OptimizerFactory.list_optimizers().get('tiling') is not None

    def test_alias_coverage_tiling(self):
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()
        from neoswga.core.optimizer_factory import OptimizerRegistry
        assert OptimizerRegistry.is_registered('coverage-tiling')

    def test_alias_interval_tiling(self):
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()
        from neoswga.core.optimizer_factory import OptimizerRegistry
        assert OptimizerRegistry.is_registered('interval-tiling')

    def test_factory_creates_instance(self):
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()
        cache = _make_mock_cache({})
        optimizer = OptimizerFactory.create(
            'tiling',
            position_cache=cache,
            fg_prefixes=['g1'],
            fg_seq_lengths=[100000],
        )
        assert optimizer.name == 'tiling'


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
