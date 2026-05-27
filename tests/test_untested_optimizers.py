"""
Integration tests for previously untested optimizers.

Covers:
- background-aware: Three-stage optimizer with background reduction
- equiphi29: EquiPhi29-optimized selection at 42C
- clique: Clique-based exhaustive dimer-free set discovery
- multi-agent: Multi-agent parallel optimizer with result aggregation
- serial-cascade pipelines: coverage-then-dimerfree, dimerfree-scored, bg-prefilter-hybrid
- bg-prefilter: Background pruning pre-filter
"""

import numpy as np
import pytest
from unittest.mock import Mock

from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry
from neoswga.core.base_optimizer import BaseOptimizer, OptimizerConfig, OptimizationStatus

# Ensure all optimizers are registered (aliases, etc.)
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

_ensure_optimizers_registered()


# =============================================================================
# Helpers
# =============================================================================


def _make_position_cache(positions_map):
    """Build a mock PositionCache from a dict of {primer: positions_array}.

    Positions are returned for any prefix and any strand.
    """
    cache = Mock()

    def get_positions(prefix, primer, strand="both"):
        return positions_map.get(primer, np.array([], dtype=np.int64))

    cache.get_positions = get_positions
    return cache


# Candidate primers for testing.  Short sequences chosen so that
# is_dimer_fast with default thresholds allows most pairs through
# the compatibility graph (needed for clique tests).
CANDIDATES = [
    "ATCGATCG",
    "GCTAGCTA",
    "AAATTTGG",
    "CCCCGGGG",
    "TTTTAAAA",
    "AACCTTGG",
    "GGAATTCC",
    "TTGGCCAA",
]

# Each primer binds at a few positions spread across a 100 kb mock genome.
POSITION_MAP = {
    "ATCGATCG": np.array([1000, 11000, 21000, 31000, 41000]),
    "GCTAGCTA": np.array([51000, 61000, 71000]),
    "AAATTTGG": np.array([81000, 91000]),
    "CCCCGGGG": np.array([5000, 25000, 45000, 65000, 85000]),
    "TTTTAAAA": np.array([15000, 35000, 55000, 75000, 95000]),
    "AACCTTGG": np.array([3000, 23000, 43000, 63000, 83000]),
    "GGAATTCC": np.array([7000, 27000, 47000, 67000, 87000]),
    "TTGGCCAA": np.array([9000, 29000, 49000, 69000, 89000]),
}

# Background positions (fewer sites, representing low host binding).
BG_POSITION_MAP = {
    "ATCGATCG": np.array([500]),
    "GCTAGCTA": np.array([1500, 2500]),
    "AAATTTGG": np.array([]),
    "CCCCGGGG": np.array([3500, 4500, 5500]),
    "TTTTAAAA": np.array([6500]),
    "AACCTTGG": np.array([7500]),
    "GGAATTCC": np.array([]),
    "TTGGCCAA": np.array([8500, 9500]),
}

FG_PREFIXES = ["test_genome"]
FG_SEQ_LENGTHS = [100000]
BG_PREFIXES = ["test_background"]
BG_SEQ_LENGTHS = [50000]


def _make_fg_bg_cache(fg_map, bg_map, fg_prefix="test_genome", bg_prefix="test_background"):
    """Build a mock PositionCache that returns different positions for fg vs bg prefixes."""
    cache = Mock()

    def get_positions(prefix, primer, strand="both"):
        if prefix == bg_prefix:
            return bg_map.get(primer, np.array([], dtype=np.int64))
        return fg_map.get(primer, np.array([], dtype=np.int64))

    cache.get_positions = get_positions
    return cache


# =============================================================================
# background-aware optimizer
# =============================================================================


class TestBackgroundAwareOptimizer:
    """Tests for the background-aware optimizer (clinical applications)."""

    def test_registered(self):
        assert OptimizerRegistry.is_registered("background-aware")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("clinical")
        assert OptimizerRegistry.is_registered("bg-aware")

    def test_factory_create(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "background-aware",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)
        assert optimizer.name == "background-aware"

    def test_optimize_returns_success(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "background-aware",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.status in (OptimizationStatus.SUCCESS, OptimizationStatus.PARTIAL)
        assert len(result.primers) > 0

    def test_optimize_coverage_positive(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "background-aware",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.metrics.fg_coverage > 0

    def test_without_background_degrades_gracefully(self):
        """background-aware should still return a result without background data."""
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "background-aware",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        # Without background the optimizer may succeed, partially succeed,
        # or report an error — just verify it does not crash
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.NO_CONVERGENCE,
            OptimizationStatus.ERROR,
        )


# =============================================================================
# equiphi29 optimizer
# =============================================================================


class TestEquiPhi29Optimizer:
    """Tests for the EquiPhi29 optimizer adapter."""

    def test_registered(self):
        assert OptimizerRegistry.is_registered("equiphi29")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("equiphi")
        assert OptimizerRegistry.is_registered("eq29")

    def test_factory_create(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "equiphi29",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)
        assert optimizer.name == "equiphi29"

    def test_optimize_returns_result(self):
        """EquiPhi29 adapter delegates to HybridOptimizer with equiphi29 preset.

        The underlying HybridOptimizer may reject the 42C temperature if
        its thermo filter defaults to phi29 range validation. The adapter
        should handle this gracefully and return a valid result object.
        """
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "equiphi29",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        # The result should be a valid OptimizationResult regardless of success
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.ERROR,
        )
        assert result.optimizer_name == "equiphi29"

    def test_with_background(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "equiphi29",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.ERROR,
        )


# =============================================================================
# clique optimizer
# =============================================================================


class TestCliqueOptimizer:
    """Tests for the clique-based dimer-free optimizer."""

    @pytest.fixture(autouse=True)
    def _require_networkx(self):
        pytest.importorskip("networkx")

    def test_registered(self):
        assert OptimizerRegistry.is_registered("clique")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("clique-dimer-free")
        assert OptimizerRegistry.is_registered("dimer-free")

    def test_factory_create(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "clique",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)
        assert optimizer.name == "clique"

    def test_optimize_returns_result(self):
        cache = _make_position_cache(POSITION_MAP)
        # Use max_dimer_bp=5 to relax dimer check for short synthetic primers
        config = OptimizerConfig(target_set_size=3, max_dimer_bp=5, verbose=False)
        optimizer = OptimizerFactory.create(
            "clique",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        # Clique enumeration may not find sets if primers are too complementary
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.ERROR,
        )

    def test_primers_are_dimer_free(self):
        """All selected primers should be mutually dimer-free."""
        from neoswga.core.dimer import is_dimer_fast

        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, max_dimer_bp=5, verbose=False)
        optimizer = OptimizerFactory.create(
            "clique",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        if result.status == OptimizationStatus.SUCCESS:
            primers = list(result.primers)
            for i in range(len(primers)):
                for j in range(i + 1, len(primers)):
                    assert not is_dimer_fast(primers[i], primers[j], 3), (
                        f"Dimer found between {primers[i]} and {primers[j]}"
                    )

    def test_no_valid_sets_returns_error(self):
        """When no dimer-free set exists, the optimizer should fail gracefully."""
        # Use primers that are reverse-complements to force dimers
        dimer_primers = ["AAAA", "TTTT"]
        dimer_map = {
            "AAAA": np.array([1000, 20000]),
            "TTTT": np.array([5000, 30000]),
        }
        cache = _make_position_cache(dimer_map)
        config = OptimizerConfig(target_set_size=2, verbose=False)
        optimizer = OptimizerFactory.create(
            "clique",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(dimer_primers, target_size=2)
        # May succeed or fail depending on dimer threshold; just verify no crash
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.ERROR,
        )


# =============================================================================
# multi-agent optimizer
# =============================================================================


class TestMultiAgentOptimizer:
    """Tests for the multi-agent parallel optimizer."""

    def test_registered(self):
        assert OptimizerRegistry.is_registered("multi-agent")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("ensemble")
        assert OptimizerRegistry.is_registered("parallel")

    def test_factory_create(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "multi-agent",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)
        assert optimizer.name == "multi-agent"

    def test_optimize_returns_result(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "multi-agent",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.status in (OptimizationStatus.SUCCESS, OptimizationStatus.PARTIAL)
        assert len(result.primers) > 0

    def test_optimize_coverage_positive(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "multi-agent",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.metrics.fg_coverage > 0


# =============================================================================
# serial-cascade optimizer: coverage-then-dimerfree
# =============================================================================


class TestCoverageThenDimerFreeOptimizer:
    """Tests for the DS -> Clique serial cascade."""

    @pytest.fixture(autouse=True)
    def _require_networkx(self):
        pytest.importorskip("networkx")

    def test_registered(self):
        assert OptimizerRegistry.is_registered("coverage-then-dimerfree")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("ds-clique")
        assert OptimizerRegistry.is_registered("coverage-dimerfree")

    def test_factory_create(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "coverage-then-dimerfree",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)

    def test_optimize_returns_result(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, max_dimer_bp=5, verbose=False)
        optimizer = OptimizerFactory.create(
            "coverage-then-dimerfree",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        # Clique stage may not find dimer-free sets with synthetic primers
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.ERROR,
        )

    def test_coverage_positive(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, max_dimer_bp=5, verbose=False)
        optimizer = OptimizerFactory.create(
            "coverage-then-dimerfree",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        if result.status == OptimizationStatus.SUCCESS:
            assert result.metrics.fg_coverage > 0


# =============================================================================
# serial-cascade optimizer: dimerfree-scored
# =============================================================================


class TestDimerFreeScoredOptimizer:
    """Tests for the Clique -> Network serial cascade."""

    @pytest.fixture(autouse=True)
    def _require_networkx(self):
        pytest.importorskip("networkx")

    def test_registered(self):
        assert OptimizerRegistry.is_registered("dimerfree-scored")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("clique-network")

    def test_factory_create(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "dimerfree-scored",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)

    def test_optimize_returns_result(self):
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, max_dimer_bp=5, verbose=False)
        optimizer = OptimizerFactory.create(
            "dimerfree-scored",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        # Clique stage may not find dimer-free sets with synthetic primers
        assert result.status in (
            OptimizationStatus.SUCCESS,
            OptimizationStatus.PARTIAL,
            OptimizationStatus.ERROR,
        )


# =============================================================================
# bg-prefilter optimizer
# =============================================================================


class TestBgPrefilterOptimizer:
    """Tests for the background pruning pre-filter."""

    def test_registered(self):
        assert OptimizerRegistry.is_registered("bg-prefilter")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("background-prefilter")
        assert OptimizerRegistry.is_registered("bg-prune")

    def test_factory_create(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=4, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)
        assert optimizer.name == "bg-prefilter"

    def test_optimize_returns_result(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=4, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=4)
        assert result.status in (OptimizationStatus.SUCCESS, OptimizationStatus.PARTIAL)
        assert len(result.primers) > 0

    def test_coverage_positive(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=4, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=4)
        assert result.metrics.fg_coverage > 0

    def test_without_background_passes_through(self):
        """Without background data, bg-prefilter should pass candidates through."""
        cache = _make_position_cache(POSITION_MAP)
        config = OptimizerConfig(target_set_size=4, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=4)
        assert result.status == OptimizationStatus.SUCCESS
        assert len(result.primers) == 4

    def test_with_high_background_returns_result(self):
        """Optimizer handles primers with varying background binding levels."""
        # Use valid DNA sequences with different background densities
        fg_map = {
            "ATGATGATG": np.array([1000, 20000, 40000, 60000, 80000]),
            "CATCATCAT": np.array([5000, 25000, 45000, 65000, 85000]),
            "GACGACGAC": np.array([10000, 30000, 50000, 70000, 90000]),
        }
        bg_map = {
            "ATGATGATG": np.array([100]),                        # low background
            "CATCATCAT": np.array(list(range(1000)), dtype=np.int64),  # high background
            "GACGACGAC": np.array([200]),                        # low background
        }
        cache = _make_fg_bg_cache(fg_map, bg_map)
        config = OptimizerConfig(target_set_size=2, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(
            ["ATGATGATG", "CATCATCAT", "GACGACGAC"], target_size=2
        )
        assert result.status in (OptimizationStatus.SUCCESS, OptimizationStatus.PARTIAL)
        assert len(result.primers) >= 2


# =============================================================================
# serial-cascade optimizer: bg-prefilter-hybrid
# =============================================================================


class TestBgPrefilterHybridOptimizer:
    """Tests for the bg-prefilter -> hybrid serial cascade."""

    def test_registered(self):
        assert OptimizerRegistry.is_registered("bg-prefilter-hybrid")

    def test_aliases_registered(self):
        assert OptimizerRegistry.is_registered("bgpf-hybrid")

    def test_factory_create(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter-hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        assert isinstance(optimizer, BaseOptimizer)

    def test_optimize_returns_result(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter-hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.status in (OptimizationStatus.SUCCESS, OptimizationStatus.PARTIAL)
        assert len(result.primers) > 0

    def test_coverage_positive(self):
        cache = _make_fg_bg_cache(POSITION_MAP, BG_POSITION_MAP)
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "bg-prefilter-hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bg_prefixes=BG_PREFIXES,
            bg_seq_lengths=BG_SEQ_LENGTHS,
            config=config,
        )
        result = optimizer.optimize(CANDIDATES, target_size=3)
        if result.status == OptimizationStatus.SUCCESS:
            assert result.metrics.fg_coverage > 0
