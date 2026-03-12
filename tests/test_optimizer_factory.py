"""
Tests for optimizer_factory and hybrid_optimizer modules.

Covers:
- OptimizerRegistry: registration, lookup, aliases, error handling
- OptimizerFactory: creation, validation, from_params
- HybridBaseOptimizer: instantiation via factory and end-to-end optimization
- HybridOptimizer: direct two-stage optimization with synthetic data
"""

import pytest
import numpy as np
from unittest.mock import Mock

from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry
from neoswga.core.base_optimizer import BaseOptimizer, OptimizerConfig, OptimizationStatus
from neoswga.core.exceptions import OptimizerNotFoundError

# Import optimizer modules to trigger @OptimizerFactory.register() decorators.
# Without these imports the registry remains empty because registration
# happens at module load time via decorators in each optimizer file.
import neoswga.core.hybrid_optimizer  # noqa: F401
import neoswga.core.greedy_optimizer  # noqa: F401
import neoswga.core.dominating_set_adapter  # noqa: F401
import neoswga.core.network_optimizer  # noqa: F401
import neoswga.core.genetic_algorithm  # noqa: F401


# =============================================================================
# Helpers
# =============================================================================

def _make_position_cache(positions_map):
    """
    Build a mock PositionCache from a dict of {primer: positions_array}.

    Positions are returned for any prefix and any strand.
    """
    cache = Mock()

    def get_positions(prefix, primer, strand='both'):
        return positions_map.get(primer, np.array([], dtype=np.int64))

    cache.get_positions = get_positions
    return cache


# Candidate primers that spell out simple DNA sequences for testing.
CANDIDATES = ["ATCGATCG", "GCTAGCTA", "AAATTTGG", "CCCCGGGG", "TTTTAAAA"]

# Each primer binds at a few positions spread across a 100 kb mock genome.
POSITION_MAP = {
    "ATCGATCG": np.array([1000, 11000, 21000, 31000, 41000]),
    "GCTAGCTA": np.array([51000, 61000, 71000]),
    "AAATTTGG": np.array([81000, 91000]),
    "CCCCGGGG": np.array([5000, 25000, 45000, 65000, 85000]),
    "TTTTAAAA": np.array([15000, 35000, 55000, 75000, 95000]),
}

FG_PREFIXES = ["test_genome"]
FG_SEQ_LENGTHS = [100000]


# =============================================================================
# OptimizerRegistry Tests
# =============================================================================

class TestOptimizerRegistry:
    """Tests for the optimizer registry singleton."""

    def test_hybrid_is_registered(self):
        """The hybrid optimizer must be registered at import time."""
        assert OptimizerRegistry.is_registered("hybrid")

    def test_hybrid_aliases_registered(self):
        """Hybrid optimizer aliases should resolve correctly."""
        assert OptimizerRegistry.is_registered("hybrid-optimizer")
        assert OptimizerRegistry.is_registered("two-stage")

    def test_greedy_is_registered(self):
        """Greedy optimizer should be registered."""
        assert OptimizerRegistry.is_registered("greedy")

    def test_dominating_set_is_registered(self):
        """Dominating-set optimizer should be registered."""
        assert OptimizerRegistry.is_registered("dominating-set")

    def test_network_is_registered(self):
        """Network optimizer should be registered."""
        assert OptimizerRegistry.is_registered("network")

    def test_genetic_is_registered(self):
        """Genetic algorithm optimizer should be registered."""
        assert OptimizerRegistry.is_registered("genetic")

    def test_get_returns_base_optimizer_subclass(self):
        """Registry.get should return a class that is a BaseOptimizer subclass."""
        cls = OptimizerRegistry.get("hybrid")
        assert issubclass(cls, BaseOptimizer)

    def test_get_unknown_raises_optimizer_not_found(self):
        """Looking up an unregistered name should raise OptimizerNotFoundError."""
        with pytest.raises(OptimizerNotFoundError):
            OptimizerRegistry.get("nonexistent-optimizer")

    def test_alias_resolves_to_same_class(self):
        """An alias should resolve to the same class as the canonical name."""
        hybrid_cls = OptimizerRegistry.get("hybrid")
        alias_cls = OptimizerRegistry.get("two-stage")
        assert hybrid_cls is alias_cls

    def test_list_all_contains_hybrid(self):
        """list_all should include 'hybrid' in the returned dict."""
        all_optimizers = OptimizerRegistry.list_all()
        assert "hybrid" in all_optimizers

    def test_list_all_returns_descriptions(self):
        """Descriptions should be non-empty strings."""
        all_optimizers = OptimizerRegistry.list_all()
        for name, desc in all_optimizers.items():
            assert isinstance(desc, str)
            assert len(desc) > 0, f"Empty description for {name}"


# =============================================================================
# OptimizerFactory.create Tests
# =============================================================================

class TestOptimizerFactoryCreate:
    """Tests for OptimizerFactory.create."""

    @pytest.fixture
    def cache(self):
        return _make_position_cache(POSITION_MAP)

    def test_create_hybrid(self, cache):
        """Factory should create a hybrid optimizer instance."""
        optimizer = OptimizerFactory.create(
            "hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        assert optimizer.name == "hybrid"
        assert isinstance(optimizer, BaseOptimizer)

    def test_create_with_alias(self, cache):
        """Factory should accept aliases."""
        optimizer = OptimizerFactory.create(
            "two-stage",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        assert optimizer.name == "hybrid"

    def test_create_greedy(self, cache):
        """Factory should create a greedy optimizer."""
        optimizer = OptimizerFactory.create(
            "greedy",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        assert optimizer.name == "greedy"

    def test_create_dominating_set(self, cache):
        """Factory should create a dominating-set optimizer."""
        optimizer = OptimizerFactory.create(
            "dominating-set",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        assert optimizer.name == "dominating-set"

    def test_create_unknown_raises(self, cache):
        """Creating with unknown name should raise OptimizerNotFoundError."""
        with pytest.raises(OptimizerNotFoundError):
            OptimizerFactory.create(
                "does-not-exist",
                position_cache=cache,
                fg_prefixes=FG_PREFIXES,
                fg_seq_lengths=FG_SEQ_LENGTHS,
            )

    def test_create_empty_name_raises(self, cache):
        """Empty optimizer name should raise ValueError."""
        with pytest.raises(ValueError, match="cannot be empty"):
            OptimizerFactory.create(
                "",
                position_cache=cache,
                fg_prefixes=FG_PREFIXES,
                fg_seq_lengths=FG_SEQ_LENGTHS,
            )

    def test_create_none_cache_raises(self, cache):
        """None position_cache should raise ValueError."""
        with pytest.raises(ValueError, match="position_cache"):
            OptimizerFactory.create(
                "hybrid",
                position_cache=None,
                fg_prefixes=FG_PREFIXES,
                fg_seq_lengths=FG_SEQ_LENGTHS,
            )

    def test_create_empty_prefixes_raises(self, cache):
        """Empty fg_prefixes should raise ValueError."""
        with pytest.raises(ValueError, match="fg_prefixes"):
            OptimizerFactory.create(
                "hybrid",
                position_cache=cache,
                fg_prefixes=[],
                fg_seq_lengths=[],
            )

    def test_create_mismatched_lengths_raises(self, cache):
        """Mismatched fg_prefixes and fg_seq_lengths should raise ValueError."""
        with pytest.raises(ValueError, match="same length"):
            OptimizerFactory.create(
                "hybrid",
                position_cache=cache,
                fg_prefixes=["g1", "g2"],
                fg_seq_lengths=[100000],
            )

    def test_create_with_config(self, cache):
        """Custom OptimizerConfig should be passed through."""
        config = OptimizerConfig(target_set_size=4, verbose=False)
        optimizer = OptimizerFactory.create(
            "hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        assert optimizer.config.target_set_size == 4
        assert optimizer.config.verbose is False


# =============================================================================
# OptimizerFactory.list_optimizers Tests
# =============================================================================

class TestOptimizerFactoryList:
    """Tests for OptimizerFactory.list_optimizers."""

    def test_list_returns_dict(self):
        """list_optimizers should return a dict."""
        result = OptimizerFactory.list_optimizers()
        assert isinstance(result, dict)

    def test_list_contains_core_optimizers(self):
        """All core optimizers should be listed."""
        result = OptimizerFactory.list_optimizers()
        expected = {"hybrid", "greedy", "dominating-set", "network", "genetic"}
        assert expected.issubset(set(result.keys()))


# =============================================================================
# OptimizerFactory.get_optimizer_info Tests
# =============================================================================

class TestOptimizerFactoryInfo:
    """Tests for OptimizerFactory.get_optimizer_info."""

    def test_info_for_hybrid(self):
        """get_optimizer_info should return metadata for hybrid."""
        info = OptimizerFactory.get_optimizer_info("hybrid")
        assert info["name"] == "hybrid"
        assert "class" in info
        assert "description" in info
        assert "module" in info

    def test_info_unknown_raises(self):
        """get_optimizer_info for unknown name should raise."""
        with pytest.raises(OptimizerNotFoundError):
            OptimizerFactory.get_optimizer_info("nonexistent")


# =============================================================================
# HybridBaseOptimizer (via BaseOptimizer interface) Tests
# =============================================================================

class TestHybridBaseOptimizer:
    """Tests for the hybrid optimizer accessed through the BaseOptimizer interface."""

    @pytest.fixture
    def cache(self):
        return _make_position_cache(POSITION_MAP)

    @pytest.fixture
    def optimizer(self, cache):
        """Create a HybridBaseOptimizer via the factory."""
        config = OptimizerConfig(target_set_size=3, verbose=False)
        return OptimizerFactory.create(
            "hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )

    def test_name(self, optimizer):
        """Optimizer should report its name."""
        assert optimizer.name == "hybrid"

    def test_description(self, optimizer):
        """Optimizer should have a non-empty description."""
        assert len(optimizer.description) > 0

    def test_optimize_returns_result(self, optimizer):
        """optimize() should return an OptimizationResult."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert hasattr(result, "primers")
        assert hasattr(result, "score")
        assert hasattr(result, "status")
        assert hasattr(result, "metrics")
        assert hasattr(result, "optimizer_name")

    def test_optimize_selects_primers(self, optimizer):
        """optimize() should select a non-empty primer set."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert len(result.primers) > 0

    def test_optimize_respects_target_size(self, optimizer):
        """Selected set size should not exceed target_size."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert len(result.primers) <= 3

    def test_optimize_status_is_success(self, optimizer):
        """Optimization on valid data should succeed."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result.status == OptimizationStatus.SUCCESS

    def test_optimize_returns_valid_metrics(self, optimizer):
        """Metrics should contain coverage and selectivity information."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        metrics = result.metrics
        assert 0.0 <= metrics.fg_coverage <= 1.0
        assert metrics.total_fg_sites >= 0
        assert metrics.mean_tm > 0

    def test_optimize_primers_are_from_candidates(self, optimizer):
        """All selected primers should come from the candidate pool."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        for primer in result.primers:
            assert primer in CANDIDATES

    def test_optimize_empty_candidates_raises(self, optimizer):
        """optimize() with empty candidates should raise ValueError."""
        with pytest.raises(ValueError):
            optimizer.optimize([], target_size=3)

    def test_optimize_with_fixed_primers(self, optimizer):
        """Fixed primers should appear in the result."""
        fixed = ["ATCGATCG"]
        result = optimizer.optimize(CANDIDATES, target_size=3, fixed_primers=fixed)
        assert "ATCGATCG" in result.primers

    def test_optimize_result_serializable(self, optimizer):
        """OptimizationResult.to_dict should produce a dict."""
        result = optimizer.optimize(CANDIDATES, target_size=3)
        d = result.to_dict()
        assert isinstance(d, dict)
        assert "primers" in d
        assert "score" in d
        assert "status" in d


# =============================================================================
# HybridOptimizer (direct, non-BaseOptimizer interface) Tests
# =============================================================================

class TestHybridOptimizerDirect:
    """Tests for the HybridOptimizer class used directly."""

    @pytest.fixture
    def cache(self):
        return _make_position_cache(POSITION_MAP)

    @pytest.fixture
    def hybrid(self, cache):
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        return HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            bin_size=10000,
        )

    def test_instantiation(self, hybrid):
        """HybridOptimizer should instantiate without errors."""
        assert hybrid.fg_prefixes == FG_PREFIXES
        assert hybrid.fg_seq_lengths == FG_SEQ_LENGTHS

    def test_optimize_returns_hybrid_result(self, hybrid):
        """optimize() should return a HybridResult dataclass."""
        from neoswga.core.hybrid_optimizer import HybridResult
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)
        assert isinstance(result, HybridResult)

    def test_hybrid_result_has_expected_fields(self, hybrid):
        """HybridResult should contain stage-level and overall metrics."""
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)

        assert hasattr(result, "primers")
        assert hasattr(result, "stage1_primers")
        assert hasattr(result, "stage1_coverage")
        assert hasattr(result, "stage2_primers")
        assert hasattr(result, "final_coverage")
        assert hasattr(result, "final_connectivity")
        assert hasattr(result, "final_predicted_amplification")
        assert hasattr(result, "total_runtime")

    def test_optimize_selects_primers(self, hybrid):
        """Should select a non-empty primer set."""
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)
        assert len(result.primers) > 0

    def test_optimize_final_count(self, hybrid):
        """Selected set should not exceed final_count."""
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)
        assert len(result.primers) <= 3

    def test_coverage_nonnegative(self, hybrid):
        """Coverage values should be between 0 and 1."""
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)
        assert 0.0 <= result.final_coverage <= 1.0
        assert 0.0 <= result.stage1_coverage <= 1.0

    def test_runtime_positive(self, hybrid):
        """Total runtime should be a positive number."""
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)
        assert result.total_runtime >= 0.0

    def test_stage1_has_more_primers_than_final(self, hybrid):
        """Stage 1 should select more primers than the final count to allow refinement."""
        result = hybrid.optimize(CANDIDATES, final_count=2, verbose=False)
        # Stage 1 should select at least as many as final (often more)
        assert len(result.stage1_primers) >= len(result.primers)

    def test_fixed_primers_included(self, hybrid):
        """Fixed primers should appear in the final set."""
        fixed = ["GCTAGCTA"]
        result = hybrid.optimize(
            CANDIDATES, final_count=3, fixed_primers=fixed, verbose=False
        )
        assert "GCTAGCTA" in result.primers

    def test_fixed_primers_exceed_target(self, hybrid):
        """When fixed primers meet target, they should be returned directly."""
        fixed = ["ATCGATCG", "GCTAGCTA", "AAATTTGG"]
        result = hybrid.optimize(
            CANDIDATES, final_count=3, fixed_primers=fixed, verbose=False
        )
        assert set(fixed).issubset(set(result.primers))

    def test_str_representation(self, hybrid):
        """HybridResult __str__ should produce readable output."""
        result = hybrid.optimize(CANDIDATES, final_count=3, verbose=False)
        text = str(result)
        assert "Hybrid Optimization Result" in text
        assert "Coverage" in text

    def test_polymerase_preset_phi29(self, cache):
        """phi29 preset should use 70000 bp max extension."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            polymerase="phi29",
        )
        assert h.max_extension == 70000

    def test_polymerase_preset_equiphi29(self, cache):
        """equiphi29 preset should use 80000 bp max extension."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            polymerase="equiphi29",
        )
        assert h.max_extension == 80000

    def test_polymerase_preset_bst(self, cache):
        """bst preset should use 10000 bp max extension."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            polymerase="bst",
        )
        assert h.max_extension == 10000

    def test_explicit_max_extension_overrides_preset(self, cache):
        """Explicitly setting max_extension should override the preset."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            polymerase="bst",
            max_extension=50000,
        )
        assert h.max_extension == 50000


# =============================================================================
# Edge Cases
# =============================================================================

class TestHybridEdgeCases:
    """Edge cases for hybrid optimization."""

    def test_single_candidate(self):
        """Optimization with a single candidate should return it."""
        cache = _make_position_cache({
            "ATCGATCG": np.array([5000, 50000]),
        })
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        result = h.optimize(["ATCGATCG"], final_count=5, verbose=False)
        assert result.primers == ["ATCGATCG"]

    def test_candidates_with_no_positions(self):
        """Primers with no binding positions should still not cause errors."""
        cache = _make_position_cache({})  # No positions for any primer
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
        # Should not raise, even if result is poor
        result = optimizer.optimize(CANDIDATES, target_size=3)
        assert result is not None

    def test_large_target_size(self):
        """Target size larger than candidate pool should return all candidates."""
        cache = _make_position_cache(POSITION_MAP)
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        result = h.optimize(CANDIDATES, final_count=100, verbose=False)
        # Should return at most len(CANDIDATES) primers
        assert len(result.primers) <= len(CANDIDATES)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
