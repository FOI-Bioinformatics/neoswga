"""
Consolidated tests for optimizer factory, registry, and hybrid optimizer.

Covers:
- OptimizerRegistry: registration, lookup, aliases, error handling
- OptimizerFactory: creation, validation, from_params
- Greedy-dropout merger and legacy-bfs removal
- Dominating-set adapter delegation
- Polymerase-aware configuration in HybridOptimizer
- Background pruning in HybridOptimizer
- Background data flow and pre-filter
- HybridBaseOptimizer: instantiation via factory and end-to-end optimization
- HybridOptimizer: direct two-stage optimization with synthetic data
- CLI integration checks
"""

import inspect

import numpy as np
import pytest
from unittest.mock import Mock, MagicMock

from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry
from neoswga.core.base_optimizer import BaseOptimizer, OptimizerConfig, OptimizationStatus
from neoswga.core.exceptions import OptimizerNotFoundError

# Import optimizer modules to trigger @OptimizerFactory.register() decorators.
import neoswga.core.hybrid_optimizer  # noqa: F401
import neoswga.core.greedy_optimizer  # noqa: F401
import neoswga.core.dominating_set_adapter  # noqa: F401
import neoswga.core.network_optimizer  # noqa: F401
import neoswga.core.genetic_algorithm  # noqa: F401

# Ensure unified_optimizer registrations are loaded (aliases, etc.)
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


# Candidate primers for testing.
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

    def test_all_primary_names(self):
        """All expected primary optimizer names should be registered."""
        expected = [
            "greedy",
            "dominating-set",
            "weighted-set-cover",
            "hybrid",
            "network",
            "genetic",
            "background-aware",
            "equiphi29",
            "milp",
            "normalized",
            "tiling",
        ]
        for name in expected:
            assert OptimizerRegistry.is_registered(name), f"'{name}' not registered"

    def test_all_aliases(self):
        """All expected optimizer aliases should be registered."""
        expected_aliases = [
            "bfs",
            "greedy-bfs",
            "greedy-dropout",
            "dropout",
            "ds",
            "set-cover",
            "wsc",
            "hybrid-optimizer",
            "two-stage",
            "clinical",
            "bg-aware",
            "equiphi",
            "eq29",
        ]
        for alias in expected_aliases:
            assert OptimizerRegistry.is_registered(alias), f"Alias '{alias}' not registered"


# =============================================================================
# Greedy-dropout merger
# =============================================================================


class TestGreedyDropoutMerged:
    """Test that greedy-dropout is merged into greedy."""

    def test_greedy_dropout_alias_resolves(self):
        """'greedy-dropout' alias should resolve to the greedy optimizer."""
        assert OptimizerRegistry.is_registered("greedy-dropout")
        greedy_cls = OptimizerRegistry.get("greedy")
        dropout_cls = OptimizerRegistry.get("greedy-dropout")
        assert greedy_cls is dropout_cls

    def test_dropout_alias_resolves(self):
        """'dropout' alias should resolve to the greedy optimizer."""
        assert OptimizerRegistry.is_registered("dropout")
        greedy_cls = OptimizerRegistry.get("greedy")
        dropout_cls = OptimizerRegistry.get("dropout")
        assert greedy_cls is dropout_cls

    def test_greedy_has_dropout_method(self):
        """GreedyOptimizer should have _apply_dropout method."""
        from neoswga.core.greedy_optimizer import GreedyOptimizer

        assert hasattr(GreedyOptimizer, "_apply_dropout")

    def test_greedy_config_has_dropout_flag(self):
        """GreedyConfig should have enable_dropout flag."""
        from neoswga.core.greedy_optimizer import GreedyConfig

        config = GreedyConfig()
        assert hasattr(config, "enable_dropout")
        assert config.enable_dropout is True
        assert hasattr(config, "drop_out_iteration")
        assert config.drop_out_iteration == 4

    def test_no_greedy_dropout_class(self):
        """GreedyDropoutOptimizer class should no longer exist."""
        import neoswga.core.greedy_optimizer as module

        assert not hasattr(module, "GreedyDropoutOptimizer")


# =============================================================================
# Legacy-bfs removed from CLI
# =============================================================================


class TestLegacyBfsRemoved:
    """Test that legacy-bfs is removed from CLI choices."""

    def test_legacy_bfs_not_in_cli(self):
        """legacy-bfs should not appear in CLI optimization method choices."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()

        for action in parser._subparsers._group_actions:
            for choice_name, subparser in action.choices.items():
                if choice_name == "optimize":
                    for sub_action in subparser._actions:
                        if (
                            hasattr(sub_action, "option_strings")
                            and "--optimization-method" in sub_action.option_strings
                        ):
                            assert "legacy-bfs" not in sub_action.choices, (
                                "'legacy-bfs' should be removed from CLI choices"
                            )
                            return
        pytest.fail("Could not find --optimization-method argument")


# =============================================================================
# Dominating-set adapter delegation
# =============================================================================


class TestDominatingSetConsolidation:
    """Test that dominating-set adapter delegates to original optimizer."""

    def test_adapter_wraps_original(self):
        """DominatingSetAdapter should wrap the original DominatingSetOptimizer."""
        from neoswga.core.dominating_set_adapter import DominatingSetAdapter
        from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        adapter = DominatingSetAdapter(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
        )

        assert hasattr(adapter, "_optimizer")
        assert isinstance(adapter._optimizer, DominatingSetOptimizer)

    def test_factory_registration_intact(self):
        """'dominating-set' factory registration should work."""
        assert OptimizerRegistry.is_registered("dominating-set")
        assert OptimizerRegistry.is_registered("ds")
        assert OptimizerRegistry.is_registered("set-cover")

    def test_weighted_set_cover_still_registered(self):
        """WeightedSetCoverOptimizer should still be registered."""
        assert OptimizerRegistry.is_registered("weighted-set-cover")
        assert OptimizerRegistry.is_registered("wsc")


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
# Polymerase-aware configuration in HybridOptimizer
# =============================================================================


class TestPolymeraseConfig:
    """Test polymerase-aware configuration in HybridOptimizer."""

    def test_polymerase_presets_exist(self):
        """All polymerase presets should be defined."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS

        assert "phi29" in POLYMERASE_PRESETS
        assert "equiphi29" in POLYMERASE_PRESETS
        assert "bst" in POLYMERASE_PRESETS
        assert "klenow" in POLYMERASE_PRESETS

    def test_phi29_defaults(self):
        """Phi29 should have standard defaults."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS

        phi29 = POLYMERASE_PRESETS["phi29"]
        assert phi29.max_extension == 70000
        assert phi29.thermo_filter is False
        assert phi29.primer_multiplier == 1.0
        assert phi29.reaction_temp == 30.0

    def test_equiphi29_preset(self):
        """EquiPhi29 should have higher temp and thermo-filtering."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS

        eq29 = POLYMERASE_PRESETS["equiphi29"]
        assert eq29.max_extension == 80000
        assert eq29.thermo_filter is True
        assert eq29.primer_multiplier == 0.85
        assert eq29.reaction_temp == 42.0

    def test_bst_preset(self):
        """BST should have short extension and high temp."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS

        bst = POLYMERASE_PRESETS["bst"]
        assert bst.max_extension == 10000
        assert bst.thermo_filter is True
        assert bst.reaction_temp == 60.0

    def test_hybrid_accepts_polymerase(self):
        """HybridOptimizer should accept polymerase parameter."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
            polymerase="equiphi29",
        )

        assert opt.polymerase == "equiphi29"
        assert opt.max_extension == 80000
        assert opt.poly_config.thermo_filter is True

    def test_hybrid_phi29_default(self):
        """HybridOptimizer should default to phi29."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
        )

        assert opt.polymerase == "phi29"
        assert opt.max_extension == 70000
        assert opt.poly_config.thermo_filter is False

    def test_hybrid_explicit_max_extension_overrides_preset(self):
        """Explicit max_extension should override polymerase preset."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
            polymerase="equiphi29",
            max_extension=50000,
        )

        assert opt.max_extension == 50000

    def test_gc_adaptive_adjustment(self):
        """GC-rich genomes should trigger wider GC acceptance."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
            polymerase="equiphi29",
            genome_gc_content=0.70,
        )

        assert opt.poly_config.max_gc == 0.80

    def test_equiphi29_adapter_delegates_to_hybrid(self):
        """EquiPhi29OptimizerAdapter should use HybridOptimizer internally."""
        from neoswga.core.equiphi29_optimizer import EquiPhi29OptimizerAdapter
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        adapter = EquiPhi29OptimizerAdapter(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
        )

        assert hasattr(adapter, "_hybrid")
        assert isinstance(adapter._hybrid, HybridOptimizer)
        assert adapter._hybrid.polymerase == "equiphi29"


# =============================================================================
# Background pruning in HybridOptimizer
# =============================================================================


class TestBackgroundPruning:
    """Test background pruning stage in HybridOptimizer."""

    def test_background_pruning_default_off(self):
        """Background pruning should be off by default."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
        )

        assert opt.background_pruning is False

    def test_background_pruning_can_be_enabled(self):
        """Background pruning should be configurable."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
            background_pruning=True,
            background_weight=3.0,
            min_coverage_threshold=0.90,
        )

        assert opt.background_pruning is True
        assert opt.background_weight == 3.0
        assert opt.min_coverage_threshold == 0.90

    def test_background_aware_adapter_uses_hybrid(self):
        """BackgroundAwareBaseOptimizer should use HybridOptimizer with bg pruning."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareBaseOptimizer
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        adapter = BackgroundAwareBaseOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[1000],
            bg_prefixes=["bg_test"],
            bg_seq_lengths=[2000],
        )

        assert hasattr(adapter, "_hybrid")
        assert isinstance(adapter._hybrid, HybridOptimizer)
        assert adapter._hybrid.background_pruning is True

    def test_hybrid_has_prune_background_method(self):
        """HybridOptimizer should have _prune_background method."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        assert hasattr(HybridOptimizer, "_prune_background")

    def test_hybrid_has_count_background_sites_method(self):
        """HybridOptimizer should have _count_background_sites method."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        assert hasattr(HybridOptimizer, "_count_background_sites")

    def test_background_pruning_flag_enables(self):
        """Background pruning should be enabled when set."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([100, 500, 1000])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["test"],
            fg_seq_lengths=[10000],
            background_pruning=True,
        )

        assert opt.background_pruning is True


# =============================================================================
# Background data flow fixes
# =============================================================================


class TestPositionCacheIncludesBgPrefixes:
    """Test that PositionCache in run_optimization includes bg_prefixes."""

    def test_bg_prefixes_in_source(self):
        """run_optimization source should include bg_prefixes in PositionCache."""
        from neoswga.core.unified_optimizer import run_optimization

        src = inspect.getsource(run_optimization)
        assert "bg_prefixes or []" in src, (
            "bg_prefixes not included in PositionCache creation"
        )

    def test_prefilter_function_exists(self):
        """_prefilter_by_background function should exist."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        assert callable(_prefilter_by_background)


class TestBackgroundAwareOptimizerParameterOrder:
    """Test that BackgroundAwareOptimizer uses correct get_positions order."""

    def test_calculate_coverage_parameter_order(self):
        """_calculate_coverage should iterate over all fg_prefixes."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer

        src = inspect.getsource(BackgroundAwareOptimizer._calculate_coverage)
        # Must loop over all prefixes, not just [0]
        assert "for prefix in self.fg_prefixes" in src, (
            "_calculate_coverage must iterate over all fg_prefixes"
        )
        assert "get_positions(prefix, primer" in src, (
            "_calculate_coverage uses wrong parameter order"
        )
        assert "positions[0]" not in src, (
            "_calculate_coverage still uses tuple indexing"
        )
        assert "positions[1]" not in src, (
            "_calculate_coverage still uses tuple indexing"
        )

    def test_count_background_sites_parameter_order(self):
        """_count_background_sites should call get_positions(prefix, primer, ...)."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer

        src = inspect.getsource(BackgroundAwareOptimizer._count_background_sites)
        assert "get_positions(bg_prefix, primer" in src, (
            "_count_background_sites uses wrong parameter order"
        )
        assert "positions[0]" not in src, (
            "_count_background_sites still uses tuple indexing"
        )
        assert "positions[1]" not in src, (
            "_count_background_sites still uses tuple indexing"
        )

    def test_calculate_coverage_with_mock(self):
        """_calculate_coverage should return correct coverage using mock cache."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([10, 20, 30, 40, 50])

        opt = BackgroundAwareOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["fg_test"],
            bg_prefixes=["bg_test"],
            fg_seq_lengths=[1000],
            bg_seq_lengths=[2000],
        )

        coverage = opt._calculate_coverage(["ATCGATCG", "GCTAGCTA"])
        assert coverage > 0
        calls = mock_cache.get_positions.call_args_list
        for call in calls:
            args = call[0]
            assert args[0] == "fg_test", f"First arg should be prefix, got {args[0]}"

    def test_count_background_sites_with_mock(self):
        """_count_background_sites should count sites using correct parameter order."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([100, 200, 300])

        opt = BackgroundAwareOptimizer(
            position_cache=mock_cache,
            fg_prefixes=["fg_test"],
            bg_prefixes=["bg_test"],
            fg_seq_lengths=[1000],
            bg_seq_lengths=[2000],
        )

        sites = opt._count_background_sites(["ATCGATCG"])
        assert sites == 3
        call_args = mock_cache.get_positions.call_args[0]
        assert call_args[0] == "bg_test", (
            f"First arg should be bg prefix, got {call_args[0]}"
        )


# =============================================================================
# Background pre-filter tests
# =============================================================================


class TestBackgroundPrefilter:
    """Test the _prefilter_by_background function."""

    def _make_cache(self, fg_counts, bg_counts):
        """Create a mock cache with specified fg/bg binding counts."""
        mock_cache = MagicMock()

        def get_positions(prefix, primer, strand="both"):
            if prefix == "fg_test":
                count = fg_counts.get(primer, 0)
            elif prefix == "bg_test":
                count = bg_counts.get(primer, 0)
            else:
                count = 0
            return np.arange(count)

        mock_cache.get_positions.side_effect = get_positions
        return mock_cache

    def test_removes_high_background_primers(self):
        """Primers with poor fg/bg ratio should be removed."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        cache = self._make_cache(
            fg_counts={"primer_a": 10, "primer_b": 10, "primer_c": 10},
            bg_counts={"primer_a": 1, "primer_b": 100, "primer_c": 2},
        )

        result = _prefilter_by_background(
            cache,
            ["primer_a", "primer_b", "primer_c"],
            fg_prefixes=["fg_test"],
            bg_prefixes=["bg_test"],
            min_ratio=1.0,
            max_removal_fraction=0.50,
        )

        assert "primer_a" in result
        assert "primer_c" in result
        assert "primer_b" not in result

    def test_respects_max_removal_fraction(self):
        """Should not remove more than max_removal_fraction of candidates."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        cache = self._make_cache(
            fg_counts={"p1": 1, "p2": 1, "p3": 1, "p4": 1, "p5": 1},
            bg_counts={"p1": 10, "p2": 20, "p3": 30, "p4": 40, "p5": 50},
        )

        candidates = ["p1", "p2", "p3", "p4", "p5"]
        result = _prefilter_by_background(
            cache,
            candidates,
            fg_prefixes=["fg_test"],
            bg_prefixes=["bg_test"],
            min_ratio=1.0,
            max_removal_fraction=0.20,
        )

        assert len(result) >= 4

    def test_never_returns_empty(self):
        """Pre-filter should never return an empty list."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        cache = self._make_cache(
            fg_counts={"p1": 0},
            bg_counts={"p1": 100},
        )

        result = _prefilter_by_background(
            cache,
            ["p1"],
            fg_prefixes=["fg_test"],
            bg_prefixes=["bg_test"],
            min_ratio=1.0,
            max_removal_fraction=0.50,
        )

        assert len(result) >= 1

    def test_no_removal_when_all_good(self):
        """No primers removed when all have good ratio."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        cache = self._make_cache(
            fg_counts={"p1": 100, "p2": 100, "p3": 100},
            bg_counts={"p1": 1, "p2": 1, "p3": 1},
        )

        candidates = ["p1", "p2", "p3"]
        result = _prefilter_by_background(
            cache,
            candidates,
            fg_prefixes=["fg_test"],
            bg_prefixes=["bg_test"],
            min_ratio=1.0,
            max_removal_fraction=0.20,
        )

        assert len(result) == 3


# =============================================================================
# CLI flag tests
# =============================================================================


class TestNoBgPrefilterCliFlag:
    """Test --no-bg-prefilter CLI flag."""

    def test_flag_exists_in_parser(self):
        """--no-bg-prefilter flag should exist in optimize subparser."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()

        for action in parser._subparsers._group_actions:
            for choice_name, subparser in action.choices.items():
                if choice_name == "optimize":
                    option_strings = []
                    for sub_action in subparser._actions:
                        if hasattr(sub_action, "option_strings"):
                            option_strings.extend(sub_action.option_strings)
                    assert "--no-bg-prefilter" in option_strings, (
                        "'--no-bg-prefilter' flag not found in optimize subparser"
                    )
                    return
        pytest.fail("Could not find optimize subparser")


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
        """Stage 1 should select at least as many primers as the final set."""
        result = hybrid.optimize(CANDIDATES, final_count=2, verbose=False)
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
        cache = _make_position_cache(
            {
                "ATCGATCG": np.array([5000, 50000]),
            }
        )
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        h = HybridOptimizer(
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
        )
        result = h.optimize(["ATCGATCG"], final_count=5, verbose=False)
        assert result.primers == ["ATCGATCG"]

    def test_candidates_with_no_positions(self):
        """Primers with no binding positions should not cause errors."""
        cache = _make_position_cache({})
        config = OptimizerConfig(target_set_size=3, verbose=False)
        optimizer = OptimizerFactory.create(
            "hybrid",
            position_cache=cache,
            fg_prefixes=FG_PREFIXES,
            fg_seq_lengths=FG_SEQ_LENGTHS,
            config=config,
        )
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
        assert len(result.primers) <= len(CANDIDATES)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
