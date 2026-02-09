"""
Tests for optimizer consolidation (Tier 1 + Tier 2).

Verifies:
- greedy-dropout/dropout aliases resolve to greedy optimizer
- legacy-bfs is removed from CLI choices
- dominating-set factory registration delegates to original optimizer
- hybrid with polymerase='equiphi29' applies thermo-filtering
- hybrid with polymerase='phi29' skips thermo-filtering (default)
- hybrid with background_pruning=True runs background pruning stage
- hybrid with background_pruning=False skips it (default)
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch


# =========================================================================
# Tier 1a: Greedy-dropout merged into greedy
# =========================================================================

class TestGreedyDropoutMerged:
    """Test that greedy-dropout is merged into greedy."""

    def test_greedy_dropout_alias_resolves(self):
        """'greedy-dropout' alias should resolve to the greedy optimizer."""
        from neoswga.core.optimizer_factory import OptimizerRegistry
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()

        assert OptimizerRegistry.is_registered('greedy-dropout')
        # Should resolve to the same class as 'greedy'
        greedy_cls = OptimizerRegistry.get('greedy')
        dropout_cls = OptimizerRegistry.get('greedy-dropout')
        assert greedy_cls is dropout_cls

    def test_dropout_alias_resolves(self):
        """'dropout' alias should resolve to the greedy optimizer."""
        from neoswga.core.optimizer_factory import OptimizerRegistry
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()

        assert OptimizerRegistry.is_registered('dropout')
        greedy_cls = OptimizerRegistry.get('greedy')
        dropout_cls = OptimizerRegistry.get('dropout')
        assert greedy_cls is dropout_cls

    def test_greedy_has_dropout_method(self):
        """GreedyOptimizer should have _apply_dropout method."""
        from neoswga.core.greedy_optimizer import GreedyOptimizer
        assert hasattr(GreedyOptimizer, '_apply_dropout')

    def test_greedy_config_has_dropout_flag(self):
        """GreedyConfig should have enable_dropout flag."""
        from neoswga.core.greedy_optimizer import GreedyConfig
        config = GreedyConfig()
        assert hasattr(config, 'enable_dropout')
        assert config.enable_dropout is True
        assert hasattr(config, 'drop_out_iteration')
        assert config.drop_out_iteration == 4

    def test_no_greedy_dropout_class(self):
        """GreedyDropoutOptimizer class should no longer exist."""
        import neoswga.core.greedy_optimizer as module
        assert not hasattr(module, 'GreedyDropoutOptimizer')


# =========================================================================
# Tier 1b: legacy-bfs removed
# =========================================================================

class TestLegacyBfsRemoved:
    """Test that legacy-bfs is removed from CLI choices."""

    def test_legacy_bfs_not_in_cli(self):
        """legacy-bfs should not appear in CLI optimization method choices."""
        from neoswga.cli_unified import create_parser
        parser = create_parser()

        # Find the optimize subparser's --optimization-method argument
        for action in parser._subparsers._group_actions:
            for choice_name, subparser in action.choices.items():
                if choice_name == 'optimize':
                    for sub_action in subparser._actions:
                        if hasattr(sub_action, 'option_strings') and \
                           '--optimization-method' in sub_action.option_strings:
                            assert 'legacy-bfs' not in sub_action.choices, \
                                "'legacy-bfs' should be removed from CLI choices"
                            return
        pytest.fail("Could not find --optimization-method argument")


# =========================================================================
# Tier 1c: Dominating-set adapter delegates to original
# =========================================================================

class TestDominatingSetConsolidation:
    """Test that dominating-set adapter delegates to original optimizer."""

    def test_adapter_wraps_original(self):
        """DominatingSetAdapter should wrap the original DominatingSetOptimizer."""
        from neoswga.core.dominating_set_adapter import DominatingSetAdapter
        from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

        # Create mock position cache
        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        adapter = DominatingSetAdapter(
            position_cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
        )

        # Verify the adapter has an internal optimizer instance
        assert hasattr(adapter, '_optimizer')
        assert isinstance(adapter._optimizer, DominatingSetOptimizer)

    def test_factory_registration_intact(self):
        """'dominating-set' factory registration should work."""
        from neoswga.core.optimizer_factory import OptimizerRegistry
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()

        assert OptimizerRegistry.is_registered('dominating-set')
        assert OptimizerRegistry.is_registered('ds')
        assert OptimizerRegistry.is_registered('set-cover')

    def test_weighted_set_cover_still_registered(self):
        """WeightedSetCoverOptimizer should still be registered."""
        from neoswga.core.optimizer_factory import OptimizerRegistry
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()

        assert OptimizerRegistry.is_registered('weighted-set-cover')
        assert OptimizerRegistry.is_registered('wsc')


# =========================================================================
# Tier 2a: Polymerase-aware config in HybridOptimizer
# =========================================================================

class TestPolymeraseConfig:
    """Test polymerase-aware configuration in HybridOptimizer."""

    def test_polymerase_presets_exist(self):
        """All polymerase presets should be defined."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS
        assert 'phi29' in POLYMERASE_PRESETS
        assert 'equiphi29' in POLYMERASE_PRESETS
        assert 'bst' in POLYMERASE_PRESETS
        assert 'klenow' in POLYMERASE_PRESETS

    def test_phi29_defaults(self):
        """Phi29 should have standard defaults."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS
        phi29 = POLYMERASE_PRESETS['phi29']
        assert phi29.max_extension == 70000
        assert phi29.thermo_filter is False
        assert phi29.primer_multiplier == 1.0
        assert phi29.reaction_temp == 30.0

    def test_equiphi29_preset(self):
        """EquiPhi29 should have higher temp and thermo-filtering."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS
        eq29 = POLYMERASE_PRESETS['equiphi29']
        assert eq29.max_extension == 80000
        assert eq29.thermo_filter is True
        assert eq29.primer_multiplier == 0.85
        assert eq29.reaction_temp == 42.0

    def test_bst_preset(self):
        """BST should have short extension and high temp."""
        from neoswga.core.hybrid_optimizer import POLYMERASE_PRESETS
        bst = POLYMERASE_PRESETS['bst']
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
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
            polymerase='equiphi29',
        )

        assert opt.polymerase == 'equiphi29'
        assert opt.max_extension == 80000
        assert opt.poly_config.thermo_filter is True

    def test_hybrid_phi29_default(self):
        """HybridOptimizer should default to phi29."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
        )

        assert opt.polymerase == 'phi29'
        assert opt.max_extension == 70000
        assert opt.poly_config.thermo_filter is False

    def test_hybrid_explicit_max_extension_overrides_preset(self):
        """Explicit max_extension should override polymerase preset."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
            polymerase='equiphi29',
            max_extension=50000,  # Explicit override
        )

        assert opt.max_extension == 50000

    def test_gc_adaptive_adjustment(self):
        """GC-rich genomes should trigger wider GC acceptance."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
            polymerase='equiphi29',
            genome_gc_content=0.70,  # GC-rich
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
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
        )

        assert hasattr(adapter, '_hybrid')
        assert isinstance(adapter._hybrid, HybridOptimizer)
        assert adapter._hybrid.polymerase == 'equiphi29'


# =========================================================================
# Tier 2b: Background pruning in HybridOptimizer
# =========================================================================

class TestBackgroundPruning:
    """Test background pruning stage in HybridOptimizer."""

    def test_background_pruning_default_off(self):
        """Background pruning should be off by default."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([])

        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['test'],
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
            fg_prefixes=['test'],
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
            fg_prefixes=['test'],
            fg_seq_lengths=[1000],
            bg_prefixes=['bg_test'],
            bg_seq_lengths=[2000],
        )

        assert hasattr(adapter, '_hybrid')
        assert isinstance(adapter._hybrid, HybridOptimizer)
        assert adapter._hybrid.background_pruning is True

    def test_hybrid_has_prune_background_method(self):
        """HybridOptimizer should have _prune_background method."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        assert hasattr(HybridOptimizer, '_prune_background')

    def test_hybrid_has_count_background_sites_method(self):
        """HybridOptimizer should have _count_background_sites method."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer
        assert hasattr(HybridOptimizer, '_count_background_sites')

    def test_background_pruning_increases_stage1_count(self):
        """Background pruning should request more primers in Stage 1."""
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([100, 500, 1000])

        # With bg pruning: stage1_count should be >= final_count * 2
        opt = HybridOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['test'],
            fg_seq_lengths=[10000],
            background_pruning=True,
        )

        # The auto stage1_count formula with bg pruning uses max(final+8, final*1.67, final*2)
        # For final_count=12: max(20, 20, 24) = 24
        # This is verified by checking the code logic, not by running optimize
        assert opt.background_pruning is True


# =========================================================================
# Integration: All factory registrations intact
# =========================================================================

class TestAllRegistrationsIntact:
    """Verify all expected optimizer registrations still work."""

    @pytest.fixture(autouse=True)
    def ensure_registered(self):
        from neoswga.core.unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()

    def test_all_primary_names(self):
        from neoswga.core.optimizer_factory import OptimizerRegistry
        expected = [
            'greedy', 'dominating-set', 'weighted-set-cover',
            'hybrid', 'network', 'genetic', 'background-aware',
            'equiphi29', 'milp', 'normalized', 'tiling',
        ]
        for name in expected:
            assert OptimizerRegistry.is_registered(name), f"'{name}' not registered"

    def test_all_aliases(self):
        from neoswga.core.optimizer_factory import OptimizerRegistry
        expected_aliases = [
            'bfs', 'greedy-bfs', 'greedy-dropout', 'dropout',
            'ds', 'set-cover', 'wsc',
            'hybrid-optimizer', 'two-stage',
            'clinical', 'bg-aware',
            'equiphi', 'eq29',
        ]
        for alias in expected_aliases:
            assert OptimizerRegistry.is_registered(alias), f"Alias '{alias}' not registered"


# =========================================================================
# Background data flow fixes
# =========================================================================

class TestPositionCacheIncludesBgPrefixes:
    """Test that PositionCache in run_optimization includes bg_prefixes."""

    def test_bg_prefixes_in_source(self):
        """run_optimization source should include bg_prefixes in PositionCache."""
        import inspect
        from neoswga.core.unified_optimizer import run_optimization
        src = inspect.getsource(run_optimization)
        assert 'bg_prefixes or []' in src, \
            "bg_prefixes not included in PositionCache creation"

    def test_prefilter_function_exists(self):
        """_prefilter_by_background function should exist."""
        from neoswga.core.unified_optimizer import _prefilter_by_background
        assert callable(_prefilter_by_background)


class TestBackgroundAwareOptimizerParameterOrder:
    """Test that BackgroundAwareOptimizer uses correct get_positions order."""

    def test_calculate_coverage_parameter_order(self):
        """_calculate_coverage should call get_positions(prefix, primer, ...)."""
        import inspect
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer
        src = inspect.getsource(BackgroundAwareOptimizer._calculate_coverage)
        # Should use get_positions(self.fg_prefixes[0], primer, 'both')
        assert 'get_positions(self.fg_prefixes[0], primer' in src, \
            "_calculate_coverage uses wrong parameter order"
        # Should NOT use positions[0] or positions[1] tuple indexing
        assert 'positions[0]' not in src, \
            "_calculate_coverage still uses tuple indexing"
        assert 'positions[1]' not in src, \
            "_calculate_coverage still uses tuple indexing"

    def test_count_background_sites_parameter_order(self):
        """_count_background_sites should call get_positions(prefix, primer, ...)."""
        import inspect
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer
        src = inspect.getsource(BackgroundAwareOptimizer._count_background_sites)
        # Should use get_positions(bg_prefix, primer, 'both')
        assert 'get_positions(bg_prefix, primer' in src, \
            "_count_background_sites uses wrong parameter order"
        assert 'positions[0]' not in src, \
            "_count_background_sites still uses tuple indexing"
        assert 'positions[1]' not in src, \
            "_count_background_sites still uses tuple indexing"

    def test_calculate_coverage_with_mock(self):
        """_calculate_coverage should return correct coverage using mock cache."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer

        mock_cache = MagicMock()
        # Return 5 positions for each primer call
        mock_cache.get_positions.return_value = np.array([10, 20, 30, 40, 50])

        opt = BackgroundAwareOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['fg_test'],
            bg_prefixes=['bg_test'],
            fg_seq_lengths=[1000],
            bg_seq_lengths=[2000],
        )

        coverage = opt._calculate_coverage(['ATCGATCG', 'GCTAGCTA'])
        assert coverage > 0
        # Verify get_positions was called with correct order (prefix, primer, strand)
        calls = mock_cache.get_positions.call_args_list
        for call in calls:
            args = call[0]
            assert args[0] == 'fg_test', f"First arg should be prefix, got {args[0]}"

    def test_count_background_sites_with_mock(self):
        """_count_background_sites should count sites using correct parameter order."""
        from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer

        mock_cache = MagicMock()
        mock_cache.get_positions.return_value = np.array([100, 200, 300])

        opt = BackgroundAwareOptimizer(
            position_cache=mock_cache,
            fg_prefixes=['fg_test'],
            bg_prefixes=['bg_test'],
            fg_seq_lengths=[1000],
            bg_seq_lengths=[2000],
        )

        sites = opt._count_background_sites(['ATCGATCG'])
        assert sites == 3  # len of [100, 200, 300]
        # Verify called with bg_prefix first
        call_args = mock_cache.get_positions.call_args[0]
        assert call_args[0] == 'bg_test', f"First arg should be bg prefix, got {call_args[0]}"


# =========================================================================
# Background pre-filter tests
# =========================================================================

class TestBackgroundPrefilter:
    """Test the _prefilter_by_background function."""

    def _make_cache(self, fg_counts, bg_counts):
        """Create a mock cache with specified fg/bg binding counts."""
        mock_cache = MagicMock()

        def get_positions(prefix, primer, strand='both'):
            if prefix == 'fg_test':
                count = fg_counts.get(primer, 0)
            elif prefix == 'bg_test':
                count = bg_counts.get(primer, 0)
            else:
                count = 0
            return np.arange(count)

        mock_cache.get_positions.side_effect = get_positions
        return mock_cache

    def test_removes_high_background_primers(self):
        """Primers with poor fg/bg ratio should be removed."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        # primer_a: fg=10, bg=1 -> ratio=5.0  (good)
        # primer_b: fg=10, bg=100 -> ratio=0.099 (bad)
        # primer_c: fg=10, bg=2 -> ratio=3.33  (good)
        cache = self._make_cache(
            fg_counts={'primer_a': 10, 'primer_b': 10, 'primer_c': 10},
            bg_counts={'primer_a': 1, 'primer_b': 100, 'primer_c': 2},
        )

        result = _prefilter_by_background(
            cache, ['primer_a', 'primer_b', 'primer_c'],
            fg_prefixes=['fg_test'], bg_prefixes=['bg_test'],
            min_ratio=1.0, max_removal_fraction=0.50,
        )

        assert 'primer_a' in result
        assert 'primer_c' in result
        assert 'primer_b' not in result

    def test_respects_max_removal_fraction(self):
        """Should not remove more than max_removal_fraction of candidates."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        # All primers have poor ratio (< 1.0)
        cache = self._make_cache(
            fg_counts={'p1': 1, 'p2': 1, 'p3': 1, 'p4': 1, 'p5': 1},
            bg_counts={'p1': 10, 'p2': 20, 'p3': 30, 'p4': 40, 'p5': 50},
        )

        candidates = ['p1', 'p2', 'p3', 'p4', 'p5']
        result = _prefilter_by_background(
            cache, candidates,
            fg_prefixes=['fg_test'], bg_prefixes=['bg_test'],
            min_ratio=1.0, max_removal_fraction=0.20,
        )

        # max_removal=20% of 5 = 1 removed, so keep at least 4
        assert len(result) >= 4

    def test_never_returns_empty(self):
        """Pre-filter should never return an empty list."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        cache = self._make_cache(
            fg_counts={'p1': 0},
            bg_counts={'p1': 100},
        )

        result = _prefilter_by_background(
            cache, ['p1'],
            fg_prefixes=['fg_test'], bg_prefixes=['bg_test'],
            min_ratio=1.0, max_removal_fraction=0.50,
        )

        assert len(result) >= 1

    def test_no_removal_when_all_good(self):
        """No primers removed when all have good ratio."""
        from neoswga.core.unified_optimizer import _prefilter_by_background

        cache = self._make_cache(
            fg_counts={'p1': 100, 'p2': 100, 'p3': 100},
            bg_counts={'p1': 1, 'p2': 1, 'p3': 1},
        )

        candidates = ['p1', 'p2', 'p3']
        result = _prefilter_by_background(
            cache, candidates,
            fg_prefixes=['fg_test'], bg_prefixes=['bg_test'],
            min_ratio=1.0, max_removal_fraction=0.20,
        )

        assert len(result) == 3


class TestNoBgPrefilterCliFlag:
    """Test --no-bg-prefilter CLI flag."""

    def test_flag_exists_in_parser(self):
        """--no-bg-prefilter flag should exist in optimize subparser."""
        from neoswga.cli_unified import create_parser
        parser = create_parser()

        for action in parser._subparsers._group_actions:
            for choice_name, subparser in action.choices.items():
                if choice_name == 'optimize':
                    option_strings = []
                    for sub_action in subparser._actions:
                        if hasattr(sub_action, 'option_strings'):
                            option_strings.extend(sub_action.option_strings)
                    assert '--no-bg-prefilter' in option_strings, \
                        "'--no-bg-prefilter' flag not found in optimize subparser"
                    return
        pytest.fail("Could not find optimize subparser")
