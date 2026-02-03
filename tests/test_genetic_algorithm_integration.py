"""
Tests for Phase 3: Genetic Algorithm Integration.

Verifies that the genetic algorithm properly evaluates individuals using
actual position data instead of random placeholders.
"""

import pytest
import numpy as np
from unittest.mock import Mock, MagicMock


class TestGeneticAlgorithmEvaluation:
    """Test GA evaluation with position cache."""

    def test_ga_accepts_position_cache(self):
        """PrimerSetGA should accept position_cache parameter."""
        from neoswga.core.genetic_algorithm import PrimerSetGA, GAConfig
        from neoswga.core.reaction_conditions import ReactionConditions

        # Create mock cache
        mock_cache = Mock()
        mock_cache.get_positions = Mock(return_value=np.array([]))

        # Create conditions
        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        config = GAConfig(
            population_size=10,
            generations=2
        )

        # Should not raise error
        ga = PrimerSetGA(
            primer_pool=['ATCGATCG', 'GCTAGCTA'],
            fg_prefixes=['test_fg'],
            bg_prefixes=['test_bg'],
            fg_lengths=[100000],
            bg_lengths=[100000],
            conditions=conditions,
            config=config,
            position_cache=mock_cache
        )

        assert ga.position_cache is mock_cache

    def test_evaluate_individual_uses_cache(self):
        """_evaluate_individual should use position cache when available."""
        from neoswga.core.genetic_algorithm import PrimerSetGA, GAConfig, Individual
        from neoswga.core.reaction_conditions import ReactionConditions

        # Create mock cache that returns positions
        mock_cache = Mock()
        def mock_get_positions(prefix, primer, strand):
            # Return some mock positions
            if 'fg' in prefix:
                return np.array([1000, 2000, 3000, 4000, 5000])
            else:
                return np.array([1000])  # Less in background

        mock_cache.get_positions = Mock(side_effect=mock_get_positions)

        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        config = GAConfig(
            population_size=10,
            generations=2
        )

        ga = PrimerSetGA(
            primer_pool=['ATCGATCG', 'GCTAGCTA'],
            fg_prefixes=['test_fg'],
            bg_prefixes=['test_bg'],
            fg_lengths=[100000],
            bg_lengths=[100000],
            conditions=conditions,
            config=config,
            position_cache=mock_cache
        )

        # Evaluate an individual
        individual = Individual(primers=['ATCGATCG'])
        fitness, metrics = ga._evaluate_individual(individual)

        # Should have called get_positions
        assert mock_cache.get_positions.called

        # Metrics should include fg_sites and bg_sites
        assert 'fg_sites' in metrics
        assert 'bg_sites' in metrics
        assert metrics['fg_sites'] > 0

    def test_evaluate_individual_without_cache_uses_fallback(self):
        """Without position cache, should use fallback (not crash)."""
        from neoswga.core.genetic_algorithm import PrimerSetGA, GAConfig, Individual
        from neoswga.core.reaction_conditions import ReactionConditions

        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        config = GAConfig(
            population_size=10,
            generations=2
        )

        ga = PrimerSetGA(
            primer_pool=['ATCGATCG', 'GCTAGCTA'],
            fg_prefixes=['test_fg'],
            bg_prefixes=['test_bg'],
            fg_lengths=[100000],
            bg_lengths=[100000],
            conditions=conditions,
            config=config,
            position_cache=None  # No cache
        )

        # Should not crash without cache
        individual = Individual(primers=['ATCGATCG'])
        fitness, metrics = ga._evaluate_individual(individual)

        # Should still return valid fitness
        assert 0 <= fitness <= 1
        assert 'coverage' in metrics


class TestGeneticAlgorithmCLI:
    """Test CLI options for genetic algorithm."""

    def test_cli_has_genetic_option(self):
        """CLI should have 'genetic' as optimization method."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['optimize', '-j', 'test.json', '--optimization-method', 'genetic'])

        assert args.optimization_method == 'genetic'


class TestGeneticAlgorithmGini:
    """Test Gini calculation."""

    def test_gini_calculation(self):
        """Gini should be calculated correctly."""
        from neoswga.core.genetic_algorithm import PrimerSetGA, GAConfig
        from neoswga.core.reaction_conditions import ReactionConditions

        conditions = ReactionConditions(temp=30.0, na_conc=50.0)
        config = GAConfig(population_size=10, generations=2)

        ga = PrimerSetGA(
            primer_pool=['ATCGATCG'],
            fg_prefixes=['test_fg'],
            bg_prefixes=['test_bg'],
            fg_lengths=[100000],
            bg_lengths=[100000],
            conditions=conditions,
            config=config
        )

        # Uniform distribution should have low Gini
        uniform = np.array([100, 100, 100, 100, 100])
        gini_uniform = ga._calculate_gini(uniform)

        # Unequal distribution should have high Gini
        unequal = np.array([1, 1, 1, 1, 1000])
        gini_unequal = ga._calculate_gini(unequal)

        assert gini_uniform < gini_unequal
        assert 0 <= gini_uniform <= 1
        assert 0 <= gini_unequal <= 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
