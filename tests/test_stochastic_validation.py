"""
Tests for Phase 5: Stochastic Validation Integration.

Verifies that the stochastic simulation (Gillespie algorithm) is properly
integrated for validating network-based predictions.
"""

import pytest
import numpy as np
from unittest.mock import Mock, patch


class TestGillespieSimulator:
    """Test GillespieSimulator class."""

    def test_simulator_initializes(self):
        """Simulator should initialize without errors."""
        from neoswga.core.stochastic_simulator import GillespieSimulator, ReactionParameters

        # Create mock networks
        fg_network = Mock()
        fg_network.binding_sites = [100, 500, 1000, 1500]
        fg_network.graph = Mock()
        fg_network.graph.edges = Mock(return_value=[])

        bg_network = Mock()
        bg_network.binding_sites = [200, 600]
        bg_network.graph = Mock()
        bg_network.graph.edges = Mock(return_value=[])

        primers = ['ATCGATCG', 'GCTAGCTA']
        params = ReactionParameters(temperature=30.0)

        simulator = GillespieSimulator(primers, fg_network, bg_network, params)

        assert simulator is not None
        assert simulator.primers == primers

    def test_reaction_parameters_defaults(self):
        """ReactionParameters should have sensible defaults."""
        from neoswga.core.stochastic_simulator import ReactionParameters

        params = ReactionParameters()

        assert params.k_on > 0
        assert params.k_off > 0
        assert params.extension_rate > 0
        assert params.processivity_step > 0 and params.processivity_step <= 1
        assert params.max_extension > 0

    def test_extension_probability_decreases_with_distance(self):
        """Extension probability should decrease with distance."""
        from neoswga.core.stochastic_simulator import ReactionParameters

        params = ReactionParameters(processivity_step=0.95)

        prob_short = params.extension_probability(1000)
        prob_long = params.extension_probability(10000)

        assert prob_long < prob_short


class TestReactionState:
    """Test ReactionState class."""

    def test_state_initializes_with_defaults(self):
        """ReactionState should initialize with sensible defaults."""
        from neoswga.core.stochastic_simulator import ReactionState

        state = ReactionState()

        assert state.time == 0.0
        assert state.fg_molecules > 0
        assert state.bg_molecules > 0
        assert state.fg_products == 0
        assert state.bg_products == 0
        assert state.primer_conc > 0
        assert state.dNTP > 0


class TestValidateNetworkPredictions:
    """Test validate_network_predictions function."""

    def test_function_exists(self):
        """validate_network_predictions should be importable."""
        from neoswga.core.stochastic_simulator import validate_network_predictions
        assert callable(validate_network_predictions)


class TestStochasticValidationCLI:
    """Test CLI options for stochastic validation."""

    def test_cli_has_validate_simulation_option(self):
        """CLI should have --validate-simulation option."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['optimize', '-j', 'test.json', '--validate-simulation'])

        assert hasattr(args, 'validate_simulation')
        assert args.validate_simulation == True

    def test_cli_has_simulation_time_option(self):
        """CLI should have --simulation-time option."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'optimize', '-j', 'test.json',
            '--simulation-time', '7200'
        ])

        assert hasattr(args, 'simulation_time')
        assert args.simulation_time == 7200.0

    def test_cli_simulation_time_default(self):
        """Default --simulation-time should be 3600 seconds."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['optimize', '-j', 'test.json'])

        assert args.simulation_time == 3600.0


class TestSimulationPlotting:
    """Test simulation result plotting functions."""

    def test_plot_function_exists(self):
        """plot_simulation_results should be importable."""
        from neoswga.core.stochastic_simulator import plot_simulation_results
        assert callable(plot_simulation_results)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
