"""
Unit tests for parameter module.

Tests:
- PipelineParameters dataclass
- get_current_config() function
- set_from_config() function
- reset_to_defaults() function
- Round-trip consistency
"""

import pytest
from dataclasses import fields

from neoswga.core.parameter import (
    PipelineParameters,
    get_current_config,
    set_from_config,
    reset_to_defaults,
)


class TestPipelineParameters:
    """Tests for PipelineParameters dataclass."""

    def test_default_values(self):
        """Test default parameter values."""
        params = PipelineParameters()

        assert params.min_k == 6
        assert params.max_k == 12
        assert params.min_fg_freq == 1e-5
        assert params.max_gini == 0.6
        assert params.polymerase == 'phi29'
        assert params.fg_circular is True
        assert params.verbose is False

    def test_custom_values(self):
        """Test custom parameter values."""
        params = PipelineParameters(
            min_k=8,
            max_k=15,
            polymerase='equiphi29',
            reaction_temp=42.0,
            verbose=True
        )

        assert params.min_k == 8
        assert params.max_k == 15
        assert params.polymerase == 'equiphi29'
        assert params.reaction_temp == 42.0
        assert params.verbose is True

    def test_list_fields_default(self):
        """Test that list fields default to empty lists."""
        params = PipelineParameters()

        assert params.fg_genomes == []
        assert params.bg_genomes == []
        assert params.fg_prefixes == []
        assert params.fg_seq_lengths == []

    def test_list_fields_custom(self):
        """Test setting custom list values."""
        params = PipelineParameters(
            fg_genomes=['/path/to/genome.fasta'],
            fg_prefixes=['/path/to/data/target']
        )

        assert params.fg_genomes == ['/path/to/genome.fasta']
        assert params.fg_prefixes == ['/path/to/data/target']

    def test_thermodynamic_additives(self):
        """Test thermodynamic additive parameters."""
        params = PipelineParameters(
            dmso_percent=5.0,
            betaine_m=1.0,
            mg_conc=3.5
        )

        assert params.dmso_percent == 5.0
        assert params.betaine_m == 1.0
        assert params.mg_conc == 3.5


class TestGetCurrentConfig:
    """Tests for get_current_config() function."""

    def test_returns_pipeline_parameters(self):
        """Test that get_current_config returns PipelineParameters."""
        config = get_current_config()

        assert isinstance(config, PipelineParameters)

    def test_has_expected_attributes(self):
        """Test that config has expected attributes."""
        config = get_current_config()

        # Check key attributes exist
        assert hasattr(config, 'min_k')
        assert hasattr(config, 'max_k')
        assert hasattr(config, 'polymerase')
        assert hasattr(config, 'reaction_temp')
        assert hasattr(config, 'fg_genomes')

    def test_returns_sensible_defaults(self):
        """Test that returned config has sensible default values."""
        config = get_current_config()

        # K-mer range should be valid
        assert config.min_k >= 4
        assert config.max_k >= config.min_k

        # GC range should be valid
        assert 0 <= config.gc_min <= 1
        assert 0 <= config.gc_max <= 1
        assert config.gc_min <= config.gc_max


# =============================================================================
# set_from_config Tests
# =============================================================================

class TestSetFromConfig:
    """Tests for set_from_config function."""

    def test_sets_global_state(self):
        """Test that set_from_config updates globals."""
        import neoswga.core.parameter as param

        original_min_k = param.min_k
        original_polymerase = param.polymerase

        try:
            config = PipelineParameters(
                min_k=15,
                max_k=25,
                polymerase='bst',
            )
            set_from_config(config)

            assert param.min_k == 15
            assert param.max_k == 25
            assert param.polymerase == 'bst'
        finally:
            # Restore
            param.min_k = original_min_k
            param.polymerase = original_polymerase

    def test_round_trip_consistency(self):
        """Test that get -> set -> get produces same values."""
        # Get current config
        config1 = get_current_config()

        # Modify and set
        config2 = PipelineParameters(
            min_k=8,
            max_k=16,
            min_fg_freq=1e-4,
            polymerase='equiphi29',
            reaction_temp=45.0,
        )
        set_from_config(config2)

        # Get again and compare
        config3 = get_current_config()

        assert config3.min_k == config2.min_k
        assert config3.max_k == config2.max_k
        assert config3.min_fg_freq == config2.min_fg_freq
        assert config3.polymerase == config2.polymerase
        assert config3.reaction_temp == config2.reaction_temp

        # Restore original
        set_from_config(config1)


# =============================================================================
# reset_to_defaults Tests
# =============================================================================

class TestResetToDefaults:
    """Tests for reset_to_defaults function."""

    def test_resets_values(self):
        """Test that reset_to_defaults restores default values."""
        import neoswga.core.parameter as param

        # Modify some values
        param.min_k = 99
        param.polymerase = 'custom'

        # Reset
        reset_to_defaults()

        # Check defaults are restored
        config = get_current_config()
        assert config.min_k == 6
        assert config.polymerase == 'phi29'

    def test_clears_lists(self):
        """Test that reset_to_defaults clears list fields."""
        import neoswga.core.parameter as param

        # Set some genome data
        param.fg_genomes = ['/path/1.fa', '/path/2.fa']
        param.fg_prefixes = ['/prefix1', '/prefix2']

        # Reset
        reset_to_defaults()

        # Check lists are empty
        config = get_current_config()
        assert config.fg_genomes == []
        assert config.fg_prefixes == []


# =============================================================================
# Field Coverage Tests
# =============================================================================

class TestFieldCoverage:
    """Tests to ensure all fields are handled by get/set."""

    def test_all_dataclass_fields_are_getable(self):
        """Test that all PipelineParameters fields are returned by get_current_config."""
        config = get_current_config()

        for f in fields(PipelineParameters):
            # All fields should be accessible
            assert hasattr(config, f.name), f"Field {f.name} not accessible"

    def test_all_dataclass_fields_are_setable(self):
        """Test that all PipelineParameters fields are set by set_from_config."""
        import neoswga.core.parameter as param

        # Create a config with specific values for each field
        config = PipelineParameters(
            min_k=10,
            max_k=20,
            min_fg_freq=2e-5,
            max_bg_freq=1e-5,
            min_tm=20.0,
            max_tm=50.0,
            max_gini=0.5,
            max_primer=1000,
            min_amp_pred=15.0,
            max_dimer_bp=4,
            max_self_dimer_bp=5,
            gc_min=0.4,
            gc_max=0.6,
            gc_tolerance=0.1,
            genome_gc=0.5,
            polymerase='equiphi29',
            reaction_temp=45.0,
            na_conc=100.0,
            mg_conc=3.0,
            primer_conc=1e-6,
            dmso_percent=5.0,
            betaine_m=1.0,
            iterations=10,
            max_sets=8,
            fg_circular=False,
            bg_circular=True,
            sample_rate=0.5,
            min_sample_count=10,
            use_bloom_filter=True,
            data_dir='/test/data',
            cpus=8,
            verbose=True,
        )

        set_from_config(config)

        # Verify each field was set
        assert param.min_k == 10
        assert param.max_k == 20
        assert param.polymerase == 'equiphi29'
        assert param.reaction_temp == 45.0
        assert param.dmso_percent == 5.0
        assert param.betaine_m == 1.0
        assert param.fg_circular is False
        assert param.bg_circular is True
        assert param.use_bloom_filter is True
        assert param.verbose is True

        # Clean up
        reset_to_defaults()


# =============================================================================
# Integration Tests
# =============================================================================

class TestParameterIntegration:
    """Integration tests for parameter management."""

    def test_pipeline_workflow(self):
        """Test typical pipeline configuration workflow."""
        # 1. Create a typed configuration
        config = PipelineParameters(
            min_k=6,
            max_k=12,
            polymerase='phi29',
            reaction_temp=30.0,
            min_fg_freq=1e-5,
            max_bg_freq=5e-6,
        )

        # 2. Apply to globals
        set_from_config(config)

        # 3. Modules can now access via getattr
        import neoswga.core.parameter as param
        assert getattr(param, 'polymerase') == 'phi29'
        assert getattr(param, 'reaction_temp') == 30.0

        # 4. Can also get typed config
        current = get_current_config()
        assert current.polymerase == 'phi29'

        # Clean up
        reset_to_defaults()

    def test_equiphi29_configuration(self):
        """Test EquiPhi29 configuration pattern."""
        config = PipelineParameters(
            min_k=12,
            max_k=18,
            polymerase='equiphi29',
            reaction_temp=45.0,
            min_tm=35.0,
            max_tm=65.0,
            dmso_percent=5.0,
            betaine_m=1.0,
        )

        set_from_config(config)

        import neoswga.core.parameter as param
        assert param.polymerase == 'equiphi29'
        assert param.min_k == 12
        assert param.max_k == 18
        assert param.dmso_percent == 5.0

        reset_to_defaults()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
