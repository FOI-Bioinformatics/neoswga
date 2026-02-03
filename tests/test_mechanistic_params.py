"""
Tests for mechanistic model parameters.

Validates the structure and values of literature-based parameters
used in the mechanistic SWGA model.
"""

import pytest
from neoswga.core.mechanistic_params import (
    MECHANISTIC_MODEL_PARAMS,
    APPLICATION_PROFILES,
    get_polymerase_params,
    get_application_profile,
    list_applications,
    list_polymerases,
)


class TestMechanisticModelParams:
    """Test the MECHANISTIC_MODEL_PARAMS structure."""

    def test_has_all_pathways(self):
        """Verify all four pathways are present."""
        assert 'tm' in MECHANISTIC_MODEL_PARAMS
        assert 'structure' in MECHANISTIC_MODEL_PARAMS
        assert 'enzyme' in MECHANISTIC_MODEL_PARAMS
        assert 'kinetics' in MECHANISTIC_MODEL_PARAMS
        assert 'interactions' in MECHANISTIC_MODEL_PARAMS

    def test_tm_pathway_params(self):
        """Verify Tm pathway has required parameters."""
        tm = MECHANISTIC_MODEL_PARAMS['tm']

        # Uniform correction coefficients
        assert 'dmso_coef' in tm
        assert 'formamide_coef' in tm
        assert 'trehalose_coef' in tm
        assert 'ethanol_coef' in tm
        assert 'betaine_uniform_coef' in tm
        assert 'urea_coef' in tm
        assert 'tmac_uniform_coef' in tm

        # GC sigmoid parameters
        assert 'betaine_gc_midpoint' in tm
        assert 'betaine_gc_steepness' in tm
        assert 'tmac_gc_midpoint' in tm
        assert 'tmac_gc_steepness' in tm

    def test_structure_pathway_params(self):
        """Verify structure pathway has required parameters."""
        structure = MECHANISTIC_MODEL_PARAMS['structure']

        assert 'dmso_melt_rate' in structure
        assert 'dmso_max_effect' in structure
        assert 'betaine_melt_rate' in structure
        assert 'betaine_max_effect' in structure
        assert 'temp_structure_coef' in structure
        assert 'gc_structure_coef' in structure

    def test_enzyme_pathway_params(self):
        """Verify enzyme pathway has required parameters."""
        enzyme = MECHANISTIC_MODEL_PARAMS['enzyme']

        # Polymerase-specific params
        assert 'phi29' in enzyme
        assert 'equiphi29' in enzyme
        assert 'bst' in enzyme
        assert 'klenow' in enzyme

        # General enzyme params
        assert 'betaine_peak' in enzyme
        assert 'betaine_enhancement' in enzyme
        assert 'formamide_coef' in enzyme
        assert 'mg_optimal' in enzyme

    def test_kinetics_pathway_params(self):
        """Verify kinetics pathway has required parameters."""
        kinetics = MECHANISTIC_MODEL_PARAMS['kinetics']

        assert 'optimal_delta_t' in kinetics
        assert 'delta_t_width' in kinetics
        assert 'betaine_kon_boost' in kinetics
        assert 'dmso_kon_boost' in kinetics
        assert 'ssb_kon_multiplier' in kinetics

    def test_interactions_params(self):
        """Verify interaction parameters are present."""
        interactions = MECHANISTIC_MODEL_PARAMS['interactions']

        assert 'dmso_mg_chelation' in interactions
        assert 'betaine_trehalose_synergy' in interactions


class TestPolymeraseParams:
    """Test polymerase-specific parameters."""

    def test_phi29_params(self):
        """Verify phi29 has correct parameters."""
        phi29 = MECHANISTIC_MODEL_PARAMS['enzyme']['phi29']

        assert phi29['optimal_temp'] == 30.0
        assert phi29['processivity'] == 70000
        assert 'extension_rate' in phi29
        assert 'dmso_threshold' in phi29

    def test_equiphi29_params(self):
        """Verify equiphi29 has correct parameters."""
        equiphi29 = MECHANISTIC_MODEL_PARAMS['enzyme']['equiphi29']

        assert equiphi29['optimal_temp'] == 42.0
        assert equiphi29['processivity'] == 80000

    def test_get_polymerase_params(self):
        """Test get_polymerase_params function."""
        phi29 = get_polymerase_params('phi29')
        assert phi29['optimal_temp'] == 30.0

        # Case insensitive
        equiphi = get_polymerase_params('EquiPhi29')
        assert equiphi['optimal_temp'] == 42.0

    def test_get_polymerase_params_invalid(self):
        """Test get_polymerase_params with invalid polymerase."""
        with pytest.raises(ValueError, match="Unknown polymerase"):
            get_polymerase_params('invalid_polymerase')

    def test_list_polymerases(self):
        """Test list_polymerases function."""
        polymerases = list_polymerases()

        assert 'phi29' in polymerases
        assert 'equiphi29' in polymerases
        assert 'bst' in polymerases
        assert 'klenow' in polymerases

        # Check format: (optimal_temp, processivity)
        assert polymerases['phi29'] == (30.0, 70000)


class TestApplicationProfiles:
    """Test application profiles for set size optimization."""

    def test_all_profiles_present(self):
        """Verify all expected profiles exist."""
        assert 'discovery' in APPLICATION_PROFILES
        assert 'clinical' in APPLICATION_PROFILES
        assert 'enrichment' in APPLICATION_PROFILES
        assert 'metagenomics' in APPLICATION_PROFILES

    def test_profile_structure(self):
        """Verify each profile has required fields."""
        for name, profile in APPLICATION_PROFILES.items():
            assert 'target_coverage' in profile, f"{name} missing target_coverage"
            assert 'min_specificity' in profile, f"{name} missing min_specificity"
            assert 'typical_size' in profile, f"{name} missing typical_size"
            assert 'description' in profile, f"{name} missing description"

    def test_coverage_values_valid(self):
        """Verify coverage values are in valid range."""
        for name, profile in APPLICATION_PROFILES.items():
            assert 0 < profile['target_coverage'] <= 1.0, \
                f"{name} target_coverage out of range"
            assert 0 < profile['min_specificity'] <= 1.0, \
                f"{name} min_specificity out of range"

    def test_typical_size_format(self):
        """Verify typical_size is a tuple of (min, max)."""
        for name, profile in APPLICATION_PROFILES.items():
            size = profile['typical_size']
            assert len(size) == 2, f"{name} typical_size wrong length"
            assert size[0] < size[1], f"{name} typical_size min >= max"
            assert size[0] >= 4, f"{name} typical_size min too small"
            assert size[1] <= 25, f"{name} typical_size max too large"

    def test_get_application_profile(self):
        """Test get_application_profile function."""
        clinical = get_application_profile('clinical')

        assert clinical['target_coverage'] == 0.70
        assert clinical['min_specificity'] == 0.90

    def test_get_application_profile_invalid(self):
        """Test get_application_profile with invalid application."""
        with pytest.raises(ValueError, match="Unknown application"):
            get_application_profile('invalid_application')

    def test_list_applications(self):
        """Test list_applications function."""
        apps = list_applications()

        assert 'discovery' in apps
        assert 'clinical' in apps
        assert isinstance(apps['discovery'], str)


class TestParameterValues:
    """Test that parameter values are scientifically reasonable."""

    def test_dmso_coefficient_negative(self):
        """DMSO should lower Tm (negative coefficient)."""
        assert MECHANISTIC_MODEL_PARAMS['tm']['dmso_coef'] < 0

    def test_formamide_coefficient_negative(self):
        """Formamide should lower Tm (negative coefficient)."""
        assert MECHANISTIC_MODEL_PARAMS['tm']['formamide_coef'] < 0

    def test_betaine_enhancement_positive(self):
        """Betaine should enhance enzyme at low concentration."""
        assert MECHANISTIC_MODEL_PARAMS['enzyme']['betaine_enhancement'] > 0

    def test_processivity_ordering(self):
        """Phi29 should have highest processivity among standard polymerases."""
        enzyme = MECHANISTIC_MODEL_PARAMS['enzyme']

        phi29_proc = enzyme['phi29']['processivity']
        bst_proc = enzyme['bst']['processivity']
        klenow_proc = enzyme['klenow']['processivity']

        assert phi29_proc > klenow_proc > bst_proc

    def test_optimal_temps_ordered(self):
        """Optimal temps should increase: phi29 < equiphi29 < bst."""
        enzyme = MECHANISTIC_MODEL_PARAMS['enzyme']

        assert enzyme['phi29']['optimal_temp'] < enzyme['equiphi29']['optimal_temp']
        assert enzyme['equiphi29']['optimal_temp'] < enzyme['bst']['optimal_temp']

    def test_clinical_more_specific_than_discovery(self):
        """Clinical should require higher specificity than discovery."""
        clinical = APPLICATION_PROFILES['clinical']
        discovery = APPLICATION_PROFILES['discovery']

        assert clinical['min_specificity'] > discovery['min_specificity']

    def test_discovery_higher_coverage_than_clinical(self):
        """Discovery should target higher coverage than clinical."""
        clinical = APPLICATION_PROFILES['clinical']
        discovery = APPLICATION_PROFILES['discovery']

        assert discovery['target_coverage'] > clinical['target_coverage']
