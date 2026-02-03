"""
Tests for automatic primer set size optimization.

Tests the set_size_optimizer module's ability to recommend
appropriate primer set sizes based on application profiles
and reaction conditions.
"""

import pytest
from neoswga.core.set_size_optimizer import (
    estimate_optimal_set_size,
    recommend_set_size,
    quick_size_estimate,
    create_baseline_effects,
    get_size_recommendation_summary,
)
from neoswga.core.mechanistic_model import MechanisticEffects


@pytest.fixture
def baseline_effects():
    """Create baseline MechanisticEffects for testing."""
    return create_baseline_effects()


@pytest.fixture
def reduced_binding_effects():
    """Create MechanisticEffects with reduced binding rate."""
    return MechanisticEffects(
        tm_correction=0.0,
        effective_tm=35.0,
        accessibility_factor=1.0,
        processivity_factor=1.0,
        speed_factor=1.0,
        stability_factor=1.0,
        kon_factor=1.0,
        koff_factor=1.0,
        effective_binding_rate=0.4,  # Reduced binding
        effective_extension_rate=1.0,
        predicted_amplification_factor=0.4,
    )


@pytest.fixture
def reduced_processivity_effects():
    """Create MechanisticEffects with reduced processivity."""
    return MechanisticEffects(
        tm_correction=0.0,
        effective_tm=35.0,
        accessibility_factor=1.0,
        processivity_factor=0.5,  # Reduced processivity
        speed_factor=1.0,
        stability_factor=1.0,
        kon_factor=1.0,
        koff_factor=1.0,
        effective_binding_rate=0.8,
        effective_extension_rate=0.5,
        predicted_amplification_factor=0.4,
    )


class TestEstimateOptimalSetSize:
    """Test the estimate_optimal_set_size function."""

    def test_returns_positive_integer(self, baseline_effects):
        """Result should be a positive integer."""
        size = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert isinstance(size, int)
        assert size > 0

    def test_bounded_range(self, baseline_effects):
        """Result should be within reasonable bounds (4-20)."""
        size = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert 4 <= size <= 20

    def test_shorter_primers_need_fewer(self, baseline_effects):
        """Shorter primers have more binding sites, need fewer primers."""
        short = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=8,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )
        long = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=15,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert short <= long

    def test_larger_genome_needs_more(self, baseline_effects):
        """Larger genomes need more primers for same coverage."""
        small = estimate_optimal_set_size(
            genome_length=500_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )
        large = estimate_optimal_set_size(
            genome_length=5_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert small <= large

    def test_higher_coverage_needs_more(self, baseline_effects):
        """Higher target coverage needs more primers."""
        low_cov = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.5,
            processivity=70000,
            mech_effects=baseline_effects,
        )
        high_cov = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.95,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert low_cov <= high_cov

    def test_reduced_binding_needs_more(self, baseline_effects, reduced_binding_effects):
        """Reduced binding efficiency means more primers needed."""
        normal = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )
        reduced = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=reduced_binding_effects,
        )

        assert normal <= reduced

    def test_reduced_processivity_needs_more(self, baseline_effects, reduced_processivity_effects):
        """Reduced processivity means more primers needed."""
        normal = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )
        reduced = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=reduced_processivity_effects,
        )

        assert normal <= reduced


class TestRecommendSetSize:
    """Test the recommend_set_size function."""

    def test_returns_dict_with_required_keys(self, baseline_effects):
        """Result should have all required keys."""
        result = recommend_set_size(
            application='clinical',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
        )

        assert 'recommended_size' in result
        assert 'size_range' in result
        assert 'target_coverage' in result
        assert 'min_specificity' in result
        assert 'rationale' in result
        assert 'base_estimate' in result
        assert 'application' in result

    def test_clinical_recommends_fewer_than_discovery(self, baseline_effects):
        """Clinical profile should recommend fewer primers for specificity."""
        clinical = recommend_set_size(
            application='clinical',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
        )
        discovery = recommend_set_size(
            application='discovery',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
        )

        assert clinical['recommended_size'] <= discovery['recommended_size']

    def test_metagenomics_recommends_most(self, baseline_effects):
        """Metagenomics should recommend most primers for diversity."""
        metagenomics = recommend_set_size(
            application='metagenomics',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
        )
        clinical = recommend_set_size(
            application='clinical',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
        )

        assert metagenomics['recommended_size'] >= clinical['recommended_size']

    def test_size_within_range(self, baseline_effects):
        """Recommended size should be within the typical range."""
        for app in ['discovery', 'clinical', 'enrichment', 'metagenomics']:
            result = recommend_set_size(
                application=app,
                genome_length=1_000_000,
                primer_length=10,
                mech_effects=baseline_effects,
            )

            min_size, max_size = result['size_range']
            assert min_size <= result['recommended_size'] <= max_size, \
                f"{app}: {result['recommended_size']} not in {result['size_range']}"

    def test_unknown_application_defaults_to_enrichment(self, baseline_effects):
        """Unknown application should default to enrichment."""
        result = recommend_set_size(
            application='unknown_app',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
        )

        assert result['application'] == 'enrichment'

    def test_different_applications_have_different_coverage_targets(self, baseline_effects):
        """Different applications should have different coverage targets."""
        clinical = recommend_set_size('clinical', 1_000_000, 10, baseline_effects)
        discovery = recommend_set_size('discovery', 1_000_000, 10, baseline_effects)

        assert clinical['target_coverage'] != discovery['target_coverage']

    def test_different_applications_have_different_specificity(self, baseline_effects):
        """Different applications should have different specificity requirements."""
        clinical = recommend_set_size('clinical', 1_000_000, 10, baseline_effects)
        discovery = recommend_set_size('discovery', 1_000_000, 10, baseline_effects)

        # Clinical requires higher specificity
        assert clinical['min_specificity'] > discovery['min_specificity']


class TestQuickSizeEstimate:
    """Test the quick_size_estimate convenience function."""

    def test_returns_positive_integer(self):
        """Result should be a positive integer."""
        size = quick_size_estimate('clinical', 1_000_000)

        assert isinstance(size, int)
        assert size > 0

    def test_clinical_vs_discovery(self):
        """Clinical should give smaller or equal size than discovery."""
        clinical = quick_size_estimate('clinical', 1_000_000)
        discovery = quick_size_estimate('discovery', 1_000_000)

        assert clinical <= discovery

    def test_accepts_all_applications(self):
        """Should accept all valid application names."""
        for app in ['discovery', 'clinical', 'enrichment', 'metagenomics']:
            size = quick_size_estimate(app, 1_000_000)
            assert size > 0


class TestCreateBaselineEffects:
    """Test the create_baseline_effects function."""

    def test_returns_mechanistic_effects(self):
        """Should return a MechanisticEffects instance."""
        effects = create_baseline_effects()

        assert isinstance(effects, MechanisticEffects)

    def test_has_optimal_values(self):
        """Baseline should have optimal (1.0) factor values."""
        effects = create_baseline_effects()

        assert effects.processivity_factor == 1.0
        assert effects.accessibility_factor == 1.0
        assert effects.speed_factor == 1.0
        assert effects.stability_factor == 1.0

    def test_has_reasonable_binding_rate(self):
        """Baseline should have reasonable binding rate."""
        effects = create_baseline_effects()

        assert 0.5 <= effects.effective_binding_rate <= 1.0


class TestGetSizeRecommendationSummary:
    """Test the get_size_recommendation_summary function."""

    def test_returns_string(self, baseline_effects):
        """Should return a formatted string."""
        summary = get_size_recommendation_summary(
            'clinical', 1_000_000, 10, baseline_effects
        )

        assert isinstance(summary, str)

    def test_contains_key_information(self, baseline_effects):
        """Summary should contain key information."""
        summary = get_size_recommendation_summary(
            'clinical', 1_000_000, 10, baseline_effects
        )

        assert 'clinical' in summary.lower()
        assert 'Recommended primers' in summary
        assert 'coverage' in summary.lower()
        assert 'specificity' in summary.lower()

    def test_contains_input_parameters(self, baseline_effects):
        """Summary should show input parameters."""
        summary = get_size_recommendation_summary(
            'clinical', 1_000_000, 10, baseline_effects
        )

        assert '1,000,000' in summary  # Genome length
        assert '10' in summary  # Primer length


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_small_genome(self, baseline_effects):
        """Should handle very small genomes."""
        size = estimate_optimal_set_size(
            genome_length=10_000,
            primer_length=8,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert 4 <= size <= 20

    def test_very_large_genome(self, baseline_effects):
        """Should handle very large genomes."""
        size = estimate_optimal_set_size(
            genome_length=100_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert 4 <= size <= 20

    def test_very_long_primers(self, baseline_effects):
        """Should handle very long primers."""
        size = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=18,
            target_coverage=0.8,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert 4 <= size <= 20

    def test_high_coverage_target(self, baseline_effects):
        """Should handle high coverage targets."""
        size = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.99,
            processivity=70000,
            mech_effects=baseline_effects,
        )

        assert size >= 4

    def test_low_processivity(self, baseline_effects):
        """Should handle low processivity (Bst-like)."""
        size = estimate_optimal_set_size(
            genome_length=1_000_000,
            primer_length=10,
            target_coverage=0.8,
            processivity=2000,  # Bst processivity
            mech_effects=baseline_effects,
        )

        assert size >= 4
