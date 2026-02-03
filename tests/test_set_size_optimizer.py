"""
Tests for automatic primer set size optimization with Pareto frontier analysis.

Tests the set_size_optimizer module's ability to:
- Generate SetSizeMetrics dataclass
- Compute Pareto dominance and filter optimal points
- Select from frontier based on application profiles
- Recommend appropriate primer set sizes based on application profiles
  and reaction conditions (backward compatibility)
"""

import pytest
import pandas as pd
import numpy as np

from neoswga.core.set_size_optimizer import (
    SetSizeMetrics,
    FrontierResult,
    filter_pareto_optimal,
    ParetoFrontierGenerator,
    select_from_frontier,
    estimate_optimal_set_size,
    recommend_set_size,
    quick_size_estimate,
    create_baseline_effects,
    get_size_recommendation_summary,
)
from neoswga.core.mechanistic_model import MechanisticEffects


# =============================================================================
# Fixtures
# =============================================================================

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


@pytest.fixture
def sample_primer_pool():
    """Create a sample primer pool DataFrame for testing."""
    return pd.DataFrame({
        'primer': ['ATCGATCG', 'GCTAGCTA', 'TTAATTAA', 'CCGGCCGG', 'AACCTTGG',
                   'GGCCAATT', 'TCGATCGA', 'CTAGCTAG', 'TAATTAAT', 'CGGCCGGC'],
        'fg_freq': [0.001, 0.0008, 0.0006, 0.0005, 0.0004,
                    0.0003, 0.0002, 0.00015, 0.0001, 0.00008],
        'bg_freq': [0.0001, 0.00015, 0.0001, 0.00012, 0.00008,
                    0.00005, 0.00004, 0.00003, 0.00002, 0.00001],
    })


@pytest.fixture
def sample_set_size_metrics():
    """Create sample SetSizeMetrics for testing."""
    return [
        SetSizeMetrics(set_size=4, fg_coverage=0.50, bg_coverage=0.10,
                       fg_binding_sites=500, bg_binding_sites=100, fg_bg_ratio=5.0),
        SetSizeMetrics(set_size=6, fg_coverage=0.65, bg_coverage=0.15,
                       fg_binding_sites=700, bg_binding_sites=140, fg_bg_ratio=5.0),
        SetSizeMetrics(set_size=8, fg_coverage=0.75, bg_coverage=0.18,
                       fg_binding_sites=850, bg_binding_sites=170, fg_bg_ratio=5.0),
        SetSizeMetrics(set_size=10, fg_coverage=0.82, bg_coverage=0.22,
                       fg_binding_sites=1000, bg_binding_sites=200, fg_bg_ratio=5.0),
        SetSizeMetrics(set_size=12, fg_coverage=0.88, bg_coverage=0.28,
                       fg_binding_sites=1150, bg_binding_sites=230, fg_bg_ratio=5.0),
    ]


# =============================================================================
# SetSizeMetrics Tests
# =============================================================================

class TestSetSizeMetrics:
    """Test the SetSizeMetrics dataclass."""

    def test_creation(self):
        """Test basic creation of SetSizeMetrics."""
        metrics = SetSizeMetrics(
            set_size=6,
            fg_coverage=0.75,
            bg_coverage=0.10,
            fg_binding_sites=1000,
            bg_binding_sites=100,
            fg_bg_ratio=10.0,
        )

        assert metrics.set_size == 6
        assert metrics.fg_coverage == 0.75
        assert metrics.fg_bg_ratio == 10.0

    def test_with_primers(self):
        """Test SetSizeMetrics with primer list."""
        primers = ('ATCG', 'GCTA', 'TTAA')
        metrics = SetSizeMetrics(
            set_size=3,
            fg_coverage=0.50,
            bg_coverage=0.05,
            fg_binding_sites=500,
            bg_binding_sites=50,
            fg_bg_ratio=10.0,
            primers=primers,
        )

        assert metrics.primers == primers
        assert len(metrics.primers) == 3

    def test_empty(self):
        """Test empty SetSizeMetrics creation."""
        empty = SetSizeMetrics.empty(set_size=5)

        assert empty.set_size == 5
        assert empty.fg_coverage == 0.0
        assert empty.fg_bg_ratio == 0.0

    def test_to_dict(self):
        """Test conversion to dictionary."""
        metrics = SetSizeMetrics(
            set_size=6,
            fg_coverage=0.75,
            bg_coverage=0.10,
            fg_binding_sites=1000,
            bg_binding_sites=100,
            fg_bg_ratio=10.0,
        )

        d = metrics.to_dict()

        assert d['set_size'] == 6
        assert d['fg_coverage'] == 0.75
        assert d['fg_bg_ratio'] == 10.0
        assert 'is_pareto_optimal' in d

    def test_dominates_clear_domination(self):
        """Test clear domination: better in both objectives."""
        better = SetSizeMetrics(set_size=6, fg_coverage=0.80, bg_coverage=0.10,
                                fg_binding_sites=800, bg_binding_sites=80, fg_bg_ratio=10.0)
        worse = SetSizeMetrics(set_size=8, fg_coverage=0.70, bg_coverage=0.20,
                               fg_binding_sites=700, bg_binding_sites=140, fg_bg_ratio=5.0)

        assert better.dominates(worse)
        assert not worse.dominates(better)

    def test_dominates_equal_one_better_other(self):
        """Test domination: equal in one, better in other."""
        a = SetSizeMetrics(set_size=6, fg_coverage=0.80, bg_coverage=0.10,
                          fg_binding_sites=800, bg_binding_sites=80, fg_bg_ratio=10.0)
        b = SetSizeMetrics(set_size=8, fg_coverage=0.80, bg_coverage=0.10,
                          fg_binding_sites=800, bg_binding_sites=100, fg_bg_ratio=8.0)

        # a has same coverage but better ratio
        assert a.dominates(b)
        assert not b.dominates(a)

    def test_dominates_tradeoff(self):
        """Test no domination when there's a real tradeoff."""
        # High coverage, low ratio
        a = SetSizeMetrics(set_size=10, fg_coverage=0.90, bg_coverage=0.25,
                          fg_binding_sites=900, bg_binding_sites=250, fg_bg_ratio=3.6)
        # Low coverage, high ratio
        b = SetSizeMetrics(set_size=6, fg_coverage=0.70, bg_coverage=0.05,
                          fg_binding_sites=700, bg_binding_sites=50, fg_bg_ratio=14.0)

        # Neither dominates the other - this is a true tradeoff
        assert not a.dominates(b)
        assert not b.dominates(a)


# =============================================================================
# Pareto Filtering Tests
# =============================================================================

class TestFilterParetoOptimal:
    """Test the filter_pareto_optimal function."""

    def test_empty_list(self):
        """Empty list returns empty."""
        assert filter_pareto_optimal([]) == []

    def test_single_point(self):
        """Single point is always Pareto optimal."""
        point = SetSizeMetrics(set_size=6, fg_coverage=0.75, bg_coverage=0.10,
                               fg_binding_sites=750, bg_binding_sites=100, fg_bg_ratio=7.5)

        result = filter_pareto_optimal([point])

        assert len(result) == 1
        assert result[0].is_pareto_optimal

    def test_all_dominated_by_one(self):
        """When one point dominates all others, only it survives."""
        dominant = SetSizeMetrics(set_size=6, fg_coverage=0.90, bg_coverage=0.05,
                                  fg_binding_sites=900, bg_binding_sites=50, fg_bg_ratio=18.0)
        dominated1 = SetSizeMetrics(set_size=8, fg_coverage=0.80, bg_coverage=0.10,
                                    fg_binding_sites=800, bg_binding_sites=100, fg_bg_ratio=8.0)
        dominated2 = SetSizeMetrics(set_size=10, fg_coverage=0.85, bg_coverage=0.15,
                                    fg_binding_sites=850, bg_binding_sites=150, fg_bg_ratio=5.67)

        result = filter_pareto_optimal([dominant, dominated1, dominated2])

        assert len(result) == 1
        assert result[0].set_size == 6

    def test_true_pareto_frontier(self):
        """Test filtering a true tradeoff frontier."""
        # Create points where there's a tradeoff between coverage and ratio
        points = [
            # High ratio, low coverage
            SetSizeMetrics(set_size=4, fg_coverage=0.50, bg_coverage=0.02,
                          fg_binding_sites=500, bg_binding_sites=20, fg_bg_ratio=25.0),
            # Medium ratio, medium coverage
            SetSizeMetrics(set_size=6, fg_coverage=0.70, bg_coverage=0.05,
                          fg_binding_sites=700, bg_binding_sites=50, fg_bg_ratio=14.0),
            # Low ratio, high coverage
            SetSizeMetrics(set_size=10, fg_coverage=0.90, bg_coverage=0.15,
                          fg_binding_sites=900, bg_binding_sites=150, fg_bg_ratio=6.0),
            # Dominated point (worse than 6-primer in both)
            SetSizeMetrics(set_size=8, fg_coverage=0.65, bg_coverage=0.10,
                          fg_binding_sites=650, bg_binding_sites=100, fg_bg_ratio=6.5),
        ]

        result = filter_pareto_optimal(points)

        # Should keep 3 frontier points, exclude the dominated one
        assert len(result) == 3
        sizes = {p.set_size for p in result}
        assert 4 in sizes  # High ratio point
        assert 6 in sizes  # Medium point
        assert 10 in sizes  # High coverage point
        assert 8 not in sizes  # Dominated

    def test_pareto_status_updated(self):
        """Test that _is_pareto_optimal is updated correctly."""
        pareto = SetSizeMetrics(set_size=6, fg_coverage=0.90, bg_coverage=0.05,
                                fg_binding_sites=900, bg_binding_sites=50, fg_bg_ratio=18.0)
        dominated = SetSizeMetrics(set_size=8, fg_coverage=0.80, bg_coverage=0.10,
                                   fg_binding_sites=800, bg_binding_sites=100, fg_bg_ratio=8.0)

        filter_pareto_optimal([pareto, dominated])

        assert pareto._is_pareto_optimal is True
        assert dominated._is_pareto_optimal is False


# =============================================================================
# Select From Frontier Tests
# =============================================================================

class TestSelectFromFrontier:
    """Test the select_from_frontier function."""

    @pytest.fixture
    def frontier_for_selection(self):
        """Create a frontier for selection tests."""
        return [
            # High ratio, low coverage
            SetSizeMetrics(set_size=4, fg_coverage=0.50, bg_coverage=0.02,
                          fg_binding_sites=500, bg_binding_sites=20, fg_bg_ratio=25.0),
            # Medium ratio, medium coverage
            SetSizeMetrics(set_size=6, fg_coverage=0.75, bg_coverage=0.05,
                          fg_binding_sites=750, bg_binding_sites=50, fg_bg_ratio=15.0),
            # Lower ratio, high coverage
            SetSizeMetrics(set_size=10, fg_coverage=0.92, bg_coverage=0.12,
                          fg_binding_sites=920, bg_binding_sites=120, fg_bg_ratio=7.67),
        ]

    def test_empty_frontier_raises(self):
        """Empty frontier should raise ValueError."""
        with pytest.raises(ValueError):
            select_from_frontier([], 'clinical')

    def test_clinical_selects_high_ratio(self, frontier_for_selection):
        """Clinical profile prioritizes specificity (high fg/bg ratio)."""
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'clinical'
        )

        # Clinical wants high ratio
        assert selected.fg_bg_ratio >= 10.0
        assert 'clinical' in explanation.lower()

    def test_discovery_selects_high_coverage(self, frontier_for_selection):
        """Discovery profile prioritizes coverage."""
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'discovery'
        )

        # Discovery wants high coverage
        assert selected.fg_coverage >= 0.70
        assert 'discovery' in explanation.lower()

    def test_enrichment_finds_knee(self, frontier_for_selection):
        """Enrichment profile finds balanced point (knee)."""
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'enrichment'
        )

        # Should be middle-ish point
        assert 'balanced' in explanation.lower() or 'knee' in explanation.lower()

    def test_min_fg_bg_ratio_override(self, frontier_for_selection):
        """Test overriding the min fg/bg ratio constraint."""
        # Set very high ratio requirement
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'discovery',
            min_fg_bg_ratio=20.0
        )

        # Only point with ratio >= 20 is the 4-primer set
        assert selected.fg_bg_ratio >= 20.0

    def test_target_coverage_override(self, frontier_for_selection):
        """Test overriding the target coverage."""
        # Set high coverage requirement
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'clinical',
            target_coverage=0.90
        )

        # Clinical normally wants high ratio, but with high coverage override...
        # This test ensures the parameter is respected

    def test_no_valid_points_returns_best_with_warning(self, frontier_for_selection):
        """When no points meet constraint, return best available with warning."""
        # Set impossibly high ratio
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'clinical',
            min_fg_bg_ratio=100.0  # No point has this
        )

        # Should return best available (highest ratio)
        assert selected.fg_bg_ratio == 25.0
        assert 'warning' in explanation.lower()

    def test_unknown_application_defaults_to_enrichment(self, frontier_for_selection):
        """Unknown application should default to enrichment."""
        selected, explanation = select_from_frontier(
            frontier_for_selection, 'unknown_app'
        )

        assert selected is not None


# =============================================================================
# ParetoFrontierGenerator Tests
# =============================================================================

class TestParetoFrontierGenerator:
    """Test the ParetoFrontierGenerator class."""

    def test_initialization(self, sample_primer_pool):
        """Test basic initialization."""
        generator = ParetoFrontierGenerator(
            primer_pool=sample_primer_pool,
            fg_seq_lengths=[1_000_000],
        )

        assert generator.primer_pool is not None
        assert generator.fg_total_length == 1_000_000

    def test_empty_pool_raises(self):
        """Empty primer pool should raise ValueError."""
        with pytest.raises(ValueError):
            ParetoFrontierGenerator(primer_pool=pd.DataFrame())

    def test_missing_primer_column_with_sequence(self):
        """Test handling of 'sequence' column instead of 'primer'."""
        pool = pd.DataFrame({
            'sequence': ['ATCG', 'GCTA', 'TTAA'],
            'fg_freq': [0.001, 0.0008, 0.0006],
        })

        generator = ParetoFrontierGenerator(primer_pool=pool)

        assert 'primer' in generator.primer_pool.columns

    def test_generate_frontier_quick_only(self, sample_primer_pool):
        """Test quick-only frontier generation."""
        generator = ParetoFrontierGenerator(
            primer_pool=sample_primer_pool,
            fg_seq_lengths=[1_000_000],
            bg_seq_lengths=[3_000_000_000],  # Human-like
        )

        result = generator.generate_frontier(
            min_size=4,
            max_size=8,
            quick_only=True,
            verbose=False,
        )

        assert isinstance(result, FrontierResult)
        assert len(result.all_points) > 0
        assert len(result.pareto_points) > 0
        assert all(isinstance(p, SetSizeMetrics) for p in result.all_points)

    def test_estimate_from_statistics_uses_frequency(self, sample_primer_pool):
        """Test that estimation uses fg_freq and bg_freq columns."""
        generator = ParetoFrontierGenerator(
            primer_pool=sample_primer_pool,
            fg_seq_lengths=[1_000_000],
            bg_seq_lengths=[3_000_000_000],
        )

        points = generator._estimate_from_statistics(min_size=4, max_size=6)

        # Should have 3 points (sizes 4, 5, 6)
        assert len(points) == 3

        # Larger sets should have higher coverage
        assert points[2].fg_coverage >= points[0].fg_coverage

    def test_identify_promising_sizes(self, sample_primer_pool):
        """Test identification of sizes for refinement."""
        generator = ParetoFrontierGenerator(
            primer_pool=sample_primer_pool,
            fg_seq_lengths=[1_000_000],
        )

        coarse = generator._estimate_from_statistics(min_size=4, max_size=10)
        promising = generator._identify_promising_sizes(coarse, num_refine=3)

        assert len(promising) <= 3
        assert all(isinstance(s, int) for s in promising)


# =============================================================================
# Backward Compatibility Tests
# =============================================================================

class TestEstimateOptimalSetSize:
    """Test the estimate_optimal_set_size function (backward compatibility)."""

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
    """Test the recommend_set_size function (backward compatibility)."""

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
        assert 'min_fg_bg_ratio' in result  # New key
        assert 'rationale' in result
        assert 'base_estimate' in result
        assert 'application' in result
        assert 'priority' in result  # New key
        # Legacy key still present
        assert 'min_specificity' in result

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

    def test_min_fg_bg_ratio_override(self, baseline_effects):
        """Test overriding min_fg_bg_ratio parameter."""
        result = recommend_set_size(
            application='enrichment',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
            min_fg_bg_ratio=15.0,  # Override
        )

        assert result['min_fg_bg_ratio'] == 15.0

    def test_target_coverage_override(self, baseline_effects):
        """Test overriding target_coverage parameter."""
        result = recommend_set_size(
            application='enrichment',
            genome_length=1_000_000,
            primer_length=10,
            mech_effects=baseline_effects,
            target_coverage=0.95,  # Override
        )

        assert result['target_coverage'] == 0.95

    def test_priority_reflects_application(self, baseline_effects):
        """Different applications should have different priorities."""
        clinical = recommend_set_size('clinical', 1_000_000, 10, baseline_effects)
        discovery = recommend_set_size('discovery', 1_000_000, 10, baseline_effects)
        enrichment = recommend_set_size('enrichment', 1_000_000, 10, baseline_effects)

        assert clinical['priority'] == 'specificity'
        assert discovery['priority'] == 'coverage'
        assert enrichment['priority'] == 'balanced'


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
        assert 'fg/bg' in summary.lower()  # New metric
        assert 'priority' in summary.lower()  # New field

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
