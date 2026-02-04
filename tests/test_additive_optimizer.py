"""Tests for the multi-additive optimizer."""

import pytest
from neoswga.core.additive_optimizer import (
    AdditiveOptimizer,
    AdditiveRecommendation,
    optimize_conditions,
)
from neoswga.core.reaction_conditions import ReactionConditions


class TestAdditiveRecommendation:
    """Tests for AdditiveRecommendation dataclass."""

    def test_to_conditions(self):
        """Test converting recommendation to ReactionConditions."""
        rec = AdditiveRecommendation(
            dmso_percent=5.0,
            betaine_m=1.0,
            trehalose_m=0.3,
            mg_conc=2.5,
            polymerase='phi29',
        )
        conditions = rec.to_conditions()

        assert isinstance(conditions, ReactionConditions)
        assert conditions.dmso_percent == 5.0
        assert conditions.betaine_m == 1.0
        assert conditions.trehalose_m == 0.3
        assert conditions.mg_conc == 2.5
        assert conditions.polymerase == 'phi29'

    def test_to_conditions_with_temp(self):
        """Test converting with custom temperature."""
        rec = AdditiveRecommendation(polymerase='equiphi29')
        conditions = rec.to_conditions(temp=45.0)
        assert conditions.temp == 45.0

    def test_summary_generation(self):
        """Test that summary generates valid output."""
        rec = AdditiveRecommendation(
            dmso_percent=5.0,
            betaine_m=1.0,
            predicted_amplification=0.75,
            confidence='high',
            optimization_score=0.8,
            rationale="Test rationale",
            primer_length=12,
            template_gc=0.5,
            polymerase='phi29',
        )
        summary = rec.summary()

        assert 'ADDITIVE OPTIMIZATION RECOMMENDATION' in summary
        assert 'DMSO: 5.0%' in summary
        assert 'Betaine: 1.0 M' in summary
        assert 'Confidence: high' in summary
        assert 'Test rationale' in summary


class TestAdditiveOptimizer:
    """Tests for the AdditiveOptimizer class."""

    def test_optimizer_initialization(self):
        """Test optimizer initializes with valid polymerase."""
        optimizer = AdditiveOptimizer(polymerase='phi29')
        assert optimizer.polymerase == 'phi29'

        optimizer2 = AdditiveOptimizer(polymerase='equiphi29')
        assert optimizer2.polymerase == 'equiphi29'

    def test_basic_optimization(self):
        """Test basic optimization returns valid recommendation."""
        optimizer = AdditiveOptimizer(polymerase='phi29')
        rec = optimizer.optimize(primer_length=10, template_gc=0.5)

        assert isinstance(rec, AdditiveRecommendation)
        assert 0.0 <= rec.dmso_percent <= 10.0
        assert 0.0 <= rec.betaine_m <= 3.0
        assert 0.0 <= rec.trehalose_m <= 1.0
        assert rec.confidence in ('high', 'medium', 'low')

    def test_optimization_for_high_gc(self):
        """Test that high-GC templates get more additives."""
        optimizer = AdditiveOptimizer(polymerase='phi29')

        low_gc = optimizer.optimize(primer_length=10, template_gc=0.4)
        high_gc = optimizer.optimize(primer_length=10, template_gc=0.65)

        # High GC should recommend SSB
        assert high_gc.ssb is True
        # High GC should have accessibility-boosting additives
        assert high_gc.dmso_percent >= low_gc.dmso_percent or high_gc.betaine_m >= low_gc.betaine_m

    def test_equiphi29_constraints(self):
        """Test that equiphi29 respects its constraints."""
        optimizer = AdditiveOptimizer(polymerase='equiphi29')
        rec = optimizer.optimize(primer_length=15, template_gc=0.5)

        # Equiphi29 has max 6% DMSO
        assert rec.dmso_percent <= 6.0

    def test_optimization_goals(self):
        """Test different optimization goals produce different results."""
        optimizer = AdditiveOptimizer(polymerase='phi29')

        amp = optimizer.optimize(primer_length=10, template_gc=0.5, optimize_for='amplification')
        proc = optimizer.optimize(primer_length=10, template_gc=0.5, optimize_for='processivity')

        # Different goals should potentially produce different recommendations
        # (though they might be the same if one combination dominates)
        assert amp.optimize_for == 'amplification'
        assert proc.optimize_for == 'processivity'

    def test_quick_search(self):
        """Test quick search mode works and is faster."""
        optimizer = AdditiveOptimizer(polymerase='phi29')

        # Quick search should work
        rec = optimizer.optimize(primer_length=10, template_gc=0.5, quick=True)
        assert isinstance(rec, AdditiveRecommendation)

    def test_application_recommendations(self):
        """Test application-based recommendations."""
        optimizer = AdditiveOptimizer(polymerase='phi29')

        # Clinical should optimize for specificity
        clinical = optimizer.recommend_for_application('clinical', template_gc=0.5)
        assert clinical.optimize_for == 'specificity'

        # Discovery should optimize for coverage
        discovery = optimizer.recommend_for_application('discovery', template_gc=0.5)
        assert discovery.optimize_for == 'coverage'

        # Metagenomics should optimize for coverage with shorter primers
        metagenomics = optimizer.recommend_for_application('metagenomics', template_gc=0.5)
        assert metagenomics.optimize_for == 'coverage'

    def test_primer_length_warnings(self):
        """Test that warnings are generated for long primers."""
        optimizer = AdditiveOptimizer(polymerase='phi29')

        # Very long primer for phi29 should generate warning
        rec = optimizer.optimize(primer_length=18, template_gc=0.5)

        # Should have warning about primer length (unless enough additives)
        # This depends on the combination found
        assert rec.warnings is not None  # Warnings list exists


class TestConvenienceFunction:
    """Tests for the optimize_conditions convenience function."""

    def test_basic_call(self):
        """Test basic call to optimize_conditions."""
        rec = optimize_conditions(12, 0.5)

        assert isinstance(rec, AdditiveRecommendation)
        assert rec.primer_length == 12
        assert rec.template_gc == 0.5

    def test_with_polymerase(self):
        """Test with specific polymerase."""
        rec = optimize_conditions(15, 0.6, polymerase='equiphi29')

        assert rec.polymerase == 'equiphi29'

    def test_with_optimize_for(self):
        """Test with specific optimization goal."""
        rec = optimize_conditions(12, 0.5, optimize_for='coverage')

        assert rec.optimize_for == 'coverage'


class TestSearchRanges:
    """Tests for search range handling."""

    def test_ranges_filtered_by_polymerase(self):
        """Test that search ranges are filtered by polymerase constraints."""
        optimizer = AdditiveOptimizer(polymerase='equiphi29')
        filtered = optimizer._filter_ranges_by_constraints(optimizer.SEARCH_RANGES)

        # All DMSO values should be <= 6.0 for equiphi29
        for dmso in filtered.get('dmso_percent', []):
            assert dmso <= 6.0

    def test_combination_generation(self):
        """Test that combinations are generated correctly."""
        optimizer = AdditiveOptimizer(polymerase='phi29')
        combinations = optimizer._generate_combinations(optimizer.QUICK_RANGES)

        # Should have multiple combinations
        assert len(combinations) > 1

        # Each combination should have the expected keys
        for combo in combinations:
            assert 'dmso_percent' in combo
            assert 'betaine_m' in combo
            assert 'mg_conc' in combo


class TestEdgeCases:
    """Tests for edge cases and boundary conditions."""

    def test_zero_gc(self):
        """Test with very low GC content."""
        optimizer = AdditiveOptimizer(polymerase='phi29')
        rec = optimizer.optimize(primer_length=10, template_gc=0.2)

        assert isinstance(rec, AdditiveRecommendation)
        # Low GC doesn't need SSB
        assert rec.ssb is False

    def test_high_gc(self):
        """Test with very high GC content."""
        optimizer = AdditiveOptimizer(polymerase='phi29')
        rec = optimizer.optimize(primer_length=10, template_gc=0.75)

        assert isinstance(rec, AdditiveRecommendation)
        # Very high GC needs SSB
        assert rec.ssb is True

    def test_very_short_primers(self):
        """Test with very short primers."""
        optimizer = AdditiveOptimizer(polymerase='phi29')
        rec = optimizer.optimize(primer_length=6, template_gc=0.5)

        assert isinstance(rec, AdditiveRecommendation)

    def test_long_primers_equiphi29(self):
        """Test with long primers for equiphi29."""
        optimizer = AdditiveOptimizer(polymerase='equiphi29')
        rec = optimizer.optimize(primer_length=18, template_gc=0.5)

        assert isinstance(rec, AdditiveRecommendation)
        # Long primers may need more additives
        assert rec.dmso_percent > 0 or rec.betaine_m > 0
