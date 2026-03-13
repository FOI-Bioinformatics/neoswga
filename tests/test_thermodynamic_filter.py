"""
Comprehensive unit tests for thermodynamic_filter module.

Tests all functionality including:
- ThermodynamicCriteria defaults and validation
- PrimerThermodynamics dataclass
- ThermodynamicFilter analysis and filtering
- Factory functions for creating filters
- Adaptive GC range calculations
- Edge cases and error handling
"""

import unittest
import logging
import warnings
from typing import List

from neoswga.core.thermodynamic_filter import (
    ThermodynamicCriteria,
    PrimerThermodynamics,
    ThermodynamicFilter,
    calculate_adaptive_gc_range,
    calculate_adaptive_dimer_threshold,
    create_filter_from_conditions,
    create_thermodynamic_filter_adaptive,
)

# Suppress logging during tests
logging.disable(logging.CRITICAL)


class TestThermodynamicCriteria(unittest.TestCase):
    """Test ThermodynamicCriteria dataclass."""

    def test_default_values(self):
        """Test default criteria values are set correctly."""
        criteria = ThermodynamicCriteria()

        self.assertEqual(criteria.min_tm, 25.0)
        self.assertEqual(criteria.max_tm, 50.0)
        self.assertEqual(criteria.target_tm, 35.0)
        self.assertEqual(criteria.na_conc, 50.0)
        self.assertEqual(criteria.mg_conc, 0.0)
        self.assertEqual(criteria.max_homodimer_dg, -9.0)
        self.assertEqual(criteria.max_heterodimer_dg, -9.0)
        self.assertEqual(criteria.max_hairpin_dg, -2.0)
        self.assertEqual(criteria.min_gc, 0.30)
        self.assertEqual(criteria.max_gc, 0.70)
        self.assertEqual(criteria.reaction_temp, 30.0)

    def test_custom_values(self):
        """Test custom criteria values."""
        criteria = ThermodynamicCriteria(
            min_tm=30.0,
            max_tm=60.0,
            na_conc=100.0,
            min_gc=0.40,
            max_gc=0.60
        )

        self.assertEqual(criteria.min_tm, 30.0)
        self.assertEqual(criteria.max_tm, 60.0)
        self.assertEqual(criteria.na_conc, 100.0)
        self.assertEqual(criteria.min_gc, 0.40)
        self.assertEqual(criteria.max_gc, 0.60)


class TestPrimerThermodynamics(unittest.TestCase):
    """Test PrimerThermodynamics dataclass."""

    def test_passing_primer(self):
        """Test primer that passes all filters."""
        primer = PrimerThermodynamics(
            sequence='ATCGATCGATCG',
            tm=35.0,
            gc=0.50,
            homodimer_dg=-5.0,
            hairpin_dg=-1.0,
            passes_filters=True,
            failure_reasons=[]
        )

        self.assertEqual(primer.sequence, 'ATCGATCGATCG')
        self.assertEqual(primer.tm, 35.0)
        self.assertEqual(primer.gc, 0.50)
        self.assertTrue(primer.passes_filters)
        self.assertEqual(len(primer.failure_reasons), 0)

    def test_failing_primer(self):
        """Test primer that fails filters."""
        primer = PrimerThermodynamics(
            sequence='GGGGGGGGGGGG',
            tm=45.0,
            gc=1.0,
            homodimer_dg=-15.0,
            hairpin_dg=-5.0,
            passes_filters=False,
            failure_reasons=['GC too high', 'Strong homodimer']
        )

        self.assertFalse(primer.passes_filters)
        self.assertEqual(len(primer.failure_reasons), 2)
        self.assertIn('GC too high', primer.failure_reasons)


class TestThermodynamicFilter(unittest.TestCase):
    """Test ThermodynamicFilter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.default_filter = ThermodynamicFilter()
        self.strict_criteria = ThermodynamicCriteria(
            min_tm=30.0,
            max_tm=45.0,
            min_gc=0.40,
            max_gc=0.60,
            max_homodimer_dg=-6.0,
            max_hairpin_dg=-1.0
        )
        self.strict_filter = ThermodynamicFilter(self.strict_criteria)

    def test_filter_initialization_with_defaults(self):
        """Test filter initializes with default criteria."""
        filter_obj = ThermodynamicFilter()
        self.assertIsNotNone(filter_obj.criteria)
        self.assertEqual(filter_obj.criteria.min_tm, 25.0)

    def test_filter_initialization_with_custom_criteria(self):
        """Test filter initializes with custom criteria."""
        self.assertEqual(self.strict_filter.criteria.min_gc, 0.40)
        self.assertEqual(self.strict_filter.criteria.max_gc, 0.60)

    def test_analyze_balanced_primer(self):
        """Test analysis of a balanced primer."""
        result = self.default_filter.analyze_primer('ATCGATCGATCG')

        self.assertEqual(result.sequence, 'ATCGATCGATCG')
        self.assertGreater(result.tm, 0)
        self.assertAlmostEqual(result.gc, 0.50, places=2)
        # Balanced primer should pass with default criteria
        # (depends on actual Tm calculation)

    def test_analyze_gc_rich_primer(self):
        """Test analysis of a GC-rich primer."""
        result = self.default_filter.analyze_primer('GCGCGCGCGCGC')

        self.assertEqual(result.gc, 1.0)
        # Should fail GC check with default criteria (max 70%)
        self.assertFalse(result.passes_filters)
        self.assertTrue(any('GC' in r for r in result.failure_reasons))

    def test_analyze_at_rich_primer(self):
        """Test analysis of an AT-rich primer."""
        result = self.default_filter.analyze_primer('ATAATATATATAT')

        self.assertLess(result.gc, 0.30)
        # Should fail GC check with default criteria (min 30%)
        self.assertFalse(result.passes_filters)
        self.assertTrue(any('GC' in r for r in result.failure_reasons))

    def test_filter_empty_list(self):
        """Test filtering empty primer list."""
        passing, stats = self.default_filter.filter_candidates([])

        self.assertEqual(len(passing), 0)
        self.assertEqual(stats['total_candidates'], 0)
        self.assertEqual(stats['final_passing'], 0)

    def test_filter_duplicates_allowed(self):
        """Test that duplicate primers are processed."""
        candidates = ['ATCGATCGATCG', 'ATCGATCGATCG', 'ATCGATCGATCG']
        passing, stats = self.default_filter.filter_candidates(
            candidates, check_heterodimers=False
        )

        self.assertEqual(stats['total_candidates'], 3)
        # Duplicates should be processed (same result for each)

    def test_filter_statistics(self):
        """Test that filtering returns correct statistics."""
        candidates = ['ATCGATCGATCG', 'GCTAGCTAGCTA', 'AAAAAAAAAAAA']
        passing, stats = self.default_filter.filter_candidates(
            candidates, check_heterodimers=False
        )

        self.assertIn('total_candidates', stats)
        self.assertIn('final_passing', stats)
        self.assertIn('mean_tm', stats)
        self.assertIn('mean_gc', stats)
        self.assertEqual(stats['total_candidates'], 3)

    def test_tm_range_enforcement(self):
        """Test that Tm range is enforced."""
        # Use criteria with narrow Tm range
        narrow_criteria = ThermodynamicCriteria(
            min_tm=30.0,
            max_tm=35.0
        )
        filter_obj = ThermodynamicFilter(narrow_criteria)

        # Very short primer should have low Tm
        result = filter_obj.analyze_primer('ATCG')
        # Tm should be outside narrow range
        if result.tm < 30.0:
            self.assertFalse(result.passes_filters)
            self.assertTrue(any('Tm' in r for r in result.failure_reasons))

    def test_gc_range_enforcement(self):
        """Test that GC range is enforced."""
        # Test min_gc enforcement
        result = self.default_filter.analyze_primer('AAAAAAAAAA')
        self.assertLess(result.gc, 0.30)
        self.assertFalse(result.passes_filters)

        # Test max_gc enforcement
        result = self.default_filter.analyze_primer('GGGGGGGGGG')
        self.assertGreater(result.gc, 0.70)
        self.assertFalse(result.passes_filters)

    def test_homodimer_threshold(self):
        """Test homodimer threshold enforcement."""
        strict_criteria = ThermodynamicCriteria(
            max_homodimer_dg=-3.0  # Very strict threshold
        )
        filter_obj = ThermodynamicFilter(strict_criteria)

        # Self-complementary sequences form strong homodimers
        result = filter_obj.analyze_primer('GCGCGCGCGC')
        # Should detect homodimer (strong self-complementarity)
        self.assertLessEqual(result.homodimer_dg, 0.0)

    def test_adjust_criteria_for_conditions(self):
        """Test criteria adjustment for different conditions."""
        filter_obj = ThermodynamicFilter()
        original_min_tm = filter_obj.criteria.min_tm

        filter_obj.adjust_criteria_for_conditions(
            temperature=42.0,
            gc_content=0.50,
            betaine_m=0.0
        )

        # Tm range should adjust for higher temperature
        self.assertEqual(filter_obj.criteria.reaction_temp, 42.0)
        self.assertEqual(filter_obj.criteria.min_tm, 37.0)  # 42 - 5
        self.assertEqual(filter_obj.criteria.max_tm, 62.0)  # 42 + 20

    def test_adjust_criteria_with_betaine(self):
        """Test criteria adjustment with betaine."""
        # Use custom criteria with higher min_gc to allow room to decrease
        custom_criteria = ThermodynamicCriteria(min_gc=0.40, max_gc=0.60)
        filter_obj = ThermodynamicFilter(custom_criteria)
        original_min_gc = filter_obj.criteria.min_gc
        original_max_gc = filter_obj.criteria.max_gc

        filter_obj.adjust_criteria_for_conditions(
            temperature=30.0,
            gc_content=0.50,
            betaine_m=1.0
        )

        # GC range should widen with betaine
        self.assertLess(filter_obj.criteria.min_gc, original_min_gc)
        self.assertGreater(filter_obj.criteria.max_gc, original_max_gc)

    def test_adjust_criteria_high_temp_relaxes_dimers(self):
        """Test that high temperature relaxes dimer constraints."""
        filter_obj = ThermodynamicFilter()
        original_homodimer = filter_obj.criteria.max_homodimer_dg

        filter_obj.adjust_criteria_for_conditions(
            temperature=45.0,  # >= 42
            gc_content=0.50,
            betaine_m=0.0
        )

        # Dimer threshold should be more lenient (more negative)
        self.assertLess(filter_obj.criteria.max_homodimer_dg, original_homodimer)


class TestAdaptiveGCRange(unittest.TestCase):
    """Test calculate_adaptive_gc_range function."""

    def test_at_rich_genome_francisella(self):
        """Test adaptive GC range for AT-rich genome (Francisella ~32%)."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.32)

        # Range should be centered around 32%
        self.assertLess(min_gc, 0.30)  # Below default minimum
        self.assertLess(max_gc, 0.70)  # Below default maximum
        # Width should be approximately 30% (±15%)
        self.assertAlmostEqual(max_gc - min_gc, 0.30, places=1)

    def test_gc_rich_genome_burkholderia(self):
        """Test adaptive GC range for GC-rich genome (Burkholderia ~67%)."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.67)

        # Range should be centered around 67%
        self.assertGreater(min_gc, 0.30)  # Above default minimum
        self.assertGreater(max_gc, 0.70)  # Above default maximum
        # Width should be approximately 30% (±15%)
        self.assertAlmostEqual(max_gc - min_gc, 0.30, places=1)

    def test_balanced_genome(self):
        """Test adaptive GC range for balanced genome (~50%)."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.50)

        # Range should be centered around 50%
        self.assertAlmostEqual(min_gc, 0.35, places=1)
        self.assertAlmostEqual(max_gc, 0.65, places=1)

    def test_extreme_at_rich_floor(self):
        """Test floor enforcement for extreme AT-rich genome."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.10)  # 10% GC

        # Should hit floor at 15%
        self.assertEqual(min_gc, 0.15)
        self.assertAlmostEqual(max_gc, 0.25, places=1)

    def test_extreme_gc_rich_ceiling(self):
        """Test ceiling enforcement for extreme GC-rich genome."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.90)  # 90% GC

        # Should hit ceiling at 85%
        self.assertAlmostEqual(min_gc, 0.75, places=1)
        self.assertEqual(max_gc, 0.85)

    def test_custom_tolerance(self):
        """Test custom tolerance parameter."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.50, tolerance=0.25)

        # Width should be 50% (±25%)
        self.assertAlmostEqual(max_gc - min_gc, 0.50, places=1)


class TestAdaptiveDimerThreshold(unittest.TestCase):
    """Test calculate_adaptive_dimer_threshold function."""

    def test_at_rich_more_lenient(self):
        """Test that AT-rich genomes get more lenient threshold."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.32)

        # AT-rich should get less negative (more lenient) threshold
        self.assertGreater(threshold, -9.0)

    def test_gc_rich_more_strict(self):
        """Test that GC-rich genomes get stricter threshold."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.67)

        # GC-rich should get more negative (stricter) threshold
        self.assertLess(threshold, -9.0)

    def test_balanced_unchanged(self):
        """Test that balanced genome keeps same threshold."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.50)

        # Balanced should be approximately unchanged
        self.assertAlmostEqual(threshold, -9.0, places=1)

    def test_floor_enforcement(self):
        """Test floor at -15 kcal/mol."""
        # Very GC-rich should hit floor
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.95)

        self.assertGreaterEqual(threshold, -15.0)

    def test_ceiling_enforcement(self):
        """Test ceiling at -5 kcal/mol."""
        # Very AT-rich should hit ceiling
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.05)

        self.assertLessEqual(threshold, -5.0)


class TestFactoryFunctions(unittest.TestCase):
    """Test factory functions for creating filters."""

    def test_create_filter_from_conditions_phi29(self):
        """Test creating filter for Phi29 conditions."""
        filter_obj = create_filter_from_conditions(
            polymerase='phi29',
            temperature=30.0,
            gc_content=0.50,
            betaine_m=0.0
        )

        self.assertIsInstance(filter_obj, ThermodynamicFilter)
        self.assertEqual(filter_obj.criteria.reaction_temp, 30.0)

    def test_create_filter_from_conditions_equiphi29(self):
        """Test creating filter for EquiPhi29 conditions."""
        filter_obj = create_filter_from_conditions(
            polymerase='equiphi29',
            temperature=42.0,
            gc_content=0.50,
            betaine_m=1.0
        )

        self.assertIsInstance(filter_obj, ThermodynamicFilter)
        self.assertEqual(filter_obj.criteria.reaction_temp, 42.0)

    def test_create_adaptive_filter_francisella(self):
        """Test creating adaptive filter for Francisella (AT-rich)."""
        filter_obj = create_thermodynamic_filter_adaptive(genome_gc=0.32)

        self.assertIsInstance(filter_obj, ThermodynamicFilter)
        # Should have wider low-GC range
        self.assertLess(filter_obj.criteria.min_gc, 0.30)
        # Should have more lenient dimer threshold
        self.assertGreater(filter_obj.criteria.max_homodimer_dg, -9.0)

    def test_create_adaptive_filter_burkholderia(self):
        """Test creating adaptive filter for Burkholderia (GC-rich)."""
        filter_obj = create_thermodynamic_filter_adaptive(genome_gc=0.67)

        self.assertIsInstance(filter_obj, ThermodynamicFilter)
        # Should have wider high-GC range
        self.assertGreater(filter_obj.criteria.max_gc, 0.70)
        # Should have stricter dimer threshold
        self.assertLess(filter_obj.criteria.max_homodimer_dg, -9.0)

    def test_create_adaptive_filter_with_betaine(self):
        """Test adaptive filter with betaine expands GC range."""
        filter_without = create_thermodynamic_filter_adaptive(
            genome_gc=0.50, betaine_m=0.0
        )
        filter_with = create_thermodynamic_filter_adaptive(
            genome_gc=0.50, betaine_m=1.0
        )

        # With betaine should have wider GC range
        range_without = filter_without.criteria.max_gc - filter_without.criteria.min_gc
        range_with = filter_with.criteria.max_gc - filter_with.criteria.min_gc
        self.assertGreater(range_with, range_without)

    def test_create_adaptive_filter_high_temp(self):
        """Test adaptive filter at high temperature relaxes dimers."""
        filter_low = create_thermodynamic_filter_adaptive(
            genome_gc=0.50, temperature=30.0
        )
        filter_high = create_thermodynamic_filter_adaptive(
            genome_gc=0.50, temperature=45.0
        )

        # High temp should have more lenient dimer threshold
        self.assertGreater(
            filter_high.criteria.max_homodimer_dg,
            filter_low.criteria.max_homodimer_dg
        )


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def test_analyze_very_short_primer(self):
        """Test analysis of very short primer (4bp)."""
        filter_obj = ThermodynamicFilter()
        result = filter_obj.analyze_primer('ATCG')

        # Should handle short sequences without error
        self.assertIsInstance(result, PrimerThermodynamics)
        self.assertEqual(len(result.sequence), 4)

    def test_analyze_single_base(self):
        """Test analysis of single base."""
        filter_obj = ThermodynamicFilter()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            result = filter_obj.analyze_primer('A')

        # Should handle without crashing
        self.assertIsInstance(result, PrimerThermodynamics)
        self.assertFalse(result.passes_filters)

    def test_analyze_empty_sequence(self):
        """Test analysis of empty sequence."""
        filter_obj = ThermodynamicFilter()
        # Should handle gracefully
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            try:
                result = filter_obj.analyze_primer('')
                self.assertFalse(result.passes_filters)
            except (ValueError, ZeroDivisionError):
                pass  # Acceptable to raise error for empty sequence

    def test_analyze_lowercase_sequence(self):
        """Test analysis of lowercase sequence."""
        filter_obj = ThermodynamicFilter()
        result = filter_obj.analyze_primer('atcgatcgatcg')

        # Should handle lowercase (thermodynamics module normalizes)
        self.assertIsInstance(result, PrimerThermodynamics)

    def test_filter_single_candidate(self):
        """Test filtering single candidate (no heterodimer check possible)."""
        filter_obj = ThermodynamicFilter()
        candidates = ['ATCGATCGATCG']
        passing, stats = filter_obj.filter_candidates(candidates)

        self.assertEqual(stats['total_candidates'], 1)
        # Heterodimer issues should be 0 with single primer
        self.assertEqual(stats['heterodimer_issues'], 0)

    def test_filter_with_heterodimers_disabled(self):
        """Test filtering with heterodimer check disabled."""
        filter_obj = ThermodynamicFilter()
        candidates = ['ATCGATCGATCG', 'GCTAGCTAGCTA', 'AATTAATTAATT']
        passing, stats = filter_obj.filter_candidates(
            candidates, check_heterodimers=False
        )

        # Should not report heterodimer issues when disabled
        self.assertEqual(stats['heterodimer_issues'], 0)


class TestLargeSetPerformance(unittest.TestCase):
    """Test performance with larger primer sets."""

    def test_filter_100_primers(self):
        """Test filtering 100 primers performs reasonably."""
        import time

        filter_obj = ThermodynamicFilter()

        # Generate 100 random-ish primers
        bases = ['A', 'T', 'C', 'G']
        import random
        random.seed(42)
        candidates = [
            ''.join(random.choices(bases, k=10))
            for _ in range(100)
        ]

        start = time.time()
        passing, stats = filter_obj.filter_candidates(
            candidates, check_heterodimers=False
        )
        elapsed = time.time() - start

        # Should complete in reasonable time (< 10 seconds)
        self.assertLess(elapsed, 10.0)
        self.assertEqual(stats['total_candidates'], 100)


if __name__ == '__main__':
    unittest.main()
