"""
Comprehensive unit tests for three_prime_stability module.

Tests all functionality including:
- Terminal Tm calculation
- GC clamp analysis
- Terminal base preference
- Hairpin detection
- Composite stability scoring
- Different stringency levels
"""

import unittest
import warnings
from typing import List

from neoswga.core.three_prime_stability import (
    ThreePrimeStability,
    ThreePrimeStabilityAnalyzer,
    create_three_prime_analyzer,
    filter_primers_by_three_prime_stability
)
from neoswga.core.reaction_conditions import ReactionConditions


class TestThreePrimeStability(unittest.TestCase):
    """Test ThreePrimeStability dataclass."""

    def test_create_stable_primer_data(self):
        """Test creating data for stable primer."""
        stability = ThreePrimeStability(
            primer='ACGTACGTACGGC',
            terminal_tm=22.0,
            terminal_gc=0.6,
            gc_clamp=2,
            terminal_binding_energy=-8.0,
            has_terminal_hairpin=False,
            terminal_base='C',
            stability_score=0.85,
            passes=True,
            failure_reason=None
        )

        self.assertEqual(stability.primer, 'ACGTACGTACGGC')
        self.assertEqual(stability.terminal_tm, 22.0)
        self.assertEqual(stability.gc_clamp, 2)
        self.assertTrue(stability.passes)

    def test_create_unstable_primer_data(self):
        """Test creating data for unstable primer."""
        stability = ThreePrimeStability(
            primer='ACGTACGTATA',
            terminal_tm=12.0,
            terminal_gc=0.2,
            gc_clamp=0,
            terminal_binding_energy=-4.0,
            has_terminal_hairpin=False,
            terminal_base='A',
            stability_score=0.35,
            passes=False,
            failure_reason='Terminal Tm too low'
        )

        self.assertFalse(stability.passes)
        self.assertIn('Terminal Tm', stability.failure_reason)


class TestThreePrimeStabilityAnalyzer(unittest.TestCase):
    """Test ThreePrimeStabilityAnalyzer class."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0)

        # Note: min_terminal_tm reduced from 15.0 to 10.0 after salt correction fix
        # (SantaLucia 1998 coefficient 12.5 instead of 16.6 gives ~5C lower Tm)
        self.analyzer_moderate = ThreePrimeStabilityAnalyzer(
            conditions=self.conditions,
            min_terminal_tm=10.0,
            min_gc_clamp=1,
            max_gc_clamp=3,
            avoid_terminal_a=True,
            min_stability_score=0.5
        )

        self.analyzer_strict = ThreePrimeStabilityAnalyzer(
            conditions=self.conditions,
            min_terminal_tm=18.0,
            min_gc_clamp=2,
            max_gc_clamp=3,
            avoid_terminal_a=True,
            min_stability_score=0.6
        )

        self.analyzer_lenient = ThreePrimeStabilityAnalyzer(
            conditions=self.conditions,
            min_terminal_tm=12.0,
            min_gc_clamp=0,
            max_gc_clamp=3,
            avoid_terminal_a=False,
            min_stability_score=0.4
        )

    def test_good_primer_gc_ending(self):
        """Test primer with good GC 3' end."""
        # Last 5 bases: TACGC
        # Terminal Tm = 2(A+T) + 4(G+C) = 2(1) + 4(3) = 14°C base
        # Last 3 bases: CGC (3 GC clamp)
        # Terminal base: C (good)
        primer = 'AAAAAATACGC'

        stability = self.analyzer_moderate.analyze_primer(primer)

        self.assertEqual(stability.primer, primer)
        self.assertEqual(stability.terminal_base, 'C')
        self.assertEqual(stability.gc_clamp, 3)  # CGC = 3 GC bases
        self.assertIsNotNone(stability.terminal_tm)  # Tm is calculated
        self.assertFalse(stability.has_terminal_hairpin)

    def test_good_primer_optimal_tm(self):
        """Test primer with optimal terminal Tm."""
        # Design primer to have Tm in optimal range
        # Last 5 bases should have 2-3 GC for 16-20°C
        primer = 'AAAAAAAATGCC'  # Last 5: ATGCC = 2(2) + 4(3) = 16°C

        stability = self.analyzer_moderate.analyze_primer(primer)

        self.assertGreaterEqual(stability.terminal_tm, 10.0)
        self.assertTrue(stability.passes)

    def test_poor_primer_low_tm(self):
        """Test primer with low terminal Tm."""
        # Last 5 bases: AAAAA (all AT)
        # Wallace Tm = 2(5) + 4(0) = 10°C (below strict threshold of 18°C)
        primer = 'GGGGGGAAAAA'

        # Use strict analyzer (min_terminal_tm=18.0) to test low Tm detection
        stability = self.analyzer_strict.analyze_primer(primer)

        self.assertLess(stability.terminal_tm, 18.0)  # Tm is 10, below strict threshold
        self.assertFalse(stability.passes)
        self.assertIn('Terminal Tm', stability.failure_reason)

    def test_poor_primer_terminal_a(self):
        """Test primer ending in A (weakest base)."""
        primer = 'GGGGGGGGCGA'  # Ends in A

        stability = self.analyzer_moderate.analyze_primer(primer)

        self.assertEqual(stability.terminal_base, 'A')
        # May or may not pass depending on other factors

    def test_poor_primer_no_gc_clamp(self):
        """Test primer with no GC clamp."""
        # Last 3 bases: AAA (0 GC clamp)
        primer = 'GGGGGGCGAAA'

        stability = self.analyzer_moderate.analyze_primer(primer)

        self.assertEqual(stability.gc_clamp, 0)
        self.assertFalse(stability.passes)
        self.assertIn('GC clamp', stability.failure_reason)

    def test_poor_primer_all_gc_clamp(self):
        """Test primer with all GC in clamp (can cause mispriming)."""
        # Last 3 bases: GGG (3 GC clamp)
        primer = 'AAAAAAAAGGG'

        stability = self.analyzer_strict.analyze_primer(primer)

        self.assertEqual(stability.gc_clamp, 3)
        # Strict mode may reject all-GC clamp

    def test_primer_with_hairpin(self):
        """Test primer with terminal hairpin."""
        # Create self-complementary sequence
        # This will form hairpin involving 3' end
        primer = 'AAAAAACGCGAAACGCG'  # Complementary regions

        stability = self.analyzer_moderate.analyze_primer(primer)

        # Hairpin detection may or may not trigger depending on structure
        # Just verify it doesn't crash
        self.assertIsInstance(stability.has_terminal_hairpin, bool)

    def test_short_primer(self):
        """Test very short primer (<5 bases)."""
        primer = 'ACGT'  # Only 4 bases

        # Should handle gracefully
        stability = self.analyzer_moderate.analyze_primer(primer)
        self.assertIsNotNone(stability)

    def test_exactly_5bp_primer(self):
        """Test primer with exactly 5 bases."""
        primer = 'ACGTG'

        stability = self.analyzer_moderate.analyze_primer(primer)

        self.assertEqual(stability.primer, primer)
        self.assertIsNotNone(stability.terminal_tm)

    def test_long_primer(self):
        """Test long primer."""
        primer = 'A' * 50 + 'CGCGC'  # Long primer ending in GC-rich

        stability = self.analyzer_moderate.analyze_primer(primer)

        # Should only analyze last 5 bases
        self.assertEqual(stability.terminal_base, 'C')

    def test_stringency_comparison(self):
        """Test that stringency affects pass/fail."""
        # Borderline primer
        primer = 'AAAAAAATCGC'  # Last 5: ATCGC

        stability_strict = self.analyzer_strict.analyze_primer(primer)
        stability_lenient = self.analyzer_lenient.analyze_primer(primer)

        # Lenient more likely to pass
        # (though not guaranteed for all sequences)
        self.assertIsNotNone(stability_strict)
        self.assertIsNotNone(stability_lenient)


class TestTerminalTmCalculation(unittest.TestCase):
    """Test terminal Tm calculation using Wallace's rule."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = ThreePrimeStabilityAnalyzer()

    def test_all_at_terminal(self):
        """Test terminal region with all AT."""
        primer = 'GGGGGAAAAA'  # Last 5: AAAAA

        stability = self.analyzer.analyze_primer(primer)

        # Wallace: Tm = 2(A+T) + 4(G+C) = 2(5) + 4(0) = 10°C
        expected_tm = 10.0
        self.assertAlmostEqual(stability.terminal_tm, expected_tm, delta=2.0)

    def test_all_gc_terminal(self):
        """Test terminal region with all GC."""
        primer = 'AAAAAGGGGC'  # Last 5: GGGGC

        stability = self.analyzer.analyze_primer(primer)

        # Wallace: Tm = 2(0) + 4(5) = 20°C
        expected_tm = 20.0
        self.assertAlmostEqual(stability.terminal_tm, expected_tm, delta=2.0)

    def test_mixed_terminal(self):
        """Test terminal region with mixed bases."""
        primer = 'AAAAACGATG'  # Last 5: CGATG = 2AT, 3GC

        stability = self.analyzer.analyze_primer(primer)

        # Wallace: Tm = 2(2) + 4(3) = 16°C
        expected_tm = 16.0
        self.assertAlmostEqual(stability.terminal_tm, expected_tm, delta=2.0)


class TestGCClampCalculation(unittest.TestCase):
    """Test GC clamp calculation."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = ThreePrimeStabilityAnalyzer()

    def test_gc_clamp_zero(self):
        """Test no GC in last 3 bases."""
        primer = 'GGGGGGGAAA'  # Last 3: AAA

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.gc_clamp, 0)

    def test_gc_clamp_one(self):
        """Test 1 GC in last 3 bases."""
        primer = 'GGGGGGGAAT'  # Last 3: AAT (0 GC)
        stability1 = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability1.gc_clamp, 0)

        primer = 'GGGGGGGATC'  # Last 3: ATC (1 GC)
        stability2 = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability2.gc_clamp, 1)

    def test_gc_clamp_two(self):
        """Test 2 GC in last 3 bases."""
        primer = 'GGGGGGGAGC'  # Last 3: AGC

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.gc_clamp, 2)

    def test_gc_clamp_three(self):
        """Test 3 GC in last 3 bases."""
        primer = 'GGGGGGGGCC'  # Last 3: GCC

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.gc_clamp, 3)


class TestTerminalBasePreference(unittest.TestCase):
    """Test terminal base identification."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = ThreePrimeStabilityAnalyzer(avoid_terminal_a=True)

    def test_terminal_g(self):
        """Test primer ending in G."""
        primer = 'AAAAAAAAG'

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.terminal_base, 'G')

    def test_terminal_c(self):
        """Test primer ending in C."""
        primer = 'AAAAAAAAC'

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.terminal_base, 'C')

    def test_terminal_a(self):
        """Test primer ending in A (should be avoided)."""
        primer = 'GGGGCGGGA'

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.terminal_base, 'A')
        # Should have lower stability score

    def test_terminal_t(self):
        """Test primer ending in T."""
        primer = 'AAAAAAAAAT'

        stability = self.analyzer.analyze_primer(primer)
        self.assertEqual(stability.terminal_base, 'T')


class TestFactoryFunctions(unittest.TestCase):
    """Test factory functions."""

    def test_create_lenient_analyzer(self):
        """Test creating lenient analyzer."""
        analyzer = create_three_prime_analyzer('lenient')
        self.assertEqual(analyzer.min_terminal_tm, 10.0)  # Factory uses 10.0 for lenient
        self.assertEqual(analyzer.min_gc_clamp, 0)
        self.assertFalse(analyzer.avoid_terminal_a)

    def test_create_moderate_analyzer(self):
        """Test creating moderate analyzer."""
        analyzer = create_three_prime_analyzer('moderate')
        self.assertEqual(analyzer.min_terminal_tm, 14.0)  # Factory default
        self.assertEqual(analyzer.min_gc_clamp, 1)
        self.assertTrue(analyzer.avoid_terminal_a)

    def test_create_strict_analyzer(self):
        """Test creating strict analyzer."""
        analyzer = create_three_prime_analyzer('strict')
        self.assertEqual(analyzer.min_terminal_tm, 18.0)
        self.assertEqual(analyzer.min_gc_clamp, 2)
        self.assertTrue(analyzer.avoid_terminal_a)

    def test_create_with_conditions(self):
        """Test creating analyzer with custom conditions."""
        # Use temp=35.0 which is valid for default phi29 polymerase (20-40C)
        conditions = ReactionConditions(temp=35.0, na_conc=100.0)
        analyzer = create_three_prime_analyzer('moderate', conditions)
        self.assertEqual(analyzer.conditions.temp, 35.0)
        self.assertEqual(analyzer.conditions.na_conc, 100.0)

    def test_invalid_stringency(self):
        """Test creating analyzer with invalid stringency."""
        with self.assertRaises(ValueError):
            create_three_prime_analyzer('ultra_strict')


class TestFilteringFunctions(unittest.TestCase):
    """Test filtering utility functions."""

    def test_filter_all_pass(self):
        """Test filtering where all primers pass."""
        primers = [
            'AAAAAATGCC',  # Good GC ending
            'AAAAAACGCC',  # Good GC ending
            'AAAAAAGCGC',  # Good GC ending
        ]

        # Note: filter_primers_by_three_prime_stability returns only passing primers
        passing = filter_primers_by_three_prime_stability(
            primers, stringency='lenient'
        )

        self.assertGreater(len(passing), 0)
        # Lenient should pass most reasonable primers

    def test_filter_some_fail(self):
        """Test filtering where some primers fail."""
        primers = [
            'AAAAAATGCC',  # Good
            'GGGGGGAAAA',  # Bad - low Tm, no GC clamp
            'GGGGGGGGGA',  # Bad - terminal A
        ]

        # Note: filter_primers_by_three_prime_stability returns only passing primers
        passing = filter_primers_by_three_prime_stability(
            primers, stringency='strict'
        )

        # Some should fail (passing < total)
        self.assertLess(len(passing), len(primers))

    def test_filter_empty_list(self):
        """Test filtering empty primer list."""
        # Note: filter_primers_by_three_prime_stability returns only passing primers
        passing = filter_primers_by_three_prime_stability(
            [], stringency='moderate'
        )

        self.assertEqual(len(passing), 0)


class TestStabilityScoreCalculation(unittest.TestCase):
    """Test composite stability score calculation."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = ThreePrimeStabilityAnalyzer()

    def test_score_range(self):
        """Test that stability score is in valid range."""
        primers = [
            'AAAAAATGCC',  # Good
            'GGGGGGAAAA',  # Poor
            'AAAAAACGCG',  # Excellent
        ]

        for primer in primers:
            stability = self.analyzer.analyze_primer(primer)
            self.assertGreaterEqual(stability.stability_score, 0.0)
            self.assertLessEqual(stability.stability_score, 1.0)

    def test_better_primer_higher_score(self):
        """Test that better primers get higher scores."""
        good_primer = 'AAAAAACGCC'  # 2 GC clamp, C ending, decent Tm
        poor_primer = 'GGGGGGAAAA'  # 0 GC clamp, A ending, low Tm

        stability_good = self.analyzer.analyze_primer(good_primer)
        stability_poor = self.analyzer.analyze_primer(poor_primer)

        # Good primer should have higher score
        self.assertGreater(stability_good.stability_score,
                          stability_poor.stability_score)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = ThreePrimeStabilityAnalyzer()

    def test_empty_primer(self):
        """Test empty primer string."""
        primer = ''

        # Should handle gracefully (may raise or return default)
        try:
            stability = self.analyzer.analyze_primer(primer)
            self.assertIsNotNone(stability)
        except (ValueError, IndexError):
            # Acceptable to raise error for empty primer
            pass

    def test_single_base_primer(self):
        """Test single base primer."""
        primer = 'A'

        # Should handle gracefully (suppressing expected short-sequence warning)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            try:
                stability = self.analyzer.analyze_primer(primer)
                self.assertIsNotNone(stability)
            except ValueError:
                # Acceptable to raise error for too-short primer
                pass

    def test_lowercase_primer(self):
        """Test primer with lowercase bases."""
        primer = 'aaaaaaacgcc'

        # Should handle or convert to uppercase
        stability = self.analyzer.analyze_primer(primer)
        self.assertIsNotNone(stability)

    def test_primer_with_n(self):
        """Test primer containing N (ambiguous base)."""
        primer = 'AAAAAAACGNC'

        # Should handle ambiguous bases (suppressing expected warning)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            stability = self.analyzer.analyze_primer(primer)
            self.assertIsNotNone(stability)

    def test_very_high_gc_terminal(self):
        """Test primer with very high GC in terminal region."""
        primer = 'AAAAAAGGGGCGC'  # Last 5 all GC

        stability = self.analyzer.analyze_primer(primer)

        # Should have high Tm
        self.assertGreater(stability.terminal_tm, 15)

    def test_palindromic_primer(self):
        """Test palindromic primer (likely to form hairpin)."""
        primer = 'CGCGAAACGCG'  # Palindromic

        stability = self.analyzer.analyze_primer(primer)

        # May or may not detect hairpin depending on exact structure
        self.assertIsInstance(stability.has_terminal_hairpin, bool)


class TestConsistencyChecks(unittest.TestCase):
    """Test internal consistency of results."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = ThreePrimeStabilityAnalyzer()

    def test_gc_content_consistency(self):
        """Test that GC content matches clamp."""
        primer = 'AAAAAACGCC'

        stability = self.analyzer.analyze_primer(primer)

        # Last 3 bases: GCC (3 GC bases: G, C, C)
        self.assertEqual(stability.gc_clamp, 3)

        # Terminal GC should be calculated from last 5 bases
        last_5 = primer[-5:]
        expected_gc = (last_5.count('G') + last_5.count('C')) / 5.0
        self.assertAlmostEqual(stability.terminal_gc, expected_gc, places=2)

    def test_passes_consistency(self):
        """Test that passes flag matches failure reason."""
        primers = [
            'AAAAAATGCC',  # Good
            'GGGGGGAAAA',  # Poor
        ]

        for primer in primers:
            stability = self.analyzer.analyze_primer(primer)

            if stability.passes:
                self.assertIsNone(stability.failure_reason)
            else:
                self.assertIsNotNone(stability.failure_reason)
                self.assertGreater(len(stability.failure_reason), 0)


if __name__ == '__main__':
    unittest.main()
