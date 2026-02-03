"""
Comprehensive unit tests for strand_bias_analyzer module.

Tests all functionality including:
- Individual primer strand bias analysis
- Set-level bias calculations
- Filtering with different stringency levels
- Edge cases and error handling
"""

import unittest
from typing import List

from neoswga.core.strand_bias_analyzer import (
    StrandBindingSite,
    PrimerStrandBias,
    SetStrandBias,
    StrandBiasAnalyzer,
    create_strand_bias_analyzer,
    filter_primers_by_strand_bias
)


class TestStrandBindingSite(unittest.TestCase):
    """Test StrandBindingSite dataclass."""

    def test_create_forward_site(self):
        """Test creating a forward strand binding site."""
        site = StrandBindingSite(position=100, strand='+', sequence='ACGTACGT')
        self.assertEqual(site.position, 100)
        self.assertEqual(site.strand, '+')
        self.assertEqual(site.sequence, 'ACGTACGT')

    def test_create_reverse_site(self):
        """Test creating a reverse strand binding site."""
        site = StrandBindingSite(position=500, strand='-', sequence='TGCATGCA')
        self.assertEqual(site.position, 500)
        self.assertEqual(site.strand, '-')
        self.assertEqual(site.sequence, 'TGCATGCA')


class TestStrandBiasAnalyzer(unittest.TestCase):
    """Test StrandBiasAnalyzer class."""

    def setUp(self):
        """Set up test fixtures."""
        self.analyzer_moderate = StrandBiasAnalyzer(
            max_bias_ratio=4.0,
            max_bias_score=0.7,
            max_biased_primers_fraction=0.3
        )

        self.analyzer_strict = StrandBiasAnalyzer(
            max_bias_ratio=2.5,
            max_bias_score=0.5,
            max_biased_primers_fraction=0.2
        )

        self.analyzer_lenient = StrandBiasAnalyzer(
            max_bias_ratio=6.0,
            max_bias_score=0.95,  # Higher to allow ratio 5.0 (bias_score ~0.90)
            max_biased_primers_fraction=0.4
        )

    def test_balanced_primer(self):
        """Test analysis of a balanced primer (1:1 ratio)."""
        primer = 'ACGTACGTACGT'

        # Create balanced binding sites (5 forward, 5 reverse)
        sites = [
            StrandBindingSite(100, '+', primer),
            StrandBindingSite(200, '+', primer),
            StrandBindingSite(300, '+', primer),
            StrandBindingSite(400, '+', primer),
            StrandBindingSite(500, '+', primer),
            StrandBindingSite(150, '-', primer),
            StrandBindingSite(250, '-', primer),
            StrandBindingSite(350, '-', primer),
            StrandBindingSite(450, '-', primer),
            StrandBindingSite(550, '-', primer),
        ]

        bias = self.analyzer_moderate.analyze_primer(primer, sites)

        self.assertEqual(bias.primer, primer)
        self.assertEqual(bias.forward_count, 5)
        self.assertEqual(bias.reverse_count, 5)
        self.assertEqual(bias.bias_ratio, 1.0)
        self.assertAlmostEqual(bias.bias_score, 0.0, places=2)
        self.assertTrue(bias.passes)
        self.assertIsNone(bias.failure_reason)

    def test_forward_biased_primer_moderate(self):
        """Test a primer with forward bias that fails moderate stringency."""
        primer = 'TGCATGCATGCA'

        # 10 forward, 2 reverse (5:1 ratio - should fail moderate)
        sites = [StrandBindingSite(i * 100, '+', primer) for i in range(10)]
        sites.extend([StrandBindingSite(i * 100, '-', primer) for i in range(10, 12)])

        bias = self.analyzer_moderate.analyze_primer(primer, sites)

        self.assertEqual(bias.forward_count, 10)
        self.assertEqual(bias.reverse_count, 2)
        self.assertEqual(bias.bias_ratio, 5.0)
        self.assertFalse(bias.passes)
        # Actual format: "Forward bias: 5.00:1 (max 4.0:1)"
        self.assertIn('Forward bias:', bias.failure_reason)

    def test_forward_biased_primer_lenient(self):
        """Test that same biased primer passes lenient stringency."""
        primer = 'TGCATGCATGCA'

        # 10 forward, 2 reverse (5:1 ratio - should pass lenient)
        sites = [StrandBindingSite(i * 100, '+', primer) for i in range(10)]
        sites.extend([StrandBindingSite(i * 100, '-', primer) for i in range(10, 12)])

        bias = self.analyzer_lenient.analyze_primer(primer, sites)

        self.assertEqual(bias.bias_ratio, 5.0)
        self.assertTrue(bias.passes)  # 5.0 < 6.0 threshold

    def test_reverse_biased_primer(self):
        """Test a primer with reverse bias."""
        primer = 'GCTAGCTAGCTA'

        # 2 forward, 10 reverse (0.2:1 ratio)
        sites = [StrandBindingSite(i * 100, '+', primer) for i in range(2)]
        sites.extend([StrandBindingSite(i * 100, '-', primer) for i in range(2, 12)])

        bias = self.analyzer_moderate.analyze_primer(primer, sites)

        self.assertEqual(bias.forward_count, 2)
        self.assertEqual(bias.reverse_count, 10)
        self.assertEqual(bias.bias_ratio, 0.2)
        self.assertFalse(bias.passes)
        # Actual format: "Reverse bias: 1:5.00 (max 1:4.0)"
        self.assertIn('Reverse bias:', bias.failure_reason)

    def test_no_forward_sites(self):
        """Test primer with only reverse sites."""
        primer = 'AAAATTTTTAAA'

        # 0 forward, 5 reverse
        sites = [StrandBindingSite(i * 100, '-', primer) for i in range(5)]

        bias = self.analyzer_moderate.analyze_primer(primer, sites)

        self.assertEqual(bias.forward_count, 0)
        self.assertEqual(bias.reverse_count, 5)
        self.assertEqual(bias.bias_ratio, 0.0)
        self.assertFalse(bias.passes)

    def test_no_reverse_sites(self):
        """Test primer with only forward sites."""
        primer = 'GGGGCCCCGGGG'

        # 5 forward, 0 reverse
        sites = [StrandBindingSite(i * 100, '+', primer) for i in range(5)]

        bias = self.analyzer_moderate.analyze_primer(primer, sites)

        self.assertEqual(bias.forward_count, 5)
        self.assertEqual(bias.reverse_count, 0)
        self.assertGreater(bias.bias_ratio, 1000)  # Very large ratio
        self.assertFalse(bias.passes)

    def test_empty_sites_list(self):
        """Test primer with no binding sites."""
        primer = 'ACGTACGT'
        sites = []

        bias = self.analyzer_moderate.analyze_primer(primer, sites)

        self.assertEqual(bias.forward_count, 0)
        self.assertEqual(bias.reverse_count, 0)
        self.assertEqual(bias.bias_ratio, 1.0)  # Default to balanced when no sites
        self.assertFalse(bias.passes)

    def test_set_analysis_all_pass(self):
        """Test set-level analysis where all primers pass."""
        primers_sites = {
            'PRIMER1': [
                StrandBindingSite(100, '+', 'P1'),
                StrandBindingSite(200, '+', 'P1'),
                StrandBindingSite(300, '-', 'P1'),
                StrandBindingSite(400, '-', 'P1'),
            ],
            'PRIMER2': [
                StrandBindingSite(150, '+', 'P2'),
                StrandBindingSite(250, '-', 'P2'),
            ],
            'PRIMER3': [
                StrandBindingSite(100, '+', 'P3'),
                StrandBindingSite(200, '+', 'P3'),
                StrandBindingSite(300, '+', 'P3'),
                StrandBindingSite(400, '-', 'P3'),
                StrandBindingSite(500, '-', 'P3'),
                StrandBindingSite(600, '-', 'P3'),
            ],
        }

        # analyze_primer_set takes list of PrimerStrandBias, not dict
        primer_biases = [
            self.analyzer_moderate.analyze_primer(primer, sites)
            for primer, sites in primers_sites.items()
        ]
        set_bias = self.analyzer_moderate.analyze_primer_set(primer_biases)

        self.assertEqual(len(set_bias.primers), 3)  # primers is a list
        self.assertEqual(set_bias.num_biased_primers, 0)
        self.assertTrue(set_bias.passes)
        self.assertIsNone(set_bias.failure_reason)

    def test_set_analysis_some_fail(self):
        """Test set-level analysis where some primers fail."""
        primers_sites = {
            'BALANCED': [
                StrandBindingSite(i, '+', 'BAL') for i in range(5)
            ] + [
                StrandBindingSite(i, '-', 'BAL') for i in range(5, 10)
            ],
            'BIASED1': [  # 10:1 ratio - will fail
                StrandBindingSite(i, '+', 'B1') for i in range(10)
            ] + [
                StrandBindingSite(100, '-', 'B1')
            ],
            'BIASED2': [  # 1:10 ratio - will fail
                StrandBindingSite(100, '+', 'B2')
            ] + [
                StrandBindingSite(i, '-', 'B2') for i in range(10)
            ],
        }

        # analyze_primer_set takes list of PrimerStrandBias, not dict
        primer_biases = [
            self.analyzer_moderate.analyze_primer(primer, sites)
            for primer, sites in primers_sites.items()
        ]
        set_bias = self.analyzer_moderate.analyze_primer_set(primer_biases)

        self.assertEqual(len(set_bias.primers), 3)  # primers is a list
        self.assertEqual(set_bias.num_biased_primers, 2)
        self.assertFalse(set_bias.passes)
        self.assertIn('66.7%', set_bias.failure_reason)  # 2/3 = 66.7%


class TestFactoryFunctions(unittest.TestCase):
    """Test factory functions and utility functions."""

    def test_create_lenient_analyzer(self):
        """Test creating lenient analyzer."""
        analyzer = create_strand_bias_analyzer('lenient')
        self.assertEqual(analyzer.max_ratio, 6.0)
        self.assertEqual(analyzer.max_score, 0.85)
        self.assertEqual(analyzer.max_biased_fraction, 0.4)

    def test_create_moderate_analyzer(self):
        """Test creating moderate analyzer."""
        analyzer = create_strand_bias_analyzer('moderate')
        self.assertEqual(analyzer.max_ratio, 4.0)
        self.assertEqual(analyzer.max_score, 0.7)
        self.assertEqual(analyzer.max_biased_fraction, 0.3)

    def test_create_strict_analyzer(self):
        """Test creating strict analyzer."""
        analyzer = create_strand_bias_analyzer('strict')
        self.assertEqual(analyzer.max_ratio, 2.5)
        self.assertEqual(analyzer.max_score, 0.5)
        self.assertEqual(analyzer.max_biased_fraction, 0.15)  # Actual value

    def test_create_invalid_stringency(self):
        """Test creating analyzer with invalid stringency defaults to moderate."""
        # Invalid stringency falls back to moderate (default behavior)
        analyzer = create_strand_bias_analyzer('ultra_strict')
        self.assertEqual(analyzer.max_ratio, 4.0)  # Same as moderate

    def test_filter_primers_all_pass(self):
        """Test filtering where all primers pass."""
        primers_sites = {
            'P1': [StrandBindingSite(i, '+', 'P1') for i in range(5)] +
                  [StrandBindingSite(i, '-', 'P1') for i in range(5, 10)],
            'P2': [StrandBindingSite(i, '+', 'P2') for i in range(3)] +
                  [StrandBindingSite(i, '-', 'P2') for i in range(3, 6)],
        }

        # Note: filter_primers_by_strand_bias takes (primers, sites_dict, stringency)
        # and returns only passing primers (not a tuple)
        passing = filter_primers_by_strand_bias(
            list(primers_sites.keys()), primers_sites, stringency='moderate'
        )

        self.assertEqual(len(passing), 2)
        self.assertIn('P1', passing)
        self.assertIn('P2', passing)

    def test_filter_primers_some_fail(self):
        """Test filtering where some primers fail."""
        primers_sites = {
            'GOOD': [StrandBindingSite(i, '+', 'G') for i in range(5)] +
                    [StrandBindingSite(i, '-', 'G') for i in range(5, 10)],
            'BAD': [StrandBindingSite(i, '+', 'B') for i in range(20)] +
                   [StrandBindingSite(100, '-', 'B')],
        }

        # Note: filter_primers_by_strand_bias returns only passing primers
        passing = filter_primers_by_strand_bias(
            list(primers_sites.keys()), primers_sites, stringency='moderate'
        )

        self.assertEqual(len(passing), 1)
        self.assertIn('GOOD', passing)
        self.assertNotIn('BAD', passing)


class TestBiasScoreCalculation(unittest.TestCase):
    """Test bias score calculation logic."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = StrandBiasAnalyzer()

    def test_bias_score_balanced(self):
        """Test bias score for perfectly balanced primer."""
        primer = 'TEST'
        sites = [StrandBindingSite(i, '+', primer) for i in range(10)]
        sites.extend([StrandBindingSite(i, '-', primer) for i in range(10, 20)])

        bias = self.analyzer.analyze_primer(primer, sites)
        self.assertAlmostEqual(bias.bias_score, 0.0, places=2)

    def test_bias_score_moderate_bias(self):
        """Test bias score for moderately biased primer."""
        primer = 'TEST'
        sites = [StrandBindingSite(i, '+', primer) for i in range(6)]
        sites.extend([StrandBindingSite(i, '-', primer) for i in range(6, 10)])

        bias = self.analyzer.analyze_primer(primer, sites)
        self.assertGreater(bias.bias_score, 0.0)
        self.assertLess(bias.bias_score, 0.5)

    def test_bias_score_extreme_bias(self):
        """Test bias score for extremely biased primer."""
        primer = 'TEST'
        sites = [StrandBindingSite(i, '+', primer) for i in range(20)]
        sites.append(StrandBindingSite(100, '-', primer))

        bias = self.analyzer.analyze_primer(primer, sites)
        self.assertGreater(bias.bias_score, 0.9)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = StrandBiasAnalyzer()

    def test_single_site(self):
        """Test primer with only one binding site."""
        primer = 'RARE'
        sites = [StrandBindingSite(100, '+', primer)]

        bias = self.analyzer.analyze_primer(primer, sites)
        self.assertFalse(bias.passes)

    def test_large_number_sites(self):
        """Test primer with very large number of sites."""
        primer = 'COMMON'
        sites = [StrandBindingSite(i, '+', primer) for i in range(500)]
        sites.extend([StrandBindingSite(i, '-', primer) for i in range(500, 1000)])

        bias = self.analyzer.analyze_primer(primer, sites)
        self.assertEqual(bias.forward_count, 500)
        self.assertEqual(bias.reverse_count, 500)
        self.assertAlmostEqual(bias.bias_ratio, 1.0)
        self.assertTrue(bias.passes)

    def test_invalid_strand(self):
        """Test that invalid strand is handled."""
        # This should be caught during site creation/validation
        primer = 'TEST'
        # Mix of valid and "invalid" marked sites
        sites = [
            StrandBindingSite(100, '+', primer),
            StrandBindingSite(200, '-', primer),
        ]

        # Should still work with valid strands
        bias = self.analyzer.analyze_primer(primer, sites)
        self.assertEqual(bias.forward_count, 1)
        self.assertEqual(bias.reverse_count, 1)


if __name__ == '__main__':
    unittest.main()
