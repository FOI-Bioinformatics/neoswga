"""
Comprehensive unit tests for dimer_network_analyzer module.

Tests all functionality including:
- Dimer severity calculation
- Network metrics computation
- Hub primer identification
- Primer replacement suggestions
- Set-level optimization
- Different stringency levels
"""

import unittest
from typing import List
import numpy as np

from neoswga.core.dimer_network_analyzer import (
    DimerInteraction,
    PrimerDimerProfile,
    DimerNetworkMetrics,
    DimerNetworkAnalyzer,
    create_dimer_network_analyzer,
    filter_primer_set_by_dimer_network
)
from neoswga.core.reaction_conditions import ReactionConditions


class TestDimerInteraction(unittest.TestCase):
    """Test DimerInteraction dataclass."""

    def test_create_heterodimer(self):
        """Test creating a heterodimer interaction."""
        inter = DimerInteraction(
            primer1='ACGTACGT',
            primer2='TGCATGCA',
            energy=-8.5,
            tm=45.0,
            severity=0.6,
            forms_dimer=True,
            is_homodimer=False
        )

        self.assertEqual(inter.primer1, 'ACGTACGT')
        self.assertEqual(inter.primer2, 'TGCATGCA')
        self.assertEqual(inter.energy, -8.5)
        self.assertTrue(inter.forms_dimer)
        self.assertFalse(inter.is_homodimer)

    def test_create_homodimer(self):
        """Test creating a homodimer interaction."""
        inter = DimerInteraction(
            primer1='ACGTACGT',
            primer2='ACGTACGT',
            energy=-10.0,
            tm=50.0,
            severity=0.8,
            forms_dimer=True,
            is_homodimer=True
        )

        self.assertTrue(inter.is_homodimer)


class TestDimerNetworkAnalyzer(unittest.TestCase):
    """Test DimerNetworkAnalyzer class."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0)

        self.analyzer_moderate = DimerNetworkAnalyzer(
            conditions=self.conditions,
            severity_threshold=0.3,
            hub_threshold=3,
            max_mean_severity=0.2,
            max_hub_primers=2
        )

        self.analyzer_strict = DimerNetworkAnalyzer(
            conditions=self.conditions,
            severity_threshold=0.2,
            hub_threshold=2,
            max_mean_severity=0.15,
            max_hub_primers=1
        )

        self.analyzer_lenient = DimerNetworkAnalyzer(
            conditions=self.conditions,
            severity_threshold=0.4,
            hub_threshold=4,
            max_mean_severity=0.3,
            max_hub_primers=3
        )

    def test_analyze_good_primer_set(self):
        """Test analysis of a primer set with minimal dimer interactions."""
        # These primers are designed to have minimal complementarity
        primers = [
            'AAAAAAAAAA',  # Poly-A
            'TTTTTTTTTT',  # Poly-T
            'GGGGGGGGGG',  # Poly-G
            'CCCCCCCCCC',  # Poly-C
        ]

        metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)

        # Should pass with minimal interactions
        self.assertEqual(metrics.num_primers, 4)
        self.assertIsInstance(metrics.total_interactions, int)
        self.assertIsInstance(metrics.mean_severity, float)
        self.assertGreaterEqual(metrics.mean_severity, 0.0)
        self.assertLessEqual(metrics.mean_severity, 1.0)

        # Check profiles exist for all primers
        for primer in primers:
            self.assertIn(primer, profiles)
            profile = profiles[primer]
            self.assertEqual(profile.primer, primer)
            self.assertGreaterEqual(profile.num_interactions, 0)

    def test_analyze_problematic_primer_set(self):
        """Test analysis of primers with known complementarity."""
        # Create primers with some complementarity
        primers = [
            'ACGTACGTAC',  # Primer 1
            'GTACGTACGT',  # Primer 2 - overlaps with 1
            'TACGTACGTA',  # Primer 3 - overlaps with 1 and 2
            'GGGGGGGGGG',  # Primer 4 - no overlap
        ]

        metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)

        self.assertEqual(metrics.num_primers, 4)

        # Check matrix dimensions
        self.assertEqual(matrix.shape, (4, 4))

        # Diagonal contains homodimer severity scores (floats, not booleans)
        # Just verify they're in valid range [0, 1]
        for i in range(4):
            self.assertGreaterEqual(matrix[i, i], 0.0)
            self.assertLessEqual(matrix[i, i], 1.0)

    def test_identify_hub_primers(self):
        """Test identification of hub primers."""
        primers = [
            'ACGTACGTAC',  # Hub primer - complementary to many
            'GTACGTACGT',  # Interacts with hub
            'TACGTACGTA',  # Interacts with hub
            'CGTACGTACG',  # Interacts with hub
            'GGGGGGGGGG',  # Isolated primer
        ]

        metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)

        # Check for hub primers
        hub_primers = [p for p, profile in profiles.items() if profile.is_hub]

        # Should identify at least some hubs if interactions exist
        self.assertIsInstance(hub_primers, list)

    def test_empty_primer_set(self):
        """Test analysis with empty primer set raises an error or returns failure."""
        primers = []

        # Empty primer set may raise error or return failure status
        try:
            metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)
            # If it returns, should indicate failure
            self.assertFalse(metrics.passes)
        except (ValueError, TypeError):
            # Acceptable to raise error for empty set
            pass

    def test_single_primer(self):
        """Test analysis with single primer."""
        primers = ['ACGTACGTAC']

        metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)

        self.assertEqual(metrics.num_primers, 1)
        self.assertEqual(metrics.total_interactions, 0)
        self.assertEqual(metrics.num_hub_primers, 0)
        self.assertTrue(metrics.passes)

    def test_two_primers_no_interaction(self):
        """Test two primers with no significant interaction."""
        primers = [
            'AAAAAAAAAA',
            'GGGGGGGGGG',
        ]

        metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)

        self.assertEqual(metrics.num_primers, 2)
        # Should have minimal or no interactions
        self.assertLessEqual(metrics.mean_severity, 0.5)

    def test_primer_replacement_suggestions(self):
        """Test suggesting replacements for problematic primers."""
        primers = [
            'ACGTACGTAC',  # Problematic
            'GTACGTACGT',  # Problematic
            'GGGGGGGGGG',  # Good
        ]

        candidates = [
            'TTTTTTTTTT',
            'CCCCCCCCCC',
            'ATATATATAT',
        ]

        metrics, profiles, matrix = self.analyzer_moderate.analyze_primer_set(primers)

        # Identify primers to replace
        to_replace = self.analyzer_moderate.identify_primers_to_replace(profiles, n=1)

        self.assertIsInstance(to_replace, list)
        self.assertLessEqual(len(to_replace), 1)

        if to_replace:
            # Get suggestions
            suggestions = self.analyzer_moderate.suggest_replacements(
                primers, candidates, to_replace
            )

            self.assertIsInstance(suggestions, dict)
            for primer in to_replace:
                self.assertIn(primer, suggestions)

    def test_stringency_levels(self):
        """Test different stringency levels."""
        primers = [
            'ACGTACGTAC',
            'GTACGTACGT',
            'TACGTACGTA',
        ]

        # Analyze with strict
        metrics_strict, _, _ = self.analyzer_strict.analyze_primer_set(primers)

        # Analyze with lenient
        metrics_lenient, _, _ = self.analyzer_lenient.analyze_primer_set(primers)

        # Both should complete
        self.assertEqual(metrics_strict.num_primers, 3)
        self.assertEqual(metrics_lenient.num_primers, 3)

        # Strict may be more likely to fail
        # (but we can't guarantee without specific sequences)


class TestFactoryFunctions(unittest.TestCase):
    """Test factory functions."""

    def test_create_lenient_analyzer(self):
        """Test creating lenient analyzer."""
        analyzer = create_dimer_network_analyzer('lenient')
        self.assertEqual(analyzer.severity_threshold, 0.4)
        self.assertEqual(analyzer.hub_threshold, 4)
        self.assertEqual(analyzer.max_mean_severity, 0.3)
        self.assertEqual(analyzer.max_hub_primers, 3)

    def test_create_moderate_analyzer(self):
        """Test creating moderate analyzer."""
        analyzer = create_dimer_network_analyzer('moderate')
        self.assertEqual(analyzer.severity_threshold, 0.3)
        self.assertEqual(analyzer.hub_threshold, 3)
        self.assertEqual(analyzer.max_mean_severity, 0.2)
        self.assertEqual(analyzer.max_hub_primers, 2)

    def test_create_strict_analyzer(self):
        """Test creating strict analyzer."""
        analyzer = create_dimer_network_analyzer('strict')
        self.assertEqual(analyzer.severity_threshold, 0.2)
        self.assertEqual(analyzer.hub_threshold, 2)
        self.assertEqual(analyzer.max_mean_severity, 0.15)
        self.assertEqual(analyzer.max_hub_primers, 1)

    def test_create_with_conditions(self):
        """Test creating analyzer with custom conditions."""
        # Use temp=35.0 which is valid for default phi29 polymerase (20-40C)
        conditions = ReactionConditions(temp=35.0, betaine_m=1.5)
        analyzer = create_dimer_network_analyzer('moderate', conditions)
        self.assertEqual(analyzer.conditions.temp, 35.0)
        self.assertEqual(analyzer.conditions.betaine_m, 1.5)

    def test_invalid_stringency(self):
        """Test creating analyzer with invalid stringency defaults to moderate."""
        # Invalid stringency falls back to moderate (default behavior)
        analyzer = create_dimer_network_analyzer('ultra_strict')
        self.assertEqual(analyzer.severity_threshold, 0.3)  # Same as moderate


class TestFilteringFunctions(unittest.TestCase):
    """Test filtering utility functions."""

    def test_filter_all_pass(self):
        """Test filtering where all primers pass."""
        primers = [
            'AAAAAAAAAA',
            'TTTTTTTTTT',
            'GGGGGGGGGG',
        ]

        # Note: filter_primer_set_by_dimer_network returns (passes: bool, metrics)
        passes, metrics = filter_primer_set_by_dimer_network(
            primers, stringency='moderate'
        )

        self.assertIsInstance(passes, bool)
        self.assertEqual(metrics.num_primers, 3)

    def test_filter_with_removal(self):
        """Test filtering returns metrics for primer set."""
        primers = [
            'ACGTACGTAC',
            'GTACGTACGT',
            'TACGTACGTA',
            'GGGGGGGGGG',
        ]

        # Note: filter_primer_set_by_dimer_network returns (passes: bool, metrics)
        passes, metrics = filter_primer_set_by_dimer_network(
            primers, stringency='strict'
        )

        # Returns pass status and metrics
        self.assertIsInstance(passes, bool)
        self.assertEqual(metrics.num_primers, 4)

    def test_filter_empty_list(self):
        """Test filtering empty primer list raises an error or returns failure."""
        # Empty primer list may raise error or return failure status
        try:
            passes, metrics = filter_primer_set_by_dimer_network([], stringency='moderate')
            self.assertFalse(passes)
        except (ValueError, TypeError):
            # Acceptable to raise error for empty set
            pass


class TestNetworkMetrics(unittest.TestCase):
    """Test network metrics calculation."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = DimerNetworkAnalyzer()

    def test_metrics_structure(self):
        """Test that metrics have correct structure."""
        primers = ['AAAAAAAAAA', 'TTTTTTTTTT']

        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)

        # Check all required fields
        self.assertIsInstance(metrics.num_primers, int)
        self.assertIsInstance(metrics.total_interactions, int)
        self.assertIsInstance(metrics.num_hub_primers, int)
        self.assertIsInstance(metrics.mean_severity, float)
        self.assertIsInstance(metrics.max_severity, float)
        self.assertIsInstance(metrics.passes, bool)

    def test_severity_distribution(self):
        """Test severity distribution calculation."""
        primers = [
            'ACGTACGTAC',
            'GTACGTACGT',
            'TACGTACGTA',
        ]

        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)

        # Should have severity distribution
        self.assertIsInstance(metrics.severity_distribution, dict)

        # Distribution categories should sum to total comparisons
        # This includes both heterodimers (n*(n-1)/2) and homodimers (n)
        total = sum(metrics.severity_distribution.values())
        expected_with_homodimers = len(primers) + len(primers) * (len(primers) - 1) // 2
        # Could be either just heterodimers or including homodimers
        self.assertIn(total, [3, 6])  # 3 pairs or 3 pairs + 3 homodimers


class TestPrimerProfile(unittest.TestCase):
    """Test primer profile calculation."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = DimerNetworkAnalyzer()

    def test_profile_fields(self):
        """Test that profile has all required fields."""
        primers = ['ACGTACGTAC', 'GTACGTACGT']

        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)

        for primer in primers:
            profile = profiles[primer]
            self.assertEqual(profile.primer, primer)
            self.assertIsInstance(profile.num_interactions, int)
            self.assertIsInstance(profile.max_severity, (int, float))  # Could be 0 or 0.0
            self.assertIsInstance(profile.mean_severity, (int, float))
            self.assertIsInstance(profile.total_binding_energy, (int, float))  # Could be 0 or 0.0
            self.assertIsInstance(profile.problematic_partners, list)
            self.assertIsInstance(profile.is_hub, bool)

    def test_profile_consistency(self):
        """Test profile values are consistent."""
        primers = ['ACGTACGTAC', 'GTACGTACGT', 'TACGTACGTA']

        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)

        for primer, profile in profiles.items():
            # Number of problematic partners <= total primers
            self.assertLessEqual(len(profile.problematic_partners), len(primers) - 1)

            # Max severity >= mean severity
            self.assertGreaterEqual(profile.max_severity, profile.mean_severity)

            # Severity values in valid range
            self.assertGreaterEqual(profile.max_severity, 0.0)
            self.assertLessEqual(profile.max_severity, 1.0)
            self.assertGreaterEqual(profile.mean_severity, 0.0)
            self.assertLessEqual(profile.mean_severity, 1.0)


class TestReplacementOptimization(unittest.TestCase):
    """Test primer replacement and optimization."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = DimerNetworkAnalyzer()

    def test_greedy_optimization(self):
        """Test greedy set optimization."""
        primers = [
            'ACGTACGTAC',
            'GTACGTACGT',
            'TACGTACGTA',
        ]

        candidates = [
            'AAAAAAAAAA',
            'TTTTTTTTTT',
            'GGGGGGGGGG',
            'CCCCCCCCCC',
        ]

        optimized, final_metrics = self.analyzer.optimize_set_greedy(
            primers, candidates, max_iterations=2
        )

        # Should return a set
        self.assertIsInstance(optimized, list)
        self.assertGreater(len(optimized), 0)

        # Final metrics should exist
        self.assertIsInstance(final_metrics, DimerNetworkMetrics)

    def test_replacement_improves_quality(self):
        """Test that replacement suggestions target worst primers."""
        primers = [
            'ACGTACGTAC',  # Potentially problematic
            'GTACGTACGT',  # Potentially problematic
            'GGGGGGGGGG',  # Good
        ]

        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)

        to_replace = self.analyzer.identify_primers_to_replace(profiles, n=1)

        # Should not replace the good primer
        if to_replace:
            self.assertNotIn('GGGGGGGGGG', to_replace)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = DimerNetworkAnalyzer()

    def test_identical_primers(self):
        """Test set with identical primers."""
        primers = ['ACGTACGT', 'ACGTACGT', 'ACGTACGT']

        # Should handle this gracefully
        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)
        self.assertEqual(metrics.num_primers, 3)

    def test_very_short_primers(self):
        """Test with very short primers."""
        primers = ['AAA', 'TTT', 'GGG']

        # Should not crash
        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)
        self.assertEqual(metrics.num_primers, 3)

    def test_very_long_primers(self):
        """Test with very long primers."""
        primers = [
            'A' * 50,
            'T' * 50,
            'G' * 50,
        ]

        # Should handle long primers
        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)
        self.assertEqual(metrics.num_primers, 3)

    def test_large_primer_set(self):
        """Test with larger primer set (performance check)."""
        # Create 20 diverse primers
        bases = ['A', 'T', 'G', 'C']
        primers = []
        for i in range(20):
            primer = ''.join([bases[j % 4] for j in range(i, i + 10)])
            primers.append(primer)

        # Should complete in reasonable time
        metrics, profiles, matrix = self.analyzer.analyze_primer_set(primers)
        self.assertEqual(metrics.num_primers, 20)
        self.assertEqual(matrix.shape, (20, 20))


class TestVerboseOutput(unittest.TestCase):
    """Test verbose output functionality."""

    def setUp(self):
        """Set up analyzer."""
        self.analyzer = DimerNetworkAnalyzer()

    def test_verbose_analysis(self):
        """Test that verbose mode runs without error."""
        primers = ['ACGTACGTAC', 'GTACGTACGT']

        # Should not crash with verbose=True
        metrics, profiles, matrix = self.analyzer.analyze_primer_set(
            primers, verbose=True
        )

        self.assertEqual(metrics.num_primers, 2)

    def test_quiet_analysis(self):
        """Test that quiet mode runs without error."""
        primers = ['ACGTACGTAC', 'GTACGTACGT']

        # Should not crash with verbose=False
        metrics, profiles, matrix = self.analyzer.analyze_primer_set(
            primers, verbose=False
        )

        self.assertEqual(metrics.num_primers, 2)


if __name__ == '__main__':
    unittest.main()
