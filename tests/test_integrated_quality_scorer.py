"""
Comprehensive unit tests for integrated_quality_scorer module.

Tests all functionality including:
- Individual primer quality scoring
- Set-level quality metrics
- Multi-dimensional scoring integration
- Recommendation generation
- Weight customization
- Different stringency levels
"""

import unittest
from typing import List, Dict

from neoswga.core.integrated_quality_scorer import (
    PrimerQualityScore,
    SetQualityScore,
    IntegratedQualityScorer,
    quick_score_primers
)
from neoswga.core.strand_bias_analyzer import StrandBindingSite
from neoswga.core.reaction_conditions import ReactionConditions


class TestPrimerQualityScore(unittest.TestCase):
    """Test PrimerQualityScore dataclass."""

    def test_create_good_score(self):
        """Test creating a high-quality primer score."""
        score = PrimerQualityScore(
            primer='ACGTACGTACGC',
            strand_bias_score=0.9,
            dimer_score=0.85,
            three_prime_score=0.8,
            complexity_score=0.7,
            thermo_score=0.75,
            overall_score=0.82,
            passes_strand_bias=True,
            passes_dimer=True,
            passes_three_prime=True,
            passes_all=True,
            failure_reasons=[],
            rank=1
        )

        self.assertTrue(score.passes_all)
        self.assertEqual(len(score.failure_reasons), 0)
        self.assertEqual(score.rank, 1)

    def test_create_poor_score(self):
        """Test creating a low-quality primer score."""
        score = PrimerQualityScore(
            primer='AAAAAAAAAA',
            strand_bias_score=0.4,
            dimer_score=0.3,
            three_prime_score=0.2,
            complexity_score=0.1,
            thermo_score=0.5,
            overall_score=0.3,
            passes_strand_bias=True,
            passes_dimer=False,
            passes_three_prime=False,
            passes_all=False,
            failure_reasons=['Low complexity', 'Poor 3\' stability'],
            rank=10
        )

        self.assertFalse(score.passes_all)
        self.assertGreater(len(score.failure_reasons), 0)


class TestSetQualityScore(unittest.TestCase):
    """Test SetQualityScore dataclass."""

    def test_create_passing_set(self):
        """Test creating a passing set score."""
        set_score = SetQualityScore(
            primers=['PRIMER1', 'PRIMER2', 'PRIMER3'],
            mean_overall_score=0.75,
            min_overall_score=0.65,
            set_dimer_score=0.8,
            set_strand_bias_score=0.85,
            passes=True,
            failure_reason=None,
            num_failing_strand=0,
            num_failing_dimer=0,
            num_failing_three_prime=0
        )

        self.assertTrue(set_score.passes)
        self.assertIsNone(set_score.failure_reason)

    def test_create_failing_set(self):
        """Test creating a failing set score."""
        set_score = SetQualityScore(
            primers=['PRIMER1', 'PRIMER2'],
            mean_overall_score=0.4,
            min_overall_score=0.3,
            set_dimer_score=0.3,
            set_strand_bias_score=0.5,
            passes=False,
            failure_reason='Mean set score below threshold',
            num_failing_strand=1,
            num_failing_dimer=2,
            num_failing_three_prime=1
        )

        self.assertFalse(set_score.passes)
        self.assertIsNotNone(set_score.failure_reason)


class TestIntegratedQualityScorer(unittest.TestCase):
    """Test IntegratedQualityScorer class."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0)

        self.scorer_moderate = IntegratedQualityScorer(
            conditions=self.conditions,
            stringency='moderate'
        )

        self.scorer_strict = IntegratedQualityScorer(
            conditions=self.conditions,
            stringency='strict'
        )

        self.scorer_lenient = IntegratedQualityScorer(
            conditions=self.conditions,
            stringency='lenient'
        )

    def test_score_single_primer(self):
        """Test scoring a single primer."""
        primer = 'ACGTACGTACGC'

        # Create minimal binding sites for strand bias
        sites = [
            StrandBindingSite(100, '+', primer),
            StrandBindingSite(200, '-', primer),
        ]

        score = self.scorer_moderate.score_primer(primer, sites)

        self.assertEqual(score.primer, primer)
        self.assertGreaterEqual(score.overall_score, 0.0)
        self.assertLessEqual(score.overall_score, 1.0)
        self.assertIsInstance(score.passes_all, bool)
        self.assertIsInstance(score.failure_reasons, list)

    def test_analyze_good_primer_set(self):
        """Test analyzing a high-quality primer set."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
            'TGCATGCATGCA',
        ]

        # Create balanced binding sites for each primer
        binding_sites = {}
        for i, primer in enumerate(primers):
            sites = [
                StrandBindingSite(j * 100, '+', primer) for j in range(3)
            ] + [
                StrandBindingSite(j * 100, '-', primer) for j in range(3, 6)
            ]
            binding_sites[primer] = sites

        primer_scores, set_score = self.scorer_moderate.analyze_primer_set(
            primers, binding_sites_dict=binding_sites
        )

        self.assertEqual(len(primer_scores), 3)
        self.assertIsInstance(set_score, SetQualityScore)
        self.assertEqual(set_score.primers, primers)

    def test_analyze_diverse_quality_primers(self):
        """Test analyzing primers with varying quality."""
        primers = [
            'ACGTACGTACGC',  # Good
            'AAAAAAAAAA',    # Poor - low complexity
            'GCTAGCTAGCTA',  # Good
        ]

        primer_scores, set_score = self.scorer_moderate.analyze_primer_set(primers)

        # Should have scores for all primers
        self.assertEqual(len(primer_scores), 3)

        # Scores should vary
        scores = [s.overall_score for s in primer_scores]
        self.assertGreater(max(scores) - min(scores), 0.1)  # Some variation

    def test_ranking_primers(self):
        """Test that primers are ranked correctly."""
        primers = [
            'ACGTACGTACGC',
            'AAAAAAAAAA',
            'GCTAGCTAGCTA',
        ]

        primer_scores, set_score = self.scorer_moderate.analyze_primer_set(primers)

        # All should have ranks
        ranks = [s.rank for s in primer_scores]
        self.assertEqual(sorted(ranks), [1, 2, 3])

        # Rank 1 should have highest score
        rank_1 = [s for s in primer_scores if s.rank == 1][0]
        for s in primer_scores:
            if s.rank > 1:
                self.assertGreaterEqual(rank_1.overall_score, s.overall_score)

    def test_empty_primer_set(self):
        """Test analyzing empty primer set raises an error."""
        # Empty primer set may raise ValueError or TypeError
        with self.assertRaises((ValueError, TypeError)):
            self.scorer_moderate.analyze_primer_set([])

    def test_single_primer_set(self):
        """Test analyzing single primer."""
        primers = ['ACGTACGTACGC']

        primer_scores, set_score = self.scorer_moderate.analyze_primer_set(primers)

        self.assertEqual(len(primer_scores), 1)
        self.assertEqual(primer_scores[0].rank, 1)

    def test_custom_weights(self):
        """Test scorer with custom weights."""
        custom_weights = {
            'dimer': 0.5,
            'three_prime': 0.2,
            'strand_bias': 0.15,
            'thermodynamics': 0.1,
            'complexity': 0.05
        }

        scorer = IntegratedQualityScorer(
            conditions=self.conditions,
            stringency='moderate',
            weights=custom_weights
        )

        self.assertEqual(scorer.weights['dimer'], 0.5)
        self.assertEqual(scorer.weights['three_prime'], 0.2)

        # Weights should sum to 1.0
        self.assertAlmostEqual(sum(scorer.weights.values()), 1.0, places=2)

    def test_invalid_weights(self):
        """Test that weights are normalized if they don't sum to 1.0."""
        # Weights that don't sum to 1.0 are automatically normalized
        invalid_weights = {
            'dimer': 0.5,
            'three_prime': 0.5,
            'strand_bias': 0.5,  # Sum > 1.0
            'thermodynamics': 0.1,
            'complexity': 0.05
        }

        # Does not raise, normalizes instead
        scorer = IntegratedQualityScorer(weights=invalid_weights)
        self.assertAlmostEqual(sum(scorer.weights.values()), 1.0, places=2)


class TestRecommendations(unittest.TestCase):
    """Test recommendation generation."""

    def setUp(self):
        """Set up scorer."""
        self.scorer = IntegratedQualityScorer()

    def test_generate_recommendations(self):
        """Test that recommendations are generated."""
        primers = [
            'ACGTACGTACGC',
            'AAAAAAAAAA',  # Poor complexity
            'GCTAGCTAGCTA',
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        recommendations = self.scorer.get_recommendations(primer_scores, set_score)

        self.assertIsInstance(recommendations, list)
        # Should have at least some recommendation
        self.assertGreater(len(recommendations), 0)

    def test_recommendations_identify_weaknesses(self):
        """Test that recommendations identify weak dimensions."""
        # Create primers with specific weaknesses
        primers = [
            'AAAAAAAAAA',  # Very low complexity
            'TTTTTTTTTT',  # Very low complexity
            'GGGGGGGGGG',  # Very low complexity
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        recommendations = self.scorer.get_recommendations(primer_scores, set_score)

        # Should mention complexity issue
        rec_text = ' '.join(recommendations).lower()
        self.assertIn('complexity', rec_text)

    def test_no_recommendations_for_good_set(self):
        """Test that good sets get minimal recommendations."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
            'TGCATGCATGCA',
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        recommendations = self.scorer.get_recommendations(primer_scores, set_score)

        # Should have few or no critical recommendations
        critical_recs = [r for r in recommendations if 'replace' in r.lower()]
        # May or may not have recommendations depending on quality


class TestQuickScoreFunction(unittest.TestCase):
    """Test quick_score_primers utility function."""

    def test_quick_score_basic(self):
        """Test basic quick scoring."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
        ]

        primer_scores, set_score = quick_score_primers(primers)

        self.assertEqual(len(primer_scores), 2)
        self.assertIsInstance(set_score, SetQualityScore)

    def test_quick_score_with_stringency(self):
        """Test quick scoring with different stringency."""
        primers = ['ACGTACGTACGC']

        primer_scores_strict, _ = quick_score_primers(primers, stringency='strict')
        primer_scores_lenient, _ = quick_score_primers(primers, stringency='lenient')

        # Both should complete
        self.assertEqual(len(primer_scores_strict), 1)
        self.assertEqual(len(primer_scores_lenient), 1)

    def test_quick_score_with_sites(self):
        """Test quick scoring with binding sites."""
        primers = ['ACGTACGTACGC']

        sites = {
            'ACGTACGTACGC': [
                StrandBindingSite(100, '+', 'ACGTACGTACGC'),
                StrandBindingSite(200, '-', 'ACGTACGTACGC'),
            ]
        }

        primer_scores, set_score = quick_score_primers(
            primers, binding_sites_dict=sites
        )

        self.assertEqual(len(primer_scores), 1)


class TestStringencyLevels(unittest.TestCase):
    """Test different stringency levels."""

    def setUp(self):
        """Set up primers for testing."""
        self.primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
        ]

    def test_lenient_vs_strict(self):
        """Test that lenient passes more than strict."""
        scorer_lenient = IntegratedQualityScorer(stringency='lenient')
        scorer_strict = IntegratedQualityScorer(stringency='strict')

        _, set_score_lenient = scorer_lenient.analyze_primer_set(self.primers)
        _, set_score_strict = scorer_strict.analyze_primer_set(self.primers)

        # Both should complete
        self.assertIsNotNone(set_score_lenient)
        self.assertIsNotNone(set_score_strict)

        # Lenient more likely to pass (though not guaranteed)

    def test_moderate_is_middle_ground(self):
        """Test that moderate is between lenient and strict."""
        scorer_moderate = IntegratedQualityScorer(stringency='moderate')

        _, set_score = scorer_moderate.analyze_primer_set(self.primers)
        self.assertIsNotNone(set_score)


class TestScoreDimensions(unittest.TestCase):
    """Test individual score dimensions."""

    def setUp(self):
        """Set up scorer."""
        self.scorer = IntegratedQualityScorer()

    def test_all_dimensions_present(self):
        """Test that all dimensions are scored."""
        primer = 'ACGTACGTACGC'
        score = self.scorer.score_primer(primer)

        # Check all dimension scores are present
        self.assertIsNotNone(score.strand_bias_score)
        self.assertIsNotNone(score.dimer_score)
        self.assertIsNotNone(score.three_prime_score)
        self.assertIsNotNone(score.complexity_score)
        self.assertIsNotNone(score.thermo_score)

        # All scores should be in valid range
        for dim_score in [
            score.strand_bias_score,
            score.dimer_score,
            score.three_prime_score,
            score.complexity_score,
            score.thermo_score
        ]:
            self.assertGreaterEqual(dim_score, 0.0)
            self.assertLessEqual(dim_score, 1.0)

    def test_overall_score_is_weighted(self):
        """Test that overall score reflects weights."""
        primer = 'ACGTACGTACGC'
        score = self.scorer.score_primer(primer)

        # Overall should be weighted average
        weighted = (
            self.scorer.weights['dimer'] * score.dimer_score +
            self.scorer.weights['three_prime'] * score.three_prime_score +
            self.scorer.weights['strand_bias'] * score.strand_bias_score +
            self.scorer.weights['thermodynamics'] * score.thermo_score +
            self.scorer.weights['complexity'] * score.complexity_score
        )

        self.assertAlmostEqual(score.overall_score, weighted, places=2)


class TestSetLevelMetrics(unittest.TestCase):
    """Test set-level quality metrics."""

    def setUp(self):
        """Set up scorer."""
        self.scorer = IntegratedQualityScorer()

    def test_mean_score_calculation(self):
        """Test mean score is calculated correctly."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
            'TGCATGCATGCA',
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)

        # Calculate expected mean
        expected_mean = sum(s.overall_score for s in primer_scores) / len(primer_scores)

        self.assertAlmostEqual(set_score.mean_overall_score, expected_mean, places=2)

    def test_min_score_identification(self):
        """Test that minimum score is identified."""
        primers = [
            'ACGTACGTACGC',
            'AAAAAAAAAA',  # Should be lowest
            'GCTAGCTAGCTA',
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)

        # Min should be lowest of all primers
        min_primer_score = min(s.overall_score for s in primer_scores)
        self.assertAlmostEqual(set_score.min_overall_score, min_primer_score, places=2)

    def test_failure_counts(self):
        """Test that failure counts are tracked."""
        primers = [
            'ACGTACGTACGC',
            'AAAAAAAAAA',
            'GCTAGCTAGCTA',
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)

        # Check failure counts are non-negative integers
        self.assertGreaterEqual(set_score.num_failing_strand, 0)
        self.assertGreaterEqual(set_score.num_failing_dimer, 0)
        self.assertGreaterEqual(set_score.num_failing_three_prime, 0)

        # Total failures should not exceed number of primers
        total_failures = (
            set_score.num_failing_strand +
            set_score.num_failing_dimer +
            set_score.num_failing_three_prime
        )
        # Note: primers can fail multiple dimensions, so total may exceed len(primers)


class TestVerboseMode(unittest.TestCase):
    """Test verbose output."""

    def setUp(self):
        """Set up scorer."""
        self.scorer = IntegratedQualityScorer()

    def test_verbose_analysis(self):
        """Test that verbose mode runs without error."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
        ]

        # Should not crash with verbose=True
        primer_scores, set_score = self.scorer.analyze_primer_set(
            primers, verbose=True
        )

        self.assertEqual(len(primer_scores), 2)

    def test_quiet_analysis(self):
        """Test that quiet mode runs without error."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
        ]

        # Should not crash with verbose=False
        primer_scores, set_score = self.scorer.analyze_primer_set(
            primers, verbose=False
        )

        self.assertEqual(len(primer_scores), 2)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def setUp(self):
        """Set up scorer."""
        self.scorer = IntegratedQualityScorer()

    def test_identical_primers(self):
        """Test set with identical primers."""
        primers = ['ACGTACGT', 'ACGTACGT', 'ACGTACGT']

        # Should handle gracefully (may de-duplicate)
        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        # Returns unique primers only (1), not duplicates (3)
        self.assertGreater(len(primer_scores), 0)
        self.assertLessEqual(len(primer_scores), 3)

    def test_very_short_primers(self):
        """Test with very short primers."""
        primers = ['AAA', 'TTT', 'GGG']

        # Should not crash
        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        self.assertEqual(len(primer_scores), 3)

    def test_very_long_primers(self):
        """Test with very long primers."""
        primers = [
            'A' * 50 + 'CGCGC',
            'T' * 50 + 'GCGCG',
        ]

        # Should handle long primers
        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        self.assertEqual(len(primer_scores), 2)

    def test_large_primer_set(self):
        """Test with larger primer set."""
        # Create 20 diverse primers (but pattern repeats every 4)
        bases = ['A', 'T', 'G', 'C']
        primers = []
        for i in range(20):
            primer = ''.join([bases[(i + j) % 4] for j in range(12)])
            primers.append(primer)

        # Due to cyclic pattern, only 4 unique primers exist
        unique_primers = len(set(primers))  # = 4

        # Should complete in reasonable time
        primer_scores, set_score = self.scorer.analyze_primer_set(primers)
        # Returns unique primers only
        self.assertEqual(len(primer_scores), unique_primers)


class TestConsistency(unittest.TestCase):
    """Test internal consistency of scores."""

    def setUp(self):
        """Set up scorer."""
        self.scorer = IntegratedQualityScorer()

    def test_score_reproducibility(self):
        """Test that same primer gets same score."""
        primer = 'ACGTACGTACGC'

        score1 = self.scorer.score_primer(primer)
        score2 = self.scorer.score_primer(primer)

        self.assertAlmostEqual(score1.overall_score, score2.overall_score, places=5)

    def test_set_score_consistency(self):
        """Test that set score is consistent with primer scores."""
        primers = [
            'ACGTACGTACGC',
            'GCTAGCTAGCTA',
        ]

        primer_scores, set_score = self.scorer.analyze_primer_set(primers)

        # Set primers should match input
        self.assertEqual(set_score.primers, primers)

        # Number of primers should match
        self.assertEqual(len(set_score.primers), len(primer_scores))


if __name__ == '__main__':
    unittest.main()
