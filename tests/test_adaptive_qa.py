#!/usr/bin/env python3
"""
Unit tests for genome-adaptive QA threshold calculations.

Tests all adaptive threshold functions across:
1. Three-prime stability (terminal Tm adaptation)
2. Thermodynamic filter (GC range and dimer threshold adaptation)
3. Genome analysis utilities
4. Integrated quality scorer

Author: NeoSWGA Development Team
Date: November 2025
"""

import unittest
import tempfile
import os
from pathlib import Path

from neoswga.core.three_prime_stability import (
    calculate_gc_deviation,
    calculate_adaptive_terminal_tm,
    create_three_prime_analyzer_adaptive
)

from neoswga.core.thermodynamic_filter import (
    calculate_adaptive_gc_range,
    calculate_adaptive_dimer_threshold,
    create_thermodynamic_filter_adaptive
)

from neoswga.core.genome_analysis import (
    calculate_genome_gc,
    get_gc_class,
    recommend_adaptive_qa,
    calculate_genome_stats
)

from neoswga.core.integrated_quality_scorer import (
    create_quality_scorer
)


class TestGCDeviation(unittest.TestCase):
    """Test GC deviation calculation."""

    def test_balanced_genome(self):
        """50% GC genome should have zero deviation."""
        deviation = calculate_gc_deviation(0.50)
        self.assertAlmostEqual(deviation, 0.0, places=6)

    def test_at_rich_genome(self):
        """AT-rich genome should have negative deviation."""
        deviation = calculate_gc_deviation(0.32)  # Francisella
        self.assertAlmostEqual(deviation, -0.18, places=6)

    def test_gc_rich_genome(self):
        """GC-rich genome should have positive deviation."""
        deviation = calculate_gc_deviation(0.67)  # Burkholderia
        self.assertAlmostEqual(deviation, 0.17, places=6)

    def test_extreme_at(self):
        """Extreme AT-rich genome."""
        deviation = calculate_gc_deviation(0.19)  # Plasmodium
        self.assertAlmostEqual(deviation, -0.31, places=6)

    def test_extreme_gc(self):
        """Extreme GC-rich genome."""
        deviation = calculate_gc_deviation(0.72)  # Streptomyces
        self.assertAlmostEqual(deviation, 0.22, places=6)


class TestAdaptiveTerminalTm(unittest.TestCase):
    """Test adaptive terminal Tm threshold calculation."""

    def test_balanced_genome_moderate(self):
        """Balanced genome should use standard threshold."""
        tm = calculate_adaptive_terminal_tm(14.0, 0.50, 'moderate')
        self.assertAlmostEqual(tm, 14.0, places=1)

    def test_francisella_moderate(self):
        """Francisella (32% GC) should reduce terminal Tm."""
        tm = calculate_adaptive_terminal_tm(14.0, 0.32, 'moderate')
        # Expected: 14 + (-0.18) * (20/0.35) = 14 - 10.3 = 3.7, but floored at 8.0
        self.assertAlmostEqual(tm, 8.0, places=1)

    def test_burkholderia_moderate(self):
        """Burkholderia (67% GC) should increase terminal Tm."""
        tm = calculate_adaptive_terminal_tm(14.0, 0.67, 'moderate')
        # Expected: 14 + 0.17 * (20/0.35) = 14 + 9.7 = 23.7
        self.assertAlmostEqual(tm, 23.7, places=1)

    def test_lenient_floor(self):
        """Lenient stringency should have 6.0°C floor."""
        tm = calculate_adaptive_terminal_tm(10.0, 0.19, 'lenient')
        self.assertGreaterEqual(tm, 6.0)

    def test_strict_floor(self):
        """Strict stringency should have 12.0°C floor."""
        tm = calculate_adaptive_terminal_tm(18.0, 0.19, 'strict')
        self.assertGreaterEqual(tm, 12.0)


class TestAdaptiveGCRange(unittest.TestCase):
    """Test adaptive GC range calculation."""

    def test_balanced_genome(self):
        """Balanced genome should center around 50%."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.50)
        self.assertAlmostEqual(min_gc, 0.35, places=2)
        self.assertAlmostEqual(max_gc, 0.65, places=2)

    def test_francisella(self):
        """Francisella (32% GC) should have AT-shifted range."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.32)
        self.assertAlmostEqual(min_gc, 0.17, places=2)
        self.assertAlmostEqual(max_gc, 0.47, places=2)

    def test_burkholderia(self):
        """Burkholderia (67% GC) should have GC-shifted range."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.67)
        self.assertAlmostEqual(min_gc, 0.52, places=2)
        self.assertAlmostEqual(max_gc, 0.82, places=2)

    def test_floor_at_15_percent(self):
        """GC range should not go below 15%."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.10)  # Very AT-rich
        self.assertGreaterEqual(min_gc, 0.15)

    def test_ceiling_at_85_percent(self):
        """GC range should not exceed 85%."""
        min_gc, max_gc = calculate_adaptive_gc_range(0.80)  # Very GC-rich
        self.assertLessEqual(max_gc, 0.85)


class TestAdaptiveDimerThreshold(unittest.TestCase):
    """Test adaptive dimer threshold calculation."""

    def test_balanced_genome(self):
        """Balanced genome should use standard threshold."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.50)
        self.assertAlmostEqual(threshold, -9.0, places=1)

    def test_francisella(self):
        """Francisella (32% GC) should have more lenient threshold."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.32)
        # Expected: -9.0 - (-0.18) * (3/0.35) = -9.0 + 1.54 = -7.46
        self.assertAlmostEqual(threshold, -7.46, places=1)
        self.assertGreater(threshold, -9.0)  # Less negative = more lenient

    def test_burkholderia(self):
        """Burkholderia (67% GC) should have stricter threshold."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.67)
        # Expected: -9.0 - 0.17 * (3/0.35) = -9.0 - 1.46 = -10.46
        self.assertAlmostEqual(threshold, -10.46, places=1)
        self.assertLess(threshold, -9.0)  # More negative = stricter

    def test_floor_at_minus_15(self):
        """Threshold should not be more negative than -15.0."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.80)
        self.assertGreaterEqual(threshold, -15.0)

    def test_ceiling_at_minus_5(self):
        """Threshold should not be less negative than -5.0."""
        threshold = calculate_adaptive_dimer_threshold(-9.0, 0.15)
        self.assertLessEqual(threshold, -5.0)


class TestGCClassification(unittest.TestCase):
    """Test genome GC classification."""

    def test_extreme_at(self):
        """Plasmodium (19% GC) should be extreme_at."""
        gc_class = get_gc_class(0.19)
        self.assertEqual(gc_class, 'extreme_at')

    def test_at_rich(self):
        """Francisella (32% GC) should be at_rich."""
        gc_class = get_gc_class(0.32)
        self.assertEqual(gc_class, 'at_rich')

    def test_balanced(self):
        """E. coli (50% GC) should be balanced."""
        gc_class = get_gc_class(0.50)
        self.assertEqual(gc_class, 'balanced')

    def test_gc_rich(self):
        """Burkholderia (67% GC) should be gc_rich."""
        gc_class = get_gc_class(0.67)
        self.assertEqual(gc_class, 'gc_rich')

    def test_extreme_gc(self):
        """Streptomyces (72% GC) should be extreme_gc."""
        gc_class = get_gc_class(0.72)
        self.assertEqual(gc_class, 'extreme_gc')

    def test_boundary_cases(self):
        """Test classification boundary cases."""
        self.assertEqual(get_gc_class(0.24), 'extreme_at')
        self.assertEqual(get_gc_class(0.25), 'at_rich')
        self.assertEqual(get_gc_class(0.40), 'balanced')
        self.assertEqual(get_gc_class(0.60), 'gc_rich')
        self.assertEqual(get_gc_class(0.70), 'extreme_gc')


class TestAdaptiveQARecommendations(unittest.TestCase):
    """Test adaptive QA recommendations."""

    def test_francisella_recommendation(self):
        """Francisella should recommend adaptive QA."""
        rec = recommend_adaptive_qa(0.32)
        self.assertTrue(rec['use_adaptive'])
        self.assertEqual(rec['gc_class'], 'at_rich')
        self.assertIn('strongly recommended', rec['reason'].lower())

    def test_balanced_recommendation(self):
        """Balanced genome should not require adaptive QA."""
        rec = recommend_adaptive_qa(0.50)
        self.assertFalse(rec['use_adaptive'])
        self.assertEqual(rec['gc_class'], 'balanced')
        self.assertIn('optional', rec['reason'].lower())

    def test_burkholderia_recommendation(self):
        """Burkholderia should recommend adaptive QA."""
        rec = recommend_adaptive_qa(0.67)
        self.assertTrue(rec['use_adaptive'])
        self.assertEqual(rec['gc_class'], 'gc_rich')
        self.assertIn('strongly recommended', rec['reason'].lower())

    def test_extreme_at_critical(self):
        """Extreme AT-rich should be marked CRITICAL."""
        rec = recommend_adaptive_qa(0.19)
        self.assertTrue(rec['use_adaptive'])
        self.assertIn('CRITICAL', rec['reason'])


class TestGenomeGCCalculation(unittest.TestCase):
    """Test genome GC content calculation from FASTA."""

    def setUp(self):
        """Create temporary FASTA files for testing."""
        self.temp_dir = tempfile.mkdtemp()

    def tearDown(self):
        """Clean up temporary files."""
        import shutil
        shutil.rmtree(self.temp_dir)

    def test_simple_fasta(self):
        """Test GC calculation on simple FASTA."""
        fasta_path = Path(self.temp_dir) / "test.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">test_seq\n")
            f.write("ATGCATGC\n")  # 50% GC (4 GC, 4 AT)

        gc = calculate_genome_gc(fasta_path)
        self.assertAlmostEqual(gc, 0.50, places=2)

    def test_multi_contig_fasta(self):
        """Test GC calculation on multi-contig FASTA."""
        fasta_path = Path(self.temp_dir) / "multi.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">contig1\n")
            f.write("GGGGCCCC\n")  # 100% GC
            f.write(">contig2\n")
            f.write("AAAATTTT\n")  # 0% GC

        gc = calculate_genome_gc(fasta_path)
        self.assertAlmostEqual(gc, 0.50, places=2)  # Average

    def test_at_rich_fasta(self):
        """Test AT-rich genome."""
        fasta_path = Path(self.temp_dir) / "at_rich.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq\n")
            f.write("AAAAATTTTTAAAAA\n")  # 0% GC

        gc = calculate_genome_gc(fasta_path)
        self.assertAlmostEqual(gc, 0.0, places=2)

    def test_gc_rich_fasta(self):
        """Test GC-rich genome."""
        fasta_path = Path(self.temp_dir) / "gc_rich.fasta"
        with open(fasta_path, 'w') as f:
            f.write(">seq\n")
            f.write("GGGGGGCCCCC\n")  # 100% GC

        gc = calculate_genome_gc(fasta_path)
        self.assertAlmostEqual(gc, 1.0, places=2)


class TestAdaptiveAnalyzerCreation(unittest.TestCase):
    """Test creation of adaptive analyzers."""

    def test_adaptive_three_prime_analyzer(self):
        """Test creation of adaptive 3' stability analyzer."""
        analyzer = create_three_prime_analyzer_adaptive(
            genome_gc=0.32,
            stringency='moderate'
        )
        self.assertIsNotNone(analyzer)
        # The terminal Tm threshold should be adjusted
        # This is internal to the analyzer, but we can test it scores correctly

    def test_adaptive_thermodynamic_filter(self):
        """Test creation of adaptive thermodynamic filter."""
        filter_obj = create_thermodynamic_filter_adaptive(
            genome_gc=0.32,
            temperature=30.0
        )
        self.assertIsNotNone(filter_obj)
        # Check that GC range was adapted
        self.assertAlmostEqual(filter_obj.criteria.min_gc, 0.17, places=2)
        self.assertAlmostEqual(filter_obj.criteria.max_gc, 0.47, places=2)

    def test_adaptive_quality_scorer(self):
        """Test creation of adaptive integrated quality scorer."""
        scorer = create_quality_scorer(
            stringency='moderate',
            genome_gc=0.32
        )
        self.assertIsNotNone(scorer)
        self.assertAlmostEqual(scorer.genome_gc, 0.32, places=2)


class TestAdaptiveQAEndToEnd(unittest.TestCase):
    """End-to-end tests of adaptive QA system."""

    def test_francisella_at_rich_primer(self):
        """Test that AT-rich primers score better with Francisella adaptive QA."""
        at_rich_primer = "ATATATATATAT"

        # Standard QA
        scorer_std = create_quality_scorer('moderate')
        score_std = scorer_std.score_primer(at_rich_primer)

        # Adaptive QA for Francisella
        scorer_adaptive = create_quality_scorer('moderate', genome_gc=0.32)
        score_adaptive = scorer_adaptive.score_primer(at_rich_primer)

        # Adaptive QA should give same or better 3' stability score for AT-rich primer
        # (because it reduces terminal Tm requirements)
        self.assertGreaterEqual(
            score_adaptive.three_prime_score,
            score_std.three_prime_score - 0.01  # Allow small numerical difference
        )

    def test_burkholderia_gc_rich_primer(self):
        """Test that GC-rich primers work with Burkholderia adaptive QA."""
        gc_rich_primer = "GCGCGCGCGCGC"

        # Standard QA
        scorer_std = create_quality_scorer('moderate')
        score_std = scorer_std.score_primer(gc_rich_primer)

        # Adaptive QA for Burkholderia
        scorer_adaptive = create_quality_scorer('moderate', genome_gc=0.67)
        score_adaptive = scorer_adaptive.score_primer(gc_rich_primer)

        # Both should work well for GC-rich primer
        self.assertGreater(score_std.three_prime_score, 0.5)
        self.assertGreater(score_adaptive.three_prime_score, 0.5)


class TestAdaptiveQAValidation(unittest.TestCase):
    """Validation tests for adaptive QA calculations."""

    def test_gc_deviation_symmetry(self):
        """Test that GC deviation is symmetric around 50%."""
        dev_at = calculate_gc_deviation(0.32)
        dev_gc = calculate_gc_deviation(0.68)
        self.assertAlmostEqual(abs(dev_at), abs(dev_gc), places=6)

    def test_terminal_tm_monotonic(self):
        """Test that terminal Tm increases monotonically with GC (when not floored)."""
        # Use strict stringency to avoid hitting floor at low GC
        tm_19 = calculate_adaptive_terminal_tm(18.0, 0.19, 'strict')
        tm_32 = calculate_adaptive_terminal_tm(18.0, 0.32, 'strict')
        tm_50 = calculate_adaptive_terminal_tm(18.0, 0.50, 'strict')
        tm_67 = calculate_adaptive_terminal_tm(18.0, 0.67, 'strict')
        tm_72 = calculate_adaptive_terminal_tm(18.0, 0.72, 'strict')

        # With strict stringency (floor = 12.0), 19% and 32% won't hit floor
        self.assertLessEqual(tm_19, tm_32)  # Use <= to allow floor
        self.assertLess(tm_32, tm_50)
        self.assertLess(tm_50, tm_67)
        self.assertLess(tm_67, tm_72)

    def test_dimer_threshold_monotonic(self):
        """Test that dimer threshold becomes stricter with higher GC."""
        dg_19 = calculate_adaptive_dimer_threshold(-9.0, 0.19)
        dg_32 = calculate_adaptive_dimer_threshold(-9.0, 0.32)
        dg_50 = calculate_adaptive_dimer_threshold(-9.0, 0.50)
        dg_67 = calculate_adaptive_dimer_threshold(-9.0, 0.67)
        dg_72 = calculate_adaptive_dimer_threshold(-9.0, 0.72)

        # More negative = stricter, should increase with GC
        self.assertGreater(dg_19, dg_32)  # Less negative for AT-rich
        self.assertGreater(dg_32, dg_50)
        self.assertGreater(dg_50, dg_67)
        self.assertGreater(dg_67, dg_72)  # More negative for GC-rich


if __name__ == '__main__':
    # Run all tests
    unittest.main(verbosity=2)
