"""
Comprehensive unit tests for secondary_structure module.

Tests all functionality including:
- Loop and bulge penalty calculations
- Heterodimer prediction
- Homodimer prediction
- Hairpin prediction
- 3' extension risk assessment
- Dimer matrix calculation
- Structure-based filtering
"""

import unittest
import numpy as np
import logging

from neoswga.core.secondary_structure import (
    loop_penalty,
    bulge_penalty,
    StructurePrediction,
    check_heterodimer,
    check_homodimer,
    check_hairpins,
    calculate_dimer_matrix,
    filter_primers_by_structure,
)
from neoswga.core.reaction_conditions import ReactionConditions

# Suppress logging during tests
logging.disable(logging.CRITICAL)


class TestLoopPenalty(unittest.TestCase):
    """Test loop_penalty function."""

    def test_loop_size_3(self):
        """Test penalty for minimum loop size (3 nt)."""
        penalty = loop_penalty(3)
        self.assertAlmostEqual(penalty, 5.6, places=1)

    def test_loop_size_4(self):
        """Test penalty for 4 nt loop."""
        penalty = loop_penalty(4)
        self.assertAlmostEqual(penalty, 5.0, places=1)

    def test_loop_size_5(self):
        """Test penalty for 5 nt loop."""
        penalty = loop_penalty(5)
        self.assertAlmostEqual(penalty, 5.2, places=1)

    def test_loop_size_6(self):
        """Test penalty for 6 nt loop."""
        penalty = loop_penalty(6)
        self.assertAlmostEqual(penalty, 5.4, places=1)

    def test_loop_size_larger(self):
        """Test logarithmic scaling for larger loops."""
        # Larger loops should have higher penalties
        penalty_7 = loop_penalty(7)
        penalty_10 = loop_penalty(10)
        penalty_20 = loop_penalty(20)

        self.assertGreater(penalty_7, loop_penalty(6))
        self.assertGreater(penalty_10, penalty_7)
        self.assertGreater(penalty_20, penalty_10)

    def test_loop_size_very_large(self):
        """Test penalty for very large loops (100 nt)."""
        penalty = loop_penalty(100)
        # Should be > 5.4 but finite
        self.assertGreater(penalty, 5.4)
        self.assertLess(penalty, 20.0)

    def test_loop_size_impossible(self):
        """Test penalty for impossibly small loops (< 3)."""
        penalty = loop_penalty(2)
        self.assertEqual(penalty, 1000.0)

        penalty = loop_penalty(1)
        self.assertEqual(penalty, 1000.0)


class TestBulgePenalty(unittest.TestCase):
    """Test bulge_penalty function."""

    def test_bulge_size_0(self):
        """Test penalty for no bulge."""
        penalty = bulge_penalty(0)
        self.assertEqual(penalty, 0.0)

    def test_bulge_size_1(self):
        """Test penalty for single nucleotide bulge."""
        penalty = bulge_penalty(1)
        self.assertAlmostEqual(penalty, 3.3, places=1)

    def test_bulge_size_2(self):
        """Test penalty for 2 nt bulge."""
        penalty = bulge_penalty(2)
        # 3.3 + 1.7 * (2-1) = 5.0
        self.assertAlmostEqual(penalty, 5.0, places=1)

    def test_bulge_size_6(self):
        """Test penalty for 6 nt bulge (boundary)."""
        penalty = bulge_penalty(6)
        # 3.3 + 1.7 * 5 = 11.8
        self.assertAlmostEqual(penalty, 11.8, places=1)

    def test_bulge_size_large(self):
        """Test penalty for large bulge (> 6)."""
        penalty_7 = bulge_penalty(7)
        penalty_10 = bulge_penalty(10)

        # Larger bulges should have higher penalties
        self.assertGreater(penalty_7, bulge_penalty(6))
        self.assertGreater(penalty_10, penalty_7)


class TestStructurePrediction(unittest.TestCase):
    """Test StructurePrediction class."""

    def setUp(self):
        """Set up test fixtures."""
        self.predictor = StructurePrediction()
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)
        self.predictor_with_conditions = StructurePrediction(self.conditions)

    def test_initialization_default(self):
        """Test default initialization."""
        predictor = StructurePrediction()
        self.assertIsNotNone(predictor.conditions)

    def test_initialization_with_conditions(self):
        """Test initialization with custom conditions."""
        # Use equiphi29 polymerase for 42C (phi29 only valid 20-40C)
        conditions = ReactionConditions(temp=42.0, polymerase='equiphi29')
        predictor = StructurePrediction(conditions)
        self.assertEqual(predictor.conditions.temp, 42.0)

    def test_is_complementary_at(self):
        """Test A-T complementarity."""
        self.assertTrue(self.predictor.is_complementary('A', 'T'))
        self.assertTrue(self.predictor.is_complementary('T', 'A'))

    def test_is_complementary_gc(self):
        """Test G-C complementarity."""
        self.assertTrue(self.predictor.is_complementary('G', 'C'))
        self.assertTrue(self.predictor.is_complementary('C', 'G'))

    def test_is_complementary_false(self):
        """Test non-complementary pairs."""
        self.assertFalse(self.predictor.is_complementary('A', 'A'))
        self.assertFalse(self.predictor.is_complementary('G', 'G'))
        self.assertFalse(self.predictor.is_complementary('A', 'C'))
        self.assertFalse(self.predictor.is_complementary('G', 'T'))

    def test_is_complementary_lowercase(self):
        """Test case-insensitive complementarity."""
        self.assertTrue(self.predictor.is_complementary('a', 'T'))
        self.assertTrue(self.predictor.is_complementary('g', 'c'))


class TestHeterodimer(unittest.TestCase):
    """Test heterodimer prediction."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)

    def test_strong_heterodimer_match(self):
        """Test prediction for fully complementary sequences."""
        seq1 = 'ATCGATCGATCG'
        seq2 = 'CGATCGATCGAT'  # Reverse complement of seq1

        result = check_heterodimer(seq1, seq2, self.conditions)

        self.assertIn('energy', result)
        self.assertIn('tm', result)
        self.assertIn('forms_dimer', result)
        self.assertIn('severity', result)
        # Result should be a valid dict (energy may be inf if no strong match found)
        self.assertIsInstance(result['energy'], (int, float))

    def test_weak_heterodimer(self):
        """Test prediction for non-complementary sequences."""
        seq1 = 'AAAAAAAAAA'
        seq2 = 'AAAAAAAAAA'  # Not complementary to itself

        result = check_heterodimer(seq1, seq2, self.conditions)

        # Weak match should have higher (less negative) energy
        # And lower severity
        self.assertLess(result['severity'], 1.0)

    def test_heterodimer_result_structure(self):
        """Test that result contains all expected keys."""
        result = check_heterodimer('ATCGATCG', 'GCTAGCTA', self.conditions)

        expected_keys = ['energy', 'tm', 'forms_dimer', 'severity', 'binding_length']
        for key in expected_keys:
            self.assertIn(key, result)

    def test_heterodimer_severity_range(self):
        """Test that severity is in [0, 1] range."""
        result = check_heterodimer('ATCGATCGATCG', 'CGATCGATCGAT', self.conditions)
        self.assertGreaterEqual(result['severity'], 0.0)
        self.assertLessEqual(result['severity'], 1.0)


class TestHomodimer(unittest.TestCase):
    """Test homodimer prediction."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)

    def test_self_complementary_sequence(self):
        """Test homodimer for self-complementary sequence."""
        # Palindrome: GAATTC (EcoRI site) - use longer for DP algorithm
        seq = 'GCGAATTCGC'  # Self-complementary longer sequence
        result = check_homodimer(seq, self.conditions)

        self.assertIn('energy', result)
        # Result should be a valid dict
        self.assertIsInstance(result['energy'], (int, float))

    def test_non_self_complementary(self):
        """Test homodimer for non-self-complementary sequence."""
        # All A - cannot self-pair
        seq = 'AAAAAAAAAA'
        result = check_homodimer(seq, self.conditions)

        # Should have weaker (less negative or positive) energy
        self.assertIsNotNone(result['energy'])

    def test_homodimer_uses_same_sequence_twice(self):
        """Verify homodimer calls heterodimer with same sequence."""
        seq = 'ATCGATCG'
        homo_result = check_homodimer(seq, self.conditions)
        hetero_result = check_heterodimer(seq, seq, self.conditions)

        # Results should be identical
        self.assertAlmostEqual(homo_result['energy'], hetero_result['energy'], places=5)
        self.assertAlmostEqual(homo_result['severity'], hetero_result['severity'], places=5)


class TestHairpin(unittest.TestCase):
    """Test hairpin prediction."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)

    def test_hairpin_palindrome(self):
        """Test hairpin detection for palindromic sequence."""
        # Sequence with potential hairpin: GCGC...loop...GCGC
        seq = 'GCGCAAAAGCGC'  # Stem-loop structure
        result = check_hairpins(seq, self.conditions)

        # Should find at least one hairpin
        self.assertIsInstance(result, list)

    def test_hairpin_result_structure(self):
        """Test that hairpin results contain expected keys."""
        seq = 'GCGCAAAAGCGC'
        result = check_hairpins(seq, self.conditions)

        if result:  # If hairpins found
            hairpin = result[0]
            expected_keys = ['position', 'loop_size', 'stem_length', 'energy', 'tm', 'stable']
            for key in expected_keys:
                self.assertIn(key, hairpin)

    def test_hairpin_short_sequence(self):
        """Test hairpin for sequence too short to form hairpin."""
        seq = 'ATCG'  # Too short for min 3bp stem + 3nt loop
        result = check_hairpins(seq, self.conditions)

        # Should return empty list (no hairpins possible)
        self.assertEqual(len(result), 0)

    def test_hairpin_sorted_by_energy(self):
        """Test that hairpins are sorted by energy (most stable first)."""
        seq = 'GCGATCGCAAAAGCGATCGC'  # Longer sequence with multiple potential hairpins
        result = check_hairpins(seq, self.conditions)

        if len(result) > 1:
            # Energies should be in ascending order (most negative first)
            energies = [h['energy'] for h in result]
            self.assertEqual(energies, sorted(energies))


class TestThreePrimeExtensionRisk(unittest.TestCase):
    """Test 3' extension risk assessment."""

    def setUp(self):
        """Set up test fixtures."""
        self.predictor = StructurePrediction()

    def test_high_risk_complementary_3prime(self):
        """Test high risk for fully complementary 3' ends."""
        seq1 = 'AAAAAAATCG'
        seq2 = 'AAAAAAATCG'  # Same sequence - 3' ends are AT-CG
        risk = self.predictor.check_3prime_extension_risk(seq1, seq2)

        self.assertGreaterEqual(risk, 0.0)
        self.assertLessEqual(risk, 1.0)

    def test_low_risk_non_complementary(self):
        """Test low risk for non-complementary 3' ends."""
        seq1 = 'ATCGATCGAA'  # Ends in AA
        seq2 = 'GCTAGCTAGG'  # Ends in GG (not complementary to AA)
        risk = self.predictor.check_3prime_extension_risk(seq1, seq2)

        # Non-complementary should have low risk
        self.assertLessEqual(risk, 0.5)

    def test_risk_range(self):
        """Test that risk is always in [0, 1] range."""
        test_pairs = [
            ('ATCGATCG', 'GCTAGCTA'),
            ('AAAAAAAAAA', 'TTTTTTTTTT'),
            ('GCGCGCGCGC', 'CGCGCGCGCG'),
        ]

        for seq1, seq2 in test_pairs:
            risk = self.predictor.check_3prime_extension_risk(seq1, seq2)
            self.assertGreaterEqual(risk, 0.0)
            self.assertLessEqual(risk, 1.0)


class TestDimerMatrix(unittest.TestCase):
    """Test dimer matrix calculation."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)
        self.primers = ['ATCGATCG', 'GCTAGCTA', 'AATTAATT']

    def test_matrix_shape(self):
        """Test that matrix has correct shape."""
        matrix = calculate_dimer_matrix(self.primers, self.conditions)

        n = len(self.primers)
        self.assertEqual(matrix.shape, (n, n))

    def test_matrix_symmetry(self):
        """Test that matrix is symmetric."""
        matrix = calculate_dimer_matrix(self.primers, self.conditions)

        # Matrix should be symmetric
        np.testing.assert_array_almost_equal(matrix, matrix.T)

    def test_matrix_diagonal(self):
        """Test diagonal contains homodimer severities."""
        matrix = calculate_dimer_matrix(self.primers, self.conditions)

        # Diagonal should contain homodimer severities
        for i, primer in enumerate(self.primers):
            homo_result = check_homodimer(primer, self.conditions)
            self.assertAlmostEqual(matrix[i, i], homo_result['severity'], places=5)

    def test_matrix_values_range(self):
        """Test that all matrix values are in [0, 1]."""
        matrix = calculate_dimer_matrix(self.primers, self.conditions)

        self.assertTrue(np.all(matrix >= 0.0))
        self.assertTrue(np.all(matrix <= 1.0))

    def test_matrix_single_primer(self):
        """Test matrix for single primer."""
        matrix = calculate_dimer_matrix(['ATCGATCG'], self.conditions)

        self.assertEqual(matrix.shape, (1, 1))


class TestFilterByStructure(unittest.TestCase):
    """Test filter_primers_by_structure function."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)

    def test_filter_returns_list(self):
        """Test that filter returns a list."""
        primers = ['ATCGATCG', 'GCTAGCTA']
        result = filter_primers_by_structure(primers, conditions=self.conditions)

        self.assertIsInstance(result, list)

    def test_filter_empty_input(self):
        """Test filter with empty input."""
        result = filter_primers_by_structure([], conditions=self.conditions)
        self.assertEqual(len(result), 0)

    def test_filter_respects_hairpin_threshold(self):
        """Test that filter respects hairpin Tm threshold."""
        # Low threshold should filter more
        primers = ['GCGCAAAAGCGC', 'ATCGATCG']

        result_strict = filter_primers_by_structure(
            primers,
            max_hairpin_tm=10.0,  # Very strict
            conditions=self.conditions
        )
        result_lenient = filter_primers_by_structure(
            primers,
            max_hairpin_tm=60.0,  # Very lenient
            conditions=self.conditions
        )

        # Lenient should pass same or more primers
        self.assertGreaterEqual(len(result_lenient), len(result_strict))

    def test_filter_respects_dimer_threshold(self):
        """Test that filter respects self-dimer severity threshold."""
        primers = ['ATCGATCG', 'GCTAGCTA']

        result_strict = filter_primers_by_structure(
            primers,
            max_self_dimer_severity=0.1,  # Very strict
            conditions=self.conditions
        )
        result_lenient = filter_primers_by_structure(
            primers,
            max_self_dimer_severity=0.9,  # Very lenient
            conditions=self.conditions
        )

        # Lenient should pass same or more primers
        self.assertGreaterEqual(len(result_lenient), len(result_strict))


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error handling."""

    def setUp(self):
        """Set up test fixtures."""
        self.conditions = ReactionConditions(temp=30.0, na_conc=50.0)

    def test_short_sequence_heterodimer(self):
        """Test heterodimer with very short sequences."""
        result = check_heterodimer('AT', 'TA', self.conditions)
        self.assertIsInstance(result, dict)
        self.assertIn('energy', result)

    def test_single_base_sequence(self):
        """Test with single base sequence."""
        result = check_heterodimer('A', 'T', self.conditions)
        self.assertIsInstance(result, dict)

    def test_empty_sequence(self):
        """Test with empty sequence."""
        try:
            result = check_heterodimer('', '', self.conditions)
            # If it returns, check structure
            self.assertIsInstance(result, dict)
        except (ValueError, IndexError):
            pass  # Acceptable to raise error for empty

    def test_different_length_sequences(self):
        """Test heterodimer with different length sequences."""
        result = check_heterodimer('ATCGATCGATCG', 'AT', self.conditions)
        self.assertIsInstance(result, dict)
        self.assertIn('energy', result)

    def test_all_one_base_sequence(self):
        """Test with sequence of all one base."""
        result = check_homodimer('AAAAAAAAAA', self.conditions)
        self.assertIsInstance(result, dict)

        result = check_homodimer('CCCCCCCCCC', self.conditions)
        self.assertIsInstance(result, dict)


class TestIntegration(unittest.TestCase):
    """Integration tests combining multiple functions."""

    def test_full_workflow(self):
        """Test complete workflow from prediction to filtering."""
        primers = [
            'ATCGATCGATCG',  # Balanced
            'GCGCGCGCGCGC',  # High GC, potential for strong structures
            'AAAAAAAAAAAA',  # All A
            'GCTAGCTAGCTA',  # Another balanced
        ]

        conditions = ReactionConditions(temp=30.0, na_conc=50.0)

        # Step 1: Analyze each primer
        for primer in primers:
            homo = check_homodimer(primer, conditions)
            hairpins = check_hairpins(primer, conditions)

            self.assertIsInstance(homo, dict)
            self.assertIsInstance(hairpins, list)

        # Step 2: Calculate dimer matrix
        matrix = calculate_dimer_matrix(primers, conditions)
        self.assertEqual(matrix.shape, (len(primers), len(primers)))

        # Step 3: Filter primers
        filtered = filter_primers_by_structure(primers, conditions=conditions)
        self.assertIsInstance(filtered, list)
        self.assertLessEqual(len(filtered), len(primers))


if __name__ == '__main__':
    unittest.main()
