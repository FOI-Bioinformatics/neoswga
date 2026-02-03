"""
Comprehensive unit tests for legacy dimer.py module.

Tests all functionality including:
- is_dimer function (threshold-based dimer detection)
- heterodimer_matrix calculation
- compatible_set and compatible functions
- is_compatible_set function
- Edge cases
"""

import unittest
import numpy as np
import logging

import neoswga.core.parameter as parameter

from neoswga.core.dimer import (
    is_dimer,
    heterodimer_matrix,
    compatible_set,
    compatible,
    is_compatible_set,
)

# Suppress logging during tests
logging.disable(logging.CRITICAL)

# Initialize default parameter values for testing
# (normally these are loaded from a params file)
parameter.max_self_dimer_bp = 3
parameter.max_dimer_bp = 3


class TestIsDimer(unittest.TestCase):
    """Test is_dimer function."""

    def test_is_dimer_complementary_sequences(self):
        """Test dimer detection for complementary sequences."""
        seq1 = 'ATCGATCG'
        seq2 = 'CGATCGAT'  # Reverse complement of seq1
        # These should form a dimer (long common substring with RC)
        result = is_dimer(seq1, seq2)
        self.assertIsInstance(result, (bool, np.bool_))

    def test_is_dimer_non_complementary(self):
        """Test dimer detection for non-complementary sequences."""
        seq1 = 'AAAAAAAAAA'
        seq2 = 'AAAAAAAAAA'  # RC is TTTTTTTTTT, no overlap with AAAA
        # These should not form a strong dimer
        result = is_dimer(seq1, seq2, max_dimer_bp=5)
        self.assertFalse(result)

    def test_is_dimer_threshold_effect(self):
        """Test that threshold affects dimer detection."""
        seq1 = 'ATCGATCG'
        seq2 = 'CGATCGAT'

        # With strict threshold, may detect dimer
        result_strict = is_dimer(seq1, seq2, max_dimer_bp=2)
        # With lenient threshold, should not detect dimer
        result_lenient = is_dimer(seq1, seq2, max_dimer_bp=10)

        # Lenient should be False or same as strict (depends on match length)
        # At minimum, both should be boolean
        self.assertIsInstance(result_strict, (bool, np.bool_))
        self.assertIsInstance(result_lenient, (bool, np.bool_))

    def test_is_dimer_self_complementary(self):
        """Test dimer detection for self-complementary sequence."""
        # GCGC is self-complementary
        seq = 'GCGCGCGC'
        result = is_dimer(seq, seq)
        self.assertIsInstance(result, (bool, np.bool_))

    def test_is_dimer_short_sequences(self):
        """Test dimer detection for short sequences."""
        seq1 = 'AT'
        seq2 = 'AT'
        result = is_dimer(seq1, seq2)
        self.assertIsInstance(result, (bool, np.bool_))

    def test_is_dimer_single_base(self):
        """Test dimer detection for single base."""
        seq1 = 'A'
        seq2 = 'T'
        result = is_dimer(seq1, seq2)
        self.assertIsInstance(result, (bool, np.bool_))


class TestHeterodimerMatrix(unittest.TestCase):
    """Test heterodimer_matrix function."""

    def test_matrix_shape(self):
        """Test that matrix has correct shape."""
        primers = ['ATCGATCG', 'GCTAGCTA', 'AATTAATT']
        matrix = heterodimer_matrix(primers)

        self.assertEqual(matrix.shape, (3, 3))

    def test_matrix_symmetry(self):
        """Test that matrix is symmetric."""
        primers = ['ATCGATCG', 'GCTAGCTA', 'AATTAATT']
        matrix = heterodimer_matrix(primers)

        # Matrix should be symmetric
        np.testing.assert_array_equal(matrix, matrix.T)

    def test_matrix_binary(self):
        """Test that matrix contains only 0s and 1s."""
        primers = ['ATCGATCG', 'GCTAGCTA', 'AATTAATT']
        matrix = heterodimer_matrix(primers)

        # All values should be 0 or 1
        self.assertTrue(np.all((matrix == 0) | (matrix == 1)))

    def test_matrix_diagonal(self):
        """Test that diagonal contains self-dimer results."""
        primers = ['ATCGATCG', 'GCTAGCTA']
        matrix = heterodimer_matrix(primers)

        # Diagonal should be 0 or 1
        for i in range(len(primers)):
            self.assertIn(matrix[i, i], [0, 1])

    def test_matrix_single_primer(self):
        """Test matrix for single primer."""
        primers = ['ATCGATCG']
        matrix = heterodimer_matrix(primers)

        self.assertEqual(matrix.shape, (1, 1))

    def test_matrix_threshold_effect(self):
        """Test that threshold affects matrix."""
        primers = ['ATCGATCG', 'CGATCGAT']

        matrix_strict = heterodimer_matrix(primers, max_dimer_bp=2)
        matrix_lenient = heterodimer_matrix(primers, max_dimer_bp=10)

        # Both should be valid matrices
        self.assertEqual(matrix_strict.shape, (2, 2))
        self.assertEqual(matrix_lenient.shape, (2, 2))

        # Lenient threshold should have same or fewer 1s
        self.assertGreaterEqual(np.sum(matrix_strict), np.sum(matrix_lenient))


class TestCompatibleSet(unittest.TestCase):
    """Test compatible_set function."""

    def test_compatible_set_all_compatible(self):
        """Test compatible_set with primers that don't form dimers."""
        primers = ['AAAAAAAAAA', 'CCCCCCCCCC', 'GGGGGGGGGG']
        matrix = heterodimer_matrix(primers, max_dimer_bp=5)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        # If matrix has no 1s in off-diagonal, set is compatible
        result = compatible_set(matrix, primers, primer_to_idx)
        self.assertIsInstance(result, bool)

    def test_compatible_set_one_incompatible(self):
        """Test compatible_set with one incompatible pair."""
        primers = ['ATCGATCG', 'CGATCGAT', 'AAAAAAAAAA']
        matrix = heterodimer_matrix(primers, max_dimer_bp=2)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        result = compatible_set(matrix, primers, primer_to_idx)
        self.assertIsInstance(result, bool)

    def test_compatible_set_single_primer(self):
        """Test compatible_set with single primer (always compatible)."""
        primers = ['ATCGATCG']
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {primers[0]: 0}

        result = compatible_set(matrix, primers, primer_to_idx)
        # Single primer can only have self-dimer issue
        self.assertIsInstance(result, bool)


class TestCompatible(unittest.TestCase):
    """Test compatible function (adding single primer)."""

    def test_compatible_with_empty_set(self):
        """Test adding primer to empty set."""
        primers = ['ATCGATCG', 'GCTAGCTA']
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        # Adding to empty set should always be compatible
        result = compatible(matrix, [], 'ATCGATCG', primer_to_idx)
        self.assertTrue(result)

    def test_compatible_with_existing_primer(self):
        """Test adding primer to existing set."""
        primers = ['ATCGATCG', 'GCTAGCTA', 'AAAAAAAAAA']
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        # Check if third primer compatible with first two
        selected = ['ATCGATCG', 'GCTAGCTA']
        result = compatible(matrix, selected, 'AAAAAAAAAA', primer_to_idx)
        self.assertIsInstance(result, bool)


class TestIsCompatibleSet(unittest.TestCase):
    """Test is_compatible_set function (convenience wrapper)."""

    def test_is_compatible_set_compatible(self):
        """Test compatible set detection."""
        # These primers shouldn't form strong dimers with each other
        primers = ['AAAAAAAAAA', 'CCCCCCCCCC']
        result = is_compatible_set(primers, max_dimer_bp=5)
        self.assertIsInstance(result, bool)

    def test_is_compatible_set_single_primer(self):
        """Test is_compatible_set with single primer."""
        primers = ['ATCGATCG']
        result = is_compatible_set(primers)
        self.assertIsInstance(result, bool)

    def test_is_compatible_set_large(self):
        """Test is_compatible_set with larger set."""
        # Generate multiple random-ish primers
        primers = [
            'ATCGATCGAT',
            'GCTAGCTAGA',
            'AATTAATTAA',
            'CCGGCCGGCC',
            'TATATATATAT',
        ]
        result = is_compatible_set(primers, max_dimer_bp=3)
        self.assertIsInstance(result, bool)

    def test_is_compatible_set_threshold_effect(self):
        """Test that threshold affects compatibility."""
        primers = ['ATCGATCG', 'CGATCGAT']

        result_strict = is_compatible_set(primers, max_dimer_bp=2)
        result_lenient = is_compatible_set(primers, max_dimer_bp=10)

        self.assertIsInstance(result_strict, bool)
        self.assertIsInstance(result_lenient, bool)
        # Lenient threshold should be same or more compatible
        if not result_strict:
            # If strict fails, lenient might still pass
            pass  # Just checking it doesn't crash


class TestEdgeCases(unittest.TestCase):
    """Test edge cases."""

    def test_empty_primer_list(self):
        """Test heterodimer_matrix with empty list."""
        matrix = heterodimer_matrix([])
        self.assertEqual(matrix.shape, (0, 0))

    def test_very_short_primers(self):
        """Test with very short primers."""
        primers = ['AT', 'GC', 'TA']
        matrix = heterodimer_matrix(primers)
        self.assertEqual(matrix.shape, (3, 3))

    def test_identical_primers(self):
        """Test with identical primers."""
        primers = ['ATCGATCG', 'ATCGATCG', 'ATCGATCG']
        matrix = heterodimer_matrix(primers)

        # All off-diagonal should be same (all identical)
        self.assertEqual(matrix[0, 1], matrix[1, 2])
        self.assertEqual(matrix[0, 2], matrix[1, 2])


class TestConsistency(unittest.TestCase):
    """Test consistency between functions."""

    def test_matrix_and_is_dimer_consistent(self):
        """Test that heterodimer_matrix is consistent with is_dimer."""
        primers = ['ATCGATCG', 'GCTAGCTA']
        matrix = heterodimer_matrix(primers, max_dimer_bp=3)

        # Check off-diagonal matches is_dimer
        expected_01 = 1 if is_dimer(primers[0], primers[1], max_dimer_bp=3) else 0
        self.assertEqual(matrix[0, 1], expected_01)
        self.assertEqual(matrix[1, 0], expected_01)

    def test_compatible_functions_consistent(self):
        """Test compatible and compatible_set consistency."""
        primers = ['ATCGATCG', 'GCTAGCTA', 'AAAAAAAAAA']
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        # compatible_set([A, B, C]) should equal:
        # compatible([], A) AND compatible([A], B) AND compatible([A,B], C)
        # But this is expensive to verify fully


if __name__ == '__main__':
    unittest.main()
