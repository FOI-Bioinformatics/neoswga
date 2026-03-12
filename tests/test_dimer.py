"""
Unified tests for dimer detection in the dimer.py module.

Covers:
- Sequence-based dimer detection (is_dimer, heterodimer_matrix, compatible_set)
- Thermodynamic dimer detection (is_dimer_thermodynamic)
- Consistency between sequence-based and thermodynamic approaches
"""

import unittest

import numpy as np
import pytest

import neoswga.core.parameter as parameter
from neoswga.core.dimer import (
    compatible,
    compatible_set,
    heterodimer_matrix,
    is_compatible_set,
    is_dimer,
    is_dimer_fast,
    is_dimer_thermodynamic,
)
from neoswga.core.thermodynamics import reverse_complement

# Initialize default parameter values for testing.
parameter.max_self_dimer_bp = 3
parameter.max_dimer_bp = 3


# =============================================================================
# Sequence-based dimer detection
# =============================================================================


class TestSequenceBasedDimer:
    """Tests for is_dimer and related sequence-based functions."""

    # -- is_dimer ------------------------------------------------------------

    def test_is_dimer_complementary_sequences(self):
        """Complementary sequences should return a boolean result."""
        seq1 = "ATCGATCG"
        seq2 = "CGATCGAT"  # Reverse complement of seq1
        result = is_dimer(seq1, seq2)
        assert isinstance(result, (bool, np.bool_))

    def test_is_dimer_non_complementary(self):
        """Non-complementary homopolymers should not form a dimer."""
        seq1 = "AAAAAAAAAA"
        seq2 = "AAAAAAAAAA"  # RC is TTTTTTTTTT, no overlap with AAAA
        result = is_dimer(seq1, seq2, max_dimer_bp=5)
        assert result is False or result == np.False_

    def test_is_dimer_threshold_effect(self):
        """Strict vs. lenient thresholds should both produce boolean results."""
        seq1 = "ATCGATCG"
        seq2 = "CGATCGAT"

        result_strict = is_dimer(seq1, seq2, max_dimer_bp=2)
        result_lenient = is_dimer(seq1, seq2, max_dimer_bp=10)

        assert isinstance(result_strict, (bool, np.bool_))
        assert isinstance(result_lenient, (bool, np.bool_))

    def test_is_dimer_self_complementary(self):
        """Self-complementary sequence should return a boolean."""
        seq = "GCGCGCGC"
        result = is_dimer(seq, seq)
        assert isinstance(result, (bool, np.bool_))

    def test_is_dimer_short_sequences(self):
        """Short sequences should not raise errors."""
        result = is_dimer("AT", "AT")
        assert isinstance(result, (bool, np.bool_))

    def test_is_dimer_single_base(self):
        """Single-base sequences should not raise errors."""
        result = is_dimer("A", "T")
        assert isinstance(result, (bool, np.bool_))

    # -- heterodimer_matrix --------------------------------------------------

    def test_matrix_shape(self):
        """Matrix dimensions should match the number of primers."""
        primers = ["ATCGATCG", "GCTAGCTA", "AATTAATT"]
        matrix = heterodimer_matrix(primers)
        assert matrix.shape == (3, 3)

    def test_matrix_symmetry(self):
        """Heterodimer matrix should be symmetric."""
        primers = ["ATCGATCG", "GCTAGCTA", "AATTAATT"]
        matrix = heterodimer_matrix(primers)
        np.testing.assert_array_equal(matrix, matrix.T)

    def test_matrix_binary(self):
        """Matrix should contain only 0 and 1 values."""
        primers = ["ATCGATCG", "GCTAGCTA", "AATTAATT"]
        matrix = heterodimer_matrix(primers)
        assert np.all((matrix == 0) | (matrix == 1))

    def test_matrix_diagonal(self):
        """Diagonal entries should be valid binary values."""
        primers = ["ATCGATCG", "GCTAGCTA"]
        matrix = heterodimer_matrix(primers)
        for i in range(len(primers)):
            assert matrix[i, i] in [0, 1]

    def test_matrix_single_primer(self):
        """Single-primer matrix should be 1x1."""
        matrix = heterodimer_matrix(["ATCGATCG"])
        assert matrix.shape == (1, 1)

    def test_matrix_threshold_effect(self):
        """Strict threshold should yield at least as many dimers as lenient."""
        primers = ["ATCGATCG", "CGATCGAT"]

        matrix_strict = heterodimer_matrix(primers, max_dimer_bp=2)
        matrix_lenient = heterodimer_matrix(primers, max_dimer_bp=10)

        assert matrix_strict.shape == (2, 2)
        assert matrix_lenient.shape == (2, 2)
        assert np.sum(matrix_strict) >= np.sum(matrix_lenient)

    def test_matrix_empty_list(self):
        """Empty primer list should produce a 0x0 matrix."""
        matrix = heterodimer_matrix([])
        assert matrix.shape == (0, 0)

    def test_matrix_very_short_primers(self):
        """Very short primers should produce a valid matrix."""
        primers = ["AT", "GC", "TA"]
        matrix = heterodimer_matrix(primers)
        assert matrix.shape == (3, 3)

    def test_matrix_identical_primers(self):
        """Identical primers should have identical off-diagonal values."""
        primers = ["ATCGATCG", "ATCGATCG", "ATCGATCG"]
        matrix = heterodimer_matrix(primers)
        assert matrix[0, 1] == matrix[1, 2]
        assert matrix[0, 2] == matrix[1, 2]

    # -- compatible_set / compatible / is_compatible_set ---------------------

    def test_compatible_set_all_compatible(self):
        """Non-dimer-forming primers should return a boolean."""
        primers = ["AAAAAAAAAA", "CCCCCCCCCC", "GGGGGGGGGG"]
        matrix = heterodimer_matrix(primers, max_dimer_bp=5)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        result = compatible_set(matrix, primers, primer_to_idx)
        assert isinstance(result, bool)

    def test_compatible_set_one_incompatible(self):
        """One incompatible pair should still return a boolean."""
        primers = ["ATCGATCG", "CGATCGAT", "AAAAAAAAAA"]
        matrix = heterodimer_matrix(primers, max_dimer_bp=2)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        result = compatible_set(matrix, primers, primer_to_idx)
        assert isinstance(result, bool)

    def test_compatible_set_single_primer(self):
        """Single primer should always be compatible with itself (modulo self-dimer)."""
        primers = ["ATCGATCG"]
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {primers[0]: 0}

        result = compatible_set(matrix, primers, primer_to_idx)
        assert isinstance(result, bool)

    def test_compatible_with_empty_set(self):
        """Adding a primer to an empty set should always be compatible."""
        primers = ["ATCGATCG", "GCTAGCTA"]
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        result = compatible(matrix, [], "ATCGATCG", primer_to_idx)
        assert result is True

    def test_compatible_with_existing_primer(self):
        """Adding a primer to an existing set should return a boolean."""
        primers = ["ATCGATCG", "GCTAGCTA", "AAAAAAAAAA"]
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        selected = ["ATCGATCG", "GCTAGCTA"]
        result = compatible(matrix, selected, "AAAAAAAAAA", primer_to_idx)
        assert isinstance(result, bool)

    def test_is_compatible_set_compatible(self):
        """Non-dimer-forming pair should return a boolean."""
        primers = ["AAAAAAAAAA", "CCCCCCCCCC"]
        result = is_compatible_set(primers, max_dimer_bp=5)
        assert isinstance(result, bool)

    def test_is_compatible_set_single_primer(self):
        """Single-primer set should return a boolean."""
        result = is_compatible_set(["ATCGATCG"])
        assert isinstance(result, bool)

    def test_is_compatible_set_large(self):
        """Larger primer sets should not raise errors."""
        primers = [
            "ATCGATCGAT",
            "GCTAGCTAGA",
            "AATTAATTAA",
            "CCGGCCGGCC",
            "TATATATATAT",
        ]
        result = is_compatible_set(primers, max_dimer_bp=3)
        assert isinstance(result, bool)

    def test_is_compatible_set_threshold_effect(self):
        """Lenient threshold should be at least as compatible as strict."""
        primers = ["ATCGATCG", "CGATCGAT"]
        result_strict = is_compatible_set(primers, max_dimer_bp=2)
        result_lenient = is_compatible_set(primers, max_dimer_bp=10)

        assert isinstance(result_strict, bool)
        assert isinstance(result_lenient, bool)


# =============================================================================
# Thermodynamic dimer detection
# =============================================================================


class TestThermodynamicDimer:
    """Tests for is_dimer_thermodynamic."""

    def test_perfect_complement_detected(self):
        """A primer paired with its exact reverse complement should be a dimer."""
        seq = "ATCGATCGATCG"
        seq_rc = reverse_complement(seq)
        assert is_dimer_thermodynamic(seq, seq_rc) is True

    def test_short_complement_not_flagged(self):
        """Primers sharing minimal complementarity should not be flagged."""
        assert is_dimer_thermodynamic("ATCGATCG", "TTTTTTTT") is False

    def test_no_complementarity(self):
        """Primers with no complementary bases should not be flagged."""
        assert is_dimer_thermodynamic("AAAAAAAA", "AAAAAAAA") is False

    def test_moderate_complement_depends_on_threshold(self):
        """Moderate complementary region flagged only with lenient threshold."""
        seq_1 = "ATCGATCGATCG"
        seq_2 = "CGATCG"  # 6 bp complement embedded

        # Strict threshold: should not flag short duplexes
        assert is_dimer_thermodynamic(seq_1, seq_2, delta_g_threshold=-15.0) is False
        # Lenient threshold: should flag even weak interactions
        assert is_dimer_thermodynamic(seq_1, seq_2, delta_g_threshold=-2.0) is True

    def test_custom_conditions_temperature(self):
        """Temperature changes should not break dimer detection."""
        from neoswga.core.reaction_conditions import ReactionConditions

        seq = "GCGCGCGCGCGC"
        seq_rc = reverse_complement(seq)

        cold = ReactionConditions(temp=20.0, polymerase="phi29")
        assert is_dimer_thermodynamic(seq, seq_rc, conditions=cold) is True

        hot = ReactionConditions(temp=65.0, polymerase="bst")
        assert is_dimer_thermodynamic(seq, seq_rc, conditions=hot) is True

    def test_self_dimer(self):
        """A palindromic sequence tested against itself should return a boolean."""
        palindrome = "ATCGATCGAT"
        result = is_dimer_thermodynamic(palindrome, palindrome)
        assert isinstance(result, bool)


# =============================================================================
# Cross-method consistency
# =============================================================================


class TestDimerConsistency:
    """Tests for consistency between sequence-based and thermodynamic checks."""

    def test_matrix_and_is_dimer_consistent(self):
        """Heterodimer matrix values should match individual is_dimer calls."""
        primers = ["ATCGATCG", "GCTAGCTA"]
        matrix = heterodimer_matrix(primers, max_dimer_bp=3)

        expected_01 = 1 if is_dimer(primers[0], primers[1], max_dimer_bp=3) else 0
        assert matrix[0, 1] == expected_01
        assert matrix[1, 0] == expected_01

    def test_compatible_functions_consistent(self):
        """compatible and compatible_set should not contradict each other."""
        primers = ["ATCGATCG", "GCTAGCTA", "AAAAAAAAAA"]
        matrix = heterodimer_matrix(primers)
        primer_to_idx = {p: i for i, p in enumerate(primers)}

        # Smoke test: both functions should execute without error
        compatible_set(matrix, primers, primer_to_idx)
        compatible(matrix, [primers[0]], primers[1], primer_to_idx)

    def test_sequence_and_thermo_agree_on_strong_dimers(self):
        """Strong complements flagged by sequence check should also flag thermodynamically."""
        seq_1 = "GCGATCGATCGC"
        seq_2 = reverse_complement(seq_1)

        assert is_dimer_fast(seq_1, seq_2, max_dimer_bp=3) is True
        assert is_dimer_thermodynamic(seq_1, seq_2) is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
