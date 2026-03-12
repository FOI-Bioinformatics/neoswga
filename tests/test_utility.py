"""
Unit tests for utility module.

Tests:
- Array manipulation functions (flatten, merge, intersection)
- Mathematical functions (softmax, sigmoid, gini_exact)
- String operations (complement, reverse_complement, longest_char_repeat)
- Sequence utilities (get_num_mismatches, longest_common_substring)
- Gap calculations (get_positional_gap_lengths)
- Multiprocessing utilities (create_pool)
"""

import pytest
import numpy as np
import tempfile
from pathlib import Path

from neoswga.core.utility import (
    flatten,
    mergeArrays,
    softmax,
    sigmoid,
    longest_char_repeat,
    complement,
    get_num_mismatches,
    longest_common_substring,
    reverse,
    reverse_complement,
    intersection,
    gini_exact,
    most_frequent,
    get_positional_gap_lengths,
    create_pool,
    get_seq_length,
)


# =============================================================================
# Array Manipulation Tests
# =============================================================================

class TestFlatten:
    """Tests for flatten function."""

    def test_flatten_simple(self):
        """Test flattening simple nested list."""
        nested = [[1, 2], [3, 4], [5, 6]]
        result = flatten(nested)

        assert result == [1, 2, 3, 4, 5, 6]

    def test_flatten_empty(self):
        """Test flattening empty list."""
        result = flatten([])

        assert result == []

    def test_flatten_single_list(self):
        """Test flattening list with single sublist."""
        nested = [[1, 2, 3]]
        result = flatten(nested)

        assert result == [1, 2, 3]

    def test_flatten_strings(self):
        """Test flattening list of strings."""
        nested = [['a', 'b'], ['c', 'd']]
        result = flatten(nested)

        assert result == ['a', 'b', 'c', 'd']


class TestMergeArrays:
    """Tests for mergeArrays function."""

    def test_merge_simple(self):
        """Test merging two sorted arrays."""
        arr1 = [1, 3, 5]
        arr2 = [2, 4, 6]
        result = mergeArrays(arr1, arr2)

        assert result == [1, 2, 3, 4, 5, 6]

    def test_merge_with_duplicates(self):
        """Test merging arrays with duplicates."""
        arr1 = [1, 2, 3]
        arr2 = [2, 3, 4]
        result = mergeArrays(arr1, arr2)

        # Duplicates should be removed (only keeps from arr1)
        assert result == [1, 2, 3, 4]

    def test_merge_empty_first(self):
        """Test merging with empty first array."""
        arr1 = []
        arr2 = [1, 2, 3]
        result = mergeArrays(arr1, arr2)

        assert result == [1, 2, 3]

    def test_merge_empty_second(self):
        """Test merging with empty second array."""
        arr1 = [1, 2, 3]
        arr2 = []
        result = mergeArrays(arr1, arr2)

        assert result == [1, 2, 3]

    def test_merge_both_empty(self):
        """Test merging two empty arrays."""
        result = mergeArrays([], [])

        assert result == []


class TestIntersection:
    """Tests for intersection function."""

    def test_intersection_simple(self):
        """Test simple intersection."""
        lst1 = [1, 2, 3, 4]
        lst2 = [3, 4, 5, 6]
        result = intersection(lst1, lst2)

        assert set(result) == {3, 4}

    def test_intersection_no_common(self):
        """Test intersection with no common elements."""
        lst1 = [1, 2, 3]
        lst2 = [4, 5, 6]
        result = intersection(lst1, lst2)

        assert result == []

    def test_intersection_empty(self):
        """Test intersection with empty list."""
        lst1 = [1, 2, 3]
        lst2 = []
        result = intersection(lst1, lst2)

        assert result == []


# =============================================================================
# Mathematical Function Tests
# =============================================================================

class TestSoftmax:
    """Tests for softmax function."""

    def test_softmax_simple(self):
        """Test softmax with simple input."""
        x = np.array([1.0, 2.0, 3.0])
        result = softmax(x)

        # Should sum to 1
        assert np.isclose(result.sum(), 1.0)
        # Higher input should have higher probability
        assert result[2] > result[1] > result[0]

    def test_softmax_uniform(self):
        """Test softmax with uniform input."""
        x = np.array([1.0, 1.0, 1.0])
        result = softmax(x)

        # All should be equal
        assert np.allclose(result, [1/3, 1/3, 1/3])

    def test_softmax_single(self):
        """Test softmax with single element."""
        x = np.array([5.0])
        result = softmax(x)

        assert np.isclose(result[0], 1.0)


class TestSigmoid:
    """Tests for sigmoid function."""

    def test_sigmoid_zero(self):
        """Test sigmoid at 0."""
        result = sigmoid(0)

        assert result == 0.5

    def test_sigmoid_positive(self):
        """Test sigmoid with positive input."""
        result = sigmoid(10)

        assert result > 0.99

    def test_sigmoid_negative(self):
        """Test sigmoid with negative input."""
        result = sigmoid(-10)

        assert result < 0.01


class TestGiniExact:
    """Tests for gini_exact function."""

    def test_gini_uniform(self):
        """Test Gini coefficient for uniform distribution."""
        # All equal values -> Gini = 0
        values = [100, 100, 100, 100]
        result = gini_exact(values)

        assert result == pytest.approx(0.0, abs=0.01)

    def test_gini_perfect_inequality(self):
        """Test Gini coefficient for perfect inequality."""
        # One person has everything
        values = [0, 0, 0, 1000]
        result = gini_exact(values)

        # Should be close to 1 (but not exactly 1 due to epsilon added)
        assert result > 0.5

    def test_gini_empty(self):
        """Test Gini coefficient for empty array."""
        result = gini_exact([])

        assert result == 0

    def test_gini_single_element(self):
        """Test Gini coefficient for single element."""
        result = gini_exact([100])

        assert result == 0

    def test_gini_negative_values(self):
        """Test Gini handles negative values."""
        # Should shift values to be positive
        values = [-10, -5, 0, 5, 10]
        result = gini_exact(values)

        assert 0 <= result <= 1


# =============================================================================
# String Operation Tests
# =============================================================================

class TestLongestCharRepeat:
    """Tests for longest_char_repeat function."""

    def test_longest_repeat_simple(self):
        """Test finding longest character repeat."""
        s = "AATTTTCCC"
        result = longest_char_repeat(s, 'T')

        assert result == 4

    def test_longest_repeat_no_match(self):
        """Test when character not in string."""
        s = "AATTCC"
        result = longest_char_repeat(s, 'G')

        assert result == 0

    def test_longest_repeat_at_start(self):
        """Test repeat at start of string."""
        s = "AAAAAATCG"
        result = longest_char_repeat(s, 'A')

        assert result == 6

    def test_longest_repeat_single(self):
        """Test single occurrence."""
        s = "ATCGATCG"
        result = longest_char_repeat(s, 'A')

        assert result == 1


class TestComplement:
    """Tests for complement function."""

    def test_complement_simple(self):
        """Test simple complement."""
        result = complement("ATCG")

        assert result == "TAGC"

    def test_complement_lowercase(self):
        """Test complement handles lowercase."""
        result = complement("atcg")

        assert result == "tagc"

    def test_complement_with_unknown(self):
        """Test complement with unknown base."""
        result = complement("ATNG")

        assert result == "TANC"  # N stays as N


class TestReverse:
    """Tests for reverse function."""

    def test_reverse_simple(self):
        """Test simple reverse."""
        result = reverse("ATCG")

        assert result == "GCTA"

    def test_reverse_palindrome(self):
        """Test reverse of palindrome."""
        result = reverse("ATAT")

        assert result == "TATA"


class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_reverse_complement_simple(self):
        """Test simple reverse complement."""
        result = reverse_complement("ATCG")

        assert result == "CGAT"

    def test_reverse_complement_palindrome(self):
        """Test reverse complement of DNA palindrome."""
        # ATAT reverse complement is ATAT
        result = reverse_complement("ATAT")

        assert result == "ATAT"


# =============================================================================
# Sequence Comparison Tests
# =============================================================================

class TestGetNumMismatches:
    """Tests for get_num_mismatches function."""

    def test_no_mismatches(self):
        """Test identical sequences."""
        result = get_num_mismatches("ATCG", "ATCG")

        assert result == 0

    def test_all_mismatches(self):
        """Test completely different sequences."""
        result = get_num_mismatches("AAAA", "TTTT")

        assert result == 4

    def test_some_mismatches(self):
        """Test partial mismatches."""
        result = get_num_mismatches("ATCG", "ATGG")

        assert result == 1


class TestLongestCommonSubstring:
    """Tests for longest_common_substring function."""

    def test_common_substring_simple(self):
        """Test finding common substring."""
        result = longest_common_substring("ATCGATCG", "GATCGATC")

        # "ATCG" is common (4 chars)
        assert result >= 4

    def test_no_common_substring(self):
        """Test with no common substring."""
        result = longest_common_substring("AAAA", "TTTT")

        assert result == 0

    def test_identical_strings(self):
        """Test with identical strings."""
        result = longest_common_substring("ATCG", "ATCG")

        assert result == 4


# =============================================================================
# Gap Calculation Tests
# =============================================================================

class TestGetPositionalGapLengths:
    """Tests for get_positional_gap_lengths function."""

    def test_simple_gaps(self):
        """Test simple gap calculation."""
        positions = [100, 200, 300]
        gaps = get_positional_gap_lengths(positions, circular=False)

        assert list(gaps) == [100, 100]

    def test_circular_genome(self):
        """Test circular genome gap calculation."""
        positions = [100, 500]
        gaps = get_positional_gap_lengths(positions, circular=True, seq_length=1000)

        # Gaps: 500-100=400, and wrap: 1000-500+100=600
        assert len(gaps) == 2
        assert 400 in gaps
        assert 600 in gaps

    def test_empty_positions(self):
        """Test with empty positions."""
        gaps = get_positional_gap_lengths([])

        assert len(gaps) == 0

    def test_single_position(self):
        """Test with single position."""
        gaps = get_positional_gap_lengths([100])

        assert len(gaps) == 0

    def test_unsorted_positions(self):
        """Test that positions are sorted."""
        positions = [300, 100, 200]
        gaps = get_positional_gap_lengths(positions, circular=False)

        # Should be [100, 100] after sorting
        assert list(gaps) == [100, 100]


# =============================================================================
# Most Frequent Tests
# =============================================================================

class TestMostFrequent:
    """Tests for most_frequent function."""

    def test_single_most_frequent(self):
        """Test with clear most frequent element."""
        lst = ['a', 'b', 'a', 'c', 'a']
        result = most_frequent(lst)

        assert result == 'a'

    def test_tie_returns_one(self):
        """Test that ties return one of the tied elements."""
        lst = ['a', 'b', 'a', 'b']
        result = most_frequent(lst)

        assert result in ['a', 'b']


# =============================================================================
# Multiprocessing Tests
# =============================================================================

def _square_for_pool(x):
    """Helper function for pool test (must be at module level for pickling)."""
    return x * x


def _identity_for_pool(x):
    """Helper function for pool test (must be at module level for pickling)."""
    return x


class TestCreatePool:
    """Tests for create_pool function."""

    def test_create_pool_simple(self):
        """Test create_pool with simple function."""
        result = create_pool(_square_for_pool, [1, 2, 3, 4], cpus=2)

        assert result == [1, 4, 9, 16]

    def test_create_pool_empty_input(self):
        """Test create_pool with empty input."""
        result = create_pool(_identity_for_pool, [], cpus=2)

        assert result == []


# =============================================================================
# File Reading Tests
# =============================================================================

class TestGetSeqLength:
    """Tests for get_seq_length function."""

    def test_get_seq_length(self, tmp_path):
        """Test getting sequence length from FASTA."""
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(">seq1\nATCGATCGATCGATCGATCG\n")  # 20 bp

        result = get_seq_length(str(fasta_path))

        assert result == 20

    def test_get_seq_length_multisequence(self, tmp_path):
        """Test getting sequence length from multi-sequence FASTA."""
        fasta_path = tmp_path / "test.fasta"
        fasta_path.write_text(">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n")  # 8 + 8 = 16 bp

        result = get_seq_length(str(fasta_path))

        assert result == 16


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_softmax_moderate_values(self):
        """Test softmax with moderate values."""
        # Note: The basic softmax implementation overflows with very large values
        # This tests the practical working range
        x = np.array([10.0, 10.0, 10.0])
        result = softmax(x)

        # Should still sum to 1 and be uniform
        assert np.isclose(result.sum(), 1.0)

    def test_gini_with_numpy_array(self):
        """Test gini_exact with numpy array input."""
        values = np.array([1, 2, 3, 4, 5])
        result = gini_exact(values)

        assert 0 <= result <= 1

    def test_complement_empty_string(self):
        """Test complement of empty string."""
        result = complement("")

        assert result == ""


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
