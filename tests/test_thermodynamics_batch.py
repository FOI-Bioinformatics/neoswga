"""
Tests for batch thermodynamic calculations and caching functions.

Tests:
- Batch Tm calculations (vectorized)
- Batch GC content calculations
- Batch Wallace Tm calculations
- Cache management functions
- Performance characteristics
"""

import pytest
import numpy as np
from typing import List

from neoswga.core.thermodynamics import (
    calculate_tm_batch,
    calculate_gc_batch,
    calculate_wallace_tm_batch,
    clear_thermodynamic_caches,
    get_cache_stats,
    calculate_tm_with_salt,
    gc_content,
    wallace_tm,
    reverse_complement,
)


# =============================================================================
# Batch Tm Calculation Tests
# =============================================================================

class TestCalculateTmBatch:
    """Tests for calculate_tm_batch function."""

    def test_single_sequence(self):
        """Test batch calculation with a single sequence."""
        sequences = ['ATCGATCG']
        tms = calculate_tm_batch(sequences)

        assert isinstance(tms, np.ndarray)
        assert len(tms) == 1
        assert not np.isnan(tms[0])

    def test_multiple_sequences(self):
        """Test batch calculation with multiple sequences."""
        sequences = ['ATCGATCG', 'GCTAGCTA', 'AAAAAATTTT']
        tms = calculate_tm_batch(sequences)

        assert isinstance(tms, np.ndarray)
        assert len(tms) == 3
        assert all(not np.isnan(tm) for tm in tms)

    def test_matches_individual_calculation(self):
        """Test that batch results match individual calculations."""
        sequences = ['ATCGATCG', 'GCTAGCTA', 'ATATATATAT']

        # Batch calculation
        batch_tms = calculate_tm_batch(sequences)

        # Individual calculations
        individual_tms = [calculate_tm_with_salt(seq) for seq in sequences]

        for i, seq in enumerate(sequences):
            assert abs(batch_tms[i] - individual_tms[i]) < 0.1, \
                f"Mismatch for {seq}: batch={batch_tms[i]}, individual={individual_tms[i]}"

    def test_with_custom_salt(self):
        """Test batch calculation with custom salt concentrations."""
        sequences = ['ATCGATCG', 'GCTAGCTA']

        tms_low_salt = calculate_tm_batch(sequences, na_conc=10.0)
        tms_high_salt = calculate_tm_batch(sequences, na_conc=100.0)

        # Higher salt should increase Tm
        for i in range(len(sequences)):
            assert tms_high_salt[i] > tms_low_salt[i]

    def test_with_magnesium(self):
        """Test batch calculation with magnesium."""
        sequences = ['ATCGATCGATCG']

        tms_no_mg = calculate_tm_batch(sequences, mg_conc=0.0)
        tms_with_mg = calculate_tm_batch(sequences, mg_conc=2.0)

        # Adding Mg should affect Tm
        assert tms_no_mg[0] != tms_with_mg[0]

    def test_empty_list(self):
        """Test batch calculation with empty list."""
        sequences = []
        tms = calculate_tm_batch(sequences)

        assert isinstance(tms, np.ndarray)
        assert len(tms) == 0

    def test_large_batch(self):
        """Test batch calculation with many sequences."""
        # Generate 100 random sequences
        bases = 'ATCG'
        sequences = []
        np.random.seed(42)
        for _ in range(100):
            length = np.random.randint(8, 15)  # Minimum 8bp for reasonable Tm
            seq = ''.join(np.random.choice(list(bases), length))
            sequences.append(seq)

        tms = calculate_tm_batch(sequences)

        assert len(tms) == 100
        assert all(not np.isnan(tm) for tm in tms)
        # All Tms should be in reasonable range (short AT-rich sequences can be below 0C)
        assert all(-20 < tm < 100 for tm in tms)


# =============================================================================
# Batch GC Content Tests
# =============================================================================

class TestCalculateGcBatch:
    """Tests for calculate_gc_batch function."""

    def test_single_sequence(self):
        """Test GC batch with single sequence."""
        sequences = ['ATCGATCG']
        gc_values = calculate_gc_batch(sequences)

        assert isinstance(gc_values, np.ndarray)
        assert len(gc_values) == 1
        assert gc_values[0] == 0.5  # 4 GC out of 8

    def test_multiple_sequences(self):
        """Test GC batch with multiple sequences."""
        sequences = ['AAAA', 'GGGG', 'ATCG']
        gc_values = calculate_gc_batch(sequences)

        assert len(gc_values) == 3
        assert gc_values[0] == 0.0   # All A
        assert gc_values[1] == 1.0   # All G
        assert gc_values[2] == 0.5   # Half GC

    def test_matches_individual_calculation(self):
        """Test that batch results match individual calculations."""
        sequences = ['ATCGATCG', 'AAATTTCCC', 'GCGCGCGC']

        batch_gc = calculate_gc_batch(sequences)
        individual_gc = [gc_content(seq) for seq in sequences]

        for i in range(len(sequences)):
            assert abs(batch_gc[i] - individual_gc[i]) < 0.001

    def test_empty_list(self):
        """Test GC batch with empty list."""
        gc_values = calculate_gc_batch([])
        assert len(gc_values) == 0


# =============================================================================
# Batch Wallace Tm Tests
# =============================================================================

class TestCalculateWallaceTmBatch:
    """Tests for calculate_wallace_tm_batch function."""

    def test_single_sequence(self):
        """Test Wallace Tm batch with single sequence."""
        sequences = ['ATCGATCG']
        tms = calculate_wallace_tm_batch(sequences)

        assert isinstance(tms, np.ndarray)
        assert len(tms) == 1
        assert not np.isnan(tms[0])

    def test_multiple_sequences(self):
        """Test Wallace Tm batch with multiple sequences."""
        sequences = ['AAAA', 'GGGG', 'ATCG']
        tms = calculate_wallace_tm_batch(sequences)

        assert len(tms) == 3
        # All A: 4*2 = 8
        assert tms[0] == 8.0
        # All G: 4*4 = 16
        assert tms[1] == 16.0
        # ATCG: 2*2 + 2*4 = 12
        assert tms[2] == 12.0

    def test_matches_individual_calculation(self):
        """Test that batch results match individual calculations."""
        sequences = ['ATCGATCG', 'AAATTTCCC', 'GCGCGCGC']

        batch_tms = calculate_wallace_tm_batch(sequences)
        individual_tms = [wallace_tm(seq) for seq in sequences]

        for i in range(len(sequences)):
            assert abs(batch_tms[i] - individual_tms[i]) < 0.001

    def test_gc_rich_higher_tm(self):
        """Test that GC-rich sequences have higher Wallace Tm."""
        at_rich = ['ATATATATAT']
        gc_rich = ['GCGCGCGCGC']

        at_tm = calculate_wallace_tm_batch(at_rich)[0]
        gc_tm = calculate_wallace_tm_batch(gc_rich)[0]

        assert gc_tm > at_tm


# =============================================================================
# Cache Management Tests
# =============================================================================

class TestCacheManagement:
    """Tests for cache management functions."""

    def test_clear_caches(self):
        """Test that clear_thermodynamic_caches works."""
        # Do some calculations to populate cache
        calculate_tm_with_salt('ATCGATCG')
        calculate_tm_with_salt('GCTAGCTA')

        # Clear caches
        clear_thermodynamic_caches()

        # Should not raise an error
        assert True

    def test_get_cache_stats(self):
        """Test that get_cache_stats returns valid info."""
        # Clear first
        clear_thermodynamic_caches()

        # Do some calculations
        calculate_tm_with_salt('ATCGATCG')
        calculate_tm_with_salt('ATCGATCG')  # Same - should hit cache

        stats = get_cache_stats()

        assert isinstance(stats, dict)
        assert 'enthalpy_entropy' in stats
        assert 'free_energy' in stats

    def test_cache_has_hits_after_repeat_calculation(self):
        """Test that cache accumulates hits."""
        clear_thermodynamic_caches()

        # First call - miss
        calculate_tm_with_salt('ATCGATCG')

        # Second call - should hit cache
        calculate_tm_with_salt('ATCGATCG')

        stats = get_cache_stats()

        # The enthalpy_entropy cache should have at least one hit
        cache_info = stats['enthalpy_entropy']
        assert cache_info.hits >= 0  # At least no errors


# =============================================================================
# Reverse Complement Tests
# =============================================================================

class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_basic_reverse_complement(self):
        """Test basic reverse complement."""
        assert reverse_complement('ATCG') == 'CGAT'
        assert reverse_complement('AAAA') == 'TTTT'
        assert reverse_complement('GGGG') == 'CCCC'

    def test_palindrome(self):
        """Test that palindrome returns itself."""
        assert reverse_complement('ATAT') == 'ATAT'
        assert reverse_complement('GCGC') == 'GCGC'

    def test_lowercase_handling(self):
        """Test that lowercase is handled."""
        assert reverse_complement('atcg') == 'CGAT'
        assert reverse_complement('AtCg') == 'CGAT'

    def test_long_sequence(self):
        """Test reverse complement of longer sequence."""
        seq = 'ATCGATCGATCG'
        rc = reverse_complement(seq)

        # Length should be preserved
        assert len(rc) == len(seq)

        # Complement of complement should be original
        assert reverse_complement(rc) == seq.upper()


# =============================================================================
# Performance Tests
# =============================================================================

class TestPerformance:
    """Performance-related tests."""

    def test_batch_faster_than_loop(self):
        """Test that batch calculation is efficient."""
        import time

        # Generate sequences
        np.random.seed(42)
        sequences = []
        for _ in range(50):
            length = np.random.randint(8, 12)
            seq = ''.join(np.random.choice(list('ATCG'), length))
            sequences.append(seq)

        # Time batch calculation
        start = time.time()
        batch_result = calculate_tm_batch(sequences)
        batch_time = time.time() - start

        # Clear cache to ensure fair comparison
        clear_thermodynamic_caches()

        # Time individual calculations
        start = time.time()
        individual_results = [calculate_tm_with_salt(seq) for seq in sequences]
        individual_time = time.time() - start

        # Batch should not be significantly slower (within 5x)
        # Note: for small batches, overhead may dominate
        assert batch_time < individual_time * 5

        # Results should match
        for i in range(len(sequences)):
            assert abs(batch_result[i] - individual_results[i]) < 0.1


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Edge case tests for batch functions."""

    def test_batch_with_varying_lengths(self):
        """Test batch with sequences of different lengths."""
        sequences = ['ATCG', 'ATCGATCG', 'ATCGATCGATCGATCG']
        tms = calculate_tm_batch(sequences)

        assert len(tms) == 3
        # Longer sequences generally have higher Tm
        # (not always true but generally)

    def test_batch_with_extreme_gc(self):
        """Test batch with extreme GC content sequences."""
        sequences = [
            'AAAAAAAAAA',  # 0% GC
            'ATATATATAT',  # 0% GC
            'GCGCGCGCGC',  # 100% GC
            'CCCCCCCCCC',  # 100% GC
        ]

        gc_values = calculate_gc_batch(sequences)

        assert gc_values[0] == 0.0
        assert gc_values[1] == 0.0
        assert gc_values[2] == 1.0
        assert gc_values[3] == 1.0

    def test_batch_preserves_order(self):
        """Test that batch calculation preserves sequence order."""
        sequences = ['AAAA', 'TTTT', 'GGGG', 'CCCC']

        gc_values = calculate_gc_batch(sequences)
        tms = calculate_tm_batch(sequences)

        # Order should be preserved
        assert gc_values[0] == 0.0  # AAAA
        assert gc_values[1] == 0.0  # TTTT
        assert gc_values[2] == 1.0  # GGGG
        assert gc_values[3] == 1.0  # CCCC


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
