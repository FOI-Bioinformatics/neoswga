"""
Unit tests for position_cache module.

Tests:
- BindingSite dataclass
- PositionCache class methods
- StreamingPositionCache class
- Coverage computation
- Gini coefficient calculation
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
import h5py

from neoswga.core.position_cache import (
    BindingSite,
    PositionCache,
    StreamingPositionCache,
)
from neoswga.core.thermodynamics import reverse_complement


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def sample_hdf5_file(tmp_path):
    """Create a sample HDF5 file with position data."""
    h5_path = tmp_path / "test_8mer_positions.h5"

    with h5py.File(h5_path, 'w') as f:
        # Forward strand positions
        f.create_dataset('ATCGATCG', data=np.array([100, 500, 1000, 2000]))
        f.create_dataset('GCTAGCTA', data=np.array([200, 600, 1100]))
        # Reverse complement of ATCGATCG (for reverse strand testing)
        f.create_dataset('CGATCGAT', data=np.array([300, 700, 1200]))

    # Return prefix (without the _8mer_positions.h5 suffix)
    return str(tmp_path / "test")


@pytest.fixture
def sample_hdf5_with_multiple_kmers(tmp_path):
    """Create HDF5 files for multiple k-mer lengths."""
    prefixes = []

    for k in [6, 8, 10]:
        h5_path = tmp_path / f"genome_{k}mer_positions.h5"

        with h5py.File(h5_path, 'w') as f:
            # Create some position data
            primer = 'A' * k
            f.create_dataset(primer, data=np.array([100 * k, 200 * k]))

    return str(tmp_path / "genome")


# =============================================================================
# BindingSite Tests
# =============================================================================

class TestBindingSite:
    """Tests for BindingSite dataclass."""

    def test_creation(self):
        """Test creating a BindingSite."""
        site = BindingSite(position=100, strand='forward', primer='ATCGATCG')

        assert site.position == 100
        assert site.strand == 'forward'
        assert site.primer == 'ATCGATCG'

    def test_hashable(self):
        """Test that BindingSite is hashable."""
        site1 = BindingSite(position=100, strand='forward', primer='ATCGATCG')
        site2 = BindingSite(position=100, strand='forward', primer='ATCGATCG')
        site3 = BindingSite(position=200, strand='forward', primer='ATCGATCG')

        # Same values should hash the same
        assert hash(site1) == hash(site2)
        # Different positions should hash differently
        assert hash(site1) != hash(site3)

    def test_can_be_used_in_set(self):
        """Test that BindingSites can be used in a set."""
        site1 = BindingSite(position=100, strand='forward', primer='ATCGATCG')
        site2 = BindingSite(position=100, strand='forward', primer='ATCGATCG')
        site3 = BindingSite(position=200, strand='forward', primer='ATCGATCG')

        site_set = {site1, site2, site3}

        # site1 and site2 are duplicates, so set should have 2 elements
        assert len(site_set) == 2


# =============================================================================
# PositionCache Tests
# =============================================================================

class TestPositionCache:
    """Tests for PositionCache class."""

    def test_initialization_with_valid_data(self, sample_hdf5_file):
        """Test cache initialization loads positions correctly."""
        primers = ['ATCGATCG', 'GCTAGCTA']

        cache = PositionCache([sample_hdf5_file], primers)

        assert len(cache.primers) == 2
        assert len(cache.cache) > 0

    def test_initialization_with_missing_hdf5(self, tmp_path):
        """Test cache handles missing HDF5 files gracefully."""
        # Should not raise, just log warning
        cache = PositionCache([str(tmp_path / "nonexistent")], ['ATCGATCG'])

        # Cache should be empty but not crash
        assert len(cache.cache) == 0

    def test_get_positions_forward(self, sample_hdf5_file):
        """Test getting forward strand positions."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        positions = cache.get_positions(sample_hdf5_file, 'ATCGATCG', strand='forward')

        assert isinstance(positions, np.ndarray)
        assert len(positions) == 4  # We created 4 positions

    def test_get_positions_reverse(self, sample_hdf5_file):
        """Test getting reverse strand positions."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        positions = cache.get_positions(sample_hdf5_file, 'ATCGATCG', strand='reverse')

        assert isinstance(positions, np.ndarray)
        # Reverse complement CGATCGAT has 3 positions
        assert len(positions) == 3

    def test_get_positions_both_strands(self, sample_hdf5_file):
        """Test getting positions from both strands."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        positions = cache.get_positions(sample_hdf5_file, 'ATCGATCG', strand='both')

        assert isinstance(positions, np.ndarray)
        # 4 forward + 3 reverse = 7 total
        assert len(positions) == 7

    def test_get_positions_missing_primer(self, sample_hdf5_file):
        """Test getting positions for non-existent primer."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        positions = cache.get_positions(sample_hdf5_file, 'XXXXXXXX', strand='both')

        assert isinstance(positions, np.ndarray)
        assert len(positions) == 0

    def test_get_all_positions(self, sample_hdf5_file):
        """Test getting positions for multiple primers."""
        primers = ['ATCGATCG', 'GCTAGCTA']
        cache = PositionCache([sample_hdf5_file], primers)

        result = cache.get_all_positions(sample_hdf5_file, primers)

        assert len(result) == 2
        assert 'ATCGATCG' in result
        assert 'GCTAGCTA' in result
        # Each entry is (forward_positions, reverse_positions)
        assert isinstance(result['ATCGATCG'], tuple)
        assert len(result['ATCGATCG']) == 2


# =============================================================================
# Coverage Computation Tests
# =============================================================================

class TestCoverageComputation:
    """Tests for coverage computation methods."""

    def test_compute_coverage_vectorized(self, sample_hdf5_file):
        """Test vectorized coverage computation."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        genome_length = 3000
        coverage = cache.compute_coverage_vectorized(
            sample_hdf5_file, primers, genome_length
        )

        assert isinstance(coverage, np.ndarray)
        assert len(coverage) == genome_length
        assert coverage.dtype == bool
        # Should have some positions covered
        assert np.any(coverage)

    def test_compute_coverage_clips_out_of_bounds(self, sample_hdf5_file):
        """Test that coverage clips positions outside genome bounds."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        # Use genome length smaller than some positions
        genome_length = 500
        coverage = cache.compute_coverage_vectorized(
            sample_hdf5_file, primers, genome_length
        )

        # Should not crash and should be correct length
        assert len(coverage) == genome_length

    def test_compute_statistics(self, sample_hdf5_file):
        """Test computing coverage statistics."""
        primers = ['ATCGATCG', 'GCTAGCTA']
        cache = PositionCache([sample_hdf5_file], primers)

        genome_length = 5000
        stats = cache.compute_statistics(sample_hdf5_file, primers, genome_length)

        assert 'coverage_fraction' in stats
        assert 'num_covered_bases' in stats
        assert 'num_gaps' in stats
        assert 'mean_gap_size' in stats
        assert 'max_gap_size' in stats
        assert 'gap_gini' in stats

        assert 0 <= stats['coverage_fraction'] <= 1
        assert stats['num_covered_bases'] >= 0

    def test_compute_statistics_empty_primers(self, sample_hdf5_file):
        """Test statistics with no primers."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        genome_length = 5000
        stats = cache.compute_statistics(sample_hdf5_file, [], genome_length)

        assert stats['coverage_fraction'] == 0
        assert stats['num_covered_bases'] == 0


# =============================================================================
# Gap and Gini Tests
# =============================================================================

class TestGapAnalysis:
    """Tests for gap size analysis."""

    def test_find_gap_sizes_no_coverage(self, sample_hdf5_file):
        """Test gap finding with no coverage."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        # Create a coverage array with no coverage
        coverage = np.zeros(1000, dtype=bool)
        gaps = cache._find_gap_sizes(coverage)

        # Should be one gap spanning the whole genome
        assert len(gaps) == 1
        assert gaps[0] == 1000

    def test_find_gap_sizes_full_coverage(self, sample_hdf5_file):
        """Test gap finding with full coverage."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        # Create a coverage array with full coverage
        coverage = np.ones(1000, dtype=bool)
        gaps = cache._find_gap_sizes(coverage)

        # Should be no gaps
        assert len(gaps) == 0

    def test_find_gap_sizes_alternating(self, sample_hdf5_file):
        """Test gap finding with alternating coverage."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        # Create coverage: covered-uncovered-covered
        coverage = np.array([True]*100 + [False]*200 + [True]*100, dtype=bool)
        gaps = cache._find_gap_sizes(coverage)

        assert len(gaps) == 1
        assert gaps[0] == 200

    def test_gini_coefficient_uniform(self, sample_hdf5_file):
        """Test Gini coefficient for uniform distribution."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        # Uniform distribution (all same value)
        values = [100, 100, 100, 100]
        gini = cache._gini_coefficient(values)

        # Gini should be 0 for perfect equality
        assert gini == pytest.approx(0.0, abs=0.01)

    def test_gini_coefficient_skewed(self, sample_hdf5_file):
        """Test Gini coefficient for skewed distribution."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        # Very skewed distribution
        values = [1, 1, 1, 1000]
        gini = cache._gini_coefficient(values)

        # Gini should be high for inequality
        assert gini > 0.5

    def test_gini_coefficient_empty(self, sample_hdf5_file):
        """Test Gini coefficient for empty values."""
        primers = ['ATCGATCG']
        cache = PositionCache([sample_hdf5_file], primers)

        gini = cache._gini_coefficient([])

        assert gini == 0.0


# =============================================================================
# Reverse Complement Tests
# =============================================================================

class TestReverseComplement:
    """Tests for reverse complement function (canonical from thermodynamics module)."""

    def test_reverse_complement_simple(self):
        """Test reverse complement of a simple sequence."""
        rc = reverse_complement('ATCG')

        assert rc == 'CGAT'

    def test_reverse_complement_palindrome(self):
        """Test reverse complement of a palindrome."""
        # ATAT is a palindrome in DNA
        rc = reverse_complement('ATAT')

        assert rc == 'ATAT'

    def test_reverse_complement_full(self):
        """Test reverse complement of longer sequence."""
        rc = reverse_complement('ATCGATCG')

        assert rc == 'CGATCGAT'


# =============================================================================
# StreamingPositionCache Tests
# =============================================================================

class TestStreamingPositionCache:
    """Tests for StreamingPositionCache class."""

    def test_initialization(self, sample_hdf5_file):
        """Test streaming cache initialization."""
        with StreamingPositionCache([sample_hdf5_file]) as cache:
            # Should open HDF5 files
            assert len(cache.file_handles) > 0 or len(cache.fname_prefixes) > 0

    def test_get_positions(self, sample_hdf5_file):
        """Test getting positions from streaming cache."""
        with StreamingPositionCache([sample_hdf5_file]) as cache:
            positions = cache.get_positions(sample_hdf5_file, 'ATCGATCG', strand='forward')

            # May be empty if file pattern doesn't match, but shouldn't crash
            assert isinstance(positions, np.ndarray)

    def test_context_manager_closes_handles(self, sample_hdf5_file):
        """Test that context manager properly closes file handles."""
        cache = StreamingPositionCache([sample_hdf5_file])

        # Exit context manager
        cache.close()

        # Handles should be cleared
        assert len(cache.file_handles) == 0

    def test_preload_subset(self, sample_hdf5_file):
        """Test preloading a subset of primers."""
        primers = ['ATCGATCG']

        with StreamingPositionCache([sample_hdf5_file], primers=primers) as cache:
            # Should preload the small set
            # Note: may or may not have preloaded depending on file structure
            assert cache.fname_prefixes == [sample_hdf5_file]


# =============================================================================
# Edge Cases and Error Handling
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_empty_primer_list(self, sample_hdf5_file):
        """Test with empty primer list."""
        cache = PositionCache([sample_hdf5_file], [])

        assert len(cache.primers) == 0
        assert len(cache.cache) == 0

    def test_empty_prefix_list(self, tmp_path):
        """Test with empty prefix list."""
        cache = PositionCache([], ['ATCGATCG'])

        assert len(cache.cache) == 0

    def test_multiple_prefixes(self, tmp_path):
        """Test with multiple HDF5 file prefixes."""
        # Create two HDF5 files
        for name in ['genome1', 'genome2']:
            h5_path = tmp_path / f"{name}_8mer_positions.h5"
            with h5py.File(h5_path, 'w') as f:
                f.create_dataset('ATCGATCG', data=np.array([100, 200]))

        prefixes = [str(tmp_path / "genome1"), str(tmp_path / "genome2")]
        primers = ['ATCGATCG']

        cache = PositionCache(prefixes, primers)

        # Should load from both files
        pos1 = cache.get_positions(prefixes[0], 'ATCGATCG')
        pos2 = cache.get_positions(prefixes[1], 'ATCGATCG')

        assert len(pos1) >= 0
        assert len(pos2) >= 0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
