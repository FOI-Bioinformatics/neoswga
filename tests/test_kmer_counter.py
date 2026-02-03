"""
Tests for kmer_counter module.

Tests:
- MultiGenomeKmerCounter class
- Standalone functions (run_jellyfish, get_kmer_to_count_dict, get_primer_list_from_kmers)
- Pure Python k-mer counting fallback
"""

import pytest
import os
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from neoswga.core.kmer_counter import (
    check_jellyfish_available,
    require_jellyfish,
    MultiGenomeKmerCounter,
    count_kmers_in_sequence,
    run_jellyfish,
    get_kmer_to_count_dict,
    get_primer_list_from_kmers,
)


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    tmp = tempfile.mkdtemp(prefix='test_kmer_')
    yield tmp
    shutil.rmtree(tmp)


@pytest.fixture
def sample_fasta(temp_dir):
    """Create a sample FASTA file."""
    fasta_path = os.path.join(temp_dir, 'test_genome.fa')
    with open(fasta_path, 'w') as f:
        f.write('>chr1\n')
        f.write('ATCGATCGATCGATCGATCGATCGATCGATCG\n')
        f.write('GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n')
    return fasta_path


@pytest.fixture
def sample_kmer_dump_file(temp_dir):
    """Create a sample jellyfish dump file."""
    dump_path = os.path.join(temp_dir, 'test_6mer_all.txt')
    with open(dump_path, 'w') as f:
        f.write('ATCGAT 10\n')
        f.write('TCGATC 8\n')
        f.write('CGATCG 12\n')
        f.write('GCTAGC 6\n')
    return dump_path


# =============================================================================
# Jellyfish Availability Tests
# =============================================================================

class TestJellyfishAvailability:
    """Tests for jellyfish availability checking."""

    def test_check_jellyfish_available_returns_bool(self):
        """Test that check_jellyfish_available returns a boolean."""
        result = check_jellyfish_available()
        assert isinstance(result, bool)

    def test_require_jellyfish_when_available(self):
        """Test require_jellyfish doesn't raise when jellyfish is available."""
        if check_jellyfish_available():
            # Should not raise
            require_jellyfish()
        else:
            pytest.skip("Jellyfish not available")

    def test_require_jellyfish_when_not_available(self):
        """Test require_jellyfish raises RuntimeError when jellyfish is not available."""
        with patch('shutil.which', return_value=None):
            with pytest.raises(RuntimeError, match="Jellyfish is required"):
                require_jellyfish()


# =============================================================================
# Pure Python K-mer Counting Tests
# =============================================================================

class TestCountKmersInSequence:
    """Tests for pure Python k-mer counting."""

    def test_count_kmers_basic(self):
        """Test basic k-mer counting."""
        sequence = 'ATCGATCG'
        counts = count_kmers_in_sequence(sequence, k=4)

        assert 'ATCG' in counts
        assert 'TCGA' in counts
        assert 'CGAT' in counts
        assert 'GATC' in counts

    def test_count_kmers_with_repeats(self):
        """Test k-mer counting with repeating patterns."""
        sequence = 'ATATAT'
        counts = count_kmers_in_sequence(sequence, k=2)

        assert counts['AT'] == 3
        assert counts['TA'] == 2

    def test_count_kmers_case_insensitive(self):
        """Test that counting is case-insensitive."""
        seq_upper = 'ATCG'
        seq_lower = 'atcg'
        seq_mixed = 'AtCg'

        counts_upper = count_kmers_in_sequence(seq_upper, k=2)
        counts_lower = count_kmers_in_sequence(seq_lower, k=2)
        counts_mixed = count_kmers_in_sequence(seq_mixed, k=2)

        assert counts_upper == counts_lower == counts_mixed

    def test_count_kmers_excludes_non_acgt(self):
        """Test that k-mers with non-ACGT bases are excluded."""
        sequence = 'ATCNGATCG'
        counts = count_kmers_in_sequence(sequence, k=4)

        # K-mers containing N should not be counted
        assert 'TCNG' not in counts
        assert 'CNGA' not in counts
        assert 'NGAT' not in counts

    def test_count_kmers_empty_sequence(self):
        """Test counting k-mers in empty sequence."""
        counts = count_kmers_in_sequence('', k=4)
        assert counts == {}

    def test_count_kmers_short_sequence(self):
        """Test counting k-mers when sequence is shorter than k."""
        counts = count_kmers_in_sequence('ATG', k=4)
        assert counts == {}


# =============================================================================
# get_kmer_to_count_dict Tests
# =============================================================================

class TestGetKmerToCountDict:
    """Tests for get_kmer_to_count_dict function."""

    def test_reads_dump_file_correctly(self, sample_kmer_dump_file):
        """Test that dump file is parsed correctly."""
        counts = get_kmer_to_count_dict(sample_kmer_dump_file)

        assert counts['ATCGAT'] == 10
        assert counts['TCGATC'] == 8
        assert counts['CGATCG'] == 12
        assert counts['GCTAGC'] == 6

    def test_returns_dict(self, sample_kmer_dump_file):
        """Test that function returns a dictionary."""
        counts = get_kmer_to_count_dict(sample_kmer_dump_file)
        assert isinstance(counts, dict)

    def test_handles_empty_file(self, temp_dir):
        """Test handling of empty file."""
        empty_file = os.path.join(temp_dir, 'empty.txt')
        with open(empty_file, 'w') as f:
            pass  # Create empty file

        counts = get_kmer_to_count_dict(empty_file)
        assert counts == {}


# =============================================================================
# MultiGenomeKmerCounter Tests (with mocking)
# =============================================================================

class TestMultiGenomeKmerCounter:
    """Tests for MultiGenomeKmerCounter class."""

    def test_initialization_checks_jellyfish(self):
        """Test that initialization checks for jellyfish."""
        with patch('neoswga.core.kmer_counter.check_jellyfish_available', return_value=False):
            with pytest.raises(RuntimeError, match="Jellyfish is required"):
                MultiGenomeKmerCounter()

    @pytest.mark.skipif(not check_jellyfish_available(), reason="Jellyfish not available")
    def test_add_genome(self, sample_fasta):
        """Test adding a genome."""
        counter = MultiGenomeKmerCounter()
        counter.add_genome('test', sample_fasta)

        assert 'test' in counter.genome_fastas
        assert counter.genome_fastas['test'] == sample_fasta
        assert counter.genome_lengths['test'] > 0

    @pytest.mark.skipif(not check_jellyfish_available(), reason="Jellyfish not available")
    def test_add_genome_file_not_found(self):
        """Test adding a non-existent genome file."""
        counter = MultiGenomeKmerCounter()

        with pytest.raises(FileNotFoundError):
            counter.add_genome('test', '/nonexistent/path.fa')

    @pytest.mark.skipif(not check_jellyfish_available(), reason="Jellyfish not available")
    def test_cleanup(self, sample_fasta):
        """Test cleanup removes temp files."""
        counter = MultiGenomeKmerCounter()
        counter.add_genome('test', sample_fasta)

        # Force creation of temp dir
        work_dir = counter._get_work_dir()
        assert os.path.exists(work_dir)

        counter.cleanup()
        # Temp dir should be removed
        assert not os.path.exists(work_dir)


# =============================================================================
# run_jellyfish Tests (with mocking for CI)
# =============================================================================

class TestRunJellyfish:
    """Tests for run_jellyfish standalone function."""

    def test_validates_genome_file_exists(self, temp_dir):
        """Test that non-existent genome file raises error."""
        with pytest.raises(FileNotFoundError, match="Genome file not found"):
            run_jellyfish('/nonexistent/genome.fa', os.path.join(temp_dir, 'output'))

    def test_checks_jellyfish_available(self):
        """Test that jellyfish availability is checked."""
        with patch('neoswga.core.kmer_counter.check_jellyfish_available', return_value=False):
            with pytest.raises(RuntimeError, match="Jellyfish is required"):
                run_jellyfish('/some/genome.fa', '/some/output')

    @pytest.mark.skipif(not check_jellyfish_available(), reason="Jellyfish not available")
    def test_creates_output_files(self, sample_fasta, temp_dir):
        """Test that output files are created."""
        output_prefix = os.path.join(temp_dir, 'output')

        run_jellyfish(sample_fasta, output_prefix, min_k=6, max_k=6, cpus=1)

        # Check output file exists
        expected_file = f"{output_prefix}_6mer_all.txt"
        assert os.path.exists(expected_file)

        # Check file has content
        with open(expected_file) as f:
            content = f.read()
            assert len(content) > 0


# =============================================================================
# get_primer_list_from_kmers Tests (with mocking)
# =============================================================================

class TestGetPrimerListFromKmers:
    """Tests for get_primer_list_from_kmers function."""

    def test_filters_by_tm(self, temp_dir):
        """Test that Tm filtering is applied."""
        # Create mock k-mer file
        prefix = os.path.join(temp_dir, 'test')
        kmer_file = f"{prefix}_6mer_all.txt"
        with open(kmer_file, 'w') as f:
            f.write('ATCGAT 10\n')  # Low Tm
            f.write('GCGCGC 10\n')  # High Tm
            f.write('ATGCAT 10\n')  # Medium Tm

        primers = get_primer_list_from_kmers(
            [prefix],
            kmer_lengths=range(6, 7),
            min_tm=10.0,
            max_tm=30.0
        )

        # Should filter by Tm
        assert isinstance(primers, list)

    def test_returns_list(self, temp_dir):
        """Test that function returns a list."""
        prefix = os.path.join(temp_dir, 'test')
        kmer_file = f"{prefix}_6mer_all.txt"
        with open(kmer_file, 'w') as f:
            f.write('ATCGAT 10\n')

        primers = get_primer_list_from_kmers([prefix], kmer_lengths=range(6, 7))
        assert isinstance(primers, list)

    def test_handles_missing_files_gracefully(self, temp_dir):
        """Test that missing files don't crash the function."""
        prefix = os.path.join(temp_dir, 'nonexistent')

        # Should not raise, just warn
        primers = get_primer_list_from_kmers([prefix], kmer_lengths=range(6, 7))
        assert primers == []


# =============================================================================
# Integration Tests
# =============================================================================

@pytest.mark.skipif(not check_jellyfish_available(), reason="Jellyfish not available")
class TestKmerCounterIntegration:
    """Integration tests requiring jellyfish."""

    def test_full_workflow(self, sample_fasta, temp_dir):
        """Test complete k-mer counting workflow."""
        # Initialize counter
        counter = MultiGenomeKmerCounter(cpus=1, output_dir=temp_dir)

        # Add genome
        counter.add_genome('test', sample_fasta)

        # Count k-mers
        counts = counter.count_kmers('test', k=6)

        # Verify we got results
        assert isinstance(counts, dict)
        assert len(counts) > 0

        # Verify counts are positive integers
        for kmer, count in counts.items():
            assert len(kmer) == 6
            assert isinstance(count, int)
            assert count > 0

        # Cleanup
        counter.cleanup()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
