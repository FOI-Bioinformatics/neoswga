"""
Unit tests for genome_io module.

Tests:
- GenomeStats dataclass
- GenomeLoader class (FASTA, gzip, zip)
- GenomeCache class
- Format detection
- Validation functions
"""

import pytest
import gzip
import zipfile
import tempfile
from pathlib import Path

from neoswga.core.genome_io import (
    GenomeStats,
    GenomeLoader,
    GenomeCache,
    load_genome,
    get_genome_stats,
    validate_genome_file,
)


# =============================================================================
# Test Fixtures
# =============================================================================

@pytest.fixture
def simple_fasta_file(tmp_path):
    """Create a simple FASTA file."""
    fasta_path = tmp_path / "test.fasta"
    content = """>seq1
ATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTA
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def gzipped_fasta_file(tmp_path):
    """Create a gzipped FASTA file."""
    fasta_path = tmp_path / "test.fasta.gz"
    content = """>seq1
ATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTA
"""
    with gzip.open(fasta_path, 'wt') as f:
        f.write(content)
    return fasta_path


@pytest.fixture
def zipped_fasta_file(tmp_path):
    """Create a zipped FASTA file."""
    zip_path = tmp_path / "test.fasta.zip"
    content = """>seq1
ATCGATCGATCGATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTA
"""
    with zipfile.ZipFile(zip_path, 'w') as zf:
        zf.writestr("genome.fasta", content)
    return zip_path


@pytest.fixture
def fasta_with_n_bases(tmp_path):
    """Create a FASTA file with N bases."""
    fasta_path = tmp_path / "with_n.fasta"
    content = """>seq1
ATCGATCGNNNNNNNATCGATCG
>seq2
GCTAGCTAGCTAGCTAGCTA
"""
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def large_fasta_file(tmp_path):
    """Create a larger FASTA file for validation testing."""
    fasta_path = tmp_path / "large.fasta"
    # Create a sequence of 200,000 bp
    seq = "ATCG" * 50000
    content = f">seq1\n{seq}\n"
    fasta_path.write_text(content)
    return fasta_path


# =============================================================================
# GenomeStats Tests
# =============================================================================

class TestGenomeStats:
    """Tests for GenomeStats dataclass."""

    def test_creation(self, tmp_path):
        """Test creating GenomeStats."""
        stats = GenomeStats(
            file_path=tmp_path / "test.fasta",
            file_format="fasta",
            total_sequences=2,
            total_length=40,
            gc_content=0.5,
            n_count=0,
            n_fraction=0.0,
            sequence_lengths=[20, 20],
            mean_length=20.0,
            max_length=20,
            min_length=20
        )

        assert stats.total_sequences == 2
        assert stats.total_length == 40
        assert stats.gc_content == 0.5

    def test_string_representation(self, tmp_path):
        """Test GenomeStats __str__ method."""
        stats = GenomeStats(
            file_path=tmp_path / "test.fasta",
            file_format="fasta",
            total_sequences=2,
            total_length=40,
            gc_content=0.5,
            n_count=5,
            n_fraction=0.125,
            sequence_lengths=[20, 20],
            mean_length=20.0,
            max_length=20,
            min_length=20
        )

        str_repr = str(stats)

        assert "Genome Statistics:" in str_repr
        assert "test.fasta" in str_repr
        assert "40" in str_repr  # total length


# =============================================================================
# GenomeLoader Format Detection Tests
# =============================================================================

class TestFormatDetection:
    """Tests for file format detection."""

    def test_detect_plain_fasta(self, simple_fasta_file):
        """Test detecting plain FASTA format."""
        loader = GenomeLoader()
        format_type = loader.detect_format(simple_fasta_file)

        assert format_type == "fasta"

    def test_detect_gzip_fasta(self, gzipped_fasta_file):
        """Test detecting gzipped FASTA format."""
        loader = GenomeLoader()
        format_type = loader.detect_format(gzipped_fasta_file)

        assert format_type == "fasta.gz"

    def test_detect_zip_fasta(self, zipped_fasta_file):
        """Test detecting zipped FASTA format."""
        loader = GenomeLoader()
        format_type = loader.detect_format(zipped_fasta_file)

        assert format_type == "fasta.zip"

    def test_detect_by_magic_bytes_gzip(self, tmp_path):
        """Test detecting gzip by magic bytes."""
        # Create a file with gzip magic bytes but wrong extension
        gzip_path = tmp_path / "test.bin"
        content = """>seq1
ATCG
"""
        with gzip.open(gzip_path, 'wt') as f:
            f.write(content)

        loader = GenomeLoader()
        format_type = loader.detect_format(gzip_path)

        assert format_type == "fasta.gz"

    def test_detect_nonexistent_file_raises(self, tmp_path):
        """Test that detecting format of nonexistent file raises."""
        loader = GenomeLoader()

        with pytest.raises(FileNotFoundError):
            loader.detect_format(tmp_path / "nonexistent.fasta")


# =============================================================================
# GenomeLoader Loading Tests
# =============================================================================

class TestGenomeLoading:
    """Tests for genome loading functionality."""

    def test_load_plain_fasta(self, simple_fasta_file):
        """Test loading plain FASTA file."""
        loader = GenomeLoader()
        sequence = loader.load_genome(simple_fasta_file)

        assert isinstance(sequence, str)
        assert len(sequence) == 40  # Two 20 bp sequences
        assert 'ATCG' in sequence

    def test_load_gzipped_fasta(self, gzipped_fasta_file):
        """Test loading gzipped FASTA file."""
        loader = GenomeLoader()
        sequence = loader.load_genome(gzipped_fasta_file)

        assert isinstance(sequence, str)
        assert len(sequence) == 40

    def test_load_zipped_fasta(self, zipped_fasta_file):
        """Test loading zipped FASTA file."""
        loader = GenomeLoader()
        sequence = loader.load_genome(zipped_fasta_file)

        assert isinstance(sequence, str)
        assert len(sequence) == 40

    def test_load_genome_uppercase(self, tmp_path):
        """Test that loaded sequence is uppercase."""
        fasta_path = tmp_path / "lowercase.fasta"
        fasta_path.write_text(">seq1\natcgatcg\n")

        loader = GenomeLoader()
        sequence = loader.load_genome(fasta_path)

        assert sequence == "ATCGATCG"

    def test_load_genome_concatenates_sequences(self, simple_fasta_file):
        """Test that multiple sequences are concatenated."""
        loader = GenomeLoader()
        sequence = loader.load_genome(simple_fasta_file)

        # Should concatenate both sequences
        assert len(sequence) == 40

    def test_load_genome_computes_stats(self, simple_fasta_file):
        """Test that loading computes statistics."""
        loader = GenomeLoader()
        sequence = loader.load_genome(simple_fasta_file, return_stats=True)
        stats = loader.get_stats()

        assert stats is not None
        assert stats.total_sequences == 2
        assert stats.total_length == 40

    def test_load_genome_without_stats(self, simple_fasta_file):
        """Test loading without computing stats."""
        loader = GenomeLoader()
        sequence = loader.load_genome(simple_fasta_file, return_stats=False)

        # Stats should still be None
        stats = loader.get_stats()
        assert stats is None

    def test_load_empty_fasta_raises(self, tmp_path):
        """Test that loading empty FASTA raises."""
        empty_fasta = tmp_path / "empty.fasta"
        empty_fasta.write_text("")

        loader = GenomeLoader()

        with pytest.raises(ValueError, match="No sequences found"):
            loader.load_genome(empty_fasta)


# =============================================================================
# Streaming Loading Tests
# =============================================================================

class TestStreamingLoading:
    """Tests for streaming genome loading."""

    def test_load_genome_streaming(self, simple_fasta_file):
        """Test streaming loading yields sequences."""
        loader = GenomeLoader()
        sequences = list(loader.load_genome_streaming(simple_fasta_file))

        assert len(sequences) == 2
        assert all(isinstance(s, str) for s in sequences)

    def test_load_genome_streaming_gzip(self, gzipped_fasta_file):
        """Test streaming loading of gzipped file."""
        loader = GenomeLoader()
        sequences = list(loader.load_genome_streaming(gzipped_fasta_file))

        assert len(sequences) == 2


# =============================================================================
# GC Content Tests
# =============================================================================

class TestGCContent:
    """Tests for GC content calculation."""

    def test_gc_content_50_percent(self, tmp_path):
        """Test GC content calculation for 50% GC."""
        fasta_path = tmp_path / "gc50.fasta"
        fasta_path.write_text(">seq1\nATCGATCGATCGATCGATCG\n")

        loader = GenomeLoader()
        loader.load_genome(fasta_path)
        stats = loader.get_stats()

        assert stats.gc_content == pytest.approx(0.5, abs=0.01)

    def test_gc_content_high(self, tmp_path):
        """Test GC content calculation for high GC."""
        fasta_path = tmp_path / "gc_high.fasta"
        # 80% GC content
        fasta_path.write_text(">seq1\nGGGGCCCCAT\n")

        loader = GenomeLoader()
        loader.load_genome(fasta_path)
        stats = loader.get_stats()

        assert stats.gc_content == pytest.approx(0.8, abs=0.01)

    def test_gc_content_low(self, tmp_path):
        """Test GC content calculation for low GC."""
        fasta_path = tmp_path / "gc_low.fasta"
        # 20% GC content (2 G + 2 C out of 20 bases)
        fasta_path.write_text(">seq1\nAAAAAAAAGCGCAAAAAAAA\n")

        loader = GenomeLoader()
        loader.load_genome(fasta_path)
        stats = loader.get_stats()

        assert stats.gc_content == pytest.approx(0.2, abs=0.01)


# =============================================================================
# N Base Content Tests
# =============================================================================

class TestNContent:
    """Tests for N base counting."""

    def test_n_count(self, fasta_with_n_bases):
        """Test counting N bases."""
        loader = GenomeLoader()
        loader.load_genome(fasta_with_n_bases)
        stats = loader.get_stats()

        assert stats.n_count == 7
        assert stats.n_fraction > 0

    def test_no_n_bases(self, simple_fasta_file):
        """Test file with no N bases."""
        loader = GenomeLoader()
        loader.load_genome(simple_fasta_file)
        stats = loader.get_stats()

        assert stats.n_count == 0
        assert stats.n_fraction == 0.0


# =============================================================================
# Validation Tests
# =============================================================================

class TestValidation:
    """Tests for genome validation."""

    def test_validate_valid_genome(self, large_fasta_file):
        """Test validating a valid genome."""
        loader = GenomeLoader()
        is_valid, issues = loader.validate_genome(large_fasta_file)

        assert is_valid is True
        assert len(issues) == 0

    def test_validate_genome_too_short(self, simple_fasta_file):
        """Test validating a genome that's too short."""
        loader = GenomeLoader()
        is_valid, issues = loader.validate_genome(
            simple_fasta_file,
            min_length=1000000  # 1 Mbp minimum
        )

        assert is_valid is False
        assert any("too short" in issue.lower() for issue in issues)

    def test_validate_too_many_n_bases(self, tmp_path):
        """Test validating genome with too many N bases."""
        fasta_path = tmp_path / "many_n.fasta"
        # 50% N bases
        fasta_path.write_text(">seq1\nATCGNNNNNNNNNNATCG\n")

        loader = GenomeLoader()
        is_valid, issues = loader.validate_genome(
            fasta_path,
            max_n_fraction=0.10,
            min_length=10
        )

        assert is_valid is False
        assert any("N bases" in issue for issue in issues)

    def test_validate_unusual_gc(self, tmp_path):
        """Test validating genome with unusual GC content."""
        fasta_path = tmp_path / "weird_gc.fasta"
        # 90% GC content
        fasta_path.write_text(">seq1\n" + "G" * 900000 + "A" * 100000 + "\n")

        loader = GenomeLoader()
        is_valid, issues = loader.validate_genome(fasta_path)

        # Should have a warning about unusual GC
        assert any("GC content" in issue for issue in issues)


# =============================================================================
# GenomeCache Tests
# =============================================================================

class TestGenomeCache:
    """Tests for GenomeCache class."""

    def test_cache_stores_genome(self, simple_fasta_file):
        """Test that cache stores loaded genomes."""
        cache = GenomeCache()

        seq1, stats1 = cache.get(simple_fasta_file)
        seq2, stats2 = cache.get(simple_fasta_file)

        # Should return same data
        assert seq1 == seq2
        assert stats1.total_length == stats2.total_length

    def test_cache_eviction(self, tmp_path):
        """Test cache eviction when full."""
        # Create cache with max size 2
        cache = GenomeCache(max_cache_size=2)

        # Create 3 FASTA files
        for i in range(3):
            fasta_path = tmp_path / f"genome{i}.fasta"
            fasta_path.write_text(f">seq{i}\nATCGATCG\n")
            cache.get(fasta_path)

        # Cache should only have 2 entries
        assert len(cache.cache) == 2

    def test_cache_clear(self, simple_fasta_file):
        """Test clearing the cache."""
        cache = GenomeCache()

        cache.get(simple_fasta_file)
        assert len(cache.cache) == 1

        cache.clear()
        assert len(cache.cache) == 0

    def test_cache_lru_order(self, tmp_path):
        """Test that cache uses LRU eviction."""
        cache = GenomeCache(max_cache_size=2)

        # Create 3 files
        files = []
        for i in range(3):
            fasta_path = tmp_path / f"genome{i}.fasta"
            fasta_path.write_text(f">seq{i}\nATCGATCG\n")
            files.append(fasta_path)

        # Load first two
        cache.get(files[0])
        cache.get(files[1])

        # Access first again (makes it most recently used)
        cache.get(files[0])

        # Load third - should evict second (least recently used)
        cache.get(files[2])

        # First and third should be in cache
        keys = list(cache.cache.keys())
        assert str(files[0].resolve()) in keys
        assert str(files[2].resolve()) in keys


# =============================================================================
# Convenience Function Tests
# =============================================================================

class TestConvenienceFunctions:
    """Tests for module-level convenience functions."""

    def test_load_genome_function(self, simple_fasta_file):
        """Test load_genome convenience function."""
        sequence = load_genome(simple_fasta_file)

        assert isinstance(sequence, str)
        assert len(sequence) == 40

    def test_get_genome_stats_function(self, simple_fasta_file):
        """Test get_genome_stats convenience function."""
        stats = get_genome_stats(simple_fasta_file)

        assert isinstance(stats, GenomeStats)
        assert stats.total_sequences == 2

    def test_validate_genome_file_valid(self, large_fasta_file):
        """Test validate_genome_file with valid file."""
        result = validate_genome_file(large_fasta_file)

        assert result is True

    def test_validate_genome_file_invalid(self, simple_fasta_file):
        """Test validate_genome_file with invalid file (too short)."""
        with pytest.raises(ValueError, match="Invalid genome"):
            validate_genome_file(simple_fasta_file)


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases."""

    def test_single_base_sequence(self, tmp_path):
        """Test loading a single-base sequence."""
        fasta_path = tmp_path / "single.fasta"
        fasta_path.write_text(">seq1\nA\n")

        loader = GenomeLoader()
        sequence = loader.load_genome(fasta_path)

        assert sequence == "A"

    def test_multiline_sequence(self, tmp_path):
        """Test loading a multiline sequence."""
        fasta_path = tmp_path / "multiline.fasta"
        fasta_path.write_text(">seq1\nATCG\nGCTA\nATCG\n")

        loader = GenomeLoader()
        sequence = loader.load_genome(fasta_path)

        assert sequence == "ATCGGCTAATCG"

    def test_fa_extension(self, tmp_path):
        """Test loading file with .fa extension."""
        fasta_path = tmp_path / "test.fa"
        fasta_path.write_text(">seq1\nATCGATCG\n")

        loader = GenomeLoader()
        format_type = loader.detect_format(fasta_path)

        assert format_type == "fasta"

    def test_fna_extension(self, tmp_path):
        """Test loading file with .fna extension."""
        fasta_path = tmp_path / "test.fna"
        fasta_path.write_text(">seq1\nATCGATCG\n")

        loader = GenomeLoader()
        format_type = loader.detect_format(fasta_path)

        assert format_type == "fasta"

    def test_zip_without_fasta(self, tmp_path):
        """Test loading zip without FASTA file raises."""
        zip_path = tmp_path / "no_fasta.zip"
        with zipfile.ZipFile(zip_path, 'w') as zf:
            zf.writestr("readme.txt", "This is not a FASTA file")

        loader = GenomeLoader()

        with pytest.raises(ValueError, match="No FASTA file found"):
            loader.load_genome(zip_path)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
