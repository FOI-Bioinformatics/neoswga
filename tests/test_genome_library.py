"""Tests for genome library module."""

import json
import os
import pytest
from unittest.mock import patch, MagicMock


class TestGenomeLibraryEntry:
    """Test GenomeLibraryEntry dataclass."""

    def test_serialization_roundtrip(self):
        from neoswga.core.genome_library import GenomeLibraryEntry

        entry = GenomeLibraryEntry(
            name="test-genome",
            species="Test species",
            fasta_path="/tmp/test.fna",
            genome_size=1000000,
            gc_content=0.45,
            role="bg",
            kmer_dir="/tmp/test-genome",
            kmer_prefix="/tmp/test-genome/test-genome",
            computed_k_ranges=[(6, 12)],
        )
        d = entry.to_dict()
        restored = GenomeLibraryEntry.from_dict(d)
        assert restored.name == entry.name
        assert restored.genome_size == entry.genome_size
        assert restored.gc_content == entry.gc_content
        assert restored.computed_k_ranges == [(6, 12)]

    def test_default_role(self):
        from neoswga.core.genome_library import GenomeLibraryEntry

        entry = GenomeLibraryEntry(
            name="test", species="", fasta_path="", genome_size=0, gc_content=0.0
        )
        assert entry.role == "bg"


class TestGenomeLibrary:
    """Test GenomeLibrary class."""

    def test_init_creates_directory(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        lib_dir = str(tmp_path / "genomes")
        library = GenomeLibrary(library_dir=lib_dir)
        assert os.path.isdir(lib_dir)

    def test_list_empty(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        assert library.list() == []

    def test_get_nonexistent(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        assert library.get("nonexistent") is None

    def test_get_kmer_prefix_nonexistent(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        assert library.get_kmer_prefix("nonexistent") is None

    def test_has_kmers_for_range_nonexistent(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        assert library.has_kmers_for_range("nonexistent", 6, 12) is False

    def test_remove_nonexistent(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        assert library.remove("nonexistent") is False

    def test_add_missing_fasta(self, tmp_path):
        from neoswga.core.genome_library import GenomeLibrary

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        with pytest.raises(FileNotFoundError):
            library.add("test", "/nonexistent/path.fna")

    @patch("neoswga.core.genome_library.GenomeLibrary._calculate_genome_stats")
    @patch("neoswga.core.kmer_counter.run_jellyfish")
    def test_add_and_get_lifecycle(self, mock_jf, mock_stats, tmp_path):
        """Test add -> get -> list -> remove lifecycle."""
        from neoswga.core.genome_library import GenomeLibrary

        mock_stats.return_value = (100000, 0.45)
        mock_jf.return_value = None

        # Create a dummy FASTA file
        fasta = tmp_path / "test.fna"
        fasta.write_text(">seq1\nATCG\n")

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        entry = library.add("test-genome", str(fasta), k_ranges=[(6, 12)], build_bloom="no")

        assert entry.name == "test-genome"
        assert entry.genome_size == 100000
        assert entry.gc_content == 0.45

        # Get
        retrieved = library.get("test-genome")
        assert retrieved is not None
        assert retrieved.name == "test-genome"

        # List
        entries = library.list()
        assert len(entries) == 1
        assert entries[0].name == "test-genome"

        # get_kmer_prefix
        prefix = library.get_kmer_prefix("test-genome")
        assert prefix is not None
        assert "test-genome" in prefix

        # Remove
        assert library.remove("test-genome") is True
        assert library.get("test-genome") is None
        assert library.list() == []

    @patch("neoswga.core.genome_library.GenomeLibrary._calculate_genome_stats")
    @patch("neoswga.core.kmer_counter.run_jellyfish")
    def test_has_kmers_for_range_with_files(self, mock_jf, mock_stats, tmp_path):
        """has_kmers_for_range returns True when k-mer files exist."""
        from neoswga.core.genome_library import GenomeLibrary

        mock_stats.return_value = (100000, 0.5)
        mock_jf.return_value = None

        fasta = tmp_path / "test.fna"
        fasta.write_text(">seq1\nATCG\n")

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        entry = library.add("test-g", str(fasta), k_ranges=[(6, 8)], build_bloom="no")

        # Create fake k-mer files
        for k in range(6, 9):
            kmer_file = f"{entry.kmer_prefix}_{k}mer_all.txt"
            with open(kmer_file, "w") as f:
                f.write("ATCG 1\n")

        assert library.has_kmers_for_range("test-g", 6, 8) is True
        assert library.has_kmers_for_range("test-g", 6, 12) is False  # Missing 9-12

    @patch("neoswga.core.genome_library.GenomeLibrary._calculate_genome_stats")
    @patch("neoswga.core.kmer_counter.run_jellyfish")
    def test_add_creates_symlink(self, mock_jf, mock_stats, tmp_path):
        """add() creates a symlink to the original FASTA."""
        from neoswga.core.genome_library import GenomeLibrary

        mock_stats.return_value = (50000, 0.4)
        mock_jf.return_value = None

        fasta = tmp_path / "original.fna"
        fasta.write_text(">seq1\nATCG\n")

        library = GenomeLibrary(library_dir=str(tmp_path / "genomes"))
        library.add("symlink-test", str(fasta), k_ranges=[(6, 8)], build_bloom="no")

        link = tmp_path / "genomes" / "symlink-test" / "genome.fna"
        assert link.is_symlink()
        assert os.path.realpath(str(link)) == str(fasta.resolve())
