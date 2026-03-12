"""Tests for neoswga.core.string_search module.

Covers the main search functions using synthetic genome sequences,
including forward/reverse complement searches, overlapping matches,
edge cases, circular genome handling, Aho-Corasick multi-k search,
per-k fallback search, and the genome caching layer.
"""

import pytest
from unittest.mock import patch, MagicMock

import neoswga.core.string_search as string_search_mod
from neoswga.core.string_search import (
    get_all_positions_per_k,
    get_all_positions_multi_k,
    get_cached_genome_sequence,
    clear_genome_cache,
    get_genome_cache_stats,
    write_to_h5py,
    AHOCORASICK_AVAILABLE,
)
from neoswga.core.thermodynamics import reverse_complement


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SIMPLE_GENOME = "ATCGATCGATCG"  # 12 bp, "ATCG" at positions 0, 4, 8
REPEAT_GENOME = "AAAAAAAAAA"    # 10 bp, all A
NO_MATCH_GENOME = "TTTTTTTTTT"  # 10 bp, no ATCG match
MIXED_GENOME = "ATCGNNNATCG"   # 11 bp with ambiguous bases


@pytest.fixture(autouse=True)
def _clear_cache():
    """Ensure a clean genome cache for every test."""
    string_search_mod._genome_cache.clear()
    yield
    string_search_mod._genome_cache.clear()


def _patch_genome(sequence: str):
    """Return a mock patch that makes get_cached_genome_sequence return *sequence*."""
    return patch(
        "neoswga.core.string_search.get_cached_genome_sequence",
        return_value=sequence,
    )


# ---------------------------------------------------------------------------
# get_all_positions_per_k -- basic behaviour
# ---------------------------------------------------------------------------

class TestGetAllPositionsPerK:
    """Tests for the per-k sliding-window search."""

    def test_find_known_positions(self):
        """ATCG should be found at positions 0, 4, 8 in ATCGATCGATCG."""
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_per_k(["ATCG"], "dummy.fna", circular=False)
        assert result["ATCG"] == [0, 4, 8]

    def test_multiple_kmers(self):
        """Search for several k-mers simultaneously."""
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_per_k(
                ["ATCG", "TCGA", "CGAT"], "dummy.fna", circular=False
            )
        assert result["ATCG"] == [0, 4, 8]
        assert result["TCGA"] == [1, 5]
        assert result["CGAT"] == [2, 6]

    def test_no_matches(self):
        """A k-mer absent from the genome should yield an empty list."""
        with _patch_genome(NO_MATCH_GENOME):
            result = get_all_positions_per_k(["ATCG"], "dummy.fna", circular=False)
        assert result["ATCG"] == []

    def test_overlapping_matches(self):
        """AA appears at every position in AAAAAAAAAA (overlapping)."""
        with _patch_genome(REPEAT_GENOME):
            result = get_all_positions_per_k(["AA"], "dummy.fna", circular=False)
        assert result["AA"] == list(range(9))  # positions 0-8

    def test_empty_kmer_list(self):
        """An empty primer list should return an empty dict."""
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_per_k([], "dummy.fna", circular=False)
        assert result == {}

    def test_kmer_longer_than_genome(self):
        """A k-mer longer than the genome should produce no matches."""
        short_genome = "ATCG"
        with _patch_genome(short_genome):
            result = get_all_positions_per_k(
                ["ATCGATCG"], "dummy.fna", circular=False
            )
        assert result["ATCGATCG"] == []

    def test_kmer_equals_genome(self):
        """A k-mer identical to the full genome should match at position 0."""
        genome = "ATCGATCG"
        with _patch_genome(genome):
            result = get_all_positions_per_k(["ATCGATCG"], "dummy.fna", circular=False)
        assert result["ATCGATCG"] == [0]

    def test_single_nucleotide_genome(self):
        """Single-base genome: only a single-base k-mer can match."""
        with _patch_genome("A"):
            result = get_all_positions_per_k(["A"], "dummy.fna", circular=False)
        assert result["A"] == [0]

        with _patch_genome("A"):
            result = get_all_positions_per_k(["AT"], "dummy.fna", circular=False)
        assert result["AT"] == []

    def test_empty_genome(self):
        """An empty genome should yield no matches for any k-mer."""
        with _patch_genome(""):
            result = get_all_positions_per_k(["A"], "dummy.fna", circular=False)
        assert result["A"] == []


# ---------------------------------------------------------------------------
# Circular genome handling
# ---------------------------------------------------------------------------

class TestCircularGenome:
    """Tests for circular (wrap-around) genome searches."""

    def test_circular_wraparound_match(self):
        """A k-mer that spans the origin of a circular genome should be found."""
        # Genome: CGATCGAT (8 bp). Circular wraps CG AT -> positions should
        # include the wrap.  "ATCG" starts at position 2 and 6 (wrap: AT from
        # end + CG from start).
        genome = "CGATCGAT"
        with _patch_genome(genome):
            result = get_all_positions_per_k(["ATCG"], "dummy.fna", circular=True)
        # Position 2 is a normal hit; position 6 wraps around (AT|CG)
        assert 2 in result["ATCG"]
        assert 6 in result["ATCG"]

    def test_linear_no_wraparound(self):
        """The wrap-around hit should NOT appear in linear mode."""
        genome = "CGATCGAT"
        with _patch_genome(genome):
            result = get_all_positions_per_k(["ATCG"], "dummy.fna", circular=False)
        assert 6 not in result["ATCG"]
        assert result["ATCG"] == [2]


# ---------------------------------------------------------------------------
# Reverse complement searches
# ---------------------------------------------------------------------------

class TestReverseComplementSearch:
    """Verify that reverse-complement primers are found correctly."""

    def test_reverse_complement_positions(self):
        """Searching the RC of ATCG (which is CGAT) should find CGAT positions."""
        rc = reverse_complement("ATCG")
        assert rc == "CGAT"
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_per_k([rc], "dummy.fna", circular=False)
        assert result["CGAT"] == [2, 6]

    def test_forward_and_rc_together(self):
        """Searching both a primer and its RC should each return independent hits."""
        primer = "ATCG"
        rc = reverse_complement(primer)
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_per_k(
                [primer, rc], "dummy.fna", circular=False
            )
        assert result[primer] == [0, 4, 8]
        assert result[rc] == [2, 6]

    def test_palindromic_kmer(self):
        """A palindromic k-mer (equal to its own RC) should still work."""
        palindrome = "AATT"
        assert reverse_complement(palindrome) == palindrome
        genome = "AATTGGAATT"
        with _patch_genome(genome):
            result = get_all_positions_per_k([palindrome], "dummy.fna", circular=False)
        assert result[palindrome] == [0, 6]


# ---------------------------------------------------------------------------
# Aho-Corasick multi-k search
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not AHOCORASICK_AVAILABLE,
    reason="pyahocorasick not installed",
)
class TestAhoCorasickMultiK:
    """Tests for the Aho-Corasick multi-k search path."""

    def test_basic_multi_k(self):
        """Multi-k search should find primers of different lengths."""
        primer_lists_by_k = {
            4: ["ATCG", "CGAT"],
            3: ["ATC", "TCG"],
        }
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_multi_k(
                primer_lists_by_k, "dummy.fna", circular=False
            )
        assert result["ATCG"] == [0, 4, 8]
        assert result["CGAT"] == [2, 6]
        assert result["ATC"] == [0, 4, 8]
        assert result["TCG"] == [1, 5, 9]

    def test_multi_k_no_matches(self):
        """Multi-k search with no hits should return empty lists."""
        primer_lists_by_k = {4: ["GGGG"]}
        with _patch_genome(SIMPLE_GENOME):
            result = get_all_positions_multi_k(
                primer_lists_by_k, "dummy.fna", circular=False
            )
        assert result["GGGG"] == []

    def test_multi_k_circular(self):
        """Multi-k with circular genome should find wrap-around hits."""
        genome = "CGATCGAT"
        primer_lists_by_k = {4: ["ATCG"]}
        with _patch_genome(genome):
            result = get_all_positions_multi_k(
                primer_lists_by_k, "dummy.fna", circular=True
            )
        assert 2 in result["ATCG"]
        assert 6 in result["ATCG"]


# ---------------------------------------------------------------------------
# Aho-Corasick vs per-k equivalence
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not AHOCORASICK_AVAILABLE,
    reason="pyahocorasick not installed",
)
class TestAhoCorasickEquivalence:
    """The Aho-Corasick path and the per-k fallback must produce identical results."""

    @pytest.mark.parametrize(
        "genome,circular",
        [
            (SIMPLE_GENOME, False),
            (SIMPLE_GENOME, True),
            (REPEAT_GENOME, False),
            ("ACGTACGTACGT", True),
            ("GATTACAGATTACA", False),
        ],
    )
    def test_equivalence(self, genome: str, circular: bool):
        """Both search paths should return exactly the same positions."""
        primers_4 = ["ATCG", "CGAT", "GATC", "TCGA"]
        primers_3 = ["ATC", "TCG", "CGA", "GAT"]

        with _patch_genome(genome):
            per_k_4 = get_all_positions_per_k(primers_4, "dummy.fna", circular=circular)
            per_k_3 = get_all_positions_per_k(primers_3, "dummy.fna", circular=circular)

        with _patch_genome(genome):
            multi_k = get_all_positions_multi_k(
                {4: primers_4, 3: primers_3}, "dummy.fna", circular=circular
            )

        for primer in primers_4:
            assert sorted(multi_k[primer]) == sorted(per_k_4[primer]), (
                f"Mismatch for {primer} in genome={genome!r}, circular={circular}"
            )
        for primer in primers_3:
            assert sorted(multi_k[primer]) == sorted(per_k_3[primer]), (
                f"Mismatch for {primer} in genome={genome!r}, circular={circular}"
            )


# ---------------------------------------------------------------------------
# Genome cache
# ---------------------------------------------------------------------------

class TestGenomeCache:
    """Tests for the module-level genome sequence cache."""

    def test_clear_cache(self):
        """clear_genome_cache should empty the cache."""
        string_search_mod._genome_cache["test_key"] = "ACGT"
        assert get_genome_cache_stats()["num_genomes"] == 1
        clear_genome_cache()
        assert get_genome_cache_stats()["num_genomes"] == 0

    def test_cache_stats(self):
        """get_genome_cache_stats should report correct totals."""
        string_search_mod._genome_cache["g1"] = "ACGT"
        string_search_mod._genome_cache["g2"] = "TTTTTTTT"
        stats = get_genome_cache_stats()
        assert stats["num_genomes"] == 2
        assert stats["total_bp"] == 12

    def test_cached_genome_returns_same_object(self):
        """A cached genome should be returned by reference (no copy)."""
        string_search_mod._genome_cache["f.fna"] = "ACGT"
        seq1 = get_cached_genome_sequence("f.fna")
        seq2 = get_cached_genome_sequence("f.fna")
        assert seq1 is seq2

    def test_get_cached_genome_reads_file_once(self):
        """The file should be read only on the first call; subsequent calls use the cache."""
        clear_genome_cache()
        fake_seq = "ACGTACGT"
        unique_file = f"_test_cache_reads_once_{id(self)}.fna"
        with patch("neoswga.core.utility.read_fasta_file") as mock_read:
            mock_read.return_value = iter([fake_seq])
            seq1 = get_cached_genome_sequence(unique_file)
            seq2 = get_cached_genome_sequence(unique_file)
        assert seq1 == fake_seq.upper()
        assert seq2 == fake_seq.upper()
        mock_read.assert_called_once_with(unique_file)


# ---------------------------------------------------------------------------
# write_to_h5py
# ---------------------------------------------------------------------------

class TestWriteToH5py:
    """Tests for HDF5 output writing."""

    def test_write_empty_dict(self, tmp_path):
        """Writing an empty dict should be a no-op (no file created)."""
        prefix = str(tmp_path / "test")
        write_to_h5py({}, prefix)
        import os
        # No file should be created
        assert not any(f.endswith(".h5") for f in os.listdir(tmp_path))

    def test_write_and_read_back(self, tmp_path):
        """Positions written to HDF5 should be recoverable."""
        import h5py
        prefix = str(tmp_path / "test")
        kmer_dict = {"ATCG": [0, 4, 8], "CGAT": [2, 6]}
        # Create the file first (write_to_h5py opens in r+ mode)
        h5_path = prefix + "_4mer_positions.h5"
        with h5py.File(h5_path, "w"):
            pass
        write_to_h5py(kmer_dict, prefix)
        with h5py.File(h5_path, "r") as f:
            assert list(f["ATCG"][:]) == [0, 4, 8]
            assert list(f["CGAT"][:]) == [2, 6]

    def test_overwrite_existing_dataset(self, tmp_path):
        """Overwriting an existing dataset should update the values."""
        import h5py
        prefix = str(tmp_path / "test")
        h5_path = prefix + "_4mer_positions.h5"
        # Write initial data
        with h5py.File(h5_path, "w") as f:
            f.create_dataset("ATCG", data=[0, 4])
        # Overwrite with new data (different length)
        write_to_h5py({"ATCG": [0, 4, 8]}, prefix)
        with h5py.File(h5_path, "r") as f:
            assert list(f["ATCG"][:]) == [0, 4, 8]
