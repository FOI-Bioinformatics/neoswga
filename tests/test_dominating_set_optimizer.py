"""
Unit tests for dominating_set_optimizer module.

Tests:
- CoverageRegion dataclass
- BipartiteGraph class
- DominatingSetOptimizer class
- Greedy set cover algorithm
"""

import pytest
import numpy as np
from unittest.mock import Mock, MagicMock, patch
from typing import Dict, Set

from neoswga.core.dominating_set_optimizer import (
    CoverageRegion,
    BipartiteGraph,
    DominatingSetOptimizer,
)


# =============================================================================
# CoverageRegion Tests
# =============================================================================

class TestCoverageRegion:
    """Tests for CoverageRegion dataclass."""

    def test_creation(self):
        """Test creating a CoverageRegion."""
        region = CoverageRegion(
            chromosome="chr1",
            start=0,
            end=10000,
            covered_by=set()
        )

        assert region.chromosome == "chr1"
        assert region.start == 0
        assert region.end == 10000
        assert region.covered_by == set()

    def test_hashable(self):
        """Test that CoverageRegion is hashable."""
        region1 = CoverageRegion("chr1", 0, 10000, set())
        region2 = CoverageRegion("chr1", 0, 10000, set())
        region3 = CoverageRegion("chr1", 10000, 20000, set())

        # Same chromosome/start/end should hash the same
        assert hash(region1) == hash(region2)
        # Different positions should hash differently
        assert hash(region1) != hash(region3)

    def test_can_be_used_in_set(self):
        """Test that CoverageRegions can be used in a set."""
        region1 = CoverageRegion("chr1", 0, 10000, set())
        region2 = CoverageRegion("chr1", 0, 10000, {"ATCG"})  # Same position, different coverage
        region3 = CoverageRegion("chr1", 10000, 20000, set())

        region_set = {region1, region2, region3}

        # region1 and region2 have same hash (same position), so set depends on equality
        # They should be considered equal by hash
        assert len(region_set) == 2 or len(region_set) == 3  # depends on __eq__

    def test_covered_by_modification(self):
        """Test that covered_by can be modified."""
        region = CoverageRegion("chr1", 0, 10000, set())

        region.covered_by.add("ATCGATCG")
        region.covered_by.add("GCTAGCTA")

        assert "ATCGATCG" in region.covered_by
        assert "GCTAGCTA" in region.covered_by
        assert len(region.covered_by) == 2


# =============================================================================
# BipartiteGraph Tests
# =============================================================================

class TestBipartiteGraph:
    """Tests for BipartiteGraph class."""

    def test_initialization(self):
        """Test graph initialization."""
        graph = BipartiteGraph(bin_size=10000)

        assert graph.bin_size == 10000
        assert len(graph.primers) == 0
        assert len(graph.regions) == 0
        assert len(graph.primer_to_regions) == 0
        assert len(graph.region_to_primers) == 0

    def test_add_primer_coverage_single_position(self):
        """Test adding a primer with a single position."""
        graph = BipartiteGraph(bin_size=10000)

        positions = np.array([5000])
        graph.add_primer_coverage("ATCGATCG", positions, "genome1", 100000)

        assert "ATCGATCG" in graph.primers
        assert len(graph.regions) == 1
        assert "ATCGATCG" in graph.primer_to_regions

    def test_add_primer_coverage_multiple_positions(self):
        """Test adding a primer with multiple positions across bins."""
        graph = BipartiteGraph(bin_size=10000)

        # Positions in different bins
        positions = np.array([5000, 15000, 25000, 35000])
        graph.add_primer_coverage("ATCGATCG", positions, "genome1", 100000)

        assert "ATCGATCG" in graph.primers
        assert len(graph.regions) == 4  # 4 different bins

        # Each region should be covered by the primer
        for region in graph.regions:
            assert "ATCGATCG" in region.covered_by

    def test_add_multiple_primers(self):
        """Test adding multiple primers with overlapping coverage."""
        graph = BipartiteGraph(bin_size=10000)

        # First primer covers bins 0-1
        positions1 = np.array([5000, 15000])
        graph.add_primer_coverage("PRIMER1", positions1, "genome1", 100000)

        # Second primer covers bins 1-2 (overlaps bin 1)
        positions2 = np.array([15000, 25000])
        graph.add_primer_coverage("PRIMER2", positions2, "genome1", 100000)

        assert len(graph.primers) == 2
        assert len(graph.regions) == 3  # bins 0, 1, 2

        # Find bin 1 region (the overlapping one)
        for region in graph.regions:
            if region.start == 10000:  # bin 1
                assert "PRIMER1" in region.covered_by
                assert "PRIMER2" in region.covered_by

    def test_get_uncovered_regions_all_covered(self):
        """Test get_uncovered_regions when all regions are covered."""
        graph = BipartiteGraph(bin_size=10000)

        positions = np.array([5000, 15000])
        graph.add_primer_coverage("PRIMER1", positions, "genome1", 100000)

        uncovered = graph.get_uncovered_regions({"PRIMER1"})
        assert len(uncovered) == 0

    def test_get_uncovered_regions_none_covered(self):
        """Test get_uncovered_regions when no primers selected."""
        graph = BipartiteGraph(bin_size=10000)

        positions = np.array([5000, 15000])
        graph.add_primer_coverage("PRIMER1", positions, "genome1", 100000)

        uncovered = graph.get_uncovered_regions(set())
        assert len(uncovered) == 2  # Both regions uncovered

    def test_get_uncovered_regions_partial_coverage(self):
        """Test get_uncovered_regions with partial coverage."""
        graph = BipartiteGraph(bin_size=10000)

        # Two primers, each covering different regions
        positions1 = np.array([5000])
        graph.add_primer_coverage("PRIMER1", positions1, "genome1", 100000)

        positions2 = np.array([15000])
        graph.add_primer_coverage("PRIMER2", positions2, "genome1", 100000)

        # Only select PRIMER1
        uncovered = graph.get_uncovered_regions({"PRIMER1"})
        assert len(uncovered) == 1  # Region covered by PRIMER2 is uncovered

    def test_get_coverage_score_empty_graph(self):
        """Test coverage score for empty graph."""
        graph = BipartiteGraph(bin_size=10000)

        score = graph.get_coverage_score(set())
        assert score == 0.0

    def test_get_coverage_score_full_coverage(self):
        """Test coverage score for full coverage."""
        graph = BipartiteGraph(bin_size=10000)

        positions = np.array([5000, 15000])
        graph.add_primer_coverage("PRIMER1", positions, "genome1", 100000)

        score = graph.get_coverage_score({"PRIMER1"})
        assert score == 1.0

    def test_get_coverage_score_partial_coverage(self):
        """Test coverage score for partial coverage."""
        graph = BipartiteGraph(bin_size=10000)

        # Two primers, each covering one region
        graph.add_primer_coverage("PRIMER1", np.array([5000]), "genome1", 100000)
        graph.add_primer_coverage("PRIMER2", np.array([15000]), "genome1", 100000)

        # Only select PRIMER1 - should get 50% coverage
        score = graph.get_coverage_score({"PRIMER1"})
        assert score == 0.5


# =============================================================================
# DominatingSetOptimizer Tests
# =============================================================================

class TestDominatingSetOptimizer:
    """Tests for DominatingSetOptimizer class."""

    @pytest.fixture
    def mock_cache(self):
        """Create a mock PositionCache."""
        cache = Mock()
        return cache

    @pytest.fixture
    def simple_optimizer(self, mock_cache):
        """Create optimizer with simple mock data."""
        # Configure mock to return positions for specific primers
        def get_positions(prefix, primer, strand):
            # Return positions based on primer
            positions_map = {
                "PRIMER1": np.array([5000, 15000, 25000]),  # Covers bins 0, 1, 2
                "PRIMER2": np.array([35000, 45000]),  # Covers bins 3, 4
                "PRIMER3": np.array([15000, 35000]),  # Overlaps bins 1 and 3
                "PRIMER4": np.array([55000]),  # Covers bin 5
            }
            return positions_map.get(primer, np.array([]))

        mock_cache.get_positions = get_positions

        optimizer = DominatingSetOptimizer(
            cache=mock_cache,
            fg_prefixes=["genome1"],
            fg_seq_lengths=[100000],
            bin_size=10000
        )
        return optimizer

    def test_initialization(self, mock_cache):
        """Test optimizer initialization."""
        optimizer = DominatingSetOptimizer(
            cache=mock_cache,
            fg_prefixes=["genome1", "genome2"],
            fg_seq_lengths=[100000, 50000],
            bin_size=5000
        )

        assert optimizer.cache == mock_cache
        assert optimizer.fg_prefixes == ["genome1", "genome2"]
        assert optimizer.fg_seq_lengths == [100000, 50000]
        assert optimizer.bin_size == 5000

    def test_optimize_greedy_basic(self, simple_optimizer):
        """Test basic greedy optimization."""
        candidates = ["PRIMER1", "PRIMER2", "PRIMER3", "PRIMER4"]

        result = simple_optimizer.optimize_greedy(
            candidates,
            max_primers=10,
            verbose=False
        )

        # Should return a valid result dict
        assert "primers" in result
        assert "n_primers" in result
        assert "coverage" in result
        assert "covered_regions" in result
        assert "total_regions" in result

        # Should have selected some primers
        assert result["n_primers"] > 0
        assert len(result["primers"]) == result["n_primers"]

    def test_optimize_greedy_achieves_coverage(self, simple_optimizer):
        """Test that greedy optimization achieves coverage."""
        candidates = ["PRIMER1", "PRIMER2", "PRIMER4"]  # Can cover bins 0-5

        result = simple_optimizer.optimize_greedy(
            candidates,
            max_primers=10,
            verbose=False
        )

        # Should achieve full coverage
        assert result["coverage"] == 1.0
        assert result["uncovered_regions"] == 0

    def test_optimize_greedy_respects_max_primers(self, simple_optimizer):
        """Test that greedy respects max_primers limit."""
        candidates = ["PRIMER1", "PRIMER2", "PRIMER3", "PRIMER4"]

        result = simple_optimizer.optimize_greedy(
            candidates,
            max_primers=2,
            verbose=False
        )

        # Should not exceed max_primers
        assert result["n_primers"] <= 2

    def test_optimize_greedy_selects_efficiently(self, simple_optimizer):
        """Test that greedy selects primers efficiently."""
        # PRIMER1 covers 3 regions, PRIMER3 only helps with overlap
        candidates = ["PRIMER1", "PRIMER3"]

        result = simple_optimizer.optimize_greedy(
            candidates,
            max_primers=5,
            verbose=False
        )

        # PRIMER1 should be selected first (covers more)
        assert "PRIMER1" in result["primers"]

    def test_optimize_greedy_empty_candidates(self, mock_cache):
        """Test greedy with empty candidates list."""
        mock_cache.get_positions = Mock(return_value=np.array([]))

        optimizer = DominatingSetOptimizer(
            cache=mock_cache,
            fg_prefixes=["genome1"],
            fg_seq_lengths=[100000],
            bin_size=10000
        )

        result = optimizer.optimize_greedy([], max_primers=10, verbose=False)

        assert result["n_primers"] == 0
        assert result["coverage"] == 0.0

    def test_optimize_greedy_no_coverage_primers(self, mock_cache):
        """Test greedy when no primers have positions."""
        mock_cache.get_positions = Mock(return_value=np.array([]))

        optimizer = DominatingSetOptimizer(
            cache=mock_cache,
            fg_prefixes=["genome1"],
            fg_seq_lengths=[100000],
            bin_size=10000
        )

        result = optimizer.optimize_greedy(
            ["PRIMER1", "PRIMER2"],
            max_primers=10,
            verbose=False
        )

        assert result["n_primers"] == 0
        assert result["total_regions"] == 0


# =============================================================================
# Integration Tests
# =============================================================================

class TestDominatingSetIntegration:
    """Integration tests for dominating set optimization."""

    def test_full_workflow_mock(self):
        """Test complete optimization workflow with mock cache."""
        # Create mock cache
        cache = Mock()

        def get_positions(prefix, primer, strand):
            # Simulate primers covering different parts of genome
            if primer == "ATCGATCG":
                return np.array([1000, 2000, 3000])
            elif primer == "GCTAGCTA":
                return np.array([4000, 5000, 6000])
            elif primer == "AAATTTCCC":
                return np.array([7000, 8000, 9000])
            return np.array([])

        cache.get_positions = get_positions

        optimizer = DominatingSetOptimizer(
            cache=cache,
            fg_prefixes=["test_genome"],
            fg_seq_lengths=[10000],
            bin_size=1000  # Small bins
        )

        candidates = ["ATCGATCG", "GCTAGCTA", "AAATTTCCC"]
        result = optimizer.optimize_greedy(
            candidates,
            max_primers=5,
            verbose=False
        )

        # Should cover all regions
        assert result["coverage"] == 1.0
        # All 3 primers needed (they cover non-overlapping regions)
        assert result["n_primers"] == 3

    def test_overlapping_coverage_efficiency(self):
        """Test that optimizer efficiently handles overlapping coverage."""
        cache = Mock()

        def get_positions(prefix, primer, strand):
            # MEGA_PRIMER covers everything
            if primer == "MEGA_PRIMER":
                return np.array([500, 1500, 2500, 3500, 4500])
            # SMALL_1 only covers bin 0
            elif primer == "SMALL_1":
                return np.array([500])
            # SMALL_2 only covers bin 1
            elif primer == "SMALL_2":
                return np.array([1500])
            return np.array([])

        cache.get_positions = get_positions

        optimizer = DominatingSetOptimizer(
            cache=cache,
            fg_prefixes=["genome"],
            fg_seq_lengths=[5000],
            bin_size=1000
        )

        # Should prefer MEGA_PRIMER since it covers everything
        result = optimizer.optimize_greedy(
            ["MEGA_PRIMER", "SMALL_1", "SMALL_2"],
            max_primers=10,
            verbose=False
        )

        # MEGA_PRIMER should be selected first and achieve full coverage
        assert "MEGA_PRIMER" in result["primers"]
        assert result["n_primers"] == 1
        assert result["coverage"] == 1.0


# =============================================================================
# Edge Cases
# =============================================================================

class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_single_bin_genome(self):
        """Test with a genome smaller than bin size."""
        cache = Mock()
        cache.get_positions = Mock(return_value=np.array([500]))

        optimizer = DominatingSetOptimizer(
            cache=cache,
            fg_prefixes=["tiny_genome"],
            fg_seq_lengths=[5000],  # Smaller than default bin
            bin_size=10000
        )

        result = optimizer.optimize_greedy(
            ["PRIMER1"],
            max_primers=5,
            verbose=False
        )

        # Should create one region and cover it
        assert result["total_regions"] == 1
        assert result["coverage"] == 1.0

    def test_large_position_values(self):
        """Test with large position values."""
        cache = Mock()
        cache.get_positions = Mock(return_value=np.array([
            1000000, 2000000, 3000000, 4000000
        ]))

        optimizer = DominatingSetOptimizer(
            cache=cache,
            fg_prefixes=["large_genome"],
            fg_seq_lengths=[5000000],
            bin_size=100000
        )

        result = optimizer.optimize_greedy(
            ["PRIMER1"],
            max_primers=5,
            verbose=False
        )

        # Should handle large positions correctly
        assert result["total_regions"] == 4
        assert result["coverage"] == 1.0

    def test_multiple_genomes(self):
        """Test with multiple target genomes."""
        cache = Mock()

        def get_positions(prefix, primer, strand):
            if prefix == "genome1":
                return np.array([1000, 2000])
            elif prefix == "genome2":
                return np.array([500, 1500])
            return np.array([])

        cache.get_positions = get_positions

        optimizer = DominatingSetOptimizer(
            cache=cache,
            fg_prefixes=["genome1", "genome2"],
            fg_seq_lengths=[5000, 3000],
            bin_size=1000
        )

        result = optimizer.optimize_greedy(
            ["PRIMER1"],
            max_primers=5,
            verbose=False
        )

        # Should cover regions in both genomes
        assert result["total_regions"] == 4  # 2 bins in each genome
        assert result["coverage"] == 1.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
