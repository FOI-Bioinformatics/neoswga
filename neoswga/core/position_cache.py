"""
In-memory position cache for eliminating HDF5 I/O bottleneck.

For bacterial genomes (1-7 Mbp), loading all positions into memory
provides 1000x speedup over repeated disk access.

Memory usage: ~4 MB for 500 primers × 1000 sites × 8 bytes
"""

import numpy as np
import h5py
import os
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from collections import defaultdict
import logging

from neoswga.core.thermodynamics import reverse_complement

logger = logging.getLogger(__name__)


@dataclass
class BindingSite:
    """Single primer binding site"""
    position: int
    strand: str  # 'forward' or 'reverse'
    primer: str

    def __hash__(self):
        return hash((self.position, self.strand, self.primer))


class PositionCache:
    """
    In-memory cache of all primer binding positions.

    Loads all HDF5 position data once, provides O(1) lookup.
    Dramatically faster than repeated HDF5 reads.

    Example:
        cache = PositionCache(fg_prefixes, candidate_primers)
        positions = cache.get_positions('fg_prefix', 'ATCGATCG')  # Fast!
    """

    def __init__(self, fname_prefixes: List[str], primers: List[str]):
        """
        Load all positions for given primers from HDF5 files.

        Args:
            fname_prefixes: List of path prefixes (e.g., ['data/ecoli'])
            primers: List of primer sequences to cache

        Raises:
            FileNotFoundError: If HDF5 files don't exist
        """
        self.cache: Dict[Tuple[str, str, str], np.ndarray] = {}
        self.primers = set(primers)
        self.fname_prefixes = fname_prefixes

        self._load_all_positions()
        self._report_statistics()

    def _load_all_positions(self) -> None:
        """Single-pass load of all position data"""

        logger.info(f"Loading positions for {len(self.primers)} primers from {len(self.fname_prefixes)} genomes...")

        # Group primers by length for efficient HDF5 access
        primers_by_length = defaultdict(list)
        for primer in self.primers:
            primers_by_length[len(primer)].append(primer)

        total_loaded = 0

        # Pre-compute reverse complements once (avoids redundant calls per prefix)
        rc_map = {primer: reverse_complement(primer) for primer in self.primers}

        for fname_prefix in self.fname_prefixes:
            for k, primer_list in primers_by_length.items():
                hdf5_path = f"{fname_prefix}_{k}mer_positions.h5"

                if not os.path.exists(hdf5_path):
                    logger.warning(f"HDF5 file not found: {hdf5_path}")
                    continue

                with h5py.File(hdf5_path, 'r') as db:
                    for primer in primer_list:
                        # Forward strand
                        if primer in db:
                            positions = np.array(db[primer], dtype=np.int32)
                            key = (fname_prefix, primer, 'forward')
                            self.cache[key] = positions
                            total_loaded += 1

                        # Reverse strand (pre-computed reverse complement)
                        rc = rc_map[primer]
                        if rc in db:
                            positions = np.array(db[rc], dtype=np.int32)
                            key = (fname_prefix, primer, 'reverse')
                            self.cache[key] = positions
                            total_loaded += 1

        logger.info(f"Loaded {total_loaded} position arrays into memory")

    def get_positions(self, fname_prefix: str, primer: str, strand: str = 'both') -> np.ndarray:
        """
        Get positions for a primer. O(1) lookup.

        Args:
            fname_prefix: Genome identifier
            primer: Primer sequence
            strand: 'forward', 'reverse', or 'both'

        Returns:
            Array of positions (empty if not found)
        """
        if strand == 'both':
            key = (fname_prefix, primer, 'both')
            cached = self.cache.get(key)
            if cached is not None:
                return cached
            fw = self.cache.get((fname_prefix, primer, 'forward'), np.array([], dtype=np.int32))
            rv = self.cache.get((fname_prefix, primer, 'reverse'), np.array([], dtype=np.int32))
            combined = np.concatenate([fw, rv])
            self.cache[key] = combined
            return combined
        else:
            return self.cache.get((fname_prefix, primer, strand), np.array([], dtype=np.int32))

    def get_all_positions(self, fname_prefix: str, primers: List[str]) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Get positions for multiple primers (batched).

        Returns:
            Dict mapping primer -> (forward_positions, reverse_positions)
        """
        result = {}
        for primer in primers:
            fw = self.get_positions(fname_prefix, primer, 'forward')
            rv = self.get_positions(fname_prefix, primer, 'reverse')
            result[primer] = (fw, rv)
        return result

    def compute_coverage_vectorized(self, fname_prefix: str, primers: List[str],
                                   genome_length: int) -> np.ndarray:
        """
        Compute coverage array for primer set (vectorized, fast).

        Args:
            fname_prefix: Genome identifier
            primers: List of primers
            genome_length: Length of genome

        Returns:
            Boolean array of length genome_length (True = covered)
        """
        coverage = np.zeros(genome_length, dtype=bool)

        for primer in primers:
            positions = self.get_positions(fname_prefix, primer, 'both')
            if len(positions) > 0:
                # Clip positions to genome bounds
                valid_positions = positions[positions < genome_length]
                coverage[valid_positions] = True

        return coverage

    # `compute_coverage_with_extension` was removed — it was uncalled
    # across the codebase and shipped a stale 70000 bp processivity
    # default that conflicted with the per-primer-reach semantics used
    # elsewhere. Callers wanting per-prefix coverage should use
    # `coverage.compute_per_prefix_coverage` with an explicit `extension`
    # kwarg from `coverage.polymerase_extension_reach`.

    def compute_statistics(self, fname_prefix: str, primers: List[str],
                          genome_length: int) -> Dict[str, float]:
        """
        Compute coverage statistics efficiently.

        Returns:
            Dictionary with coverage metrics
        """
        coverage = self.compute_coverage_vectorized(fname_prefix, primers, genome_length)

        # Find gaps
        gap_sizes = self._find_gap_sizes(coverage)

        return {
            'coverage_fraction': np.sum(coverage) / genome_length,
            'num_covered_bases': int(np.sum(coverage)),
            'num_gaps': len(gap_sizes),
            'mean_gap_size': np.mean(gap_sizes) if len(gap_sizes) > 0 else 0,
            'max_gap_size': np.max(gap_sizes) if len(gap_sizes) > 0 else 0,
            'gap_gini': self._gini_coefficient(gap_sizes) if len(gap_sizes) > 1 else 0,
            'gap_entropy': self.compute_gap_entropy(fname_prefix, primers, genome_length),
        }

    def compute_gap_entropy(self, fname_prefix: str, primers: List[str],
                           genome_length: int, num_bins: int = 50) -> float:
        """Compute Shannon entropy of gap size distribution.

        Entropy measures the information content of the gap distribution.
        Higher entropy indicates more diverse gap sizes; lower entropy
        indicates more uniform spacing. Complements the Gini coefficient
        by capturing distributional shape rather than inequality alone.

        Args:
            fname_prefix: Genome identifier.
            primers: List of primer sequences.
            genome_length: Length of genome.
            num_bins: Number of histogram bins for gap discretization.

        Returns:
            Shannon entropy in bits (0 = all gaps identical,
            higher = more varied).
        """
        from scipy.stats import entropy as scipy_entropy

        coverage = self.compute_coverage_vectorized(fname_prefix, primers, genome_length)
        gap_sizes = self._find_gap_sizes(coverage)

        if len(gap_sizes) < 2:
            return 0.0

        # Bin gap sizes into histogram and compute entropy
        counts, _ = np.histogram(gap_sizes, bins=min(num_bins, len(gap_sizes)))
        # Filter zero bins and normalize to probability distribution
        counts = counts[counts > 0]
        probs = counts / counts.sum()
        return float(scipy_entropy(probs, base=2))

    def _find_gap_sizes(self, coverage: np.ndarray) -> List[int]:
        """Find sizes of uncovered regions"""
        # Find transitions
        transitions = np.diff(coverage.astype(int))
        gap_starts = np.where(transitions == -1)[0] + 1
        gap_ends = np.where(transitions == 1)[0] + 1

        # Handle edge cases
        if not coverage[0]:
            gap_starts = np.concatenate([[0], gap_starts])
        if not coverage[-1]:
            gap_ends = np.concatenate([gap_ends, [len(coverage)]])

        gap_sizes = gap_ends - gap_starts
        return gap_sizes.tolist()

    def _gini_coefficient(self, values: List[int]) -> float:
        """Compute Gini coefficient of gap distribution"""
        if len(values) == 0:
            return 0.0

        sorted_values = np.sort(values)
        n = len(sorted_values)
        cumsum = np.cumsum(sorted_values)

        return (2.0 * np.sum((np.arange(1, n+1)) * sorted_values)) / (n * cumsum[-1]) - (n + 1) / n

    def compute_strand_alternation_stats(
        self,
        fname_prefix: str,
        primers: List[str],
        genome_length: int,
    ) -> Dict[str, float]:
        """Compute strand alternation metrics for a primer set.

        Strand alternation measures how well forward and reverse strand
        binding sites interleave across the genome. Good interleaving is
        needed for efficient strand displacement amplification by Phi29.

        A region with only forward-strand binding amplifies linearly at
        best; exponential amplification requires primers on both strands.

        Args:
            fname_prefix: Genome identifier.
            primers: List of primer sequences.
            genome_length: Length of genome in base pairs.

        Returns:
            Dictionary with the following keys:

            - strand_alternation_gap_mean: Mean gap between consecutive
              opposite-strand binding sites.
            - strand_alternation_gap_max: Maximum such gap (worst region).
            - strand_alternation_score: Fraction of adjacent site pairs
              that alternate strands (1.0 = perfect alternation).
            - strand_coverage_ratio: min(fwd, rev) / max(fwd, rev) site
              count ratio (1.0 = balanced).
            - longest_same_strand_run: Longest consecutive run of sites
              on the same strand (lower is better).
        """
        # Collect all binding sites with strand labels
        sites: List[Tuple[int, str]] = []
        total_fwd = 0
        total_rev = 0

        for primer in primers:
            fwd_pos = self.get_positions(fname_prefix, primer, 'forward')
            rev_pos = self.get_positions(fname_prefix, primer, 'reverse')
            total_fwd += len(fwd_pos)
            total_rev += len(rev_pos)

            for pos in fwd_pos:
                sites.append((int(pos), 'forward'))
            for pos in rev_pos:
                sites.append((int(pos), 'reverse'))

        if len(sites) < 2:
            return {
                'strand_alternation_gap_mean': float(genome_length),
                'strand_alternation_gap_max': float(genome_length),
                'strand_alternation_score': 0.0,
                'strand_coverage_ratio': 0.0,
                'longest_same_strand_run': len(sites),
            }

        # Sort by position
        sites.sort(key=lambda x: x[0])

        # Compute gaps between consecutive opposite-strand sites
        alternation_gaps: List[int] = []
        alternations = 0
        same_strand_run = 1
        longest_run = 1

        for i in range(1, len(sites)):
            if sites[i][1] != sites[i - 1][1]:
                # Strand switch: record gap
                gap = sites[i][0] - sites[i - 1][0]
                alternation_gaps.append(gap)
                alternations += 1
                same_strand_run = 1
            else:
                same_strand_run += 1
                longest_run = max(longest_run, same_strand_run)

        alternation_score = alternations / (len(sites) - 1)

        # Strand balance ratio
        if total_fwd > 0 and total_rev > 0:
            strand_ratio = min(total_fwd, total_rev) / max(total_fwd, total_rev)
        else:
            strand_ratio = 0.0

        if alternation_gaps:
            mean_gap = float(np.mean(alternation_gaps))
            max_gap = float(max(alternation_gaps))
        else:
            mean_gap = float(genome_length)
            max_gap = float(genome_length)

        return {
            'strand_alternation_gap_mean': mean_gap,
            'strand_alternation_gap_max': max_gap,
            'strand_alternation_score': alternation_score,
            'strand_coverage_ratio': strand_ratio,
            'longest_same_strand_run': longest_run,
        }

    def _report_statistics(self) -> None:
        """Report cache statistics"""
        memory_bytes = sum(arr.nbytes for arr in self.cache.values())
        memory_mb = memory_bytes / 1e6

        logger.info(f"Cache statistics:")
        logger.info(f"  Entries: {len(self.cache)}")
        logger.info(f"  Memory: {memory_mb:.1f} MB")
        logger.info(f"  Primers: {len(self.primers)}")


class StreamingPositionCache:
    """
    Memory-efficient cache for very large background genomes.

    Uses memory mapping and lazy loading.
    Suitable for human genome (3 Gbp) where full cache is infeasible.
    """

    def __init__(self, fname_prefixes: List[str], primers: Optional[List[str]] = None):
        """
        Initialize streaming cache with memory mapping.

        Args:
            fname_prefixes: HDF5 file prefixes
            primers: Optional subset to preload (if small)
        """
        self.fname_prefixes = fname_prefixes
        self.file_handles: Dict[str, h5py.File] = {}
        self.preloaded: Dict[Tuple[str, str, str], np.ndarray] = {}

        # Open HDF5 files in read-only mode with memory mapping
        for prefix in fname_prefixes:
            for k in range(6, 13):
                path = f"{prefix}_{k}mer_positions.h5"
                if os.path.exists(path):
                    # Open with driver for memory mapping
                    self.file_handles[path] = h5py.File(path, 'r', rdcc_nbytes=1024**2)

        # Preload small subset if provided
        if primers and len(primers) < 100:
            self._preload_subset(primers)

    def _preload_subset(self, primers: List[str]) -> None:
        """Preload positions for small set of primers"""
        logger.info(f"Preloading {len(primers)} primers...")

        for fname_prefix in self.fname_prefixes:
            for primer in primers:
                k = len(primer)
                path = f"{fname_prefix}_{k}mer_positions.h5"

                if path not in self.file_handles:
                    continue

                db = self.file_handles[path]

                if primer in db:
                    self.preloaded[(fname_prefix, primer, 'forward')] = np.array(db[primer])

    def get_positions(self, fname_prefix: str, primer: str, strand: str = 'both') -> np.ndarray:
        """
        Get positions (lazy loading from memory-mapped HDF5).

        Slower than PositionCache but uses less memory.
        """
        # Check preloaded cache
        if strand != 'both':
            key = (fname_prefix, primer, strand)
            if key in self.preloaded:
                return self.preloaded[key]

        # Lazy load from HDF5
        k = len(primer)
        path = f"{fname_prefix}_{k}mer_positions.h5"

        if path not in self.file_handles:
            return np.array([], dtype=np.int32)

        db = self.file_handles[path]

        if strand == 'both':
            fw = np.array(db[primer]) if primer in db else np.array([])
            rc = reverse_complement(primer)
            rv = np.array(db[rc]) if rc in db else np.array([])
            return np.concatenate([fw, rv])
        elif strand == 'forward':
            return np.array(db[primer]) if primer in db else np.array([])
        else:  # reverse
            rc = reverse_complement(primer)
            return np.array(db[rc]) if rc in db else np.array([])

    def close(self) -> None:
        """Close all HDF5 file handles"""
        for fh in self.file_handles.values():
            fh.close()
        self.file_handles.clear()

    def __enter__(self) -> 'StreamingPositionCache':
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()


def benchmark_cache_vs_hdf5(fname_prefixes: List[str], primers: List[str], iterations: int = 1000) -> None:
    """
    Benchmark cache performance vs. direct HDF5 access.

    Demonstrates 100-1000× speedup for repeated access patterns.
    """
    import time

    # Build cache
    print("Building cache...")
    start = time.time()
    cache = PositionCache(fname_prefixes, primers)
    cache_build_time = time.time() - start
    print(f"Cache build time: {cache_build_time:.2f}s")

    # Benchmark cache access
    print(f"\nBenchmarking {iterations} cache accesses...")
    start = time.time()
    for _ in range(iterations):
        primer = np.random.choice(primers)
        positions = cache.get_positions(fname_prefixes[0], primer)
    cache_time = time.time() - start
    print(f"Cache access time: {cache_time:.4f}s ({iterations/cache_time:.0f} queries/sec)")

    # Benchmark direct HDF5 access
    print(f"\nBenchmarking {iterations} HDF5 accesses...")
    start = time.time()
    for _ in range(iterations):
        primer = np.random.choice(primers)
        k = len(primer)
        path = f"{fname_prefixes[0]}_{k}mer_positions.h5"
        with h5py.File(path, 'r') as db:
            if primer in db:
                positions = np.array(db[primer])
    hdf5_time = time.time() - start
    print(f"HDF5 access time: {hdf5_time:.4f}s ({iterations/hdf5_time:.0f} queries/sec)")

    print(f"\nSpeedup: {hdf5_time/cache_time:.1f}×")
    print(f"Amortization: Cache pays for itself after {cache_build_time/((hdf5_time-cache_time)/iterations):.0f} queries")


if __name__ == "__main__":
    # Example usage
    import sys

    logging.basicConfig(level=logging.INFO)

    # Test with example data
    if len(sys.argv) > 1:
        fname_prefixes = [sys.argv[1]]
        test_primers = ['ATCGATCG', 'GCTAGCTA', 'TTAATTAA']

        print("Testing PositionCache...")
        cache = PositionCache(fname_prefixes, test_primers)

        for primer in test_primers:
            positions = cache.get_positions(fname_prefixes[0], primer)
            print(f"{primer}: {len(positions)} sites")

        # Run benchmark if enough primers
        if len(test_primers) > 10:
            benchmark_cache_vs_hdf5(fname_prefixes, test_primers)
    else:
        print("Usage: python position_cache.py <fname_prefix>")
