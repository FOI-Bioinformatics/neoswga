"""
Efficient background genome filtering using probabilistic data structures.

For massive genomes (human 3 Gbp, tick 2.1 Gbp), exact counting is infeasible.
Use Bloom filters and sampling for fast negative selection.
"""

import logging
import os
import pickle
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

try:
    from pybloom_live import BloomFilter
except ImportError:
    print("WARNING: pybloom_live not installed. Install with: pip install pybloom-live")
    BloomFilter = None

logger = logging.getLogger(__name__)


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _reverse_complement(seq: str) -> str:
    """Reverse complement of an ACGT string."""
    return seq.translate(_COMPLEMENT)[::-1]


def distinct_kmer_capacity(
    genome_size: int, min_k: int = 6, max_k: int = 12, margin: float = 1.1
) -> int:
    """Upper bound on the DISTINCT k-mers a genome holds, over min_k..max_k.

    pybloom allocates its bit array upfront from `capacity`, and its `add`
    increments the count only for an item the filter did not already hold. So
    this is the quantity capacity must bound.

    The previous heuristic, ten times the genome's base count, bounded
    INSERTIONS instead: every position contributes one k-mer per length, and
    the multiplier was chosen for those seven lengths. On hg38 that asked for
    33e9 capacity, about 39.6 GB of bits, where the distinct count is 22.4
    million and about 27 MB. It could not be allocated on the one background
    this module exists for, while every small-genome test passed.

    Each term saturates at 4**k, which is why a 144 Mb and a 3.3 Gb genome
    return the same bound: above the k-mer space a longer genome cannot hold
    more distinct k-mers. The margin covers pybloom raising AT capacity rather
    than above it, and the floor keeps the result usable, since a capacity of
    zero is rejected.
    """
    total = 0
    for k in range(min_k, max_k + 1):
        positions = max(0, genome_size - k + 1)
        total += min(4**k, positions)
    return max(1, int(total * margin))


@dataclass
class BackgroundFilterConfig:
    """Configuration for background filtering"""

    max_exact_matches: int = 10  # Max perfect matches in background
    # Bounded by 1 + 3k neighbours REPORTED PRESENT, not by a site count: 37
    # for a 12-mer. See BackgroundBloomFilter.count_present_neighbours.
    max_1mm_matches: int = 100
    bloom_fp_rate: float = 0.01  # Bloom filter false positive rate
    sample_rate: int = 100  # For sampled suffix array
    use_repeat_filter: bool = True
    # The oligo lengths the filter must answer for. A design using a length
    # outside this range reads as absent from the background and clears the
    # gate unscreened, so these must cover the design's own min_k..max_k.
    min_k: int = 6
    max_k: int = 12


class BackgroundBloomFilter:
    """
    Bloom filter for fast background genome screening.

    Memory efficient: 3 Gbp genome → ~4 GB Bloom filter
    vs. ~170 GB exact HDF5 index.

    Query time: O(1) per primer (vs. O(genome_size) for exact search)
    """

    def __init__(self, capacity: int = 3e9, error_rate: float = 0.01):
        """
        Initialize Bloom filter.

        Args:
            capacity: Expected number of k-mers (e.g., 3e9 for human)
            error_rate: False positive rate (default 1%)
        """
        if BloomFilter is None:
            raise ImportError("pybloom_live required. Install: pip install pybloom-live")

        self.bloom = BloomFilter(capacity=int(capacity), error_rate=error_rate)
        self.kmer_count = 0
        self.genome_size = 0
        # The oligo lengths this filter was built over. None means the range
        # was never recorded, which is what an artifact written before this
        # field looks like -- unknown, not empty.
        self.min_k: Optional[int] = None
        self.max_k: Optional[int] = None

    def _record_range(self, min_k: int, max_k: int) -> None:
        """Widen the recorded range; two builds on one filter give the union."""
        self.min_k = min_k if self.min_k is None else min(self.min_k, min_k)
        self.max_k = max_k if self.max_k is None else max(self.max_k, max_k)

    def covers_length(self, k: int) -> Optional[bool]:
        """Whether this filter can answer for a k-mer of length k.

        None when the range was never recorded, which a caller must treat as
        unknown rather than as either answer.
        """
        if self.min_k is None or self.max_k is None:
            return None
        return self.min_k <= k <= self.max_k

    def add(self, kmer: str):
        """
        Add a single k-mer to the Bloom filter.

        Args:
            kmer: K-mer sequence to add (must contain only ATCG)
        """
        if self._is_valid_kmer(kmer):
            self.bloom.add(kmer)
            self.kmer_count += 1

    def add_genome(
        self,
        fasta_path: str,
        include_mismatches: bool = False,
        min_k: int = 6,
        max_k: int = 12,
        chunk_size: int = 1000000,
    ):
        """
        Add all k-mers from genome to Bloom filter.

        Optimized for large genomes with:
        - Batch k-mer generation (reduces Python overhead)
        - Progress reporting with tqdm
        - Single-pass sequence processing per k-mer length
        - Optional mismatch variants (disabled by default for speed)

        Args:
            fasta_path: Path to FASTA file
            include_mismatches: Add 1-mismatch variants (slow, default False)
            min_k: Minimum k-mer length (default 6)
            max_k: Maximum k-mer length (default 12)
            chunk_size: Positions to process per progress update
        """
        logger.info(f"Adding genome to Bloom filter: {fasta_path}")
        logger.info(f"  K-mer range: {min_k}-{max_k}bp")
        self._record_range(min_k, max_k)

        from Bio import SeqIO

        try:
            from tqdm import tqdm

            use_tqdm = True
        except ImportError:
            use_tqdm = False

        # First pass: count total sequence length for progress bar
        total_length = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            total_length += len(record.seq)

        logger.info(f"  Total genome length: {total_length:,} bp")

        # Precompile the valid bases set for faster lookup
        valid_bases = set("ATCG")

        # Second pass: add k-mers
        processed = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            seq_len = len(seq)
            self.genome_size += seq_len

            logger.info(f"  Processing {record.id}: {seq_len:,} bp")

            # Process each k-mer length separately (more cache-friendly)
            for k in range(min_k, max_k + 1):
                n_positions = seq_len - k + 1
                if n_positions <= 0:
                    continue

                # Use tqdm for progress on this k value
                desc = f"    {k}bp k-mers"
                iterator = range(n_positions)
                if use_tqdm:
                    iterator = tqdm(
                        iterator, desc=desc, unit=" pos", mininterval=1.0, disable=False
                    )

                batch_kmers = []
                for i in iterator:
                    kmer = seq[i : i + k]

                    # Fast validity check using set membership
                    if all(b in valid_bases for b in kmer):
                        batch_kmers.append(kmer)

                        # Batch add to reduce per-item overhead
                        if len(batch_kmers) >= 10000:
                            for km in batch_kmers:
                                self.bloom.add(km)
                            self.kmer_count += len(batch_kmers)
                            batch_kmers = []

                # Add remaining batch
                if batch_kmers:
                    for km in batch_kmers:
                        self.bloom.add(km)
                    self.kmer_count += len(batch_kmers)

                # Optional: add 1-mismatch variants (very slow, typically skip)
                if include_mismatches:
                    logger.warning("  Adding 1-mismatch variants (slow)...")
                    # Only process a sample for mismatches
                    for i in range(0, n_positions, 100):  # Every 100th position
                        kmer = seq[i : i + k]
                        if all(b in valid_bases for b in kmer):
                            for variant in self._generate_1mm_variants(kmer):
                                self.bloom.add(variant)

        logger.info(f"Bloom filter built: {self.kmer_count:,} k-mers, {self.genome_size:,} bp")

    def add_from_kmer_files(self, kmer_prefix: str, min_k: int = 6, max_k: int = 12):
        """
        Build Bloom filter from pre-computed jellyfish k-mer files.

        MUCH faster than add_genome() because:
        - Only processes unique k-mers (not every position)
        - K-mer files already validated by jellyfish
        - No FASTA parsing overhead

        Args:
            kmer_prefix: Path prefix for k-mer files (e.g., 'data/human_chr1')
                        Will look for {prefix}_{k}mer_all.txt files
            min_k: Minimum k-mer length (default 6)
            max_k: Maximum k-mer length (default 12)
        """
        logger.info(f"Building Bloom filter from k-mer files: {kmer_prefix}")
        logger.info(f"  K-mer range: {min_k}-{max_k}bp")
        self._record_range(min_k, max_k)

        try:
            from tqdm import tqdm

            use_tqdm = True
        except ImportError:
            use_tqdm = False

        for k in range(min_k, max_k + 1):
            fpath = f"{kmer_prefix}_{k}mer_all.txt"
            if not os.path.exists(fpath):
                logger.warning(f"  K-mer file not found: {fpath}")
                continue

            # Count lines for progress bar
            n_lines = sum(1 for _ in open(fpath))
            logger.info(f"  Loading {k}bp k-mers: {n_lines:,} entries")

            with open(fpath, "r") as f:
                iterator = f
                if use_tqdm:
                    iterator = tqdm(
                        f, total=n_lines, desc=f"    {k}bp", unit=" kmers", mininterval=1.0
                    )

                batch = []
                skipped = 0
                for line in iterator:
                    parts = line.split()
                    if not parts:
                        continue
                    kmer = parts[0].upper()
                    # `add_genome` skips a k-mer holding a base outside ACGT,
                    # and this route applied no check at all, so the two
                    # disagreed about what the filter contains and kmer_count
                    # counted tokens that are not k-mers. The length check is
                    # the same claim: an entry in the k-mer table for length k
                    # that is not k bases long did not come from that table's
                    # counting run.
                    if len(kmer) != k or not self._is_valid_kmer(kmer):
                        skipped += 1
                        continue
                    batch.append(kmer)

                    if len(batch) >= 50000:
                        for km in batch:
                            self.bloom.add(km)
                        self.kmer_count += len(batch)
                        batch = []

                # Add remaining batch
                if batch:
                    for km in batch:
                        self.bloom.add(km)
                    self.kmer_count += len(batch)

                if skipped:
                    logger.warning(
                        "  %s: skipped %d entries that are not %d-mers over ACGT",
                        fpath,
                        skipped,
                        k,
                    )

        logger.info(f"Bloom filter built from k-mer files: {self.kmer_count:,} unique k-mers")

    def contains(self, kmer: str) -> bool:
        """
        Check if k-mer likely exists in background genome, on either strand.

        The reverse complement is checked too, because the background is
        double-stranded: a primer whose revcomp occurs in the host genome binds
        the host just as surely as one matching the forward strand.

        This used to test the forward strand only, while `add_genome` also
        indexed only the forward strand -- so roughly half of the primers that
        bind the background were reported absent and survived filtering. It
        also disagreed with the non-Bloom path, which counts CANONICAL k-mers
        (`jellyfish -C`) and so has always been strand-symmetric.

        Checking both at query time rather than indexing both keeps filters
        built by earlier versions correct without a rebuild, and does not grow
        the filter.

        Returns:
            True if kmer or its reverse complement is in the genome (or a false
            positive), False if both are definitely absent.
        """
        if kmer in self.bloom:
            return True
        return _reverse_complement(kmer) in self.bloom

    def estimate_match_count(self, primer: str, max_mismatches: int = 1) -> int:
        """
        Estimate number of matches (approximate).

        Returns lower bound on match count.
        """
        matches = 0

        # Exact matches
        if self.contains(primer):
            matches += 1

        # Mismatch variants
        if max_mismatches >= 1:
            for variant in self._generate_1mm_variants(primer):
                if self.contains(variant):
                    matches += 1

        return matches

    def _is_valid_kmer(self, kmer: str) -> bool:
        """Check if k-mer contains only ATCG"""
        return all(base in "ATCG" for base in kmer)

    def _generate_1mm_variants(self, seq: str) -> Set[str]:
        """All 1-mismatch variants of a sequence.

        Delegates so this and `mismatch_counts` cannot come to disagree about
        what a mismatch neighbour is -- they are counting the same population
        for the same reason, one probabilistically and one exactly.
        """
        from neoswga.core.mismatch_counts import one_mismatch_variants

        return one_mismatch_variants(seq)

    def save(self, path: str):
        """Save Bloom filter to disk"""
        logger.info(f"Saving Bloom filter to {path}")
        with open(path, "wb") as f:
            pickle.dump(
                {
                    "bloom": self.bloom,
                    "kmer_count": self.kmer_count,
                    "genome_size": self.genome_size,
                    "min_k": self.min_k,
                    "max_k": self.max_k,
                },
                f,
            )

    @classmethod
    def load(cls, path: str) -> "BackgroundBloomFilter":
        """Load Bloom filter from disk"""
        logger.info(f"Loading Bloom filter from {path}")
        from neoswga.core.safe_pickle import safe_load

        data = safe_load(path, context="bloom_filter")

        instance = cls.__new__(cls)
        instance.bloom = data["bloom"]
        instance.kmer_count = data["kmer_count"]
        instance.genome_size = data["genome_size"]
        # `.get` rather than `[...]`: a filter written before the range was
        # recorded must load and report the range as unknown, not raise.
        instance.min_k = data.get("min_k")
        instance.max_k = data.get("max_k")

        return instance

    def memory_usage_mb(self) -> float:
        """Estimate memory usage in MB"""
        import sys

        return sys.getsizeof(pickle.dumps(self.bloom)) / 1e6


class SampledGenomeIndex:
    """
    Sampled suffix array for approximate counting.

    Indexes every Nth position of genome (e.g., N=100).
    Provides ~1% accuracy estimates with 1% memory.

    For human genome:
    - Full index: 170 GB
    - Sampled (1/100): 1.7 GB
    - Accuracy: ±10% count estimate
    """

    def __init__(self, sample_rate: int = 100):
        """
        Initialize sampled index.

        Args:
            sample_rate: Index every Nth position (e.g., 100 = 1% sample)
        """
        self.sample_rate = sample_rate
        self.kmers: Dict[str, int] = defaultdict(int)
        self.genome_size = 0
        self.min_k: Optional[int] = None
        self.max_k: Optional[int] = None
        # Which quantity the counts are. "sampled_positions" means every
        # sample_rate-th position was stored and `estimate_count` extrapolates;
        # "kmer_counts" means exact jellyfish counts at sample_rate 1, where
        # that extrapolation is the identity and there is no sparsity to
        # assess. Both are written to `bg_sampled.pkl`, so without this a
        # reader cannot tell which one they hold. None means unrecorded.
        self.source: Optional[str] = None

    def _record_range(self, min_k: int, max_k: int) -> None:
        """Widen the recorded range; two builds on one index give the union."""
        self.min_k = min_k if self.min_k is None else min(self.min_k, min_k)
        self.max_k = max_k if self.max_k is None else max(self.max_k, max_k)

    def add_genome(self, fasta_path: str, min_k: int = 6, max_k: int = 12):
        """
        Add genome with sampling.

        Optimized for large genomes with progress reporting.

        Args:
            fasta_path: Path to FASTA file
            min_k: Minimum k-mer length (default 6)
            max_k: Maximum k-mer length (default 12)
        """
        logger.info(f"Building sampled index (rate=1/{self.sample_rate}): {fasta_path}")
        self._record_range(min_k, max_k)
        self.source = "sampled_positions"

        from Bio import SeqIO

        try:
            from tqdm import tqdm

            use_tqdm = True
        except ImportError:
            use_tqdm = False

        valid_bases = set("ATCG")

        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            seq_len = len(seq)
            self.genome_size += seq_len

            logger.info(
                f"  Processing {record.id}: {seq_len:,} bp (sampling every {self.sample_rate}th position)"
            )

            # Sample positions for each k-mer length
            for k in range(min_k, max_k + 1):
                n_sampled = (seq_len - k + 1) // self.sample_rate + 1

                desc = f"    {k}bp sampled"
                positions = range(0, seq_len - k + 1, self.sample_rate)

                if use_tqdm:
                    positions = tqdm(
                        positions, desc=desc, total=n_sampled, unit=" pos", mininterval=1.0
                    )

                for i in positions:
                    kmer = seq[i : i + k]
                    if all(b in valid_bases for b in kmer):
                        self.kmers[kmer] += 1

        logger.info(f"Sampled index built: {len(self.kmers):,} unique k-mers")

    def add_from_kmer_files(self, kmer_prefix: str, min_k: int = 6, max_k: int = 12) -> None:
        """Load exact counts from jellyfish tables, rather than sampling.

        This index then holds a different quantity from the sampled one: exact
        counts at `sample_rate` 1, where the extrapolation in `estimate_count`
        is the identity. Both are written to `bg_sampled.pkl`, so `source`
        records which a reader holds; without it the only thing separating them
        was that this route happened to leave `genome_size` at 0.

        jellyfish counts canonical k-mers, so one of a k-mer and its reverse
        complement carries the combined double-stranded count and the other is
        absent. That is what `estimate_count` adding both is correct for.
        """
        if self.sample_rate != 1:
            raise ValueError(
                f"exact k-mer counts must be stored at sample_rate 1, not "
                f"{self.sample_rate}: estimate_count multiplies by the rate, "
                f"so any other value would scale counts that need no scaling"
            )
        self._record_range(min_k, max_k)
        self.source = "kmer_counts"

        for k in range(min_k, max_k + 1):
            fpath = f"{kmer_prefix}_{k}mer_all.txt"
            if not os.path.exists(fpath):
                logger.warning(f"  K-mer file not found: {fpath}")
                continue
            skipped = 0
            with open(fpath) as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 2:
                        continue
                    kmer = parts[0].upper()
                    if len(kmer) != k or not self._is_valid_kmer(kmer):
                        skipped += 1
                        continue
                    self.kmers[kmer] = int(parts[1])
            if skipped:
                logger.warning(
                    "  %s: skipped %d entries that are not %d-mers over ACGT",
                    fpath,
                    skipped,
                    k,
                )

        logger.info(f"Sampled index built from k-mer files: {len(self.kmers):,} k-mers")

    def estimate_count(self, kmer: str) -> int:
        """
        Estimate count (extrapolate from sample).

        Returns:
            Estimated count in full genome
        """
        # Both strands, for the same reason as BackgroundBloomFilter.contains:
        # the background is double-stranded and the non-Bloom path counts
        # canonical k-mers.
        sampled_count = self.kmers.get(kmer, 0)
        if kmer != _reverse_complement(kmer):
            sampled_count += self.kmers.get(_reverse_complement(kmer), 0)
        return sampled_count * self.sample_rate

    def _is_valid_kmer(self, kmer: str) -> bool:
        return all(base in "ATCG" for base in kmer)

    def save(self, path: str):
        """Save index to disk"""
        with open(path, "wb") as f:
            pickle.dump(
                {
                    "kmers": dict(self.kmers),
                    "sample_rate": self.sample_rate,
                    "genome_size": self.genome_size,
                    "min_k": self.min_k,
                    "max_k": self.max_k,
                    "source": self.source,
                },
                f,
            )

    @classmethod
    def load(cls, path: str) -> "SampledGenomeIndex":
        """Load index from disk"""
        from neoswga.core.safe_pickle import safe_load

        data = safe_load(path, context="bloom_filter")

        instance = cls(sample_rate=data["sample_rate"])
        instance.kmers = defaultdict(int, data["kmers"])
        instance.genome_size = data["genome_size"]
        # `.get` for the same reason as the Bloom filter's: an index written
        # before these fields existed is unknown, not malformed.
        instance.min_k = data.get("min_k")
        instance.max_k = data.get("max_k")
        instance.source = data.get("source")
        return instance


class BackgroundFilter:
    """
    Complete background filtering pipeline.

    Combines Bloom filter (fast rejection) + sampled index (count estimation).
    """

    def __init__(
        self,
        bloom_filter: Optional[BackgroundBloomFilter] = None,
        sampled_index: Optional[SampledGenomeIndex] = None,
        config: Optional[BackgroundFilterConfig] = None,
    ):
        """
        Initialize filter.

        Args:
            bloom_filter: Pre-built Bloom filter (or None to build)
            sampled_index: Pre-built sampled index (or None to build)
            config: Filter configuration
        """
        self.bloom = bloom_filter
        self.sampled_index = sampled_index
        self.config = config or BackgroundFilterConfig()

    def build_from_genome(self, fasta_path: str):
        """
        Build both Bloom filter and sampled index from genome.

        This is a one-time cost (30-60 minutes for human genome).
        """
        genome_size = self._estimate_genome_size(fasta_path)

        # Capacity bounds DISTINCT k-mers, which is what pybloom counts and
        # allocates for. Sizing it to the base count -- once to the count
        # itself, then to ten times it -- overflowed on small genomes and
        # asked for 39.6 GB on hg38. See `distinct_kmer_capacity`.
        #
        # Mismatch variants are not indexed here for the same reason the CLI
        # path and build_background_filter skip them: each adds a further 3*k
        # entries. `count_present_neighbours` still generates them at query
        # time, so the capability is not lost, only the cost of pre-computing
        # it.
        capacity = distinct_kmer_capacity(genome_size, self.config.min_k, self.config.max_k)

        logger.info("Building Bloom filter (capacity=%s)...", f"{capacity:,}")
        self.bloom = BackgroundBloomFilter(capacity=capacity, error_rate=self.config.bloom_fp_rate)
        self.bloom.add_genome(
            fasta_path,
            include_mismatches=False,
            min_k=self.config.min_k,
            max_k=self.config.max_k,
        )

        # Build sampled index
        logger.info("Building sampled index...")
        self.sampled_index = SampledGenomeIndex(sample_rate=self.config.sample_rate)
        self.sampled_index.add_genome(
            fasta_path, min_k=self.config.min_k, max_k=self.config.max_k
        )

    def filter_primers(self, candidates: List[str]) -> List[str]:
        """
        Filter primers against background genome.

        Returns:
            List of primers passing all filters
        """
        if self.bloom is None:
            raise ValueError("Bloom filter not initialized. Call build_from_genome() first.")

        passed = []
        stats = {"total": len(candidates), "bloom_rejected": 0, "count_rejected": 0, "passed": 0}

        for primer in candidates:
            # Fast rejection via Bloom filter
            if self.bloom.contains(primer):
                # Likely exists in background (or false positive)
                # Get more accurate count estimate
                if self.sampled_index:
                    estimated_count = self.sampled_index.estimate_count(primer)
                    if estimated_count > self.config.max_exact_matches:
                        stats["count_rejected"] += 1
                        continue
                else:
                    stats["bloom_rejected"] += 1
                    continue

            # Check 1-mismatch matches
            if self.config.max_1mm_matches > 0:
                mm_matches = self._count_mismatch_matches(primer)
                if mm_matches > self.config.max_1mm_matches:
                    stats["count_rejected"] += 1
                    continue

            passed.append(primer)
            stats["passed"] += 1

        logger.info(
            f"Background filtering: {stats['total']} → {stats['passed']} "
            f"({100*stats['passed']/stats['total']:.1f}% passed)"
        )
        logger.info(f"  Bloom rejected: {stats['bloom_rejected']}")
        logger.info(f"  Count rejected: {stats['count_rejected']}")

        return passed

    def _count_mismatch_matches(self, primer: str, max_mismatches: int = 1) -> int:
        """Count approximate mismatch matches"""
        if not self.bloom:
            return 0

        return self.bloom.estimate_match_count(primer, max_mismatches)

    def _estimate_genome_size(self, fasta_path: str) -> int:
        """Quick estimate of genome size"""
        from Bio import SeqIO

        total = 0
        for record in SeqIO.parse(fasta_path, "fasta"):
            total += len(record.seq)
        return total

    def save(self, bloom_path: str, index_path: str):
        """Save both components"""
        if self.bloom:
            self.bloom.save(bloom_path)
        if self.sampled_index:
            self.sampled_index.save(index_path)

    @classmethod
    def load(
        cls, bloom_path: str, index_path: str, config: Optional[BackgroundFilterConfig] = None
    ) -> "BackgroundFilter":
        """Load pre-built filter"""
        bloom = BackgroundBloomFilter.load(bloom_path) if os.path.exists(bloom_path) else None
        index = SampledGenomeIndex.load(index_path) if os.path.exists(index_path) else None

        return cls(bloom_filter=bloom, sampled_index=index, config=config)


def build_human_genome_filter(human_fasta: str, output_dir: str):
    """
    Build and save human genome filter (one-time setup).

    Usage:
        build_human_genome_filter('human_genome.fasta', 'data/')
        # Creates: data/human_bloom.pkl, data/human_sampled.pkl
    """
    os.makedirs(output_dir, exist_ok=True)

    filter = BackgroundFilter()
    filter.build_from_genome(human_fasta)

    bloom_path = os.path.join(output_dir, "human_bloom.pkl")
    index_path = os.path.join(output_dir, "human_sampled.pkl")

    filter.save(bloom_path, index_path)

    logger.info(f"Human genome filter saved:")
    logger.info(f"  Bloom filter: {bloom_path} ({filter.bloom.memory_usage_mb():.1f} MB)")
    logger.info(f"  Sampled index: {index_path}")


def build_background_filter(
    genome_fasta: str,
    output_dir: str,
    capacity: int = None,
    error_rate: float = 0.01,
    verbose: bool = True,
    min_k: int = 6,
    max_k: int = 12,
):
    """
    Build and save background genome filter (CLI entry point).

    Creates a Bloom filter and optionally a sampled index for efficient
    background filtering of large genomes.

    Args:
        genome_fasta: Path to background genome FASTA file
        output_dir: Directory to save filter files
        capacity: Bloom filter capacity (auto-detected from genome if None)
        error_rate: Bloom filter false positive rate (default: 0.01 = 1%)
        verbose: Print progress messages
        min_k: Shortest oligo length to index (default: 6)
        max_k: Longest oligo length to index (default: 12)

    Returns:
        Tuple of (bloom_path, sampled_path) for the created filter files
    """
    os.makedirs(output_dir, exist_ok=True)

    if verbose:
        logger.info(f"Building Bloom filter for: {genome_fasta}")

    # Auto-detect capacity from genome size if not specified
    if capacity is None:
        from Bio import SeqIO

        total_size = 0
        for record in SeqIO.parse(genome_fasta, "fasta"):
            total_size += len(record.seq)
        # Bound the DISTINCT k-mers, which is what pybloom allocates for.
        capacity = distinct_kmer_capacity(total_size, min_k, max_k)
        if verbose:
            logger.info(f"Auto-detected genome size: {total_size:,} bp")
            logger.info(f"Using capacity: {capacity:,} (distinct k-mers over k={min_k}-{max_k})")

    # Build Bloom filter
    if verbose:
        logger.info(f"Building Bloom filter (capacity={capacity:,}, error_rate={error_rate})")

    bloom = BackgroundBloomFilter(capacity=capacity, error_rate=error_rate)
    # Faster without mismatch variants; they are generated at query time.
    bloom.add_genome(genome_fasta, include_mismatches=False, min_k=min_k, max_k=max_k)

    # Build sampled index for count estimation
    if verbose:
        logger.info("Building sampled index for count estimation...")
    sampled_index = SampledGenomeIndex(sample_rate=100)
    sampled_index.add_genome(genome_fasta, min_k=min_k, max_k=max_k)

    # Save filter files
    bloom_path = os.path.join(output_dir, "bg_bloom.pkl")
    sampled_path = os.path.join(output_dir, "bg_sampled.pkl")

    bloom.save(bloom_path)
    sampled_index.save(sampled_path)

    if verbose:
        logger.info(f"Filter built successfully!")
        logger.info(f"  Bloom filter: {bloom_path} ({bloom.memory_usage_mb():.1f} MB)")
        logger.info(f"  Sampled index: {sampled_path}")
        logger.info(f"  Total k-mers indexed: {bloom.kmer_count:,}")

    return bloom_path, sampled_path


if __name__ == "__main__":
    import sys

    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) > 1:
        command = sys.argv[1]

        if command == "build" and len(sys.argv) == 4:
            # Build filter from genome
            fasta_path = sys.argv[2]
            output_dir = sys.argv[3]
            build_human_genome_filter(fasta_path, output_dir)

        elif command == "test" and len(sys.argv) == 4:
            # Test filter on primers
            bloom_path = sys.argv[2]
            index_path = sys.argv[3]

            filter = BackgroundFilter.load(bloom_path, index_path)

            test_primers = ["ATCGATCG", "GCGCGCGC", "AAAAAAAA", "TTTTTTTT"]
            passed = filter.filter_primers(test_primers)

            print(f"Test primers: {test_primers}")
            print(f"Passed: {passed}")

        else:
            print("Usage:")
            print("  Build: python background_filter.py build <genome.fasta> <output_dir>")
            print("  Test:  python background_filter.py test <bloom.pkl> <index.pkl>")
    else:
        print("Usage: python background_filter.py <command> [args]")
