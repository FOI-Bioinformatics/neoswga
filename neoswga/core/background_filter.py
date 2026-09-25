"""
Efficient background genome filtering using probabilistic data structures.

For massive genomes (human 3 Gbp, tick 2.1 Gbp), exact counting is infeasible.
Use Bloom filters and sampling for fast negative selection.
"""

import logging
import os
import pickle
import re
from collections import defaultdict
from dataclasses import dataclass

import numpy as np

logger = logging.getLogger(__name__)

try:
    from pybloom_live import BloomFilter
except ImportError:
    # A library module must not write to stdout on import. Every other optional
    # dependency in this package reports through the logger, and a print here
    # reached anyone importing neoswga for an unrelated reason.
    logger.warning(
        "pybloom_live not installed, so the Bloom background path is "
        "unavailable. Install with: pip install pybloom-live"
    )
    BloomFilter = None


_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Maximal runs of unambiguous bases. A k-mer straddling anything else was
# skipped by the position-wise validity check this replaces, so scanning inside
# runs reproduces exactly the set that check kept.
_ACGT_RUN = re.compile(r"[ACGT]+")


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


# Measured with ru_maxrss, one size per process, on a dict of 12-mer string
# keys (scripts/benchmarking/sampled_index_rss.py): 124.8 B/entry at 1 million
# entries and 125.5 at 4 million. 200,000 entries reads 137.8 because fixed
# process overhead is a larger share of a small delta, the caveat
# count_coverage_rss.py records for the same method.
SAMPLED_INDEX_BYTES_PER_ENTRY = 125

# Above this the index is worth a warning. It is a resource threshold, not a
# scientific one: the point is that the structure the Bloom filter was chosen
# to avoid has reappeared beside it, at roughly a hundred times the filter's
# own size.
_LARGE_SAMPLED_INDEX_BYTES = 1e9


def projected_sampled_entries(genome_size: int, min_k: int, max_k: int, sample_rate: int) -> int:
    """How many entries a sampled index over this genome will hold.

    Each term saturates at 4**k, so above the k-mer space a longer genome adds
    nothing: the count is the same for a 144 Mb and a 3.3 Gb background at
    k 6-12.
    """
    total = 0
    for k in range(min_k, max_k + 1):
        sampled_positions = max(0, genome_size - k + 1) // max(1, sample_rate)
        total += min(4**k, sampled_positions)
    return total


def warn_if_sampled_index_is_large(
    genome_size: int, min_k: int, max_k: int, sample_rate: int
) -> None:
    """Say so when the companion index dwarfs the filter it accompanies.

    The Bloom filter exists so a host-sized background need not be held as an
    exact index. The sampled index beside it is a plain dict of k-mer strings,
    and at host scale it is by far the larger of the two: hg38 at sample rate
    100 over k 6-12 projects to 22.4 million entries and about 2.8 GB, against
    26.8 MB for the filter. A user who reached for Bloom to save memory should
    be told where the memory went.

    Extrapolated from a measured per-entry constant, not from a host-scale
    build: no Bloom filter has ever been built against a host genome here.
    """
    entries = projected_sampled_entries(genome_size, min_k, max_k, sample_rate)
    projected = entries * SAMPLED_INDEX_BYTES_PER_ENTRY
    if projected < _LARGE_SAMPLED_INDEX_BYTES:
        return

    logger.warning(
        "Sampled index over %s bp at k=%d-%d, sample rate %d, projects to "
        "%s entries and about %.1f GB, which is far larger than the Bloom "
        "filter it accompanies. Building it from pre-counted k-mer tables "
        "instead stores unique k-mers with exact counts rather than sampled "
        "positions: 'neoswga build-filter --genome <kmer_prefix> -o <dir> "
        "--from-kmers'. Figure extrapolated from %d bytes per entry measured "
        "on this structure, not from a build at this scale.",
        f"{genome_size:,}",
        min_k,
        max_k,
        sample_rate,
        f"{entries:,}",
        projected / 1e9,
        SAMPLED_INDEX_BYTES_PER_ENTRY,
    )


@dataclass
class BackgroundFilterConfig:
    """Configuration for background filtering"""

    max_exact_matches: int = 10  # Max perfect matches in background
    # A ceiling on how many of a primer's 1 + 3k mismatch neighbours the
    # background holds, NOT on a site count: 37 for a 12-mer, 91 at k=30. See
    # BackgroundBloomFilter.count_present_neighbours.
    #
    # Unset by default. The old default of 100 sat above the largest value the
    # quantity can take, so the gate could not fire at any oligo length, and no
    # measurement here supports a particular ceiling. Shipping it unset says
    # that, where a smaller number would assert a threshold nothing validates.
    max_1mm_matches: int | None = None
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

    def __init__(self, capacity: int, error_rate: float = 0.01):
        """
        Initialize Bloom filter.

        Args:
            capacity: Expected number of DISTINCT k-mers. pybloom allocates its
                bit array from this upfront, at about 9.59 bits per item at a
                1% error rate, and raises IndexError once it is exceeded, so it
                is a real allocation and a real ceiling rather than a hint. Use
                `distinct_kmer_capacity` to derive it from a genome length and
                a k range. It has no default: the former 3e9 was a 3.6 GB
                allocation nobody chose, and the one caller that took it needed
                a different number in both directions.
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
        self.min_k: int | None = None
        self.max_k: int | None = None

    def _record_range(self, min_k: int, max_k: int) -> None:
        """Widen the recorded range; two builds on one filter give the union."""
        self.min_k = min_k if self.min_k is None else min(self.min_k, min_k)
        self.max_k = max_k if self.max_k is None else max(self.max_k, max_k)

    def covers_length(self, k: int) -> bool | None:
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

        # Second pass: add k-mers
        for record in SeqIO.parse(fasta_path, "fasta"):
            seq = str(record.seq).upper()
            seq_len = len(seq)
            self.genome_size += seq_len

            logger.info(f"  Processing {record.id}: {seq_len:,} bp")

            # Locate the unambiguous stretches once per record, as spans rather
            # than substrings so nothing is copied. Every k-mer the old loop
            # kept lies wholly inside one of these, and every k-mer it skipped
            # straddles or sits outside them, so the set is unchanged.
            spans = [match.span() for match in _ACGT_RUN.finditer(seq)]

            # Process each k-mer length separately (more cache-friendly)
            for k in range(min_k, max_k + 1):
                n_positions = seq_len - k + 1
                if n_positions <= 0:
                    continue

                # The old loop re-checked all k bases at every position, an
                # O(k) test per position per k value, and collected k-mers into
                # a list only to loop over that list calling add() one at a
                # time -- pybloom has no bulk insert, so the batching bought
                # nothing. Sliding inside a known-good span makes the validity
                # test free.
                add = self.bloom.add
                added = 0
                progress = None
                if use_tqdm:
                    progress = tqdm(
                        total=n_positions,
                        desc=f"    {k}bp k-mers",
                        unit=" pos",
                        mininterval=1.0,
                    )

                for start, end in spans:
                    in_run = end - k + 1
                    if in_run <= start:
                        continue
                    for i in range(start, in_run):
                        add(seq[i : i + k])
                    added += in_run - start
                    if progress is not None:
                        progress.update(in_run - start)

                if progress is not None:
                    progress.close()
                self.kmer_count += added

                # Optional: add 1-mismatch variants (very slow, typically skip)
                if include_mismatches:
                    logger.warning("  Adding 1-mismatch variants (slow)...")
                    # Only process a sample for mismatches
                    for i in range(0, n_positions, 100):  # Every 100th position
                        kmer = seq[i : i + k]
                        if self._is_valid_kmer(kmer):
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

            with open(fpath) as f:
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

    def count_present_neighbours(self, primer: str, max_mismatches: int = 1) -> int:
        """How many of the primer and its 1-mismatch neighbours the filter holds.

        This is NOT a match count. Each neighbour contributes at most 1
        regardless of how often it occurs in the background, so the value is
        bounded by 1 + 3k -- 37 for a 12-mer, 91 at the longest oligo this tool
        supports.

        The former name, `estimate_match_count`, described it as a lower bound
        on the match count. True, but weak enough to mislead: it was compared
        against `BackgroundFilterConfig.max_1mm_matches`, whose default of 100
        no oligo length could reach, so that gate never fired whatever the
        background held.

        A Bloom filter holds presence. A count would have to come from the
        sampled index beside it, and manufacturing one is the mistake the
        sentinel removed from `get_bg_rates_via_bloom` made.
        """
        present = 1 if self.contains(primer) else 0

        if max_mismatches >= 1:
            for variant in self._generate_1mm_variants(primer):
                if self.contains(variant):
                    present += 1

        return present

    def _is_valid_kmer(self, kmer: str) -> bool:
        """Check if k-mer contains only ATCG"""
        return all(base in "ATCG" for base in kmer)

    def _generate_1mm_variants(self, seq: str) -> set[str]:
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
        """Size of the filter's bit array, in MB.

        This used to pickle the whole filter and take `sys.getsizeof` of the
        resulting bytes: a second copy of a structure sized in tens of MB for
        a host background, allocated to produce one log line, and measuring
        the SERIALISED form rather than the live allocation. `num_bits` is the
        allocation, and it is fixed at construction, so the figure does not
        move with how many k-mers have been added.
        """
        return self.bloom.num_bits / 8 / 1e6


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
        self.kmers: dict[str, int] = defaultdict(int)
        self.genome_size = 0
        self.min_k: int | None = None
        self.max_k: int | None = None
        # Which quantity the counts are. "sampled_positions" means every
        # sample_rate-th position was stored and `estimate_count` extrapolates;
        # "kmer_counts" means exact jellyfish counts at sample_rate 1, where
        # that extrapolation is the identity and there is no sparsity to
        # assess. Both are written to `bg_sampled.pkl`, so without this a
        # reader cannot tell which one they hold. None means unrecorded.
        self.source: str | None = None

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
        bloom_filter: BackgroundBloomFilter | None = None,
        sampled_index: SampledGenomeIndex | None = None,
        config: BackgroundFilterConfig | None = None,
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
        warn_if_sampled_index_is_large(
            genome_size, self.config.min_k, self.config.max_k, self.config.sample_rate
        )

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
        self.sampled_index.add_genome(fasta_path, min_k=self.config.min_k, max_k=self.config.max_k)

    def filter_primers(self, candidates: list[str]) -> list[str]:
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

            # Check 1-mismatch neighbours, when a ceiling is configured.
            if self.config.max_1mm_matches:
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
        """How many 1-mismatch neighbours the background holds; see
        BackgroundBloomFilter.count_present_neighbours for what it is not."""
        if not self.bloom:
            return 0

        return self.bloom.count_present_neighbours(primer, max_mismatches)

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
        cls, bloom_path: str, index_path: str, config: BackgroundFilterConfig | None = None
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

    bg_filter = BackgroundFilter()
    bg_filter.build_from_genome(human_fasta)

    bloom_path = os.path.join(output_dir, "human_bloom.pkl")
    index_path = os.path.join(output_dir, "human_sampled.pkl")

    bg_filter.save(bloom_path, index_path)

    logger.info("Human genome filter saved:")
    logger.info(f"  Bloom filter: {bloom_path} ({bg_filter.bloom.memory_usage_mb():.1f} MB)")
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
    warn_if_sampled_index_is_large(bloom.genome_size, min_k, max_k, 100)

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
