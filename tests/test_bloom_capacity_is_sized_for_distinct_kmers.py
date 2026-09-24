"""Capacity must bound DISTINCT k-mers, not the genome's base count.

pybloom allocates its bit array upfront from `capacity`, at about 9.59 bits per
item for a 1% error rate, and its `add` increments the count only for an item
the filter did not already hold. So capacity is a bound on the DISTINCT k-mers
inserted, not on the number of insertions attempted.

`capacity = genome_size * 10` treated it as the second. The multiplier was
chosen for the seven k-mer lengths each position contributes, which is an
insertion count. Measured against the installed pybloom_live:

    genome                current capacity    allocation    distinct k-mers
    plasmid 5.4 kb                  53,860       0.06 MB             36,361
    E. coli 4.64 Mb             46,416,520      55.60 MB         10,232,681
    Drosophila 144 Mb        1,440,000,000       1.73 GB         22,368,256
    hg38 3.3 Gb             33,000,000,000      39.56 GB         22,368,256

The last two agree because each term saturates at 4**k: above the k-mer space a
longer genome cannot hold more distinct k-mers. hg38 is the documented reason
this module exists, and it is the row that cannot be allocated.
"""

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import distinct_kmer_capacity

HG38 = 3_300_000_000
# Measured: BloomFilter(capacity=1000, error_rate=0.01).num_bits == 9590
BITS_PER_ITEM = 9.59


def test_capacity_saturates_at_the_kmer_space():
    """Above 4**max_k a longer genome cannot hold more distinct k-mers."""
    assert distinct_kmer_capacity(HG38, 6, 12) == distinct_kmer_capacity(10 * HG38, 6, 12)


def test_a_host_sized_filter_fits_in_memory():
    capacity = distinct_kmer_capacity(HG38, 6, 12)
    implied_bytes = capacity * BITS_PER_ITEM / 8
    assert implied_bytes < 100e6, (
        f"a host-sized filter needs {implied_bytes / 1e9:.2f} GB; the "
        f"genome_size*10 heuristic needed 39.56 GB and could not be allocated"
    )


def test_capacity_covers_every_kmer_a_small_genome_holds():
    """The bound must not be so tight that a real build overflows it."""
    size = 5_000
    capacity = distinct_kmer_capacity(size, 6, 12)
    upper = sum(min(4**k, size) for k in range(6, 13))
    assert capacity >= upper


def test_a_narrow_range_costs_less_than_a_wide_one():
    assert distinct_kmer_capacity(HG38, 12, 12) < distinct_kmer_capacity(HG38, 6, 12)


def test_an_empty_genome_still_returns_a_usable_capacity():
    """pybloom rejects a capacity of zero, so the floor is not cosmetic."""
    assert distinct_kmer_capacity(0, 6, 12) >= 1
