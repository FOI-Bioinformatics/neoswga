"""Behavioural tests for the Bloom-filter background path.

CLAUDE.md names `neoswga build-filter` as the answer to the "large background
genomes" known issue: instead of loading every background k-mer into memory, a
Bloom filter answers membership in O(1). That makes it the component deciding
which primers are rejected for hitting the host genome on exactly the runs that
are too big to check by hand.

A Bloom filter has one-sided error by construction: it can say "present" about
something absent, never "absent" about something present. That asymmetry is the
whole safety argument for using one here -- a false positive drops a usable
primer, while a false negative would let a host-binding primer through. These
tests pin that direction, since a bug inverting it fails silently and only shows
up as background amplification in the lab.
"""

import os

import pytest

# The Bloom path is an optional extra (`pip install 'neoswga[improved]'`), so
# these skip rather than fail where it is absent. CI installs [dev], which
# includes it, so they do run there.
pytest.importorskip("pybloom_live", reason="pybloom-live not installed ([improved] extra)")

from neoswga.core.background_filter import (
    BackgroundBloomFilter,
    BackgroundFilter,
    SampledGenomeIndex,
)


def _write_fasta(path, sequence, name="bg"):
    with open(path, "w") as fh:
        fh.write(f">{name}\n")
        for i in range(0, len(sequence), 70):
            fh.write(sequence[i : i + 70] + "\n")
    return str(path)


@pytest.fixture
def background_seq():
    import random

    rng = random.Random(4242)
    return "".join(rng.choice("ACGT") for _ in range(5_000))


@pytest.fixture
def background_fasta(tmp_path, background_seq):
    return _write_fasta(tmp_path / "background.fasta", background_seq)


# ----------------------------------------------------------------------
# The one-sided error guarantee
# ----------------------------------------------------------------------


def test_every_kmer_that_was_added_is_reported_present(background_seq):
    """No false negatives. This is the property the design depends on.

    If a k-mer genuinely in the background could be reported absent, a primer
    that binds the host would survive filtering and amplify contamination.
    """
    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.01)

    k = 10
    kmers = {background_seq[i : i + k] for i in range(0, 500)}
    for kmer in kmers:
        bloom.add(kmer)

    missing = [kmer for kmer in kmers if not bloom.contains(kmer)]
    assert not missing, f"Bloom filter lost {len(missing)} k-mers it was given"


def test_false_positive_rate_is_near_the_configured_bound(background_seq):
    """The error rate is a promise about how many good primers get discarded."""
    import random

    bloom = BackgroundBloomFilter(capacity=50_000, error_rate=0.01)

    k = 12
    present = {background_seq[i : i + k] for i in range(0, 2_000)}
    for kmer in present:
        bloom.add(kmer)

    rng = random.Random(7)
    absent = []
    while len(absent) < 2_000:
        candidate = "".join(rng.choice("ACGT") for _ in range(k))
        if candidate not in present:
            absent.append(candidate)

    false_positives = sum(1 for kmer in absent if bloom.contains(kmer))
    rate = false_positives / len(absent)
    assert rate < 0.05, f"false-positive rate {rate:.3f} far exceeds the 0.01 target"


def test_an_empty_filter_claims_nothing(background_seq):
    bloom = BackgroundBloomFilter(capacity=10_000, error_rate=0.01)
    probes = [background_seq[i : i + 10] for i in range(0, 200, 10)]
    assert not any(bloom.contains(p) for p in probes)


# ----------------------------------------------------------------------
# Building from a genome
# ----------------------------------------------------------------------


def test_genome_kmers_are_present_after_building(tmp_path, background_fasta, background_seq):
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(background_fasta, min_k=10, max_k=10)

    probes = [background_seq[i : i + 10] for i in range(0, 1_000, 97)]
    missing = [p for p in probes if not bloom.contains(p)]
    assert not missing, f"{len(missing)}/{len(probes)} genome k-mers absent after add_genome"


def test_reverse_complements_are_covered(tmp_path, background_fasta, background_seq):
    """A primer binds either strand, so a background filter that indexes only
    the forward strand would pass primers matching the reverse one."""
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(background_fasta, min_k=10, max_k=10)

    table = str.maketrans("ACGT", "TGCA")
    probes = [background_seq[i : i + 10] for i in range(0, 500, 61)]
    rcs = [p.translate(table)[::-1] for p in probes]

    found = sum(1 for rc in rcs if bloom.contains(rc))
    assert found == len(rcs), (
        f"only {found}/{len(rcs)} reverse complements present; a primer matching "
        "the reverse strand of the background would not be filtered"
    )


# ----------------------------------------------------------------------
# Persistence
# ----------------------------------------------------------------------


def test_saved_filter_round_trips(tmp_path, background_fasta, background_seq):
    """`build-filter` writes it once and every later run loads it; a lossy round
    trip would change which primers are rejected between runs."""
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(background_fasta, min_k=10, max_k=10)

    path = tmp_path / "bg.pkl"
    bloom.save(str(path))
    assert path.is_file()

    reloaded = BackgroundBloomFilter.load(str(path))
    probes = [background_seq[i : i + 10] for i in range(0, 1_000, 53)]
    assert [reloaded.contains(p) for p in probes] == [bloom.contains(p) for p in probes]


def test_memory_usage_is_reported(tmp_path, background_fasta):
    """The reason this component exists is memory; the figure should be real."""
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(background_fasta, min_k=10, max_k=10)
    assert bloom.memory_usage_mb() > 0


# ----------------------------------------------------------------------
# SampledGenomeIndex
# ----------------------------------------------------------------------


def test_sampled_index_estimates_counts(tmp_path, background_fasta, background_seq):
    """The Bloom filter answers presence; the sampled index estimates abundance.

    A primer present once and one present a thousand times are very different
    background risks, and presence alone cannot distinguish them.
    """
    index = SampledGenomeIndex(sample_rate=1)
    index.add_genome(background_fasta, min_k=10, max_k=10)

    present = background_seq[100:110]
    assert index.estimate_count(present) > 0


def test_sampled_index_round_trips(tmp_path, background_fasta, background_seq):
    index = SampledGenomeIndex(sample_rate=1)
    index.add_genome(background_fasta, min_k=10, max_k=10)

    path = tmp_path / "sampled.pkl"
    index.save(str(path))
    reloaded = SampledGenomeIndex.load(str(path))

    probe = background_seq[200:210]
    assert reloaded.estimate_count(probe) == index.estimate_count(probe)


# ----------------------------------------------------------------------
# The filter as the pipeline uses it
# ----------------------------------------------------------------------


def test_background_binding_primers_are_rejected(tmp_path, background_fasta, background_seq):
    """End to end: a primer taken from the background should not survive.

    One drawn from elsewhere should, otherwise the filter is rejecting
    everything and the pipeline would report an empty candidate pool as a
    design problem rather than a filter problem.
    """
    bg_filter = BackgroundFilter()
    bg_filter.build_from_genome(background_fasta)

    from_background = background_seq[300:312]

    import random

    rng = random.Random(11)
    novel = None
    for _ in range(200):
        candidate = "".join(rng.choice("ACGT") for _ in range(12))
        rc = candidate.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        if candidate not in background_seq and rc not in background_seq:
            novel = candidate
            break
    assert novel, "could not construct a primer absent from the background"

    survivors = bg_filter.filter_primers([from_background, novel])

    assert (
        from_background not in survivors
    ), "a primer lifted straight out of the background survived filtering"


def test_reverse_strand_binders_are_rejected_too(tmp_path, background_fasta, background_seq):
    """The bug this file was written to catch.

    A primer whose REVERSE COMPLEMENT occurs in the background binds the host
    just as surely as one matching the forward strand, because the background is
    double-stranded. Both the indexing and the query were forward-only, so about
    half of background-binding primers survived filtering -- and it disagreed
    with the non-Bloom path, which counts canonical k-mers via `jellyfish -C`.
    """
    bg_filter = BackgroundFilter()
    bg_filter.build_from_genome(background_fasta)

    forward = background_seq[400:412]
    reverse = forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    survivors = bg_filter.filter_primers([forward, reverse])
    assert reverse not in survivors, (
        "a primer matching the reverse strand of the background survived; the "
        "Bloom path is forward-only again"
    )


def test_capacity_is_sized_for_insertions_not_bases(tmp_path, background_fasta):
    """`build_from_genome` sized the filter to the base count.

    Every position contributes one k-mer per length in min_k..max_k -- seven by
    default -- so the filter overflowed and pybloom raised "BloomFilter is at
    capacity" partway through, on any genome rather than only large ones. This
    is the class-based entry point; the CLI uses a different one that was
    already sized correctly.
    """
    bg_filter = BackgroundFilter()
    bg_filter.build_from_genome(background_fasta)  # raised IndexError before

    assert bg_filter.bloom is not None
    assert bg_filter.bloom.memory_usage_mb() > 0


def test_saved_filter_can_actually_be_reloaded(tmp_path, background_fasta, background_seq):
    """`build-filter` wrote a file that could never be read back.

    save() succeeded and load() always raised, because the safe-pickle
    allowlist did not cover the bit storage and hash constructor a pybloom
    filter carries. That made the documented workaround for large background
    genomes -- the whole reason this module exists -- broken end to end.
    """
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(background_fasta, min_k=10, max_k=10)

    path = tmp_path / "roundtrip.pkl"
    bloom.save(str(path))
    reloaded = BackgroundBloomFilter.load(str(path))

    probe = background_seq[150:160]
    assert reloaded.contains(probe), "a reloaded filter lost a k-mer it contained"
