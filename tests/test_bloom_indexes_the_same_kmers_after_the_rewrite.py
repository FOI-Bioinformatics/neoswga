"""The scan rewrite must index exactly the k-mer set the old loop did.

The old loop walked every position and validated the k-mer base by base,
skipping any that held a base outside ACGT. Splitting the record into maximal
ACGT runs and sliding within each must reproduce that set exactly, including
at run boundaries, where an off-by-one drops or invents the k-1 k-mers on
either side.

The reference below is a brute-force sliding window written in this file. It
calls neither scanner nor any helper they call, so agreement is evidence
rather than a restatement.
"""

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import BackgroundBloomFilter


def _reference_kmers(seq, k):
    """Brute force, deliberately the slow implementation."""
    found = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if all(base in "ACGT" for base in kmer):
            found.add(kmer)
    return found


def _reference_positions(seq, k):
    return sum(
        1
        for i in range(max(0, len(seq) - k + 1))
        if all(base in "ACGT" for base in seq[i : i + k])
    )


SEQUENCES = [
    "ACGTACGTACGT",
    "ACGTNACGTACGTN",
    "NNNNACGTACGTACGTNNNN",
    "ACGT" * 50 + "N" + "TGCA" * 50,
    "ACGTAC",  # exactly k
    "ACGTA",  # shorter than k
    "NNNNNN",
    "",
    "ACGTNNNNNNACGTACGTAC",
    "acgtacgtacgt",  # lower case, upper-cased on read
]


@pytest.mark.parametrize("seq", SEQUENCES)
def test_every_valid_kmer_is_present_and_the_count_matches(tmp_path, seq):
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + seq + "\n")

    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.001)
    bloom.add_genome(str(fasta), min_k=6, max_k=6)

    upper = seq.upper()
    for kmer in _reference_kmers(upper, 6):
        assert bloom.contains(kmer), f"{kmer} was indexed before and is absent now"
    assert bloom.kmer_count == _reference_positions(upper, 6)


def test_a_kmer_spanning_an_ambiguous_base_is_not_indexed(tmp_path):
    """The sharpest case: the filter must not answer for something it skipped.

    A false positive can make this pass by luck, so the error rate is set low
    and the probe is checked against the reference set rather than assumed.
    """
    seq = "TTTTTTTTTTTTNGGGGGGGGGGGG"
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + seq + "\n")

    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.0001)
    bloom.add_genome(str(fasta), min_k=6, max_k=6)

    spanning = "TTTTTN"
    assert spanning not in _reference_kmers(seq, 6)
    assert not bloom.contains(spanning)


def test_multiple_records_are_scanned_independently(tmp_path):
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">a\nACGTACGTAC\n>b\nTGCATGCATG\n")

    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.001)
    bloom.add_genome(str(fasta), min_k=6, max_k=6)

    expected = _reference_positions("ACGTACGTAC", 6) + _reference_positions("TGCATGCATG", 6)
    assert bloom.kmer_count == expected
