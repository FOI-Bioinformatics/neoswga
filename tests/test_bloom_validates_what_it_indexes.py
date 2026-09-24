"""The k-mer-file route indexed whatever the first column held.

`add_genome` skips a k-mer containing a base outside ACGT.
`add_from_kmer_files` applied no check at all: no bases, no length. So the two
routes disagreed about what a filter contains, and `kmer_count` counted tokens
that are not k-mers.

A jellyfish table should not contain such a line. But the prefix is a path the
user supplies, and `--from-kmers` is the route this work steers host-sized
backgrounds towards, so what it accepts is worth pinning.
"""

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import BackgroundBloomFilter, SampledGenomeIndex


@pytest.fixture
def table(tmp_path):
    (tmp_path / "bg_8mer_all.txt").write_text(
        "ACGTACGT 5\n"      # valid
        "ACGTNCGT 3\n"      # ambiguous base
        "ACGTAC 2\n"        # too short for this table
        "ACGTACGTAC 7\n"    # too long for this table
        "acgtacga 4\n"      # valid, lower case
        "\n"                # blank
    )
    return str(tmp_path / "bg")


def test_only_real_kmers_reach_the_filter(table):
    bloom = BackgroundBloomFilter(capacity=10_000, error_rate=0.01)
    bloom.add_from_kmer_files(table, min_k=8, max_k=8)

    assert bloom.contains("ACGTACGT")
    assert bloom.contains("ACGTACGA"), "a lower-case entry is still a k-mer"
    assert not bloom.contains("ACGTNCGT")
    assert bloom.kmer_count == 2


def test_the_sampled_index_applies_the_same_rule(table):
    index = SampledGenomeIndex(sample_rate=1)
    index.add_from_kmer_files(table, min_k=8, max_k=8)
    assert set(index.kmers) == {"ACGTACGT", "ACGTACGA"}


def test_what_was_skipped_is_reported(table, caplog):
    bloom = BackgroundBloomFilter(capacity=10_000, error_rate=0.01)
    with caplog.at_level("WARNING"):
        bloom.add_from_kmer_files(table, min_k=8, max_k=8)
    assert "skipped 3" in caplog.text


def test_exact_counts_must_not_be_scaled():
    """`estimate_count` multiplies by `sample_rate`, so storing exact counts at
    any other rate would silently inflate every background count."""
    index = SampledGenomeIndex(sample_rate=100)
    with pytest.raises(ValueError) as excinfo:
        index.add_from_kmer_files("unused", min_k=8, max_k=8)
    assert "sample_rate 1" in str(excinfo.value)
