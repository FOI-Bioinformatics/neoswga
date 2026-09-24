"""The mismatch estimate counts NEIGHBOURS PRESENT, not matches.

Each neighbour contributes at most 1 regardless of how often it occurs in the
background, so the value is bounded by 1 + 3k -- 37 for a 12-mer. The old name,
`estimate_match_count`, was read as a site count by
`BackgroundFilterConfig.max_1mm_matches`, whose default of 100 no 12-mer can
reach, so that gate could not fire whatever the background held.

A Bloom filter holds presence. A count would have to come from the sampled
index beside it, and inventing one is what the sentinel removed from
`get_bg_rates_via_bloom` did. So the repair is to name the quantity.
"""

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import BackgroundBloomFilter


@pytest.fixture
def bloom(tmp_path):
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + "ACGT" * 5000 + "\n")
    filt = BackgroundBloomFilter(capacity=100_000, error_rate=0.01)
    filt.add_genome(str(fasta), min_k=10, max_k=10)
    return filt


def test_the_estimate_is_bounded_by_the_neighbourhood(bloom):
    primer = "ACGTACGTAC"
    assert bloom.count_present_neighbours(primer) <= 1 + 3 * len(primer)


def test_a_primer_absent_with_all_its_neighbours_counts_zero(bloom):
    """Every 1-mismatch neighbour of a poly-A 10-mer is absent from an ACGT
    repeat, so the honest answer is zero rather than a floor."""
    assert bloom.count_present_neighbours("AAAAAAAAAA") == 0


def test_the_old_name_is_gone(bloom):
    assert not hasattr(bloom, "estimate_match_count"), (
        "the old name promised a match count the filter cannot supply"
    )


def test_the_configured_ceiling_is_reachable_or_unset():
    """`max_1mm_matches` must be comparable with what the method returns.

    At the old default of 100 the gate could not fire for any oligo length
    this tool supports, since 1 + 3k is 91 even at k=30. It ships unset rather
    than lowered, because no measurement here supports a particular ceiling
    and a smaller number would assert one.
    """
    from neoswga.core.background_filter import BackgroundFilterConfig

    longest_supported_k = 30
    configured = BackgroundFilterConfig().max_1mm_matches
    assert configured is None or configured <= 1 + 3 * longest_supported_k
