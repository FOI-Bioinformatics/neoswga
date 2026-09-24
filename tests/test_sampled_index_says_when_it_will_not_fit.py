"""The sampled index must say when its projected size is large.

The Bloom filter exists so a host-sized background need not be held as an exact
index. The sampled index beside it is a plain Python dict keyed by k-mer
string, and at host scale it is by far the larger of the two, so a user who
reached for Bloom to save memory should be told where the memory went.

Measured with ru_maxrss, one size per process
(scripts/benchmarking/sampled_index_rss.py):

    200,000 entries      27.6 MB     137.8 B/entry
    1,000,000 entries   124.8 MB     124.8 B/entry
    4,000,000 entries   502.1 MB     125.5 B/entry

The smallest reads high because fixed process overhead is a larger share of a
small delta, which is the caveat `count_coverage_rss.py` already records for
the same method. Taking 125 B/entry, hg38 at sample rate 100 over k 6-12
projects to 22.4 million entries and about 2.8 GB, against 26.8 MB for the
Bloom filter it accompanies.

That is an extrapolation from a measured per-entry constant, not a measured
host-scale build. No Bloom filter has ever been built against a host genome in
this repository.
"""

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import (
    SAMPLED_INDEX_BYTES_PER_ENTRY,
    projected_sampled_entries,
    warn_if_sampled_index_is_large,
)

HG38 = 3_300_000_000


def test_the_projection_saturates_at_the_kmer_space():
    assert projected_sampled_entries(HG38, 6, 12, 100) == projected_sampled_entries(
        10 * HG38, 6, 12, 100
    )


def test_a_coarser_sample_projects_fewer_entries():
    assert projected_sampled_entries(HG38, 6, 11, 1000) < projected_sampled_entries(
        HG38, 6, 11, 100
    )


def test_the_measured_constant_is_the_one_the_benchmark_found():
    """If this is retuned, the docstring above and the script must move with it."""
    assert SAMPLED_INDEX_BYTES_PER_ENTRY == 125


def test_a_host_sized_projection_warns_and_names_the_cheaper_route(caplog):
    with caplog.at_level("WARNING"):
        warn_if_sampled_index_is_large(HG38, 6, 12, 100)
    assert "--from-kmers" in caplog.text
    assert "GB" in caplog.text


def test_a_bacterial_projection_is_silent(caplog):
    with caplog.at_level("WARNING"):
        warn_if_sampled_index_is_large(4_641_652, 6, 12, 100)
    assert caplog.text == ""


def test_building_a_small_index_warns_about_nothing(tmp_path, caplog):
    from neoswga.core.background_filter import BackgroundFilter, BackgroundFilterConfig

    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + "ACGTTGCA" * 500 + "\n")
    with caplog.at_level("WARNING"):
        BackgroundFilter(config=BackgroundFilterConfig(min_k=8, max_k=10)).build_from_genome(
            str(fasta)
        )
    assert "sampled index" not in caplog.text.lower()


def test_a_real_build_reaches_the_warning(tmp_path, caplog, monkeypatch):
    """Pins the PATH, not just the helper.

    A helper that warns correctly and is called by nothing is the defect this
    repository names Known Issue 8. The per-entry constant is inflated so a
    small genome crosses the threshold, rather than building a 3 Gb index.
    """
    import neoswga.core.background_filter as bf
    from neoswga.core.background_filter import BackgroundFilter, BackgroundFilterConfig

    monkeypatch.setattr(bf, "SAMPLED_INDEX_BYTES_PER_ENTRY", 10**9)

    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + "ACGTTGCA" * 500 + "\n")
    with caplog.at_level("WARNING"):
        BackgroundFilter(config=BackgroundFilterConfig(min_k=8, max_k=10)).build_from_genome(
            str(fasta)
        )
    assert "--from-kmers" in caplog.text
