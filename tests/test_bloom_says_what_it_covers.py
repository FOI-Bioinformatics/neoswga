"""An index must say what it covers, and a design outside that must refuse.

A Bloom filter built over k 6-12 answers False for every 13-mer, because no
13-mer was ever inserted. The background count is then zero, and zero clears
any frequency gate. So a design at k 13-18 screened against a phi29-range
filter passes its whole candidate pool unscreened, silently.

The filter is not at fault: absence is the honest answer to a question outside
the domain it was built over. What was missing is any record of that domain,
and any check against it.

A filter carrying no recorded range predates this field. That is UNKNOWN rather
than wrong, so it warns and proceeds, which is the rule `digest_algorithm`
established for a k-mer table written before its provenance sidecar existed.
"""

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core import filter as filter_mod
from neoswga.core.background_filter import (
    BackgroundBloomFilter,
    BackgroundFilter,
    BackgroundFilterConfig,
    SampledGenomeIndex,
)
from neoswga.core.exceptions import ReferenceDataError


@pytest.fixture
def bg_fasta(tmp_path):
    fasta = tmp_path / "bg.fasta"
    fasta.write_text(">bg\n" + "ACGTTGCA" * 500 + "\n")
    return str(fasta)


@pytest.fixture
def built(tmp_path, bg_fasta):
    """A filter and a sampled index over k 8-10, as `build-filter` writes them."""
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(bg_fasta, min_k=8, max_k=10)
    bloom.save(str(tmp_path / "bg_bloom.pkl"))

    index = SampledGenomeIndex(sample_rate=1)
    index.add_genome(bg_fasta, min_k=8, max_k=10)
    index.save(str(tmp_path / "bg_sampled.pkl"))
    return str(tmp_path / "bg_bloom.pkl")


def test_a_filter_records_the_lengths_it_indexed(bg_fasta):
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(bg_fasta, min_k=8, max_k=10)
    assert (bloom.min_k, bloom.max_k) == (8, 10)


def test_the_range_survives_a_round_trip(tmp_path, bg_fasta):
    bloom = BackgroundBloomFilter(capacity=200_000, error_rate=0.01)
    bloom.add_genome(bg_fasta, min_k=8, max_k=10)
    path = tmp_path / "rt.pkl"
    bloom.save(str(path))
    reloaded = BackgroundBloomFilter.load(str(path))
    assert (reloaded.min_k, reloaded.max_k) == (8, 10)


def test_two_builds_on_one_filter_describe_the_union(bg_fasta):
    bloom = BackgroundBloomFilter(capacity=400_000, error_rate=0.01)
    bloom.add_genome(bg_fasta, min_k=8, max_k=10)
    bloom.add_genome(bg_fasta, min_k=11, max_k=12)
    assert (bloom.min_k, bloom.max_k) == (8, 12)


def test_an_older_filter_reports_an_unknown_range(tmp_path):
    """Unknown is not stale: an artifact written before this field must load."""
    bloom = BackgroundBloomFilter(capacity=1000, error_rate=0.01)
    bloom.add("ACGTACGTAC")
    path = tmp_path / "old.pkl"
    bloom.save(str(path))

    import pickle

    with open(path, "rb") as fh:
        payload = pickle.load(fh)
    del payload["min_k"], payload["max_k"]
    with open(path, "wb") as fh:
        pickle.dump(payload, fh)

    reloaded = BackgroundBloomFilter.load(str(path))
    assert reloaded.min_k is None and reloaded.max_k is None


def test_a_primer_longer_than_the_filter_refuses(built):
    with pytest.raises(ReferenceDataError) as excinfo:
        filter_mod.get_bg_rates_via_bloom(["ACGTTGCAACGTT"], built)
    message = str(excinfo.value)
    assert "13" in message, "the refusal must name the length it cannot answer for"
    assert "build-filter" in message, "the refusal must name the command that fixes it"


def test_a_primer_shorter_than_the_filter_refuses_too(built):
    with pytest.raises(ReferenceDataError):
        filter_mod.get_bg_rates_via_bloom(["ACGTTG"], built)


def test_a_primer_inside_the_range_is_answered(built):
    counts = filter_mod.get_bg_rates_via_bloom(["ACGTTGCAA"], built)
    assert counts["ACGTTGCAA"] > 0


def test_an_unrecorded_range_warns_rather_than_refusing(tmp_path, built, caplog):
    import pickle

    with open(built, "rb") as fh:
        payload = pickle.load(fh)
    payload["min_k"] = payload["max_k"] = None
    with open(built, "wb") as fh:
        pickle.dump(payload, fh)

    with caplog.at_level("WARNING"):
        counts = filter_mod.get_bg_rates_via_bloom(["ACGTTGCAACGTT"], built)
    assert counts == {"ACGTTGCAACGTT": 0}
    assert "records no k-mer range" in caplog.text


def test_the_sampled_index_says_which_quantity_it_holds(bg_fasta):
    index = SampledGenomeIndex(sample_rate=100)
    index.add_genome(bg_fasta, min_k=8, max_k=10)
    assert index.source == "sampled_positions"


def test_the_class_entry_point_honours_the_configured_range(bg_fasta):
    """build_from_genome took add_genome's 6-12 defaults regardless of config,
    so a design at another length indexed the wrong one and could not know."""
    bg = BackgroundFilter(config=BackgroundFilterConfig(min_k=9, max_k=11))
    bg.build_from_genome(bg_fasta)
    assert (bg.bloom.min_k, bg.bloom.max_k) == (9, 11)
    assert (bg.sampled_index.min_k, bg.sampled_index.max_k) == (9, 11)
