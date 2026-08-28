"""Position scanning on genomes past the 2 Gb mark.

`get_all_positions_multi_k` hands the whole genome to pyahocorasick as one
string. Its `Automaton.iter()` indexes with a 32-bit int, so for any sequence
longer than 2**31 - 1 the scan loop never executes and it yields **nothing** --
no exception, no warning, an empty iterator. Measured directly, with the needle
planted at offset 1000:

    length 2**30      -> 1 match
    length 2**31 - 10 -> 1 match
    length 2**31 + 10 -> 0 matches

Human is 3.1 Gb, mouse 2.7 Gb. So for exactly the host this tool exists to
discriminate against, every primer got an empty position list, and the metrics
built on those positions -- `total_bg_sites`, `bg_coverage` -- reported 0.
A silent zero from a background is indistinguishable from perfect specificity,
which is the most expensive way this codebase can be wrong.

Measured on a real run: a 27-primer panel scored `total_bg_sites: 0` against
hg38 where the jellyfish count files put the true figure at 860.

The fix scans in chunks below the limit, overlapping by `k - 1` so a match
straddling a boundary is still found, and attributing each match to the chunk
its start falls in so the overlap does not double-count.
"""

import pytest

from neoswga.core import string_search


@pytest.fixture(autouse=True)
def _clear_cache():
    string_search.clear_genome_cache()
    yield
    string_search.clear_genome_cache()


@pytest.fixture
def genome(tmp_path):
    """A sequence long enough to span several test-sized chunks."""
    import random

    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(600))
    path = tmp_path / "g.fna"
    path.write_text(">g\n" + seq + "\n")
    return str(path), seq


def _brute(seq, primer):
    hits, at = [], seq.find(primer)
    while at != -1:
        hits.append(at)
        at = seq.find(primer, at + 1)
    return hits


def test_chunk_size_stays_under_the_pyahocorasick_limit():
    """The constant is the whole fix; pin it."""
    assert string_search.MAX_SCAN_CHUNK <= 2**31 - 1


def test_chunked_scan_finds_every_occurrence(genome):
    """Chunking must not lose matches -- including across chunk boundaries."""
    path, seq = genome
    primers = [seq[100:110], seq[297:307], seq[450:460]]

    found = string_search.get_all_positions_multi_k(
        {10: primers}, path, circular=False, chunk_size=64
    )

    for primer in primers:
        assert sorted(found[primer]) == _brute(seq, primer), primer


def test_chunked_scan_matches_an_unchunked_scan(genome):
    """A chunk larger than the genome is the old behaviour; they must agree."""
    path, seq = genome
    primers = [seq[100:110], seq[297:307]]

    chunked = string_search.get_all_positions_multi_k(
        {10: primers}, path, circular=False, chunk_size=64
    )
    string_search.clear_genome_cache()
    whole = string_search.get_all_positions_multi_k(
        {10: primers}, path, circular=False, chunk_size=10_000
    )

    assert {p: sorted(v) for p, v in chunked.items()} == {
        p: sorted(v) for p, v in whole.items()
    }


def test_overlap_does_not_duplicate_a_boundary_match(genome):
    """A match inside the overlap belongs to exactly one chunk."""
    path, seq = genome
    # Straddles the 64/128 boundary.
    primer = seq[60:70]

    found = string_search.get_all_positions_multi_k(
        {10: [primer]}, path, circular=False, chunk_size=64
    )

    assert sorted(found[primer]) == _brute(seq, primer)
    assert len(found[primer]) == len(set(found[primer]))


def test_circular_wrap_still_works_when_chunked(genome):
    """The wrap-around primer spans the join and must survive chunking."""
    path, seq = genome
    primer = seq[-5:] + seq[:5]

    found = string_search.get_all_positions_multi_k(
        {10: [primer]}, path, circular=True, chunk_size=64
    )

    assert found[primer] == [len(seq) - 5]
