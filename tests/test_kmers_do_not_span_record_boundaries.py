"""A k-mer that exists in no record must not be found.

Audit finding F3. `string_search` concatenates every FASTA record into one
string with no separator and scans the result, so a k-mer straddling the join
between two records is reported as a hit. The audit's probe plants
`CCCCCCGGGGGG` across the boundary of a two-record file in which neither record
contains it.

The concatenation is also the coordinate system every position is expressed in,
so it is kept. What changes is that a match spanning a record boundary is
rejected rather than recorded.

This matters for fragmented and multi-chromosome references. The shipped wMel
target is a single record, but its Drosophila background is not, so background
site counts -- the denominator of every specificity claim -- were exposed.

**Still open after this fix:** a coverage window anchored near the end of one
record still extends into the next, because windows are marked on the
concatenated array. That needs per-record lengths threaded into the coverage
helpers and is not addressed here.
"""

import pytest

from neoswga.core import string_search


@pytest.fixture
def two_records(tmp_path):
    """Neither record contains CCCCCCGGGGGG; their junction does."""
    path = tmp_path / "two.fasta"
    path.write_text(">a\nAAAAAACCCCCC\n>b\nGGGGGGTTTTTT\n")
    string_search.clear_genome_cache()
    yield str(path)
    string_search.clear_genome_cache()


def test_a_kmer_across_the_join_is_not_reported(two_records):
    spanning = "CCCCCCGGGGGG"
    present = "AAAAAACCCCCC"

    found = string_search.get_all_positions_multi_k(
        {12: [spanning, present]}, two_records, circular=False
    )

    assert (
        found[spanning] == []
    ), "a 12-mer formed only by concatenating two records was reported as a hit"
    assert found[present] == [0], "a genuine hit inside one record must survive"


def test_a_kmer_wholly_inside_the_second_record_keeps_its_offset(two_records):
    """Coordinates stay in the concatenated frame; only spanning hits go."""
    found = string_search.get_all_positions_multi_k(
        {12: ["GGGGGGTTTTTT"]}, two_records, circular=False
    )

    assert found["GGGGGGTTTTTT"] == [12]


def test_a_single_record_file_is_unaffected(tmp_path):
    """The common case must not change."""
    path = tmp_path / "one.fasta"
    path.write_text(">a\nAAAAAACCCCCCGGGGGGTTTTTT\n")
    string_search.clear_genome_cache()

    found = string_search.get_all_positions_multi_k(
        {12: ["CCCCCCGGGGGG"]}, str(path), circular=False
    )

    assert found["CCCCCCGGGGGG"] == [6]
