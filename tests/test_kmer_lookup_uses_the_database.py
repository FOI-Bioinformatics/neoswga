"""Counting a known candidate list is a set operation, not a scan.

`filter` asks one question of a k-mer table: what are the counts of these
specific primers. Answering it by streaming every line into Python and testing
set membership costs 8.01 s on Drosophila at k=18 (116,702,442 lines), on top
of 9.89 s to write the 2.3 GB text table in the first place. Asking KMC to
intersect the table with a database of the candidates costs 2.5 s and produces
a result the size of the candidate list.

Measured 25 September 2026; see
docs/validation/kmer_counter_comparison_2026-09-25.md.

The tests that need a counter skip without one. The text-fallback tests do
not, and they are the ones that matter most, because every data directory that
exists today holds text tables and no database.
"""

import pytest

from neoswga.core import kmer_tables
from neoswga.core.kmer_backend import select_backend


@pytest.fixture
def counted(tmp_path):
    """A small reference counted by the configured backend."""
    backend = select_backend()
    if not backend.available():
        pytest.skip(f"the {backend.name} counter is not installed")
    fasta = tmp_path / "g.fna"
    fasta.write_text(">g\n" + "ACGTTGCAAGGCTTAC" * 60 + "\n")
    prefix = str(tmp_path / "g")
    backend.count(str(fasta), 8, prefix)
    return prefix


# --------------------------------------------------------------------------
# These need no counter: they are the path every existing directory takes.
# --------------------------------------------------------------------------


def test_a_text_only_directory_still_answers(tmp_path):
    (tmp_path / "old_8mer_all.txt").write_text("ACGTTGCA 7\nTGCAACGT 3\n")
    counts = kmer_tables.counts_for(str(tmp_path / "old"), 8, ["ACGTTGCA", "TTTTTTTT"])
    assert counts == {"ACGTTGCA": 7, "TTTTTTTT": 0}


def test_an_absent_kmer_reports_zero_rather_than_going_missing(tmp_path):
    """Zero is a measurement. A missing key reads as unknown, and an unknown
    background count PASSES the frequency gate -- the silent-zero family."""
    (tmp_path / "old_8mer_all.txt").write_text("ACGTTGCA 7\n")
    counts = kmer_tables.counts_for(str(tmp_path / "old"), 8, ["TTTTTTTT"])
    assert counts["TTTTTTTT"] == 0


def test_an_empty_request_does_no_work(tmp_path):
    assert kmer_tables.counts_for(str(tmp_path / "nothing"), 8, []) == {}


def test_no_table_at_all_is_refused_not_answered_with_zeros(tmp_path):
    """A directory with neither form has not been counted. Reporting zeros
    would say every candidate is absent from the background, which is the
    most permissive answer available."""
    with pytest.raises(Exception) as excinfo:
        kmer_tables.counts_for(str(tmp_path / "never_counted"), 8, ["ACGTTGCA"])
    assert "count-kmers" in str(excinfo.value)


# --------------------------------------------------------------------------
# These exercise the database path.
# --------------------------------------------------------------------------


def test_a_present_kmer_reports_a_positive_count(counted):
    assert kmer_tables.counts_for(counted, 8, ["ACGTTGC" + "A"])["ACGTTGCA"] > 0


def test_every_requested_kmer_appears_in_the_result(counted):
    asked = ["ACGTTGCA", "TTTTTTTT", "GGGGGGGG"]
    assert set(kmer_tables.counts_for(counted, 8, asked)) == set(asked)


def test_the_database_and_the_text_scan_agree(counted):
    """The fast path must answer what the slow path answers, or the speed is
    worthless. The counter benchmark applies the same rule to the two tools."""
    every = dict(kmer_tables._counts_from_stream(counted, 8))
    asked = sorted(every)[:20] + ["TTTTTTTT", "GGGGGGGG"]
    via_database = kmer_tables.counts_for(counted, 8, asked)
    expected = {kmer: every.get(kmer, 0) for kmer in asked}
    assert via_database == expected


def test_a_wrong_length_request_is_zero_not_an_error(counted):
    """KMC cannot hold a 5-mer in an 8-mer database. The caller asked a
    question the table cannot answer about a k-mer it cannot contain."""
    assert kmer_tables.counts_for(counted, 8, ["ACGTT"]) == {"ACGTT": 0}


def test_the_candidate_database_does_not_survive_the_call(counted, tmp_path):
    """A 500,000-primer candidate file left behind would be a real cost."""
    before = set(p.name for p in tmp_path.iterdir())
    kmer_tables.counts_for(counted, 8, ["ACGTTGCA", "TTTTTTTT"])
    after = set(p.name for p in tmp_path.iterdir())
    assert after == before
