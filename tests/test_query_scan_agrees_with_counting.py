"""The query scan counts what a counted table counts.

The scan exists so a host genome need not be counted at all, which means its
answer is never checked against a table in the runs that use it. So it is
checked here twice: against a brute-force count written in this file, and
against `kmer_tables.counts_for` driven by a real counter.

The oracle calls nothing from `query_scan` -- not the encoder, not the record
reader -- so agreement is evidence rather than a restatement. It is the
obvious slow implementation: build every window with Python string slicing,
take the lexicographically smaller of the window and its reverse complement,
and count. An error in the 2-bit rolling encoding, the ambiguity mask or the
chunk overlap is invisible from inside that formulation and plain from here.
"""

import random
import shutil

import pytest

from neoswga.core import kmer_tables, query_scan
from neoswga.core.kmer_backend import JellyfishBackend, KmcBackend

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def oracle(records, k):
    """Canonical k-mer -> count, by brute force over each record separately."""
    counts = {}
    for record in records:
        record = record.upper()
        for i in range(len(record) - k + 1):
            window = record[i : i + k]
            if any(base not in "ACGT" for base in window):
                continue
            reverse = window.translate(_COMPLEMENT)[::-1]
            key = min(window, reverse)
            counts[key] = counts.get(key, 0) + 1
    return counts


def canonical(kmer):
    return min(kmer, kmer.translate(_COMPLEMENT)[::-1])


def write_fasta(path, records):
    path.write_text("".join(f">r{i}\n{seq}\n" for i, seq in enumerate(records)))
    return str(path)


def every_kmer(records, k):
    seen = []
    for record in records:
        for i in range(len(record.upper()) - k + 1):
            window = record.upper()[i : i + k]
            if window not in seen:
                seen.append(window)
    return seen


def check(tmp_path, records, k, queries=None, chunk=query_scan.DEFAULT_CHUNK):
    """Scan agrees with the oracle for `queries`, defaulting to every window."""
    path = write_fasta(tmp_path / "ref.fasta", records)
    if queries is None:
        queries = every_kmer(records, k)
    expected = oracle(records, k)
    found = query_scan.count_kmers(path, k, queries, chunk=chunk)
    assert set(found) == set(queries), "every requested k-mer must appear in the result"
    for kmer in queries:
        assert found[kmer] == expected.get(canonical(kmer), 0), kmer
    return found


def test_a_single_record_is_counted_window_by_window(tmp_path):
    check(tmp_path, ["ACGTACGTTGCA"], 4)


def test_overlapping_occurrences_all_count(tmp_path):
    found = check(tmp_path, ["AAAAAAA"], 3, queries=["AAA"])
    assert found["AAA"] == 5


def test_a_palindrome_is_counted_once_per_position(tmp_path):
    # ACGT is its own reverse complement, so forward and reverse agree and the
    # canonical minimum must not double count it: two occurrences, not four.
    found = check(tmp_path, ["ACGTACGT"], 4, queries=["ACGT"])
    assert found["ACGT"] == 2


def test_a_kmer_and_its_reverse_complement_are_one_quantity(tmp_path):
    found = check(tmp_path, ["AAAC" + "G" * 6 + "GTTT"], 4, queries=["AAAC", "GTTT"])
    assert found["AAAC"] == found["GTTT"] == 2


def test_no_kmer_spans_a_record_join(tmp_path):
    """The defect this package already fixed for binding positions.

    `ACGGTA` occurs in neither record and is formed only by concatenating
    them, which is exactly the fixture `string_search` uses.
    """
    found = check(tmp_path, ["TTACG", "GTACC"], 6, queries=["ACGGTA", "TACCGT"])
    assert found["ACGGTA"] == 0
    assert found["TACCGT"] == 0


def test_an_ambiguous_base_is_not_read_as_an_a(tmp_path):
    found = check(tmp_path, ["AANAA"], 3, queries=["AAA", "AAN"])
    assert found["AAA"] == 0, "a window containing N must not be counted"
    assert found["AAN"] == 0, "a query containing N occurs nowhere"


def test_lower_case_sequence_counts(tmp_path):
    found = check(tmp_path, ["acgtacgt"], 4, queries=["ACGT"])
    assert found["ACGT"] == 2


def test_a_window_spanning_a_chunk_edge_is_counted_exactly_once(tmp_path):
    """The overlap is k-1, so an off-by-one drops or doubles edge windows."""
    rng = random.Random(11)
    record = "".join(rng.choice("ACGT") for _ in range(4000))
    for chunk in (7, 8, 9, 64, 997):
        check(tmp_path, [record], 8, chunk=chunk)


def test_many_records_and_ambiguity_together(tmp_path):
    rng = random.Random(23)
    records = []
    for _ in range(12):
        seq = "".join(rng.choice("ACGTN") for _ in range(300))
        records.append(seq)
    check(tmp_path, records, 6, chunk=101)


def test_an_absent_kmer_is_reported_as_a_measured_zero(tmp_path):
    path = write_fasta(tmp_path / "ref.fasta", ["ACGTACGT"])
    found = query_scan.count_kmers(path, 4, ["TTTT", "ACGT"])
    assert found == {"TTTT": 0, "ACGT": 2}


def test_a_kmer_of_another_length_is_zero_rather_than_an_error(tmp_path):
    path = write_fasta(tmp_path / "ref.fasta", ["ACGTACGT"])
    assert query_scan.count_kmers(path, 4, ["ACG", "ACGTA"]) == {"ACG": 0, "ACGTA": 0}


def test_no_queries_is_an_empty_result(tmp_path):
    path = write_fasta(tmp_path / "ref.fasta", ["ACGT"])
    assert query_scan.count_kmers(path, 4, []) == {}


def test_a_missing_reference_raises_rather_than_counting_zero(tmp_path):
    with pytest.raises(FileNotFoundError):
        query_scan.count_kmers(str(tmp_path / "absent.fasta"), 4, ["ACGT"])


def test_a_k_too_long_to_encode_raises(tmp_path):
    path = write_fasta(tmp_path / "ref.fasta", ["ACGT"])
    with pytest.raises(ValueError, match="at most"):
        query_scan.count_kmers(path, 40, ["A" * 40])


# ---------------------------------------------------------------------------
# Against a real counter
# ---------------------------------------------------------------------------


def _a_counter_is_available():
    return JellyfishBackend().available() or KmcBackend().available()


@pytest.mark.skipif(not _a_counter_is_available(), reason="no k-mer counter installed")
def test_the_scan_agrees_with_the_counted_table(tmp_path):
    """The claim that matters: same reference, same k, same counts.

    Without this the scan could be self-consistently wrong in the same way as
    its own oracle, since both were written here. A counter is a third party.
    """
    from neoswga.core import kmer_counter

    rng = random.Random(5)
    records = ["".join(rng.choice("ACGT") for _ in range(2000)) for _ in range(3)]
    genome = write_fasta(tmp_path / "ref.fasta", records)

    prefix = str(tmp_path / "ref")
    k = 10
    kmer_counter.run_jellyfish(genome, prefix, min_k=k, max_k=k, cpus=1)
    assert kmer_tables.table_exists(prefix, k)

    queries = every_kmer(records, k)[:500]
    # Only canonical spellings: `counts_for` reports 0 for the other spelling
    # of a pair, since a canonical table stores one of the two. Every caller in
    # this package passes k-mers read from such a table, so this is the
    # comparison the pipeline actually makes.
    queries = [kmer for kmer in queries if canonical(kmer) == kmer]
    assert queries, "the fixture produced no canonical queries to compare"

    from_table = kmer_tables.counts_for(prefix, k, queries)
    from_scan = query_scan.count_kmers(genome, k, queries)
    assert from_scan == from_table
    assert sum(from_scan.values()) > 0, "the fixture counted nothing, so this proves nothing"


@pytest.mark.skipif(shutil.which("kmc") is None, reason="KMC is not installed")
def test_the_scan_agrees_with_kmc_specifically(tmp_path):
    """Both counters are canonical, and the scan must match either."""
    from neoswga.core import kmer_backend

    rng = random.Random(7)
    records = ["".join(rng.choice("ACGTN") for _ in range(3000)) for _ in range(2)]
    genome = write_fasta(tmp_path / "ref.fasta", records)
    prefix = str(tmp_path / "kmcref")
    k = 12
    KmcBackend().count(genome, k, prefix)
    assert kmer_backend.KmcBackend().database(prefix, k)

    queries = [kmer for kmer in every_kmer(records, k)[:300] if canonical(kmer) == kmer]
    assert queries
    assert query_scan.count_kmers(genome, k, queries) == kmer_tables.counts_for(prefix, k, queries)


@pytest.mark.skipif(not _a_counter_is_available(), reason="no k-mer counter installed")
def test_the_one_query_the_two_routes_answer_differently(tmp_path):
    """A non-canonical spelling: the table says 0, the scan says the real count.

    A canonical table stores one spelling of each reverse-complement pair, so
    `counts_for` cannot answer for the other one and returns 0 -- a zero that
    is not a measurement. The scan counts the pair.

    Both behaviours are asserted rather than reconciled. Every caller in this
    package reads its k-mers from a canonical table, so no query the pipeline
    makes is affected: measured on Drosophila at k=18, all 608 canonical
    queries agreed and all 13 disagreements were non-canonical spellings.
    Pinning both is what stops the gap widening unnoticed, and what would make
    a deliberate fix to `counts_for` visible as a change here.
    """
    from neoswga.core import kmer_counter

    rng = random.Random(17)
    record = "".join(rng.choice("ACGT") for _ in range(3000))
    genome = write_fasta(tmp_path / "ref.fasta", [record])
    prefix = str(tmp_path / "ref")
    k = 10
    kmer_counter.run_jellyfish(genome, prefix, min_k=k, max_k=k, cpus=1)

    spellings = [w for w in every_kmer([record], k) if canonical(w) != w]
    assert spellings, "the fixture produced no non-canonical window"
    probe = spellings[0]

    assert kmer_tables.counts_for(prefix, k, [probe])[probe] == 0
    assert query_scan.count_kmers(genome, k, [probe])[probe] >= 1
    # And the canonical spelling of the same k-mer agrees on both routes.
    mate = canonical(probe)
    assert kmer_tables.counts_for(prefix, k, [mate]) == query_scan.count_kmers(genome, k, [mate])
