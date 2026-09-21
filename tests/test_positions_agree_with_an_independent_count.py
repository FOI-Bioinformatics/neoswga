"""The position scan is checked against a count this file computes itself.

Task 3 of the 2026-09-21 valid-design plan. Coverage and specificity are both
computed from stored positions, so a scan that finds too few sites understates
host binding and overstates how specific a panel is. Both of this project's
2**31 defects had exactly that shape, and both sat behind a passing suite
because every test ran on a reference small enough that neither limit applied.

The oracle here is a brute-force sliding window written in this file. It does
not call the scanner, and it does not call any helper the scanner calls, so
agreement between the two is evidence rather than a tautology. The cases are
the ones where a scanner plausibly goes wrong:

- a k-mer that overlaps itself, where a naive scan advances by k and finds one
  occurrence where there are several
- a palindrome, which is its own reverse complement, so a scan that stores both
  strands separately can double it
- a reverse complement that also occurs in the forward direction
- ambiguous bases, which match nothing and must not match everything
- a circular origin, where an occurrence spans the join
- several FASTA records, where a match must not be formed across the joins

`get_all_positions_per_k` stores forward matches only; the reverse strand is
recovered by looking up the reverse complement, which is what `PositionCache`
does. The oracle follows that convention rather than inventing another one.
"""

import pytest

from neoswga.core import string_search

K = 6


def occurrences(sequence, pattern, circular=False):
    """Every start offset of `pattern` in `sequence`, overlaps included.

    A sliding window, one base at a time. `str.find` in a loop advancing by
    `len(pattern)` would miss the overlapping occurrences that make this check
    worth running at all.
    """
    if circular:
        # An occurrence may span the origin, so the search space is the
        # sequence with its own first k-1 bases appended. Offsets stay in the
        # original coordinate system.
        extended = sequence + sequence[: len(pattern) - 1]
    else:
        extended = sequence
    return [
        index
        for index in range(len(extended) - len(pattern) + 1)
        if extended[index : index + len(pattern)] == pattern and index < len(sequence)
    ]


def write_fasta(path, records):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    string_search.clear_genome_cache()
    return str(path)


def scanned(fasta, kmers, circular=False):
    found = string_search.get_all_positions_per_k(list(kmers), fasta, circular)
    return {kmer: sorted(int(p) for p in found.get(kmer, [])) for kmer in kmers}


def concatenated(records):
    """The coordinate system every stored position uses."""
    return "".join(seq for _name, seq in records)


# ---------------------------------------------------------------------------


def test_overlapping_occurrences_are_all_found(tmp_path):
    """`AAAAAA` in a run of ten A's occurs five times, not once."""
    records = [("one", "GC" + "A" * 10 + "GC")]
    fasta = write_fasta(tmp_path / "overlap.fna", records)
    sequence = concatenated(records)
    kmer = "A" * K

    assert scanned(fasta, [kmer])[kmer] == occurrences(sequence, kmer)
    assert len(occurrences(sequence, kmer)) == 5


def test_a_palindrome_is_not_counted_twice(tmp_path):
    """`GAATTC` is its own reverse complement; forward storage must hold once."""
    kmer = "GAATTC"
    assert string_search_reverse_complement(kmer) == kmer

    records = [("one", "TTAC" + kmer + "CAGG" + kmer + "AT")]
    fasta = write_fasta(tmp_path / "palindrome.fna", records)
    sequence = concatenated(records)

    assert scanned(fasta, [kmer])[kmer] == occurrences(sequence, kmer)
    assert len(scanned(fasta, [kmer])[kmer]) == 2


def test_a_kmer_and_its_reverse_complement_are_separate_answers(tmp_path):
    """Both strands are recovered by looking up two sequences, not one.

    A scanner that stored them together would make a primer's forward and
    reverse sites indistinguishable, and `PositionCache` unions them itself.
    """
    forward = "ACGGTA"
    reverse = string_search_reverse_complement(forward)
    records = [("one", "TT" + forward + "CC" + reverse + "AA" + forward + "G")]
    fasta = write_fasta(tmp_path / "rc.fna", records)
    sequence = concatenated(records)

    found = scanned(fasta, [forward, reverse])
    assert found[forward] == occurrences(sequence, forward)
    assert found[reverse] == occurrences(sequence, reverse)
    assert len(found[forward]) == 2
    assert len(found[reverse]) == 1


def test_ambiguous_bases_match_nothing_rather_than_everything(tmp_path):
    """An N in the reference must not create a site for a primer over it."""
    records = [("one", "ACGTTNNNNNGGCCAA" + "ACGGTA")]
    fasta = write_fasta(tmp_path / "ambiguous.fna", records)
    sequence = concatenated(records)

    for kmer in ("NNNNNN", "TNNNNN", "ACGGTA"):
        assert scanned(fasta, [kmer])[kmer] == occurrences(sequence, kmer), kmer
    assert scanned(fasta, ["ACGGTA"])["ACGGTA"], "a real k-mer beside the Ns is still found"


def test_a_circular_origin_is_crossed_exactly_once(tmp_path):
    """An occurrence spanning the join is one site, at its start offset."""
    tail, head = "ACG", "GTA"
    records = [("one", head + "TTTTCCCCGGGG" + tail)]
    fasta = write_fasta(tmp_path / "circular.fna", records)
    sequence = concatenated(records)
    kmer = tail + head

    assert occurrences(sequence, kmer, circular=True) == [len(sequence) - len(tail)]
    assert scanned(fasta, [kmer], circular=True)[kmer] == occurrences(sequence, kmer, circular=True)
    assert scanned(fasta, [kmer], circular=False)[kmer] == []


def test_no_match_is_formed_across_a_record_join(tmp_path):
    """Two records concatenated must not manufacture a k-mer at their seam.

    This is why the index stores record starts. A site invented at a join is
    credited to a contig it is not on, and its coverage window then extends
    into a neighbour.
    """
    records = [("one", "TTTTACG"), ("two", "GTACCCC")]
    fasta = write_fasta(tmp_path / "join.fna", records)
    seam = "ACGGTA"

    found = scanned(fasta, [seam])[seam]
    boundaries = string_search.get_cached_record_boundaries(fasta)

    # The first record starts at 0 implicitly; only interior joins are stored.
    assert boundaries == [7], boundaries
    assert found == [], (
        "a k-mer spanning the join between two records was reported as a site; " f"got {found}"
    )


def test_every_kmer_of_the_reference_is_accounted_for(tmp_path):
    """Exhaustive: every window of the reference, checked against the oracle.

    A spot check can miss a scanner that is wrong only near an edge. This asks
    about every position there is.
    """
    records = [("one", "ACGTTGCAAGGCTTACCGATGCATGGCTAACGTCAGTCCAAGTTGCACTGA")]
    fasta = write_fasta(tmp_path / "exhaustive.fna", records)
    sequence = concatenated(records)

    kmers = sorted({sequence[i : i + K] for i in range(len(sequence) - K + 1)})
    found = scanned(fasta, kmers)

    mismatched = {
        kmer: (found[kmer], occurrences(sequence, kmer))
        for kmer in kmers
        if found[kmer] != occurrences(sequence, kmer)
    }
    assert not mismatched, mismatched
    assert sum(len(v) for v in found.values()) == len(sequence) - K + 1


def test_the_total_count_matches_the_number_of_windows(tmp_path):
    """Counts and positions are two views of one fact and must agree.

    Every window of the reference is exactly one occurrence of exactly one
    k-mer, so the stored positions must sum to the window count. A scan that
    dropped sites silently would fail here without anyone naming which ones.
    """
    records = [("one", "ACGTTGCAAGGCTTACCGATG"), ("two", "CATGGCTAACGTCAGTCCAAG")]
    fasta = write_fasta(tmp_path / "counts.fna", records)
    sequence = concatenated(records)

    kmers = sorted({sequence[i : i + K] for i in range(len(sequence) - K + 1)})
    found = scanned(fasta, kmers)

    # Windows that straddle the join are not real k-mers of either record, so
    # they are excluded from the expected total rather than from the check.
    join = len(records[0][1])
    straddling = sum(1 for i in range(len(sequence) - K + 1) if i < join < i + K)
    assert sum(len(v) for v in found.values()) == len(sequence) - K + 1 - straddling


def string_search_reverse_complement(sequence):
    """Local helper, so the oracle does not borrow the scanner's own."""
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[b] for b in reversed(sequence))


def test_the_local_reverse_complement_helper_is_right():
    """The oracle's helper is itself checked, against a hand-written case."""
    assert string_search_reverse_complement("ACGGTA") == "TACCGT"
    assert string_search_reverse_complement("GAATTC") == "GAATTC"


# ---------------------------------------------------------------------------
# The two scanners must agree
# ---------------------------------------------------------------------------


def test_both_scanners_report_the_same_sites(tmp_path):
    """The Aho-Corasick path and the sliding-window path are one quantity.

    They disagreed until 2026-09-21. `get_all_positions_multi_k` rejected a
    match formed only by the join between two records; `get_all_positions_per_k`
    did not, so the same reference produced different site sets depending on
    whether pyahocorasick was installed. Two implementations of one quantity is
    how this codebase has produced disagreeing coverage numbers before.
    """
    ahocorasick = pytest.importorskip("ahocorasick")
    assert ahocorasick

    records = [("one", "ACGTTGCAAGGCTTACCGATG"), ("two", "CATGGCTAACGTCAGTCCAAG")]
    fasta = write_fasta(tmp_path / "both.fna", records)
    sequence = concatenated(records)
    kmers = sorted({sequence[i : i + K] for i in range(len(sequence) - K + 1)})

    per_k = scanned(fasta, kmers)
    multi = string_search.get_all_positions_multi_k({K: list(kmers)}, fasta, False)
    multi = {kmer: sorted(int(p) for p in multi.get(kmer, [])) for kmer in kmers}

    assert per_k == multi


def test_the_shared_join_rule_keeps_a_match_that_starts_a_record():
    """A match beginning exactly on a boundary starts a record and is real.

    Off by one here would silently delete the first k-1 sites of every contig,
    which on a draft assembly is a large and invisible loss.
    """
    from neoswga.core.string_search import spans_a_record_join

    boundaries = [7, 20]

    assert spans_a_record_join(4, 6, boundaries) is True, "4..9 contains the join at 7"
    assert spans_a_record_join(7, 6, boundaries) is False, "starts the second record"
    assert spans_a_record_join(1, 6, boundaries) is False, "1..6 ends before the join"
    assert spans_a_record_join(0, 6, []) is False, "a single record has no joins"
