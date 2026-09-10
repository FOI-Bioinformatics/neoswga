"""The genome cache must join whole records, never single characters.

`get_cached_genome_sequence` used to build its string with
`"".join(utility.read_fasta_file(path))`, and `read_fasta_file` yields one
character at a time. `str.join` materialises its argument into a list before
concatenating, so that list holds one pointer per base: for hg38 (3.3 Gbp)
about 3.3e9 pointers, roughly 26 GB, on a machine with 18 GB of RAM. The
process was killed with SIGKILL while "Loading genome sequence from ..." was
the last line in the log, which reads like a slow load rather than a defect.

Measured on 2026-09-02: a `neoswga filter` run against hg38 was killed with
exit 137 while 5.7 GB of RAM was free and no other process held more than
0.5 GB. After the fix a single run completes.

The genome is only ever needed whole, so the fix is to join the per-record
strings the streaming loader already produces (it uppercases them there).
Peak memory becomes the result plus the largest chromosome.

This is invisible on a small genome: joining 4.4 million single characters is
wasteful, not fatal, so every existing test passed. Only a background above a
couple of gigabases reaches it -- the same "test against a whole genome, not a
chromosome" lesson as CLAUDE.md Known Issues #5 and #7.
"""

import pytest

from neoswga.core import string_search


@pytest.fixture
def multi_record_fasta(tmp_path):
    """Three records, mixed case, so soft-masked input is covered too."""
    path = tmp_path / "genome.fna"
    path.write_text(
        ">chr1\n" + "acgt" * 8 + "\n" ">chr2\n" + "GGCC" * 8 + "\n" ">chr3\n" + "TtAa" * 8 + "\n"
    )
    return path


@pytest.fixture(autouse=True)
def _clear_cache():
    string_search._genome_cache.clear()
    yield
    string_search._genome_cache.clear()


def test_does_not_use_the_character_wise_reader(multi_record_fasta, monkeypatch):
    """`utility.read_fasta_file` yields per character and must not be joined.

    Asserted structurally rather than by measuring memory: a memory assertion
    would need a multi-gigabase fixture to separate the two cases, and the
    defect is exactly that small fixtures cannot tell them apart.
    """

    from neoswga.core import utility

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "get_cached_genome_sequence called utility.read_fasta_file, which "
            "yields one character per base; joining it allocates one pointer "
            "per base (~26 GB for hg38)"
        )

    monkeypatch.setattr(utility, "read_fasta_file", forbidden)

    sequence = string_search.get_cached_genome_sequence(str(multi_record_fasta))
    assert len(sequence) == 96


def test_joins_one_piece_per_record(multi_record_fasta, monkeypatch):
    """The piece count must scale with records, not with bases."""
    import neoswga.core.genome_io as genome_io

    seen = []
    real_streaming = genome_io.GenomeLoader.load_genome_streaming

    def spy(self, file_path):
        for record in real_streaming(self, file_path):
            seen.append(len(record))
            yield record

    monkeypatch.setattr(genome_io.GenomeLoader, "load_genome_streaming", spy)

    string_search.get_cached_genome_sequence(str(multi_record_fasta))

    assert seen == [32, 32, 32], f"expected three whole records, got {seen}"


def test_uppercases_soft_masked_bases(multi_record_fasta):
    """Soft-masked (lowercase) input must still be matched.

    The old code applied `.upper()` to the joined string, a second full-length
    copy. Dropping it is safe only because the streaming loader uppercases each
    record, which is what this pins.
    """
    sequence = string_search.get_cached_genome_sequence(str(multi_record_fasta))
    assert sequence == ("ACGT" * 8) + ("GGCC" * 8) + ("TTAA" * 8)


def test_result_is_cached_by_path(multi_record_fasta):
    first = string_search.get_cached_genome_sequence(str(multi_record_fasta))
    second = string_search.get_cached_genome_sequence(str(multi_record_fasta))
    assert first is second


def test_matches_the_character_wise_result(multi_record_fasta):
    """Equivalence with the old implementation, on a genome small enough to run it."""
    from neoswga.core import utility

    old = "".join(utility.read_fasta_file(str(multi_record_fasta))).upper()
    new = string_search.get_cached_genome_sequence(str(multi_record_fasta))
    assert new == old


def test_holds_at_most_one_record_at_a_time(tmp_path, monkeypatch):
    """`str.join` materialises its argument, so every record is alive at once.

    `str.join` calls `PySequence_Fast` on its argument, so a generator is turned
    into a list of every record before any concatenation begins. For hg38 that
    is about 3.3 GB of records alive at the moment the 3.3 GB result is
    allocated, which is the recorded 8.5 GB peak that forces hg38 filter runs to
    be serialised one at a time.

    Asserted by liveness rather than by resident memory: an RSS assertion would
    need a multi-gigabase fixture and would still be at the mercy of the
    platform's memory compression. Records are yielded as a `str` subclass,
    which unlike `str` itself can be weak-referenced, so their release is
    observable. A correct accumulator holds exactly one previously yielded
    record at any moment -- the consumer's loop variable -- so the count of live
    predecessors never exceeds one.

    Measured on a 960 MB synthetic genome in 24 records: 1927 MB peak for the
    join against 1168 MB for the loop, five deterministic replicates each in
    separate processes. Those are the figures recorded above
    `get_cached_genome_sequence`; an earlier unattributed pair, 1842 against
    850, appeared here and did not come from that run.
    """
    import weakref

    import neoswga.core.genome_io as genome_io

    class _Record(str):
        """A str that supports weak references, so release is observable."""

    refs = []
    live_at_yield = []

    def fake_streaming(self, file_path):
        for _ in range(6):
            record = _Record("ACGT" * 8)
            live_at_yield.append(sum(1 for ref in refs if ref() is not None))
            refs.append(weakref.ref(record))
            yield record

    monkeypatch.setattr(genome_io.GenomeLoader, "load_genome_streaming", fake_streaming)

    path = tmp_path / "g.fna"
    path.write_text(">g\nACGT\n")
    sequence = string_search.get_cached_genome_sequence(str(path))

    assert sequence == "ACGT" * 48
    assert max(live_at_yield) <= 1, (
        "more than one previously yielded record was still alive; the "
        "accumulator is holding every record before concatenating "
        f"(live predecessor counts {live_at_yield})"
    )


def test_a_bytearray_accumulator_is_not_used(tmp_path, monkeypatch):
    """The result must be built as `str`, not decoded from a byte buffer.

    A bytearray accumulates in place, which fixes the liveness problem above,
    but its final `.decode()` is a second full-length copy: measured peak
    1956 MB on the 960 MB fixture, worse than the join it would replace, which
    peaks at 1927 MB on the same fixture. `Automaton.iter()` also raises
    `TypeError: string required` for `bytes` and `bytearray`, so a byte buffer
    cannot be handed to the scanner to avoid that decode.
    """
    import neoswga.core.genome_io as genome_io

    def fake_streaming(self, file_path):
        yield "ACGT"
        yield "TTTT"

    monkeypatch.setattr(genome_io.GenomeLoader, "load_genome_streaming", fake_streaming)

    path = tmp_path / "g.fna"
    path.write_text(">g\nACGT\n")
    sequence = string_search.get_cached_genome_sequence(str(path))

    assert type(sequence) is str
    assert sequence == "ACGTTTTT"
