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
        ">chr1\n" + "acgt" * 8 + "\n"
        ">chr2\n" + "GGCC" * 8 + "\n"
        ">chr3\n" + "TtAa" * 8 + "\n"
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
