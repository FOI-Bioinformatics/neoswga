"""A prefix is a FASTA file in a concatenated space, not one BAM record.

Finding F8 of the 2026-09-16 pipeline audit. `bam_coverage.match_contigs` maps
one foreground PREFIX to one BAM reference, matching on name or on total
length. A prefix is one FASTA file laid out as a single concatenated
coordinate space, so a two-record FASTA matches no BAM contig by either route
and the BAM contributes nothing. A name match is also accepted without checking
the record length, so a partial match reports a whole record as zero depth.

Both are the silent-zero shape of Known Issues 5, 6 and 13: absent evidence
presented as measured absence. A genome with no observed depth reads as a
region that amplified poorly.

`core/reference_layout.py` is the coordinate contract the rest of Phase 6
depends on. It reads the record names and spans from the FASTA and refuses to
proceed when they disagree with what the position index was built over, rather
than matching on a name and hoping.
"""

import pytest

from neoswga.core.reference_layout import (
    LayoutMismatch,
    ReferenceLayout,
    read_layout,
    verify_layout,
)


def _fasta(tmp_path, name, records):
    path = tmp_path / name
    with open(path, "w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i : i + 60] + "\n")
    return str(path)


class TestReadingTheLayout:
    def test_a_single_record_is_one_span_at_zero(self, tmp_path):
        path = _fasta(tmp_path, "one.fna", [("chr1", "ACGT" * 25)])

        layout = read_layout(path, prefix="one")

        assert [r.name for r in layout.records] == ["chr1"]
        assert layout.records[0].start == 0
        assert layout.records[0].length == 100
        assert layout.total_length == 100

    def test_records_are_laid_out_end_to_end(self, tmp_path):
        """The concatenated space the position index and the caches use."""
        path = _fasta(
            tmp_path, "two.fna", [("chr1", "A" * 100), ("chr2", "C" * 50), ("chr3", "G" * 7)]
        )

        layout = read_layout(path, prefix="two")

        assert [(r.name, r.start, r.length) for r in layout.records] == [
            ("chr1", 0, 100),
            ("chr2", 100, 50),
            ("chr3", 150, 7),
        ]
        assert layout.total_length == 157

    def test_the_header_name_is_the_first_token(self, tmp_path):
        """BAM `@SQ` names carry no description, so the layout must not either,
        or every record binding would fail on a descriptive FASTA."""
        path = _fasta(tmp_path, "d.fna", [("chr1 Homo sapiens chromosome 1", "A" * 10)])

        layout = read_layout(path, prefix="d")

        assert layout.records[0].name == "chr1"

    def test_line_breaks_and_case_do_not_change_a_span(self, tmp_path):
        path = _fasta(tmp_path, "w.fna", [("r", "acgtACGTnn" * 13)])

        layout = read_layout(path, prefix="w")

        assert layout.total_length == 130

    def test_an_empty_file_has_no_records(self, tmp_path):
        path = tmp_path / "empty.fna"
        path.write_text("")

        layout = read_layout(str(path), prefix="e")

        assert layout.records == ()
        assert layout.total_length == 0

    def test_a_missing_file_is_named_in_the_error(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="absent.fna"):
            read_layout(str(tmp_path / "absent.fna"), prefix="x")


class TestLookups:
    def test_a_record_is_found_by_name(self, tmp_path):
        path = _fasta(tmp_path, "t.fna", [("a", "A" * 5), ("b", "C" * 9)])
        layout = read_layout(path, prefix="t")

        assert layout.offset_of("b") == 5
        assert layout.record("b").length == 9

    def test_an_unknown_name_raises_rather_than_returning_zero(self, tmp_path):
        """Returning 0 would write another record's depth at the first
        record's offset, which is worse than refusing."""
        path = _fasta(tmp_path, "t.fna", [("a", "A" * 5)])
        layout = read_layout(path, prefix="t")

        with pytest.raises(KeyError, match="nope"):
            layout.offset_of("nope")


class TestVerification:
    """The cross-check that makes the layout trustworthy."""

    def _layout(self, tmp_path):
        path = _fasta(tmp_path, "v.fna", [("a", "A" * 100), ("b", "C" * 50)])
        return read_layout(path, prefix="v")

    def test_matching_starts_and_length_pass(self, tmp_path):
        layout = self._layout(tmp_path)

        verify_layout(layout, record_starts=[0, 100], configured_length=150)

    def test_disagreeing_record_starts_raise(self, tmp_path):
        """A stale index against a re-fetched assembly. Every depth written
        afterwards would land at the wrong offset."""
        layout = self._layout(tmp_path)

        with pytest.raises(LayoutMismatch, match="record starts"):
            verify_layout(layout, record_starts=[0, 90], configured_length=150)

    def test_a_different_record_count_raises(self, tmp_path):
        layout = self._layout(tmp_path)

        with pytest.raises(LayoutMismatch, match="record starts"):
            verify_layout(layout, record_starts=[0], configured_length=150)

    def test_a_disagreeing_configured_length_raises(self, tmp_path):
        layout = self._layout(tmp_path)

        with pytest.raises(LayoutMismatch, match="length"):
            verify_layout(layout, record_starts=[0, 100], configured_length=151)

    def test_an_index_predating_record_starts_is_tolerated(self, tmp_path):
        """`PositionCache.get_record_starts` documents an empty list for an
        index built before they were stored. That is an absent check, not a
        failed one, so it must not raise -- but it must be visible."""
        layout = self._layout(tmp_path)

        result = verify_layout(layout, record_starts=[], configured_length=150)

        assert result.record_starts_checked is False

    def test_a_checked_layout_says_so(self, tmp_path):
        layout = self._layout(tmp_path)

        result = verify_layout(layout, record_starts=[0, 100], configured_length=150)

        assert result.record_starts_checked is True

    def test_an_absent_configured_length_is_tolerated(self, tmp_path):
        layout = self._layout(tmp_path)

        verify_layout(layout, record_starts=[0, 100], configured_length=None)


class TestTheMultiRecordCaseFindingF8Describes:
    def test_a_two_record_fasta_binds_both_records(self, tmp_path):
        """The case that matched nothing before: neither record's name nor the
        file's total length equals any single BAM contig."""
        path = _fasta(tmp_path, "g.fna", [("chrI", "A" * 1000), ("chrII", "C" * 2000)])
        layout = read_layout(path, prefix="g")

        bound = layout.bind({"chrI": 1000, "chrII": 2000})

        assert bound.matched == ("chrI", "chrII")
        assert bound.unmatched_records == ()
        assert bound.offsets == {"chrI": 0, "chrII": 1000}

    def test_a_record_whose_length_disagrees_is_refused(self, tmp_path):
        """A name match with the wrong length reported a whole record as zero
        depth. Absent evidence must not read as measured absence."""
        path = _fasta(tmp_path, "g.fna", [("chrI", "A" * 1000)])
        layout = read_layout(path, prefix="g")

        bound = layout.bind({"chrI": 999})

        assert bound.matched == ()
        assert "chrI" in bound.length_mismatches

    def test_records_absent_from_the_bam_are_reported_not_ignored(self, tmp_path):
        """They become non-evaluable, so a denominator can exclude them rather
        than counting them as zero depth."""
        path = _fasta(tmp_path, "g.fna", [("chrI", "A" * 10), ("chrII", "C" * 20)])
        layout = read_layout(path, prefix="g")

        bound = layout.bind({"chrI": 10})

        assert bound.matched == ("chrI",)
        assert bound.unmatched_records == ("chrII",)

    def test_binding_nothing_at_all_is_visible(self, tmp_path):
        path = _fasta(tmp_path, "g.fna", [("chrI", "A" * 10)])
        layout = read_layout(path, prefix="g")

        bound = layout.bind({"totally_different": 10})

        assert bound.matched == ()
        assert bound.unmatched_records == ("chrI",)

    def test_a_chr_prefix_difference_still_binds(self, tmp_path):
        """Ensembl against UCSC naming, which is the ordinary case."""
        path = _fasta(tmp_path, "g.fna", [("chr1", "A" * 10)])
        layout = read_layout(path, prefix="g")

        bound = layout.bind({"1": 10})

        assert bound.matched == ("chr1",)
        assert bound.offsets == {"chr1": 0}

    def test_an_explicit_alias_wins(self, tmp_path):
        path = _fasta(tmp_path, "g.fna", [("weird_name", "A" * 10)])
        layout = read_layout(path, prefix="g")

        bound = layout.bind({"NC_000001.11": 10}, aliases={"weird_name": "NC_000001.11"})

        assert bound.matched == ("weird_name",)
        assert bound.bam_name_for == {"weird_name": "NC_000001.11"}


class TestTheLayoutIsFrozen:
    def test_it_cannot_be_edited_under_a_running_analysis(self, tmp_path):
        path = _fasta(tmp_path, "f.fna", [("a", "A" * 4)])
        layout = read_layout(path, prefix="f")

        with pytest.raises(Exception):
            layout.total_length = 5

    def test_it_is_a_reference_layout(self, tmp_path):
        path = _fasta(tmp_path, "f.fna", [("a", "A" * 4)])

        assert isinstance(read_layout(path, prefix="f"), ReferenceLayout)
