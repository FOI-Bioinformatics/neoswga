"""Depth belongs at the record's offset, and unobserved is not zero.

Finding F8 of the 2026-09-16 pipeline audit, the half that shows up in numbers.

`bam_gaps` matched one foreground PREFIX to one BAM reference. A prefix is one
FASTA laid out as a concatenated coordinate space, so a two-record FASTA
matched nothing by name and nothing by total length, the BAM contributed
nothing, and the function returned no gaps at all. Gap lists on multi-record
references were therefore always empty, which reads as a genome with no
coverage holes.

And a name match was accepted without checking the record length, so a BAM
carrying a truncated or differently-versioned record reported the whole record
as zero depth. Zero depth is exactly what BAM-guided expansion goes and
targets, so absent evidence became a design decision.

Both are the silent-zero shape of Known Issues 5, 6 and 13.
"""

import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from neoswga.core.bam_coverage import bam_depth_profile, bam_gaps


def _fasta(tmp_path, name, records):
    path = tmp_path / name
    with open(path, "w") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")
    return str(path)


def _bam(tmp_path, name, refs, reads):
    """`refs` is [(name, length)]; `reads` is [(ref_index, start, length)]."""
    path = tmp_path / name
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": ref_name, "LN": length} for ref_name, length in refs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for index, (ref_index, start, length) in enumerate(reads):
            segment = pysam.AlignedSegment()
            segment.query_name = f"r{index}"
            segment.query_sequence = "A" * length
            segment.flag = 0
            segment.reference_id = ref_index
            segment.reference_start = start
            segment.mapping_quality = 60
            segment.cigar = [(0, length)]
            out.write(segment)
    pysam.index(str(path))
    return str(path)


@pytest.fixture
def two_records(tmp_path):
    """A 400 bp file: `chrA` 0-199, `chrB` 200-399."""
    fasta = _fasta(tmp_path, "g.fna", [("chrA", "A" * 200), ("chrB", "C" * 200)])
    return fasta, 400


class TestDepthLandsAtTheRecordOffset:
    def test_the_second_record_is_written_after_the_first(self, tmp_path, two_records):
        """The whole point of the concatenated space. Depth on `chrB` at 0 is
        depth at 200 in the prefix's coordinates, and writing it at 0 would put
        it on top of `chrA`.
        """
        fasta, total = two_records
        bam = _bam(
            tmp_path,
            "d.bam",
            [("chrA", 200), ("chrB", 200)],
            [(1, 50, 20)],  # only chrB, at its own offset 50
        )

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert profile.depth[250] == 1, "chrB's depth did not land at 200 + 50"
        assert profile.depth[50] == 0, "chrB's depth was written over chrA"

    def test_both_records_contribute(self, tmp_path, two_records):
        fasta, total = two_records
        bam = _bam(
            tmp_path,
            "d.bam",
            [("chrA", 200), ("chrB", 200)],
            [(0, 10, 20), (1, 10, 20)],
        )

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert profile.depth[15] == 1
        assert profile.depth[215] == 1

    def test_the_array_spans_the_whole_prefix(self, tmp_path, two_records):
        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200), ("chrB", 200)], [(0, 0, 10)])

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert profile.depth.shape == (total,)


class TestUnobservedIsNotZero:
    def test_a_record_absent_from_the_bam_is_non_evaluable(self, tmp_path, two_records):
        """Not zero depth. The BAM says nothing about it, and saying nothing is
        not the same as reporting an amplification failure."""
        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200)], [(0, 0, 10)])

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert profile.evaluable[:200].all(), "chrA was observed and must be evaluable"
        assert not profile.evaluable[200:].any(), "chrB is unobserved, not zero"

    def test_a_length_mismatch_makes_the_record_non_evaluable(self, tmp_path, two_records):
        """A name match with the wrong length reported the whole record as zero
        depth. It must be refused instead."""
        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200), ("chrB", 199)], [(0, 0, 10)])

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert not profile.evaluable[200:].any()
        assert "chrB" in profile.bound.length_mismatches

    def test_both_denominators_are_reported(self, tmp_path, two_records):
        """Masking a record must not be able to inflate apparent recovery, so
        the evaluable and the total denominators are both carried."""
        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200)], [(0, 0, 10)])

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert profile.evaluable_bases == 200
        assert profile.total_bases == 400

    def test_a_bam_matching_nothing_is_evaluable_nowhere(self, tmp_path, two_records):
        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("unrelated", 400)], [(0, 0, 10)])

        profile = bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total)

        assert profile.evaluable_bases == 0
        assert not profile.evaluable.any()


class TestGapsOnAMultiRecordReference:
    def test_gaps_are_found_at_all(self, tmp_path, two_records):
        """They were always empty, because nothing ever matched."""
        fasta, total = two_records
        bam = _bam(
            tmp_path,
            "d.bam",
            [("chrA", 200), ("chrB", 200)],
            [(0, 0, 40), (1, 0, 40)],  # both records covered only at their start
        )

        gaps = bam_gaps(
            bam,
            ["g"],
            [total],
            min_depth=1,
            min_gap_size=50,
            fg_genomes=[fasta],
        )

        assert gaps, "no gaps on a two-record reference, which is the F8 symptom"

    def test_a_gap_never_straddles_a_record_boundary(self, tmp_path, two_records):
        """`chrA`'s tail and `chrB`'s head are adjacent in the concatenated
        space and on different molecules. A polymerase cannot travel between
        them, so one gap spanning both is not a distance that exists."""
        fasta, total = two_records
        bam = _bam(
            tmp_path,
            "d.bam",
            [("chrA", 200), ("chrB", 200)],
            [(0, 0, 20), (1, 180, 20)],  # start of chrA, end of chrB
        )

        gaps = bam_gaps(
            bam, ["g"], [total], min_depth=1, min_gap_size=20, fg_genomes=[fasta]
        )

        for gap in gaps:
            assert not (gap.start < 200 < gap.end), (
                f"gap {gap.start}-{gap.end} spans the boundary at 200, which is "
                "a join between two molecules"
            )

    def test_no_gap_is_reported_inside_an_unobserved_record(self, tmp_path, two_records):
        """Inventing a coverage hole from absent data is the defect, not the
        fix. `chrB` has no BAM record, so nothing may be claimed about it."""
        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200)], [(0, 0, 20)])

        gaps = bam_gaps(
            bam, ["g"], [total], min_depth=1, min_gap_size=20, fg_genomes=[fasta]
        )

        for gap in gaps:
            assert gap.start < 200, (
                f"gap {gap.start}-{gap.end} lies in chrB, which the BAM does "
                "not cover; absent evidence is not a low-depth region"
            )


class TestTheSingleRecordCaseIsUnchanged:
    def test_a_one_record_reference_still_finds_its_gap(self, tmp_path):
        """The case that worked before must keep working."""
        fasta = _fasta(tmp_path, "one.fna", [("chr1", "A" * 400)])
        bam = _bam(tmp_path, "d.bam", [("chr1", 400)], [(0, 0, 50), (0, 350, 50)])

        gaps = bam_gaps(
            bam, ["one"], [400], min_depth=1, min_gap_size=50, fg_genomes=[fasta]
        )

        assert len(gaps) == 1
        assert 40 <= gaps[0].start <= 60
        assert 340 <= gaps[0].end <= 360


class TestTheLayoutIsCrossChecked:
    def test_a_configured_length_that_disagrees_raises(self, tmp_path, two_records):
        from neoswga.core.reference_layout import LayoutMismatch

        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200)], [(0, 0, 10)])

        with pytest.raises(LayoutMismatch):
            bam_depth_profile(bam, prefix="g", fasta_path=fasta, configured_length=total + 1)

    def test_record_starts_are_checked_when_supplied(self, tmp_path, two_records):
        from neoswga.core.reference_layout import LayoutMismatch

        fasta, total = two_records
        bam = _bam(tmp_path, "d.bam", [("chrA", 200)], [(0, 0, 10)])

        with pytest.raises(LayoutMismatch):
            bam_depth_profile(
                bam,
                prefix="g",
                fasta_path=fasta,
                configured_length=total,
                record_starts=[0, 190],
            )


class TestWithoutGenomesNothingSilentlyChanges:
    def test_omitting_the_fasta_keeps_the_old_prefix_matching(self, tmp_path):
        """Three callers pass no genomes today. They must keep working, and a
        single-record reference is what they have always been able to handle.
        """
        bam = _bam(tmp_path, "d.bam", [("solo", 400)], [(0, 0, 50), (0, 350, 50)])

        gaps = bam_gaps(bam, ["solo"], [400], min_depth=1, min_gap_size=50)

        assert len(gaps) == 1

    def test_a_multi_record_reference_without_a_fasta_says_so(self, tmp_path, caplog):
        """It cannot be handled without the layout, and it used to return an
        empty list that read as "no gaps"."""
        import logging

        bam = _bam(tmp_path, "d.bam", [("chrA", 200), ("chrB", 200)], [(0, 0, 10)])

        with caplog.at_level(logging.WARNING):
            bam_gaps(bam, ["g"], [400], min_depth=1, min_gap_size=20)

        assert any("fg_genomes" in record.message for record in caplog.records), (
            "a reference that cannot be matched must name the missing input "
            "rather than returning an empty gap list"
        )
