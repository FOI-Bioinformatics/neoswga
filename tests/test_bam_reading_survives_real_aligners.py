"""Depth counting, against the record shapes real aligners emit.

`DepthPolicy` is unit-tested in `test_the_depth_policy_is_explicit.py` against a
stub `_Read` with the four boolean attributes it reads. That pins the RULE and
cannot pin the PATH: `compute_bam_depth` hands the policy to
`pysam.AlignmentFile.count_coverage` as a `read_callback`, and whether pysam
honours it -- for which record kinds, under which of its own default filters --
is a fact about pysam rather than about the stub. Both ends existed and nothing
walked between them.

So these build real BAMs and read them back. Everything here passes today; it is
a ratchet rather than a repair, and the reason it is worth having is that a
pysam upgrade changing `count_coverage`'s filtering would move every depth
figure in this repository with no test to notice.

The shapes are the ones the three aligners this project expects actually
produce:

- **minimap2**: supplementary alignments with HARD clips, `=`/`X` CIGAR ops
  under `--eqx`, reference skips on spliced input, and records with no base
  qualities at all where the input was FASTA.
- **bwa mem**: supplementary alignments with hard clips, MAPQ up to 60.
- **bowtie2**: soft clips under `--local`, MAPQ up to 42, no supplementary
  records.

Duplicates are COUNTED and that is deliberate: `core/depth_policy.py` records
that identical start coordinates are independent priming events under
hyperbranched amplification, so marking them duplicates is an artefact of the
protocol rather than of the library.
"""

import pathlib

import numpy as np
import pytest

pysam = pytest.importorskip("pysam")

from neoswga.core.bam_coverage import compute_bam_depth  # noqa: E402
from neoswga.core.depth_policy import DepthPolicy  # noqa: E402

LENGTH = 1_000
SEQ = "ACGT" * 5  # 20 bp
HEADER = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": LENGTH}]}

# CIGAR operation codes, by name, because `(5, 10)` reads as nothing.
MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, EQUAL, DIFF = 0, 1, 2, 3, 4, 5, 7, 8


def record(name="r", pos=100, seq=SEQ, flag=0, cigar=None, mapq=60, qualities=True):
    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = seq
    read.flag = flag
    read.reference_id = 0
    read.reference_start = pos
    read.mapping_quality = mapq
    read.cigartuples = cigar or [(MATCH, len(seq))]
    if qualities:
        read.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return read


@pytest.fixture
def bam(tmp_path):
    """Build an indexed BAM from records and return its path."""

    def build(records, name="aln.bam", index=True):
        path = str(tmp_path / name)
        with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
            for read in records:
                out.write(read)
        if index:
            pysam.index(path)
        return path

    return build


def covered(path, policy=None):
    return int((compute_bam_depth(path, "chr1", LENGTH, policy=policy) > 0).sum())


# ---------------------------------------------------------------------------
# Which records count
# ---------------------------------------------------------------------------


def test_a_primary_alignment_counts(bam):
    assert covered(bam([record()])) == 20


@pytest.mark.parametrize(
    "flag,kind",
    [(0x100, "secondary"), (0x800, "supplementary"), (0x200, "QC-failed")],
)
def test_the_records_the_policy_excludes_do_not_reach_the_depth(bam, flag, kind):
    """The integration the stub test cannot do. A supplementary alignment is
    one chimeric molecule counted in several places, and minimap2 and bwa mem
    both emit them routinely; counting them inflates depth exactly where a
    long read spans a junction."""
    path = bam([record(pos=100), record(name="x", pos=500, flag=flag)])

    assert covered(path) == 20, f"{kind} reached the depth"


def test_duplicates_do_count_and_that_is_the_decision(bam):
    """Not an oversight. `core/depth_policy.py` records that identical start
    coordinates are independent priming events under hyperbranched
    amplification."""
    path = bam([record(pos=100), record(name="d", pos=500, flag=0x400)])

    assert covered(path) == 40


def test_the_policy_can_be_asked_to_count_supplementary(bam):
    """The knob exists and reaches pysam, which is the half a stub cannot
    show."""
    path = bam([record(pos=100), record(name="x", pos=500, flag=0x800)])

    assert covered(path, DepthPolicy(count_supplementary=True)) == 40


def test_a_mapq_floor_reaches_pysam(bam):
    """bowtie2 tops out at 42 and bwa at 60, so a floor written for one scale
    silently discards everything on the other."""
    path = bam([record(mapq=30)])

    assert covered(path, DepthPolicy(min_mapping_quality=20)) == 20
    assert covered(path, DepthPolicy(min_mapping_quality=42)) == 0


# ---------------------------------------------------------------------------
# CIGAR shapes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "cigar,expected,why",
    [
        ([(MATCH, 20)], 20, "plain match"),
        ([(SOFT_CLIP, 10), (MATCH, 20)], 20, "bowtie2 --local soft clip"),
        ([(HARD_CLIP, 10), (MATCH, 20)], 20, "hard clip, as on a supplementary"),
        ([(EQUAL, 10), (DIFF, 2), (EQUAL, 8)], 20, "minimap2 --eqx"),
        ([(MATCH, 10), (DEL, 5), (MATCH, 10)], 20, "deletion spans no read base"),
        ([(MATCH, 10), (REF_SKIP, 5), (MATCH, 10)], 20, "spliced reference skip"),
        ([(MATCH, 10), (INS, 5), (MATCH, 5)], 15, "insertion consumes query only"),
    ],
)
def test_every_cigar_shape_counts_the_reference_bases_it_should(bam, cigar, expected, why):
    consumed = sum(n for op, n in cigar if op in (MATCH, EQUAL, DIFF, INS, SOFT_CLIP))
    path = bam([record(seq="A" * consumed, cigar=cigar)])

    assert covered(path) == expected, why


def test_a_record_with_no_base_qualities_still_counts(bam):
    """minimap2 on a FASTA input writes `*` for QUAL. A depth path that needs
    base qualities would read the whole run as empty."""
    path = bam([record(qualities=False)])

    assert covered(path) == 20


def test_a_base_quality_floor_is_applied_when_one_is_asked_for(bam):
    read = record()
    read.query_qualities = pysam.qualitystring_to_array("!" * len(SEQ))  # Q0
    path = bam([read])

    assert covered(path) == 20
    assert covered(path, DepthPolicy(min_base_quality=20)) == 0


# ---------------------------------------------------------------------------
# What a user hits first
# ---------------------------------------------------------------------------


def test_a_bam_without_an_index_says_to_index_it(bam):
    """The first wall, and pysam's own message is `fetch called on bamfile
    without index`, which names neither the file nor the remedy."""
    path = bam([record()], index=False)

    with pytest.raises(Exception) as caught:
        compute_bam_depth(path, "chr1", LENGTH)

    message = str(caught.value)
    assert "samtools index" in message, message
    assert pathlib.Path(path).name in message, message


def test_a_missing_file_says_so_and_names_it(tmp_path):
    with pytest.raises(Exception) as caught:
        compute_bam_depth(str(tmp_path / "absent.bam"), "chr1", LENGTH)

    assert "absent.bam" in str(caught.value)


def test_a_file_that_is_not_an_alignment_says_which_file(tmp_path):
    """pysam answers `file has no sequences defined ... Consider opening with
    check_sq=False`, which is advice for a different problem."""
    path = tmp_path / "reads.fastq"
    path.write_text("@r\nACGT\n+\nIIII\n")

    with pytest.raises(Exception) as caught:
        compute_bam_depth(str(path), "chr1", LENGTH)

    message = str(caught.value)
    assert "reads.fastq" in message, message
    assert "check_sq" not in message, "pysam's own advice leaked through"


def test_a_contig_absent_from_the_bam_says_which_contigs_there_are(bam):
    """The commonest real failure: `chr1` against `1`, or an assembly whose
    names differ. `match_contigs` handles the design path; this is the raw
    one."""
    path = bam([record()])

    with pytest.raises(Exception) as caught:
        compute_bam_depth(path, "chrX", LENGTH)

    assert "chrX" in str(caught.value)


# ---------------------------------------------------------------------------
# Guard the guards
# ---------------------------------------------------------------------------


def test_an_empty_bam_is_zero_depth_rather_than_an_error(bam):
    """No reads is a measurement. It is the answer a failed amplification
    gives, and it must not be confused with a broken file."""
    path = bam([])

    assert covered(path) == 0


def test_depth_is_counted_not_merely_flagged(bam):
    """Two overlapping reads are depth 2, so the array is a count and the
    `> 0` used elsewhere in this file is not hiding a boolean."""
    path = bam([record(name="a", pos=100), record(name="b", pos=100)])

    depth = compute_bam_depth(path, "chr1", LENGTH)
    assert int(depth.max()) == 2
    assert depth.dtype == np.int32
