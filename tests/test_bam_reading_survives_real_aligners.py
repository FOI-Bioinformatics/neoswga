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


# ---------------------------------------------------------------------------
# CRAM, which cannot be read without the reference it was compressed against
# ---------------------------------------------------------------------------


@pytest.fixture
def cram(tmp_path):
    """A CRAM and the FASTA it was written against, as separate paths."""
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + "ACGT" * (LENGTH // 4) + "\n")
    pysam.faidx(str(reference))
    path = str(tmp_path / "aln.cram")
    with pysam.AlignmentFile(path, "wc", header=HEADER, reference_filename=str(reference)) as out:
        out.write(record())
    pysam.index(path)
    return path, reference


def test_a_cram_reads_when_its_reference_is_where_it_says(cram):
    """The mode stays "rb" and htslib detects the format from the magic bytes,
    so this works without a CRAM-specific branch."""
    path, _reference = cram

    assert int((compute_bam_depth(path, "chr1", LENGTH) > 0).sum()) == 20


def test_a_cram_whose_reference_moved_says_so_rather_than_truncated(cram, tmp_path, monkeypatch):
    """CRAM stores no sequence, only differences from a reference, so reading
    one without that reference is not possible. htslib says `truncated file`,
    which is a claim about the CRAM and is wrong -- the file is intact and the
    FASTA is what is missing.

    `REF_PATH` and `REF_CACHE` are cleared because htslib will otherwise fetch
    the reference by MD5 from a cache or the EBI, so on some machines this
    would pass for a reason that has nothing to do with the code.
    """
    path, reference = cram
    monkeypatch.setenv("REF_PATH", "/nonexistent")
    monkeypatch.setenv("REF_CACHE", str(tmp_path / "no-cache"))
    reference.rename(tmp_path / "moved.fa")

    with pytest.raises(Exception) as caught:
        compute_bam_depth(path, "chr1", LENGTH)

    message = str(caught.value)
    assert "reference" in message.lower(), message
    assert "--reference" in message, message
    assert "truncated" not in message.lower(), "htslib's wrong diagnosis leaked through"


def test_a_cram_reads_when_the_reference_is_supplied_explicitly(cram, tmp_path, monkeypatch):
    """The remedy the message names has to exist."""
    path, reference = cram
    monkeypatch.setenv("REF_PATH", "/nonexistent")
    monkeypatch.setenv("REF_CACHE", str(tmp_path / "no-cache"))
    moved = tmp_path / "moved.fa"
    reference.rename(moved)
    pysam.faidx(str(moved))

    depth = compute_bam_depth(path, "chr1", LENGTH, reference=str(moved))

    assert int((depth > 0).sum()) == 20


def test_a_reference_that_does_not_exist_is_refused_before_the_read(tmp_path, bam):
    """A path nobody can open is a configuration mistake, and saying so at the
    point it was given beats an htslib decode failure later."""
    path = bam([record()])

    with pytest.raises(Exception) as caught:
        compute_bam_depth(path, "chr1", LENGTH, reference=str(tmp_path / "absent.fa"))

    assert "absent.fa" in str(caught.value)


# ---------------------------------------------------------------------------
# What `count_coverage` cannot do, measured so nobody re-derives it
# ---------------------------------------------------------------------------


def test_an_n_in_a_read_contributes_no_depth(bam):
    """`count_coverage` tallies A/C/G/T, so an ambiguous base is a coverage
    hole rather than a covered base.

    This is the one place `core/depth_policy.py`'s reasoning does not carry
    through. It declines a mapping-quality floor because a gap is what
    expansion then designs primers for, and an ambiguous BASE call is the same
    situation -- the region did amplify. Closing it needs the pileup API, so it
    is pinned here rather than fixed.
    """
    path = bam([record(seq="ACGT" + "N" * 10 + "ACGT")])

    assert covered(path) == 8, "18 reference bases spanned, 10 of them N"


def test_count_secondary_cannot_count_a_record_with_no_sequence(bam):
    """bwa mem writes secondary alignments with `SEQ` set to `*`. There are no
    bases to tally, so the knob is honest for an aligner that repeats the
    sequence and inert for one that does not.

    No production path sets it -- every `DepthPolicy` built in `neoswga/` takes
    the defaults -- so this is a limit on a library knob, not on a run.
    """
    secondary = record(name="s", pos=500, flag=0x100, qualities=False)
    secondary.query_sequence = None
    secondary.cigartuples = [(MATCH, 20)]
    path = bam([record(pos=100), secondary])

    assert covered(path) == 20
    assert covered(path, DepthPolicy(count_secondary=True)) == 20, "no bases to count"


# ---------------------------------------------------------------------------
# The remedy the CRAM message names has to be reachable
# ---------------------------------------------------------------------------


def test_every_command_taking_a_bam_also_takes_a_reference():
    """The CRAM failure says "pass --reference". A message naming a flag that
    does not exist is this repository's Known Issue 8 class pointed the other
    way: not an option nobody reads, but advice nobody can follow.

    `tests/test_every_cli_option_has_an_effect.py` checks the other half, that
    the flag is read once declared.
    """
    from neoswga.cli_unified import create_parser

    actions = {}
    for group in create_parser()._subparsers._group_actions:
        for name, sub in (group.choices or {}).items():
            actions[name] = {opt for act in sub._actions for opt in act.option_strings}

    takes_bam = {name for name, flags in actions.items() if "--bam" in flags}
    assert takes_bam, "no command declares --bam; this test has gone stale"

    missing = sorted(name for name in takes_bam if "--reference" not in actions[name])
    assert not missing, f"these take --bam but cannot be told where the reference is: {missing}"


def test_open_alignment_is_the_only_door_into_an_alignment_file():
    """Every actionable message above comes from `open_alignment`. A second
    `pysam.AlignmentFile(...)` anywhere in the package gets none of them: a
    missing file, a FASTQ handed to `--bam` and an unresolvable CRAM reference
    all revert to pysam's own wording, and `--reference` is ignored.

    Two such sites existed when this audit started -- `calibrate-reach`'s
    header read and `expand-primers`' contig-name hint -- and both were reached
    only after a successful open elsewhere, which is why neither was visible in
    any failing run.

    This is a source check rather than a behaviour one on purpose: a new raw
    open is invisible to every behavioural test until someone hits the failure
    it mishandles.
    """
    import pathlib

    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    owner = package / "core" / "bam_coverage.py"

    offenders = []
    for path in package.rglob("*.py"):
        if path == owner:
            continue
        if "AlignmentFile" in path.read_text():
            offenders.append(str(path.relative_to(package.parent)))

    assert not offenders, (
        "open an alignment file through bam_coverage.open_alignment, which "
        f"names the file and the remedy and honours --reference: {offenders}"
    )


# ---------------------------------------------------------------------------
# Depth is never reported for bases the BAM cannot answer for
# ---------------------------------------------------------------------------


@pytest.fixture
def fully_covered(bam):
    """A BAM covering its whole 1,000 bp contig, so any zero in a depth array
    longer than the contig is fabricated rather than measured."""
    return bam([record(name=f"r{i}", pos=i) for i in range(0, LENGTH - 20, 10)])


def test_a_length_past_the_contig_end_is_refused_not_padded(fully_covered):
    """`count_coverage` CLAMPS its `stop` to the contig length, and the caller
    wrote the shorter result into a `length`-sized array of zeros. A zero
    there is not "no reads" -- it is no sequence to have reads on.

    Measured on a fully covered 2 kb contig asked for 5,000 bases, `bam_gaps`
    reported a gap of (1950, 5000): 3,050 bp, of which 3,000 bp is invented.
    `expand-primers` designs oligos AT gaps, so they would target a region the
    BAM says nothing whatever about, and `calibrate-reach` would fit a reach
    against the same invented zeros.
    """
    with pytest.raises(Exception) as caught:
        compute_bam_depth(fully_covered, "chr1", LENGTH * 5)

    message = str(caught.value)
    assert "chr1" in message
    assert str(LENGTH) in message.replace(",", ""), "the real contig length must be named"
    assert "--contig-alias" in message, "the remedy must be named"


def test_a_length_within_the_contig_is_measured_silently(fully_covered):
    """A contig LONGER than the configured length reads a prefix of it, and
    every base reported was observed. Refusing that would break the ordinary
    case to fix the fabricating one."""
    assert int((compute_bam_depth(fully_covered, "chr1", LENGTH // 2) > 0).sum()) == LENGTH // 2
    assert len(compute_bam_depth(fully_covered, "chr1", LENGTH)) == LENGTH


def test_the_record_path_was_already_immune(tmp_path, bam):
    """`ReferenceLayout.bind` refuses a length mismatch outright, so every path
    carrying `fg_genomes` could not reach this. Pinned so that a future
    loosening there does not quietly reopen the fabrication."""
    from neoswga.core.bam_coverage import bam_gaps

    fasta = tmp_path / "g.fasta"
    fasta.write_text(">chr1\n" + "ACGT" * (LENGTH * 5 // 4) + "\n")
    path = bam([record(name=f"r{i}", pos=i) for i in range(0, LENGTH - 20, 10)])

    gaps = bam_gaps(
        path,
        ["g"],
        [LENGTH * 5],
        min_depth=1,
        min_gap_size=100,
        fg_genomes=[str(fasta)],
        contig_aliases={"chr1": "chr1"},
    )

    assert gaps == [], "an unbound record contributes no gap, rather than a whole-record hole"


def test_a_bwa_mem_repeat_is_credited_to_one_copy_only(bam):
    """bwa mem reports a multi-mapping read at ONE locus and names the others
    in an `XA` tag on that same primary record. Nothing in this package reads
    that tag, so the other copies get no depth.

    This is the failure `core/depth_policy.py` declines a mapping-quality
    floor to prevent -- "an unmappable repeat does not become a coverage gap"
    -- reached by a route the floor cannot block. Those reads are MAPQ 0 and
    ARE counted; counting them once does not cover the other copies.

    Pinned rather than fixed: reading `XA` would credit depth to a placement
    the aligner declined to make, which is a decision and not a repair.
    """
    copies = [100, 300, 500, 700]
    reads = []
    for n in range(20):
        read = record(name=f"rep{n}", pos=copies[0], seq="A" * 50, mapq=0)
        read.set_tag("XA", "".join(f"chr1,+{p},50M,0;" for p in copies[1:]))
        reads.append(read)
    path = bam(reads)

    depth = compute_bam_depth(path, "chr1", LENGTH)

    assert int(depth[copies[0] + 25]) == 20, "the primary placement is counted"
    for other in copies[1:]:
        assert int(depth[other + 25]) == 0, f"copy at {other} got no depth from XA"
