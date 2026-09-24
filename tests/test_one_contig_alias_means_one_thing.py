"""`--contig-alias` means the same thing on both binding rules.

There are two, they key their aliases differently, and the CLI documents only
one of them:

- `match_contigs` keys on the foreground PREFIX or its basename. It is what
  `calibrate-reach` uses, and it is what the flag's help describes: "Map a
  foreground prefix/basename to a BAM contig."
- `ReferenceLayout.bind` keys on the FASTA RECORD name. Every path carrying
  `fg_genomes` uses it -- `expand-primers`, `analyze-coverage`, `iterate` --
  which makes it the production default.

So the documented form was the one that failed on the commoner path, and it
failed by binding nothing: the record was reported NOT EVALUABLE, no BAM gaps
were produced, and the run SUCCEEDED having quietly ignored the sequencing
data `--bam` exists to use. Expansion then proceeded on in-silico gaps alone.

Measured before the fix, on a single-record `mygenome.fasta` whose record is
`contig_A`, against a BAM contig `BAMNAME`:

    prefix-keyed {"mygenome": "BAMNAME"} -> 0 gaps
    record-keyed {"contig_A": "BAMNAME"} -> 1 gap

`_record_keyed_aliases` translates the prefix form where it is unambiguous,
which is a reference holding exactly one record. On a multi-record reference a
prefix names a FILE and a BAM contig names a molecule, so there is no sound
translation and it warns rather than guessing -- guessing there would bind a
whole file's depth to one contig.
"""

import pytest

pysam = pytest.importorskip("pysam")

from neoswga.core.bam_coverage import bam_gaps, match_contigs  # noqa: E402

LENGTH = 2_000


def _write_bam(path, contigs, covered_to):
    header = {
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for index, (_name, _length) in enumerate(contigs):
            for start in range(0, covered_to - 100, 50):
                read = pysam.AlignedSegment()
                read.query_name = f"r{index}_{start}"
                read.query_sequence = "A" * 100
                read.flag = 0
                read.reference_id = index
                read.reference_start = start
                read.mapping_quality = 60
                read.cigartuples = [(0, 100)]
                read.query_qualities = pysam.qualitystring_to_array("I" * 100)
                out.write(read)
    pysam.index(str(path))
    return str(path)


@pytest.fixture
def single_record(tmp_path):
    """A FASTA whose file stem and record header deliberately differ, which is
    the ordinary case: `mygenome.fasta` holding `>contig_A`."""
    fasta = tmp_path / "mygenome.fasta"
    fasta.write_text(">contig_A\n" + "ACGT" * (LENGTH // 4) + "\n")
    bam = _write_bam(tmp_path / "a.bam", [("BAMNAME", LENGTH)], covered_to=1_000)
    return str(fasta), bam


@pytest.mark.parametrize("key", ["mygenome", "contig_A"])
def test_either_spelling_binds_on_the_record_path(single_record, key):
    """The point. The prefix form is what the help documents and what a user
    writes; the record form is what `bind` documents. Both must work."""
    fasta, bam = single_record

    gaps = bam_gaps(
        bam,
        ["mygenome"],
        [LENGTH],
        min_depth=1,
        min_gap_size=10,
        fg_genomes=[fasta],
        contig_aliases={key: "BAMNAME"},
    )

    assert len(gaps) == 1, f"alias keyed on {key!r} bound nothing"


def test_the_prefix_form_still_binds_on_the_prefix_path():
    """`calibrate-reach`'s rule is unchanged, and this is what pins that the
    fix did not move it."""
    assert match_contigs(
        ["BAMNAME"], [LENGTH], ["mygenome"], [LENGTH], {"mygenome": "BAMNAME"}
    ) == {"mygenome": "BAMNAME"}


def test_a_record_keyed_alias_wins_over_a_prefix_keyed_one(tmp_path):
    """The more specific claim, and the one `bind` documents. A translation
    that overrode it would let the vaguer key decide."""
    fasta = tmp_path / "mygenome.fasta"
    fasta.write_text(">contig_A\n" + "ACGT" * (LENGTH // 4) + "\n")
    bam = _write_bam(tmp_path / "a.bam", [("RIGHT", LENGTH), ("WRONG", LENGTH)], covered_to=1_000)

    gaps = bam_gaps(
        bam,
        ["mygenome"],
        [LENGTH],
        min_depth=1,
        min_gap_size=10,
        fg_genomes=[str(fasta)],
        contig_aliases={"contig_A": "RIGHT", "mygenome": "WRONG"},
    )

    assert len(gaps) == 1


def test_a_prefix_alias_on_a_multi_record_reference_is_refused_not_guessed(tmp_path, caplog):
    """A prefix names a FILE and a BAM contig names a molecule, so there is no
    sound translation for two records. Binding the file's whole depth to one
    contig is the kind of guess this repository keeps removing."""
    fasta = tmp_path / "two.fasta"
    half = LENGTH // 2
    fasta.write_text(">recA\n" + "ACGT" * (half // 4) + "\n>recB\n" + "ACGT" * (half // 4) + "\n")
    bam = _write_bam(tmp_path / "a.bam", [("BAMNAME", half)], covered_to=400)

    with caplog.at_level("WARNING"):
        bam_gaps(
            bam,
            ["two"],
            [LENGTH],
            min_depth=1,
            min_gap_size=10,
            fg_genomes=[str(fasta)],
            contig_aliases={"two": "BAMNAME"},
        )

    text = caplog.text
    assert "--contig-alias" in text
    assert "recA" in text, "the warning must name the records that CAN be aliased"
