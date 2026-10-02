"""A variant file is read in one place, and a table it cannot trust blocks.

Two halves, both from rules this repository has paid for.

**One door.** `pysam.VariantFile` is opened in `core/variant_table.py` and
nowhere else, the same arrangement `bam_coverage.open_alignment` holds for
alignment files. A second raw open gets none of the refusals below: a table
made against another assembly reads as a table about this one, and the figures
that follow are plausible rather than obviously wrong.

**An artifact that cannot be parsed blocks.** A contig the reference does not
have, a REF allele that disagrees with the FASTA, an unsorted file: each is a
`ReferenceDataError` rather than a skipped row. The REF check is also what
catches a coordinate convention error, since both accepted formats are 1-based
and a 0-based table disagrees with the FASTA almost everywhere.
"""

import pathlib

import numpy as np
import pytest

from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.variant_sites import evaluate_strain_panel
from neoswga.core.variant_table import UNNAMED_STRAIN, open_variants
from tests.variant_invariants import assert_site_counts_add_up


def _panel(*args, **kwargs):
    """`evaluate_strain_panel`, with the site-count invariants asserted on the result."""
    block = evaluate_strain_panel(*args, **kwargs)
    assert_site_counts_add_up(block)
    return block


SEQUENCE = "ACGTTGCAGGTA" * 200  # 2400 bp, no run of N
SECOND = "TTACGGCATTCA" * 100  # 1200 bp


@pytest.fixture
def reference(tmp_path):
    path = tmp_path / "reference.fna"
    path.write_text(f">chr1\n{SEQUENCE}\n")
    return str(path)


@pytest.fixture
def two_record_reference(tmp_path):
    path = tmp_path / "two.fna"
    path.write_text(f">chrA\n{SEQUENCE}\n>chrB\n{SECOND}\n")
    return str(path)


def _vcf(tmp_path, rows, samples=(), name="variants.vcf"):
    lines = ["##fileformat=VCFv4.2", '##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">']
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if samples:
        header += "\tFORMAT\t" + "\t".join(samples)
    lines.append(header)
    for row in rows:
        lines.append("\t".join(str(field) for field in row))
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _tsv(tmp_path, header, rows, name="variants.tsv"):
    body = ["\t".join(header)] + ["\t".join(str(field) for field in row) for row in rows]
    path = tmp_path / name
    path.write_text("\n".join(body) + "\n")
    return str(path)


def _base(position):
    """The reference base at a 1-based position of `SEQUENCE`."""
    return SEQUENCE[position - 1]


def _other(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


# ----------------------------------------------------------------------
# Refusals
# ----------------------------------------------------------------------


def test_a_contig_the_reference_does_not_have_is_refused(tmp_path, reference):
    path = _vcf(tmp_path, [("chrZ", 10, ".", _base(10), "A", ".", "PASS", ".")])
    with pytest.raises(ReferenceDataError, match="chrZ"):
        open_variants(path, reference)


def test_a_ref_allele_that_disagrees_with_the_fasta_is_refused(tmp_path, reference):
    wrong = _other(_base(100))
    path = _vcf(tmp_path, [("chr1", 100, ".", wrong, "A", ".", "PASS", ".")])
    with pytest.raises(ReferenceDataError, match="another assembly|REF"):
        open_variants(path, reference)


def test_a_zero_based_tsv_is_refused_by_the_ref_check(tmp_path, reference):
    """The off-by-one case, which is why the REF check exists.

    `SEQUENCE` is a repeat of a 12-mer with no two adjacent bases equal, so a
    shift of one base changes the expected REF at every row.
    """
    rows = [
        ("chr1", position - 1, _base(position), _other(_base(position))) for position in (50, 90)
    ]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt"], rows)
    with pytest.raises(ReferenceDataError, match="1-based|another assembly"):
        open_variants(path, reference)


def test_the_same_tsv_written_one_based_is_accepted(tmp_path, reference):
    """The companion to the test above: the refusal is about the convention,
    not about the rows."""
    rows = [("chr1", position, _base(position), _other(_base(position))) for position in (50, 90)]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt"], rows)
    table = open_variants(path, reference)
    assert list(table.strains[0].starts) == [49, 89]


def test_an_unsorted_table_is_refused(tmp_path, reference):
    rows = [
        ("chr1", 200, ".", _base(200), "A", ".", "PASS", "."),
        ("chr1", 100, ".", _base(100), "A", ".", "PASS", "."),
    ]
    path = _vcf(tmp_path, rows)
    with pytest.raises(ReferenceDataError, match="sort"):
        open_variants(path, reference)


def test_an_interleaved_contig_is_refused(tmp_path, two_record_reference):
    rows = [
        ("chrA", 10, _base(10), _other(_base(10))),
        ("chrB", 10, SECOND[9], _other(SECOND[9])),
        ("chrA", 20, _base(20), _other(_base(20))),
    ]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt"], rows)
    with pytest.raises(ReferenceDataError, match="appears again"):
        open_variants(path, two_record_reference)


def test_a_tsv_without_a_header_is_refused(tmp_path, reference):
    path = tmp_path / "headerless.tsv"
    path.write_text(f"chr1\t50\t{_base(50)}\tA\n")
    with pytest.raises(ReferenceDataError, match="header"):
        open_variants(str(path), reference)


def test_a_genotype_that_is_neither_zero_nor_one_is_refused(tmp_path, reference):
    rows = [("chr1", 50, _base(50), _other(_base(50)), "2")]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt", "strain_a"], rows)
    with pytest.raises(ReferenceDataError, match="genotype"):
        open_variants(path, reference)


def test_a_short_row_is_refused_rather_than_read_as_missing_genotypes(tmp_path, reference):
    path = tmp_path / "short.tsv"
    path.write_text(f"chrom\tpos\tref\talt\tstrain_a\tstrain_b\nchr1\t50\t{_base(50)}\tA\t1\n")
    with pytest.raises(ReferenceDataError, match="column"):
        open_variants(str(path), reference)


def test_a_missing_file_is_refused(tmp_path, reference):
    with pytest.raises(ReferenceDataError, match="does not exist"):
        open_variants(str(tmp_path / "nothing.vcf"), reference)


def test_a_file_that_is_not_a_vcf_is_refused(tmp_path, reference):
    path = tmp_path / "nonsense.vcf"
    path.write_text("this is not a variant file at all\n")
    with pytest.raises(ReferenceDataError):
        open_variants(str(path), reference)


# ----------------------------------------------------------------------
# Unknown is not absent
# ----------------------------------------------------------------------


def test_a_strain_with_no_genotype_at_a_row_is_unavailable_with_a_reason(tmp_path, reference):
    rows = [
        ("chr1", 50, ".", _base(50), "A", ".", "PASS", ".", "GT", "1/1", "./."),
        ("chr1", 90, ".", _base(90), "C", ".", "PASS", ".", "GT", "0/0", "0/0"),
    ]
    path = _vcf(tmp_path, rows, samples=("strain_a", "strain_b"))
    table = open_variants(path, reference)

    measured = table.strain("strain_a")
    blocked = table.strain("strain_b")
    assert measured.status == "measured"
    assert measured.count == 1
    assert blocked.status == "unavailable"
    assert blocked.count is None
    assert "unknown rather than no" in blocked.unavailable


def test_a_tsv_dot_genotype_is_unknown_too(tmp_path, reference):
    rows = [(("chr1"), 50, _base(50), "A", "1", ".")]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt", "strain_a", "strain_b"], rows)
    table = open_variants(path, reference)
    assert table.strain("strain_a").status == "measured"
    assert table.strain("strain_b").status == "unavailable"


def test_a_worst_strain_figure_is_unknown_when_one_strain_is(tmp_path, reference):
    """Not the worst of the rest. A reduction over a set with one unmeasured
    member is unmeasured -- the rule `host_profile.aggregate_loads` breaks by
    returning 0.0 for an empty list."""
    from neoswga.core.position_cache import PositionCache

    rows = [
        ("chr1", 50, ".", _base(50), "A", ".", "PASS", ".", "GT", "1/1", "./."),
    ]
    path = _vcf(tmp_path, rows, samples=("strain_a", "strain_b"))
    table = open_variants(path, reference)

    primer = SEQUENCE[40:52]
    prefix = str(pathlib.Path(reference).with_suffix(""))
    cache = PositionCache(
        [prefix], [primer], genome_paths=[reference], circular=False, on_missing="scan"
    )
    block = _panel(cache, prefix, [primer], len(SEQUENCE), table, extension=300, circular=False)

    assert block["per_strain"]["strain_a"]["status"] == "measured"
    assert block["per_strain"]["strain_b"]["status"] == "unavailable"
    for key in ("worst_strain_coverage", "worst_strain_intact_fraction"):
        assert block[key]["value"] is None
        assert "unknown rather than the lowest of the rest" in block[key]["unavailable"]
    assert block["per_primer_intact_fraction_every_strain"][primer] is None


def test_a_file_with_no_samples_is_one_unnamed_strain(tmp_path, reference):
    rows = [("chr1", 50, ".", _base(50), "A", ".", "PASS", ".")]
    table = open_variants(_vcf(tmp_path, rows), reference)
    assert table.names == (UNNAMED_STRAIN,)
    assert table.strains[0].count == 1


def test_a_tsv_and_a_vcf_stating_the_same_variants_agree(tmp_path, reference):
    positions = (50, 90, 130)
    vcf = _vcf(
        tmp_path,
        [
            ("chr1", p, ".", _base(p), _other(_base(p)), ".", "PASS", ".", "GT", "1/1")
            for p in positions
        ],
        samples=("strain_a",),
    )
    tsv = _tsv(
        tmp_path,
        ["chrom", "pos", "ref", "alt", "strain_a"],
        [("chr1", p, _base(p), _other(_base(p)), 1) for p in positions],
    )
    from_vcf = open_variants(vcf, reference).strain("strain_a")
    from_tsv = open_variants(tsv, reference).strain("strain_a")
    assert np.array_equal(from_vcf.starts, from_tsv.starts)
    assert np.array_equal(from_vcf.ends, from_tsv.ends)


def test_an_indel_covers_every_base_of_its_reference_allele(tmp_path, reference):
    """And the table says so in its notes: the coordinate shift downstream is
    not modelled, so a gap length in a carrying strain is approximate."""
    deletion = SEQUENCE[49:53]
    rows = [("chr1", 50, ".", deletion, deletion[0], ".", "PASS", ".")]
    table = open_variants(_vcf(tmp_path, rows), reference)
    strain = table.strains[0]
    assert list(strain.starts) == [49]
    assert list(strain.ends) == [53]
    assert any("indel" in note or "longer than one base" in note for note in table.notes)


# ----------------------------------------------------------------------
# The door itself
# ----------------------------------------------------------------------

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "neoswga"


def test_variant_file_is_opened_in_one_module_only():
    """A source check rather than a behavioural one, for the reason the
    equivalent alignment test gives: a new raw open is invisible to every
    behavioural test until someone hits the failure it mishandles."""
    owner = PACKAGE / "core" / "variant_table.py"
    offenders = [
        str(path.relative_to(PACKAGE.parent))
        for path in PACKAGE.rglob("*.py")
        if path != owner and "VariantFile" in path.read_text()
    ]
    assert not offenders, (
        "open a variant file through variant_table.open_variants, which places "
        "every coordinate in concatenated space and refuses a table made "
        f"against another assembly: {offenders}"
    )


def test_importing_the_variant_modules_does_not_import_pysam():
    """pysam is the optional `[bam]` extra, and a TSV table needs none of it.

    Run in a subprocess: `sys.modules` in the test process is already polluted
    by every other test in the suite.
    """
    import subprocess
    import sys

    code = (
        "import neoswga.core.variant_table, neoswga.core.variant_sites, sys, json\n"
        "print(json.dumps('pysam' in sys.modules))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip().splitlines()[-1] == "false"


def test_no_affected_site_is_scored_with_the_uniform_mismatch_model():
    """The guard the plan states: an affected site is counted as lost and
    nothing else. `occupancy.mismatch_tm` is a uniform 4.0 C per mismatch whose
    record in `core/registry/model_evidence.json` has status `assumed`, and
    weighting a variant-affected site with it would be a discrimination claim
    with no measurement behind it."""
    import ast

    for name in ("variant_table.py", "variant_sites.py"):
        text = (PACKAGE / "core" / name).read_text()
        assert "mismatch_tm" not in text, name
        imported = set()
        for node in ast.walk(ast.parse(text)):
            if isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module)
                imported |= {f"{node.module}.{alias.name}" for alias in node.names}
            elif isinstance(node, ast.Import):
                imported |= {alias.name for alias in node.names}
        offenders = sorted(
            module
            for module in imported
            if any(part in module for part in ("occupancy", "mismatch", "three_prime_stability"))
        )
        assert not offenders, (
            f"{name} reaches a mismatch or occupancy model ({offenders}); the "
            f"variant route counts an affected site as lost and scores nothing"
        )


# ----------------------------------------------------------------------
# Half-calls, symbolic alleles, compression, and inconsistent rows
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "genotype,status,carries",
    [
        ("0/.", "unavailable", None),
        ("./0", "unavailable", None),
        (".|0", "unavailable", None),
        ("./1", "measured", 1),
        ("1|.", "measured", 1),
        ("0/0", "measured", 0),
        ("0", "measured", 0),
    ],
)
def test_a_half_call_is_unknown_unless_it_names_the_variant(
    tmp_path, reference, genotype, status, carries
):
    """Any allele that names an ALT makes a carrier. Otherwise any uncalled
    allele makes the genotype unknown: `0/.` does not establish the reference
    allele, and reading it as one reports a site intact in a strain nobody
    measured there. Only a fully stated reference genotype is a non-carrier."""
    rows = [("chr1", 50, ".", _base(50), _other(_base(50)), ".", "PASS", ".", "GT", genotype)]
    strain = open_variants(_vcf(tmp_path, rows, samples=("s",)), reference).strain("s")
    assert strain.status == status
    assert strain.count == carries
    if status == "unavailable":
        assert "no complete genotype" in strain.unavailable


def test_a_half_called_strain_does_not_report_its_site_intact(tmp_path, reference):
    from neoswga.core.position_cache import PositionCache

    rows = [("chr1", 50, ".", _base(50), _other(_base(50)), ".", "PASS", ".", "GT", "0/.")]
    table = open_variants(_vcf(tmp_path, rows, samples=("half",)), reference)
    primer = SEQUENCE[44:56]
    prefix = str(pathlib.Path(reference).with_suffix(""))
    cache = PositionCache(
        [prefix], [primer], genome_paths=[reference], circular=False, on_missing="scan"
    )
    block = _panel(cache, prefix, [primer], len(SEQUENCE), table, extension=300)
    assert block["per_strain"]["half"]["status"] == "unavailable"
    assert block["per_strain"]["half"]["intact_sites"]["value"] is None


def _symbolic_vcf(tmp_path, info, alt="<DEL>", end_header=True, name="sv.vcf"):
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=2400>"]
    if end_header:
        lines.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">')
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts")
    lines.append(f"chr1\t100\t.\t{_base(100)}\t{alt}\t.\tPASS\t{info}\tGT\t1")
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_a_symbolic_deletion_with_end_covers_its_whole_span(tmp_path, reference):
    """`<DEL>` at POS 100 with END=400 affects [99, 400). Reading only the REF
    padding base made every site inside a deletion intact."""
    from neoswga.core.position_cache import PositionCache

    table = open_variants(_symbolic_vcf(tmp_path, "END=400"), reference)
    strain = table.strain("s")
    assert list(strain.starts) == [99]
    assert list(strain.ends) == [400]
    assert any("symbolic" in note for note in table.notes)

    primer = SEQUENCE[250:262]  # inside the deletion, far from its padding base
    prefix = str(pathlib.Path(reference).with_suffix(""))
    cache = PositionCache(
        [prefix], [primer], genome_paths=[reference], circular=False, on_missing="scan"
    )
    block = _panel(cache, prefix, [primer], len(SEQUENCE), table, extension=300)
    record = block["per_strain"]["s"]
    inside = [
        int(s) for s in cache.get_positions(prefix, primer) if int(s) >= 99 and int(s) + 12 <= 400
    ]
    assert inside, "the test primer must have a site inside the deletion"
    assert record["affected_sites"]["value"] >= len(inside)


def test_a_symbolic_allele_without_end_is_refused(tmp_path, reference):
    with pytest.raises(ReferenceDataError, match="INFO/END"):
        open_variants(_symbolic_vcf(tmp_path, ".", name="no_end.vcf"), reference)


def test_a_symbolic_allele_whose_end_header_is_missing_is_refused(tmp_path, reference):
    """Without the header line pysam does not read END as the span, so the
    span is unknown and the row is refused rather than read as one base."""
    path = _symbolic_vcf(tmp_path, "END=400", end_header=False, name="no_header.vcf")
    with pytest.raises(ReferenceDataError, match="INFO/END"):
        open_variants(path, reference)


def test_a_breakend_is_refused(tmp_path, reference):
    path = _symbolic_vcf(tmp_path, ".", alt=f"{_base(100)}[chr1:900[", name="bnd.vcf")
    with pytest.raises(ReferenceDataError, match="breakend"):
        open_variants(path, reference)


def test_a_symbolic_alt_in_a_tsv_is_refused(tmp_path, reference):
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt"], [("chr1", 100, _base(100), "<DEL>")])
    with pytest.raises(ReferenceDataError, match="symbolic or breakend"):
        open_variants(path, reference)


def test_a_plain_gzipped_vcf_is_refused_with_the_bgzip_advice(tmp_path, reference):
    """pysam raises NotImplementedError while iterating a gzip (not bgzip)
    file. It used to escape as an unexpected error."""
    import gzip

    rows = [("chr1", 50, ".", _base(50), _other(_base(50)), ".", "PASS", ".")]
    plain = pathlib.Path(_vcf(tmp_path, rows))
    zipped = tmp_path / "variants.vcf.gz"
    with gzip.open(zipped, "wt") as handle:
        handle.write(plain.read_text())
    with pytest.raises(ReferenceDataError, match="bgzip"):
        open_variants(str(zipped), reference)


def test_a_tsv_row_longer_than_its_header_is_refused(tmp_path, reference):
    path = tmp_path / "long.tsv"
    path.write_text(f"chrom\tpos\tref\talt\tstrain_a\nchr1\t50\t{_base(50)}\tA\t1\t0\n")
    with pytest.raises(ReferenceDataError, match="column"):
        open_variants(str(path), reference)


def test_a_genotype_naming_an_allele_the_row_lacks_is_refused(tmp_path, reference):
    """ALT `.` with GT 1/1. pysam reports the allele as missing, which would
    have made the strain unavailable for a reason that misdescribes the row."""
    rows = [("chr1", 50, ".", _base(50), ".", ".", "PASS", ".", "GT", "1/1")]
    with pytest.raises(ReferenceDataError, match="inconsistent"):
        open_variants(_vcf(tmp_path, rows, samples=("s",)), reference)


def test_a_second_alt_allele_beyond_the_row_is_refused(tmp_path, reference):
    rows = [("chr1", 50, ".", _base(50), _other(_base(50)), ".", "PASS", ".", "GT", "0/2")]
    with pytest.raises(ReferenceDataError, match="inconsistent"):
        open_variants(_vcf(tmp_path, rows, samples=("s",)), reference)


# ----------------------------------------------------------------------
# A stale index against the FASTA the table is placed on
# ----------------------------------------------------------------------


class _StaleIndexCache:
    """Answers like a PositionCache whose index was built over another FASTA."""

    def __init__(self, starts):
        self._starts = starts

    def get_record_starts(self, prefix):
        return list(self._starts)

    def get_positions(self, prefix, primer, strand="both"):
        return np.array([], dtype=np.int64)


def _variant_args(path, linear=True, named=None):
    import argparse

    return argparse.Namespace(variants=path, variants_reference=named, linear=linear)


def test_record_starts_from_a_stale_index_are_refused(tmp_path, two_record_reference):
    """`verify_layout` compares the record starts the index stored with the
    FASTA's own. A FASTA edited since `filter` puts every variant at an offset
    the positions no longer use."""
    from neoswga.cli.evaluate import _variant_blocks

    rows = [("chrA", 10, _base(10), _other(_base(10)))]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt"], rows)
    stale = _StaleIndexCache([0, 2000])  # the FASTA's second record starts at 2400
    prefix = str(tmp_path / "two")
    with pytest.raises(ReferenceDataError, match="filter"):
        _variant_blocks(
            _variant_args(path),
            ["ACGTTGCAGGTA"],
            [],
            stale,
            ([prefix], [two_record_reference], [len(SEQUENCE) + len(SECOND)]),
            300,
        )


def test_matching_record_starts_are_accepted(tmp_path, two_record_reference):
    from neoswga.cli.evaluate import _variant_blocks

    rows = [("chrA", 10, _base(10), _other(_base(10)))]
    path = _tsv(tmp_path, ["chrom", "pos", "ref", "alt"], rows)
    current = _StaleIndexCache([0, len(SEQUENCE)])
    prefix = str(tmp_path / "two")
    block = _variant_blocks(
        _variant_args(path),
        ["ACGTTGCAGGTA"],
        [],
        current,
        ([prefix], [two_record_reference], [len(SEQUENCE) + len(SECOND)]),
        300,
    )
    assert block["reference_fasta"] == two_record_reference
