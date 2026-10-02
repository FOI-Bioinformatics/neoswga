"""The one door into a variant file, in concatenated reference coordinates.

Target diversity can be stated two ways: one FASTA per strain, which Phase 2
already evaluates reference by reference, or a table of variants against one
reference, which is what a canonical-SNP matrix or a VCF from a mapping
pipeline holds. This module reads the second form and nothing else: it decides
nothing about primers and computes no penalty.

What it returns is deliberately narrow. Per strain, two sorted int64 arrays
giving the start and the end of every variant that strain carries, in the
CONCATENATED coordinate space the position index, `PositionCache` and every
coverage helper work in (a 0-based forward-strand offset in the concatenation
of all FASTA records, no separator). `variant_sites` turns those intervals into
intact and affected binding sites; this module only locates them.

Three rules it keeps.

**One door.** `pysam.VariantFile` is opened here and nowhere else in the
package, the same convention `bam_coverage.open_alignment` holds for alignment
files, and for the same reason: the three things a user hits first -- a file
that is not a variant file, a contig the reference does not have, and a table
made against another assembly -- have actionable answers only if there is one
place to put them. pysam is the optional `[bam]` extra, so the import is lazy
and importing this module stays cheap.

**An artifact that cannot be parsed blocks.** Every refusal is a
`ReferenceDataError`. A row whose REF allele disagrees with the FASTA is the
check that catches a table made against another assembly and, because both
accepted inputs are 1-based, an off-by-one in the coordinates as well. Nothing
is skipped and nothing is repaired.

**Unknown is not absent.** A strain whose genotype a row does not state does
not thereby lack the variant. Such a strain is marked `unavailable` with the
row that did it as the reason, and `variant_sites` reports no figure for it
rather than a figure computed as though it carried nothing.

Coordinate conventions, stated because the two input formats could differ and
here they do not:

- VCF/BCF `POS` is 1-based, as the specification says.
- The TSV reader expects `pos` to be 1-based as well, NOT BED-like 0-based.
  A table written 0-based is refused by the REF check, which is what that check
  is for.
- A variant occupies `[start, start + len(REF))` in concatenated coordinates.
  An insertion therefore occupies the one reference base it is anchored at. A
  symbolic allele occupies `[POS - 1, END)` from `INFO/END`, and is refused
  without one; a breakend is refused.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.reference_layout import ReferenceLayout, read_layout

logger = logging.getLogger(__name__)

__all__ = [
    "UNNAMED_STRAIN",
    "VARIANT_DTYPE",
    "StrainVariants",
    "VariantTable",
    "open_variants",
]

#: Coordinates are int64 everywhere in this package. int32 saturated on hg38
#: and hid host sites (Known Issue 7); a variant interval is the same kind of
#: coordinate as a binding position and carries the same dtype.
VARIANT_DTYPE = np.int64

#: The strain name used when a file names no samples: a plain variant list is
#: one unnamed strain rather than no strains, so the same code path reports it.
UNNAMED_STRAIN = "all_variants"

_TSV_REQUIRED = ("chrom", "pos", "ref", "alt")


def _require_pysam():
    """Import pysam or raise an actionable error.

    The same lazy shape as `bam_coverage._require_pysam`, kept separate so a
    TSV table needs no optional dependency at all.
    """
    try:
        import pysam  # noqa: F401

        return pysam
    except ImportError as e:  # pragma: no cover - exercised via monkeypatch
        raise ReferenceDataError(
            "variant file reader",
            "reading a VCF or BCF requires pysam, which is not installed",
            "install it with `pip install 'neoswga[bam]'`, or supply the "
            "variants as a TSV with columns chrom/pos/ref/alt, which needs no "
            "optional dependency",
        ) from e


@dataclass(frozen=True)
class StrainVariants:
    """One strain's variants, or why that strain could not be read.

    `starts` and `ends` are sorted int64 arrays of the same length, in
    concatenated reference coordinates. For an unavailable strain they are
    empty, which is NOT a statement that the strain carries no variants: the
    status is what says whether the arrays mean anything, and
    `variant_sites` reports no figure for a strain that is unavailable.
    """

    name: str
    starts: np.ndarray
    ends: np.ndarray
    status: str = "measured"
    unavailable: str = ""

    @property
    def measured(self) -> bool:
        return self.status == "measured"

    @property
    def count(self) -> int | None:
        """How many variants this strain carries, or None if unknown."""
        return int(len(self.starts)) if self.measured else None

    def as_dict(self) -> dict[str, object]:
        return {
            "strain": self.name,
            "status": self.status,
            "unavailable": self.unavailable or None,
            "variants": self.count,
        }


@dataclass(frozen=True)
class VariantTable:
    """Every strain a variant file describes, against one reference.

    `records_read` is carried so a caller can say how large the table was,
    which is the figure that distinguishes "no variant falls in any binding
    site" from "the file held nothing".
    """

    path: str
    reference: str
    source: str
    layout: ReferenceLayout
    strains: tuple[StrainVariants, ...]
    records_read: int = 0
    contigs_seen: tuple[str, ...] = ()
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(strain.name for strain in self.strains)

    def strain(self, name: str) -> StrainVariants:
        for strain in self.strains:
            if strain.name == name:
                return strain
        raise KeyError(f"{name!r} is not a strain of {self.path!r}")

    @property
    def measured_strains(self) -> tuple[StrainVariants, ...]:
        return tuple(strain for strain in self.strains if strain.measured)

    def as_dict(self) -> dict[str, object]:
        return {
            "path": self.path,
            "reference": self.reference,
            "source": self.source,
            "records_read": self.records_read,
            "contigs": list(self.contigs_seen),
            "strains": [strain.as_dict() for strain in self.strains],
            "notes": list(self.notes),
        }


@dataclass
class _Row:
    """One parsed variant row, before it is assigned to strains."""

    contig: str
    pos: int  # 1-based, as both accepted formats write it
    ref: str
    carriers: tuple[str, ...]
    unknown: tuple[str, ...]
    line: int
    #: Record-local 0-based exclusive end of the reference span the variant
    #: affects. `pos - 1 + len(ref)` for a sequence allele; `INFO/END` for a
    #: symbolic one, whose REF is only the padding base.
    end: int = 0
    symbolic: bool = False

    def __post_init__(self):
        if self.end <= 0:
            self.end = self.pos - 1 + len(self.ref)


def open_variants(path: str, reference_fasta: str, *, prefix: str | None = None) -> VariantTable:
    """Read a variant file against one reference FASTA.

    Args:
        path: a VCF (plain or bgzipped), a BCF, or a TSV with a header whose
            first four columns are chrom/pos/ref/alt and whose remaining
            columns, if any, are one 0/1 genotype column per strain.
        reference_fasta: the FASTA the table was made against. Its record
            layout supplies the contig-to-offset mapping, and its sequence is
            what each REF allele is checked against.
        prefix: the position-index prefix this reference is known by, carried
            on the layout for reporting. Defaults to the FASTA's basename.

    Returns:
        A `VariantTable`.

    Raises:
        ReferenceDataError: the file does not exist or cannot be parsed; a
            contig is absent from the reference; a REF allele disagrees with
            the FASTA; the records are not sorted. None of these is repaired
            or skipped.
    """
    if not os.path.isfile(path):
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"{path} does not exist",
            "check the path passed to --variants",
        )
    if not os.path.isfile(reference_fasta):
        raise ReferenceDataError(
            f"reference FASTA {os.path.basename(str(reference_fasta))}",
            f"{reference_fasta} does not exist, so no variant coordinate can be "
            f"placed and no REF allele can be checked",
            "pass the FASTA the variant table was made against",
        )

    layout = read_layout(reference_fasta, prefix or os.path.basename(reference_fasta))
    if not layout.records:
        raise ReferenceDataError(
            f"reference FASTA {os.path.basename(str(reference_fasta))}",
            "the FASTA holds no records, so it cannot be the reference a "
            "variant table was made against",
            "check that the file is a FASTA and holds sequence",
        )

    source = "vcf" if _looks_like_vcf(path) else "tsv"
    rows, names = _read_vcf(path) if source == "vcf" else _read_tsv(path)

    _check_sorted(rows, path)
    offsets = _contig_offsets(rows, layout, path)
    _check_reference_alleles(rows, layout, reference_fasta, path)

    strains = _assign_strains(rows, names, offsets)
    return VariantTable(
        path=str(path),
        reference=str(reference_fasta),
        source=source,
        layout=layout,
        strains=strains,
        records_read=len(rows),
        contigs_seen=tuple(dict.fromkeys(row.contig for row in rows)),
        notes=_notes(rows, strains),
    )


# ----------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------


def _looks_like_vcf(path: str) -> bool:
    lowered = str(path).lower()
    return lowered.endswith((".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"))


def _read_vcf(path: str) -> tuple[list[_Row], tuple[str, ...]]:
    """Rows and sample names from a VCF or BCF. The only `VariantFile` call."""
    pysam = _require_pysam()
    rows: list[_Row] = []
    # NotImplementedError is in both guards because it is what pysam raises for
    # a `.vcf.gz` that was gzipped rather than bgzipped ("seek not implemented
    # in files compressed by method 1"), and it surfaces during iteration, not
    # when the file is opened.
    unreadable = (OSError, ValueError, NotImplementedError)
    try:
        handle = pysam.VariantFile(str(path))
    except unreadable as exc:
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"pysam could not open it as a VCF or BCF: {exc}",
            _BGZIP_ADVICE,
        ) from exc

    with handle:
        names = tuple(str(name) for name in handle.header.samples)
        try:
            records: Iterator = iter(handle)
            for index, record in enumerate(records, start=1):
                rows.append(_vcf_row(record, names, index))
        except unreadable as exc:
            raise ReferenceDataError(
                f"variant table {os.path.basename(str(path))}",
                f"the file could not be read past record {len(rows) + 1}: {exc}",
                "a truncated or corrupt variant file blocks rather than "
                "contributing the records it did hold; " + _BGZIP_ADVICE,
            ) from exc
    return rows, names


_BGZIP_ADVICE = (
    "check that the file is a VCF/BCF, and that a .gz file is compressed with "
    "bgzip (`bgzip file.vcf`) rather than gzip"
)


def _raw_genotypes(record, names: Sequence[str]) -> dict[str, tuple[int | None, ...] | None]:
    """Each sample's GT as written, allele by allele, or None when absent.

    Read from the record's own text rather than from `samples[name]["GT"]`,
    because pysam reports an allele index the row does not have (a `1` on a
    row whose ALT is `.`) as None, which is indistinguishable from a missing
    call and is the opposite answer. `.` is None; any other token must be an
    integer.
    """
    fields = str(record).rstrip("\n").split("\t")
    if len(fields) < 9 or not names:
        return {name: None for name in names}
    keys = fields[8].split(":")
    if "GT" not in keys:
        return {name: None for name in names}
    where = keys.index("GT")
    out: dict[str, tuple[int | None, ...] | None] = {}
    for name, value in zip(names, fields[9 : 9 + len(names)], strict=True):
        parts = value.split(":")
        if where >= len(parts) or parts[where] in ("", "."):
            out[name] = None
            continue
        alleles: list[int | None] = []
        for token in parts[where].replace("|", "/").split("/"):
            if token == ".":
                alleles.append(None)
                continue
            try:
                alleles.append(int(token))
            except ValueError as exc:
                raise ReferenceDataError(
                    "variant record",
                    f"record at {record.contig}:{record.pos} has genotype "
                    f"{parts[where]!r} for {name}, which is not a VCF genotype",
                    "a GT field holds allele indices separated by / or |, with . "
                    "for an allele that was not called",
                ) from exc
        out[name] = tuple(alleles)
    return out


def _vcf_row(record, names: Sequence[str], index: int) -> _Row:
    """One VCF record, with the strains that carry it and those that cannot say.

    The genotype rule, per sample:

    - any allele that names an ALT: the strain carries the variant, whatever
      its other alleles say (`./1` is a carrier);
    - otherwise, any allele not called (`.`), or no GT at all: the genotype is
      UNKNOWN and the strain is unavailable, because a half-call such as `0/.`
      does not establish the reference allele and reading it as one would
      report a site as intact in a strain nobody measured there;
    - only a genotype whose every allele is stated and is the reference allele
      is a non-carrier.

    A genotype naming an allele the row does not have is an inconsistent row
    and is refused.
    """
    where = f"record {index} ({record.contig}:{record.pos})"
    reference = str(record.ref or "")
    if not reference:
        raise ReferenceDataError(
            "variant record",
            f"{where} has no REF allele",
            "every row needs a reference allele; it is what places the "
            "variant and what is checked against the FASTA",
        )
    alts = tuple(str(alt) for alt in (record.alts or ()))
    end, symbolic = _vcf_span(record, reference, alts, where)

    carriers: list[str] = []
    unknown: list[str] = []
    for name, alleles in _raw_genotypes(record, names).items():
        if alleles is None:
            unknown.append(name)
            continue
        named = [allele for allele in alleles if allele is not None]
        if any(allele < 0 or allele > len(alts) for allele in named):
            raise ReferenceDataError(
                "variant record",
                f"{where} gives {name} allele(s) {named} but the row has "
                f"{len(alts)} ALT allele(s) ({','.join(alts) or '.'})",
                "a genotype may name only the alleles its row states; the row "
                "is inconsistent and nothing it says can be placed",
            )
        if any(allele > 0 for allele in named):
            carriers.append(name)
        elif len(named) < len(alleles):
            unknown.append(name)
    return _Row(
        contig=str(record.contig),
        pos=int(record.pos),
        ref=reference.upper(),
        carriers=tuple(carriers),
        unknown=tuple(unknown),
        line=index,
        end=end,
        symbolic=symbolic,
    )


def _is_breakend(alt: str) -> bool:
    return "[" in alt or "]" in alt or (len(alt) > 1 and (alt.startswith(".") or alt.endswith(".")))


def _vcf_span(record, reference: str, alts: Sequence[str], where: str) -> tuple[int, bool]:
    """The record-local 0-based exclusive end of the span a row affects.

    A sequence allele affects its REF bases. A symbolic allele (`<DEL>`,
    `<DUP>`, `<INV>`, ...) carries only a padding base in REF, so its span is
    `[POS - 1, END)` from `INFO/END`; without END the span is unknown and the
    row is refused rather than read as one base. A breakend joins this position
    to another and has no span on this reference at all, so it is refused.
    """
    breakends = [alt for alt in alts if _is_breakend(alt)]
    if breakends:
        raise ReferenceDataError(
            "variant record",
            f"{where} has a breakend ALT ({breakends[0]}), which joins this "
            f"position to another rather than changing a span of this reference",
            "remove breakend records, or state the rearrangement as sequence "
            "or as a symbolic allele with INFO/END",
        )
    sequence_end = int(record.pos) - 1 + len(reference)
    if not any(alt.startswith("<") for alt in alts):
        return sequence_end, False
    # pysam folds INFO/END into `stop` when the header declares END, and leaves
    # `stop` at the REF span otherwise, so a stop beyond the padding base is an
    # END that was stated and read.
    stop = int(record.stop)
    if stop <= int(record.pos):
        raise ReferenceDataError(
            "variant record",
            f"{where} has a symbolic ALT ({','.join(alts)}) and no usable "
            f"INFO/END, so the span it affects is unknown",
            "give the record INFO/END (declared in the header), or write the variant as sequence",
        )
    return max(stop, sequence_end), True


def _read_tsv(path: str) -> tuple[list[_Row], tuple[str, ...]]:
    """Rows and strain names from a chrom/pos/ref/alt table.

    `pos` is 1-based, like VCF. A 0-based table is refused by the REF check
    rather than silently shifted by one, which is the whole reason that check
    exists.
    """
    rows: list[_Row] = []
    names: tuple[str, ...] = ()
    header_seen = False
    with open(path) as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            fields = line.lstrip("#").split("\t") if not header_seen else line.split("\t")
            if not header_seen:
                header = [field.strip().lower() for field in fields]
                if tuple(header[:4]) != _TSV_REQUIRED:
                    raise ReferenceDataError(
                        f"variant table {os.path.basename(str(path))}",
                        f"the first line is not a header of chrom/pos/ref/alt: {line!r}",
                        "write a tab-separated header line with the columns "
                        "chrom, pos, ref, alt, followed by one optional 0/1 "
                        "column per strain; pos is 1-based, as in VCF",
                    )
                names = tuple(field.strip() for field in fields[4:])
                header_seen = True
                continue
            rows.append(_tsv_row(fields, names, number, path))
    if not header_seen:
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            "the file holds no header line",
            "write a tab-separated header line with the columns chrom, "
            "pos, ref, alt, followed by one optional 0/1 column per strain",
        )
    return rows, names


def _tsv_row(fields: Sequence[str], names: Sequence[str], number: int, path: str) -> _Row:
    if len(fields) != 4 + len(names):
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"line {number} has {len(fields)} column(s) where the header declares {4 + len(names)}",
            "every row needs exactly the header's columns: a short row cannot "
            "be read as missing genotypes, and a long one carries genotypes "
            "for strains the header does not name",
        )
    alt = fields[3].strip()
    if alt.startswith("<") or _is_breakend(alt):
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"line {number} has a symbolic or breakend ALT {alt!r}, whose span "
            f"a chrom/pos/ref/alt table cannot state",
            "write the variant as sequence, or give it as a VCF record with INFO/END",
        )
    try:
        position = int(fields[1])
    except ValueError as exc:
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"line {number} has a non-integer position {fields[1]!r}",
            "pos is a 1-based integer coordinate, as in VCF",
        ) from exc
    if position < 1:
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"line {number} has position {position}, which is not a 1-based coordinate",
            "this reader expects 1-based positions, as VCF writes them, not BED-like 0-based ones",
        )

    reference = fields[2].strip().upper()
    if not reference:
        raise ReferenceDataError(
            f"variant table {os.path.basename(str(path))}",
            f"line {number} has an empty REF allele",
            "every row needs a reference allele; it is what places the "
            "variant and what is checked against the FASTA",
        )

    carriers: list[str] = []
    unknown: list[str] = []
    for name, value in zip(names, fields[4:], strict=True):
        token = value.strip()
        if token in ("", ".", "NA", "na", "?"):
            unknown.append(name)
        elif token == "0":
            continue
        elif token == "1":
            carriers.append(name)
        else:
            raise ReferenceDataError(
                f"variant table {os.path.basename(str(path))}",
                f"line {number} has genotype {token!r} for strain {name!r}",
                "a genotype column holds 0 (reference), 1 (carries the "
                "variant), or `.` for a site that strain was not genotyped at",
            )
    return _Row(
        contig=fields[0].strip(),
        pos=position,
        ref=reference,
        carriers=tuple(carriers),
        unknown=tuple(unknown),
        line=number,
    )


# ----------------------------------------------------------------------
# Checking
# ----------------------------------------------------------------------


def _check_sorted(rows: Sequence[_Row], path: str) -> None:
    """Refuse records that are not grouped by contig and ascending within one.

    Sorted input is what makes the `searchsorted` intact mask correct, and
    sorting the rows here instead would hide a table that is not what it claims
    to be -- a VCF out of order is usually a concatenation accident, and the
    rows after the break may be against another reference.
    """
    seen: list[str] = []
    previous: _Row | None = None
    for row in rows:
        if previous is not None and row.contig == previous.contig:
            if row.pos < previous.pos:
                raise ReferenceDataError(
                    f"variant table {os.path.basename(str(path))}",
                    f"record {row.line} is at {row.contig}:{row.pos}, before "
                    f"record {previous.line} at {previous.contig}:{previous.pos}",
                    "sort the table by contig and position; an unsorted table "
                    "is usually two files concatenated, and the second may be "
                    "against another reference",
                )
        elif row.contig in seen:
            raise ReferenceDataError(
                f"variant table {os.path.basename(str(path))}",
                f"contig {row.contig!r} appears again at record {row.line} "
                f"after another contig's records",
                "group every contig's records together and sort them by position",
            )
        else:
            seen.append(row.contig)
        previous = row


def _contig_offsets(
    rows: Sequence[_Row], layout: ReferenceLayout, path: str
) -> Mapping[str, tuple[int, int]]:
    """Contig name to (concatenated start, length), refusing an unknown name.

    The mapping comes from `reference_layout`, which already reads record names
    and spans from a FASTA; this does not write a second one.
    """
    known = {record.name: (record.start, record.length) for record in layout.records}
    offsets: dict[str, tuple[int, int]] = {}
    for row in rows:
        if row.contig in known and row.end > known[row.contig][1]:
            # A span that runs off its record would be placed on the next
            # record in concatenated space, which no variant can reach.
            raise ReferenceDataError(
                f"variant table {os.path.basename(str(path))}",
                f"record {row.line} at {row.contig}:{row.pos} spans to base "
                f"{row.end}, past the end of {row.contig} "
                f"({known[row.contig][1]} bp)",
                "the table was made against another assembly, or its END is wrong",
            )
        if row.contig in offsets:
            continue
        if row.contig not in known:
            raise ReferenceDataError(
                f"variant table {os.path.basename(str(path))}",
                f"contig {row.contig!r} (record {row.line}) is not a record of "
                f"{os.path.basename(layout.path)}, which holds "
                f"{', '.join(layout.names[:5])}"
                f"{' and others' if len(layout.names) > 5 else ''}",
                "the variant table and the reference must be the same "
                "assembly, with the same record names",
            )
        offsets[row.contig] = known[row.contig]
    return offsets


def _check_reference_alleles(
    rows: Sequence[_Row], layout: ReferenceLayout, reference_fasta: str, path: str
) -> None:
    """Every REF allele against the FASTA, one record held at a time.

    This is the check that catches a table made against another assembly, and,
    because both accepted formats are 1-based, an off-by-one in the
    coordinates: a SNP table shifted by one base disagrees with the FASTA at
    almost every row.

    Records are streamed rather than concatenated. A host-sized reference would
    not fit otherwise, and the rows for one record are all that is needed while
    that record is in hand.
    """
    from neoswga.core import genome_io

    wanted: dict[str, list[_Row]] = {}
    for row in rows:
        wanted.setdefault(row.contig, []).append(row)
    if not wanted:
        return

    loader = genome_io.GenomeLoader()
    names = layout.names
    sequences = loader.load_genome_streaming(reference_fasta)
    for name, sequence in zip(names, sequences, strict=True):
        for row in wanted.pop(name, []):
            start = row.pos - 1
            observed = sequence[start : start + len(row.ref)].upper()
            if observed == row.ref:
                continue
            raise ReferenceDataError(
                f"variant table {os.path.basename(str(path))}",
                f"record {row.line} states REF {row.ref!r} at {row.contig}:"
                f"{row.pos} but {os.path.basename(layout.path)} has "
                f"{observed or '(past the end of the record)'!r} there",
                "the table was made against another assembly, or its "
                "positions are not 1-based; this reader expects 1-based "
                "positions, as VCF writes them",
            )
        if not wanted:
            break


# ----------------------------------------------------------------------
# Assigning rows to strains
# ----------------------------------------------------------------------


def _assign_strains(
    rows: Sequence[_Row],
    names: Sequence[str],
    offsets: Mapping[str, tuple[int, int]],
) -> tuple[StrainVariants, ...]:
    """One `StrainVariants` per strain, in concatenated coordinates.

    With no genotype columns every row belongs to one unnamed strain, which is
    what a plain variant list is. A strain whose genotype some row does not
    state is unavailable with that row as the reason: an ungenotyped site is
    not a reference allele.
    """
    if not names:
        starts, ends = _intervals(rows, offsets)
        return (StrainVariants(name=UNNAMED_STRAIN, starts=starts, ends=ends),)

    carried: dict[str, list[_Row]] = {name: [] for name in names}
    blocked: dict[str, _Row] = {}
    for row in rows:
        for name in row.carriers:
            carried[name].append(row)
        for name in row.unknown:
            blocked.setdefault(name, row)

    out: list[StrainVariants] = []
    for name in names:
        row = blocked.get(name)
        if row is not None:
            reason = (
                f"strain {name} has no complete genotype at {row.contig}:{row.pos} "
                f"(record {row.line}): an allele was not called and none names "
                f"the variant, so whether it carries that variant is unknown "
                f"rather than no"
            )
            out.append(
                StrainVariants(
                    name=name,
                    starts=np.array([], dtype=VARIANT_DTYPE),
                    ends=np.array([], dtype=VARIANT_DTYPE),
                    status="unavailable",
                    unavailable=reason,
                )
            )
            continue
        starts, ends = _intervals(carried[name], offsets)
        out.append(StrainVariants(name=name, starts=starts, ends=ends))
    return tuple(out)


def _intervals(
    rows: Sequence[_Row], offsets: Mapping[str, tuple[int, int]]
) -> tuple[np.ndarray, np.ndarray]:
    """Sorted int64 `[start, end)` intervals in concatenated coordinates."""
    starts = np.empty(len(rows), dtype=VARIANT_DTYPE)
    ends = np.empty(len(rows), dtype=VARIANT_DTYPE)
    for index, row in enumerate(rows):
        offset = offsets[row.contig][0]
        start = offset + row.pos - 1
        starts[index] = start
        ends[index] = offset + row.end
    order = np.argsort(starts, kind="stable")
    return starts[order], ends[order]


def _notes(rows: Sequence[_Row], strains: Sequence[StrainVariants]) -> tuple[str, ...]:
    notes: list[str] = []
    symbolic = [row for row in rows if row.symbolic]
    if symbolic:
        notes.append(
            f"{len(symbolic)} record(s) carry a symbolic ALT; each affects "
            f"[POS - 1, END) as INFO/END states it, padding base included."
        )
    indels = [row for row in rows if len(row.ref) > 1 and not row.symbolic]
    if indels:
        notes.append(
            f"{len(indels)} record(s) have a reference allele longer than one "
            f"base; each affects every site it overlaps, and the coordinate "
            f"shift it causes downstream is not modelled, so gap lengths in "
            f"the carrying strains are approximate."
        )
    unavailable = [strain.name for strain in strains if not strain.measured]
    if unavailable:
        notes.append(
            "No variant set could be read for " + ", ".join(unavailable) + "; see each reason."
        )
    return tuple(notes)
