"""Where each FASTA record sits in a prefix's concatenated coordinate space.

Finding F8 of the 2026-09-16 pipeline audit. `bam_coverage.match_contigs` maps
one foreground PREFIX to one BAM reference, matching on name or on total
length. A prefix is not a contig: it is one FASTA file laid out as a single
concatenated space, which is what the position index, `PositionCache` and every
coverage helper work in. So a two-record FASTA matches no BAM contig by either
route, and the BAM contributes nothing at all. A name match is also accepted
without checking the record length, so a partial match reports a whole record
as zero depth.

Both are the silent-zero shape of Known Issues 5, 6 and 13: absent evidence
presented as measured absence. Here that is worse than usual, because the
quantity being measured IS absence -- a region with no observed depth is read
as a region that amplified poorly, which is exactly what BAM-guided expansion
then goes and targets.

This module is the coordinate contract the rest of Phase 6 rests on. It reads
the record names and spans from the FASTA, refuses to proceed when they
disagree with what the position index was built over, and binds BAM references
to records rather than to prefixes.

Two things it does NOT do. It does not read sequence into memory: a foreground
target is usually small but a background is not, and a layout is a few hundred
integers whatever the genome. And it does not guess: a record whose length
disagrees with the BAM's is REFUSED rather than accepted on the strength of its
name, because that mismatch is the one that silently zeroed a whole record.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

logger = logging.getLogger(__name__)


class LayoutMismatch(RuntimeError):
    """The FASTA on disk is not the one the position index was built over.

    Raised rather than warned. Every depth written after a mismatch lands at
    the wrong offset, and a wrong offset produces a plausible coverage figure
    rather than an obviously broken one.
    """


@dataclass(frozen=True)
class Record:
    """One FASTA record, and where it starts in the concatenated space."""

    name: str
    start: int
    length: int

    @property
    def end(self) -> int:
        return self.start + self.length


@dataclass(frozen=True)
class BoundRecords:
    """Which records a BAM actually carries, and at what offsets.

    `unmatched_records` and `length_mismatches` are carried rather than logged
    away, because they are what a non-evaluable mask is built from: a record
    the BAM does not cover has UNKNOWN depth, not zero depth.
    """

    matched: tuple[str, ...]
    offsets: dict[str, int]
    bam_name_for: dict[str, str]
    unmatched_records: tuple[str, ...]
    length_mismatches: dict[str, tuple[int, int]]
    lengths: dict[str, int]

    @property
    def evaluable_length(self) -> int:
        """Bases a BAM can speak for. The denominator that is not a guess."""
        return sum(self.lengths.get(name, 0) for name in self.matched)

    @property
    def total_length(self) -> int:
        """Every base of the prefix, evaluable or not."""
        return sum(self.lengths.values())


@dataclass(frozen=True)
class LayoutCheck:
    """What was verified, so a caller can say which checks actually ran."""

    record_starts_checked: bool
    configured_length_checked: bool


@dataclass(frozen=True)
class ReferenceLayout:
    """The record layout of one prefix's FASTA."""

    prefix: str
    path: str
    records: tuple[Record, ...]
    total_length: int

    def record(self, name: str) -> Record:
        for record in self.records:
            if record.name == name:
                return record
        raise KeyError(f"{name!r} is not a record of {self.path!r}")

    def offset_of(self, name: str) -> int:
        """The record's start. Raises rather than returning 0 for an unknown
        name, since 0 would write one record's depth over another's."""
        return self.record(name).start

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(record.name for record in self.records)

    @property
    def record_starts(self) -> tuple[int, ...]:
        return tuple(record.start for record in self.records)

    def bind(
        self,
        bam_lengths: Mapping[str, int],
        aliases: Mapping[str, str] | None = None,
    ) -> BoundRecords:
        """Match this layout's records to BAM references, by RECORD.

        Order: an explicit alias, then an exact name, then a `chr` prefix
        difference, which is the ordinary Ensembl-against-UCSC case. A length
        that disagrees refuses the match at every step.
        """
        aliases = dict(aliases or {})
        stripped = {_strip_chr(name): name for name in bam_lengths}

        matched, offsets, bam_name_for = [], {}, {}
        unmatched, mismatches = [], {}
        for record in self.records:
            bam_name = _candidate_bam_name(record.name, bam_lengths, stripped, aliases)
            if bam_name is None:
                unmatched.append(record.name)
                continue
            bam_length = int(bam_lengths[bam_name])
            if bam_length != record.length:
                mismatches[record.name] = (record.length, bam_length)
                unmatched.append(record.name)
                logger.warning(
                    "Record %r is %d bp in %s and %d bp in the BAM; refusing the "
                    "match rather than reporting the record as unobserved.",
                    record.name,
                    record.length,
                    os.path.basename(self.path),
                    bam_length,
                )
                continue
            matched.append(record.name)
            offsets[record.name] = record.start
            bam_name_for[record.name] = bam_name

        return BoundRecords(
            matched=tuple(matched),
            offsets=offsets,
            bam_name_for=bam_name_for,
            unmatched_records=tuple(unmatched),
            length_mismatches=mismatches,
            lengths={record.name: record.length for record in self.records},
        )


def _strip_chr(name: str) -> str:
    return name[3:] if name.lower().startswith("chr") else name


def _candidate_bam_name(
    record_name: str,
    bam_lengths: Mapping[str, int],
    stripped: Mapping[str, str],
    aliases: Mapping[str, str],
) -> str | None:
    alias = aliases.get(record_name)
    if alias is not None and alias in bam_lengths:
        return alias
    if record_name in bam_lengths:
        return record_name
    return stripped.get(_strip_chr(record_name))


def read_layout(fasta_path: str, prefix: str) -> ReferenceLayout:
    """Record names and spans, streamed rather than loaded.

    The name is the header's first whitespace-delimited token, which is the
    convention a BAM `@SQ` line follows; keeping the description would make
    every record binding fail on a descriptive FASTA.
    """
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(f"No such FASTA: {fasta_path}")

    records: list = []
    name: str | None = None
    start = 0
    length = 0

    with open(fasta_path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append(Record(name, start, length))
                    start += length
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                length = 0
            elif name is not None:
                length += len(line.strip())
    if name is not None:
        records.append(Record(name, start, length))
        start += length

    return ReferenceLayout(
        prefix=prefix,
        path=fasta_path,
        records=tuple(records),
        total_length=start,
    )


def verify_layout(
    layout: ReferenceLayout,
    record_starts: Sequence[int],
    configured_length: int | None,
) -> LayoutCheck:
    """Refuse a layout that is not what the index and the config describe.

    An EMPTY `record_starts` is an absent check rather than a failed one:
    `PositionCache.get_record_starts` documents it for an index built before
    they were stored. The result says which checks ran, so a caller can report
    "not checked" instead of implying it passed.
    """
    starts_checked = bool(record_starts)
    if starts_checked and tuple(int(s) for s in record_starts) != layout.record_starts:
        raise LayoutMismatch(
            f"{layout.path} has record starts {layout.record_starts} but the "
            f"position index was built over {tuple(record_starts)}. The FASTA "
            f"has changed since the index was written; re-run `neoswga filter` "
            f"rather than reading depth at offsets that no longer mean what "
            f"they meant."
        )

    length_checked = configured_length is not None
    if length_checked and int(configured_length) != layout.total_length:
        raise LayoutMismatch(
            f"{layout.path} has total length {layout.total_length} bp but the run "
            f"is configured for a length of {int(configured_length)} bp. Every "
            f"coverage denominator would be wrong by the difference."
        )

    if not starts_checked:
        logger.info(
            "The position index for %r predates stored record starts, so the "
            "layout could not be cross-checked against it.",
            layout.prefix,
        )
    return LayoutCheck(
        record_starts_checked=starts_checked,
        configured_length_checked=length_checked,
    )
