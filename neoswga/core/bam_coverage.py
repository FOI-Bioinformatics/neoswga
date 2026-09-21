"""Coverage-gap detection from real sequencing depth (BAM).

After a SWGA primer set is synthesized and used, reads are mapped back to the
target genome. Regions that amplified poorly show up as low sequencing depth.
This module turns a mapped BAM into the same ``CoverageGap`` objects that
``primer_expansion.identify_gaps`` produces from in-silico binding positions,
so the two gap sources can be merged and used to add oligos that fill the
real holes.

pysam is an optional dependency (the ``[bam]`` extra). Importing this module
is cheap; the pysam import is deferred to the functions that need it.
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np

from neoswga.core.depth_policy import DepthPolicy
from neoswga.core.primer_expansion import CoverageGap
from neoswga.core.reference_layout import BoundRecords, read_layout, verify_layout

logger = logging.getLogger(__name__)


def _require_pysam():
    """Import pysam or raise an actionable error."""
    try:
        import pysam  # noqa: F401

        return pysam
    except ImportError as e:  # pragma: no cover - exercised via monkeypatch
        raise RuntimeError(
            "BAM coverage support requires pysam. Install it with:\n"
            "    pip install 'neoswga[bam]'\n"
            "(or `pip install pysam`)."
        ) from e


def _strip_chr(name: str) -> str:
    return name[3:] if name.lower().startswith("chr") else name


def match_contigs(
    bam_refs: Sequence[str],
    bam_ref_lengths: Sequence[int],
    fg_prefixes: Sequence[str],
    fg_seq_lengths: Sequence[int],
    aliases: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """Map foreground prefixes to BAM reference (contig) names.

    Strategy, in order: explicit alias, exact match, basename match, then
    chr-prefix normalization. Every one of those is a NAME agreeing with a
    name, which is a claim somebody made. A prefix with no such match is
    omitted and logged; BAM references matching nothing are ignored.

    **Length is not identity, and is no longer a fallback.** A "unique
    sequence-length match" bound a BAM contig to a foreground reference purely
    because the two were the same size and nothing else in the BAM was. Equal
    length happens: this repository ships two plasmids of 5,386 bp each, and
    two chromosomes from different assemblies routinely agree. When it
    happens every coordinate lines up, so the depth profile reads cleanly
    against a sequence the design was not made for -- and sequencing feedback
    drives redesign, so the low-depth regions it reports become targeted
    additions aimed at gaps in the wrong genome.

    The remedy the warning names is `--contig-alias`, which is the same claim
    made by someone who can check it.

    A name match whose lengths DISAGREE still binds, because the names are an
    assertion and this code should not overrule it, but it warns: that
    combination means the BAM was aligned against a different version of the
    sequence this design used.

    Args:
        bam_refs: BAM @SQ reference names (header order).
        bam_ref_lengths: BAM @SQ reference lengths, parallel to bam_refs.
        fg_prefixes: foreground HDF5 prefixes (the genome coordinate space).
        fg_seq_lengths: foreground lengths, parallel to fg_prefixes.
        aliases: optional explicit {fg_prefix_or_basename: bam_ref} overrides.

    Returns:
        dict mapping ``fg_prefix -> bam_ref`` (only confident matches).
    """
    aliases = aliases or {}
    ref_set = set(bam_refs)
    ref_by_stripped = {_strip_chr(r): r for r in bam_refs}
    # Kept to report a disagreement on a NAME match, not to make one.
    length_by_ref: Dict[str, int] = {str(r): int(ln) for r, ln in zip(bam_refs, bam_ref_lengths)}

    mapping: Dict[str, str] = {}
    for prefix, length in zip(fg_prefixes, fg_seq_lengths):
        base = os.path.basename(prefix)

        matched = None

        # 1. explicit alias (by full prefix or basename)
        if prefix in aliases and aliases[prefix] in ref_set:
            matched = aliases[prefix]
        elif base in aliases and aliases[base] in ref_set:
            matched = aliases[base]
        # 2. exact match on prefix or basename
        elif prefix in ref_set:
            matched = prefix
        elif base in ref_set:
            matched = base
        # 3. chr-prefix normalization (chrI <-> I)
        elif _strip_chr(base) in ref_by_stripped:
            matched = ref_by_stripped[_strip_chr(base)]

        if matched is not None:
            mapping[prefix] = matched
            recorded = length_by_ref.get(matched)
            if recorded is not None and int(recorded) != int(length):
                # The names agree and the sequences cannot both be right. Bind
                # it, because the name is somebody's assertion, but say so: a
                # BAM aligned against another version of this reference puts
                # every downstream coordinate slightly elsewhere.
                logger.warning(
                    "BAM contig '%s' is %d bp and foreground prefix '%s' is %d bp. "
                    "The names match, so they are treated as the same reference, "
                    "but one of them is a different version of the sequence.",
                    matched,
                    int(recorded),
                    prefix,
                    int(length),
                )
            continue

        logger.warning(
            "Could not match foreground prefix '%s' (%d bp) to any BAM contig "
            "by name (%s). It will be skipped for BAM-gap detection. Matching "
            "on sequence length alone is deliberately not attempted, because "
            "equal length is not identity; pass --contig-alias to map it "
            "explicitly.",
            prefix,
            length,
            ", ".join(bam_refs) or "<none>",
        )

    return mapping


def compute_bam_depth(
    bam_path: str,
    contig: str,
    length: int,
    policy: Optional[DepthPolicy] = None,
) -> np.ndarray:
    """Return a per-base depth array (int32, len ``length``) for ``contig``.

    Uses ``pysam.AlignmentFile.count_coverage`` (sum of A/C/G/T per base).
    Positions beyond the array are ignored; missing positions are 0.

    `policy` states which reads count. It used to be pysam's `'all'` callback
    with `quality_threshold=0`, which excluded duplicates -- wrong for
    hyperbranched amplification, where identical start coordinates are
    independent priming events -- and counted supplementary alignments, which
    counts one chimeric molecule in several places. See `core/depth_policy.py`.
    """
    policy = policy or DepthPolicy()
    pysam = _require_pysam()
    depth = np.zeros(length, dtype=np.int32)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # count_coverage returns 4 arrays (A,C,G,T) of length (stop-start).
        cov = bam.count_coverage(
            contig,
            start=0,
            stop=length,
            quality_threshold=policy.min_base_quality,
            read_callback=policy.accepts,
        )
        per_base = np.asarray(cov, dtype=np.int64).sum(axis=0)
        n = min(len(per_base), length)
        depth[:n] = per_base[:n].astype(np.int32)
    return depth


@dataclass(frozen=True)
class DepthProfile:
    """Observed depth for one prefix, and where it could be observed at all.

    `evaluable` is the half that did not exist. A base is evaluable when a BAM
    record bound to its FASTA record at a matching length, so a zero in `depth`
    means measured zero there and nothing at all elsewhere. Reporting the two
    denominators together is what stops masking a hard region inflating
    apparent recovery.
    """

    prefix: str
    depth: np.ndarray
    evaluable: np.ndarray
    bound: "BoundRecords"
    policy: DepthPolicy

    @property
    def evaluable_bases(self) -> int:
        return int(self.evaluable.sum())

    @property
    def total_bases(self) -> int:
        return int(self.evaluable.size)


def bam_depth_profile(
    bam_path: str,
    prefix: str,
    fasta_path: str,
    configured_length: Optional[int] = None,
    record_starts: Optional[Sequence[int]] = None,
    aliases: Optional[Dict[str, str]] = None,
    policy: Optional[DepthPolicy] = None,
) -> DepthProfile:
    """Depth in the prefix's concatenated space, with a non-evaluable mask.

    Binds BAM references to FASTA RECORDS rather than mapping a prefix to a
    contig, which is finding F8: a prefix is a file, not a contig, so a
    multi-record reference matched nothing and the BAM contributed nothing.

    A record whose length disagrees with the BAM's is refused rather than
    accepted on its name, because that mismatch reported a whole record as zero
    depth -- and zero depth is what BAM-guided expansion targets.
    """
    policy = policy or DepthPolicy()
    pysam = _require_pysam()
    layout = read_layout(fasta_path, prefix=prefix)
    verify_layout(layout, record_starts or [], configured_length)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        bam_lengths = {name: int(length) for name, length in zip(bam.references, bam.lengths)}

    bound = layout.bind(bam_lengths, aliases=aliases)
    depth = np.zeros(layout.total_length, dtype=np.int32)
    evaluable = np.zeros(layout.total_length, dtype=bool)

    for name in bound.matched:
        record = layout.record(name)
        record_depth = compute_bam_depth(
            bam_path, bound.bam_name_for[name], record.length, policy=policy
        )
        depth[record.start : record.end] = record_depth
        evaluable[record.start : record.end] = True

    if bound.unmatched_records:
        logger.warning(
            "%d of %d record(s) in %s are not covered by the BAM (%s); their "
            "bases are reported as NOT EVALUABLE rather than as zero depth.",
            len(bound.unmatched_records),
            len(layout.records),
            os.path.basename(fasta_path),
            ", ".join(bound.unmatched_records[:5]),
        )
    logger.info("%s", policy.describe())
    return DepthProfile(prefix=prefix, depth=depth, evaluable=evaluable, bound=bound, policy=policy)


def find_low_depth_gaps(
    depth: np.ndarray,
    prefix: str,
    min_depth: int,
    min_gap_size: int,
    circular: bool = False,
) -> List[CoverageGap]:
    """Find runs of low sequencing depth and return them as CoverageGaps.

    A position is "covered" when ``depth >= min_depth``. Contiguous runs of
    uncovered positions at least ``min_gap_size`` long become gaps, in the same
    coordinate space as ``primer_expansion.identify_gaps``.

    On a circular genome a low-depth run touching the last base is merged with
    one touching the first base (single origin-spanning gap).
    """
    length = len(depth)
    if length == 0:
        return []

    low = depth < min_depth  # True where uncovered
    if not low.any():
        return []

    # Find run boundaries via diff on the int view.
    padded = np.concatenate(([0], low.view(np.int8), [0]))
    diffs = np.diff(padded)
    starts = np.where(diffs == 1)[0]
    ends = np.where(diffs == -1)[0]  # exclusive

    runs = list(zip(starts.tolist(), ends.tolist()))

    # Circular merge: a run ending at `length` joins a run starting at 0.
    if circular and len(runs) >= 2:
        first_s, first_e = runs[0]
        last_s, last_e = runs[-1]
        if first_s == 0 and last_e == length:
            merged_size = (first_e - first_s) + (last_e - last_s)
            # Represent as the wrap gap: start at last_s, end at first_e+length.
            wrap = (last_s, first_e + length, merged_size)
            inner = runs[1:-1]
            gaps = [
                CoverageGap(chromosome=prefix, start=s, end=e, size=e - s)
                for s, e in inner
                if (e - s) >= min_gap_size
            ]
            if wrap[2] >= min_gap_size:
                gaps.append(
                    CoverageGap(chromosome=prefix, start=wrap[0], end=wrap[1], size=wrap[2])
                )
            gaps.sort(key=lambda g: g.size, reverse=True)
            return gaps

    gaps = [
        CoverageGap(chromosome=prefix, start=s, end=e, size=e - s)
        for s, e in runs
        if (e - s) >= min_gap_size
    ]
    gaps.sort(key=lambda g: g.size, reverse=True)
    return gaps


def _bam_gaps_by_record(
    bam_path,
    fg_prefixes,
    fg_seq_lengths,
    fg_genomes,
    min_depth,
    min_gap_size,
    circular,
    contig_aliases,
    record_starts_by_prefix,
):
    """Gaps found per RECORD, so none can span a join between two molecules.

    A gap is sought inside each bound record's own span. Unbound records
    contribute nothing: the BAM says nothing about them, and inventing a
    coverage hole from absent evidence is the defect rather than the fix.

    `circular` applies only to a genuinely single-record prefix. Wrapping the
    first and last records of a multi-record file would join two molecules.
    """
    all_gaps: List[CoverageGap] = []
    for prefix, length, genome in zip(fg_prefixes, fg_seq_lengths, fg_genomes):
        profile = bam_depth_profile(
            bam_path,
            prefix=prefix,
            fasta_path=genome,
            configured_length=length,
            record_starts=record_starts_by_prefix.get(prefix),
            aliases=contig_aliases,
        )
        layout_records = {r.name: r for r in read_layout(genome, prefix=prefix).records}
        single = len(layout_records) == 1
        for name in profile.bound.matched:
            record = layout_records[name]
            gaps = find_low_depth_gaps(
                profile.depth[record.start : record.end],
                prefix,
                min_depth,
                min_gap_size,
                circular=circular and single,
            )
            for gap in gaps:
                all_gaps.append(
                    CoverageGap(
                        chromosome=prefix,
                        start=gap.start + record.start,
                        end=gap.end + record.start,
                        size=gap.size,
                    )
                )
        logger.info(
            "BAM depth for %r: %d of %d bp evaluable across %d bound record(s); "
            "%d gap(s) at min_depth=%d, min_gap_size=%d.",
            prefix,
            profile.evaluable_bases,
            profile.total_bases,
            len(profile.bound.matched),
            len(all_gaps),
            min_depth,
            min_gap_size,
        )

    all_gaps.sort(key=lambda g: g.size, reverse=True)
    return all_gaps


def bam_gaps(
    bam_path: str,
    fg_prefixes: Sequence[str],
    fg_seq_lengths: Sequence[int],
    min_depth: int = 5,
    min_gap_size: int = 10000,
    circular: bool = False,
    contig_aliases: Optional[Dict[str, str]] = None,
    fg_genomes: Optional[Sequence[str]] = None,
    record_starts_by_prefix: Optional[Dict[str, Sequence[int]]] = None,
) -> List[CoverageGap]:
    """Compute low-depth coverage gaps across all foreground prefixes.

    Returns ``CoverageGap`` objects keyed by ``fg_prefix`` (the same
    coordinate space as ``identify_gaps``), sorted largest-first.

    With ``fg_genomes``, depth is bound per RECORD through
    `reference_layout`, so a multi-record reference works, a gap cannot span a
    join between two molecules, and a record the BAM does not cover yields no
    gaps rather than a whole-record hole. Without it the old prefix-level
    matching runs, which can only handle a single-record reference; a
    multi-record one is named rather than returned as an empty list.
    """
    if fg_genomes:
        return _bam_gaps_by_record(
            bam_path,
            fg_prefixes,
            fg_seq_lengths,
            fg_genomes,
            min_depth,
            min_gap_size,
            circular,
            contig_aliases,
            record_starts_by_prefix or {},
        )

    pysam = _require_pysam()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        bam_refs = list(bam.references)
        bam_ref_lengths = list(bam.lengths)

    mapping = match_contigs(bam_refs, bam_ref_lengths, fg_prefixes, fg_seq_lengths, contig_aliases)
    if not mapping:
        logger.warning(
            "No foreground prefix matched the BAM header, so no BAM gaps were "
            "produced. A prefix is a FASTA FILE and not a contig, so a "
            "multi-record reference cannot match this way: pass fg_genomes so "
            "the records can be bound individually."
        )
        return []

    all_gaps: List[CoverageGap] = []
    length_by_prefix = dict(zip(fg_prefixes, fg_seq_lengths))
    for prefix, contig in mapping.items():
        length = length_by_prefix[prefix]
        depth = compute_bam_depth(bam_path, contig, length)
        gaps = find_low_depth_gaps(depth, prefix, min_depth, min_gap_size, circular=circular)
        logger.info(
            "BAM contig '%s' -> %d low-depth gap(s) (min_depth=%d, min_gap_size=%d)",
            contig,
            len(gaps),
            min_depth,
            min_gap_size,
        )
        all_gaps.extend(gaps)

    all_gaps.sort(key=lambda g: g.size, reverse=True)
    return all_gaps
