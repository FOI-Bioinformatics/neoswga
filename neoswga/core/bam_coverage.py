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
from typing import Dict, List, Optional, Sequence

import numpy as np

from neoswga.core.primer_expansion import CoverageGap

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

    Strategy, in order: explicit alias, exact match, basename match,
    chr-prefix normalization, then unique sequence-length match. Foreground
    prefixes with no confident match are omitted (and logged); BAM references
    that match nothing are ignored.

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
    # length -> bam_refs (for unique-length fallback)
    refs_by_len: Dict[int, List[str]] = {}
    for r, ln in zip(bam_refs, bam_ref_lengths):
        refs_by_len.setdefault(int(ln), []).append(r)

    mapping: Dict[str, str] = {}
    for prefix, length in zip(fg_prefixes, fg_seq_lengths):
        base = os.path.basename(prefix)

        # 1. explicit alias (by full prefix or basename)
        if prefix in aliases and aliases[prefix] in ref_set:
            mapping[prefix] = aliases[prefix]
            continue
        if base in aliases and aliases[base] in ref_set:
            mapping[prefix] = aliases[base]
            continue

        # 2. exact match on prefix or basename
        if prefix in ref_set:
            mapping[prefix] = prefix
            continue
        if base in ref_set:
            mapping[prefix] = base
            continue

        # 3. chr-prefix normalization (chrI <-> I)
        stripped = _strip_chr(base)
        if stripped in ref_by_stripped:
            mapping[prefix] = ref_by_stripped[stripped]
            continue

        # 4. unique sequence-length match
        same_len = refs_by_len.get(int(length), [])
        if len(same_len) == 1:
            mapping[prefix] = same_len[0]
            logger.info(
                "Matched fg prefix '%s' to BAM contig '%s' by length (%d bp)",
                prefix, same_len[0], length,
            )
            continue

        logger.warning(
            "Could not match foreground prefix '%s' (%d bp) to any BAM contig "
            "(%s). It will be skipped for BAM-gap detection; pass "
            "--contig-alias to map it explicitly.",
            prefix, length, ", ".join(bam_refs) or "<none>",
        )

    return mapping


def compute_bam_depth(bam_path: str, contig: str, length: int) -> np.ndarray:
    """Return a per-base depth array (int32, len ``length``) for ``contig``.

    Uses ``pysam.AlignmentFile.count_coverage`` (sum of A/C/G/T per base).
    Positions beyond the array are ignored; missing positions are 0.
    """
    pysam = _require_pysam()
    depth = np.zeros(length, dtype=np.int32)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        # count_coverage returns 4 arrays (A,C,G,T) of length (stop-start).
        cov = bam.count_coverage(
            contig, start=0, stop=length, quality_threshold=0
        )
        per_base = np.asarray(cov, dtype=np.int64).sum(axis=0)
        n = min(len(per_base), length)
        depth[:n] = per_base[:n].astype(np.int32)
    return depth


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
                    CoverageGap(chromosome=prefix, start=wrap[0],
                                end=wrap[1], size=wrap[2])
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


def bam_gaps(
    bam_path: str,
    fg_prefixes: Sequence[str],
    fg_seq_lengths: Sequence[int],
    min_depth: int = 5,
    min_gap_size: int = 10000,
    circular: bool = False,
    contig_aliases: Optional[Dict[str, str]] = None,
) -> List[CoverageGap]:
    """Compute low-depth coverage gaps across all foreground contigs.

    Returns ``CoverageGap`` objects keyed by ``fg_prefix`` (the same
    coordinate space as ``identify_gaps``), sorted largest-first.
    """
    pysam = _require_pysam()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        bam_refs = list(bam.references)
        bam_ref_lengths = list(bam.lengths)

    mapping = match_contigs(
        bam_refs, bam_ref_lengths, fg_prefixes, fg_seq_lengths, contig_aliases
    )
    if not mapping:
        logger.warning(
            "No foreground contigs matched the BAM header; no BAM gaps produced."
        )
        return []

    all_gaps: List[CoverageGap] = []
    length_by_prefix = dict(zip(fg_prefixes, fg_seq_lengths))
    for prefix, contig in mapping.items():
        length = length_by_prefix[prefix]
        depth = compute_bam_depth(bam_path, contig, length)
        gaps = find_low_depth_gaps(
            depth, prefix, min_depth, min_gap_size, circular=circular
        )
        logger.info(
            "BAM contig '%s' -> %d low-depth gap(s) (min_depth=%d, min_gap_size=%d)",
            contig, len(gaps), min_depth, min_gap_size,
        )
        all_gaps.extend(gaps)

    all_gaps.sort(key=lambda g: g.size, reverse=True)
    return all_gaps
