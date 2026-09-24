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
from neoswga.core.exceptions import ReferenceDataError
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


def open_alignment(bam_path, require_index: bool = True, reference: Optional[str] = None):
    """Open a BAM or CRAM, or refuse in a way that names the remedy.

    pysam's own answers to the three things a user hits first are a pysam
    internal's answers. A file that is not an alignment gives "file has no
    sequences defined (mode='rb') ... Consider opening with check_sq=False",
    which is advice for a different problem and names neither the file nor what
    is wrong with it. A BAM without an index gives "fetch called on bamfile
    without index", which names neither the file nor `samtools index`.

    The mode stays "rb" and that is not a bug: htslib detects the format from
    the file's magic bytes, so a CRAM opens through it. Verified by reading one
    back.

    `require_index` is a parameter because reading the HEADER -- contig names
    and lengths, which is all `bind_bam_depth` and `bam_low_depth_gaps` want --
    needs no index, and demanding one there would refuse a file this code can
    read perfectly well.

    `reference` is the FASTA a CRAM was compressed against. CRAM stores
    differences from a reference rather than sequence, so without it the
    records cannot be decoded at all. htslib will look in the `UR` header
    field, then `REF_PATH` and `REF_CACHE`, then the EBI -- so a CRAM may read
    on one machine and not on another with no change to the file, and passing
    the path explicitly is the only way to make it deterministic.
    """
    pysam = _require_pysam()
    name = os.path.basename(str(bam_path))

    if not os.path.exists(bam_path):
        raise ReferenceDataError(
            f"alignment file {name}",
            f"{bam_path} does not exist",
            "check the path passed to --bam",
        )
    if reference is not None and not os.path.exists(reference):
        raise ReferenceDataError(
            f"reference FASTA {os.path.basename(str(reference))}",
            f"{reference} does not exist",
            "check the path passed to --reference",
        )
    try:
        handle = pysam.AlignmentFile(bam_path, "rb", reference_filename=reference)
    except ValueError as exc:
        raise ReferenceDataError(
            f"alignment file {name}",
            f"{bam_path} is not readable as BAM or CRAM",
            "check that it is an aligner's output and not a FASTQ, SAM or "
            "truncated file; `samtools quickcheck` reports the same thing",
        ) from exc  # pysam's own text is kept on __cause__ rather than inlined:
        # it ends "Consider opening with check_sq=False", which is advice for a
        # headerless SAM and would read here as a remedy for the wrong problem.

    if require_index and not handle.has_index():
        handle.close()
        raise ReferenceDataError(
            f"alignment index for {name}",
            f"{bam_path} has no .bai, .csi or .crai beside it, and depth is " f"read by region",
            f"samtools index {bam_path}",
        )
    return handle


def _strip_chr(name: str) -> str:
    return name[3:] if name.lower().startswith("chr") else name


def require_matching_targets(fg_prefixes, fg_seq_lengths) -> None:
    """Refuse a prefix list and a length list that do not correspond.

    These come from two separate params.json keys and are paired with `zip`,
    which stops at the shorter one. A file naming three targets and two lengths
    therefore analysed two of them, with no error and no warning, and every gap
    in the third was reported as absent -- so `expand-primers`, whose job is to
    add oligos for the gaps, would never add one for that target.

    `optimizer_factory.create` already refuses this, which is why the optimizer
    path was safe and this one was not: the coverage commands do not build an
    optimizer. The rule belongs in both places rather than in a shared module
    neither would naturally import.
    """
    if len(fg_prefixes) != len(fg_seq_lengths):
        raise ValueError(
            f"fg_prefixes names {len(fg_prefixes)} target(s) and fg_seq_lengths "
            f"gives {len(fg_seq_lengths)} length(s). They are paired positionally, "
            f"so a mismatch silently drops a target from the analysis. Check "
            f"params.json: the two lists must correspond one to one."
        )


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
    require_matching_targets(fg_prefixes, fg_seq_lengths)
    aliases = aliases or {}
    ref_set = set(bam_refs)
    ref_by_stripped = {_strip_chr(r): r for r in bam_refs}
    # Kept to report a disagreement on a NAME match, not to make one.
    length_by_ref: Dict[str, int] = {
        str(r): int(ln) for r, ln in zip(bam_refs, bam_ref_lengths, strict=True)
    }

    mapping: Dict[str, str] = {}
    for prefix, length in zip(fg_prefixes, fg_seq_lengths, strict=True):
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


def _require_contig_covers(bam, bam_path, contig, length):
    """Refuse to report depth for bases the BAM cannot answer for.

    `count_coverage` CLAMPS its `stop` to the contig's length, and the caller
    then wrote the shorter result into a `length`-sized array of zeros. So
    asking for 5,000 bases of a 2,000 bp contig returned 3,000 fabricated
    zeros -- and a zero here is not "no reads", it is no sequence to have
    reads on. Measured on a fully covered 2 kb contig, `bam_gaps` reported a
    gap of (1950, 5000): 3,050 bp, of which 3,000 bp is invented.
    `expand-primers` designs oligos AT gaps, so those oligos would target a
    region the BAM says nothing whatever about, and `calibrate-reach` would
    fit a polymerase reach against the same invented zeros.

    Reachable only through `match_contigs`, which binds on a NAME and warns
    rather than refusing when the lengths disagree. That is deliberate -- a
    name is an assertion this code should not overrule -- and it stays.
    Binding on the name and inventing the depth are separate decisions, and
    only the second is wrong. `ReferenceLayout.bind` refuses a length mismatch
    outright, so every path carrying `fg_genomes` was already immune.

    A contig LONGER than the configured length is fine and stays silent: that
    reads a prefix of it, and every base reported was observed.
    """
    lengths = dict(zip(bam.references, bam.lengths, strict=True))
    contig_length = lengths.get(contig)
    if contig_length is None:
        return  # pysam raises its own error next, and it names the contig
    if length <= int(contig_length):
        return
    raise ReferenceDataError(
        f"contig {contig} in {os.path.basename(str(bam_path))}",
        f"depth was requested for {length:,} bases but the contig is only "
        f"{int(contig_length):,} bp, so the remaining "
        f"{length - int(contig_length):,} would be reported as zero depth "
        f"rather than as unobserved",
        "align against the same assembly this design uses, or map the right "
        "contig with --contig-alias",
    )


def compute_bam_depth(
    bam_path: str,
    contig: str,
    length: int,
    policy: Optional[DepthPolicy] = None,
    reference: Optional[str] = None,
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
    depth = np.zeros(length, dtype=np.int32)
    with open_alignment(bam_path, reference=reference) as bam:
        _require_contig_covers(bam, bam_path, contig, length)
        # count_coverage returns 4 arrays (A,C,G,T) of length (stop-start).
        try:
            cov = bam.count_coverage(
                contig,
                start=0,
                stop=length,
                quality_threshold=policy.min_base_quality,
                read_callback=policy.accepts,
            )
        except OSError as exc:
            # htslib says "truncated file" when a CRAM's reference cannot be
            # resolved, which is a claim about the CRAM and is wrong: the file
            # is intact and the FASTA is what is missing. A BAM carries its own
            # sequence, so it cannot reach here for this reason.
            if not bam.is_cram:
                raise
            raise ReferenceDataError(
                f"reference for CRAM {os.path.basename(str(bam_path))}",
                "its records could not be decoded, which for a CRAM means the "
                "reference FASTA it was compressed against was not found",
                "pass --reference with the FASTA the reads were aligned to",
            ) from exc
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


def _record_keyed_aliases(layout, aliases):
    """Translate a PREFIX-keyed `--contig-alias` into a RECORD-keyed one.

    The two binding rules key their aliases differently and the CLI documents
    only one of them. `match_contigs`, which `calibrate-reach` uses, keys on
    the foreground PREFIX or its basename -- which is what
    `--contig-alias FG=BAMCONTIG` says. `ReferenceLayout.bind`, which every
    path carrying `fg_genomes` uses and which is therefore the production
    default, keys on the FASTA RECORD name.

    So the documented form was the one that failed on the commoner path.
    Measured on a single-record `mygenome.fasta` whose record is `contig_A`
    against a BAM contig `BAMNAME`: `mygenome=BAMNAME` produced 0 gaps and
    `contig_A=BAMNAME` produced 1. The run still succeeded, having quietly
    ignored the sequencing data `--bam` exists to use.

    A prefix alias is unambiguous only when the reference holds ONE record,
    since then there is exactly one thing it can mean. On a multi-record
    reference a prefix names a file and a BAM contig names a molecule, so
    there is no sound translation; that is warned about rather than guessed,
    because guessing here binds a whole file's depth to one contig.

    An alias already keyed on a record name always wins: it is the more
    specific claim and the one `bind` documents.
    """
    if not aliases:
        return aliases

    record_names = set(layout.names)
    prefix_keys = {layout.prefix, os.path.basename(layout.prefix)}
    translated = dict(aliases)

    for key, value in aliases.items():
        if key in record_names or key not in prefix_keys:
            continue
        if len(layout.records) == 1:
            only = layout.records[0].name
            translated.setdefault(only, value)
        else:
            logger.warning(
                "--contig-alias %r names the foreground prefix, but %s holds %d "
                "records and a prefix is a FILE rather than a contig. Alias the "
                "record instead, by its FASTA header name: %s.",
                key,
                os.path.basename(layout.path),
                len(layout.records),
                ", ".join(layout.names[:5]),
            )
    return translated


def bam_depth_profile(
    bam_path: str,
    prefix: str,
    fasta_path: str,
    configured_length: Optional[int] = None,
    record_starts: Optional[Sequence[int]] = None,
    aliases: Optional[Dict[str, str]] = None,
    policy: Optional[DepthPolicy] = None,
    reference: Optional[str] = None,
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
    layout = read_layout(fasta_path, prefix=prefix)
    verify_layout(layout, record_starts or [], configured_length)

    # Header only, so no index is needed to answer it.
    with open_alignment(bam_path, require_index=False, reference=reference) as bam:
        bam_lengths = {
            name: int(length) for name, length in zip(bam.references, bam.lengths, strict=True)
        }

    bound = layout.bind(bam_lengths, aliases=_record_keyed_aliases(layout, aliases))
    depth = np.zeros(layout.total_length, dtype=np.int32)
    evaluable = np.zeros(layout.total_length, dtype=bool)

    for name in bound.matched:
        record = layout.record(name)
        record_depth = compute_bam_depth(
            bam_path,
            bound.bam_name_for[name],
            record.length,
            policy=policy,
            reference=reference,
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

    runs = list(zip(starts.tolist(), ends.tolist(), strict=True))

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
    reference,
):
    """Gaps found per RECORD, so none can span a join between two molecules.

    A gap is sought inside each bound record's own span. Unbound records
    contribute nothing: the BAM says nothing about them, and inventing a
    coverage hole from absent evidence is the defect rather than the fix.

    `circular` applies only to a genuinely single-record prefix. Wrapping the
    first and last records of a multi-record file would join two molecules.
    """
    all_gaps: List[CoverageGap] = []
    for prefix, length, genome in zip(fg_prefixes, fg_seq_lengths, fg_genomes, strict=True):
        profile = bam_depth_profile(
            bam_path,
            prefix=prefix,
            fasta_path=genome,
            configured_length=length,
            record_starts=record_starts_by_prefix.get(prefix),
            aliases=contig_aliases,
            reference=reference,
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
    reference: Optional[str] = None,
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
    require_matching_targets(fg_prefixes, fg_seq_lengths)
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
            reference,
        )

    with open_alignment(bam_path, require_index=False, reference=reference) as bam:
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
    length_by_prefix = dict(zip(fg_prefixes, fg_seq_lengths, strict=True))
    for prefix, contig in mapping.items():
        length = length_by_prefix[prefix]
        depth = compute_bam_depth(bam_path, contig, length, reference=reference)
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
