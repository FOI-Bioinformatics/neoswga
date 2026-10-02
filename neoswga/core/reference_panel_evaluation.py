"""One record per reference genome, and the target-against-host cross table.

Everything the pipeline reports about several references today is POOLED. The
foreground frequency gate sums counts over summed length, the background gate
does the same, and `total_bg_sites`, `selectivity_ratio` and
`selectivity_density` are flat sums over prefixes against one pooled background
length (`base_optimizer` ~1124-1201). Two consequences were measured by reading
those call sites:

- a primer frequent in one target strain and absent from another passes the
  foreground gate exactly as one present in both;
- a primer clean against a 3.3 Gb host and frequent in a 2 Mb host passes the
  background gate, because the small host is a small share of the pooled bases.
  That is Known Issue 6 -- a figure that moves with background size -- in a
  second place.

This module answers the per-reference question after the fact. It measures, and
decides nothing: no threshold is defined here, and no selection stage reads it.

Three rules it keeps, each from a defect this repository has shipped.

**Unknown is not zero.** Every quantity is a `panel_evaluation.Measurement`, so
a reference whose index was never opened carries `value=None` with a reason
rather than a 0.0 that reads as "binds the host nowhere".

**A reduction over a set with one unavailable member is unavailable.** Not the
best of the rest. `host_profile.aggregate_loads` returns 0.0 for an empty list,
which is the same shape as reporting the most permissive answer available, and
is deliberately not copied.

**Prefix and genome travel together.** A reference is one `ReferenceSpec`
carrying both, and the counts route names the genome at the call site, as
`kmer_tables.counts_for` requires: pairing a prefix with a genome from a module
global is how a design once scored itself against another run's reference.

Sites are measured by one of two routes, and the record says which.

  `positions`  the binding positions, from the position index and, where the
               caller allows it, by scanning the FASTA. This is the only route
               that can also yield coverage and gap statistics.
  `counts`     `kmer_tables.counts_for`, which reads a k-mer table where one
               exists and otherwise scans the reference through `query_scan`.
               Coverage is then `unavailable` with that as the reason, which is
               the honest answer for a host too large to hold in memory:
               `docs/validation/query_scan_2026-09-25.md` measures the scan at
               332 MB peak on a 144 Mb reference and 2.3 GB on hg38.

The two agree on the site count. A canonical table stores one spelling of each
reverse-complement pair and answers 0 for the other, so the counts route asks
for both spellings and keeps the larger -- that 0 is not a measurement. The same
record documents the difference, measured at 13 of 1,272 queries at k=18.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from neoswga.core.panel_evaluation import Measurement
from neoswga.core.selectivity import selectivity_density_from_loads, selectivity_from_loads

logger = logging.getLogger(__name__)

__all__ = [
    "ROLE_HOST",
    "ROLE_TARGET",
    "PairRecord",
    "ReferencePanelAssessment",
    "ReferenceRecord",
    "ReferenceSpec",
    "evaluate_reference_panel",
]

ROLE_TARGET = "target"
ROLE_HOST = "host"

#: A reference that was measured, against one that could not be.
STATUS_MEASURED = "measured"
STATUS_UNAVAILABLE = "unavailable"

#: How sites were obtained. Reported, because the two routes answer the same
#: question at very different cost and only one of them also yields coverage.
SOURCE_POSITIONS = "positions"
SOURCE_COUNTS = "counts"


@dataclass(frozen=True)
class ReferenceSpec:
    """One genome to evaluate against, named once.

    `scan` is the cost decision, and it belongs to the caller: finding
    positions in a reference with no index holds that reference in memory, so a
    command evaluating a 3.3 Gb host asks for the counts route instead and gets
    `coverage` as unavailable rather than as zero.
    """

    prefix: str
    genome: str | None
    length: int
    role: str
    circular: bool = False
    scan: bool = True

    def __post_init__(self):
        if self.role not in (ROLE_TARGET, ROLE_HOST):
            raise ValueError(f"role must be {ROLE_TARGET!r} or {ROLE_HOST!r}, got {self.role!r}")


@dataclass(frozen=True)
class ReferenceRecord:
    """What one reference measured, or why it did not."""

    prefix: str
    role: str
    length: int
    status: str
    source: str
    sites: Measurement
    site_density: Measurement
    coverage: Measurement
    mean_gap: Measurement
    max_gap: Measurement
    gap_gini: Measurement
    weighted_load: Measurement
    per_primer_sites: Mapping[str, int | None] = field(default_factory=dict)
    per_primer_status: Mapping[str, str] = field(default_factory=dict)
    unavailable: str = ""
    genome: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return {
            "prefix": self.prefix,
            "genome": self.genome,
            "role": self.role,
            "length_bp": self.length,
            "status": self.status,
            "unavailable": self.unavailable or None,
            "site_source": self.source,
            "sites": self.sites.as_dict(),
            "sites_per_mb": self.site_density.as_dict(),
            "coverage": self.coverage.as_dict(),
            "mean_gap_bp": self.mean_gap.as_dict(),
            "max_gap_bp": self.max_gap.as_dict(),
            "gap_gini": self.gap_gini.as_dict(),
            "weighted_site_load": self.weighted_load.as_dict(),
            "per_primer_sites": dict(self.per_primer_sites),
            "per_primer_status": dict(self.per_primer_status),
        }


@dataclass(frozen=True)
class PairRecord:
    """One target against one host, each measured against its own length.

    `density` is the comparable figure. `ratio` is kept beside it because a
    reader comparing against a published fg/bg number needs it, and because the
    two disagree whenever the hosts differ in size -- which is the whole reason
    this table exists rather than one pooled ratio (Known Issue 6).
    """

    target: str
    host: str
    density: Measurement
    ratio: Measurement
    #: True when the host carries no sites at all. The ratio is then unbounded
    #: rather than large, and `MAX_SELECTIVITY` stands in for it as everywhere
    #: else in this package.
    zero_host_sites: bool = False

    def as_dict(self) -> dict[str, Any]:
        return {
            "target": self.target,
            "host": self.host,
            "selectivity_density": self.density.as_dict(),
            "selectivity_ratio": self.ratio.as_dict(),
            "zero_host_sites": self.zero_host_sites,
            "comparable": "selectivity_density",
        }


@dataclass(frozen=True)
class ReferencePanelAssessment:
    """Every reference, every pair, and the two reductions over them."""

    primers: tuple[str, ...]
    references: tuple[ReferenceRecord, ...]
    pairs: tuple[PairRecord, ...]
    worst_target_coverage: Measurement
    worst_host_density: Measurement
    extension_reach_bp: int
    notes: tuple[str, ...] = ()

    @property
    def targets(self) -> tuple[ReferenceRecord, ...]:
        return tuple(r for r in self.references if r.role == ROLE_TARGET)

    @property
    def hosts(self) -> tuple[ReferenceRecord, ...]:
        return tuple(r for r in self.references if r.role == ROLE_HOST)

    def record(self, prefix: str) -> ReferenceRecord | None:
        for reference in self.references:
            if reference.prefix == prefix:
                return reference
        return None

    def as_dict(self) -> dict[str, Any]:
        return {
            "primers": list(self.primers),
            "extension_reach_bp": self.extension_reach_bp,
            "per_target": {r.prefix: r.as_dict() for r in self.targets},
            "per_host": {r.prefix: r.as_dict() for r in self.hosts},
            "target_host_pairs": [pair.as_dict() for pair in self.pairs],
            "worst_target_coverage": self.worst_target_coverage.as_dict(),
            "worst_host_selectivity_density": self.worst_host_density.as_dict(),
            "notes": list(self.notes),
        }


def evaluate_reference_panel(
    primers: Sequence[str],
    references: Sequence[ReferenceSpec],
    conditions: Any = None,
    *,
    extension: int = 3000,
    max_mismatches: int = 1,
) -> ReferencePanelAssessment:
    """Measure one primer set against every reference, separately.

    Args:
        primers: the panel. Deduplicated, order preserved.
        references: one `ReferenceSpec` per genome, each carrying its own
            prefix, FASTA path, length and role. Lengths are per reference on
            purpose: a host's sites spread over the host's own bases, and
            dividing them by a pooled length is the defect this module exists
            to avoid.
        conditions: `ReactionConditions`, or None. Needed only for the
            occupancy-weighted load, which is `unavailable` without it.
        extension: per-primer reach for coverage, in bp. The caller resolves it
            from the polymerase (`coverage.polymerase_extension_reach`); no
            default polymerase is assumed here.
        max_mismatches: mismatch classes the weighted load sums over.

    Returns:
        A `ReferencePanelAssessment`. Nothing in it is a substituted value: a
        quantity that could not be measured is a `Measurement` with `value=None`
        and a reason, and a reduction over any such member is itself
        unavailable.
    """
    panel = _unique(primers)
    specs = list(references)
    if not panel:
        raise ValueError("evaluate_reference_panel needs at least one primer")

    seen: set[str] = set()
    for spec in specs:
        if spec.prefix in seen:
            raise ValueError(f"reference prefix {spec.prefix!r} appears twice")
        seen.add(spec.prefix)

    records: dict[str, ReferenceRecord] = {}
    # A reference with no length is not a reference that measured zero. Its
    # length is the denominator of every density, coverage figure and gap, and a
    # site count reported beside a zero length reads as a clean genome. An empty
    # or unreadable FASTA arrives here exactly this way.
    for spec in specs:
        if int(spec.length) <= 0:
            records[spec.prefix] = _unavailable_record(
                spec,
                f"the reference has no length ({spec.length} bp), so nothing about "
                f"it can be measured; check that the FASTA holds sequence",
            )

    for group in _position_groups(specs):
        if any(spec.prefix in records for spec in group):
            group = [spec for spec in group if spec.prefix not in records]
        if not group:
            continue
        records.update(_measure_positions(group, panel, extension))
    for spec in specs:
        if spec.prefix not in records:
            records[spec.prefix] = _measure_counts(spec, panel)

    # The weighted load is a separate question with its own prerequisite (a
    # table for the mismatch classes), so a reference can have sites and no
    # weighted load. It is attached per reference rather than folded into the
    # site measurement.
    for prefix, record in list(records.items()):
        records[prefix] = _with_weighted_load(record, panel, conditions, max_mismatches)

    ordered = tuple(records[spec.prefix] for spec in specs)
    pairs = _cross_table(ordered)
    return ReferencePanelAssessment(
        primers=tuple(panel),
        references=ordered,
        pairs=pairs,
        worst_target_coverage=_worst_coverage(ordered),
        worst_host_density=_worst_density(pairs),
        extension_reach_bp=int(extension),
        notes=_notes(ordered),
    )


# ----------------------------------------------------------------------
# The positions route
# ----------------------------------------------------------------------


def _position_groups(specs: Sequence[ReferenceSpec]) -> list[list[ReferenceSpec]]:
    """References whose positions can be looked up, grouped for one call each.

    `compute_per_prefix_coverage` takes one reach and one `circular`, and a
    cache is built either to scan or not to scan, so the group key is
    (role, circular, scanning). Grouping is what lets the per-prefix loop be
    CALLED once per role group rather than reimplemented here, which is how the
    duplicate in `cli/iterate.py` came to swallow an exception per primer and
    report 0.0 coverage.
    """
    groups: dict[tuple[str, bool, bool], list[ReferenceSpec]] = {}
    for spec in specs:
        scanning = _scan_possible(spec)
        if not scanning and not _index_exists(spec.prefix):
            continue
        groups.setdefault((spec.role, bool(spec.circular), scanning), []).append(spec)
    return list(groups.values())


def _scan_possible(spec: ReferenceSpec) -> bool:
    """Whether this reference's FASTA may and can be scanned for positions.

    A missing file is decided here rather than by trying and interpreting the
    failure, so one absent genome does not take its group's cache down with it.
    """
    return bool(spec.scan and spec.genome and os.path.exists(spec.genome))


def _index_exists(prefix: str) -> bool:
    """Whether any position index file exists for this prefix, at any k.

    A directory listing rather than `open_index`: the question is whether a
    file is there, and the k is not known until a primer length is. The file
    itself is only ever read through `position_index`, via `PositionCache`.
    """
    directory = os.path.dirname(prefix) or "."
    stem = os.path.basename(prefix)
    try:
        names = os.listdir(directory)
    except OSError:
        return False
    return any(name.startswith(stem + "_") and name.endswith("mer_positions.h5") for name in names)


def _measure_positions(
    group: Sequence[ReferenceSpec], panel: Sequence[str], extension: int
) -> dict[str, ReferenceRecord]:
    """Sites, coverage and gaps for one group, from one cache and one call.

    A group whose cache or coverage call raises is retried one reference at a
    time, so one unreadable genome leaves the others measured instead of taking
    the group down with it. A single reference that still raises is unavailable
    with the error as its reason.
    """
    specs = list(group)
    try:
        return _measure_group(specs, panel, extension)
    except Exception as exc:
        if len(specs) == 1:
            spec = specs[0]
            logger.debug("Positions unavailable for %s: %s", spec.prefix, exc)
            return {spec.prefix: _unavailable_record(spec, f"position lookup failed: {exc}")}
        out: dict[str, ReferenceRecord] = {}
        for spec in specs:
            out.update(_measure_positions([spec], panel, extension))
        return out


def _measure_group(
    specs: Sequence[ReferenceSpec], panel: Sequence[str], extension: int
) -> dict[str, ReferenceRecord]:
    from neoswga.core.coverage import compute_per_prefix_coverage
    from neoswga.core.position_cache import PositionCache

    prefixes = [spec.prefix for spec in specs]
    lengths = [int(spec.length) for spec in specs]
    circular = bool(specs[0].circular)

    # Scanning resolves a primer the index does not hold, which is the only way
    # an outside oligo set gets a real answer. It needs a genome per prefix, in
    # the same order, and the pairing is made here rather than read from a
    # global.
    genomes = [spec.genome for spec in specs]
    scan = all(_scan_possible(spec) for spec in specs)
    cache = PositionCache(
        prefixes,
        list(panel),
        genome_paths=list(genomes) if scan else None,
        circular=circular,
        on_missing="scan" if scan else "warn",
    )

    _overall, per_prefix = compute_per_prefix_coverage(
        cache,
        list(panel),
        prefixes,
        lengths,
        extension=extension,
        circular=circular,
    )

    # Per prefix, not pooled: a primer the index holds for one reference and not
    # for another is measured on the first and unknown on the second, and one
    # shared set of names would report both as unknown.
    unresolved: dict[str, set[str]] = {}
    for prefix, primer in cache.missing_primers:
        unresolved.setdefault(prefix, set()).add(primer)
    zero_site = set(cache.zero_site_primers)

    out: dict[str, ReferenceRecord] = {}
    for spec in specs:
        absent = unresolved.get(spec.prefix, set())
        per_primer: dict[str, int | None] = {}
        status: dict[str, str] = {}
        for primer in panel:
            if primer in absent:
                per_primer[primer] = None
                status[primer] = "not_indexed"
                continue
            count = int(len(cache.get_positions(spec.prefix, primer)))
            per_primer[primer] = count
            status[primer] = "no_sites" if count == 0 or primer in zero_site else "ok"

        if any(value is None for value in per_primer.values()):
            # Some primer's sites on this reference are unknown, so every figure
            # derived from the union of sites is a lower bound and not a
            # measurement. The per-primer block still says which ones answered.
            reason = (
                f"{sum(1 for v in per_primer.values() if v is None)} primer(s) are "
                f"absent from the position index for {spec.prefix} and no FASTA "
                f"was scanned; their sites are unknown rather than zero"
            )
            out[spec.prefix] = _unavailable_record(
                spec,
                reason,
                per_primer_sites=per_primer,
                per_primer_status=status,
                source=SOURCE_POSITIONS,
            )
            continue

        sites = sum(int(value) for value in per_primer.values())
        mean_gap, max_gap, gini = _gap_statistics(cache, panel, spec)
        out[spec.prefix] = ReferenceRecord(
            prefix=spec.prefix,
            genome=spec.genome,
            role=spec.role,
            length=int(spec.length),
            status=STATUS_MEASURED,
            source=SOURCE_POSITIONS,
            sites=Measurement("sites", float(sites), "sites", basis="exact matches, both strands"),
            site_density=_density(sites, spec.length),
            coverage=Measurement(
                "coverage",
                float(per_prefix.get(spec.prefix, 0.0)),
                "fraction",
                basis=f"union of {int(extension)} bp windows, "
                f"{'circular' if circular else 'linear'}",
            ),
            mean_gap=mean_gap,
            max_gap=max_gap,
            gap_gini=gini,
            weighted_load=Measurement(
                "weighted_site_load", None, "sites", unavailable="not attempted yet"
            ),
            per_primer_sites=per_primer,
            per_primer_status=status,
        )
    return out


def _gap_statistics(
    cache: Any, panel: Sequence[str], spec: ReferenceSpec
) -> tuple[Measurement, Measurement, Measurement]:
    """Mean, maximum and Gini of the distances between sites, or unavailable.

    Fewer than two sites gives no gap at all. `cli/evaluate.py` returns
    (0.0, 0.0, 0.0) there, which reads as perfectly even spacing; here the three
    are unavailable with that as the reason.
    """
    import numpy as np

    unmeasured = "fewer than two binding sites on this reference, so there is no gap"
    positions = sorted(
        {int(x) for primer in panel for x in cache.get_positions(spec.prefix, primer)}
    )
    if len(positions) < 2:
        return (
            Measurement("mean_gap", None, "bp", unavailable=unmeasured),
            Measurement("max_gap", None, "bp", unavailable=unmeasured),
            Measurement("gap_gini", None, "fraction", unavailable=unmeasured),
        )

    gaps = np.diff(positions).astype(float)
    if spec.circular:
        gaps = np.append(gaps, float(spec.length - positions[-1] + positions[0]))
    ordered = np.sort(gaps)
    n = len(ordered)
    total = float(ordered.sum())
    basis = f"{n} gap(s) between {len(positions)} site(s)"
    if total <= 0:
        gini = Measurement(
            "gap_gini",
            None,
            "fraction",
            unavailable="every gap measured zero, so the Gini coefficient is undefined",
        )
    else:
        value = 2.0 * float(np.sum(np.arange(1, n + 1) * ordered)) / (n * total) - (n + 1) / n
        gini = Measurement("gap_gini", float(value), "fraction", basis=basis)
    return (
        Measurement("mean_gap", float(ordered.mean()), "bp", basis=basis),
        Measurement("max_gap", float(ordered.max()), "bp", basis=basis),
        gini,
    )


# ----------------------------------------------------------------------
# The counts route
# ----------------------------------------------------------------------


def _measure_counts(spec: ReferenceSpec, panel: Sequence[str]) -> ReferenceRecord:
    """Sites from `kmer_tables.counts_for`, with coverage unavailable.

    This is the route for a reference nobody wants held in memory. It answers
    the site question exactly -- a canonical table and a scan both count every
    occurrence of a k-mer and of its reverse complement -- and it cannot answer
    a positional question at all, so coverage and the gap figures carry that as
    their reason rather than a zero.
    """
    from neoswga.core import kmer_tables
    from neoswga.core.thermodynamics import reverse_complement

    no_positions = (
        "sites were counted rather than located, so no positional figure exists "
        "for this reference; scan it to measure coverage"
    )

    by_k: dict[int, list[str]] = {}
    for primer in panel:
        by_k.setdefault(len(primer), []).append(primer)

    per_primer: dict[str, int | None] = {}
    for k, group in sorted(by_k.items()):
        # Both spellings, because a canonical table holds the pair under one of
        # them and answers 0 for the other, and that 0 is not a measurement.
        # docs/validation/query_scan_2026-09-25.md measures the difference.
        wanted = sorted({*group, *(reverse_complement(p) for p in group)})
        try:
            counts = kmer_tables.counts_for(spec.prefix, k, wanted, genome=spec.genome)
        except (FileNotFoundError, OSError, RuntimeError, ValueError) as exc:
            logger.debug("Counts unavailable for %s at k=%d: %s", spec.prefix, k, exc)
            return _unavailable_record(spec, f"no {k}-mer counts available: {exc}")
        for primer in group:
            forward = int(counts.get(primer, 0) or 0)
            reverse = int(counts.get(reverse_complement(primer), 0) or 0)
            per_primer[primer] = max(forward, reverse)

    sites = sum(int(value or 0) for value in per_primer.values())
    return ReferenceRecord(
        prefix=spec.prefix,
        genome=spec.genome,
        role=spec.role,
        length=int(spec.length),
        status=STATUS_MEASURED,
        source=SOURCE_COUNTS,
        sites=Measurement(
            "sites", float(sites), "sites", basis="k-mer counts, canonical, both strands"
        ),
        site_density=_density(sites, spec.length),
        coverage=Measurement("coverage", None, "fraction", unavailable=no_positions),
        mean_gap=Measurement("mean_gap", None, "bp", unavailable=no_positions),
        max_gap=Measurement("max_gap", None, "bp", unavailable=no_positions),
        gap_gini=Measurement("gap_gini", None, "fraction", unavailable=no_positions),
        weighted_load=Measurement(
            "weighted_site_load", None, "sites", unavailable="not attempted yet"
        ),
        per_primer_sites=per_primer,
        per_primer_status={p: ("no_sites" if not per_primer[p] else "ok") for p in panel},
    )


# ----------------------------------------------------------------------
# The weighted load, the cross table and the reductions
# ----------------------------------------------------------------------


def _with_weighted_load(
    record: ReferenceRecord, panel: Sequence[str], conditions: Any, max_mismatches: int
) -> ReferenceRecord:
    """Occupancy-weighted load for one reference, when its tables allow it."""
    import dataclasses

    if record.status == STATUS_UNAVAILABLE:
        return dataclasses.replace(
            record,
            weighted_load=Measurement(
                "weighted_site_load",
                None,
                "sites",
                unavailable="this reference was not measured at all",
            ),
        )

    if conditions is None:
        reason = "no reaction conditions were supplied, so occupancy is undefined"
        return dataclasses.replace(
            record,
            weighted_load=Measurement("weighted_site_load", None, "sites", unavailable=reason),
        )

    from neoswga.core.occupancy import weighted_site_load

    try:
        load = weighted_site_load(list(panel), [record.prefix], conditions, max_mismatches)
    except (FileNotFoundError, OSError, RuntimeError, ValueError, KeyError) as exc:
        # Absent rather than substituted by the exact count: the two are
        # different numbers, and on a set with no exact host matches they
        # disagree completely.
        logger.debug("Weighted load unavailable for %s: %s", record.prefix, exc)
        return dataclasses.replace(
            record,
            weighted_load=Measurement(
                "weighted_site_load",
                None,
                "sites",
                unavailable=f"mismatch-class counts unavailable: {exc}",
            ),
        )
    return dataclasses.replace(
        record,
        weighted_load=Measurement(
            "weighted_site_load",
            float(load),
            "sites",
            basis=f"<= {max_mismatches} mismatch(es) at {getattr(conditions, 'temp', '?')} C",
        ),
    )


def _cross_table(records: Sequence[ReferenceRecord]) -> tuple[PairRecord, ...]:
    """Every (target, host) pair, each against its own two lengths."""
    targets = [r for r in records if r.role == ROLE_TARGET]
    hosts = [r for r in records if r.role == ROLE_HOST]
    pairs: list[PairRecord] = []
    for target in targets:
        for host in hosts:
            pairs.append(_pair(target, host))
    return tuple(pairs)


def _pair(target: ReferenceRecord, host: ReferenceRecord) -> PairRecord:
    reason = _pair_unavailable(target, host)
    if reason:
        return PairRecord(
            target=target.prefix,
            host=host.prefix,
            density=Measurement("selectivity_density", None, "ratio", unavailable=reason),
            ratio=Measurement("selectivity_ratio", None, "ratio", unavailable=reason),
        )

    fg = float(target.sites.value or 0.0)
    bg = float(host.sites.value or 0.0)
    density = selectivity_density_from_loads(fg, float(target.length), bg, float(host.length))
    ratio = selectivity_from_loads(fg, bg)
    basis = (
        f"{int(fg)} target site(s) per {target.length} bp against "
        f"{int(bg)} host site(s) per {host.length} bp"
    )
    note = " (unbounded: no host sites detected)" if bg <= 0 else ""
    return PairRecord(
        target=target.prefix,
        host=host.prefix,
        density=Measurement("selectivity_density", float(density), "ratio", basis=basis + note),
        ratio=Measurement(
            "selectivity_ratio",
            float(ratio),
            "ratio",
            basis=basis + note + "; moves with host size, so not comparable across hosts",
        ),
        zero_host_sites=bg <= 0,
    )


def _pair_unavailable(target: ReferenceRecord, host: ReferenceRecord) -> str:
    for record, role in ((target, "target"), (host, "host")):
        if record.sites.value is None:
            return f"{role} {record.prefix} has no measured site count: {record.sites.unavailable}"
        if record.length <= 0:
            return f"{role} {record.prefix} has no length, so a density cannot be formed"
    return ""


def _worst_coverage(records: Sequence[ReferenceRecord]) -> Measurement:
    """The lowest target coverage, or unavailable if any target is.

    The reduction is the point of the phase. A worst case computed over the
    members that happened to answer is not a worst case, and it reads as one.
    """
    targets = [r for r in records if r.role == ROLE_TARGET]
    if not targets:
        return Measurement(
            "worst_target_coverage", None, "fraction", unavailable="no target reference was given"
        )
    missing = [r for r in targets if r.coverage.value is None]
    if missing:
        names = ", ".join(r.prefix for r in missing)
        return Measurement(
            "worst_target_coverage",
            None,
            "fraction",
            unavailable=(
                f"coverage is unavailable for {len(missing)} of {len(targets)} targets "
                f"({names}), so the worst over them is unknown rather than the lowest "
                f"of the rest"
            ),
        )
    worst = min(targets, key=lambda r: r.coverage.value)
    return Measurement(
        "worst_target_coverage",
        float(worst.coverage.value),
        "fraction",
        basis=f"lowest of {len(targets)} target(s); {worst.prefix}",
    )


def _worst_density(pairs: Sequence[PairRecord]) -> Measurement:
    """The lowest target-against-host density, or unavailable if any pair is.

    A pair with no host sites carries `MAX_SELECTIVITY`, which is larger than
    any measured value, so including it in the minimum cannot hide a real one.
    """
    if not pairs:
        return Measurement(
            "worst_host_selectivity_density",
            None,
            "ratio",
            unavailable="no target-host pair was formed",
        )
    missing = [pair for pair in pairs if pair.density.value is None]
    if missing:
        names = ", ".join(f"{pair.target}/{pair.host}" for pair in missing)
        return Measurement(
            "worst_host_selectivity_density",
            None,
            "ratio",
            unavailable=(
                f"{len(missing)} of {len(pairs)} pair(s) could not be measured ({names}), "
                f"so the worst over them is unknown rather than the lowest of the rest"
            ),
        )
    worst = min(pairs, key=lambda pair: pair.density.value)
    suffix = " (every host carried no sites)" if all(p.zero_host_sites for p in pairs) else ""
    return Measurement(
        "worst_host_selectivity_density",
        float(worst.density.value),
        "ratio",
        basis=f"lowest of {len(pairs)} pair(s); {worst.target} against {worst.host}{suffix}",
    )


def _notes(records: Sequence[ReferenceRecord]) -> tuple[str, ...]:
    notes: list[str] = []
    counted = [r.prefix for r in records if r.source == SOURCE_COUNTS]
    if counted:
        notes.append(
            "Sites were counted rather than located for "
            + ", ".join(counted)
            + "; coverage and the gap figures are unavailable for those references."
        )
    unavailable = [r.prefix for r in records if r.status == STATUS_UNAVAILABLE]
    if unavailable:
        notes.append(
            "No measurement at all for " + ", ".join(unavailable) + "; see each record's reason."
        )
    return tuple(notes)


# ----------------------------------------------------------------------
# Small shared pieces
# ----------------------------------------------------------------------


def _unique(primers: Sequence[str]) -> list[str]:
    out: list[str] = []
    for primer in primers:
        sequence = str(primer).strip().upper()
        if sequence and sequence not in out:
            out.append(sequence)
    return out


def _density(sites: int, length: int) -> Measurement:
    if length <= 0:
        return Measurement(
            "sites_per_mb",
            None,
            "sites/Mb",
            unavailable="the reference length is unknown, so a density cannot be formed",
        )
    return Measurement(
        "sites_per_mb",
        float(sites) * 1e6 / float(length),
        "sites/Mb",
        basis=f"{sites} site(s) over {length} bp",
    )


def _unavailable_record(
    spec: ReferenceSpec,
    reason: str,
    *,
    per_primer_sites: Mapping[str, int | None] | None = None,
    per_primer_status: Mapping[str, str] | None = None,
    source: str = SOURCE_POSITIONS,
) -> ReferenceRecord:
    """A reference that was not measured, saying why, with no substituted zero."""
    return ReferenceRecord(
        prefix=spec.prefix,
        genome=spec.genome,
        role=spec.role,
        length=int(spec.length),
        status=STATUS_UNAVAILABLE,
        source=source,
        sites=Measurement("sites", None, "sites", unavailable=reason),
        site_density=Measurement("sites_per_mb", None, "sites/Mb", unavailable=reason),
        coverage=Measurement("coverage", None, "fraction", unavailable=reason),
        mean_gap=Measurement("mean_gap", None, "bp", unavailable=reason),
        max_gap=Measurement("max_gap", None, "bp", unavailable=reason),
        gap_gini=Measurement("gap_gini", None, "fraction", unavailable=reason),
        weighted_load=Measurement(
            "weighted_site_load", None, "sites", unavailable="not attempted yet"
        ),
        per_primer_sites=dict(per_primer_sites or {}),
        per_primer_status=dict(per_primer_status or {}),
        unavailable=reason,
    )
