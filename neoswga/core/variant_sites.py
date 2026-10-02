"""Which binding sites of a panel survive in which strain, and what that is not.

A binding site on the reference occupies `[pos, pos + k)`, which is the
exact-match convention all three scanners in this package follow. A site is
INTACT in a strain when no variant that strain carries falls inside that
interval, and AFFECTED otherwise. That reading is binary and needs no model,
which is why it is the one reported: the shipped mismatch model is a uniform
4.0 C per mismatch with status `assumed` in `core/registry/model_evidence.json`,
and it says nothing about where in the site the mismatch sits.

Affected sites are additionally split by the distance of the nearest variant
from the primer's 3' end, on the strand the primer binds, because a mismatch
near the 3' terminus is expected to matter more than a distal one. The split is
REPORTED and nothing is computed from it. No affected site is scored, weighted,
or called tolerated.

Three limits of this route, which `evaluate-set` writes into its output:

1. Sites GAINED in a strain through a variant are invisible. The reference is
   where sites are found, so a k-mer that a variant creates in a strain is
   never looked for.
2. Indels are treated as affecting any site they overlap, and the coordinate
   shift they cause downstream is ignored. Positions in a carrying strain are
   therefore reference positions, and gap lengths in that strain are
   approximate.
3. The table says nothing about sequence absent from the reference. A region
   present in a strain and not in the reference holds no sites here, which is
   not a measurement that it holds none.

Every figure is a `panel_evaluation.Measurement`, and a reduction over a set
with one unavailable strain is unavailable rather than the worst of the rest.
"""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.panel_evaluation import Measurement
from neoswga.core.position_cache import POSITION_DTYPE
from neoswga.core.variant_table import StrainVariants, VariantTable

logger = logging.getLogger(__name__)

__all__ = [
    "LIMITS",
    "THREE_PRIME_WINDOW_NT",
    "IntactPositions",
    "SiteGeometry",
    "StrainSites",
    "classify_sites",
    "evaluate_strain_panel",
    "intact_mask",
    "nearest_variant_distances",
    "split_sites_by_strand",
    "three_prime_offset",
]

#: How close to the 3' end a variant counts as proximal, in bases. This is a
#: REPORTING split, not a threshold anything is compared against: no penalty,
#: weight or rejection follows from which side of it a site falls. Phase 4b of
#: the diversity plan is where a position-dependent model, and a width with a
#: measurement behind it, would come from.
THREE_PRIME_WINDOW_NT = 5

#: The three limits of the variant route, written into every output that uses
#: it so a reader never has the per-strain figures without them.
LIMITS = (
    "Sites gained in a strain through a variant are invisible to this route: "
    "sites are found on the reference, so a binding site a variant creates is "
    "never looked for.",
    "An indel is treated as affecting any site it overlaps, and the coordinate "
    "shift it causes downstream is ignored, so gap lengths in a strain "
    "carrying one are approximate.",
    "The table says nothing about sequence absent from the reference. A region "
    "a strain has and the reference does not holds no sites here, which is not "
    "a measurement that it holds none.",
)

#: Strand names are `PositionCache`'s: 'forward', 'reverse' or 'both'. Not
#: '+'/'-' (Known Issue 4).
_FORWARD = "forward"
_REVERSE = "reverse"


# ----------------------------------------------------------------------
# Site geometry and the mask
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class SiteGeometry:
    """How the scanner laid a reference out, which decides what a site covers.

    `length` is the concatenated length the scanner wraps on and `circular`
    whether it wrapped. Both scanners in `string_search` treat a circular
    reference as ONE ring of the concatenated sequence, not one ring per
    record: they search `sequence + sequence[:k - 1]` and keep a match whose
    start lies in `[0, length)`, so a site starting within `k - 1` of the end
    covers `[pos, length)` and then `[0, pos + k - length)`. On a multi-record
    reference that ring joins the last record to the first, and
    `spans_a_record_join` does not stop it, because no record boundary lies
    inside such a match. This module follows that definition exactly rather
    than a per-record one the scanner does not use.
    """

    length: int
    circular: bool = False


def _site_bases(sites: np.ndarray, k: int, geometry: SiteGeometry) -> tuple[np.ndarray, np.ndarray]:
    """Each site's k reference bases, and whether its geometry is certain.

    Row `i` holds the concatenated coordinates of site `i`'s bases in the order
    the k-mer is read on the forward strand, so column `j` is the k-mer's base
    `j`. A site is ASSESSABLE when the scanner could have produced it under
    `geometry`: it starts inside the reference and either fits before the end,
    or the reference is circular and the k-mer is no longer than the ring. A
    site that fits neither -- a wrap position in a linear run, an offset past
    the end, a k-mer longer than a circular reference -- has no definition this
    module can follow with certainty, so it is not assessed. It is never
    reported intact.
    """
    width = int(k)
    length = int(geometry.length)
    offsets = np.arange(width, dtype=POSITION_DTYPE)
    bases = sites[:, None] + offsets[None, :]
    inside = (sites >= 0) & (sites < length)
    fits = sites + width <= length
    wraps = inside & ~fits & bool(geometry.circular) & (width <= length)
    assessable = inside & (fits | wraps)
    if np.any(wraps):
        bases = np.where(wraps[:, None], bases % length, bases)
    return bases, assessable


def _covered(bases: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Whether each base lies inside any variant interval.

    `starts` sorted ascending, `ends` in the same order. A base `b` is covered
    when some interval begins at or before it and the LARGEST end among those
    intervals reaches past it. The running maximum is what makes this exact
    when a long interval begins before a short one: the short one's end says
    nothing about the bases the long one still covers.
    """
    if starts.size == 0:
        return np.zeros(bases.shape, dtype=bool)
    running_max_end = np.maximum.accumulate(ends)
    before = np.searchsorted(starts, bases, side="right")
    covered = np.zeros(bases.shape, dtype=bool)
    nonzero = before > 0
    covered[nonzero] = running_max_end[before[nonzero] - 1] > bases[nonzero]
    return covered


def _as_intervals(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    starts = np.asarray(starts, dtype=POSITION_DTYPE)
    ends = np.asarray(ends, dtype=POSITION_DTYPE)
    if starts.size != ends.size:
        raise ValueError(
            f"the variant interval arrays disagree: {starts.size} start(s) "
            f"against {ends.size} end(s)"
        )
    return starts, ends


def classify_sites(
    positions: np.ndarray,
    k: int,
    starts: np.ndarray,
    ends: np.ndarray,
    geometry: SiteGeometry,
) -> tuple[np.ndarray, np.ndarray]:
    """`(intact, assessable)` boolean arrays, one entry per site.

    A site is intact when it is assessable and no variant interval covers any
    of its k bases; it is affected when it is assessable and one does. A site
    that is not assessable is neither, and `intact` is False for it, so it can
    never be counted as intact by accident.
    """
    sites = np.asarray(positions, dtype=POSITION_DTYPE)
    if int(k) <= 0:
        raise ValueError(f"a binding site cannot be {k} bases long")
    if sites.size == 0:
        return np.zeros(0, dtype=bool), np.zeros(0, dtype=bool)
    starts, ends = _as_intervals(starts, ends)
    bases, assessable = _site_bases(sites, k, geometry)
    hit = _covered(bases, starts, ends).any(axis=1)
    return assessable & ~hit, assessable


def intact_mask(
    positions: np.ndarray,
    k: int,
    starts: np.ndarray,
    ends: np.ndarray,
    geometry: SiteGeometry,
) -> np.ndarray:
    """True where the site is assessable and no variant falls inside it.

    `starts` must be sorted ascending, which `variant_table` guarantees. The
    test is a `searchsorted` per site base plus a running maximum of the
    interval ends (see `_covered`). An empty variant set leaves every
    assessable site intact, which is a measurement: the strain carries no
    variants. A strain whose variants could not be READ never reaches here --
    see `StrainVariants.status`.
    """
    intact, _assessable = classify_sites(positions, k, starts, ends, geometry)
    return intact


def three_prime_offset(position: int, k: int, strand: str, geometry: SiteGeometry) -> int:
    """Where the primer's 3' terminal base sits, in concatenated coordinates.

    A site is stored as the forward-strand start offset of the k-mer, for
    either strand. A primer matching the forward strand runs 5'->3' left to
    right, so its 3' base is the k-mer's last base, `position + k - 1`; on a
    circular reference that is taken modulo the length the scanner wraps on. A
    primer matching the reverse strand is stored under its reverse complement's
    key at the same forward-strand offset and anneals antiparallel, so its 3'
    base is at `position`.
    """
    if strand == _FORWARD:
        anchor = int(position) + int(k) - 1
        if geometry.circular and anchor >= int(geometry.length):
            anchor -= int(geometry.length)
        return anchor
    if strand == _REVERSE:
        return int(position)
    raise ValueError(
        f"strand must be {_FORWARD!r} or {_REVERSE!r}, got {strand!r}; "
        f"PositionCache does not use '+'/'-' (Known Issue 4)"
    )


def nearest_variant_distances(
    positions: np.ndarray,
    k: int,
    strand: str,
    starts: np.ndarray,
    ends: np.ndarray,
    geometry: SiteGeometry,
) -> np.ndarray:
    """Bases from the primer's 3' end to the nearest variant inside each site.

    Entry `i` is -1 where no variant falls inside site `i`, or where the site
    is not assessable; the mask, not this, says whether a site is affected.
    Distances are counted along the primer, so 0 is the 3' terminal base: base
    `j` of the forward-strand k-mer is `k - 1 - j` bases from the 3' end of a
    forward-strand primer and `j` bases from that of a reverse-strand one. That
    holds for a wrapped site too, because it is counted along the k-mer rather
    than along the coordinates.
    """
    if strand not in (_FORWARD, _REVERSE):
        three_prime_offset(0, k, strand, geometry)  # raises the strand error
    sites = np.asarray(positions, dtype=POSITION_DTYPE)
    out = np.full(sites.size, -1, dtype=POSITION_DTYPE)
    starts, ends = _as_intervals(starts, ends)
    if sites.size == 0 or starts.size == 0:
        return out

    width = int(k)
    bases, assessable = _site_bases(sites, width, geometry)
    covered = _covered(bases, starts, ends) & assessable[:, None]
    offsets = np.arange(width, dtype=POSITION_DTYPE)
    distance = (width - 1 - offsets) if strand == _FORWARD else offsets
    sentinel = np.iinfo(POSITION_DTYPE).max
    per_base = np.where(covered, distance[None, :], sentinel)
    nearest = per_base.min(axis=1)
    hit = nearest != sentinel
    out[hit] = nearest[hit]
    return out


def split_sites_by_strand(cache: Any, prefix: str, primer: str) -> dict[str, np.ndarray]:
    """One primer's sites on one reference, kept per strand.

    Per strand because the 3'-end distance depends on it: the same stored
    offset means a 3' terminus at the right-hand end of the site on the forward
    strand and at the left-hand end on the reverse strand.
    """
    return {
        _FORWARD: np.asarray(cache.get_positions(prefix, primer, _FORWARD), dtype=POSITION_DTYPE),
        _REVERSE: np.asarray(cache.get_positions(prefix, primer, _REVERSE), dtype=POSITION_DTYPE),
    }


# ----------------------------------------------------------------------
# A cache-shaped view over the intact sites only
# ----------------------------------------------------------------------


class IntactPositions:
    """The sites intact in one strain, answering as a `PositionCache` does.

    This exists so the per-strain coverage and gap figures come out of the SAME
    functions the reference figures do (`coverage.compute_per_prefix_coverage`
    and `reference_panel_evaluation.gap_statistics`) rather than a second
    implementation beside them. It answers only for the prefix and the primers
    it was built over, and raises for anything else, because an empty array for
    a question nobody asked is the silent zero this package keeps meeting.
    """

    def __init__(self, cache: Any, prefix: str, sites: Mapping[str, Mapping[str, np.ndarray]]):
        self._cache = cache
        self._prefix = prefix
        self._sites = {primer: dict(per_strand) for primer, per_strand in sites.items()}

    def get_positions(self, fname_prefix: str, primer: str, strand: str = "both") -> np.ndarray:
        if fname_prefix != self._prefix:
            raise KeyError(
                f"{fname_prefix!r} is not the reference this intact-site view "
                f"was built over ({self._prefix!r})"
            )
        per_strand = self._sites.get(primer)
        if per_strand is None:
            raise KeyError(f"no intact sites were computed for {primer!r}")
        if strand == "both":
            return np.unique(np.concatenate([per_strand[_FORWARD], per_strand[_REVERSE]]))
        if strand in (_FORWARD, _REVERSE):
            return per_strand[strand]
        raise ValueError(f"strand must be {_FORWARD!r}, {_REVERSE!r} or 'both', got {strand!r}")

    def get_record_starts(self, fname_prefix: str) -> list[int]:
        """Delegated: record boundaries are the reference's, not a strain's."""
        getter = getattr(self._cache, "get_record_starts", None)
        if getter is None:
            return []
        return list(getter(fname_prefix) or [])


# ----------------------------------------------------------------------
# Per strain
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class StrainSites:
    """One strain's intact and affected sites, and the figures over them.

    `intact_sites + affected_sites + not_assessed_sites` is the panel's
    reference site count, deduplicated over the two strands, for every
    measured strain.
    """

    name: str
    status: str
    intact_sites: Measurement
    affected_sites: Measurement
    not_assessed_sites: Measurement
    proximal: Measurement
    distal: Measurement
    coverage: Measurement
    mean_gap: Measurement
    max_gap: Measurement
    gap_gini: Measurement
    intact_fraction: Measurement
    variants: int | None = None
    per_primer_intact: Mapping[str, float | None] = field(default_factory=dict)
    unavailable: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "strain": self.name,
            "status": self.status,
            "unavailable": self.unavailable or None,
            "variants": self.variants,
            "intact_sites": self.intact_sites.as_dict(),
            "affected_sites": self.affected_sites.as_dict(),
            "not_assessed_sites": self.not_assessed_sites.as_dict(),
            "affected_three_prime_proximal": self.proximal.as_dict(),
            "affected_distal": self.distal.as_dict(),
            "three_prime_window_nt": THREE_PRIME_WINDOW_NT,
            "coverage_on_intact_sites": self.coverage.as_dict(),
            "mean_gap_bp": self.mean_gap.as_dict(),
            "max_gap_bp": self.max_gap.as_dict(),
            "gap_gini": self.gap_gini.as_dict(),
            "intact_site_fraction": self.intact_fraction.as_dict(),
            "per_primer_intact_fraction": dict(self.per_primer_intact or {}),
        }


def _union(per_strand: Mapping[str, np.ndarray]) -> np.ndarray:
    """One primer's sites over both strands, each site once.

    A palindromic primer's site is stored under both strand keys at the same
    offset and is ONE site. Every count in this module -- the reference total
    as well as the intact and affected counts -- is taken over this union, so
    the figure agrees with `reference_panel_evaluation`'s `sites`.
    """
    return np.unique(
        np.concatenate(
            [
                np.asarray(per_strand[_FORWARD], dtype=POSITION_DTYPE),
                np.asarray(per_strand[_REVERSE], dtype=POSITION_DTYPE),
            ]
        )
    )


def evaluate_strain_panel(
    cache: Any,
    prefix: str,
    primers: Sequence[str],
    length: int,
    table: VariantTable,
    *,
    extension: int = 3000,
    circular: bool = False,
    window: int = THREE_PRIME_WINDOW_NT,
) -> dict[str, Any]:
    """The `per_strain` block: intact sites, their coverage, and the reductions.

    The reference figures are computed through the same two functions over the
    same view, so a strain carrying no variants reproduces them exactly rather
    than approximately.

    Args:
        cache: anything answering `get_positions(prefix, primer, strand)` for
            this reference -- normally the `PositionCache` the caller already
            built, so no reference is scanned twice.
        prefix: the one reference the variant table was read against.
        primers: the panel.
        length: that reference's length in bases, its own denominator. It must
            equal the length of the FASTA the table was placed against, since
            that is the length a circular scan wraps on.
        table: the `VariantTable` from `variant_table.open_variants`.
        extension: per-primer reach, resolved by the caller.
        circular: as the caller scanned and treats this reference.
        window: the 3'-proximal reporting width, in bases.

    Returns:
        A JSON-ready dict. Every figure is a measurement with a reason when it
        is absent, and the two reductions are None when any strain is
        unavailable.

    Raises:
        ReferenceDataError: `length` disagrees with the FASTA's own length, so
            neither the wrap point nor any coverage denominator is known.
    """
    if int(length) != int(table.layout.total_length):
        raise ReferenceDataError(
            f"reference {table.reference}",
            f"the run gives this reference a length of {int(length)} bp and the "
            f"FASTA the variant table was placed against holds "
            f"{int(table.layout.total_length)} bp",
            "the configured length and the FASTA must describe the same "
            "sequence; correct fg_seq_lengths in params.json, or re-run filter "
            "if the FASTA has changed",
        )
    geometry = SiteGeometry(length=int(table.layout.total_length), circular=bool(circular))
    panel = [str(primer).strip().upper() for primer in primers if str(primer).strip()]
    reference_sites = {primer: split_sites_by_strand(cache, prefix, primer) for primer in panel}
    reference_total = sum(int(_union(per_strand).size) for per_strand in reference_sites.values())

    reference_block = {
        name: figure.as_dict()
        for name, figure in _figures(
            cache, prefix, panel, length, reference_sites, extension, circular
        ).items()
    }
    strains = [
        _strain_record(
            cache,
            prefix,
            panel,
            length,
            reference_sites,
            reference_total,
            strain,
            geometry,
            extension,
            window,
        )
        for strain in table.strains
    ]

    return {
        "source": table.as_dict(),
        "reference_fasta": table.reference,
        "reference_prefix": prefix,
        "geometry": {"length_bp": geometry.length, "circular": geometry.circular},
        "three_prime_window_nt": int(window),
        "three_prime_window_basis": (
            "a reporting split only: no penalty, weight or rejection follows "
            "from which side of it an affected site falls"
        ),
        "reference_sites": reference_total,
        "reference": reference_block,
        "per_strain": {record.name: record.as_dict() for record in strains},
        "worst_strain_coverage": _worst(
            strains, "coverage", "worst_strain_coverage", "fraction"
        ).as_dict(),
        "worst_strain_intact_fraction": _worst(
            strains, "intact_fraction", "worst_strain_intact_fraction", "fraction"
        ).as_dict(),
        "per_primer_intact_fraction_every_strain": _per_primer_every_strain(strains, panel),
        "limits": list(LIMITS),
    }


@dataclass
class _Tally:
    intact: int = 0
    affected: int = 0
    not_assessed: int = 0
    proximal: int = 0
    distal: int = 0


def _tally_primer(
    primer: str,
    per_strand: Mapping[str, np.ndarray],
    strain: StrainVariants,
    geometry: SiteGeometry,
    window: int,
    tally: _Tally,
) -> tuple[dict[str, np.ndarray], float | None]:
    """One primer's contribution, and the intact sites per strand for the view.

    Returns the per-strand intact arrays and the primer's intact fraction,
    which is None when the primer has no site, or has a site that could not be
    assessed: in neither case is there a fraction of its sites to report.
    """
    k = len(primer)
    union = _union(per_strand)
    intact, assessable = classify_sites(union, k, strain.starts, strain.ends, geometry)
    n_intact = int(np.count_nonzero(intact))
    n_assessable = int(np.count_nonzero(assessable))
    tally.intact += n_intact
    tally.affected += n_assessable - n_intact
    tally.not_assessed += int(union.size) - n_assessable

    affected = union[assessable & ~intact]
    if affected.size:
        nearest = _nearest_over_strands(affected, k, per_strand, strain, geometry)
        measured = nearest[nearest >= 0]
        tally.proximal += int(np.count_nonzero(measured < window))
        tally.distal += int(np.count_nonzero(measured >= window))

    per_strand_intact = {}
    for name in (_FORWARD, _REVERSE):
        sites = np.asarray(per_strand[name], dtype=POSITION_DTYPE)
        per_strand_intact[name] = sites[intact_mask(sites, k, strain.starts, strain.ends, geometry)]

    fraction = None
    if union.size and n_assessable == union.size:
        fraction = n_intact / int(union.size)
    return per_strand_intact, fraction


def _strain_record(
    cache: Any,
    prefix: str,
    panel: Sequence[str],
    length: int,
    reference_sites: Mapping[str, Mapping[str, np.ndarray]],
    reference_total: int,
    strain: StrainVariants,
    geometry: SiteGeometry,
    extension: int,
    window: int,
) -> StrainSites:
    if not strain.measured:
        return _unavailable_strain(strain, panel)

    tally = _Tally()
    intact: dict[str, dict[str, np.ndarray]] = {}
    per_primer: dict[str, float | None] = {}
    for primer in panel:
        intact[primer], per_primer[primer] = _tally_primer(
            primer, reference_sites[primer], strain, geometry, window, tally
        )

    if tally.not_assessed:
        # Every figure over the intact sites would silently leave out the
        # sites nobody could assess, and read as a figure over all of them.
        reason = (
            f"{tally.not_assessed} reference site(s) have a geometry the scanner "
            f"definition does not cover at length {geometry.length} "
            f"({'circular' if geometry.circular else 'linear'}), so they were "
            f"counted neither intact nor affected; a figure over the rest would "
            f"read as a figure over all of them"
        )
        figures = {
            "coverage": Measurement("coverage", None, "fraction", unavailable=reason),
            "mean_gap": Measurement("mean_gap", None, "bp", unavailable=reason),
            "max_gap": Measurement("max_gap", None, "bp", unavailable=reason),
            "gap_gini": Measurement("gap_gini", None, "fraction", unavailable=reason),
        }
        fraction = Measurement("intact_site_fraction", None, "fraction", unavailable=reason)
    else:
        figures = _figures(cache, prefix, panel, length, intact, extension, geometry.circular)
        fraction = _intact_fraction(tally.intact, reference_total)

    return StrainSites(
        name=strain.name,
        status="measured",
        variants=strain.count,
        intact_sites=Measurement(
            "intact_sites",
            float(tally.intact),
            "sites",
            basis="reference sites with no variant inside them, both strands, each once",
        ),
        affected_sites=Measurement(
            "affected_sites",
            float(tally.affected),
            "sites",
            basis="reference sites with at least one variant inside them",
        ),
        not_assessed_sites=Measurement(
            "not_assessed_sites",
            float(tally.not_assessed),
            "sites",
            basis="sites whose geometry this route cannot place with certainty; "
            "excluded from intact and affected alike",
        ),
        proximal=Measurement(
            "affected_three_prime_proximal",
            float(tally.proximal),
            "sites",
            basis=f"nearest variant within {int(window)} base(s) of the 3' end",
        ),
        distal=Measurement(
            "affected_distal",
            float(tally.distal),
            "sites",
            basis=f"nearest variant {int(window)} or more bases from the 3' end",
        ),
        coverage=figures["coverage"],
        mean_gap=figures["mean_gap"],
        max_gap=figures["max_gap"],
        gap_gini=figures["gap_gini"],
        intact_fraction=fraction,
        per_primer_intact=per_primer,
    )


def _intact_fraction(intact: int, reference_total: int) -> Measurement:
    if not reference_total:
        return Measurement(
            "intact_site_fraction",
            None,
            "fraction",
            unavailable="the panel has no binding site on the reference, so no "
            "fraction of its sites can be formed",
        )
    return Measurement(
        "intact_site_fraction",
        intact / reference_total,
        "fraction",
        basis=f"{intact} of {reference_total} reference site(s)",
    )


def _nearest_over_strands(
    affected: np.ndarray,
    k: int,
    per_strand: Mapping[str, np.ndarray],
    strain: StrainVariants,
    geometry: SiteGeometry,
) -> np.ndarray:
    """For each affected site, the shortest 3'-end distance over its strands.

    A site occurs on one strand for almost every primer, and on both for a
    palindromic one. Taking the shorter of the two distances keeps the site
    counted once and reports the stricter reading, which is the 3'-proximal
    one.
    """
    best = np.full(affected.size, -1, dtype=POSITION_DTYPE)
    for name in (_FORWARD, _REVERSE):
        sites = np.asarray(per_strand[name], dtype=POSITION_DTYPE)
        if sites.size == 0:
            continue
        on_strand = np.isin(affected, sites)
        if not np.any(on_strand):
            continue
        distances = nearest_variant_distances(
            affected[on_strand], k, name, strain.starts, strain.ends, geometry
        )
        current = best[on_strand]
        improved = np.where(
            (distances >= 0) & ((current < 0) | (distances < current)), distances, current
        )
        best[on_strand] = improved
    return best


def _unavailable_strain(strain: StrainVariants, panel: Sequence[str]) -> StrainSites:
    reason = strain.unavailable or "this strain's variants could not be read"
    return StrainSites(
        name=strain.name,
        status="unavailable",
        variants=None,
        intact_sites=Measurement("intact_sites", None, "sites", unavailable=reason),
        affected_sites=Measurement("affected_sites", None, "sites", unavailable=reason),
        not_assessed_sites=Measurement("not_assessed_sites", None, "sites", unavailable=reason),
        proximal=Measurement("affected_three_prime_proximal", None, "sites", unavailable=reason),
        distal=Measurement("affected_distal", None, "sites", unavailable=reason),
        coverage=Measurement("coverage", None, "fraction", unavailable=reason),
        mean_gap=Measurement("mean_gap", None, "bp", unavailable=reason),
        max_gap=Measurement("max_gap", None, "bp", unavailable=reason),
        gap_gini=Measurement("gap_gini", None, "fraction", unavailable=reason),
        intact_fraction=Measurement("intact_site_fraction", None, "fraction", unavailable=reason),
        per_primer_intact={primer: None for primer in panel},
        unavailable=reason,
    )


def _figures(
    cache: Any,
    prefix: str,
    panel: Sequence[str],
    length: int,
    sites: Mapping[str, Mapping[str, np.ndarray]],
    extension: int,
    circular: bool,
) -> dict[str, Measurement]:
    """Coverage and the three gap figures over one site set.

    Both come from the functions the per-reference evaluation uses, over a view
    that answers with this site set, so the reference and every strain are
    measured the one way.
    """
    from neoswga.core.coverage import compute_per_prefix_coverage
    from neoswga.core.reference_panel_evaluation import (
        ROLE_TARGET,
        ReferenceSpec,
        gap_statistics,
    )

    view = IntactPositions(cache, prefix, sites)
    _overall, per_prefix = compute_per_prefix_coverage(
        view,
        list(panel),
        [prefix],
        [int(length)],
        extension=int(extension),
        circular=bool(circular),
    )
    spec = ReferenceSpec(
        prefix=prefix,
        genome=None,
        length=int(length),
        role=ROLE_TARGET,
        circular=bool(circular),
        scan=False,
    )
    mean_gap, max_gap, gini = gap_statistics(view, list(panel), spec)
    if prefix in per_prefix:
        coverage = Measurement(
            "coverage",
            float(per_prefix[prefix]),
            "fraction",
            basis=f"union of {int(extension)} bp windows over the sites counted here, "
            f"{'circular' if circular else 'linear'}",
        )
    else:
        coverage = Measurement(
            "coverage",
            None,
            "fraction",
            unavailable=f"no coverage was returned for {prefix}, so it is unknown rather than zero",
        )
    return {
        "coverage": coverage,
        "mean_gap": mean_gap,
        "max_gap": max_gap,
        "gap_gini": gini,
    }


def _worst(strains: Sequence[StrainSites], attribute: str, name: str, units: str) -> Measurement:
    """The lowest value over the strains, or unavailable if any strain is.

    Not the lowest of the ones that answered. A worst case computed over the
    measured members reads as a worst case over all of them, which is how an
    unavailable strain comes to look like a passing one.
    """
    if not strains:
        return Measurement(name, None, units, unavailable="no strain was read from the table")
    missing = [strain for strain in strains if getattr(strain, attribute).value is None]
    if missing:
        names = ", ".join(strain.name for strain in missing)
        return Measurement(
            name,
            None,
            units,
            unavailable=(
                f"{len(missing)} of {len(strains)} strain(s) could not be measured "
                f"({names}), so the worst over them is unknown rather than the "
                f"lowest of the rest"
            ),
        )
    worst = min(strains, key=lambda strain: getattr(strain, attribute).value)
    return Measurement(
        name,
        float(getattr(worst, attribute).value),
        units,
        basis=f"lowest of {len(strains)} strain(s); {worst.name}",
    )


def _per_primer_every_strain(
    strains: Sequence[StrainSites], panel: Sequence[str]
) -> dict[str, float | None]:
    """Per primer, the lowest intact fraction across the strains.

    None for a primer whose fraction is unknown in any strain, and None for a
    primer with no site on the reference: in neither case is there a fraction
    to report, and 0.0 would read as a primer whose sites the variants
    destroyed.
    """
    out: dict[str, float | None] = {}
    for primer in panel:
        values = [(strain.per_primer_intact or {}).get(primer) for strain in strains]
        if not values or any(value is None for value in values):
            out[primer] = None
            continue
        out[primer] = min(float(value) for value in values if value is not None)
    return out
