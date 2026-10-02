"""Start from an oligo set somebody already has, and propose small edits to it.

Phase 2b of the genomic-diversity plan. `evaluate-set` says how a set does on
each target and each host; this says which oligo is responsible for what, and
which single edit would help. It proposes and never applies: nothing here
writes a primer set, and no selection stage reads it.

Four steps, in order.

1. **Evaluate** the set through `reference_panel_evaluation`, the same function
   `evaluate-set` calls. There is no second evaluation here.
2. **Attribute, per oligo**: its sites on every reference, its marginal
   contribution to each target's coverage (the set against the set without it),
   its dimer partners inside the set, and its Tm against the reaction window.
3. **Propose** single edits, reported in four sections: drops that cost no
   coverage, adds from the candidate pool, swaps for a candidate that could not
   be added beside the oligo it replaces, and trade-off drops.
4. **Check** every proposal before reporting it. A proposed add is screened
   against every RETAINED oligo at the configured `max_dimer_bp` and against the
   Tm window, and its sites must be known on every reference; one that fails any
   of these is not an option. The configured panel limits are advisory here: a
   proposal that misses one is reported with the limit named.

Rules this module keeps.

**The prediction is the measurement.** The figures reported for a proposal are
`evaluate_reference_panel` run on the edited set. Only the reading of positions
and counts is shared between panels (`MemoisedSources`), so applying a proposal
and running `evaluate-set` on the result, with the same params file, references
and geometry flags, gives the same numbers. Both commands resolve the reach
through `coverage.resolve_coverage_reach`. That makes the prediction a
statement about the model, not an independent check of it, and the report says
so.

**Unknown is not zero, and not success.** An oligo whose sites on a target
could not be established is reported as unavailable on that target. It is never
proposed for dropping on the grounds that it binds nothing, and while any
target is unmeasured no edit is proposed at all, because the worst target is
then unknown and every ranking below is by the worst target. A candidate whose
sites on any reference, target or host, could not be established is not
proposed either: it is counted as unmeasured with the reason.

**No threshold is invented.** A drop that costs nothing is a marginal coverage
of exactly zero on every target. Nothing here says when a host cost is too
high. An oligo whose removal raises the worst target-against-host selectivity
density is listed as a TRADE-OFF, with the density before and after and the
coverage it costs on each target. That list is not a list of improvements: the
member of any set with the lowest target-to-host ratio always qualifies.

**Every entry was evaluated alone.** A proposal is one edit to the set as it
stands. Two entries are not jointly applicable in general: two oligos that
cover the same bases each have zero marginal coverage, and dropping both was
not evaluated.

**A fixed oligo stays.** `fixed_oligos` on the design request marks an oligo
as already ordered or validated. It is attributed like any other and is never
proposed for dropping, alone or as half of a swap.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

from neoswga.core.panel_evaluation import Measurement
from neoswga.core.reference_panel_evaluation import (
    PanelSources,
    ReferencePanelAssessment,
    ReferenceRecord,
    ReferenceSpec,
    evaluate_reference_panel,
)

logger = logging.getLogger(__name__)

__all__ = [
    "KIND_ADD",
    "KIND_DROP",
    "KIND_SWAP",
    "KIND_TRADE_OFF",
    "SECTION_ORDER",
    "ImprovementReport",
    "ImprovementSettings",
    "Proposal",
    "ProposalSection",
    "improve_set",
    "report_lines",
]

KIND_DROP = "drop"
KIND_ADD = "add"
KIND_SWAP = "swap"
KIND_TRADE_OFF = "trade_off_drop"

#: The order sections are reported in. `max_edits` applies to each on its own,
#: so a pool full of helpful candidates cannot push the oligo that binds
#: nothing off the report.
SECTION_ORDER = (KIND_DROP, KIND_ADD, KIND_SWAP, KIND_TRADE_OFF)

_TITLES = {
    KIND_DROP: "Drops that cost no coverage",
    KIND_ADD: "Adds that raise the worst target",
    KIND_SWAP: "Swaps that raise the worst target",
    KIND_TRADE_OFF: "Trade-offs: drops that raise selectivity density and cost coverage",
}

_NOTES = {
    KIND_DROP: (
        "Each entry is true alone. Two oligos that cover the same bases each "
        "have zero marginal coverage, and dropping both was not evaluated."
    ),
    KIND_ADD: "Each entry is one add to the set as it stands; they were not evaluated together.",
    KIND_SWAP: "Each entry is one swap in the set as it stands; they were not evaluated together.",
    KIND_TRADE_OFF: (
        "These are not improvements. Removal raises the worst target-against-host "
        "selectivity density and lowers coverage, by the amounts listed. In any "
        "set the member with the lowest target-to-host ratio qualifies, so an "
        "entry here does not show that an oligo is unusually costly. Each entry "
        "was evaluated alone."
    ),
}

#: Facts about a section, carried as data so a reader (or a test) does not have
#: to parse the prose above.
_FACTS: dict[str, dict[str, bool]] = {
    KIND_DROP: {"each_entry_evaluated_alone": True, "jointly_applicable": False},
    KIND_ADD: {"each_entry_evaluated_alone": True, "jointly_applicable": False},
    KIND_SWAP: {"each_entry_evaluated_alone": True, "jointly_applicable": False},
    KIND_TRADE_OFF: {
        "each_entry_evaluated_alone": True,
        "jointly_applicable": False,
        "is_improvement": False,
        "lowest_ratio_member_always_qualifies": True,
    },
}

#: Said once in the report, because a reader comparing the predicted figure with
#: a later `evaluate-set` run would otherwise take the agreement as evidence.
PREDICTION_NOTE = (
    "Predicted by the same evaluation that would measure the edited set, so "
    "agreement with a later evaluate-set run on the same params file, references "
    "and geometry flags is expected and is not an independent check."
)


@dataclass(frozen=True)
class ImprovementSettings:
    """What a proposal is held to. Resolved by the caller, never defaulted here.

    `extension` is the per-primer reach coverage is measured at. The Tm window
    and the dimer limits are the configured ones; a proposed add must clear
    them against every oligo that stays. `max_edits` bounds each section of the
    report on its own.
    """

    extension: int
    max_dimer_bp: int
    min_tm: float
    max_tm: float
    max_dimer_dg: float | None = None
    fixed_oligos: tuple[str, ...] = ()
    excluded_oligos: tuple[str, ...] = ()
    max_edits: int = 5

    def as_dict(self) -> dict[str, Any]:
        return {
            "extension_reach_bp": int(self.extension),
            "max_dimer_bp": int(self.max_dimer_bp),
            "max_dimer_dg": self.max_dimer_dg,
            "min_tm": float(self.min_tm),
            "max_tm": float(self.max_tm),
            "fixed_oligos": list(self.fixed_oligos),
            "excluded_oligos": list(self.excluded_oligos),
            "max_edits_per_section": int(self.max_edits),
        }


class MemoisedSources(PanelSources):
    """Positions, counts and loads read once and shared across evaluations.

    `improve_set` evaluates the set, the set without each oligo and the set
    with each candidate. Read per call, a host given by FASTA would be scanned
    once per panel. Kept here, each reference is read once for every oligo and
    candidate together, and every panel is then evaluated by the unchanged
    evaluation code.

    Not part of the module's public interface: it exists for `improve_set`.
    Every key names the reference in full (prefix, genome path and length for
    counts; prefixes and genome paths for positions), so one reference's answer
    cannot be returned for another that shares a prefix or a k-mer.

    A position cache built over more primers answers a subset identically: the
    evaluation asks it per primer. A primer the held cache was not built over
    causes a rebuild over the union, never an answer from the cache that lacks
    it. A failed read is not kept, so a reference that could not be read is
    asked again and fails the same way.
    """

    def __init__(self) -> None:
        self._caches: dict[tuple[Any, ...], Any] = {}
        self._counts: dict[tuple[str, str | None, int, str], int] = {}
        self._loads: dict[tuple[str, tuple[str, ...], str, int], float] = {}

    def position_cache(self, prefixes, panel, genomes, circular, scan):
        key = (tuple(prefixes), tuple(genomes), bool(circular), bool(scan))
        cache = self._caches.get(key)
        if cache is not None and set(panel) <= set(cache.primers):
            return cache
        held = set(cache.primers) if cache is not None else set()
        # Rebuilt over the union through the ordinary constructor, so a primer
        # that arrives late is looked up by exactly the route the first ones
        # were. `improve_set` reads every oligo and candidate in one evaluation
        # first, so within it this runs once per reference group.
        union = sorted(held | set(panel))
        cache = super().position_cache(prefixes, union, genomes, circular, scan)
        self._caches[key] = cache
        return cache

    def counts(self, spec, k, kmers):
        reference = (spec.prefix, spec.genome, int(spec.length))
        wanted = [kmer for kmer in kmers if (*reference, kmer) not in self._counts]
        if wanted:
            fetched = super().counts(spec, k, wanted)
            for kmer in wanted:
                # A k-mer the table does not list was counted and occurs zero
                # times; that is the table's answer, not a gap in it.
                self._counts[(*reference, kmer)] = int(fetched.get(kmer, 0) or 0)
        return {kmer: self._counts[(*reference, kmer)] for kmer in kmers}

    def weighted_load(self, panel, prefix, conditions, max_mismatches):
        fingerprint = getattr(conditions, "fingerprint", None)
        if not callable(fingerprint):
            # Nothing stable identifies this reaction, so nothing is kept. An
            # object identity is not one: it can be reused after collection.
            return super().weighted_load(panel, prefix, conditions, max_mismatches)
        key = (prefix, tuple(panel), str(fingerprint()), int(max_mismatches))
        if key not in self._loads:
            self._loads[key] = super().weighted_load(panel, prefix, conditions, max_mismatches)
        return self._loads[key]


@dataclass(frozen=True)
class Proposal:
    """One edit, what it is predicted to do, and what it was checked against."""

    kind: str
    drop: tuple[str, ...]
    add: tuple[str, ...]
    reason: str
    resulting_set: tuple[str, ...]
    predicted: ReferencePanelAssessment
    worst_target_gain: float
    host_load: Measurement
    coverage_change: Mapping[str, Mapping[str, float | None]]
    density_before: Measurement
    density_after: Measurement
    depends_on: str = ""
    advisory: Mapping[str, Any] | None = None

    @property
    def oligos_changed(self) -> int:
        return len(self.drop) + len(self.add)

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "drop": list(self.drop),
            "add": list(self.add),
            "reason": self.reason,
            "depends_on": self.depends_on or None,
            "oligos_changed": self.oligos_changed,
            "resulting_set": list(self.resulting_set),
            "worst_target_coverage_gain": self.worst_target_gain,
            "per_target_coverage": {k: dict(v) for k, v in self.coverage_change.items()},
            "worst_host_selectivity_density": {
                "before": self.density_before.as_dict(),
                "after": self.density_after.as_dict(),
            },
            "worst_host_sites_per_mb": self.host_load.as_dict(),
            "predicted": self.predicted.as_dict(),
            "panel_limits": dict(self.advisory) if self.advisory is not None else None,
        }


@dataclass(frozen=True)
class ProposalSection:
    """One kind of edit: what is shown, of how many considered, and how to read it."""

    kind: str
    proposals: tuple[Proposal, ...]
    considered: int

    @property
    def title(self) -> str:
        return _TITLES[self.kind]

    @property
    def note(self) -> str:
        return _NOTES[self.kind]

    @property
    def facts(self) -> Mapping[str, bool]:
        return dict(_FACTS[self.kind])

    def as_dict(self) -> dict[str, Any]:
        return {
            "kind": self.kind,
            "title": self.title,
            "shown": len(self.proposals),
            "considered": int(self.considered),
            "note": self.note,
            **self.facts,
            "proposals": [
                dict(proposal.as_dict(), rank=rank)
                for rank, proposal in enumerate(self.proposals, start=1)
            ],
        }


@dataclass(frozen=True)
class ImprovementReport:
    """The current set, the diagnosis per oligo, and the edits by section."""

    primers: tuple[str, ...]
    settings: ImprovementSettings
    current: ReferencePanelAssessment
    current_host_load: Measurement
    attribution: tuple[Mapping[str, Any], ...]
    sections: tuple[ProposalSection, ...]
    candidate_pool: Mapping[str, Any]
    proposals_unavailable: str = ""
    current_limits: Mapping[str, Any] | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)

    def section(self, kind: str) -> ProposalSection:
        return next(section for section in self.sections if section.kind == kind)

    @property
    def proposals(self) -> tuple[Proposal, ...]:
        """Every proposal shown, in section order."""
        return tuple(p for section in self.sections for p in section.proposals)

    def as_dict(self) -> dict[str, Any]:
        return {
            "primers": list(self.primers),
            "settings": self.settings.as_dict(),
            "current": self.current.as_dict(),
            "current_worst_host_sites_per_mb": self.current_host_load.as_dict(),
            "current_panel_limits": (
                dict(self.current_limits) if self.current_limits is not None else None
            ),
            "attribution": [dict(row) for row in self.attribution],
            "sections": {section.kind: section.as_dict() for section in self.sections},
            "proposals_unavailable": self.proposals_unavailable or None,
            "candidate_pool": dict(self.candidate_pool),
            "ranking_within_a_section": [
                "gain in the worst target's coverage",
                "worst host site density (lower first; unmeasured last)",
                "number of oligos changed",
            ],
            "prediction_note": PREDICTION_NOTE,
            "notes": list(self.notes),
        }


# ----------------------------------------------------------------------
# The entry point
# ----------------------------------------------------------------------


def improve_set(
    primers: Sequence[str],
    references: Sequence[ReferenceSpec],
    conditions: Any,
    settings: ImprovementSettings,
    candidates: Sequence[str] | None = None,
    *,
    limit_check: Callable[[Sequence[str]], Mapping[str, Any]] | None = None,
    pool_unavailable: str = "",
) -> ImprovementReport:
    """Diagnose one oligo set and propose the smallest edits that help it.

    Args:
        primers: the set as the user holds it. Any mix of lengths.
        references: one `ReferenceSpec` per target and host, as `evaluate-set`
            builds them.
        conditions: `ReactionConditions`. Tm is read from it, so without it no
            add can be checked against the Tm window and none is proposed.
        settings: the limits a proposal is held to.
        candidates: the pool an add is drawn from, or None when there is none.
        limit_check: evaluates the configured panel limits on a set and returns
            a JSON-ready mapping. Advisory: its outcome is attached, and it
            removes no proposal.
        pool_unavailable: why there is no candidate pool, for the report.

    Returns:
        An `ImprovementReport`. It describes; nothing is applied.
    """
    panel = _normalise(primers)
    if not panel:
        raise ValueError("improve_set needs at least one oligo")
    specs = list(references)
    evaluator = _Evaluator(specs, conditions, int(settings.extension))
    fixed = {str(p).upper() for p in settings.fixed_oligos}

    pool, pool_record = _screen_pool(panel, candidates, conditions, settings, fixed)
    pool_record["unavailable"] = pool_unavailable or pool_record.get("unavailable")

    # One evaluation over the set and every candidate that survived the hard
    # screens, so each reference is read once. Its figures are not used.
    eligible = [c for c, _partners in pool]
    if eligible:
        evaluator.quick(panel + eligible)

    current = evaluator.full(panel)
    without = {o: evaluator.quick([p for p in panel if p != o]) for o in panel if len(panel) > 1}
    screen = _configured_screen(panel, settings, conditions)
    attribution = tuple(
        _attribute(o, panel, current, without.get(o), screen, conditions, settings, fixed)
        for o in panel
    )

    worst_now = current.worst_target_coverage
    drafts: list[_Draft] = []
    reason = ""
    if worst_now.value is None:
        reason = f"no edit is proposed because the worst target is unknown: {worst_now.unavailable}"
    else:
        drafts.extend(_draft_drops(panel, current, without, fixed))
        drafts.extend(_draft_adds(panel, pool, current, evaluator, pool_record))

    limit = max(0, int(settings.max_edits))
    sections = []
    for kind in SECTION_ORDER:
        ranked = sorted((d for d in drafts if d.kind == kind), key=_rank_key)
        shown = tuple(_finish(d, panel, current, evaluator, limit_check) for d in ranked[:limit])
        sections.append(ProposalSection(kind=kind, proposals=shown, considered=len(ranked)))

    return ImprovementReport(
        primers=tuple(panel),
        settings=settings,
        current=current,
        current_host_load=_host_load(current),
        attribution=attribution,
        sections=tuple(sections),
        candidate_pool=pool_record,
        proposals_unavailable=reason,
        current_limits=limit_check(panel) if limit_check is not None else None,
        notes=tuple(current.notes),
    )


# ----------------------------------------------------------------------
# Evaluation, shared reading
# ----------------------------------------------------------------------


class _Evaluator:
    """`evaluate_reference_panel` on many panels, reading each reference once.

    `quick` leaves the occupancy-weighted load out, because that load reads the
    mismatch-class tables per panel and no ranking below uses it. `full` is the
    evaluation `evaluate-set` reports, and is what a proposal's predicted
    figures come from.
    """

    def __init__(self, specs: Sequence[ReferenceSpec], conditions: Any, extension: int):
        self.specs = list(specs)
        self.conditions = conditions
        self.extension = int(extension)
        self.sources = MemoisedSources()
        self._quick: dict[tuple[str, ...], ReferencePanelAssessment] = {}

    def quick(self, panel: Sequence[str]) -> ReferencePanelAssessment:
        key = tuple(panel)
        if key not in self._quick:
            self._quick[key] = evaluate_reference_panel(
                list(panel), self.specs, None, extension=self.extension, sources=self.sources
            )
        return self._quick[key]

    def full(self, panel: Sequence[str]) -> ReferencePanelAssessment:
        return evaluate_reference_panel(
            list(panel),
            self.specs,
            self.conditions,
            extension=self.extension,
            sources=self.sources,
        )


def _normalise(primers: Sequence[str] | None) -> list[str]:
    out: list[str] = []
    for primer in primers or ():
        sequence = str(primer).strip().upper()
        if sequence and sequence not in out:
            out.append(sequence)
    return out


def _host_load(assessment: ReferencePanelAssessment) -> Measurement:
    """The highest host site density, each host against its own length.

    Unavailable when any host is, as every reduction in the evaluation is: the
    worst over the hosts that happened to answer is not the worst.
    """
    name, units = "worst_host_sites_per_mb", "sites/Mb"
    hosts = assessment.hosts
    if not hosts:
        return Measurement(name, None, units, unavailable="no host reference was given")
    missing = [h.prefix for h in hosts if h.site_density.value is None]
    if missing:
        return Measurement(
            name,
            None,
            units,
            unavailable=(
                f"site density is unavailable for {len(missing)} of {len(hosts)} host(s) "
                f"({', '.join(missing)}), so the highest over them is unknown"
            ),
        )
    worst = max(hosts, key=lambda h: h.site_density.value)
    return Measurement(
        name,
        float(worst.site_density.value),
        units,
        basis=f"highest of {len(hosts)} host(s); {worst.prefix}",
    )


# ----------------------------------------------------------------------
# Step 2: attribution
# ----------------------------------------------------------------------


def _configured_screen(pool: Sequence[str], settings: ImprovementSettings, conditions: Any):
    """The configured dimer screen over `pool`, through the one door.

    One construction, and it always forwards the stability floor, as every
    site that builds a screen must
    (`tests/test_the_dimer_screen_can_carry_a_stability_floor.py`): a site that
    forwarded only the run limit would drop a configured floor in silence.

    The floor (`max_dimer_dg`) is a free energy at the reaction temperature.
    With a floor configured and no temperature to evaluate it at, this raises:
    screening at an assumed temperature would report pairs as compatible under
    a reaction nobody specified. With no floor the temperature is not used, so
    an absent one is simply not passed.
    """
    from neoswga.core.lazy_dimer import dimer_screen

    floor = settings.max_dimer_dg
    temp = getattr(conditions, "temp", None)
    if floor is not None and temp is None:
        raise ValueError(
            f"max_dimer_dg is set ({floor} kcal/mol) and the reaction "
            "conditions carry no temperature, so the dimer stability floor cannot be "
            "evaluated. Supply reaction conditions, or remove max_dimer_dg."
        )
    reaction = {} if temp is None else {"temp": float(temp)}
    return dimer_screen(list(pool), int(settings.max_dimer_bp), max_dimer_dg=floor, **reaction)


def _tm(oligo: str, conditions: Any, settings: ImprovementSettings) -> dict[str, Any]:
    window = [float(settings.min_tm), float(settings.max_tm)]
    if conditions is None:
        return {
            "tm_c": None,
            "window_c": window,
            "in_window": None,
            "unavailable": "no reaction conditions were supplied, so Tm is undefined",
        }
    tm = float(conditions.calculate_effective_tm(oligo))
    return {
        "tm_c": tm,
        "window_c": window,
        "in_window": bool(window[0] <= tm <= window[1]),
        "unavailable": None,
    }


def _sites_on(record: ReferenceRecord, oligo: str) -> dict[str, Any]:
    """This oligo's sites on one reference: a count, or why there is none."""
    sites = record.per_primer_sites.get(oligo)
    status = record.per_primer_status.get(oligo)
    if sites is None:
        if status == "not_indexed":
            why = (
                "absent from the position index and no FASTA was scanned; its "
                "sites here are unknown, which is not the same as none"
            )
        else:
            why = record.unavailable or "this reference was not measured"
        return {"sites": None, "status": "unavailable", "unavailable": why}
    out: dict[str, Any] = {"sites": int(sites), "status": status or "ok", "unavailable": None}
    if record.length > 0:
        out["sites_per_mb"] = float(sites) * 1e6 / float(record.length)
    return out


def _marginal(
    oligo: str,
    record: ReferenceRecord,
    without: ReferencePanelAssessment | None,
) -> Measurement:
    """Coverage of one target with the set, minus coverage without this oligo."""
    name, units = "marginal_coverage", "fraction"
    if without is None:
        return Measurement(
            name, None, units, unavailable="the set holds one oligo, so there is no set without it"
        )
    if record.coverage.value is None:
        return Measurement(name, None, units, unavailable=record.coverage.unavailable)
    other = without.record(record.prefix)
    if other is None or other.coverage.value is None:
        why = other.coverage.unavailable if other is not None else "reference not evaluated"
        return Measurement(name, None, units, unavailable=f"without {oligo}: {why}")
    return Measurement(
        name,
        float(record.coverage.value) - float(other.coverage.value),
        units,
        basis="coverage with the set minus coverage with this oligo left out",
    )


def _attribute(
    oligo: str,
    panel: Sequence[str],
    current: ReferencePanelAssessment,
    without: ReferencePanelAssessment | None,
    screen: Any,
    conditions: Any,
    settings: ImprovementSettings,
    fixed: set[str],
) -> dict[str, Any]:
    """The diagnosis for one oligo. Worth reporting with no edit proposed."""
    targets: dict[str, Any] = {}
    for record in current.targets:
        row = _sites_on(record, oligo)
        marginal = _marginal(oligo, record, without)
        row["marginal_coverage"] = marginal.as_dict()
        # The only oligo covering some region of this target, which is what a
        # non-zero marginal coverage means. Unknown when the marginal is.
        row["sole_cover_of_some_region"] = None if marginal.value is None else marginal.value > 0
        targets[record.prefix] = row
    hosts = {record.prefix: _sites_on(record, oligo) for record in current.hosts}
    partners = [other for other in panel if other != oligo and screen.dimerises(oligo, [other])]
    return {
        "oligo": oligo,
        "length": len(oligo),
        "fixed": oligo in fixed,
        "tm": _tm(oligo, conditions, settings),
        "dimer_partners": partners,
        "dimer_limit_bp": int(settings.max_dimer_bp),
        "targets": targets,
        "hosts": hosts,
    }


# ----------------------------------------------------------------------
# Step 3: drafts, and the hard screens of step 4
# ----------------------------------------------------------------------


@dataclass(frozen=True)
class _Draft:
    """An edit and the evaluation it was ranked on, before the full one."""

    kind: str
    drop: tuple[str, ...]
    add: tuple[str, ...]
    reason: str
    quick: ReferencePanelAssessment
    gain: float
    host_load: Measurement
    depends_on: str = ""


def _rank_key(draft: _Draft):
    """Worst-target gain, then host load, then how much of the set changes.

    An unmeasured host load sorts after every measured one at the same gain
    rather than as a load of zero. The trailing names make the order stable.
    """
    load = draft.host_load.value
    return (
        -draft.gain,
        load is None,
        load if load is not None else 0.0,
        len(draft.drop) + len(draft.add),
        draft.drop,
        draft.add,
    )


def _gain(current: ReferencePanelAssessment, after: ReferencePanelAssessment) -> float | None:
    before, value = current.worst_target_coverage.value, after.worst_target_coverage.value
    if before is None or value is None:
        return None
    return float(value) - float(before)


def _density_text(assessment: ReferencePanelAssessment) -> str:
    """The worst density as text. With no host site left it is unbounded.

    The evaluation stands a fixed large number in for that case, and printing
    the number would read as a measured ratio.
    """
    if assessment.worst_host_density.value is None:
        return "unavailable"
    if assessment.pairs and all(pair.zero_host_sites for pair in assessment.pairs):
        return "unbounded (no host site)"
    return f"{assessment.worst_host_density.value:.4g}"


def _draft_drops(
    panel: Sequence[str],
    current: ReferencePanelAssessment,
    without: Mapping[str, ReferencePanelAssessment],
    fixed: set[str],
) -> list[_Draft]:
    """Drops that cost no coverage, and drops that trade coverage for density.

    The two are different kinds and are reported apart. The first is an oligo
    the set does as well without, on every target. The second is a comparison
    of two measured densities and nothing more: it is satisfied by the member
    of any set with the lowest target-to-host ratio, so it is listed as a
    trade-off and never as a recommendation.
    """
    drafts: list[_Draft] = []
    for oligo in panel:
        after = without.get(oligo)
        if oligo in fixed or after is None:
            continue
        gain = _gain(current, after)
        if gain is None:
            continue
        marginals = [_marginal(oligo, record, after).value for record in current.targets]
        target_sites = [record.per_primer_sites.get(oligo) for record in current.targets]
        if all(value is not None and value == 0 for value in marginals):
            if all(sites == 0 for sites in target_sites):
                reason = "binds no target: no site on any target reference"
            else:
                reason = (
                    "no marginal contribution: every base it covers on every "
                    "target is covered by another oligo"
                )
            drafts.append(_Draft(KIND_DROP, (oligo,), (), reason, after, gain, _host_load(after)))
            continue
        before_density = current.worst_host_density.value
        after_density = after.worst_host_density.value
        if before_density is None or after_density is None or after_density <= before_density:
            continue
        reason = (
            "removal raises the worst target-against-host selectivity density from "
            f"{_density_text(current)} to {_density_text(after)}; the coverage it "
            "costs is listed per target"
        )
        drafts.append(_Draft(KIND_TRADE_OFF, (oligo,), (), reason, after, gain, _host_load(after)))
    return drafts


def _screen_pool(
    panel: Sequence[str],
    candidates: Sequence[str] | None,
    conditions: Any,
    settings: ImprovementSettings,
    fixed: set[str],
) -> tuple[list[tuple[str, tuple[str, ...]]], dict[str, Any]]:
    """The candidates that clear the hard constraints, with what each displaces.

    Returns `(candidate, partners)` pairs. An empty `partners` is a candidate
    that can be added beside the whole set. One partner that is not fixed is a
    candidate that can only replace that oligo. Anything else conflicts with an
    oligo that stays and is not an option.
    """
    record: dict[str, Any] = {
        "size": 0,
        "eligible_add": 0,
        "eligible_swap": 0,
        "rejected_dimer": 0,
        "rejected_tm": 0,
        "rejected_excluded": 0,
        "unmeasured": 0,
        "unmeasured_examples": [],
        "unavailable": None,
    }
    pool = [c for c in _normalise(candidates) if c not in panel]
    record["size"] = len(pool)
    if not pool:
        if candidates is not None:
            record["unavailable"] = "the candidate pool holds no oligo outside the set"
        else:
            record["unavailable"] = "no candidate pool was supplied"
        return [], record
    if conditions is None:
        record["unavailable"] = (
            "no reaction conditions were supplied, so no candidate can be "
            "checked against the Tm window and none is proposed"
        )
        return [], record

    excluded = {str(p).upper() for p in settings.excluded_oligos}
    in_window: list[str] = []
    for candidate in pool:
        if candidate in excluded or set(candidate) - set("ACGT"):
            record["rejected_excluded"] += 1
        elif not _tm(candidate, conditions, settings)["in_window"]:
            record["rejected_tm"] += 1
        else:
            in_window.append(candidate)

    # Built over the set and every candidate it will be asked about: a dense
    # matrix answers "no dimer" for a sequence it was not built over.
    screen = _configured_screen([*panel, *in_window], settings, conditions)
    kept: list[tuple[str, tuple[str, ...]]] = []
    for candidate in in_window:
        partners = tuple(p for p in panel if screen.dimerises(candidate, [p]))
        if not partners:
            record["eligible_add"] += 1
            kept.append((candidate, ()))
        elif len(partners) == 1 and partners[0] not in fixed and len(panel) > 1:
            record["eligible_swap"] += 1
            kept.append((candidate, partners))
        else:
            record["rejected_dimer"] += 1
    return kept, record


def _unknown_sites(after: ReferencePanelAssessment, candidate: str) -> list[tuple[str, str]]:
    """References, target or host, on which this candidate's sites are unknown.

    A candidate absent from a host index that may not be scanned has an unknown
    host load. That is not a low one, and a candidate carrying it is not an
    option: proposing it would rank an unmeasured host cost beside measured
    ones.
    """
    unknown = []
    for record in after.references:
        if record.per_primer_sites.get(candidate) is None:
            unknown.append((record.prefix, _sites_on(record, candidate)["unavailable"]))
    return unknown


def _draft_adds(
    panel: Sequence[str],
    pool: Sequence[tuple[str, tuple[str, ...]]],
    current: ReferencePanelAssessment,
    evaluator: _Evaluator,
    pool_record: dict[str, Any],
) -> list[_Draft]:
    """Adds and swaps that raise the worst target's coverage.

    Judged on the worst target, not on the pooled figure: a candidate that adds
    coverage only where the set already has it does not appear, however much it
    would raise the mean.
    """
    worst = current.worst_target_coverage
    drafts: list[_Draft] = []
    if not pool:
        return drafts
    # Only reached with every target measured: the caller proposes nothing
    # while the worst target is unknown.
    worst_target = min(current.targets, key=lambda record: record.coverage.value)
    for candidate, partners in pool:
        if partners:
            kind, drop = KIND_SWAP, partners
            edited = [p for p in panel if p not in partners] + [candidate]
            depends = (
                f"{candidate} dimerises with {partners[0]} above the configured "
                f"limit, so it can only replace it"
            )
        else:
            kind, drop, edited, depends = KIND_ADD, (), [*panel, candidate], ""
        after = evaluator.quick(edited)
        gain = _gain(current, after)
        unknown = _unknown_sites(after, candidate)
        if unknown or gain is None:
            # What the candidate would do on some reference is unknown. Counted
            # with the reason, not proposed.
            pool_record["unmeasured"] += 1
            if len(pool_record["unmeasured_examples"]) < 5:
                pool_record["unmeasured_examples"].append(
                    {
                        "candidate": candidate,
                        "references": [prefix for prefix, _why in unknown],
                        "reason": (
                            unknown[0][1] if unknown else after.worst_target_coverage.unavailable
                        ),
                    }
                )
            continue
        if gain <= 0:
            continue
        reason = (
            f"raises the worst target's coverage ({_short(worst_target.prefix)}) from "
            f"{worst.value:.4f}; the lowest over the targets is then "
            f"{after.worst_target_coverage.value:.4f}"
        )
        drafts.append(
            _Draft(kind, tuple(drop), (candidate,), reason, after, gain, _host_load(after), depends)
        )
    return drafts


def _finish(
    draft: _Draft,
    panel: Sequence[str],
    current: ReferencePanelAssessment,
    evaluator: _Evaluator,
    limit_check: Callable[[Sequence[str]], Mapping[str, Any]] | None,
) -> Proposal:
    """The reported proposal: the full evaluation, and the advisory limits."""
    resulting = [p for p in panel if p not in draft.drop] + list(draft.add)
    predicted = evaluator.full(resulting)
    change: dict[str, dict[str, float | None]] = {}
    for record in current.targets:
        after = predicted.record(record.prefix)
        before_value = record.coverage.value
        after_value = after.coverage.value if after is not None else None
        change[record.prefix] = {
            "before": before_value,
            "after": after_value,
            "change": (
                None
                if before_value is None or after_value is None
                else float(after_value) - float(before_value)
            ),
        }
    return Proposal(
        kind=draft.kind,
        drop=draft.drop,
        add=draft.add,
        reason=draft.reason,
        depends_on=draft.depends_on,
        resulting_set=tuple(resulting),
        predicted=predicted,
        worst_target_gain=float(draft.gain),
        host_load=_host_load(predicted),
        coverage_change=change,
        density_before=current.worst_host_density,
        density_after=predicted.worst_host_density,
        advisory=limit_check(resulting) if limit_check is not None else None,
    )


# ----------------------------------------------------------------------
# The printed table
# ----------------------------------------------------------------------


def _short(prefix: str) -> str:
    return str(prefix).replace("\\", "/").rsplit("/", 1)[-1]


def _show(entry: Mapping[str, Any] | None, fmt: str = "{:.4g}") -> str:
    value = (entry or {}).get("value")
    if value is None:
        return "unavailable"
    return fmt.format(value)


def _reference_lines(current: Mapping[str, Any]) -> list[str]:
    lines = ["", "Current set, per reference (each against its own length):"]
    rows = list((current.get("per_target") or {}).items())
    rows += list((current.get("per_host") or {}).items())
    for name, record in rows:
        lines.append(
            f"  {record.get('role', '?'):<6s} {_short(name):<24s} "
            f"sites {_show(record.get('sites'), '{:.0f}'):>9s}  "
            f"{_show(record.get('sites_per_mb'), '{:.1f}'):>9s}/Mb  "
            f"coverage {_show(record.get('coverage'), '{:.2%}')}"
        )
        if record.get("status") == "unavailable":
            lines.append(f"         why: {record.get('unavailable')}")
    lines.append(
        f"  coverage reach: {current.get('extension_reach_bp')} bp; "
        f"worst target coverage: {_show(current.get('worst_target_coverage'), '{:.2%}')}"
    )
    worst = current.get("worst_target_coverage") or {}
    if worst.get("value") is None and worst.get("unavailable"):
        lines.append(f"         why: {worst['unavailable']}")
    lines.append(
        "  worst target-against-host density: "
        f"{_show(current.get('worst_host_selectivity_density'))}"
    )
    return lines


def _attribution_lines(rows: Sequence[Mapping[str, Any]]) -> list[str]:
    lines = [
        "",
        "Per oligo (marginal = coverage lost on that target without it):",
        f"  {'oligo':<22s} {'Tm':>6s} {'win':>4s} {'fix':>4s} {'dimers':>6s}  per reference",
    ]
    for row in rows:
        tm = row["tm"]
        tm_text = "n/a" if tm["tm_c"] is None else f"{tm['tm_c']:.1f}"
        window = "n/a" if tm["in_window"] is None else ("yes" if tm["in_window"] else "NO")
        parts = []
        for name, target in row["targets"].items():
            if target["sites"] is None:
                parts.append(f"{_short(name)}: unavailable")
                continue
            marginal = target["marginal_coverage"]["value"]
            shown = "n/a" if marginal is None else f"{marginal:+.2%}"
            parts.append(f"{_short(name)}: {target['sites']} sites, marginal {shown}")
        for name, host in row["hosts"].items():
            shown = "unavailable" if host["sites"] is None else f"{host['sites']} sites"
            parts.append(f"host {_short(name)}: {shown}")
        lines.append(
            f"  {row['oligo']:<22s} {tm_text:>6s} {window:>4s} "
            f"{'yes' if row['fixed'] else '-':>4s} {len(row['dimer_partners']):>6d}  "
            + "; ".join(parts)
        )
    return lines


def _one_proposal_lines(proposal: Mapping[str, Any]) -> list[str]:
    edit = []
    if proposal["drop"]:
        edit.append("drop " + ", ".join(proposal["drop"]))
    if proposal["add"]:
        edit.append("add " + ", ".join(proposal["add"]))
    lines = [f"    {proposal['rank']}. {'; '.join(edit)}", f"         {proposal['reason']}"]
    if proposal.get("depends_on"):
        lines.append(f"         depends on: {proposal['depends_on']}")
    for name, change in proposal["per_target_coverage"].items():
        if change["change"] is None:
            lines.append(f"         {_short(name)}: coverage unavailable")
        else:
            lines.append(
                f"         {_short(name)}: coverage {change['before']:.2%} -> "
                f"{change['after']:.2%} ({change['change']:+.2%})"
            )
    lines.append(
        "         worst host site density: "
        f"{_show(proposal.get('worst_host_sites_per_mb'), '{:.1f}')} sites/Mb"
    )
    limits = proposal.get("panel_limits")
    if limits and not limits.get("evaluated", True):
        lines.append(f"         panel limits not evaluated: {limits.get('unavailable')}")
    elif limits and limits.get("violations"):
        lines.append("         misses a configured limit: " + ", ".join(limits["violations"]))
    return lines


def _section_lines(section: Mapping[str, Any]) -> list[str]:
    """One section: its title, how many are shown of how many, and its entries.

    "none found" is said only when nothing was considered. A section cut to
    zero by `--max-edits 0` still says how many there were.
    """
    considered = int(section.get("considered") or 0)
    lines = ["", f"  {section['title']}"]
    if considered == 0:
        lines.append("    none found")
        return lines
    lines.append(f"    shown {section['shown']} of {considered} considered")
    lines.append(f"    {section['note']}")
    for proposal in section.get("proposals") or []:
        lines += _one_proposal_lines(proposal)
    return lines


def _proposal_lines(report: Mapping[str, Any]) -> list[str]:
    lines = ["", "Proposed edits, by kind (single edits; none is applied):"]
    if report.get("proposals_unavailable"):
        lines.append(f"  {report['proposals_unavailable']}")
    else:
        for section in (report.get("sections") or {}).values():
            lines += _section_lines(section)
    pool = report.get("candidate_pool") or {}
    lines.append("")
    if pool.get("unavailable"):
        lines.append(f"  adds and swaps: {pool['unavailable']}")
    else:
        lines.append(
            f"  candidate pool: {pool.get('size', 0)} examined; "
            f"{pool.get('rejected_dimer', 0)} conflict with an oligo that stays, "
            f"{pool.get('rejected_tm', 0)} outside the Tm window, "
            f"{pool.get('unmeasured', 0)} with sites that could not be established"
        )
        for example in pool.get("unmeasured_examples") or []:
            lines.append(f"    not proposed, {example['candidate']}: {example['reason']}")
    lines.append(f"  {report.get('prediction_note') or PREDICTION_NOTE}")
    return lines


def report_lines(report: Mapping[str, Any]) -> list[str]:
    """The printed table, from the JSON form of an `ImprovementReport`."""
    lines = ["", "=" * 72, "SET IMPROVEMENT (report only; nothing is applied)", "=" * 72]
    lines.append(f"  Oligos: {len(report.get('primers') or [])}")
    lines += _reference_lines(report.get("current") or {})
    lines += _attribution_lines(report.get("attribution") or [])
    lines += _proposal_lines(report)
    for note in report.get("notes") or []:
        lines.append(f"  note: {note}")
    lines.append("=" * 72)
    return lines
