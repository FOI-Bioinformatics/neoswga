"""One complete acceptance record for a panel, produced once and read everywhere.

Task 5 of the 2026-09-21 valid-design plan. A panel is currently measured in
several places and the answers are assembled differently at each: the optimizer
holds `PrimerSetMetrics`, the acceptance path builds an `AcceptanceReport`, the
summary JSON is written by a third piece of code and the report renders a
fourth. Nothing makes them agree, and this project has already shipped a run
printing two coverage figures that differ by construction.

`PanelAssessment` is the one record. It carries:

- the panel and the request that produced it, by hash, so a saved result names
  its own configuration;
- named metrics WITH their units and the reach and denominator they were
  computed at, because a coverage figure without its reach carries almost no
  information: one saved 26-oligo panel reads 41.3% at 1 kb and 93.5% at 5 kb;
- which quantities were unavailable and why, separately from quantities that
  measured zero;
- every hard-constraint violation, and qualification as a single boolean that
  is true exactly when there are none.

Three rules it enforces.

**A required quantity that is NaN or infinite fails.** Not "is reported as
NaN": a non-finite coverage compares False against every threshold, so a panel
carrying one passes no limit and fails no limit, and the run reports a number
that arithmetic cannot use.

**A verified zero background is represented, not divided by.** A panel with no
host sites has infinite selectivity, which is not a number JSON has and not a
claim about enrichment. It is recorded as a zero denominator with the site
count beside it.

**Measurement is separated from policy inside, and returned together.** What
the panel is and whether it is acceptable are different questions with
different evidence, but a caller that can get one without the other will
eventually report a metric alongside a verdict computed from something else.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, Mapping, Optional, Tuple

from .exceptions import ModelEvaluationError

__all__ = [
    "Measurement",
    "PanelAssessment",
    "evaluate_panel",
]

#: Quantities a design cannot proceed without. A non-finite value for any of
#: them fails the run rather than being reported.
REQUIRED_METRICS = ("fg_coverage",)


@dataclass(frozen=True)
class Measurement:
    """One number, its unit, and what it was computed against.

    `value` is None when the quantity could not be measured, and `unavailable`
    then says why. That is not the same as a measured zero, and the two must
    not share a representation: a panel that binds the host nowhere and a panel
    whose host index was never opened both read zero otherwise.
    """

    name: str
    value: Optional[float]
    units: str
    basis: str = ""
    unavailable: str = ""

    def __post_init__(self):
        if self.value is None and not self.unavailable:
            raise ValueError(f"{self.name} has no value and no reason for not having one")
        if self.value is not None and self.unavailable:
            raise ValueError(f"{self.name} has both a value and a reason for not having one")

    def as_dict(self) -> Dict[str, Any]:
        return {
            "value": self.value,
            "units": self.units,
            "basis": self.basis,
            "unavailable": self.unavailable or None,
        }


@dataclass(frozen=True)
class PanelAssessment:
    """What a panel is, and whether it is acceptable. One record, one answer."""

    primers: Tuple[str, ...]
    request_hash: str
    metrics: Mapping[str, Measurement]
    per_target: Mapping[str, Measurement]
    violations: Tuple[str, ...]
    qualified: bool
    model_versions: Tuple[Tuple[str, str], ...] = ()
    notes: Tuple[str, ...] = ()
    #: Present only when the background genome carries no sites at all. The
    #: selectivity ratio is then undefined rather than infinite.
    zero_background: bool = False
    evidence: Mapping[str, str] = field(default_factory=dict)

    def value(self, name: str) -> Optional[float]:
        measurement = self.metrics.get(name)
        return measurement.value if measurement else None

    def as_dict(self) -> Dict[str, Any]:
        """A JSON-serializable form. No infinities, no NaN.

        The report and the saved result both read this, so that a rendered
        figure and a stored one cannot come from different arithmetic.
        """
        return {
            "primers": list(self.primers),
            "request_hash": self.request_hash,
            "qualified": self.qualified,
            "violations": list(self.violations),
            "zero_background": self.zero_background,
            "metrics": {name: m.as_dict() for name, m in sorted(self.metrics.items())},
            "per_target": {name: m.as_dict() for name, m in sorted(self.per_target.items())},
            "model_versions": [list(pair) for pair in self.model_versions],
            "evidence": dict(sorted(self.evidence.items())),
            "notes": list(self.notes),
        }


def _finite_or_fail(name: str, value, subject) -> Optional[float]:
    """Refuse a non-finite required quantity rather than reporting it.

    NaN compares False against every threshold, so a panel carrying one passes
    no limit and fails no limit. Infinity is worse: it is not representable in
    JSON and, on a selectivity ratio, reads as guaranteed enrichment.
    """
    if value is None:
        return None
    number = float(value)
    if not math.isfinite(number):
        raise ModelEvaluationError(name, subject, f"non-finite value ({value})")
    return number


def _measure(name, value, units, basis, subject, *, required=False) -> Measurement:
    if value is None:
        if required:
            raise ModelEvaluationError(name, subject, "required quantity was not computed")
        return Measurement(name, None, units, basis, unavailable="not computed for this panel")
    return Measurement(name, _finite_or_fail(name, value, subject), units, basis)


def _panel_violations(request, primers) -> list:
    """Hard properties of a delivered panel, independent of configured limits.

    Both were enforced already and neither reached this record. Requested size
    is asked three different ways elsewhere and they are NOT interchangeable:
    `result_validation.validate_result` compares `!=`, so it reports a panel
    larger than requested too; `optimization_service`'s `assess` closure
    compares `>=` on the search's own stopping rule, where changing it moves
    delivered panels. This one asks only whether the delivered panel is short.
    The dimer half is the shared rule rather than a third reading of it.

    Candidate QC is deliberately absent -- GC window, homopolymer run, Gini.
    In this codebase those are ADMISSION rules applied during `filter`, not
    properties of a delivered panel, and no panel-level path re-checks them.

    Composition is not QC and is checked: a duplicated oligo and an oligo the
    request excluded are both facts about the panel that came back, whatever
    admitted it. An earlier version of this docstring lumped the two together
    and the distinction is the reason only one of them belongs here.
    """
    found = []

    requested = getattr(request, "target_size", None)
    if requested and len(primers) < requested:
        found.append(f"panel size {len(primers)} is below the requested {requested}")

    # Composition. Both are properties of the DELIVERED panel under the request
    # rather than admission rules applied during `filter`, which is what
    # separates them from the QC gates deliberately absent below.
    # `validate_result` checked both and this record checked neither, so two
    # records described the same panel differently.
    counts: Dict[str, int] = {}
    for oligo in primers:
        counts[oligo] = counts.get(oligo, 0) + 1
    repeated = sorted(oligo for oligo, count in counts.items() if count > 1)
    if repeated:
        found.append(f"duplicate oligo in the delivered panel: {', '.join(repeated)}")

    # The request names what must not be selected. A candidate pool that was
    # never filtered against it can reach the optimizer through expand-primers
    # or swap-primer, which is the case `validate_result`'s `forbidden_primers`
    # argument exists for.
    excluded = set(getattr(request, "excluded_oligos", ()) or ())
    reinjected = sorted(set(primers) & excluded)
    if reinjected:
        found.append(f"excluded oligo in the delivered panel: {', '.join(reinjected)}")

    # The delivered-heterodimer rule is `dimer.dimer_validation_issue`, asked
    # for rather than rewritten. This measured complementary runs itself and
    # agreed with it on every panel a resolved request can produce, which made
    # the second implementation invisible until the two drifted --
    # `string_search`'s two scanners are what that costs.
    # `optimization_service.panel_violations` is the SEARCH-time screen and a
    # wider rule: self-dimers against `max_self_dimer_bp` and the optional
    # `max_dimer_dg` floor, reported as an unattributed "dimer constraint".
    # Adopting it here would fault a delivered pool on an admission threshold.
    from neoswga.core.dimer import dimer_validation_issue

    issue = dimer_validation_issue(list(primers), getattr(request, "max_dimer_bp", None))
    if issue is not None:
        found.append(issue["detail"])
    return found


def evaluate_panel(request, oligos, metrics, *, objective=None) -> PanelAssessment:
    """One assessment of one panel under one request.

    `metrics` is the evaluator's own `PrimerSetMetrics`, passed in rather than
    recomputed here: this function's job is to give one shape to what was
    measured and to decide acceptance from it, not to add a fourth place where
    coverage is calculated.

    `objective` supplies the configured hard limits. With none, the panel has
    no limits to fail and qualifies on having been measured at all, which is
    the documented behaviour when no panel limit is configured: setting none
    must leave the delivered panel byte-identical.
    """
    primers = tuple(oligos)
    subject = f"panel of {len(primers)}"
    reach = getattr(request, "coverage_reach", None)
    basis = f"reach {reach} bp" if reach else "reach not recorded"

    named: Dict[str, Measurement] = {}
    named["fg_coverage"] = _measure(
        "fg_coverage",
        getattr(metrics, "fg_coverage", None),
        "fraction of target bases",
        f"{basis}; denominator is total target length",
        subject,
        required=True,
    )
    # The record reported the geometric figure while the live acceptance path
    # selects on this one, so the two answered different questions about the
    # same panel. Both are kept and both say which they are.
    named["effective_fg_coverage"] = _measure(
        "effective_fg_coverage",
        getattr(metrics, "effective_fg_coverage", None),
        "fraction of target bases",
        f"{basis}; occupancy-weighted, the quantity selection uses",
        subject,
    )
    named["bg_coverage"] = _measure(
        "bg_coverage",
        getattr(metrics, "bg_coverage", None),
        "fraction of host bases",
        f"{basis}; denominator is total host length",
        subject,
    )
    for name, units, note in (
        ("total_fg_sites", "exact-match sites", "exact matches only"),
        ("total_bg_sites", "exact-match sites", "exact matches only"),
        ("mean_gap", "bp", "between adjacent target sites"),
        ("max_gap", "bp", "largest target hole"),
        ("gap_gini", "dimensionless", "evenness of target gaps"),
        ("mean_tm", "C", "under the resolved reaction"),
    ):
        named[name] = _measure(name, getattr(metrics, name, None), units, note, subject)

    # A verified zero background is a measurement and a useful one. It is not a
    # selectivity of infinity, which is neither a number JSON carries nor a
    # claim anybody has evidence for.
    zero_background = named["total_bg_sites"].value == 0
    if zero_background:
        named["selectivity_ratio"] = Measurement(
            "selectivity_ratio",
            None,
            "dimensionless",
            "undefined: the host carries no exact-match site for this panel",
            unavailable="zero denominator",
        )
    else:
        named["selectivity_ratio"] = _measure(
            "selectivity_ratio",
            getattr(metrics, "selectivity_ratio", None),
            "dimensionless",
            "effective foreground load over effective background load",
            subject,
        )

    per_target: Dict[str, Measurement] = {}
    for prefix, value in (getattr(metrics, "per_target_coverage", None) or {}).items():
        per_target[str(prefix)] = _measure(
            f"coverage[{prefix}]", value, "fraction of target bases", basis, subject
        )

    # Hard properties of the delivered panel, each already enforced elsewhere
    # and brought here so one record can be believed. `objective` keeps owning
    # the CONFIGURED limits; these three are not configurable and never were.
    violations = list(objective.violations(primers)) if objective is not None else []
    violations.extend(_panel_violations(request, primers))

    notes = []
    if zero_background:
        notes.append(
            "The host carries no exact-match site for this panel. That is a "
            "measurement, not an enrichment estimate: mismatched binding is "
            "not counted here."
        )
    if named["fg_coverage"].value is not None:
        notes.append(
            "Coverage is a geometric proxy at the declared reach. It is not a "
            "predicted sequencing breadth and has not been calibrated against one."
        )

    return PanelAssessment(
        primers=primers,
        request_hash=getattr(request, "request_hash", ""),
        metrics=named,
        per_target=per_target,
        violations=tuple(violations),
        qualified=not violations,
        model_versions=tuple(getattr(request, "model_versions", ()) or ()),
        notes=tuple(notes),
        zero_background=zero_background,
        evidence={
            "coverage_reach": "assumed (design-density convention; never measured)",
            "mismatch_penalty": "assumed (uniform in identity and position)",
        },
    )
