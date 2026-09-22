"""What each chemistry constant is, and what the model built on it supports.

A number in this codebase can be one of five things, and until now they all
reached a design with equal authority:

- **measured**: the cited work reports this value for a case the model applies
  it to;
- **estimated**: extrapolated from cited data taken at another temperature, on
  longer DNA, or in another buffer;
- **empirical**: chosen so the model behaves plausibly, with no source;
- **assumed**: a modelling decision with a stated reason and no measurement;
- **absent**: no model computes this effect, and the code must not report zero
  for it.

The distinction that matters is between the first two. A coefficient with a
primary citation beside it reads as measured whatever it is, so the registry
assigns `status` from whether the cited work covers THIS case rather than from
whether a citation exists. Most of the additive coefficients are 37 C figures
for PCR-length duplexes, extrapolated to a 30 C isothermal reaction on 12-mers.
That extrapolation may well be fine. It is not a measurement of it.

**On provenance.** Every `source` in the registry is the attribution this
repository already carried. The primary literature was not re-read while
compiling it, so a source records who a value is attributed to, not that the
attribution was checked. Saying so is the point: the plan's rule is that an
assumption is not promoted to a measurement because it has a citation, and a
registry that quietly implied verification would break that rule in the act of
recording it.

**Why machine-readable.** `docs/SCIENCE_CITATIONS.md` already held most of this
in prose, and had drifted: it states Klenow processivity as 10,000 bp citing
Bambara (1978) while the shipped registry says 40 bp. A ledger nothing checks
drifts from the code it describes. `tests/test_model_evidence_contract.py`
pins the registry against the constants the code actually uses.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Tuple

from .exceptions import ReferenceDataError, UnsupportedModelError

__all__ = [
    "EVIDENCE_STATUSES",
    "EvidenceRecord",
    "load_evidence",
    "require_model_support",
]

#: The five things a constant can be. Ordered from strongest to weakest claim.
EVIDENCE_STATUSES = ("measured", "estimated", "empirical", "assumed", "absent")

_EVIDENCE_FILENAME = "model_evidence.json"
_CACHE: Optional[Dict[str, "EvidenceRecord"]] = None


@dataclass(frozen=True)
class EvidenceRecord:
    """One constant, what it is, and the domain the model supports it over."""

    quantity: str
    units: str
    model: str
    value: str
    status: str
    source: str
    sequence_domain: str
    buffer_domain: str
    temperature_domain: str
    uncertainty: str
    notes: str

    #: The temperature span the coefficient is recorded for, in Celsius, or
    #: None. Optional because exactly one of the 23 records states a range:
    #: ten give a single reference point and twelve give no number. A record
    #: without one is NOT refused -- absence of a domain is not a domain of
    #: zero, and refusing on it would turn a gap in the evidence into a gap in
    #: the tool.
    temperature_range_c: Optional[Tuple[float, float]] = None
    temperature_range_note: str = ""

    @property
    def is_supported(self) -> bool:
        """Whether a computation may rest on this at all.

        `absent` is the only status that is not a model. The other four are
        models of varying evidential strength, and refusing to compute with an
        estimate would refuse most of the chemistry in this field.
        """
        return self.status != "absent"


def _span(value) -> Optional[Tuple[float, float]]:
    """A recorded [low, high] temperature span, or None.

    Anything malformed is None rather than an error: a registry that refused
    to load over a bad optional field would take the whole package down for a
    quantity nothing may be consulting.
    """
    try:
        low, high = value  # type: ignore[misc]
        return (float(low), float(high))
    except (TypeError, ValueError):
        return None


def load_evidence() -> Dict[str, EvidenceRecord]:
    """The registry, by quantity. Raises if it is missing or malformed.

    Absent evidence is not evidence of no constraint, so a missing artifact
    fails rather than yielding an empty registry that every check then passes.
    Read through `importlib.resources` so an installed package finds it without
    a repository checkout.
    """
    global _CACHE
    if _CACHE is not None:
        return _CACHE

    import importlib.resources

    try:
        handle = importlib.resources.files("neoswga.core.registry") / _EVIDENCE_FILENAME
        payload = json.loads(handle.read_text())
    except (OSError, FileNotFoundError, ModuleNotFoundError, ValueError) as exc:
        raise ReferenceDataError(
            f"model evidence registry ({_EVIDENCE_FILENAME})",
            f"could not be loaded: {exc}",
            "Reinstall the package; this artifact ships as package data.",
        ) from exc

    if not isinstance(payload, dict) or "records" not in payload:
        raise ReferenceDataError(
            f"model evidence registry ({_EVIDENCE_FILENAME})",
            "has no 'records' list; its schema is not the one this version reads",
            "Reinstall the package.",
        )

    records: Dict[str, EvidenceRecord] = {}
    for entry in payload["records"]:
        try:
            record = EvidenceRecord(
                quantity=str(entry["quantity"]),
                units=str(entry.get("units", "")),
                model=str(entry.get("model", "")),
                value=str(entry.get("value", "")),
                status=str(entry["status"]),
                source=str(entry.get("source", "")),
                sequence_domain=str(entry.get("sequence_domain", "")),
                buffer_domain=str(entry.get("buffer_domain", "")),
                temperature_domain=str(entry.get("temperature_domain", "")),
                uncertainty=str(entry.get("uncertainty", "")),
                notes=str(entry.get("notes", "")),
                temperature_range_c=_span(entry.get("temperature_range_c")),
                temperature_range_note=str(entry.get("temperature_range_note", "")),
            )
        except KeyError as exc:
            raise ReferenceDataError(
                f"model evidence registry ({_EVIDENCE_FILENAME})",
                f"a record is missing the required field {exc}",
                "Reinstall the package.",
            ) from exc
        if record.status not in EVIDENCE_STATUSES:
            raise ReferenceDataError(
                f"model evidence registry ({_EVIDENCE_FILENAME})",
                f"{record.quantity} declares status {record.status!r}, which is not "
                f"one of {', '.join(EVIDENCE_STATUSES)}",
                "Reinstall the package.",
            )
        records[record.quantity] = record

    _CACHE = records
    return records


# Additive request fields, paired with the evidence record covering their effect
# on duplex stability. The pairing is explicit rather than derived from the
# field name so that adding an additive without an evidence record fails here
# instead of passing unexamined.
_ADDITIVE_EVIDENCE = {
    "dmso_percent": "tm_dmso",
    "betaine_m": "tm_betaine_uniform",
    "trehalose_m": "tm_trehalose",
    "formamide_percent": "tm_formamide",
    "ethanol_percent": "tm_ethanol",
    "urea_m": "tm_urea",
    "tmac_m": "tm_tmac_uniform",
    "propanediol_m": "tm_propanediol",
    "glycerol_percent": "tm_glycerol",
    "bsa_ug_ml": "tm_bsa",
    "peg_percent": "tm_peg",
}


def require_model_support(request) -> None:
    """Refuse a requested computation that falls outside a recorded domain.

    Called once, before any search. Everything it checks is knowable from the
    resolved request, so refusing here costs nothing; the alternative is a
    number produced with no evidence behind it, which looks exactly like one
    produced with evidence.

    It does NOT refuse estimates. Most of the additive chemistry in this field
    is extrapolated from PCR-length duplexes at 37 C, and a model that declined
    to run on an estimate would decline to run. What it refuses is a
    computation nothing performs: an enzyme with no parameter set, an oligo
    length outside the one the enzyme's parameters cover, and an additive whose
    duplex effect no term computes while the literature expects one.
    """
    evidence = load_evidence()

    _require_polymerase_support(request)
    _require_length_support(request)
    _require_additive_support(request, evidence)
    _require_additive_temperature_support(request, evidence)


def _require_polymerase_support(request) -> None:
    from .registry import views

    characteristics = views.as_characteristics()
    if request.polymerase not in characteristics:
        raise UnsupportedModelError(
            "polymerase parameters",
            request.polymerase,
            ", ".join(sorted(characteristics)),
        )


def _require_length_support(request) -> None:
    """Oligo lengths outside the range the enzyme's parameters were set for.

    A Bst design filtered through phi29's 6-12 bp window is a recorded defect
    in this project. The Tm computed for a 6-mer under Bst at 63 C is not
    wrong in a way anyone can see: it is a number with no evidence behind it,
    because the parameter set covers 15-25 nt for that enzyme.

    An empty length range passes: a request that names none has not asked for
    anything outside one.
    """
    from .registry import polymerases

    entry = polymerases.POLYMERASES.get(request.polymerase)
    supported = getattr(entry, "primer_length_range", None)
    if not supported or not request.primer_lengths:
        return

    low, high = int(supported[0]), int(supported[1])
    outside = sorted(k for k in request.primer_lengths if k < low or k > high)
    if outside:
        raise UnsupportedModelError(
            f"oligo length for {request.polymerase}",
            outside,
            f"{low}-{high} nt",
        )


def _require_additive_support(request, evidence: Mapping[str, EvidenceRecord]) -> None:
    """An additive set above zero whose duplex effect nothing computes.

    Three outcomes, and the registry decides which:

    - a record with a model: the effect is computed, whatever its evidential
      strength, and the request passes;
    - `assumed` with no model, meaning the absence is deliberate because the
      additive does not act on duplex stability (BSA, PEG): passes;
    - `absent`, meaning the literature expects an effect and nothing computes
      it: refused, because reporting the same Tm with and without it is the
      "known zero effect" the contract forbids.
    """
    conditions = request.conditions
    for field, quantity in _ADDITIVE_EVIDENCE.items():
        value = getattr(conditions, field, 0.0) or 0.0
        if not value:
            continue
        record = evidence.get(quantity)
        if record is None:
            raise UnsupportedModelError(
                "additive", field, "no evidence record; add one before using it"
            )
        if record.status == "absent":
            raise UnsupportedModelError(
                f"melting-temperature effect of {field}",
                value,
                "no additive with a recorded duplex model",
            )


def _require_additive_temperature_support(request, evidence: Mapping[str, EvidenceRecord]) -> None:
    """An additive used outside the temperature its coefficient covers.

    A Bst design at 63 C with DMSO returned an effective Tm of 58.72 C,
    computed from a coefficient this registry records as a 37 C reference
    extrapolated to 30-45 C. Forty-five degrees above the reference, eighteen
    above the top of the range, and nothing said so.

    Only records carrying `temperature_range_c` are checked, and exactly one
    does. Of the 23 records, one states a range, ten state a single reference
    point and twelve state no number, so there is nothing else here to check
    without re-reading the primary literature -- which the registry says
    plainly was not done when it was compiled.

    A record with no range therefore passes. Absence of a domain is not a
    domain of zero, and refusing on it would turn a gap in the evidence into a
    gap in the tool.
    """
    temp = getattr(request.conditions, "temp", None)
    if temp is None:
        return

    for field, quantity in _ADDITIVE_EVIDENCE.items():
        value = getattr(request.conditions, field, 0.0) or 0.0
        if not value:
            continue
        record = evidence.get(quantity)
        span = getattr(record, "temperature_range_c", None) if record else None
        if not span:
            continue
        low, high = float(span[0]), float(span[1])
        if low <= float(temp) <= high:
            continue
        raise UnsupportedModelError(
            f"melting-temperature effect of {field} at {temp} C",
            value,
            f"the recorded coefficient covers {low:g}-{high:g} C; "
            f"run at a temperature inside it, or remove the additive",
        )
