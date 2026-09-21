"""One immutable design request, resolved once, carrying where each value came from.

Every design command used to resolve params.json for itself. `plan-pool` read
the reaction, the reach and the dimer limits; `expand-primers` read none of
them and ran at a hard-coded 3 kb with no additive correction. `design_context`
fixed the worst of that by putting the chemistry in one place, and this module
carries the idea to the rest of the configuration: the references, the
candidate source, the fixed and excluded oligos, the panel limits, the size
policy, the search budgets, the seed, and the concentration policy.

Three properties make it worth having as a type rather than as a dict.

**It is frozen, including its nested content.** A stage that could edit the
request could change the chemistry or a threshold after the search had already
used the old one, and nothing in the output would show it. Sequences arrive as
tuples for the same reason a list field in a frozen dataclass is not frozen.

**Every value records its source.** `default_sources` says, for each resolved
setting, whether the request supplied it or a default did and which default.
"What was this designed under" then has one answer to read back, rather than
being reconstructed from whichever fallback each reader happened to apply.

**It hashes canonically.** `request_hash` is a SHA-256 over a sorted JSON form,
not Python's `hash()`, which is salted per process. That is what lets a saved
result name the configuration that produced it, and what lets an equivalent
request made through the CLI and through the library be shown to be the same
request.

Configuration only. No cache, no candidate list, no optimizer. A request is
safe to build early, to log, and to write into a result file.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field, fields
from types import SimpleNamespace
from typing import Any, Dict, Mapping, Optional, Tuple

from .exceptions import InvalidDesignRequest, UnsupportedModelError

__all__ = [
    "ConcentrationPolicy",
    "DesignRequest",
    "resolve_design_request",
]


# Keys the schema no longer accepts, with what to do instead. Refused by name,
# because a retired key that is merely ignored applies a default the user did
# not choose and believes they overrode.
RETIRED_SETTINGS: Dict[str, str] = {
    "candidate_retention:legacy": (
        "candidate_retention='legacy' gave a background index to the shortlisted "
        "candidates only, so everything else scored against an empty background "
        "and read as perfectly specific. Use 'all_qc' or 'post_gini'."
    ),
    "bg_max_removal": (
        "retired with the background prefilter that deleted candidates; "
        "ordering replaced it and deletes nothing."
    ),
}

#: Accepted concentration allocation modes.
CONCENTRATION_MODES = ("per_oligo", "fixed_total")

#: Default per-oligo concentration when nothing says otherwise, in molar.
DEFAULT_PRIMER_MOLAR = 0.5e-6


@dataclass(frozen=True)
class ConcentrationPolicy:
    """How a pool's total oligo concentration is divided among its members.

    Two declared modes, because the two make different predictions and the
    difference shows up exactly when the panel size changes.

    `per_oligo` holds each oligo at `molar` and lets the total rise with the
    panel. `fixed_total` holds the total at `molar` and divides it, so adding
    an oligo lowers every other one's concentration and therefore its occupancy.

    Which one a reaction uses is a fact about the protocol, not something
    software can infer, so it is recorded rather than assumed.
    """

    mode: str
    molar: float

    def concentrations(self, oligos: Tuple[str, ...]) -> Tuple[float, ...]:
        count = len(oligos)
        if count == 0:
            return ()
        if self.mode == "fixed_total":
            share = self.molar / count
            return tuple(share for _ in range(count))
        return tuple(self.molar for _ in range(count))

    def total_molar(self, count: int) -> float:
        if self.mode == "fixed_total":
            return self.molar
        return self.molar * count


@dataclass(frozen=True)
class DesignRequest:
    """Everything a design run is configured by, and nothing it computes."""

    # References
    fg_prefixes: Tuple[str, ...]
    fg_seq_lengths: Tuple[int, ...]
    bg_prefixes: Tuple[str, ...]
    bg_seq_lengths: Tuple[int, ...]
    fg_genomes: Tuple[str, ...]
    bg_genomes: Tuple[str, ...]
    fg_circular: bool

    # Chemistry
    polymerase: str
    conditions: Any
    coverage_reach: int
    min_tm: float
    max_tm: float

    # Candidate source identity
    data_dir: str
    primer_lengths: Tuple[int, ...]
    candidate_retention: str

    # Panel composition
    fixed_oligos: Tuple[str, ...]
    excluded_oligos: Tuple[str, ...]
    target_size: int
    max_sets: int

    # Hard limits
    max_dimer_bp: int
    max_self_dimer_bp: int
    max_dimer_dg: Optional[float]
    constraints: Any

    # Search
    total_search_evaluations: Optional[int]
    total_search_seconds: Optional[float]
    max_frontier_refills: int
    optimization_method: Optional[str]
    seed: Optional[int]

    # Models and chemistry policy
    concentration_policy: ConcentrationPolicy
    model_versions: Tuple[Tuple[str, str], ...]

    #: Where each resolved setting came from: "request", or a named default.
    default_sources: Mapping[str, str] = field(default_factory=dict)

    # ---------------------------------------------------------------- hashing

    def to_dict(self) -> Dict[str, Any]:
        """A JSON-serializable form, including the hash it produces."""
        payload = {name: _plain(getattr(self, name)) for name in _HASHED_FIELDS}
        payload["default_sources"] = dict(sorted(self.default_sources.items()))
        payload["request_hash"] = self.request_hash
        return payload

    @property
    def request_hash(self) -> str:
        """SHA-256 over a canonical form of every setting that affects a design.

        `default_sources` is excluded: it records how a value was reached, not
        what it is, so two requests reaching the same chemistry by different
        routes describe the same design and must compare equal.
        """
        payload = {name: _plain(getattr(self, name)) for name in _HASHED_FIELDS}
        canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    # ----------------------------------------------------------- concentration

    @property
    def total_primer_molar(self) -> float:
        """Total oligo concentration for the panel size this request asks for."""
        return self.concentration_policy.total_molar(self.target_size)

    def concentrations_molar(self, oligos: Tuple[str, ...]) -> Tuple[float, ...]:
        """Each oligo's concentration under the declared allocation mode."""
        return self.concentration_policy.concentrations(tuple(oligos))

    # ----------------------------------------------------------------- helpers

    def reference_manifest(self) -> Dict[str, str]:
        """Each prefix paired with the genome its index must have been built from.

        Paired positionally, which is how every other part of this codebase
        relates the two lists, and only where both are present. A prefix with no
        genome is absent from the manifest rather than paired with a guess: the
        identity check then does not run for it, which is honest, while pairing
        it with whatever happened to be in a mutable global is not. That pairing
        is what made a design refuse its own index under `pytest -n 8`.
        """
        manifest: Dict[str, str] = {}
        for prefixes, genomes in (
            (self.fg_prefixes, self.fg_genomes),
            (self.bg_prefixes, self.bg_genomes),
        ):
            for prefix, genome in zip(prefixes, genomes):
                if prefix and genome:
                    manifest[str(prefix)] = str(genome)
        return manifest

    def design_context(self):
        """The chemistry subset, for the code that already takes a context."""
        from .design_context import DesignContext

        return DesignContext(
            conditions=self.conditions,
            coverage_reach=self.coverage_reach,
            max_dimer_bp=self.max_dimer_bp,
            max_self_dimer_bp=self.max_self_dimer_bp,
            min_tm=self.min_tm,
            max_tm=self.max_tm,
            fg_circular=self.fg_circular,
            polymerase=self.polymerase,
            max_dimer_dg=self.max_dimer_dg,
            constraints=self.constraints,
        )


#: Fields that identify the design. `conditions` is folded in through its own
#: fingerprint rather than by serialising the object, so an unrelated change to
#: that class does not move every stored hash.
_HASHED_FIELDS = tuple(
    name for name in (f.name for f in fields(DesignRequest)) if name not in {"default_sources"}
)


def _plain(value):
    """A JSON-safe, order-stable form of one field."""
    if hasattr(value, "fingerprint"):
        return {"fingerprint": value.fingerprint()}
    if isinstance(value, ConcentrationPolicy):
        return {"mode": value.mode, "molar": value.molar}
    if isinstance(value, tuple):
        return [_plain(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "__dict__"):
        return {key: _plain(val) for key, val in sorted(vars(value).items())}
    return str(value)


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------


def _schema_properties() -> set:
    from .schema import load_schema

    return set(load_schema().get("properties", {}))


#: Keys a design accepts that the params schema does not declare, because they
#: describe the request rather than the pipeline configuration.
_REQUEST_ONLY_KEYS = frozenset(
    {
        "fixed_oligos",
        "excluded_oligos",
        "concentration_mode",
        "total_primer_molar",
        "seed",
    }
)


def _require_finite(params: Mapping[str, Any]) -> None:
    """A NaN threshold compares False against everything and rejects nothing."""
    for key, value in params.items():
        if isinstance(value, float) and not math.isfinite(value):
            raise InvalidDesignRequest(key, "must be finite", value)
        if isinstance(value, (list, tuple)):
            for item in value:
                if isinstance(item, float) and not math.isfinite(item):
                    raise InvalidDesignRequest(key, "must contain only finite values", item)


def _reject_unknown_and_retired(params: Mapping[str, Any]) -> None:
    for key, value in params.items():
        retired = RETIRED_SETTINGS.get(f"{key}:{value}") or RETIRED_SETTINGS.get(key)
        if retired:
            raise InvalidDesignRequest(key, retired, value)

    declared = _schema_properties() | _REQUEST_ONLY_KEYS
    # A leading underscore marks a comment. Three shipped example configs use
    # `_comment`, and a convention the repository's own files rely on is part
    # of the accepted input rather than a typo to refuse.
    supplied = {key for key in params if not str(key).startswith("_")}
    unknown = sorted(supplied - declared)
    if unknown:
        raise InvalidDesignRequest(
            unknown[0],
            "is not a setting this version accepts"
            + (f" (also: {', '.join(unknown[1:])})" if len(unknown) > 1 else ""),
        )


def _resolve_reach(params: Mapping[str, Any], polymerase: str, sources: Dict[str, str]) -> int:
    """The selection reach, refusing an explicit zero rather than replacing it.

    `design_context_from_params` used `override or params.get(...)`, so an
    explicit 0 was falsy and silently became the polymerase default. Zero is
    not a reach and the user who wrote it meant something; either way, applying
    3 kb instead and reporting coverage at it is the outcome to avoid.
    """
    from .coverage import polymerase_extension_reach

    if "coverage_reach" in params and params["coverage_reach"] is not None:
        value = params["coverage_reach"]
        if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
            raise InvalidDesignRequest(
                "coverage_reach", "must be a positive whole number of base pairs", value
            )
        sources["coverage_reach"] = "request"
        return int(value)

    sources["coverage_reach"] = f"polymerase:{polymerase}"
    return int(polymerase_extension_reach(polymerase, coverage_metric="realistic"))


def _resolve_concentration(
    params: Mapping[str, Any], sources: Dict[str, str]
) -> ConcentrationPolicy:
    mode = params.get("concentration_mode")
    if mode is None:
        mode = "per_oligo"
        sources["concentration_mode"] = "default:per_oligo"
    else:
        sources["concentration_mode"] = "request"
    if mode not in CONCENTRATION_MODES:
        raise InvalidDesignRequest(
            "concentration_mode", f"must be one of {', '.join(CONCENTRATION_MODES)}", mode
        )

    key = "total_primer_molar" if mode == "fixed_total" else "primer_conc"
    value = params.get(key)
    if value is None:
        value = DEFAULT_PRIMER_MOLAR * (
            params.get("num_primers", 6) if mode == "fixed_total" else 1
        )
        sources[key] = f"default:{value}"
    else:
        sources[key] = "request"
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise InvalidDesignRequest(key, "must be a positive concentration in molar", value)
    return ConcentrationPolicy(mode=mode, molar=value)


def _resolve_budgets(params: Mapping[str, Any], sources: Dict[str, str]) -> Dict[str, Any]:
    from .search_control import resolve_search_settings

    for key in ("total_search_evaluations", "total_search_seconds", "max_frontier_refills"):
        value = params.get(key)
        if value is not None and (not isinstance(value, (int, float)) or value < 0):
            raise InvalidDesignRequest(key, "must be non-negative", value)
        sources[key] = "request" if value is not None else "default"
    return resolve_search_settings(dict(params))


def _resolve_oligos(params: Mapping[str, Any]) -> Tuple[Tuple[str, ...], Tuple[str, ...]]:
    fixed = tuple(str(p).upper() for p in (params.get("fixed_oligos") or ()))
    excluded = tuple(str(p).upper() for p in (params.get("excluded_oligos") or ()))
    both = sorted(set(fixed) & set(excluded))
    if both:
        raise InvalidDesignRequest(
            "fixed_oligos/excluded_oligos",
            f"the same oligo cannot be both required and forbidden: {', '.join(both)}",
        )
    return fixed, excluded


def _resolve_constraints(params: Mapping[str, Any], has_background: bool):
    """Panel limits, refusing a background-measured one with no background.

    Reported as satisfied is the failure this prevents: the quantity cannot be
    measured at all without a host reference, and a limit on an unmeasurable
    quantity either passes vacuously or fails every panel. Neither is what the
    user asked for, and the remedy is to supply the reference.
    """
    from .panel_acceptance import constraints_from_parameter
    from .pool_objective import _LIMITS

    if not has_background:
        for field_name, _metric, _sense, _message, needs_background in _LIMITS:
            if needs_background and params.get(field_name) is not None:
                raise InvalidDesignRequest(
                    field_name,
                    "is measured against a background genome, and this request "
                    "names none. Supply bg_genomes/bg_prefixes, or remove the limit.",
                    params[field_name],
                )
    return constraints_from_parameter(SimpleNamespace(**dict(params)))


def _require_references(params: Mapping[str, Any]) -> None:
    """A design needs a target, named either as a genome or as a k-mer prefix.

    Both are accepted because both name the same reference at different stages:
    `fg_genomes` is what a user writes, and `fg_prefixes` is what step 1 derives
    from it. Requiring the derived form would refuse a valid params file purely
    for being read before the step that derives it.
    """
    if not (params.get("fg_prefixes") or params.get("fg_genomes")):
        raise InvalidDesignRequest(
            "fg_prefixes",
            "a design needs at least one target reference (fg_genomes or fg_prefixes)",
        )


def resolve_design_request(params: Mapping[str, Any]) -> DesignRequest:
    """Resolve a params mapping into one frozen request, or refuse it by name.

    Every refusal names the field, so the remedy is in the message rather than
    in a traceback. Nothing here reads a `parameter` module global: a request
    is a pure function of the mapping it is given, which is what lets the same
    request be built from a CLI invocation and from a library call and shown to
    be the same request.
    """
    params = dict(params)
    _reject_unknown_and_retired(params)
    _require_finite(params)

    sources: Dict[str, str] = {}
    polymerase = str(params.get("polymerase") or "phi29")
    sources["polymerase"] = "request" if params.get("polymerase") else "default:phi29"

    # Before the reference check, so a nonsense value in the request is named
    # ahead of something the request merely omits. Raises `UnsupportedModelError`
    # for an enzyme with no recorded reach, which is a different remedy from a
    # bad field: the name may be spelled correctly and simply not be modelled.
    reach = _resolve_reach(params, polymerase, sources)
    _require_references(params)

    from .reaction_conditions import build_reaction_conditions

    try:
        conditions = build_reaction_conditions(SimpleNamespace(**params))
    except (KeyError, LookupError) as exc:
        raise UnsupportedModelError("reaction conditions", polymerase) from exc
    except ValueError as exc:
        raise InvalidDesignRequest("reaction conditions", str(exc)) from exc

    fixed, excluded = _resolve_oligos(params)
    budgets = _resolve_budgets(params, sources)
    has_background = bool(params.get("bg_prefixes"))
    constraints = _resolve_constraints(params, has_background)
    policy = _resolve_concentration(params, sources)

    min_k = int(params.get("min_k", 6))
    max_k = int(params.get("max_k", 12))
    if min_k > max_k:
        raise InvalidDesignRequest("min_k", f"is above max_k ({max_k})", min_k)

    target_size = int(params.get("target_set_size") or params.get("num_primers") or 6)

    request = DesignRequest(
        fg_prefixes=tuple(str(p) for p in params.get("fg_prefixes") or ()),
        fg_seq_lengths=tuple(int(n) for n in params.get("fg_seq_lengths") or ()),
        bg_prefixes=tuple(str(p) for p in params.get("bg_prefixes") or ()),
        bg_seq_lengths=tuple(int(n) for n in params.get("bg_seq_lengths") or ()),
        fg_genomes=tuple(str(g) for g in params.get("fg_genomes") or ()),
        bg_genomes=tuple(str(g) for g in params.get("bg_genomes") or ()),
        fg_circular=bool(params.get("fg_circular", False)),
        polymerase=polymerase,
        conditions=conditions,
        coverage_reach=reach,
        min_tm=float(params.get("min_tm", 20)),
        max_tm=float(params.get("max_tm", 50)),
        data_dir=str(params.get("data_dir") or "."),
        primer_lengths=tuple(range(min_k, max_k + 1)),
        candidate_retention=str(params.get("candidate_retention") or "all_qc"),
        fixed_oligos=fixed,
        excluded_oligos=excluded,
        target_size=target_size,
        max_sets=int(params.get("max_sets", 5)),
        max_dimer_bp=int(params.get("max_dimer_bp", 3)),
        max_self_dimer_bp=int(params.get("max_self_dimer_bp", 4)),
        max_dimer_dg=params.get("max_dimer_dg"),
        constraints=constraints,
        total_search_evaluations=budgets["total_search_evaluations"],
        total_search_seconds=budgets["total_search_seconds"],
        max_frontier_refills=budgets["max_frontier_refills"],
        optimization_method=params.get("optimization_method"),
        seed=params.get("seed"),
        concentration_policy=policy,
        model_versions=_model_versions(),
        default_sources=dict(sorted(sources.items())),
    )

    # Last, because it reads the resolved request rather than the mapping. A
    # computation outside a recorded domain is refused here, before any index
    # is opened: the alternative is a number produced with no evidence behind
    # it, which looks exactly like one produced with evidence.
    from .model_evidence import require_model_support

    require_model_support(request)
    return request


def _model_versions() -> Tuple[Tuple[str, str], ...]:
    """Identifiers for the models a design's numbers depend on.

    Part of the request so a saved result names the code that produced it.
    Task 4 replaces the thermodynamic entry with a versioned evidence record;
    until then this states the parameter set by name rather than implying that
    its domain has been checked.
    """
    from neoswga import __version__

    return (
        ("neoswga", str(__version__)),
        ("nearest_neighbour", "santalucia-1998-unified"),
        ("salt_correction", "owczarzy-2008"),
    )


def params_from_parameter_module(parameter) -> Dict[str, Any]:
    """Project the resolved `parameter` module onto the settings a request takes.

    The module carries far more than configuration -- loaders, caches, helpers
    and the globals `_apply_params_only_keys` assigns -- so handing `vars()`
    straight to the resolver would report every one of them as an unknown
    setting. Projecting onto the declared keys is what makes the resolver's
    strictness usable against a module that was never a params file.

    This is the seam, not the destination. The plan's intent is that commands
    hold a request rather than reading module globals at all; until each reader
    moves, this builds the same request from what the module resolved to, so
    the hash and the recorded sources describe the run that actually happened.
    """
    accepted = _schema_properties() | _REQUEST_ONLY_KEYS
    params = {}
    for key in sorted(accepted):
        value = getattr(parameter, key, None)
        if value is not None:
            params[key] = value
    return params


def design_request_from_parameter_module(parameter) -> DesignRequest:
    """Resolve and validate the request this run is about to execute."""
    return resolve_design_request(params_from_parameter_module(parameter))


def design_request_for_run(args, parameter) -> Optional[DesignRequest]:
    """The request a command is about to execute, resolved before it starts.

    Read from the params FILE rather than from the `parameter` module, because
    `run_step4` resolves that module inside `optimize_step4` -- the same
    ordering trap `warn_on_condition_drift` documents, where every reaction
    global still holds its default at this point. The file is the request; the
    module is one resolution of it.

    Returns None when the command was given no params file, which is the
    library and `--candidates` path: the caller supplied the configuration
    directly and there is no file to validate. That is an absence, not a
    silently accepted request, and the caller can still build one itself.
    """
    path = getattr(args, "json_file", None)
    if not path:
        return None
    with open(path) as handle:
        supplied = json.load(handle)
    return resolve_design_request(supplied)
