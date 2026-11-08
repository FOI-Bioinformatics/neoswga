"""What a panel does to each host genome on its own, and the reduction over them.

`core/background_aggregate.py` says what `background_aggregate` MEANS;
`core/background_gate.py` applies it in step 2. This module is step 4's half:
it turns the per-prefix background positions `BaseOptimizer.compute_metrics`
already holds into one record per host, and reduces those records to the
figures the configured panel limits are evaluated against.

Why it is a separate module. `base_optimizer.py` is at its size budget, and the
per-prefix loop there produces ONE pooled background figure -- a flat sum of
site counts over a single pooled background length -- which is the `sum`
reading and the only one it could express. Writing the second reading into that
loop is both out of budget and the wrong shape: the quantity is per host, and
reductions over hosts are what differ between the two readings.

**The figures, per host.** Sites, occupancy-weighted load, coverage, and
`selectivity_density` against that host's own length. The names and the
arithmetic are `core/reference_panel_evaluation.py`'s -- `sites`,
`sites_per_mb`, `coverage`, `selectivity_density`, and
`worst_host_selectivity_density` for the reduction -- because that module
already measures exactly this after the fact, and two definitions of one
quantity is how they come to disagree. `selectivity_ratio` is deliberately NOT
among them: it carries no length, so it moves with how much background sequence
a host happens to be, and a per-host table of it would invite precisely the
comparison Known Issue 6 records as invalid. Density is the comparable figure.

**The reduction, per limit.** Under `worst-case` the panel faces ONE host out
of the listed set and which one is not known in advance, so a limit holds only
if it holds for every host. Each limit therefore reduces independently:
`min_selectivity_density` against the LOWEST per-host density,
`max_background_sites` against the HIGHEST per-host site count,
`max_host_coverage` against the HIGHEST per-host coverage. There is no new limit
key and no per-host threshold list -- one limit, applied per host, which is the
Phase 7 guard.

**One dispatch.** `limit_figure` is the only place anything decides which figure
a limit on one of those three quantities reads, and it is asked by everything
that DECIDES one (`pool_objective.violations` and `.shortfall`) and everything
that DISPLAYS one (`panel_acceptance`'s configured-limits table,
`panel_regime`'s criteria, `pool_planner`'s plan row, `panel_beam`'s ranking).
That second half was missed once, and the symptom was a printed table reading
`min_selectivity_density 51 10 NO` -- a value above its floor marked not met --
and a saved `panel_regime` slack wrong by 79x. A value and a verdict that come
from different figures is the defect this module exists to prevent, one layer up
from the one `PoolObjective` prevents.

**Unknown is not the worst of the rest.** A host whose positions were never
scanned has no answer, and a reduction over a set with an unmeasured member is
unmeasured: the figure is `None`, never the minimum or maximum over the hosts
that happened to answer. `pool_objective._unmet` treats `None` as missing its
limit, so a background-measured limit under `worst-case` with an unscanned host
is REFUSED rather than satisfied. `position_cache.has_entry` is the question
asked, because an empty entry (scanned, binds nowhere) and an absent one (never
scanned) are different answers and only the first is a measurement. It is asked
at every host count, one included.

**A figure nothing reads is not computed, and says so.** The per-host
occupancy-weighted load costs one k-mer table pass per host per evaluation and
only `worst-case` reads it. `per_host_load_is_in_force` is the one reader of
that question; a row whose load was skipped carries `load_unavailable` beside
real measured sites and coverage. Three states, never two: measured, not
measurable, and measured-with-one-figure-withheld.

**The search-time host term is pooled, and that is refused rather than
silent.** `host_site_charge` has one reading, and
`refuse_uncalibrated_host_term` explains why: the per-host reduction multiplies
the charge by `sum(host lengths) / worst host length`, and `background_weight`
is calibrated against the pooled magnitude. `worst-case` with
`background-aware` and two or more hosts is refused up front, naming the
measurement that would settle it.

**Nothing moves under `sum`.** `limit_figure` returns the pooled attribute it
has always returned, the pooled fields on `PrimerSetMetrics` keep their pooled
meaning under every aggregate, and with fewer than two hosts the per-host
reading IS the pooled one, so `worst-case` and `sum` agree exactly rather than
approximately. That last point is why the one-host case short-circuits instead
of being recomputed: `_coverage_over_prefixes` returns `(fraction * length) /
length`, which is not bit-identical to `fraction`, and "exactly" is the claim
Phase 7c's verification makes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .background_aggregate import SUM, WORST_CASE, resolve_background_aggregate
from .selectivity import selectivity_density_from_loads

__all__ = [
    "AGGREGATE_REDUCED",
    "BackgroundReading",
    "HostFigures",
    "host_site_charge",
    "limit_figure",
    "measure_background",
    "per_host_load_is_in_force",
    "refuse_uncalibrated_host_term",
]

#: Why a per-host occupancy load was not computed, when nothing reads it. A
#: DIFFERENT thing from a host that could not be measured, and kept apart from
#: it: "scanned, and this run does not need the figure" is not "never scanned",
#: and neither is zero.
LOAD_NOT_IN_FORCE = (
    "the per-host occupancy load is computed only under "
    "background_aggregate='worst-case', where a panel limit reads it. The "
    "pooled load for this panel is in effective_bg_sites"
)

#: The three metrics whose value a limit is read against changes with
#: `background_aggregate`, mapped to the `PrimerSetMetrics` field carrying the
#: per-host reduction. Everything not here reads its pooled attribute under
#: every aggregate, which is what keeps `sum` byte-identical.
#:
#: These are exactly the three `pool_objective._LIMITS` rows marked as needing a
#: background. Adding a fourth background-measured limit means adding its
#: reduction here, and `tests/test_per_host_figures.py` pins that correspondence
#: so a new row cannot quietly be evaluated on the pooled figure.
AGGREGATE_REDUCED: dict[str, str] = {
    "selectivity_density": "worst_host_selectivity_density",
    "total_bg_sites": "worst_host_background_sites",
    "bg_coverage": "worst_host_bg_coverage",
}


@dataclass(frozen=True)
class HostFigures:
    """What one host genome measured for one panel, or why it did not.

    Three states, not two, and they must not share a representation.

    - Measured. `unavailable` and `load_unavailable` are empty. A host that was
      scanned and carries no site has `sites == 0`: that zero is a result.
    - Not measurable. `unavailable` says why and EVERY figure is `None`. The
      host's positions were never scanned, so nothing about it is known.
    - Measured, but the occupancy-weighted load was not computed because
      nothing reads it. `load_unavailable` says so, `effective_load` and
      `selectivity_density` are `None`, and sites, sites per Mb and coverage are
      real measurements -- they come free from positions already in hand. This is
      the default-path case; see `per_host_load_is_in_force`.
    """

    prefix: str
    length: int
    sites: int | None
    sites_per_mb: float | None
    effective_load: float | None
    coverage: float | None
    selectivity_density: float | None
    unavailable: str = ""
    #: Why the occupancy-weighted load and the density derived from it are
    #: absent although the host WAS measured. Empty when they are present.
    load_unavailable: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "length_bp": self.length,
            "sites": self.sites,
            "sites_per_mb": self.sites_per_mb,
            "effective_site_load": self.effective_load,
            "coverage": self.coverage,
            # The comparable figure. Named the same as in
            # `reference_panel_evaluation.PairRecord`, and computed the same way.
            "selectivity_density": self.selectivity_density,
            "unavailable": self.unavailable or None,
            # Separate from the above: this host WAS measured and the load alone
            # was not computed. A reader that collapsed the two would report a
            # scanned host as unscanned.
            "load_unavailable": self.load_unavailable or None,
            # Said in every record rather than inferred from its absence: a
            # per-host ratio is not comparable across hosts of different size.
            "comparable": "selectivity_density",
        }


@dataclass(frozen=True)
class BackgroundReading:
    """The per-host table and the reductions a limit is read against."""

    mode: str
    per_host: tuple[HostFigures, ...] = ()
    #: Lowest per-host density, and the host that holds it. `None` when any host
    #: went unmeasured, never the lowest of the rest.
    worst_host_selectivity_density: float | None = None
    worst_host: str | None = None
    #: Highest per-host site count and coverage. Reduced per limit, because a
    #: panel clean on the host with the worst density may still be the densest
    #: on another, and `worst-case` means every host must clear every limit.
    worst_host_background_sites: int | None = None
    worst_host_bg_coverage: float | None = None
    unavailable: str = ""

    def as_dict(self) -> dict[str, Any]:
        return {
            "background_aggregate": self.mode,
            "per_host": {h.prefix: h.as_dict() for h in self.per_host},
            "worst_host": self.worst_host,
            "worst_host_selectivity_density": self.worst_host_selectivity_density,
            "worst_host_background_sites": self.worst_host_background_sites,
            "worst_host_bg_coverage": self.worst_host_bg_coverage,
            "unavailable": self.unavailable or None,
        }

    def kwargs(self) -> dict[str, Any]:
        """The fields `PrimerSetMetrics` and `PoolMetrics` carry, as kwargs.

        One place, so the full evaluator and the focused one cannot be given
        different subsets of the same reading.
        """
        return {
            "background_aggregate": self.mode,
            "per_host": {h.prefix: h.as_dict() for h in self.per_host},
            "worst_host": self.worst_host,
            "worst_host_selectivity_density": self.worst_host_selectivity_density,
            "worst_host_background_sites": self.worst_host_background_sites,
            "worst_host_bg_coverage": self.worst_host_bg_coverage,
        }


def limit_figure(metrics: Any, metric: str) -> Any:
    """The value a panel limit on `metric` is evaluated against.

    This is the one dispatch on `background_aggregate` in step 4's acceptance
    path, and `pool_objective.violations` / `pool_objective.shortfall` are its
    only callers -- so the two cannot disagree about which figure a limit bound
    on, which is the defect shape `PoolObjective` exists to prevent.

    Under `sum`, and for any metric that is not background-measured, this is the
    plain attribute lookup those two performed before: the pooled figure, with
    `None` for an evaluator that does not carry the quantity at all. The pooled
    fields never change meaning, so no output has to be re-read.

    Under `worst-case` the three background-measured metrics read their per-host
    reduction instead. An evaluator with no reduction recorded answers `None`,
    which `_unmet` treats as missing the limit: a limit must not pass because
    nothing supplied the figure it is about.
    """
    aggregate = getattr(metrics, "background_aggregate", None) or SUM
    if aggregate == WORST_CASE and metric in AGGREGATE_REDUCED:
        return getattr(metrics, AGGREGATE_REDUCED[metric], None)
    return getattr(metrics, metric, None)


def per_host_load_is_in_force(aggregate: str | None) -> bool:
    """Whether anything will READ a per-host occupancy-weighted load.

    The one reader of that question, so the module that computes the load and
    the rows that report its absence cannot disagree about whether it was in
    force. True exactly under `worst-case`, where `limit_figure` routes the
    density floor to the per-host reduction. Under `sum` every configured limit
    reads a pooled field, so computing the per-host loads would cost one k-mer
    table pass per host per evaluation for a figure nothing consults.

    Measured on a real two-host design with `selectivity_mode = occupancy`,
    computing them unconditionally cost +58 percent of `compute_metrics`
    (2.19 ms to 3.46 ms, two `weighted_site_load` calls to four) on the DEFAULT
    path. `compute_metrics` is the search's per-panel evaluator, so that is a
    multiplier on step-4 wall time rather than a one-off.

    The sites, the sites per Mb and the coverage are not gated on this: they
    come from the per-prefix positions the caller already holds and cost
    nothing, so the per-host table is still a table under `sum`. What it does
    not carry there is the load and the density, each marked
    `load_unavailable` with `LOAD_NOT_IN_FORCE` rather than zero.

    One exception, and it is free: under `selectivity_mode == "exact"` the
    pooled load IS the raw site count, so each host's load is its own count and
    no extra work exists to skip.
    """
    return resolve_background_aggregate(aggregate) == WORST_CASE


def refuse_uncalibrated_host_term(
    method: Any, aggregate: str | None, host_count: int, *, ensemble_methods: Any = None
) -> None:
    """Refuse `worst-case` with the one optimizer whose host term it miscalibrates.

    `background-aware` ranks removals on
    `background_weight * host_sites - coverage_loss * 1000`. The host term is a
    site COUNT, and `background_weight` (default 2.0) was chosen against the
    magnitude a POOLED count has. A per-host reduction cannot keep both the
    ranking and the magnitude:

    - rank on the worst host's DENSITY, carried back onto the pooled length, and
      the ranking is right -- counts are not comparable across hosts of
      different size (Known Issue 6) -- but the magnitude is multiplied by
      `sum(host lengths) / worst host length`. On this project's own three-host
      panel (3.35 Mb, 144 Mb, 3300 Mb) that factor is 1029, 24 or 1.05
      depending on which host is worst, and 44 in the Phase 7c test world. A
      `background_weight` tuned under `sum` then buys up to a thousand times
      more host avoidance than it did, and coverage is spent paying for it.
    - rank on the worst host's own count and the magnitude is bounded by the
      pooled sum, but the ranking compares counts across unequal hosts, which is
      the comparison Known Issue 6 says does not hold.

    Phase 7c shipped the first and claimed the magnitude was preserved. It is
    not, past one host. Rather than leave a term mis-scaled by three orders of
    magnitude reachable behind a note, the combination is refused.

    What would settle it: a measurement of `background_weight` under
    `worst-case` on a host panel spanning at least two orders of magnitude --
    Phase 3's D3 design rerun with `background_aggregate: worst-case` and
    `optimization_method: background-aware` across a range of
    `background_weight`, reporting delivered coverage against worst-host
    density. `scripts/benchmarking/diversity_baseline.py` is where that belongs.

    Every other method is unaffected: none of them has a host term in its
    selection, and the per-host figures and the configured limits read
    `limit_figure`, which carries no weight to calibrate. With ONE host the two
    readings are integer-identical, so there is nothing to refuse and nothing is.

    `run_optimization` calls this before any optimizer is built, because
    `_dispatch_optimizer` converts a constructor exception into a failure
    result, and a refusal must not become one.
    """
    from .exceptions import InvalidDesignRequest

    if resolve_background_aggregate(aggregate) != WORST_CASE or int(host_count) < 2:
        return
    named = {str(method or "").lower()}
    named.update(str(name).lower() for name in (ensemble_methods or ()))
    # An ensemble with no explicit method list runs `_dispatch_optimizer`'s
    # default, which includes `background-aware`. Naming the aliases here rather
    # than importing that list keeps this module free of a cycle; the ensemble
    # aliases are pinned by `tests/test_per_host_figures.py` against the
    # dispatch's own default so the two cannot drift apart.
    if ensemble_methods is None and named & {"ensemble", "auto", "all"}:
        named.add("background-aware")
    if not named & {"background-aware", "bg-aware", "clinical"}:
        return
    raise InvalidDesignRequest(
        "background_aggregate",
        "cannot be 'worst-case' with optimization_method 'background-aware' and "
        f"{int(host_count)} host genomes. That optimizer ranks removals on a host "
        "SITE COUNT scaled by background_weight, and the per-host reduction "
        "multiplies that count by sum(host lengths) / worst host length -- up to "
        "1029x on this project's own three-host panel -- so a background_weight "
        "tuned under 'sum' no longer means what it meant. No measurement "
        "establishes the weight under the per-host reading. Use 'sum' with "
        "'background-aware', or 'worst-case' with another method (the per-host "
        "figures and every configured panel limit still bind per host there). "
        "What would settle it: Phase 3's D3 design rerun under 'worst-case' "
        "across a range of background_weight, reporting delivered coverage "
        "against worst-host density",
        WORST_CASE,
    )


def host_site_charge(sites_by_prefix: dict[str, int], lengths: Any, aggregate: str | None) -> int:
    """The host load a SEARCH charges a panel, in site-count units.

    `hybrid_optimizer._count_background_sites` is the caller, and through it both
    background-aware stages: `_prune_background`'s ranking and
    `_network_refine`'s `bg_by_primer` weighting. The reduction over hosts lives
    here rather than there so the stage that selects and the metrics that report
    cannot disagree about what "the background" is -- the two halves of Known
    Issue 14.

    **There is one reading, and it is the pooled sum: every listed host present
    at once.** Phase 7c shipped a `worst-case` branch here that charged the worst
    host's site density carried back onto the pooled length. The ranking was
    right and the magnitude was not: it multiplies the charge by
    `sum(host lengths) / worst host length`, which is 1029 on this project's own
    three-host panel, and `background_weight` scales this term against a coverage
    term. `refuse_uncalibrated_host_term` carries the whole argument and the
    measurement that would settle it.

    So `aggregate` is not ignored and this is not a silent pooling: `worst-case`
    with two or more hosts RAISES here, behind the up-front refusal, so the
    combination cannot be reached by a route that skipped the gate. With one host
    the two readings are integer-identical and there is nothing to decide.
    """
    prefixes = list(sites_by_prefix)
    total = sum(int(value) for value in sites_by_prefix.values())
    if resolve_background_aggregate(aggregate) == WORST_CASE and len(prefixes) >= 2:
        refuse_uncalibrated_host_term("background-aware", aggregate, len(prefixes))
    return total


def measure_background(
    optimizer: Any,
    primers: list[str],
    positions_by_prefix: dict[str, Any],
    *,
    pooled_sites: int,
    pooled_coverage: float | None,
    pooled_density: float | None,
    effective_fg: float,
    effective_bg: float,
    selectivity_mode: str,
) -> BackgroundReading:
    """One record per host, plus the reductions, from positions already gathered.

    `positions_by_prefix` is `compute_metrics`' own per-prefix background map,
    passed in rather than re-fetched: `get_positions` is the hottest method in
    the package and a second pass over it would double the cost of every
    evaluation to recompute numbers that are already in hand.

    With fewer than two hosts the pooled figures ARE the per-host figures, so
    they are carried over rather than recomputed. That is what makes
    `worst-case` with one host equal `sum` with one host exactly; see the module
    docstring.

    The occupancy-weighted per-host load is the one quantity positions cannot
    supply. When `selectivity_mode` is "exact" the pooled load is the raw site
    count, so each host's load is its own count and nothing extra is computed.
    Otherwise it is computed only when `per_host_load_is_in_force` says something
    will read it, and a row whose load was skipped says so in
    `load_unavailable` -- a measured host with one figure withheld, which is not
    the same as a host nothing measured and not the same as zero.
    """
    mode = resolve_background_aggregate(getattr(optimizer, "background_aggregate", None))
    prefixes = list(getattr(optimizer, "bg_prefixes", None) or [])
    lengths = [int(value) for value in (getattr(optimizer, "bg_seq_lengths", None) or [])]

    if not prefixes or len(prefixes) != len(lengths):
        # No host, or prefixes and lengths that cannot be paired one to one.
        # Pairing a prefix with the wrong length is the defect
        # `tests/test_pairing_a_prefix_with_its_length_cannot_truncate.py` was
        # written for; refusing to guess is the only safe answer here, and the
        # pooled figures remain exactly what they were.
        reason = "" if not prefixes else "bg_prefixes and bg_seq_lengths do not pair one to one"
        return BackgroundReading(mode=mode, unavailable=reason)

    # Asked for EVERY host count including one. The short-circuit below used to
    # come first, so at n=1 the module's own rule -- an absent entry is not a
    # zero -- went unasked, and an indexed host holding no entry for any panel
    # oligo reported `sites: 0` with `unavailable: null` and the most favourable
    # density there is. `sum` answers the same way and always did, but
    # `worst-case` promises that unknown is refused, and a promise void at one
    # host is not the promise.
    unscanned = _unscanned_hosts(optimizer, primers, prefixes)

    if len(prefixes) == 1:
        if prefixes[0] in unscanned:
            return _reduce(mode, (_unmeasured(prefixes[0], lengths[0], unscanned[prefixes[0]]),))
        # The pooled figures ARE this host's, so they are carried over rather
        # than recomputed. `pooled_coverage` is None only for an evaluator that
        # does not measure background coverage at all (`pool_metrics`), and then
        # this host's own coverage is measured here rather than left absent.
        coverage = pooled_coverage
        if coverage is None:
            coverage = _coverage_of(optimizer, positions_by_prefix, prefixes[0], lengths[0])
        one = HostFigures(
            prefix=prefixes[0],
            length=lengths[0],
            sites=int(pooled_sites),
            sites_per_mb=_per_mb(pooled_sites, lengths[0]),
            effective_load=float(effective_bg),
            coverage=coverage,
            selectivity_density=pooled_density,
        )
        return BackgroundReading(
            mode=mode,
            per_host=(one,),
            worst_host_selectivity_density=pooled_density,
            worst_host=prefixes[0],
            worst_host_background_sites=int(pooled_sites),
            worst_host_bg_coverage=coverage,
        )

    loads, load_skipped = _per_host_loads(
        optimizer, primers, prefixes, positions_by_prefix, selectivity_mode, mode
    )
    fg_length = int(sum(getattr(optimizer, "fg_seq_lengths", None) or []))

    records = []
    for prefix, length in zip(prefixes, lengths, strict=True):
        if prefix in unscanned:
            records.append(_unmeasured(prefix, length, unscanned[prefix]))
            continue
        sites = len(positions_by_prefix.get(prefix, ()) or ())
        load = loads.get(prefix)
        if load is None and not load_skipped:
            records.append(
                _unmeasured(
                    prefix,
                    length,
                    "the occupancy-weighted load for this host could not be computed",
                )
            )
            continue
        records.append(
            HostFigures(
                prefix=prefix,
                length=length,
                sites=sites,
                sites_per_mb=_per_mb(sites, length),
                effective_load=load,
                coverage=_coverage_of(optimizer, positions_by_prefix, prefix, length),
                selectivity_density=(
                    None
                    if load is None
                    else selectivity_density_from_loads(effective_fg, fg_length, load, length)
                ),
                load_unavailable=LOAD_NOT_IN_FORCE if load is None else "",
            )
        )

    return _reduce(mode, tuple(records))


def _reduce(mode: str, records: tuple[HostFigures, ...]) -> BackgroundReading:
    """The three reductions, each `None` when any host went unmeasured.

    Per quantity rather than per host: `worst-case` says the panel must satisfy
    every limit whichever host turns up, and the host with the lowest density
    need not be the one with the most sites.
    """
    missing = [r.prefix for r in records if r.unavailable]
    if missing:
        return BackgroundReading(
            mode=mode,
            per_host=records,
            unavailable=(
                f"{len(missing)} of {len(records)} host(s) could not be measured "
                f"({', '.join(missing)}), so the worst over them is unknown rather "
                f"than the worst of the rest"
            ),
        )

    # Each reduction needs EVERY host to have answered for that quantity, not
    # just some of them. The density reduction had no such guard: it took the
    # minimum over the hosts carrying a density, which is the "worst of the
    # rest" this whole module exists to refuse. Unreachable while
    # `selectivity_density_from_loads` returns 0.0 and never None for a measured
    # host -- but reachable the moment the per-host load is withheld, which
    # `per_host_load_is_in_force` now does on the default path, so the guard is
    # load-bearing rather than defensive.
    densities = [r for r in records if r.selectivity_density is not None]
    complete = len(densities) == len(records)
    worst_density = min(densities, key=lambda r: r.selectivity_density) if complete else None
    coverages = [r.coverage for r in records if r.coverage is not None]
    sites = [r.sites for r in records if r.sites is not None]
    return BackgroundReading(
        mode=mode,
        per_host=records,
        worst_host_selectivity_density=(
            worst_density.selectivity_density if worst_density else None
        ),
        worst_host=worst_density.prefix if worst_density else None,
        worst_host_background_sites=(
            max(int(value) for value in sites) if len(sites) == len(records) else None
        ),
        worst_host_bg_coverage=max(coverages) if len(coverages) == len(records) else None,
    )


def _unmeasured(prefix: str, length: int, reason: str) -> HostFigures:
    return HostFigures(
        prefix=prefix,
        length=length,
        sites=None,
        sites_per_mb=None,
        effective_load=None,
        coverage=None,
        selectivity_density=None,
        unavailable=reason,
    )


def _per_mb(sites: int | None, length: int) -> float | None:
    if sites is None or length <= 0:
        return None
    return float(sites) * 1e6 / float(length)


def _coverage_of(optimizer, positions_by_prefix, prefix, length) -> float | None:
    """This host's own coverage fraction, in its own coordinate space.

    `optimizer._compute_coverage` rather than a second implementation: the reach
    and the `fg_circular` geometry come off the run's config, and a host
    measured by different arithmetic from the target would not be comparable
    with it.
    """
    if length <= 0:
        return None
    return optimizer._compute_coverage(positions_by_prefix.get(prefix, ()) or (), length)


def _unscanned_hosts(optimizer, primers, prefixes) -> dict[str, str]:
    """Hosts this cache holds no ANSWER for, per host, with the reason.

    `has_entry` and not a hit count (Known Issue 15, and `require_entries`'
    rule): a primer indexed against a host it binds nowhere has an entry and a
    site count of zero, which is a measurement; a primer with no entry has no
    measurement, and the zero `get_positions` returns for it would read as
    perfect specificity.

    A cache that cannot be asked is itself an absent answer. Returning "all
    measured" there would let `worst-case` report a worst host computed from
    empty arrays.
    """
    cache = getattr(optimizer, "cache", None)
    if cache is None:
        return {prefix: "no position cache was attached" for prefix in prefixes}
    has_entry = getattr(cache, "has_entry", None)
    if not callable(has_entry):
        return {
            prefix: (
                "this position cache cannot say whether it holds an answer for "
                "a primer, so a per-host figure cannot be distinguished from a zero"
            )
            for prefix in prefixes
        }
    out: dict[str, str] = {}
    for prefix in prefixes:
        absent = [primer for primer in primers if not has_entry(prefix, primer)]
        if absent:
            out[prefix] = (
                f"{len(absent)} of {len(primers)} panel oligo(s) were never scanned "
                f"against this host ({absent[0]} is the first)"
            )
    return out


def _per_host_loads(
    optimizer, primers, prefixes, positions_by_prefix, selectivity_mode, aggregate
) -> tuple[dict[str, float | None], bool]:
    """Each host's binding load, on the same footing as the pooled one.

    Returns `(loads, skipped)`. `skipped` is True when the loads were not
    attempted at all because nothing will read them, which the caller turns into
    `load_unavailable` on a row whose other figures are real measurements. A
    `None` load with `skipped` False is a host whose load was attempted and
    could not be computed, and that makes the host unmeasured.

    Under mode "exact" the pooled load IS the raw site count -- that is what
    `_effective_site_load` falls back to when the mismatch-class tables are
    absent -- so each host's load is its own count, free and exactly consistent
    with the pooled figure. Nothing is skipped there because there is no cost.

    Otherwise `occupancy.weighted_site_load` is asked once per host, and only
    when `per_host_load_is_in_force`. The pooled call is NOT replaced by the sum
    of these: `(a + b) * theta` and `a * theta + b * theta` are not the same
    float, and the pooled figures must not move under `sum`. That is also why
    this cost is additive rather than free, and why it is not paid by a run that
    reads none of it.
    """
    if selectivity_mode == "exact":
        return {
            prefix: float(len(positions_by_prefix.get(prefix, ()) or ())) for prefix in prefixes
        }, False

    if not per_host_load_is_in_force(aggregate):
        return {prefix: None for prefix in prefixes}, True

    from neoswga.core.occupancy import weighted_site_load

    conditions = getattr(optimizer, "conditions", None)
    if conditions is None:
        return {prefix: None for prefix in prefixes}, False
    max_mismatches = int(getattr(optimizer.config, "max_mismatches", 1) or 0)
    model = getattr(optimizer, "mismatch_model", None)
    penalty = getattr(optimizer, "mismatch_penalty", None)

    loads: dict[str, float | None] = {}
    for prefix in prefixes:
        try:
            loads[prefix] = weighted_site_load(
                primers, [prefix], conditions, max_mismatches, penalty, model
            )
        except (FileNotFoundError, OSError):
            # The pooled figure fell back to exact counting for the same reason,
            # but it did so for the whole background at once. One host missing
            # its tables leaves THAT host unmeasured; substituting its exact
            # count here would put two definitions of one number in one table.
            loads[prefix] = None
    return loads, False
