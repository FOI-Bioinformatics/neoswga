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

**Unknown is not the worst of the rest.** A host whose positions were never
scanned has no answer, and a reduction over a set with an unmeasured member is
unmeasured: the figure is `None`, never the minimum or maximum over the hosts
that happened to answer. `pool_objective._unmet` treats `None` as missing its
limit, so a background-measured limit under `worst-case` with an unscanned host
is REFUSED rather than satisfied. `position_cache.has_entry` is the question
asked, because an empty entry (scanned, binds nowhere) and an absent one (never
scanned) are different answers and only the first is a measurement.

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
    "limit_figure",
    "measure_background",
]

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

    `unavailable` is non-empty whenever the host could not be measured, and then
    every figure is `None`. A host that was scanned and carries no site has
    `sites == 0` and an empty `unavailable`: that zero is a result.
    """

    prefix: str
    length: int
    sites: int | None
    sites_per_mb: float | None
    effective_load: float | None
    coverage: float | None
    selectivity_density: float | None
    unavailable: str = ""

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
    Otherwise one `occupancy.weighted_site_load` call per host is made, and a
    host whose k-mer tables are missing is unmeasured rather than zero.
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

    if len(prefixes) == 1:
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

    unscanned = _unscanned_hosts(optimizer, primers, prefixes)
    loads = _per_host_loads(optimizer, primers, prefixes, positions_by_prefix, selectivity_mode)
    fg_length = int(sum(getattr(optimizer, "fg_seq_lengths", None) or []))

    records = []
    for prefix, length in zip(prefixes, lengths, strict=True):
        if prefix in unscanned:
            records.append(_unmeasured(prefix, length, unscanned[prefix]))
            continue
        sites = len(positions_by_prefix.get(prefix, ()) or ())
        load = loads.get(prefix)
        if load is None:
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
                selectivity_density=selectivity_density_from_loads(
                    effective_fg, fg_length, load, length
                ),
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

    densities = [r for r in records if r.selectivity_density is not None]
    worst_density = min(densities, key=lambda r: r.selectivity_density) if densities else None
    coverages = [r.coverage for r in records if r.coverage is not None]
    return BackgroundReading(
        mode=mode,
        per_host=records,
        worst_host_selectivity_density=(
            worst_density.selectivity_density if worst_density else None
        ),
        worst_host=worst_density.prefix if worst_density else None,
        worst_host_background_sites=max(int(r.sites) for r in records),
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
    optimizer, primers, prefixes, positions_by_prefix, selectivity_mode
) -> dict[str, float | None]:
    """Each host's binding load, on the same footing as the pooled one.

    Under mode "exact" the pooled load IS the raw site count -- that is what
    `_effective_site_load` falls back to when the mismatch-class tables are
    absent -- so each host's load is its own count, free and exactly consistent
    with the pooled figure.

    Otherwise `occupancy.weighted_site_load` is asked once per host. The pooled
    call is NOT replaced by the sum of these: `(a + b) * theta` and
    `a * theta + b * theta` are not the same float, and the pooled figures must
    not move under `sum`.
    """
    if selectivity_mode == "exact":
        return {
            prefix: float(len(positions_by_prefix.get(prefix, ()) or ())) for prefix in prefixes
        }

    from neoswga.core.occupancy import weighted_site_load

    conditions = getattr(optimizer, "conditions", None)
    if conditions is None:
        return {prefix: None for prefix in prefixes}
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
    return loads
