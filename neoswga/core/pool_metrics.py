"""The five quantities a pool design is searched and accepted on.

Task 6 of the condition-aware pool design plan, after profiling for step 204.

`BaseOptimizer.compute_metrics` answers every question anything has ever asked
of a primer set: both coverage figures, a coverage sweep across reaches,
background coverage, gap mean, maximum, Gini and Shannon entropy, melting
temperatures, dimer risk and strand alternation. That is the right shape for a
report, which is computed once.

It is the wrong shape for a search. `PoolObjective` and `plan_pool` between them
read exactly five of those fields, and the repair path calls the evaluator
thousands of times per design. Profiled on a 2 Mb synthetic target with 300
candidates and a density floor tight enough to force repairs, `plan_pool` spent
29.7 s, of which 27.6 s was inside `compute_metrics` across 4,040 calls, and a
large part of each call produced quantities nothing downstream read.

This computes the five that are read, through the same helpers, and nothing
else. `tests/test_pool_metrics_agree_with_the_full_evaluation.py` pins it
against `compute_metrics` on randomised panels, which is what stops the two
definitions drifting while each stays self-consistent.

The returned object carries only the fields something reads on purpose. A
caller reaching for one nothing computes gets an `AttributeError` rather than a
default that would be indistinguishable from a measurement.
`worst_target_coverage` joined it on 2026-10-02 for that reason and not as a
convenience: `min_per_target_coverage` became an enforced limit, so
`PoolObjective.violations` now reads that field, and the objective's `getattr`
default would have made the floor unmeetable on every panel.

**`bg_coverage` joined them on 2026-11-02, and the cost objection did not apply
to it.** This docstring used to say a caller reaching for a further field would
take that `AttributeError`, and `stage2_background_weight` was refused on those
grounds. Measured, the field costs +48% per evaluation (3.56 ms to 5.26 ms on a
4 Mb host) and is **2.05x cheaper** than the full `compute_metrics` that
`panel_refinement` already falls back to whenever a host-coverage limit is set
-- so for a host-limited design it is a speedup, not a tax. The reason it is
cheap is that this function ALREADY gathers and sorts every background position
for `total_bg_sites` and then discarded them: the marginal work is one
`_coverage_over_prefixes` call over a list in hand, not the position lookups the
27.6 s profile was about.
`docs/validation/2026-11-02-what-host-coverage-would-cost-the-search.md`.

It stays `None` with no background genome, where `compute_metrics` reports 0.0.
That divergence is deliberate: a host coverage measured against no host is not
zero, and a weighted host axis must refuse it rather than trade against a
fabricated zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class PoolMetrics:
    """What the pool objective and the plan row read, and nothing more."""

    fg_coverage: float
    effective_fg_coverage: float | None
    selectivity_density: float
    total_bg_sites: int
    max_gap: float
    #: Fraction of the background the panel reaches. `None` with no background
    #: genome -- not 0.0, which is what a weighted host axis would trade
    #: against. Defaulted so the five-field constructor calls still read.
    bg_coverage: float | None = None

    # Added 2026-10-02, when `min_per_target_coverage` became an enforced
    # limit and `PoolObjective.violations` began reading this field. It is here
    # because it is READ, which is the rule this module keeps: without it the
    # objective's `getattr(metrics, metric, None)` would hand `_unmet` a None on
    # every panel, so the floor would be unmeetable while reading as a violated
    # constraint and the bounded repair would chase it to no end.
    #
    # The MINIMUM over the per-target fractions, and `None` when any target
    # could not be measured -- never the lowest of the ones that answered. The
    # loop is `core/target_coverage.py`, the same one `compute_metrics` uses.
    worst_target_coverage: float | None = None

    # The per-host reading of the two background figures above, added
    # 2026-10-02. Not a sixth design quantity: it is the SAME two limits read
    # under `background_aggregate: worst-case`, and the search has to see what
    # the acceptance path will hold the panel to or a repair would chase the
    # pooled figure and be rejected on the worst host's. `per_host_metrics`
    # produces all of it in one call, and `background_aggregate` is here so
    # `limit_figure` can dispatch without reaching for a global.
    #
    # `worst_host_bg_coverage` is deliberately absent: `max_host_coverage`
    # routes to the full evaluator (`panel_refinement.objective_for_optimizer`),
    # and a limit this object cannot answer for reads `None` and is refused.
    # The pooled `bg_coverage` above joined on main after this was written; the
    # per-host reading does not use it yet.
    background_aggregate: str = "sum"
    per_host: dict[str, dict] = field(default_factory=dict)
    worst_host: str | None = None
    worst_host_selectivity_density: float | None = None
    worst_host_background_sites: int | None = None


def compute_pool_metrics(optimizer, primers: list[str]) -> PoolMetrics:
    """The five design quantities for one panel.

    Positions are kept PER PREFIX because each prefix is a separate sequence
    with its own 0-based coordinate space; pooling them makes position 5000 in
    two different targets the same integer. They are kept per primer as well
    because occupancy is a per-primer property and a pooled list cannot say
    which primer put a site where.
    """
    from .base_optimizer import _selectivity_density_from_loads
    from .per_host_metrics import measure_background

    if not primers:
        # An empty panel covers 0.0 of a background that EXISTS, which is a
        # measurement like the 0.0 foreground coverage beside it. Absence is
        # reserved for having no background to measure against at all.
        # `worst_target_coverage` stays None, as in `PrimerSetMetrics.empty()`:
        # no panel is not a measured zero on any target.
        empty_bg = (
            0.0 if optimizer.bg_prefixes and int(sum(optimizer.bg_seq_lengths or [])) > 0 else None
        )
        return PoolMetrics(0.0, None, 0.0, 0, float("inf"), empty_bg, None)

    fg_by_prefix = {prefix: set() for prefix in optimizer.fg_prefixes}
    bg_by_prefix = {prefix: set() for prefix in optimizer.bg_prefixes}
    fg_positions_by_primer = {prefix: {} for prefix in optimizer.fg_prefixes}

    for primer in primers:
        for prefix in optimizer.fg_prefixes:
            positions = optimizer.get_primer_positions(primer, prefix, "both")
            fg_by_prefix[prefix].update(positions.tolist())
            fg_positions_by_primer[prefix][primer] = positions.tolist()
        for prefix in optimizer.bg_prefixes:
            positions = optimizer.get_primer_positions(primer, prefix, "both")
            bg_by_prefix[prefix].update(positions.tolist())

    fg_by_prefix = {p: sorted(v) for p, v in fg_by_prefix.items()}
    bg_by_prefix = {p: sorted(v) for p, v in bg_by_prefix.items()}

    # Both through the same helpers `compute_metrics` uses, so the two
    # evaluators cannot disagree about either figure -- including about whether
    # the per-target one was computed at all, which depends on a configured
    # floor. Two different geometries on purpose: the aggregate unconfined, the
    # per-target ones confined to each record (audit F3). See
    # `core/target_coverage.py`.
    fg_coverage = optimizer._coverage_over_prefixes(
        fg_by_prefix, optimizer.fg_prefixes, optimizer.fg_seq_lengths
    )
    target_coverage = optimizer._confined_coverage_by_target(fg_by_prefix)

    effective_fg_coverage = None
    if optimizer.conditions is not None:
        total = optimizer.fg_total_length
        effective_fg_coverage = (
            sum(
                optimizer._compute_effective_coverage(fg_positions_by_primer[prefix], length)
                * length
                for prefix, length in zip(
                    optimizer.fg_prefixes, optimizer.fg_seq_lengths, strict=True
                )
                if length > 0
            )
            / total
            if total > 0
            else 0.0
        )

    # Flat lists for the site COUNTS, where summing across sequences is exactly
    # right and de-duplicating across them would be wrong.
    total_fg = sum(len(fg_by_prefix[p]) for p in optimizer.fg_prefixes)
    total_bg = sum(len(bg_by_prefix[p]) for p in optimizer.bg_prefixes)
    effective_fg, effective_bg, mode = optimizer._effective_site_load(primers)
    if mode == "exact":
        effective_fg, effective_bg = float(total_fg), float(total_bg)
    selectivity_density = _selectivity_density_from_loads(
        effective_fg,
        int(sum(optimizer.fg_seq_lengths or [])),
        effective_bg,
        int(sum(optimizer.bg_seq_lengths or [])),
    )

    # One call over `bg_by_prefix`, which is already gathered and sorted above
    # for `total_bg_sites`. `None` rather than 0.0 with no background: see the
    # module docstring.
    bg_coverage = None
    if optimizer.bg_prefixes and int(sum(optimizer.bg_seq_lengths or [])) > 0:
        bg_coverage = optimizer._coverage_over_prefixes(
            bg_by_prefix, optimizer.bg_prefixes, optimizer.bg_seq_lengths
        )

    gaps = []
    for prefix, length in zip(optimizer.fg_prefixes, optimizer.fg_seq_lengths, strict=True):
        gaps.extend(optimizer._compute_gaps(fg_by_prefix.get(prefix, []), length))

    # The same per-host reading the full evaluator records, from the per-prefix
    # positions already gathered above. The pooled coverage is passed as None,
    # so no worst-host coverage is reduced here although `bg_coverage` is now
    # computed above; `max_host_coverage` routes to `compute_metrics`, which
    # has both.
    background = measure_background(
        optimizer,
        primers,
        bg_by_prefix,
        pooled_sites=total_bg,
        pooled_coverage=None,
        pooled_density=selectivity_density,
        effective_fg=effective_fg,
        effective_bg=effective_bg,
        selectivity_mode=mode,
    )
    reading = background.kwargs()
    reading.pop("worst_host_bg_coverage")

    return PoolMetrics(
        fg_coverage=fg_coverage,
        effective_fg_coverage=effective_fg_coverage,
        selectivity_density=selectivity_density,
        total_bg_sites=total_bg,
        max_gap=max(gaps) if gaps else float("inf"),
        bg_coverage=bg_coverage,
        worst_target_coverage=target_coverage.worst,
        **reading,
    )
