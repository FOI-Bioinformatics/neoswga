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

The returned object carries only those five fields on purpose. A caller that
reaches for a sixth gets an `AttributeError` rather than a default that would be
indistinguishable from a measurement.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional


@dataclass(frozen=True)
class PoolMetrics:
    """What the pool objective and the plan row read, and nothing more."""

    fg_coverage: float
    effective_fg_coverage: Optional[float]
    selectivity_density: float
    total_bg_sites: int
    max_gap: float


def compute_pool_metrics(optimizer, primers: List[str]) -> PoolMetrics:
    """The five design quantities for one panel.

    Positions are kept PER PREFIX because each prefix is a separate sequence
    with its own 0-based coordinate space; pooling them makes position 5000 in
    two different targets the same integer. They are kept per primer as well
    because occupancy is a per-primer property and a pooled list cannot say
    which primer put a site where.
    """
    from .base_optimizer import _selectivity_density_from_loads

    if not primers:
        return PoolMetrics(0.0, None, 0.0, 0, float("inf"))

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

    fg_coverage = optimizer._coverage_over_prefixes(
        fg_by_prefix, optimizer.fg_prefixes, optimizer.fg_seq_lengths
    )

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

    gaps = []
    for prefix, length in zip(optimizer.fg_prefixes, optimizer.fg_seq_lengths, strict=True):
        gaps.extend(optimizer._compute_gaps(fg_by_prefix.get(prefix, []), length))

    return PoolMetrics(
        fg_coverage=fg_coverage,
        effective_fg_coverage=effective_fg_coverage,
        selectivity_density=selectivity_density,
        total_bg_sites=total_bg,
        max_gap=max(gaps) if gaps else float("inf"),
    )
