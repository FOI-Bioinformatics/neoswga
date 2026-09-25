"""Collecting the strand quantities, all of them, for every genome.

`PositionCache.compute_strand_alternation_stats` returns five figures per
genome. The call site in `base_optimizer._compute_metrics` read two, dropped
three, and `break`ed after the first foreground prefix, so a pan-target design
reported one target's strand structure as the panel's and the host's was never
computed although the same method would have produced it.

`strand_alternation_gap_max` is the discarded figure that matters most.
Exponential amplification needs two sites in convergent orientation within the
polymerase's reach; one site primes linearly at best. So the widest gap between
opposite-strand sites is the closest quantity in this codebase to the mechanism
SWGA runs on, and on the HOST it is the term swga 2.0 approximates with
`within_mean_gap_ratio` and fits against measured sequencing breadth. Item 3 of
`docs/validation/getting_ahead_on_spacing_2026-09-18.md`.

Nothing here scores or constrains. The quantities are collected, carried on
`PrimerSetMetrics`, serialised, and reported by `panel_regime`; no spacing
threshold derived from the reach separates the published wet-lab winners, so
turning one into a default would be the scoring change that evidence refuses.

**A prefix the cache cannot answer for is absent from the result, not zero.**
The old call site initialised both scalars to 0.0 and left them there when the
cache could not answer, so a zero meant either "measured zero" or "never
asked". A one-site panel genuinely scores 0.0 for alternation, which is why the
two cannot be told apart by value alone -- the shape of Known Issues 5, 6
and 13.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

logger = logging.getLogger(__name__)

# The five figures `compute_strand_alternation_stats` returns, named here so a
# consumer can iterate them and a dropped one is visible.
STRAND_KEYS: tuple[str, ...] = (
    "strand_alternation_gap_mean",
    "strand_alternation_gap_max",
    "strand_alternation_score",
    "strand_coverage_ratio",
    "longest_same_strand_run",
)


def collect_strand_stats(
    cache: Any,
    prefixes: Sequence[str],
    seq_lengths: Sequence[int],
    primers: Sequence[str],
) -> dict[str, dict[str, float]]:
    """The five strand figures for every genome the cache can answer for.

    Args:
        cache: A `PositionCache`. Anything without
            `compute_strand_alternation_stats` yields `{}`, because
            `StreamingPositionCache` has no such method and metrics must stay
            computable rather than raising.
        prefixes: Foreground and background prefixes together. Each is measured
            against its own length, so the host's gaps are not scaled by the
            target's size.
        seq_lengths: Aligned with `prefixes`. A mismatched pair yields `{}`
            rather than zipping short, which would silently measure the wrong
            genomes.
        primers: The delivered panel.

    Returns:
        `prefix -> {key: value}`, containing only the prefixes that answered.
        A missing prefix means the quantity was not measured, which is not the
        same as measuring zero.
    """
    if not primers or not prefixes:
        return {}
    if len(prefixes) != len(seq_lengths):
        logger.debug(
            "Strand stats skipped: %d prefixes against %d lengths",
            len(prefixes),
            len(seq_lengths),
        )
        return {}
    if not hasattr(cache, "compute_strand_alternation_stats"):
        return {}

    collected: dict[str, dict[str, float]] = {}
    for prefix, length in zip(prefixes, seq_lengths, strict=True):
        try:
            stats = cache.compute_strand_alternation_stats(prefix, list(primers), length)
        except Exception as exc:
            # Absent rather than zero: a prefix nobody could measure must not
            # read as a measurement.
            logger.debug(f"Strand stats unavailable for {prefix}: {exc}")
            continue
        collected[prefix] = {key: stats[key] for key in STRAND_KEYS if key in stats}
    return collected


def headline_strand_scalars(
    stats: dict[str, dict[str, float]],
    fg_prefixes: Sequence[str],
) -> tuple[float | None, float | None]:
    """The two scalars `PrimerSetMetrics` has always carried.

    They describe the FIRST foreground genome, which is what the previous call
    site intended. `(None, None)` when no foreground genome was measured, so a
    host-only result cannot masquerade as the target's balance and an
    uncomputed value stays distinguishable from a measured zero.
    """
    for prefix in fg_prefixes:
        entry = stats.get(prefix)
        if entry:
            return (
                entry.get("strand_alternation_score"),
                entry.get("strand_coverage_ratio"),
            )
    return None, None


def worst_convergent_gap(
    stats: dict[str, dict[str, float]], prefixes: Sequence[str]
) -> float | None:
    """The widest gap between opposite-strand sites across these genomes.

    On the foreground this is the largest stretch a panel cannot amplify
    exponentially, because it holds no convergent pair. On the background it is
    the opposite reading: a large value means the host's sites are too far apart
    to face each other, which is the structure `bg_coverage` also sees and
    `total_bg_sites` cannot.
    """
    values: list[float] = [
        stats[prefix]["strand_alternation_gap_max"]
        for prefix in prefixes
        if prefix in stats and "strand_alternation_gap_max" in stats[prefix]
    ]
    return max(values) if values else None


def panel_occupancy(primers: Sequence[str], conditions: Any) -> dict[str, float]:
    """Fraction of the time each primer in the panel is bound, per primer.

    The ingredient item 6 of
    `docs/validation/getting_ahead_on_spacing_2026-09-18.md` asked for, without
    the recipe it proposed. An occupancy-WEIGHTED gap statistic is not offered:
    the 18 published sets with wet-lab outcomes carry published gap figures
    rather than binding positions, so no weighting rule can be validated
    against them, and the unweighted statistic already fails to separate their
    winners at every reach-derived threshold.

    What this does support is telling a user how weakly their weakest primer is
    bound. Measured, occupancy spans 7 to 9 fold within a panel at equiphi29
    42 C and only 1.7 to 3.0 fold at phi29 30 C, so it is a real quantity on
    the platform where additives work and nearly constant on the one where they
    do not.

    `{}` when there are no conditions, because without a temperature there is
    no occupancy to evaluate and a fabricated zero would be indistinguishable
    from a measurement.
    """
    if not primers or conditions is None:
        return {}
    from .occupancy import site_occupancy
    from .thermodynamics import calculate_enthalpy_entropy

    temp = getattr(conditions, "temp", None)
    if temp is None:
        return {}

    out: dict[str, float] = {}
    for primer in primers:
        sequence = str(primer).upper()
        if sequence in out:
            continue
        try:
            dh, _ds = calculate_enthalpy_entropy(sequence)
            tm = conditions.calculate_effective_tm(sequence)
        except Exception as exc:
            logger.debug(f"Occupancy unavailable for {sequence}: {exc}")
            continue
        out[sequence] = site_occupancy(dh, tm, temp)
    return out
