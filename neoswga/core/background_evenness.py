"""Gap evenness per strand, and the background average swga 2.0 fits.

swga 2.0's set-level model has five terms, fitted by ridge regression with
10-fold cross-validation against the proportion of the target genome reaching
1x sequencing coverage, over 46 published sets from *M. tuberculosis* and
*H. sapiens*. `off_gap_gini` carries the second-largest recorded weight,
+0.281 behind `freq_ratio` at +0.321: a Gini index of the distances between
background binding sites, computed per strand, averaged across strands and
then across genomes. A POSITIVE weight means an uneven background is better,
which is what the amplicon geometry predicts, since clustered host sites leave
most of the host unreachable.
`docs/validation/pool_selection_audit_2026-09-18.md` records that the term
list, the definitions, the response variable, the sample size and the method
are confirmed from the paper and **the coefficient values are not**: they sit
in a CAPTCHA-gated table, and two attempts disagreed about the sign of a
different term. Treat the numbers as recorded, not verified.

This module computes the quantity and nothing else. It is REPORTED, never
scored:

- The sample behind the weight is 46 sets with one response variable, read
  from a table nobody here could open. That does not support importing a
  coefficient, and `docs/validation/tool_comparison.md` already states that
  position.
- COATswga (2025) computes no Gini and no gap statistic at all, on the stated
  ground that a per-primer Gini cannot speak for a whole set. So the field
  does not agree that gap evenness belongs in an objective, which is Known
  Issue 18's standing reason for keeping `max_gap` out of every score.

**Unknown is not zero, and here the zero is the unfavourable end.** A strand
with fewer than two sites has no distance to measure, and a Gini of 0.0 means
perfectly even spacing -- the value swga 2.0's positive weight calls WORSE.
So an unmeasurable strand must not read as a measured 0.0 in either direction:
it is left out of the average, and a genome where neither strand can be
measured is absent from the result rather than present at zero.

**Inter-site distances only.** `PrimerSetMetrics.gap_gini` pools both strands
and adds the leading and trailing segments of the sequence (and a wrap under
`fg_circular`). This one takes only the distances between consecutive
same-strand sites, which is the quantity the fitted model used, and the two are
therefore not comparable.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

#: The per-genome figure `PositionCache.compute_strand_alternation_stats`
#: returns and `strand_metrics` collects.
GINI_KEY = "gap_gini_per_strand"


def gini(values: Sequence[float]) -> float | None:
    """Gini index of a list of distances, or `None` when there are none.

    `0.0` for a single distance and for any set of equal distances, which is
    "perfectly even" and a measurement. An empty list is not a measurement.
    """
    if not values:
        return None
    ordered = sorted(float(v) for v in values)
    n = len(ordered)
    total = sum(ordered)
    if total == 0:
        return 0.0
    weighted = sum((2 * (index + 1) - n - 1) * value for index, value in enumerate(ordered))
    return weighted / (n * total)


def _distances(positions: Sequence[int]) -> list[float]:
    """Distances between consecutive sites, de-duplicated and sorted first.

    `strict=False` DELIBERATELY: pairing each site with its successor makes
    `ordered[1:]` one shorter by construction, so `strict=True` would raise on
    every call. Same pairing as `base_optimizer._compute_gaps`.
    """
    ordered = sorted(set(int(p) for p in positions))
    return [float(right - left) for left, right in zip(ordered, ordered[1:], strict=False)]


def per_strand_gap_gini(by_strand: Mapping[str, Sequence[int]]) -> float | None:
    """Mean over the strands that have a distance at all, or `None`.

    Averaging over the MEASURABLE strands rather than over both is the whole
    point: a genome bound on one strand only would otherwise be averaged with
    a fabricated zero, and zero is the end of this scale that reads as bad.
    """
    measured = [
        value
        for positions in by_strand.values()
        if (value := gini(_distances(positions))) is not None
    ]
    if not measured:
        return None
    return sum(measured) / len(measured)


def off_gap_gini(
    strand_stats: Mapping[str, Mapping[str, float | None]],
    bg_prefixes: Sequence[str],
) -> float | None:
    """swga 2.0's `off_gap_gini`: the per-strand figure, averaged over hosts.

    `None` when there is no background genome, when none of them was measured,
    or when every measured one had no same-strand distance to offer. A design
    with no host has no host evenness, which is absence and not zero.
    """
    values = [
        value
        for prefix in bg_prefixes or ()
        if (value := (strand_stats.get(prefix) or {}).get(GINI_KEY)) is not None
    ]
    if not values:
        return None
    return sum(float(v) for v in values) / len(values)
