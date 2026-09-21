"""How much of the missing depth a panel would recover.

Finding F7 of the 2026-09-16 pipeline audit. BAM-guided expansion did not
target the observed gaps, and three defects shared one cause: the objective
asked whether a binding SITE fell inside a gap rather than how much of the
missing depth a candidate's WINDOW would recover.

**The reach was never consulted.** A candidate binding just outside a 50 kb gap
whose modelled window would blanket it was discarded, while one binding at the
gap's last base and extending away was kept.

**The filter was all-or-nothing.** When fewer than `target_new` candidates
survived it, the whole filter was abandoned for the unfiltered pool, so the run
silently stopped targeting gaps.

**`gap_coverage` counted gaps rather than bases.** `1 - len(after)/len(before)`
goes NEGATIVE when one long gap splits into two, reporting progress as
regression.

One weight per base replaces all three:

    weight[r] = max(0, 1 - observed_depth[r] / desired_depth)

and a panel's worth is the weight lying under the union of its binding windows.
That is a number of bases, so it cannot go negative, it cannot double-count an
overlap, and it ranks a candidate by what it recovers rather than by where it
happens to bind.

**A base whose depth was never evaluable carries no deficit.** That is the link
to `reference_layout`: a record the BAM does not cover has UNKNOWN depth, and
treating unknown as zero would design primers for a region nothing measured.
Absent evidence is not a deficit, and this is the same distinction Known Issues
5, 6 and 13 record.

Windows are marked through `coverage._mark_window`, so they are confined to the
record holding the site when record starts are supplied. A primer near the end
of one molecule cannot recover missing depth on another.

`DeficitObjective` composes `PoolObjective` so the swap refinement can drive
selection by this quantity with no change to the refinement itself: only
`coverage` is redefined, and the constraints a design is accepted on survive
untouched. Measuring a deficit without selecting on it changes nothing useful,
which is the lesson Known Issue 14 records.
"""

from __future__ import annotations

import logging
from typing import List, Optional, Sequence, Tuple

import numpy as np

from neoswga.core.coverage import _mark_window

logger = logging.getLogger(__name__)


def deficit_weights(
    depth: Sequence[float],
    desired_depth: float,
    evaluable: Optional[Sequence[bool]] = None,
) -> np.ndarray:
    """Per-base deficit in [0, 1]: 1 where nothing was seen, 0 where enough was.

    Clamped at zero, so a base with more than the desired depth is not a
    negative deficit that another base's shortfall could cancel.

    `evaluable` zeroes every base a BAM could not speak for. Without it, an
    unobserved record reads as a total deficit and attracts the whole
    expansion budget.
    """
    if desired_depth <= 0:
        raise ValueError(f"desired_depth must be positive, not {desired_depth!r}")

    observed = np.asarray(depth, dtype=np.float64)
    weights = np.clip(1.0 - observed / float(desired_depth), 0.0, 1.0)
    if evaluable is not None:
        weights = np.where(np.asarray(evaluable, dtype=bool), weights, 0.0)
    return weights


def covered_mask(
    positions: Sequence[int],
    extension: int,
    length: int,
    circular: bool,
    record_starts: Optional[Sequence[int]] = None,
) -> np.ndarray:
    """The union of these binding windows, as a boolean array.

    A union rather than a sum, which is what makes a recovery figure a count of
    bases: two windows over the same deficit recover it once.
    """
    occupied = np.zeros(int(length), dtype=bool)
    for position in positions:
        _mark_window(
            occupied,
            int(position),
            int(extension),
            int(length),
            bool(circular),
            record_starts=record_starts,
        )
    return occupied


def recovered_deficit(
    positions: Sequence[int],
    weights: Sequence[float],
    extension: int,
    length: int,
    circular: bool,
    record_starts: Optional[Sequence[int]] = None,
) -> float:
    """Deficit bases lying under the union of these windows."""
    occupied = covered_mask(positions, extension, length, circular, record_starts)
    return float(np.asarray(weights, dtype=np.float64)[occupied].sum())


def deficit_gain(
    positions: Sequence[int],
    already_covered: np.ndarray,
    weights: Sequence[float],
    extension: int,
    length: int,
    circular: bool,
    record_starts: Optional[Sequence[int]] = None,
) -> float:
    """What these windows add over what the pool already reaches.

    Marginal by construction, so it is never negative: a candidate whose window
    falls entirely inside the panel's reach gains nothing rather than
    subtracting.
    """
    occupied = covered_mask(positions, extension, length, circular, record_starts)
    fresh = occupied & ~np.asarray(already_covered, dtype=bool)
    return float(np.asarray(weights, dtype=np.float64)[fresh].sum())


def dilate_intervals(
    intervals: Sequence[Tuple[int, int]],
    reach: int,
    length: int,
) -> List[Tuple[int, int]]:
    """Gap intervals widened by one reach, merged, clamped to the sequence.

    The prescreen the gap filter becomes. A candidate binding within one reach
    of a gap can blanket part of it, so membership was the wrong test: it threw
    away exactly the candidates whose windows would have helped.

    This is a BUDGET prescreen and nothing more. It narrows what is scored; it
    does not decide, and it has no fallback, because abandoning it when it
    leaves too few candidates is what made the old filter silently stop
    targeting gaps at all.
    """
    if not intervals:
        return []

    widened = sorted(
        (max(0, int(start) - int(reach)), min(int(length), int(end) + int(reach)))
        for start, end in intervals
    )

    merged: List[Tuple[int, int]] = []
    for start, end in widened:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


class DeficitObjective:
    """A `PoolObjective` whose coverage is recovered deficit, not genome breadth.

    Composed rather than substituted. `metrics`, `violations` and `shortfall`
    are the wrapped objective's, so the specificity floor, the constraint
    ordering and everything else a design is accepted on survive untouched;
    only `coverage` is redefined. That is what lets `refine_by_swaps` drive it
    with no change at all, since its score is
    `(-shortfall, coverage, -metrics(panel).total_bg_sites)`.

    Replacing the objective instead would silently drop every constraint the
    expansion was configured with, and redefining `shortfall` would let a
    coverage figure outbid a constraint -- the ordering `shortfall` exists to
    protect.

    Measuring the deficit is not targeting it. Known Issue 14 records the
    lesson: a stage that reads a quantity without choosing the panel changes
    nothing useful, so this exists to be attached to the search rather than
    consulted afterwards.
    """

    def __init__(
        self,
        inner,
        cache,
        weights_by_prefix,
        lengths_by_prefix,
        extension: int,
        circular: bool = False,
    ):
        self._inner = inner
        self._cache = cache
        self._weights = dict(weights_by_prefix)
        self._lengths = dict(lengths_by_prefix)
        self._extension = int(extension)
        self._circular = bool(circular)
        self._total = sum(float(np.asarray(w).sum()) for w in self._weights.values())

    metric_name = "predicted_deficit_recovery"

    # -- delegated, deliberately ------------------------------------------
    def metrics(self, primers):
        return self._inner.metrics(primers)

    def violations(self, primers):
        return self._inner.violations(primers)

    def shortfall(self, primers):
        return self._inner.shortfall(primers)

    @property
    def constraints(self):
        return getattr(self._inner, "constraints", None)

    # -- the one redefinition ---------------------------------------------
    def coverage(self, primers):
        """Fraction of the targeted deficit this panel's windows recover.

        Zero when there is no deficit anywhere, rather than a division by
        zero: a panel cannot recover what is not missing, and reporting 1.0
        would make every panel look perfect on a genome with full depth.
        """
        if self._total <= 0:
            return 0.0

        recovered = 0.0
        for prefix, weights in self._weights.items():
            length = int(self._lengths.get(prefix, len(weights)))
            positions = []
            for primer in primers:
                positions.extend(int(p) for p in self._cache.get_positions(prefix, primer, "both"))
            if not positions:
                continue
            recovered += recovered_deficit(
                positions,
                weights,
                extension=self._extension,
                length=length,
                circular=self._circular,
                record_starts=self._cache.get_record_starts(prefix) or None,
            )
        return float(min(1.0, recovered / self._total))
