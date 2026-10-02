"""Coverage per foreground target, and the worst over them, from one loop.

Extracted from `base_optimizer` on 2026-10-02 rather than added to it: that
module is at its size budget, and the answer to a budget is extraction. Same
shape as `occupancy_coverage` and `selectivity` -- `TargetCoverage` and
`coverage_by_target` are re-exported from `base_optimizer`, and
`BaseOptimizer._coverage_over_prefixes` delegates here, so there is one loop and
one definition of what a target's coverage is.

Why it exists. `per_target_coverage` was filled in after the search by
`unified_optimizer`, from `coverage.compute_per_prefix_coverage`. Search,
acceptance and the bounded repair all score a panel through
`BaseOptimizer.compute_metrics`, which did not populate the field, so a floor on
the per-target figure could be reported and never acted on. That was recorded as
a deliberate decision whose stated reason was exactly this gap. Phase 6c of the
2026-10-02 diversity plan closes the gap; the loss it acts on is measured in
`docs/validation/2026-10-02-diversity-baseline.md` (Q2: under a pooled
five-strain foreground the worst strain stays 0.24 to 0.31 below the reference;
Q3: a held-out strain loses 0.004 to 0.32, and supergroup membership does not
predict which).

Three conventions, each of which has cost this project a defect.

**The aggregate is unchanged.** `aggregate` is the same length-weighted sum of
covered bases over total bases `_coverage_over_prefixes` has always returned,
computed by the same arithmetic in the same order, including the cases where it
returned 0.0. `fg_coverage` is not redefined here; a figure is added beside it.

**An unmeasured target is absent, not zero.** A prefix with no usable length,
or one the caller supplied no positions entry for, is named in `unmeasured` and
left out of `per_target`. `compute_per_prefix_coverage` writes 0.0 for a
zero-length prefix, which is indistinguishable from a target a panel reaches
nowhere in.

**A worst case over the members that answered is not a worst case.** `worst` is
`None` whenever any target is unmeasured, and `None` whenever there are no
targets at all. `pool_objective._unmet` treats `None` as missing its limit, so a
limit resting on this figure is refused rather than satisfied. Same rule, for the
same reason, as `reference_panel_evaluation._worst_coverage`.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

#: One prefix's coverage fraction from its positions, its length and its
#: reverse-strand positions. `BaseOptimizer` supplies the reach and the geometry
#: from its own config; this module supplies the loop and the bookkeeping.
FractionFor = Callable[[Sequence[int], int, Sequence[int]], float]


@dataclass(frozen=True)
class TargetCoverage:
    """Coverage over several targets: the aggregate, each one, and the worst.

    `per_target` holds only the targets that could be measured. `unmeasured`
    names the rest, so a caller can say which ones rather than only that some
    figure is missing.
    """

    aggregate: float
    per_target: dict[str, float]
    unmeasured: tuple[str, ...] = ()

    @property
    def worst(self) -> float | None:
        """The lowest per-target coverage, or `None` when that is unknown.

        `None` with any target unmeasured, because the minimum over the subset
        that answered reads as a worst case and is not one. `None` with no
        targets, because a minimum over nothing is not 0.0 either.
        """
        if self.unmeasured or not self.per_target:
            return None
        return min(self.per_target.values())

    @property
    def worst_target(self) -> str | None:
        """Which target holds `worst`, or `None` when `worst` is `None`."""
        if self.worst is None:
            return None
        return min(self.per_target, key=lambda prefix: (self.per_target[prefix], prefix))


def coverage_by_target(
    fraction_for: FractionFor,
    positions_by_prefix: Mapping[str, Sequence[int]],
    prefixes: Sequence[str],
    seq_lengths: Sequence[int],
    *,
    reverse_by_prefix: Mapping[str, Sequence[int]] | None = None,
) -> TargetCoverage:
    """Each sequence's coverage fraction, and the length-weighted aggregate.

    Each sequence is measured in its own coordinate space and the covered BASES
    are summed, which is the only aggregate that means anything when the parts
    are separate molecules. Equivalent to one call of `fraction_for` when there
    is one prefix, so single-target behaviour is unchanged.
    """
    total = sum(int(length) for length in seq_lengths)
    reverse = reverse_by_prefix or {}

    per_target: dict[str, float] = {}
    unmeasured: list[str] = []
    covered = 0.0
    for prefix, raw_length in zip(prefixes, seq_lengths, strict=True):
        length = int(raw_length)
        if length <= 0 or prefix not in positions_by_prefix:
            # No length is no denominator, and no entry is no answer. Neither is
            # a coverage of zero, and both used to be recorded as one.
            unmeasured.append(prefix)
            continue
        fraction = fraction_for(
            positions_by_prefix[prefix],
            length,
            reverse.get(prefix, ()),
        )
        per_target[prefix] = fraction
        covered += fraction * length

    if total <= 0:
        # The aggregate `_coverage_over_prefixes` has always returned for this
        # case. The per-target record is still carried, so a caller can see that
        # every target was unmeasured rather than that the panel covers nothing.
        return TargetCoverage(0.0, per_target, tuple(unmeasured))
    return TargetCoverage(covered / total, per_target, tuple(unmeasured))
