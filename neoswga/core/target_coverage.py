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

Four conventions, each of which has cost this project a defect.

**The per-target figure is RECORD-CONFINED, and the aggregate is not.** This
reverses the first version of Phase 6c, which took the per-target dict from the
loop that produces `fg_coverage` because the plan said to. That loop is
`coverage.union_coverage`, which does not confine a window to the record holding
its site, and audit finding F3 of
`docs/validation/additives_coverage_audit_2026-09-14/` records exactly that as a
P1 error for target coverage on a fragmented or multi-chromosome reference, fixed
and pinned by `tests/test_coverage_windows_stop_at_record_edges.py`. On a
500-contig 2 Mb draft with a site 100 bp inside each contig end the two figures
read 0.775 confined and 0.9996 unconfined, so `min_per_target_coverage: 0.9`
read against the unconfined figure passes a panel that does not meet it. A floor
must be compared against the right quantity, so `confined_coverage_by_target`
uses `coverage.confined_covered_bases` -- the same function
`compute_per_prefix_coverage` uses, so `optimize` and `evaluate-set` report one
number under one key name.

**`fg_coverage` is untouched.** `_coverage_over_prefixes` still delegates here
with the unconfined fraction, by the same arithmetic in the same order,
including the cases where it returned 0.0. It is therefore NOT the
length-weighted mean of `per_target_coverage` on a multi-record reference, and
that is deliberate: redefining the aggregate would move every recorded figure
and every selection stage at once, which is a decision with its own
measurement, not a repair. The two differ only where a record join lies within
one reach of a site.

**An unmeasured target is absent, not zero.** A prefix with no usable length,
or one the caller supplied no positions entry for, is named in `unmeasured` and
left out of `per_target`. `compute_per_prefix_coverage` writes 0.0 for a
zero-length prefix, which is indistinguishable from a target a panel reaches
nowhere in. Through `compute_metrics` only the non-positive length can reach
this, because `_gather_sites` seeds every prefix; see `require_per_target`.

**A worst case over the members that answered is not a worst case.** `worst` is
`None` whenever any target is unmeasured, and `None` whenever there are no
targets at all. `pool_objective._unmet` treats `None` as missing its limit, so a
limit resting on this figure is refused rather than satisfied. Same rule, same
reason AND now the same quantity as `reference_panel_evaluation._worst_coverage`.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass

from .exceptions import ModelEvaluationError

#: One prefix's coverage fraction from its name, its positions, its length and
#: its reverse-strand positions. The prefix is passed because the confined
#: geometry needs that prefix's record starts. `BaseOptimizer` supplies the
#: reach and the geometry from its own config; this module supplies the loop and
#: the bookkeeping.
FractionFor = Callable[[str, Sequence[int], int, Sequence[int]], float]


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
            prefix,
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


def confined_coverage_by_target(
    cache,
    positions_by_prefix: Mapping[str, Sequence[int]],
    prefixes: Sequence[str],
    seq_lengths: Sequence[int],
    *,
    extension: int,
    circular: bool,
) -> TargetCoverage:
    """Each target's coverage, with every window confined to its own record.

    The figure `min_per_target_coverage` is read against, and the one
    `evaluate-set` reports: both go through `coverage.confined_covered_bases`.
    It takes the positions the caller has ALREADY gathered, so it costs one
    bitmap and one marking pass per target and no extra position lookups --
    `cache` is consulted only for the record starts, through
    `coverage._record_starts_for`, which raises rather than guessing when a
    cache that has the getter fails.

    `aggregate` here is the length-weighted mean of the confined fractions. It
    is NOT `fg_coverage`, which stays unconfined; nothing reads it today and it
    is returned only because `TargetCoverage` carries it.
    """
    from .coverage import _record_starts_for, confined_covered_bases

    def fraction_for(prefix, positions, length, _reverse):
        bases = confined_covered_bases(
            positions,
            length,
            extension,
            circular,
            record_starts=_record_starts_for(cache, prefix),
        )
        return bases / length

    return coverage_by_target(fraction_for, positions_by_prefix, prefixes, seq_lengths)


def require_per_target(per_target, worst, fg_prefixes, floor) -> None:
    """Refuse a requested per-target floor that nothing measured.

    `panel_acceptance.check_per_target_coverage` answers "nothing to say" on an
    empty dict and `result_validation` gates its warning on a non-empty one, so
    an empty dict made a REQUESTED floor pass vacuously: the panel that starves
    one target was reported as satisfying the limit set to stop exactly that,
    and aggregate coverage cannot reveal it either, which is why the floor
    exists. `docs/development/design-contracts.md` records that defect.

    The same applies one target at a time: `worst` is `None` when any target
    went unmeasured, and the minimum over the targets that answered is not a
    worst case. Both are refused here rather than reported as met.

    Nothing is refused when no floor was asked for. A single-genome run and a
    run with an unusable length are then unchanged, which is what keeps the
    per-target work from moving any delivered panel.

    **What can actually reach these two raises, through the pipeline.** Only a
    non-positive `fg_seq_lengths` entry, and an empty dict only when every entry
    is non-positive. `base_optimizer._gather_sites` seeds `out["fg"]` with every
    foreground prefix, so the "no positions entry" branch of
    `coverage_by_target` cannot fire from `compute_metrics`; and
    `base_optimizer.get_primer_positions` catches bare `Exception` and returns
    an empty array with a warning, so an unindexed prefix or an HDF5 failure
    arrives here as 0.0 coverage -- measured-looking, not unmeasured. That fails
    SAFE for a floor, since 0.0 misses any positive one, and it is why these
    raises are rarely the thing that fires. It is also a product defect of its
    own: `compute_per_prefix_coverage` catches only `(KeyError, ValueError)` and
    says "Do NOT swallow I/O / HDF5 errors here". The broad swallow predates
    this limit and is not narrowed here, because narrowing it changes what every
    existing design does with a damaged index.
    """
    if not floor:
        return
    if not per_target:
        raise ModelEvaluationError(
            "per_target_coverage",
            sorted(fg_prefixes),
            f"a per-target floor of {floor} was requested and no target's "
            "coverage could be measured, so the floor would pass vacuously",
        )
    if worst is None:
        unmeasured = sorted(set(map(str, fg_prefixes)) - set(per_target))
        raise ModelEvaluationError(
            "worst_target_coverage",
            unmeasured or sorted(fg_prefixes),
            f"a per-target floor of {floor} was requested and coverage is "
            f"unavailable for {len(unmeasured)} of {len(fg_prefixes)} targets, "
            "so the worst over them is unknown rather than the lowest of the rest",
        )


def saturation_issues(per_target, fg_prefixes, fg_seq_lengths, n_primers, extension):
    """Targets small enough that any panel covers them, as warning rows.

    When `n_primers * 2 * extension >= genome_len` a set covers the genome
    regardless of design quality, so the coverage figure is saturation-bounded
    rather than a property of the panel. Said out loud so a plasmid-scale run is
    not read as a perfect design. Warnings only: saturation is not a defect.

    A target absent from `per_target` is skipped rather than read as 0.0, which
    would be the one value that can never trip this warning and so would hide
    an unmeasured target instead of leaving it to `require_per_target`.
    """
    issues = []
    if n_primers <= 0 or extension <= 0:
        return issues
    for prefix, length in zip(fg_prefixes, fg_seq_lengths, strict=True):
        if length <= 0 or prefix not in per_target:
            continue
        window = n_primers * 2 * int(extension)
        coverage = per_target[prefix]
        if length / window < 1.0 and coverage >= 0.95:
            issues.append(
                {
                    "level": "warning",
                    "code": "coverage_saturated_on_small_genome",
                    "detail": (
                        f"{prefix}: genome={length} bp, {n_primers} primers "
                        f"x 2x {int(extension)} bp reach = {window} bp; "
                        f"coverage={coverage:.1%} is saturation-bounded "
                        "(metric unreliable)"
                    ),
                }
            )
    return issues
