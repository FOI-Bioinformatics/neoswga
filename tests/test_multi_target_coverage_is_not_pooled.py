"""Coverage across several targets must be measured per target.

`compute_metrics` collected binding positions into ONE flat list across every
foreground prefix, deduplicated it, and divided by the SUM of the genome
lengths. Each prefix is a separate HDF5 with its own 0-based coordinate space,
so position 5000 in target A and position 5000 in target B are the same integer
and `set()` collapsed them, while no position could ever exceed the longest
single genome and the denominator was the total of all of them.

Two fully covered 100 kb targets therefore reported:

    fg_coverage   1.0000 on one target, 0.5145 on two
    max_gap       102,000 bp on a 100,000 bp target

A gap longer than the sequence it is measured in is the tell. `compute_per_prefix_coverage`
on the same inputs returns 1.0 for each target, so the summary carried two
contradictory figures and the wrong one was authoritative: the report headline
reads `fg_coverage`.

It reached the delivered pool as well as the report. `--minimize-primers` gates
on `metrics.fg_coverage >= target_coverage`, and with a ceiling near 1/N the
0.70 default is unreachable, so the trim silently stopped doing anything on
multi-target runs.

The same pooling applied to background positions, and in `exact` selectivity
mode the deduplicated counts feed `selectivity_ratio`.
"""

import numpy as np
import pytest

from neoswga.core.base_optimizer import BaseOptimizer, OptimizerConfig, OptimizationResult

TARGET_LENGTH = 100_000
REACH = 3_000
#: Dense enough that each target is covered end to end at the reach above.
#: Windows are 2*REACH wide, so a stride below that leaves no hole, and the
#: explicit tail site closes the last stretch the stride cannot reach.
SITES = sorted(set(range(REACH, TARGET_LENGTH, 2 * REACH - 500)) | {TARGET_LENGTH - REACH})


class _Cache:
    """Every primer binds at `SITES` in every prefix, in that prefix's own
    coordinate space -- which is precisely the situation that collapsed."""

    def get_positions(self, prefix, primer, strand="both"):
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.array(SITES, dtype=np.int64)


class _Optimizer(BaseOptimizer):
    """Minimal concrete subclass; only `compute_metrics` is under test."""

    @property
    def name(self):
        return "probe"

    @property
    def description(self):
        return "test double"

    def optimize(self, candidates, target_size=None, **kwargs):  # pragma: no cover
        raise NotImplementedError


def _optimizer(n_targets):
    prefixes = [f"fg{i}" for i in range(n_targets)]
    return _Optimizer(
        position_cache=_Cache(),
        fg_prefixes=prefixes,
        fg_seq_lengths=[TARGET_LENGTH] * n_targets,
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=2, extension_reach=REACH),
    )


PRIMERS = ["ACCACAGATAGC", "GTTGTAGATGGA"]


def test_one_target_is_fully_covered():
    """Guard the guard. If a single target is not fully covered, the
    multi-target comparison below proves nothing."""
    metrics = _optimizer(1).compute_metrics(PRIMERS)
    assert metrics.fg_coverage == pytest.approx(1.0, abs=0.01)


@pytest.mark.parametrize("n_targets", [2, 3, 4])
def test_covering_every_target_reports_full_coverage(n_targets):
    """N copies of a fully covered target is still full coverage.

    Pooling produced roughly 1/N instead, because the identical coordinates
    deduplicated against each other while the denominator kept growing.
    """
    metrics = _optimizer(n_targets).compute_metrics(PRIMERS)

    assert metrics.fg_coverage == pytest.approx(1.0, abs=0.01), (
        f"{n_targets} identical, fully covered targets reported "
        f"fg_coverage={metrics.fg_coverage:.4f}; positions from different "
        "genomes are being treated as one coordinate space"
    )


def test_a_gap_cannot_exceed_the_sequence_it_is_measured_in():
    """`max_gap` read 102,000 bp on a 100,000 bp target."""
    metrics = _optimizer(2).compute_metrics(PRIMERS)

    assert metrics.max_gap <= TARGET_LENGTH, (
        f"max_gap is {metrics.max_gap} on targets of {TARGET_LENGTH} bp; "
        "gaps are being measured across a genome boundary"
    )


def test_the_aggregate_agrees_with_the_per_target_figures():
    """Two numbers for one quantity, in one file, must not disagree.

    `per_target_coverage` is populated from `compute_per_prefix_coverage`, which
    was always right. `fg_coverage` is what the report reads.
    """
    from neoswga.core.coverage import compute_per_prefix_coverage

    optimizer = _optimizer(3)
    metrics = optimizer.compute_metrics(PRIMERS)

    _overall, per_target = compute_per_prefix_coverage(
        cache=optimizer.cache,
        primers=PRIMERS,
        prefixes=optimizer.fg_prefixes,
        seq_lengths=optimizer.fg_seq_lengths,
        extension=REACH,
        circular=False,
    )
    expected = sum(per_target.values()) / len(per_target)

    assert metrics.fg_coverage == pytest.approx(expected, abs=0.01), (
        f"fg_coverage={metrics.fg_coverage:.4f} against a per-target mean of "
        f"{expected:.4f} ({per_target})"
    )
