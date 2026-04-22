"""Coverage wrap-around on closed-circular sequences.

Regression test for the Phase 16/17 reach fix: with realistic per-primer
reach (3 kb phi29, 4 kb equiphi29) a single binding site on a small
circular plasmid or bacterial chromosome must cover the region on both
sides of the site even when the window crosses the origin. Without
wrap-around a site at position 500 on a 6.2 kb circular plasmid covers
only [0, 3500] = 56% instead of the ~97% that circular topology implies.

Exercises ``coverage._mark_window``, ``coverage.compute_per_prefix_coverage``,
and ``base_optimizer.BaseOptimizer._compute_coverage`` — the three places
the wrap is implemented.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# _mark_window — the primitive
# ---------------------------------------------------------------------------


def _fresh(length: int) -> np.ndarray:
    return np.zeros(length, dtype=bool)


def test_mark_window_linear_mid_sequence_is_clipped_but_contiguous():
    from neoswga.core.coverage import _mark_window

    occ = _fresh(6200)
    _mark_window(occ, pos=3000, extension=3000, length=6200, circular=False)
    # Fully contained inside the genome, no clipping triggered.
    assert occ[0:6000].all()
    assert not occ[6000:].any()
    assert occ.sum() == 6000


def test_mark_window_linear_near_origin_is_clipped_no_wrap():
    from neoswga.core.coverage import _mark_window

    occ = _fresh(6200)
    _mark_window(occ, pos=500, extension=3000, length=6200, circular=False)
    # Linear: left side clipped to 0, no wrap-around at the far end.
    assert occ[0:3500].all()
    assert not occ[3500:].any()


def test_mark_window_circular_near_origin_wraps():
    from neoswga.core.coverage import _mark_window

    occ = _fresh(6200)
    _mark_window(occ, pos=500, extension=3000, length=6200, circular=True)
    # Circular: left side wraps to [6200 - 2500, 6200] = [3700, 6200]
    # and right side is [0, 3500]. The gap in [3500, 3700] stays empty.
    assert occ[0:3500].all()
    assert not occ[3500:3700].any()
    assert occ[3700:6200].all()
    assert occ.sum() == 6000  # 3500 + 2500


def test_mark_window_circular_near_end_wraps():
    from neoswga.core.coverage import _mark_window

    occ = _fresh(6200)
    _mark_window(occ, pos=5700, extension=3000, length=6200, circular=True)
    # Right side wraps to [0, (5700+3000)-6200] = [0, 2500]
    # Left side is [2700, 6200].
    assert occ[0:2500].all()
    assert not occ[2500:2700].any()
    assert occ[2700:6200].all()


def test_mark_window_circular_window_exceeds_genome_covers_all():
    from neoswga.core.coverage import _mark_window

    occ = _fresh(4000)
    _mark_window(occ, pos=2000, extension=3000, length=4000, circular=True)
    # Window diameter 6000 > genome 4000, so any position covers everything.
    assert occ.all()


# ---------------------------------------------------------------------------
# compute_per_prefix_coverage — coverage helper used for per_target_coverage
# ---------------------------------------------------------------------------


class _StubCache:
    """Minimal PositionCache-shaped object for the coverage helper.

    Returns positions for a single primer on a single prefix.
    """

    def __init__(self, positions: List[int]):
        self._positions = np.asarray(positions, dtype=int)

    def get_positions(self, prefix, primer, strand):  # noqa: ARG002
        return self._positions


def test_compute_per_prefix_coverage_linear_caps_at_clip_edge():
    from neoswga.core.coverage import compute_per_prefix_coverage

    cache = _StubCache([500])
    agg, per = compute_per_prefix_coverage(
        cache=cache,
        primers=["ACGTACGTACGTACGTA"],
        prefixes=["pl"],
        seq_lengths=[6200],
        extension=3000,
        circular=False,
    )
    # Single site at 500 on a linear 6.2 kb: [0, 3500] = 3500 bp = 56.5%.
    assert per["pl"] == pytest.approx(3500 / 6200, rel=1e-6)
    assert agg == pytest.approx(3500 / 6200, rel=1e-6)


def test_compute_per_prefix_coverage_circular_wraps_to_high_coverage():
    from neoswga.core.coverage import compute_per_prefix_coverage

    cache = _StubCache([500])
    agg, per = compute_per_prefix_coverage(
        cache=cache,
        primers=["ACGTACGTACGTACGTA"],
        prefixes=["pl"],
        seq_lengths=[6200],
        extension=3000,
        circular=True,
    )
    # Circular: [0, 3500] U [3700, 6200] = 3500 + 2500 = 6000 bp = 96.8%.
    assert per["pl"] == pytest.approx(6000 / 6200, rel=1e-6)


# ---------------------------------------------------------------------------
# base_optimizer._compute_coverage — the optimizer objective
# ---------------------------------------------------------------------------


@dataclass
class _FakeConfig:
    extension_reach: int = 3000
    fg_circular: bool = False


class _BareOptimizer:
    """Constructs just enough of a BaseOptimizer-alike to exercise the
    coverage computation without instantiating the full subclass tree.
    """

    def __init__(self, circular: bool, reach: int = 3000):
        self.config = _FakeConfig(extension_reach=reach, fg_circular=circular)


def _compute(circular: bool, positions: List[int], length: int,
             reach: int = 3000) -> float:
    from neoswga.core.base_optimizer import BaseOptimizer

    opt = _BareOptimizer(circular=circular, reach=reach)
    return BaseOptimizer._compute_coverage(opt, positions, length)


def test_base_optimizer_linear_matches_pre_fix_behavior():
    # Pre-fix / linear: single mid-sequence site fully contained.
    cov = _compute(circular=False, positions=[3000], length=6200)
    assert cov == pytest.approx(6000 / 6200, rel=1e-6)


def test_base_optimizer_linear_clips_near_edge():
    cov = _compute(circular=False, positions=[500], length=6200)
    # 0..3500 = 3500 bp
    assert cov == pytest.approx(3500 / 6200, rel=1e-6)


def test_base_optimizer_circular_wraps_near_edge():
    cov = _compute(circular=True, positions=[500], length=6200)
    # Wrap: [0, 3500] U [3700, 6200] = 6000
    assert cov == pytest.approx(6000 / 6200, rel=1e-6)


def test_base_optimizer_circular_window_exceeds_genome_is_full():
    cov = _compute(circular=True, positions=[2000], length=4000, reach=3000)
    assert cov == 1.0


def test_base_optimizer_circular_matches_linear_when_no_window_wraps():
    # Both sites are >= reach bp from either end, so the wrap branch
    # never triggers and circular/linear paths must agree exactly.
    # pos=4000 gives window [1000, 7000]; pos=6000 gives [3000, 9000].
    cov_lin = _compute(circular=False, positions=[4000, 6000], length=10000)
    cov_cir = _compute(circular=True, positions=[4000, 6000], length=10000)
    assert cov_lin == pytest.approx(cov_cir, rel=1e-9)
    # Union [1000, 9000] = 8000 bp = 80%.
    assert cov_lin == pytest.approx(0.80, rel=1e-9)


def test_base_optimizer_circular_recovers_bases_that_linear_clips():
    # Two sites near opposite edges. Linear clips; circular wraps and
    # recovers the bases that clipping would have discarded.
    # pos=2000 reach=3000 on 10000 linear: [0, 5000].
    # pos=2000 reach=3000 on 10000 circular: [0, 5000] U [9000, 10000].
    # pos=5000 reach=3000: [2000, 8000] (no wrap either way).
    # Linear union: [0, 8000] = 8000 bp = 0.80
    # Circular union: [0, 8000] U [9000, 10000] = 9000 bp = 0.90
    cov_lin = _compute(circular=False, positions=[2000, 5000], length=10000)
    cov_cir = _compute(circular=True, positions=[2000, 5000], length=10000)
    assert cov_lin == pytest.approx(0.80, rel=1e-9)
    assert cov_cir == pytest.approx(0.90, rel=1e-9)
