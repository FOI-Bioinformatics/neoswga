"""Regression tests: BackgroundAwareOptimizer._calculate_coverage must apply
strand-direction-aware extension reach.

The previous implementation extended every binding site bidirectionally
([pos-reach, pos+reach]) regardless of strand. phi29 polymerase extends
from the primer 3' end in one direction along the template, so a forward
primer extends downstream only and a reverse primer extends upstream
only. The bidirectional model over-estimated coverage by up to 2x in
strand-imbalanced primer sets.
"""

from unittest.mock import MagicMock

import pytest

from neoswga.core.background_aware_optimizer import BackgroundAwareOptimizer


def _build_optimizer(positions_for, fg_seq_length=100_000):
    """positions_for: dict (prefix, primer, strand) -> list[int]."""
    cache = MagicMock()
    cache.get_positions.side_effect = lambda prefix, primer, strand: (
        list(positions_for.get((prefix, primer, strand), []))
    )
    return BackgroundAwareOptimizer(
        position_cache=cache,
        fg_prefixes=["fg"],
        bg_prefixes=[],
        fg_seq_lengths=[fg_seq_length],
        bg_seq_lengths=[],
    )


def test_forward_only_site_extends_downstream_only():
    """A forward-strand primer at position 5000 with reach 3000 should
    cover [5000, 8000) only (downstream), not [2000, 8000) bidirectionally.
    """
    opt = _build_optimizer({
        ("fg", "AAAACCCC", "forward"): [5000],
        ("fg", "AAAACCCC", "reverse"): [],
        ("fg", "AAAACCCC", "both"): [5000],
    })
    coverage = opt._calculate_coverage(["AAAACCCC"], extension_reach=3000)
    # Bidirectional (buggy): 6000 / 100000 = 0.06
    # Strand-aware (correct): 3000 / 100000 = 0.03
    assert coverage == pytest.approx(0.03, abs=1e-3)


def test_reverse_only_site_extends_upstream_only():
    """A reverse-strand primer at position 5000 with reach 3000 should
    cover [2000, 5000) only (upstream).
    """
    opt = _build_optimizer({
        ("fg", "GGGGTTTT", "forward"): [],
        ("fg", "GGGGTTTT", "reverse"): [5000],
        ("fg", "GGGGTTTT", "both"): [5000],
    })
    coverage = opt._calculate_coverage(["GGGGTTTT"], extension_reach=3000)
    assert coverage == pytest.approx(0.03, abs=1e-3)


def test_forward_and_reverse_sites_at_same_position_cover_both_directions():
    """A site with both forward and reverse primers at position 5000
    should cover [2000, 8000), the union of upstream and downstream.
    """
    opt = _build_optimizer({
        ("fg", "AAAACCCC", "forward"): [5000],
        ("fg", "AAAACCCC", "reverse"): [5000],
        ("fg", "AAAACCCC", "both"): [5000, 5000],
    })
    coverage = opt._calculate_coverage(["AAAACCCC"], extension_reach=3000)
    assert coverage == pytest.approx(0.06, abs=1e-3)


def test_no_primers_returns_zero():
    opt = _build_optimizer({})
    assert opt._calculate_coverage([], extension_reach=3000) == 0.0
