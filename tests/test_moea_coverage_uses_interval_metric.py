"""Regression test: MOEA _calculate_coverage must use interval coverage.

The previous implementation counted unique binding positions divided by
genome length, treating each primer site as 1 nt. This produced coverage
fractions ~10000x smaller than the interval-based metric used elsewhere
in the codebase, so the multi-objective optimizer effectively ignored
the coverage objective and selected primers on other criteria only.
"""

from unittest.mock import MagicMock

import pytest

pytest.importorskip("pymoo")  # MOEA module gates on pymoo availability

from neoswga.core.moea_optimizer import PrimerSetProblem


def _build_problem(positions_for):
    """positions_for: dict (prefix, primer, strand) -> list[int]."""
    cache = MagicMock()
    cache.get_positions.side_effect = lambda prefix, primer, strand: (
        list(positions_for.get((prefix, primer, strand), []))
    )
    return PrimerSetProblem(
        candidates=["AAAACCCC", "GGGGTTTT"],
        cache=cache,
        fg_prefixes=["fg"],
        bg_prefixes=[],
        fg_seq_lengths=[100_000],
        bg_seq_lengths=[],
        max_primers=10,
    )


def test_single_site_coverage_is_interval_not_point():
    """One primer with one binding site at position 5000 should cover
    [5000-3000, 5000+3000) = 6000 bp out of 100 kb = 0.06.
    """
    problem = _build_problem({
        ("fg", "AAAACCCC", "both"): [5000],
        ("fg", "AAAACCCC", "forward"): [5000],
        ("fg", "AAAACCCC", "reverse"): [],
    })
    coverage = problem._calculate_coverage(
        ["AAAACCCC"], ["fg"], [100_000]
    )
    # Old (buggy) behaviour: 1 / 100_000 = 1e-5
    # New (interval) behaviour: 6000 / 100_000 = 0.06
    assert coverage == pytest.approx(0.06, abs=1e-3)


def test_overlapping_sites_are_unioned_not_summed():
    """Two binding sites within 2*extension of each other must not
    double-count the overlapping region.
    """
    problem = _build_problem({
        ("fg", "AAAACCCC", "both"): [5000, 7000],
        ("fg", "AAAACCCC", "forward"): [5000, 7000],
        ("fg", "AAAACCCC", "reverse"): [],
    })
    coverage = problem._calculate_coverage(
        ["AAAACCCC"], ["fg"], [100_000]
    )
    # Naive sum: 2 * 6000 = 12000 -> 0.12
    # Correct union: [2000, 8000) U [4000, 10000) = [2000, 10000) = 8000 -> 0.08
    assert coverage == pytest.approx(0.08, abs=1e-3)


def test_no_primers_returns_zero():
    problem = _build_problem({})
    assert problem._calculate_coverage([], ["fg"], [100_000]) == 0.0


def test_no_positions_returns_zero():
    problem = _build_problem({})
    coverage = problem._calculate_coverage(
        ["AAAACCCC"], ["fg"], [100_000]
    )
    assert coverage == 0.0
