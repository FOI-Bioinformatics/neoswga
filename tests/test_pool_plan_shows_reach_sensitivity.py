"""A coverage figure quoted without its reach carries almost no information.

Audit finding F2. Effective coverage marks symmetric windows around exact
binding sites, so its value is roughly linear in the window radius over the
range that matters. On one saved 26-oligo wMel panel under identical conditions
the same design reads 41.3% at 1 kb, 80.1% at 3 kb, 93.5% at 5 kb and 99.6% at
10 kb.

The 3 kb default is a design-density convention taken from sets with measured
wet-lab success, not a measured extension distribution. So a pool-size
recommendation that prints one coverage number invites the reader to treat a
modelling choice as a result.

`plan-pool` now computes the same panel at several reaches and saves the sweep,
and both reports render it. The audit's action was "always show reach
sensitivity"; this is that.
"""

import pytest

from neoswga.core.pool_planner import REACH_SENSITIVITY_FACTORS, reach_sensitivity


class _Cache:
    """Every primer binds once, in the middle of a 10 kb target."""

    def get_positions(self, prefix, primer, strand):  # noqa: ARG002
        import numpy as np

        return np.array([5_000], dtype=int)


def test_the_sweep_spans_the_configured_reach():
    rows = reach_sensitivity(
        cache=_Cache(),
        primers=["ACGTACGTACGT"],
        prefixes=["t"],
        seq_lengths=[10_000],
        reach=1_000,
    )

    reaches = [r["reach"] for r in rows]
    assert 1_000 in reaches, "the configured reach must appear in its own sweep"
    assert reaches == sorted(reaches)
    assert len(reaches) == len(set(reaches))


def test_coverage_rises_with_reach():
    rows = reach_sensitivity(
        cache=_Cache(),
        primers=["ACGTACGTACGT"],
        prefixes=["t"],
        seq_lengths=[10_000],
        reach=1_000,
    )

    coverages = [r["coverage"] for r in rows]
    assert coverages == sorted(coverages)
    assert coverages[0] < coverages[-1], "a flat sweep would tell the reader nothing"


def test_a_single_site_at_one_kb_covers_a_fifth_of_ten_kb():
    """Anchored arithmetic, so the sweep is not merely self-consistent."""
    rows = reach_sensitivity(
        cache=_Cache(),
        primers=["ACGTACGTACGT"],
        prefixes=["t"],
        seq_lengths=[10_000],
        reach=1_000,
    )
    at_1kb = next(r for r in rows if r["reach"] == 1_000)

    # One site at 5000, window [4000, 6000) on a 10 kb target.
    assert at_1kb["coverage"] == pytest.approx(0.2)


def test_the_factors_bracket_the_configured_value():
    assert min(REACH_SENSITIVITY_FACTORS) < 1.0 < max(REACH_SENSITIVITY_FACTORS)
    assert 1.0 in REACH_SENSITIVITY_FACTORS
