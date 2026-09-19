"""A strand score of 0.0 meant two different things, and only one was a number.

Known Issue 18, the half left open. `PrimerSetMetrics.strand_alternation_score`
and `strand_coverage_ratio` were made `Optional` and the aggregation in
`core/strand_metrics.py` returns `(None, None)` when no foreground genome was
measured. That closed the metrics side.

The SOURCE still conflated them. `PositionCache.compute_strand_alternation_stats`
returns 0.0 for both whenever the panel has fewer than two recorded sites --
and for alternation that is not a measured zero, it is the absence of one.
Alternation is the fraction of ADJACENT site pairs that sit on opposite
strands. With fewer than two sites there is no adjacent pair, so the fraction
is 0/0. A panel whose two sites are both forward genuinely scores 0.0, and
those two cases were indistinguishable by value.

That ambiguity is the stated reason these two are not constrainable: a limit on
them would reject a panel for a missing measurement while reporting a violated
constraint. This file pins the distinction at the source.

`strand_coverage_ratio` is treated differently on purpose, because the
arithmetic differs. It is min(fwd, rev) / max(fwd, rev), which needs one site,
not two. A single forward site really is maximally unbalanced, so 0.0 there is
a measurement. Only the no-sites case is undefined.

The gap figures are deliberately left alone. `strand_alternation_gap_mean` and
`_max` return the genome length when there is no opposite-strand pair, which is
a conservative encoding of "no convergent pair anywhere" rather than a missing
measurement, and `worst_convergent_gap` reads them that way.
"""

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.position_cache import PositionCache

GENOME = 50_000
K = 10


def _cache(tmp_path, sites_by_primer):
    """A cache whose primers bind exactly where this says."""
    prefix = str(tmp_path / "t")
    primers = list(sites_by_primer)
    with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as db:
        for primer, positions in sites_by_primer.items():
            db.create_dataset(primer, data=np.array(positions, dtype=np.int64))
        db.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))
    return PositionCache([prefix], primers), prefix


class TestAlternationSaysWhenItCannotBeMeasured:
    def test_no_sites_at_all_is_none(self, tmp_path):
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": []})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_alternation_score"] is None

    def test_one_site_is_none_not_zero(self, tmp_path):
        """The case the whole entry turns on. One site has no adjacent pair, so
        the fraction is 0/0 -- not a panel that failed to alternate."""
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": [1000]})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_alternation_score"] is None

    def test_two_sites_on_one_strand_is_a_measured_zero(self, tmp_path):
        """The other half: this panel really does not alternate, and that must
        stay distinguishable from the case above."""
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": [1000, 2000]})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_alternation_score"] == 0.0

    def test_a_measured_zero_is_not_none(self, tmp_path):
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": [1000, 2000]})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_alternation_score"] is not None


class TestTheCoverageRatioNeedsOnlyOneSite:
    """Different arithmetic, so a different rule. min/max needs one site."""

    def test_no_sites_is_none(self, tmp_path):
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": []})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_coverage_ratio"] is None

    def test_one_site_is_a_measured_zero(self, tmp_path):
        """A single forward site really is maximally unbalanced. That is a
        number, and reporting None would lose it."""
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": [1000]})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_coverage_ratio"] == 0.0


class TestTheGapsAreLeftAsTheyWere:
    def test_a_panel_with_no_convergent_pair_still_reports_the_genome_length(self, tmp_path):
        """Not a missing measurement: it says no convergent pair exists
        anywhere, which is what `worst_convergent_gap` reads it as."""
        cache, prefix = _cache(tmp_path, {"ACGTACGTAC": [1000]})

        stats = cache.compute_strand_alternation_stats(prefix, ["ACGTACGTAC"], GENOME)

        assert stats["strand_alternation_gap_max"] == float(GENOME)
        assert stats["strand_alternation_gap_mean"] == float(GENOME)


class TestNothingDownstreamBreaksOnNone:
    def test_the_aggregate_passes_none_through(self, tmp_path):
        from neoswga.core.strand_metrics import headline_strand_scalars

        stats = {"p": {"strand_alternation_score": None, "strand_coverage_ratio": None}}

        assert headline_strand_scalars(stats, ["p"]) == (None, None)

    def test_the_regime_report_does_not_crash_on_none(self):
        from types import SimpleNamespace

        from neoswga.core.panel_regime import assess_panel

        metrics = SimpleNamespace(
            fg_coverage=0.5,
            effective_fg_coverage=0.5,
            selectivity_density=10.0,
            total_bg_sites=5,
            max_gap=1000,
            bg_coverage=0.01,
            strand_coverage_ratio=None,
            strand_alternation_score=None,
        )

        regime = assess_panel(
            metrics,
            coverage_target=0.5,
            min_fg_bg_ratio=1.0,
            requested_size=6,
            delivered_size=6,
            coverage_reach=3000,
        )

        assert regime is not None
