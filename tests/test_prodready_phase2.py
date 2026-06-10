"""Phase 2 (production-readiness v2): coverage semantics & counting correctness."""

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Hybrid: Stage-1 coverage selection uses the realistic coverage_reach,
# distinct from the processivity max_extension used for connectivity.
# ---------------------------------------------------------------------------


def test_hybrid_splits_coverage_reach_from_processivity():
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    cache = _StubCache()
    opt = HybridOptimizer(
        position_cache=cache,
        fg_prefixes=["fg"],
        fg_seq_lengths=[1_000_000],
        max_extension=70000,  # processivity (Stage-2 connectivity)
        coverage_reach=3000,  # realistic per-primer reach (Stage-1 coverage)
        polymerase="phi29",
    )
    assert opt.coverage_reach == 3000
    assert opt.max_extension == 70000
    # Stage-1 set-cover optimizer must select at the realistic coverage reach.
    assert opt.dominating_optimizer.extension_reach == 3000
    # Stage-2 connectivity network keeps processivity.
    assert opt.network_optimizer.max_extension == 70000


def test_hybrid_coverage_reach_defaults_to_max_extension():
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    opt = HybridOptimizer(
        position_cache=_StubCache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[1000],
        max_extension=5000,  # no coverage_reach -> backward-compatible fallback
    )
    assert opt.coverage_reach == 5000


class _StubCache:
    def get_positions(self, prefix, primer, strand="both"):
        return np.array([], dtype=np.int64)


# ---------------------------------------------------------------------------
# PositionCache.get_positions('both') dedups palindromic double-counts.
# ---------------------------------------------------------------------------


def test_position_cache_both_dedups_palindrome():
    from neoswga.core.position_cache import PositionCache

    pc = PositionCache.__new__(PositionCache)  # bypass HDF5 loading
    pc.cache = {}
    # A palindromic primer: forward and reverse keys hold the SAME positions.
    sites = np.array([10, 50, 90], dtype=np.int32)
    pc.cache[("fg", "AATT", "forward")] = sites
    pc.cache[("fg", "AATT", "reverse")] = sites

    both = pc.get_positions("fg", "AATT", "both")
    assert sorted(both.tolist()) == [10, 50, 90]  # not 6 (deduped)


def test_position_cache_both_keeps_distinct_sites():
    from neoswga.core.position_cache import PositionCache

    pc = PositionCache.__new__(PositionCache)
    pc.cache = {}
    pc.cache[("fg", "ATCG", "forward")] = np.array([10, 20], dtype=np.int32)
    pc.cache[("fg", "ATCG", "reverse")] = np.array([30, 40], dtype=np.int32)

    both = pc.get_positions("fg", "ATCG", "both")
    assert sorted(both.tolist()) == [10, 20, 30, 40]


# ---------------------------------------------------------------------------
# NetworkOptimizer Tm scoring engages when conditions carry a temperature.
# ---------------------------------------------------------------------------


def test_network_reaction_temp_falls_back_to_conditions():
    from neoswga.core.network_optimizer import NetworkOptimizer
    from neoswga.core.reaction_conditions import ReactionConditions

    cond = ReactionConditions(temp=42.0, polymerase="equiphi29")
    opt = NetworkOptimizer(
        position_cache=_StubCache(),
        fg_prefixes=["fg"],
        bg_prefixes=[],
        fg_seq_lengths=[1000],
        bg_seq_lengths=[],
        conditions=cond,
    )
    assert opt.reaction_temp == 42.0
    # Tm score is no longer a constant 1.0 no-op when a hot primer is off-target.
    s_hot = opt._calculate_tm_score("GGGGCCCCGGGGCCCC")
    s_balanced = opt._calculate_tm_score("ATCGATCGATCG")
    assert s_hot != s_balanced


def test_network_reaction_temp_explicit_wins():
    from neoswga.core.network_optimizer import NetworkOptimizer
    from neoswga.core.reaction_conditions import ReactionConditions

    cond = ReactionConditions(temp=42.0, polymerase="equiphi29")
    opt = NetworkOptimizer(
        position_cache=_StubCache(),
        fg_prefixes=["fg"],
        bg_prefixes=[],
        fg_seq_lengths=[1000],
        bg_seq_lengths=[],
        conditions=cond,
        reaction_temp=30.0,
    )
    assert opt.reaction_temp == 30.0


# ---------------------------------------------------------------------------
# filter.get_all_rates zero-length guard (no ZeroDivisionError).
# ---------------------------------------------------------------------------


def test_get_all_rates_zero_length_no_crash(monkeypatch):
    from neoswga.core import filter as filt

    monkeypatch.setattr(
        filt,
        "get_rates_for_one_species",
        lambda primer_list, prefixes, **kw: {p: 5 for p in primer_list},
    )

    class _P:
        min_fg_freq = 1e-5
        max_bg_freq = 5e-6

    monkeypatch.setattr(filt, "parameter", _P)

    df = filt.get_all_rates(["ATCGATCG"], ["fg"], [], fg_total_length=0, bg_total_length=0)
    assert len(df) == 1  # no ZeroDivisionError


# ---------------------------------------------------------------------------
# coverage.compute_per_prefix_coverage: narrow exception handling.
# ---------------------------------------------------------------------------


def test_coverage_skips_absent_primer_but_propagates_io_error():
    from neoswga.core.coverage import compute_per_prefix_coverage

    class _KeyErrCache:
        def get_positions(self, prefix, primer, strand):
            raise KeyError(primer)  # genuinely-absent primer -> skipped

    agg, per = compute_per_prefix_coverage(_KeyErrCache(), ["ATCG"], ["fg"], [1000], extension=100)
    assert agg == 0.0

    class _IOErrCache:
        def get_positions(self, prefix, primer, strand):
            raise OSError("corrupt HDF5")  # systemic failure must surface

    with pytest.raises(OSError):
        compute_per_prefix_coverage(_IOErrCache(), ["ATCG"], ["fg"], [1000], extension=100)
