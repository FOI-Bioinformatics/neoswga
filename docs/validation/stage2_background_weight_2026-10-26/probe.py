"""Probe for the swap Stage-2 host term as a weighted axis (Known Issue 14).

Run from the repo root:

    PYTHONPATH=. python docs/validation/stage2_background_weight_2026-10-26/probe.py

It shows three things, with no external tool:

1. Off by default (weight 0): the swap score is the lexicographic tuple it has
   always been, so a lower-coverage swap is never taken however much host
   binding it would shed.
2. On (weight > 0): the same swap IS taken, because foreground coverage is now
   traded against host coverage on one axis.
3. The weights reach the inner HybridOptimizer through the factory/CLI path --
   the gap that had left `stage2_tm_weight` inert on every CLI hybrid run
   (Known Issue 8) before the wrapper was taught to forward both.
"""

from types import SimpleNamespace

from neoswga.core.swap_refinement import refine_by_swaps


class _Compat:
    def dimerises(self, primer, others):
        return False


class _TradeObjective:
    # primer -> (fg_coverage, bg_coverage, total_bg_sites)
    _T = {"A": (0.60, 0.90, 9), "B": (0.50, 0.10, 1)}

    def shortfall(self, panel):
        return 0.0

    def coverage(self, panel):
        return max(self._T[p][0] for p in panel)

    def metrics(self, panel):
        return SimpleNamespace(
            bg_coverage=max(self._T[p][1] for p in panel),
            total_bg_sites=sum(self._T[p][2] for p in panel),
        )


def _swap(weight):
    return refine_by_swaps(
        ["A"], ["A", "B"], None, None, _Compat(),
        objective=_TradeObjective(), background_weight=weight,
    ).primers


def main():
    print("Scenario: panel starts as ('A',), candidate B is available.")
    print("  A: fg_coverage=0.60  bg_coverage=0.90")
    print("  B: fg_coverage=0.50  bg_coverage=0.10  (covers less, binds the host far less)")
    print()
    off = _swap(0.0)
    print(f"weight 0.0 (default): delivered {off}  -- lexicographic coverage keeps A")
    assert off == ("A",), off
    on = _swap(0.8)
    print(f"weight 0.8          : delivered {on}  -- weighted axis prefers B")
    assert on == ("B",), on
    print("    A: 0.60 - 0.8*0.90 = -0.12")
    print("    B: 0.50 - 0.8*0.10 =  0.42  (wins)")
    print()

    # The forwarding closure. Build a hybrid through the real factory with both
    # weights set on the config and read them off the inner engine.
    import numpy as np

    from neoswga.core import unified_optimizer as uo
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.optimizer_factory import OptimizerFactory

    uo._ensure_optimizers_registered()

    class _Cache:
        def get_positions(self, *args, **kwargs):
            return np.array([], dtype=np.int64)

    config = OptimizerConfig(target_set_size=4)
    config.stage2_tm_weight = 1.0
    config.stage2_background_weight = 0.5
    opt = OptimizerFactory.create(
        name="hybrid",
        position_cache=_Cache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[10_000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=config,
        conditions=None,
    )
    inner = opt._hybrid
    print(f"factory hybrid inner.stage2_tm_weight         = {inner.stage2_tm_weight}")
    print(f"factory hybrid inner.stage2_background_weight = {inner.stage2_background_weight}")
    assert inner.stage2_tm_weight == 1.0, inner.stage2_tm_weight
    assert inner.stage2_background_weight == 0.5, inner.stage2_background_weight
    print()
    print("OK: off is unchanged, on trades coverage for host, and both weights reach the inner.")


if __name__ == "__main__":
    main()
