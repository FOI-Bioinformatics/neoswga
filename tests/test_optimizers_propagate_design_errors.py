"""A `DesignError` raised inside an optimizer reaches the command boundary.

`hybrid`, `background-aware` and `network` wrapped their whole `optimize()` in
`except Exception` and returned `OptimizationResult.failure`. A missing
reference answer or a failed model evaluation therefore became "this method
lost", and an ensemble returned another method's panel built on the same
unmeasured data (2026-11-23 sweep, tier 1 item 12). `dominating-set`, `clique`
and the ensemble member loop already re-raise.
"""

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.background_aware_optimizer import BackgroundAwareBaseOptimizer
from neoswga.core.base_optimizer import (
    CompositeOptimizer,
    OptimizationResult,
    OptimizationStatus,
    PrimerSetMetrics,
)
from neoswga.core.exceptions import ModelEvaluationError, ReferenceDataError
from neoswga.core.hybrid_optimizer import HybridBaseOptimizer
from neoswga.core.network_optimizer import NetworkBaseOptimizer
from neoswga.core.position_cache import PositionCache

PRIMERS = ["ACGTACGTAC", "TTGACCATGA", "GGCATTACGA"]
LENGTH = 10_000


@pytest.fixture
def cache_and_prefix(tmp_path):
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
        for i, primer in enumerate(PRIMERS):
            f.create_dataset(primer, data=np.array([100 + i * 3000], dtype=np.int64))
    return PositionCache([prefix], PRIMERS), prefix


def _build(cls, cache_and_prefix):
    cache, prefix = cache_and_prefix
    return cls(
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[LENGTH],
        bg_prefixes=[prefix],
        bg_seq_lengths=[LENGTH],
    )


def _raise(*_args, **_kwargs):
    raise ReferenceDataError("target", "the index was never asked")


@pytest.mark.parametrize(
    "cls,inner,method",
    [
        (HybridBaseOptimizer, "_hybrid", "optimize"),
        (BackgroundAwareBaseOptimizer, "_hybrid", "optimize"),
        (NetworkBaseOptimizer, "_network", "optimize_greedy"),
    ],
)
def test_optimizer_reraises_design_error(cls, inner, method, cache_and_prefix, monkeypatch):
    optimizer = _build(cls, cache_and_prefix)
    monkeypatch.setattr(getattr(optimizer, inner), method, _raise)
    with pytest.raises(ReferenceDataError):
        optimizer.optimize(PRIMERS, target_size=2)


def test_composite_reraises_a_member_design_error(cache_and_prefix):
    class _Fails:
        name = "fails"

        def optimize(self, *_args, **_kwargs):
            raise ModelEvaluationError("fg_coverage", ["target"], "no answer")

    class _Succeeds:
        name = "succeeds"

        def optimize(self, *_args, **_kwargs):
            return OptimizationResult(
                primers=tuple(PRIMERS[:2]),
                score=1.0,
                status=OptimizationStatus.SUCCESS,
                metrics=PrimerSetMetrics.empty(),
                iterations=1,
                optimizer_name="succeeds",
            )

    cache, prefix = cache_and_prefix
    composite = CompositeOptimizer(
        [_Fails(), _Succeeds()], cache, [prefix], [LENGTH]
    )
    with pytest.raises(ModelEvaluationError):
        composite.optimize(PRIMERS, target_size=2)


def test_a_design_error_building_the_optimizer_reaches_the_boundary(monkeypatch):
    """`_dispatch_optimizer` turned any construction failure into a failure
    result, so a DesignError lost its type and fields and the record that
    followed said only "No primer sets found" (2026-11-23 review, M5)."""
    from neoswga.core import unified_optimizer

    def _refuse(*_args, **_kwargs):
        raise ReferenceDataError("host", "never indexed")

    monkeypatch.setattr(unified_optimizer.OptimizerFactory, "create", staticmethod(_refuse))
    with pytest.raises(ReferenceDataError):
        unified_optimizer._dispatch_optimizer(
            "hybrid", None, PRIMERS, ["t"], [LENGTH], [], [], 2,
            None, None, None, False, "balanced", {},
        )
