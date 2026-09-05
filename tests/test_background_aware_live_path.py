"""The `background-aware` method the CLI actually runs, tested on that path.

Two existing files -- `test_background_aware_reach_is_configured.py` and
`test_background_aware_strand_extension.py` -- construct
`BackgroundAwareOptimizer`, the standalone three-stage class. No CLI run reaches
it. `--optimization-method background-aware` resolves to
`BackgroundAwareBaseOptimizer`, which delegates to `HybridOptimizer`, and the
standalone class is referenced only from tests. So those files assert on a
`_calculate_coverage` that no run executes, and passed throughout the period
when the live path was dropping the configured Tm window, the extension reach
and every selection weight.

This file asserts the same properties of the object the factory returns. It is
the difference between "the algorithm is right" and "the algorithm the user gets
is right", and this repository's recurring defect lives in the gap between them.
"""

import pytest

pytest.importorskip("networkx")

import numpy as np

from neoswga.core import unified_optimizer as uo
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory

GENOME = 200_000


class _Cache:
    def get_positions(self, prefix, primer, strand="both"):
        return np.array([], dtype=np.int64)


@pytest.fixture(scope="module", autouse=True)
def _registered():
    uo._ensure_optimizers_registered()


def _live(**kwargs):
    """The optimizer `--optimization-method background-aware` produces."""
    config = kwargs.pop("config", None) or OptimizerConfig(target_set_size=6)
    return OptimizerFactory.create(
        name="background-aware",
        position_cache=_Cache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[GENOME],
        bg_prefixes=["bg"],
        bg_seq_lengths=[GENOME],
        config=config,
        conditions=None,
        **kwargs,
    )


def test_the_registered_method_delegates_to_the_hybrid_optimizer():
    """Naming the live structure, so a future change to it fails here first."""
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    optimizer = _live()
    assert isinstance(optimizer._hybrid, HybridOptimizer), (
        "background-aware no longer delegates to HybridOptimizer; the "
        "properties below are asserted of the wrong object"
    )


def test_background_pruning_is_switched_on():
    """It is the stage the method exists for."""
    assert _live()._hybrid.background_pruning is True


@pytest.mark.parametrize("reach", [1_000, 3_000, 8_000])
def test_the_configured_coverage_reach_reaches_the_live_optimizer(reach):
    optimizer = _live(config=OptimizerConfig(target_set_size=6, extension_reach=reach))
    assert optimizer._hybrid.coverage_reach == reach


def test_the_configured_tm_window_reaches_the_live_optimizer():
    """The sharpest of the dropped values: on an equiphi29 pool configured for
    20-70 C, this wrapper screened at the preset 37-62 and discarded two thirds
    of the candidates before selection."""
    optimizer = _live(config=OptimizerConfig(target_set_size=6, min_tm=20.0, max_tm=70.0))
    assert optimizer._hybrid.min_tm == 20.0
    assert optimizer._hybrid.max_tm == 70.0


def test_the_pruning_floor_is_relative_on_the_live_optimizer():
    """An absolute 0.95 floor made the pruning stage a no-op at realistic
    coverage, which is to say on every ordinary run."""
    optimizer = _live()
    assert optimizer._hybrid._absolute_coverage_floor is False
    assert 0.0 < optimizer._hybrid.min_coverage_drop < 1.0
