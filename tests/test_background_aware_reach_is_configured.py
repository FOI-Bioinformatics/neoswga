"""The background-aware optimizer's reach was a literal, not a setting.

`BackgroundAwareOptimizer` hardcoded 3000 bp in two places -- the
`DominatingSetOptimizer` it builds for Stage 1, and the default argument of
`_calculate_coverage`, which every caller then took. So the coverage this
optimizer reported, and the 0.95 floor it held during background pruning, were
computed at phi29's reach whatever the polymerase or `--coverage-reach` said.
For bst that is 3x too generous (its realistic reach is 1000 bp) and for
equiphi29 it is short.

Two smaller things go with it.

`BackgroundAwareBaseOptimizer` -- the class the factory actually registers --
delegates to `HybridOptimizer` and then builds a second, unused
`BackgroundAwareOptimizer` "for backward compat". Nothing reads it, and
constructing it builds a `DominatingSetOptimizer` and a `NetworkOptimizer`, so
every background-aware run paid for two optimizers it never called.

And the module header advertises "Maintains excellent target coverage (>95%)".
That is a claim about outcomes, made on the strength of the floor above -- which
is measured on a third definition of coverage from the one reported. The floor
is a *setting*, and one the caller can lower.

NOTE ON SCOPE: this file exercises `BackgroundAwareOptimizer`, the standalone
three-stage class. No CLI run reaches it -- `--optimization-method
background-aware` resolves to `BackgroundAwareBaseOptimizer`, which delegates to
`HybridOptimizer` -- and nothing outside tests references it. So these tests
say the algorithm is right; they cannot say the algorithm the user gets is
right, and they passed throughout the period when the live wrapper was dropping
the configured Tm window, the extension reach and every selection weight.
`tests/test_background_aware_live_path.py` covers the object the factory
actually returns. Keep both: this one pins the algorithm, that one pins the
wiring.
"""

import numpy as np
import pytest

from neoswga.core.background_aware_optimizer import (
    BackgroundAwareBaseOptimizer,
    BackgroundAwareOptimizer,
)


class _Cache:
    def __init__(self, positions):
        self._positions = positions

    def get_positions(self, prefix, primer, strand="both"):
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.asarray(self._positions.get(primer, []), dtype=np.int64)


GENOME = 500_000
PRIMERS = {"AAAACCCCGGGG": [50_000, 250_000], "TTTTGGGGCCCC": [150_000, 350_000]}


def _optimizer(**kwargs):
    return BackgroundAwareOptimizer(
        position_cache=_Cache(PRIMERS),
        fg_prefixes=["fg"],
        bg_prefixes=[],
        fg_seq_lengths=[GENOME],
        bg_seq_lengths=[],
        **kwargs,
    )


@pytest.mark.parametrize("reach", [1000, 3000, 8000])
def test_coverage_uses_the_configured_reach(reach):
    """The regression: every reach gave the 3 kb answer."""
    optimizer = _optimizer(coverage_reach=reach)

    coverage = optimizer._calculate_coverage(list(PRIMERS))

    # Four sites, forward-only extension, none overlapping at these reaches.
    assert coverage == pytest.approx(4 * reach / GENOME, abs=0.01)


def test_the_reach_defaults_to_the_realistic_per_primer_value():
    """Unset must keep the previous behaviour rather than becoming zero."""
    from neoswga.core.coverage import resolve_coverage_reach

    assert _optimizer().coverage_reach == resolve_coverage_reach("phi29")


def test_the_stage_one_optimizer_gets_the_same_reach():
    """Selection and the reported figure have to agree; that is the whole
    reason `coverage_reach` is threaded rather than defaulted per site."""
    optimizer = _optimizer(coverage_reach=8000)

    assert optimizer.coverage_optimizer.extension_reach == 8000


def test_no_unused_second_optimizer_is_constructed():
    """Building it cost two optimizer constructions per run for nothing."""
    adapter = BackgroundAwareBaseOptimizer(
        position_cache=_Cache(PRIMERS),
        fg_prefixes=["fg"],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[],
        bg_seq_lengths=[],
    )

    assert not hasattr(adapter, "_direct_optimizer")


def test_the_module_does_not_advertise_a_coverage_outcome():
    """ ">95% coverage" is a claim about results, made on the strength of a
    configurable floor measured on a different definition from the reported
    one. Describing the floor is honest; promising the outcome is not."""
    import neoswga.core.background_aware_optimizer as module

    assert ">95%" not in (module.__doc__ or "")
