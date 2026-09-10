"""What the REGISTERED background-aware optimizer is wired to.

Three tests in this file exercised `BackgroundAwareOptimizer`, the standalone
three-stage class that no dispatch path reached. That class was deleted on
2026-09-10 and its cases went two ways: the reach assertions duplicated
`tests/test_background_aware_live_path.py`, which pins the same property on the
object the factory returns, and the strand-extension cases from its sibling file
were ported to `tests/test_strand_aware_extension_reach.py`, against the
implementation the shipping optimizers use.

What remains are the two assertions about the live object. Both were written
here because the deleted class was the thing being compared against.
"""

import numpy as np

from neoswga.core.background_aware_optimizer import BackgroundAwareBaseOptimizer


class _Cache:
    def __init__(self, positions):
        self._positions = positions

    def get_positions(self, prefix, primer, strand="both"):
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.asarray(self._positions.get(primer, []), dtype=np.int64)


GENOME = 500_000
PRIMERS = {"AAAACCCCGGGG": [50_000, 250_000], "TTTTGGGGCCCC": [150_000, 350_000]}


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
