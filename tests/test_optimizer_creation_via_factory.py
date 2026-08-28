"""Regression test for the Critical bug: every registered optimizer must
construct cleanly through OptimizerFactory.create and must attach the
passed ReactionConditions to its `.conditions` attribute.

Before the fix, NetworkBaseOptimizer raised NameError for `conditions`
because the identifier was not declared in its __init__ signature. The
unit test that exercised NetworkOptimizer directly passed the factory
layer entirely. This test goes through the factory so the full plumbing
chain is exercised for every registered optimizer.
"""

import pytest

from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core import unified_optimizer as _uo


class _DummyCache:
    def get_positions(self, *a, **kw):
        return []


@pytest.fixture(scope="module")
def conditions():
    return ReactionConditions(temp=30.0, polymerase="phi29", betaine_m=1.0)


@pytest.fixture(scope="module", autouse=True)
def _register_optimizers():
    _uo._ensure_optimizers_registered()


# Supported optimizers — kept after the optimizer-zoo trim.
FACTORY_NAMES_TO_TEST = [
    "dominating-set",
    "hybrid",
    "network",
    "background-aware",
]


@pytest.mark.parametrize("name", FACTORY_NAMES_TO_TEST)
def test_optimizer_factory_create_accepts_conditions(name, conditions):
    """Every registered optimizer must build via the factory without a
    NameError / TypeError and must attach conditions to self.conditions."""
    optimizer = OptimizerFactory.create(
        name=name,
        position_cache=_DummyCache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        conditions=conditions,
    )

    assert optimizer.conditions is conditions, (
        f"OptimizerFactory.create('{name}', conditions=...) did not attach "
        f"the ReactionConditions to the optimizer. Got: {optimizer.conditions!r}"
    )
