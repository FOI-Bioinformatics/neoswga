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


# Optimizers confirmed to register at import time on every supported
# platform. Optional-dep-gated entries (moea, milp when the solver is
# missing) are handled by try/except in the test.
FACTORY_NAMES_TO_TEST = [
    "greedy", "dominating-set", "weighted-set-cover",
    "hybrid", "network", "background-aware", "equiphi29",
    "tiling", "clique", "bg-prefilter", "normalized",
]


@pytest.mark.parametrize("name", FACTORY_NAMES_TO_TEST)
def test_optimizer_factory_create_accepts_conditions(name, conditions):
    """Every registered optimizer must build via the factory without a
    NameError / TypeError and must attach conditions to self.conditions."""
    try:
        optimizer = OptimizerFactory.create(
            name=name,
            position_cache=_DummyCache(),
            fg_prefixes=["fg"],
            fg_seq_lengths=[1000],
            bg_prefixes=[],
            bg_seq_lengths=[],
            conditions=conditions,
        )
    except ImportError:
        pytest.skip(f"{name} requires an optional dependency not installed")

    assert optimizer.conditions is conditions, (
        f"OptimizerFactory.create('{name}', conditions=...) did not attach "
        f"the ReactionConditions to the optimizer. Got: {optimizer.conditions!r}"
    )


@pytest.mark.parametrize("name", ["moea", "milp"])
def test_optional_optimizers_tolerate_missing_deps(name, conditions):
    """MOEA (pymoo) and MILP (mip) are optional; creation should either
    succeed with conditions or raise ImportError, not NameError."""
    try:
        optimizer = OptimizerFactory.create(
            name=name,
            position_cache=_DummyCache(),
            fg_prefixes=["fg"],
            fg_seq_lengths=[1000],
            bg_prefixes=[],
            bg_seq_lengths=[],
            conditions=conditions,
        )
    except ImportError:
        pytest.skip(f"{name}'s optional dep not installed")
    except Exception as e:
        pytest.fail(f"{name} raised unexpected {type(e).__name__}: {e}")
    else:
        assert optimizer.conditions is conditions
