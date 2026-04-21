"""Standalone tests for the MILP optimizer.

Gated on the optional `mip` dependency; skipped when it is not installed.
Closes the TODO at `milp_optimizer.py:447` (Phase 15G).
"""

import pytest


def _mip_available() -> bool:
    try:
        import mip  # noqa: F401
        return True
    except ImportError:
        return False


pytestmark = pytest.mark.skipif(
    not _mip_available(),
    reason="python-mip not installed; install with `pip install mip`",
)


def test_milp_factory_accepts_conditions():
    """The factory wrapper must accept `conditions` without NameError; this
    is the Phase 7A / code-review regression (the network wrapper had this
    bug; milp is its sibling)."""
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.reaction_conditions import ReactionConditions
    from neoswga.core import unified_optimizer as _uo
    _uo._ensure_optimizers_registered()

    class DummyCache:
        def get_positions(self, *a, **kw):
            return []

    optimizer = OptimizerFactory.create(
        name="milp",
        position_cache=DummyCache(),
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        conditions=ReactionConditions(temp=30.0, polymerase="phi29"),
    )
    assert optimizer.conditions is not None


def test_milp_reports_coverage_only_in_doctor():
    """Phase 15E: MILP is coverage-only by design. ADDITIVE_AWARE = False."""
    from neoswga.core.optimizer_factory import OptimizerRegistry
    from neoswga.core import unified_optimizer as _uo
    _uo._ensure_optimizers_registered()
    cls = OptimizerRegistry.get("milp")
    assert getattr(cls, "ADDITIVE_AWARE", False) is False, (
        "MILP is coverage-only by design; ADDITIVE_AWARE must be False"
    )
