"""The position check has to run on the real path, or it is not a check.

Phase 4 increment 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`.

`ensure_positions` existed, was tested, and never ran: nothing in production
attached a cache, and it returned quietly when none was attached. That is the
shape of the whole audit, so the wiring is pinned separately from the behaviour.

`plan_pool` holds both halves. It has the source that hands out the frontier and
the optimizer that holds the cache the frontier will be scored on, so it is the
place where the two are introduced to each other.
"""

import h5py
import numpy as np
import pytest

from neoswga.core.candidate_source import ListCandidateSource
from neoswga.core.position_cache import MissingPositionsError, PositionCache

TARGET = 6000
FOREGROUND = {
    "GCTAAAGACAAT": [100, 1500, 3000, 4500],
    "TACATAACATAC": [700, 2200, 3700, 5200],
    "ACGTCAGCACGA": [400, 2900],
}
# Deliberately short of the foreground: the third candidate has no host entry.
BACKGROUND = {
    "GCTAAAGACAAT": [11000],
    "TACATAACATAC": [],
}


def _write_index(path, entries):
    with h5py.File(path, "w") as handle:
        for primer, positions in entries.items():
            handle.create_dataset(primer, data=np.array(positions, dtype=np.int64))
        handle.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))


@pytest.fixture
def design(tmp_path):
    """A cache, an optimizer over it, and the prefixes both were built on."""
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.reaction_conditions import ReactionConditions
    from neoswga.core.unified_optimizer import _ensure_optimizers_registered

    _write_index(tmp_path / "fg_12mer_positions.h5", FOREGROUND)
    _write_index(tmp_path / "bg_12mer_positions.h5", BACKGROUND)
    fg, bg = str(tmp_path / "fg"), str(tmp_path / "bg")

    _ensure_optimizers_registered()
    cache = PositionCache([fg, bg], list(FOREGROUND), on_missing="warn")
    optimizer = OptimizerFactory.create(
        "hybrid",
        cache,
        [fg],
        [TARGET],
        [bg],
        [20000],
        config=OptimizerConfig(
            max_dimer_bp=3,
            max_self_dimer_bp=4,
            extension_reach=1000,
            refinement_method="swap",
            allow_dimer_relaxation=False,
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    return optimizer, cache


def _plan(optimizer, source):
    from neoswga.core.pool_planner import plan_pool

    return plan_pool(
        optimizer,
        source,
        [2],
        [0.5],
        primer_length=12,
        max_background_sites=1000,
        repair=False,
    )


def test_plan_pool_attaches_the_cache_to_the_source(design):
    """The source cannot vouch for a batch it has no cache for."""
    optimizer, cache = design
    source = ListCandidateSource(["GCTAAAGACAAT", "TACATAACATAC"])

    _plan(optimizer, source)

    assert getattr(source, "position_cache", None) is cache


def test_plan_pool_refuses_a_candidate_with_no_entry_on_the_background(design):
    """Unknown specificity must not be designed with as though it were zero."""
    optimizer, cache = design
    source = ListCandidateSource(["GCTAAAGACAAT", "TACATAACATAC", "ACGTCAGCACGA"])

    with pytest.raises(MissingPositionsError) as raised:
        _plan(optimizer, source)

    assert "ACGTCAGCACGA" in str(raised.value)


def test_plan_pool_accepts_a_frontier_that_is_indexed_throughout(design):
    """A measured zero on the host is a good candidate, not a missing one."""
    optimizer, cache = design
    source = ListCandidateSource(["GCTAAAGACAAT", "TACATAACATAC"])

    plan = _plan(optimizer, source)

    assert plan["candidate_source"]["examined"] == 2
