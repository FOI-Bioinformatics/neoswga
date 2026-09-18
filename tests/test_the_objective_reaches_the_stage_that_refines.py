"""The objective has to reach the object that runs the refinement.

`tests/test_stage_two_refines_on_the_shared_objective.py` pins that `plan_pool`
assigns `pool_objective` and that `refine_hybrid_stage2` reads it. Both are true
and the two were never connected.

`plan_pool` assigns it to the optimizer it was handed, which for every
command-line path is a wrapper: `OptimizerFactory` returns
`HybridBaseOptimizer` or `BackgroundAwareBaseOptimizer`, and each delegates the
search to an inner `HybridOptimizer`. `_swap_refine` is a method of the INNER
one, so `getattr(optimizer, "pool_objective", None)` inside
`refine_hybrid_stage2` read the inner object's attribute, which nobody set.
Measured through a real `plan-pool` design: the refinement ran once and received
None.

So Stage 2 refined on raw covered bases while the row was accepted on
occupancy-weighted coverage under a specificity floor, which is the two-rule
split Phase 2 set out to remove.

The existing tests could not catch it because one checks the assigning end by
AST and the other checks the reading end by source text. Neither asserts that
the object assigned to is the object read from. This file asserts the path.
"""

import json

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import reverse_complement
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

GENOME = 60_000
K = 10


@pytest.fixture(scope="module")
def designed(tmp_path_factory):
    """A target and a host with real sites, so a real optimizer can run."""
    import random

    rng = random.Random(20260918)
    target = list("".join(rng.choice("ACGT") for _ in range(GENOME)))
    host = list("".join(rng.choice("ACGT") for _ in range(GENOME)))

    primers, seen = [], set()
    while len(primers) < 24:
        candidate = "".join(rng.choice("ACGT") for _ in range(K))
        if candidate in seen or reverse_complement(candidate) in seen:
            continue
        seen.add(candidate)
        primers.append(candidate)
        for _ in range(rng.randint(4, 12)):
            pos = rng.randrange(0, GENOME - K)
            target[pos : pos + K] = list(candidate)
        for _ in range(rng.randint(6, 18) if len(primers) % 3 == 0 else rng.randint(0, 2)):
            pos = rng.randrange(0, GENOME - K)
            host[pos : pos + K] = list(candidate)

    tmp = tmp_path_factory.mktemp("objective_reaches")
    prefixes = {}
    for name, sequence in (("target", "".join(target)), ("host", "".join(host))):
        prefix = str(tmp / name)
        prefixes[name] = prefix
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as handle:
            for candidate in primers:
                for key in {candidate, reverse_complement(candidate)}:
                    hits, i = [], sequence.find(key)
                    while i != -1:
                        hits.append(i)
                        i = sequence.find(key, i + 1)
                    handle.create_dataset(key, data=np.array(hits, dtype=np.int64))
            handle.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))

    _ensure_optimizers_registered()
    return prefixes, primers


def _run(prefixes, primers, method, spy):
    """One design through the real factory, watching what the refinement gets."""
    from neoswga.core import hybrid_optimizer as ho

    cache = PositionCache([prefixes["target"], prefixes["host"]], primers)
    optimizer = OptimizerFactory.create(
        method,
        cache,
        [prefixes["target"]],
        [GENOME],
        [prefixes["host"]],
        [GENOME],
        config=OptimizerConfig(
            max_dimer_bp=3,
            max_self_dimer_bp=4,
            extension_reach=3_000,
            refinement_method="swap",
            allow_dimer_relaxation=False,
            max_frontier_refills=0,
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    original = ho.HybridOptimizer._swap_refine

    def watched(self, panel, candidates, fixed):
        spy.append(getattr(self, "pool_objective", None))
        return original(self, panel, candidates, fixed)

    ho.HybridOptimizer._swap_refine = watched
    try:
        plan_pool(
            optimizer,
            primers,
            [8],
            [0.4],
            primer_length=K,
            coverage_metric="effective",
            min_selectivity_density=20.0,
        )
    finally:
        ho.HybridOptimizer._swap_refine = original
    return optimizer


@pytest.mark.parametrize("method", ["hybrid", "background-aware"])
def test_the_refinement_receives_the_objective(designed, method):
    """The path, not the two ends.

    `background-aware` matters as much as `hybrid`: it is what `plan-pool`
    chooses by default whenever a background is supplied, and it wraps the same
    inner optimizer.
    """
    prefixes, primers = designed
    spy = []

    _run(prefixes, primers, method, spy)

    assert spy, "stage 2 refinement never ran, so this test proves nothing"
    assert all(seen is not None for seen in spy), (
        f"stage 2 refinement ran {len(spy)} time(s) under {method} and received "
        f"no objective, so it refined on raw covered bases while the row was "
        f"accepted on occupancy-weighted coverage under a specificity floor"
    )


@pytest.mark.parametrize("method", ["hybrid", "background-aware"])
def test_the_inner_optimizer_carries_the_objective_too(designed, method):
    """Asserted on the object, not on the source text.

    The wrapper is what `plan_pool` is handed; the inner one is what searches.
    Both must carry it, since either could be the one a future caller reads.
    """
    prefixes, primers = designed

    optimizer = _run(prefixes, primers, method, [])

    assert getattr(optimizer, "pool_objective", None) is not None
    inner = getattr(optimizer, "_hybrid", None)
    assert inner is not None, "the wrapper stopped delegating; update this test"
    assert getattr(inner, "pool_objective", None) is not None, (
        "the inner optimizer has no objective, so the attribute was attached to "
        "the wrapper and read from the delegate"
    )


def test_the_objective_the_refinement_gets_is_the_one_the_row_is_accepted_on(designed):
    """Not merely non-None: the same object, or the two-rule split returns."""
    prefixes, primers = designed
    spy = []

    optimizer = _run(prefixes, primers, "hybrid", spy)

    assert spy[0] is optimizer.pool_objective
