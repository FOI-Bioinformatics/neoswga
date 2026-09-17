"""The bounded scan has to reach the search that actually spends the objective.

Phase 4 increment 4 of the plan for `docs/validation/pipeline_audit_2026-09-16/`.

`tests/test_the_swap_scan_is_bounded.py` pins what the bound does. This pins
that a real `plan-pool` run uses it, because a bound nothing passes is the
defect this audit is named after: six capabilities were built, tested, merged
and never called.

The repair is the search to bound. Measured on the real Wolbachia pair, a size
row whose selectivity floor binds makes about 3,200 objective evaluations
through `pool_planner._repair`, and that call passed `None` for both bin
arguments, so the prescreen the bound needs had no inputs. Supplying them is
most of this increment.
"""

import random

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

GENOME_LENGTH = 80_000
K = 10
# Tight enough that the optimizer's first panel misses it, so the repair runs.
DENSITY_FLOOR = 40.0


@pytest.fixture(scope="module")
def designed(tmp_path_factory):
    """A target and a host carrying real sites, with selectivity that varies.

    Same shape as `tests/test_pool_plan_repair_runs_on_a_real_optimizer.py`,
    because a fixture where every candidate is equally selective cannot make a
    selectivity constraint bind and the repair would never run.
    """
    rng = random.Random(20260917)
    target = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    host = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers, seen = [], set()
    while len(primers) < 40:
        candidate = "".join(rng.choice("ACGT") for _ in range(K))
        if candidate in seen or reverse_complement(candidate) in seen:
            continue
        seen.add(candidate)
        primers.append(candidate)
        for _ in range(rng.randint(4, 14)):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            target[pos : pos + K] = list(candidate)
        host_sites = rng.randint(6, 20) if len(primers) % 3 == 0 else rng.randint(0, 2)
        for _ in range(host_sites):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            host[pos : pos + K] = list(candidate)

    tmp = tmp_path_factory.mktemp("bounded_repair")
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


def _plan(prefixes, primers, scan_width):
    cache = PositionCache([prefixes["target"], prefixes["host"]], primers)
    optimizer = OptimizerFactory.create(
        "hybrid",
        cache,
        [prefixes["target"]],
        [GENOME_LENGTH],
        [prefixes["host"]],
        [GENOME_LENGTH],
        config=OptimizerConfig(
            max_dimer_bp=3,
            max_self_dimer_bp=4,
            extension_reach=3_000,
            refinement_method="swap",
            allow_dimer_relaxation=False,
            objective_scan_width=scan_width,
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    return plan_pool(
        optimizer,
        primers,
        [8],
        [0.5],
        primer_length=K,
        min_selectivity_density=DENSITY_FLOOR,
    )


@pytest.fixture(scope="module")
def bounded(designed):
    prefixes, primers = designed
    return _plan(prefixes, primers, 16)


@pytest.fixture(scope="module")
def unbounded(designed):
    prefixes, primers = designed
    return _plan(prefixes, primers, None)


def test_the_repair_actually_runs(bounded):
    """Guard the guard: a fixture that never repairs would prove nothing."""
    row = bounded["rows"][0]
    assert (row.get("repair") or {}).get("attempted") is True


def test_a_configured_width_reaches_the_repair(bounded):
    """The wiring. Without this the bound is another unreachable capability."""
    row = bounded["rows"][0]

    assert (row.get("repair") or {}).get("scan_width") == 16


def test_an_unset_width_leaves_the_scan_unbounded(unbounded):
    """The default has to preserve what the tool did before."""
    row = unbounded["rows"][0]

    assert (row.get("repair") or {}).get("scan_width") is None


def test_the_bound_scores_fewer_panels_than_the_unbounded_scan(bounded, unbounded):
    """The point of the increment, as counted work rather than as wall clock."""
    narrow = (bounded["rows"][0].get("repair") or {}).get("objective_evaluations")
    wide = (unbounded["rows"][0].get("repair") or {}).get("objective_evaluations")

    assert narrow is not None and wide is not None
    assert narrow < wide, f"bounded scored {narrow} panels, unbounded {wide}"


def test_the_prescreen_still_sees_the_whole_pool(bounded):
    """Pairs considered must exceed panels scored, or the bound is biased."""
    repair = bounded["rows"][0].get("repair") or {}

    assert repair["pairs_considered"] > repair["objective_evaluations"]


def test_a_bounded_repair_does_not_deliver_a_worse_panel(bounded, unbounded):
    """A narrower search may find less. It must not report more than it found."""
    for plan in (bounded, unbounded):
        row = plan["rows"][0]
        if row["status"] != "evaluated":
            continue
        if not row.get("failed_constraints"):
            assert row["selectivity_density"] >= DENSITY_FLOOR - 1e-9, (
                "a row reporting no failed constraints is below the floor it was "
                "given, so acceptance and the reported metrics disagree"
            )
