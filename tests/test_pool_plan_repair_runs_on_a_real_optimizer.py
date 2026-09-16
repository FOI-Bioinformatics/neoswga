"""The repair path, through the real objects rather than stubs.

`tests/test_pool_plan_repair.py` drives `plan_pool` with a stub optimizer whose
metrics come from a table. That is the right way to pin the decision rules, and
it never touches `OptimizerFactory`, `PositionCache`, `LazyDimerCompatibility`,
`refine_by_swaps` with a real objective, or the beam over a real candidate pool.
An integration break in any of those would pass every test in that file.

So this runs `plan_pool` the way the command line does, over a real HDF5 index
written by the real scanner, with a selectivity floor tight enough that the
optimizer's first answer does not qualify.

The assertions are deliberately not "the repair succeeded". Whether a particular
synthetic pool contains a qualifying panel is a property of the fixture, and a
test asserting it would be pinning the fixture rather than the code. What is
asserted is that the machinery runs, that it reports what it did, and that every
number in a row describes the panel that row actually delivered.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_metrics import compute_pool_metrics
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import reverse_complement
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

GENOME_LENGTH = 80_000
K = 10
# Tight enough that the optimizer's first panel misses it on this fixture, so
# the repair path is entered rather than skipped.
DENSITY_FLOOR = 40.0


@pytest.fixture(scope="module")
def designed(tmp_path_factory):
    """A target and a host, both carrying real sites for every candidate.

    The host gets fewer sites per candidate than the target, and a handful of
    candidates get many host sites, so selectivity density actually varies
    across panels. A fixture where every candidate is equally selective cannot
    exercise a selectivity constraint at all.
    """
    rng = random.Random(20260916)
    target = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    host = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers, seen = [], set()
    while len(primers) < 30:
        candidate = "".join(rng.choice("ACGT") for _ in range(K))
        if candidate in seen or reverse_complement(candidate) in seen:
            continue
        seen.add(candidate)
        primers.append(candidate)
        for _ in range(rng.randint(4, 14)):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            target[pos : pos + K] = list(candidate)
        # Every third candidate is a host binder, so the panels differ in
        # selectivity and not only in coverage.
        host_sites = rng.randint(6, 20) if len(primers) % 3 == 0 else rng.randint(0, 2)
        for _ in range(host_sites):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            host[pos : pos + K] = list(candidate)

    tmp = tmp_path_factory.mktemp("plan_pool_repair")
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

    _ensure_optimizers_registered()
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
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    plan = plan_pool(
        optimizer,
        primers,
        [4, 8],
        [0.5],
        primer_length=K,
        min_selectivity_density=DENSITY_FLOOR,
    )
    return optimizer, plan


def test_the_plan_completes_through_the_real_objects(designed):
    _, plan = designed

    assert plan["rows"], "plan_pool returned no rows"
    assert plan["background_assessed"] is True


def test_every_row_reports_whether_it_was_repaired(designed):
    """Silence would leave a reader unable to tell a first answer from a second."""
    _, plan = designed

    for row in plan["rows"]:
        if row["status"] != "evaluated":
            continue
        assert "repair" in row, f"row {row['requested_size']} carries no repair record"
        assert isinstance(row["repair"]["attempted"], bool)


def test_the_repair_path_is_actually_entered(designed):
    """Guard the guard: a fixture every panel passes would prove nothing."""
    _, plan = designed
    attempted = [r for r in plan["rows"] if (r.get("repair") or {}).get("attempted")]

    assert attempted, (
        "no row missed the selectivity floor, so this fixture never reached the "
        "repair path; raise DENSITY_FLOOR until one does"
    )


def test_a_row_describes_the_panel_it_delivered(designed):
    """Recomputed from the primers in the row, not read back from the search."""
    optimizer, plan = designed

    for row in plan["rows"]:
        if row["status"] != "evaluated":
            continue
        fresh = compute_pool_metrics(optimizer, row["primers"])
        assert row["raw_coverage"] == pytest.approx(fresh.fg_coverage)
        assert row["effective_coverage"] == pytest.approx(fresh.effective_fg_coverage)
        assert row["selectivity_density"] == pytest.approx(fresh.selectivity_density)
        assert row["background_sites"] == fresh.total_bg_sites


def test_no_delivered_panel_breaks_the_dimer_limit(designed):
    """The repair may swap and rebuild; it may not relax this."""
    _, plan = designed

    for row in plan["rows"]:
        assert row.get("violating_pairs", 0) == 0, row["primers"]
        assert row.get("self_dimers", 0) == 0, row["primers"]


def test_an_eligible_row_really_clears_the_floor(designed):
    _, plan = designed

    for row in plan["rows"]:
        if row.get("eligible"):
            assert row["selectivity_density"] >= DENSITY_FLOOR
            assert not row["failed_constraints"]


def test_a_repaired_panel_is_no_worse_than_the_alternative_being_ineligible(designed):
    """A failed repair must leave a panel, not nothing."""
    _, plan = designed

    for row in plan["rows"]:
        if (row.get("repair") or {}).get("attempted") and row["status"] == "evaluated":
            assert row["primers"], "the repair returned an empty panel"
            assert len(row["primers"]) <= row["requested_size"]
