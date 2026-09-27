"""The counted allowance is never exceeded by a real search, on real stages.

`tests/test_search_budget_contract.py` pins the ledger's own arithmetic and
`tests/test_the_allowance_holds_against_a_fake_clock.py` pins its limits. Neither
runs a search. Task 6 asks for the production-path half: proposal selection,
repair, swaps, deletion, refills and alternatives, under one declared allowance.

That half matters because every defect this plan found in the ledger was a PATH
defect rather than an arithmetic one. The clique loop evaluated outside it. The
alternatives ran with it unbound. The reduction stage threw away work when it
raised. In each case the ledger was correct and something did not go through it,
which is invisible to a test that never starts a search.

Built on a synthetic genome written straight to HDF5, following
`test_hybrid_optimizer_run.py`, so it needs no external counter and no prepared
example directory.
"""

import random

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.optimization_service import OptimizationRequest, run_panel_search
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.panel_refinement import objective_for_optimizer
from neoswga.core.pool_objective import PoolConstraints
from neoswga.core.position_cache import PositionCache
from neoswga.core.search_control import SearchBudget, collect_alternative_sets
from neoswga.core.thermodynamics import reverse_complement
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

GENOME_LENGTH = 30_000


@pytest.fixture(scope="module")
def genome():
    """Twenty primers with differing coverage, so selection has real choices."""
    rng = random.Random(4242)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    primers: list[str] = []
    for i in range(20):
        primer = "".join(rng.choice("ACGT") for _ in range(10))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(2 + (i % 4)):
            pos = (i * 1_400 + j * 600) % (GENOME_LENGTH - 20)
            seq[pos : pos + 10] = list(primer)
    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture
def built(tmp_path, genome):
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as handle:
        for primer in genome["primers"]:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], genome["seq"].find(key)
                while i != -1:
                    positions.append(i)
                    i = genome["seq"].find(key, i + 1)
                if positions:
                    handle.create_dataset(key, data=np.array(positions, dtype=np.int64))
    _ensure_optimizers_registered()
    return PositionCache([prefix], genome["primers"]), prefix


def _optimizer(built, method):
    cache, prefix = built
    return OptimizerFactory.create(
        name=method,
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        config=None,
        coverage_reach=3_000,
    )


def _search(built, genome, method, allowance, **request_kwargs):
    optimizer = _optimizer(built, method)
    objective_for_optimizer(optimizer, request_kwargs.pop("constraints", None))
    budget = SearchBudget(max_evaluations=allowance)
    result = run_panel_search(
        OptimizationRequest(
            optimizer,
            tuple(genome["primers"]),
            request_kwargs.pop("target_size", 6),
            budget=budget,
            **request_kwargs,
        )
    )
    return result, budget, optimizer


@pytest.mark.parametrize("method", ["dominating-set", "hybrid"])
@pytest.mark.parametrize("allowance", [1, 5, 40])
def test_a_real_search_never_spends_past_its_allowance(built, genome, method, allowance):
    _result, budget, _optimizer = _search(built, genome, method, allowance)

    assert budget.evaluations <= allowance, (
        f"{method} spent {budget.evaluations} of an allowance of {allowance}; a "
        f"declared total that the search can exceed is not a limit"
    )


@pytest.mark.parametrize("method", ["dominating-set", "hybrid"])
def test_a_real_search_returns_a_panel_even_on_a_tiny_allowance(built, genome, method):
    """Exhaustion is a stopping point, not a failure: something usable comes back."""
    result, _budget, _optimizer = _search(built, genome, method, allowance=1)

    assert result is not None
    assert set(result.primers) <= set(genome["primers"]), "a primer was invented"


def test_deletion_under_an_allowance_stays_within_it(built, genome):
    """The `minimize` path, which is where an accepted removal used to be lost."""
    _result, budget, _optimizer = _search(
        built,
        genome,
        "dominating-set",
        allowance=25,
        minimize=True,
        target_coverage=0.1,
        constraints=PoolConstraints(),
    )

    assert budget.evaluations <= 25


def test_alternatives_never_push_the_ledger_past_the_allowance(built, genome):
    """The ledger is bound across alternatives and cannot be exceeded there."""
    result, budget, optimizer = _search(built, genome, "dominating-set", allowance=2000)
    spent_by_primary = budget.evaluations

    collect_alternative_sets(
        primary=result,
        optimizer=optimizer,
        candidates=genome["primers"],
        target_size=4,
        max_sets=5,
        max_iterations=5,
        budget=budget,
    )

    assert budget.evaluations >= spent_by_primary, "the ledger went backwards"
    assert budget.evaluations <= 2000


@pytest.mark.parametrize("method", ["dominating-set", "hybrid"])
def test_an_alternative_search_consults_no_objective_so_the_bound_is_inert(built, genome, method):
    """Measured, and recorded because it contradicts a claim I made.

    `collect_alternative_sets` was given the run's ledger on 2026-09-27, and the
    commit said alternatives then "spend the run's allowance instead of none".
    The BINDING is real; the spending is not. An alternative search calls
    `optimizer.optimize` and nothing else, and neither shipped method evaluates
    the shared objective inside `optimize`: measured at 0 objective evaluations
    for both, with ample allowance.

    So binding the ledger closed a hole without bounding anything in production
    yet, and that function's `except SearchBudgetExhausted` clause remains
    unreachable on these methods. This is the same shape one layer deeper --
    both ends exist and the path does not -- and it is what routing alternatives
    through `run_panel_search`, the remaining Task 6 item, would close.

    This test exists to fail when that happens, so the claim gets corrected in
    the same change that earns it.
    """
    from neoswga.core.search_control import budgeted_objective

    optimizer = _optimizer(built, method)
    objective = objective_for_optimizer(optimizer, None)
    budget = SearchBudget(max_evaluations=5_000)

    with budgeted_objective(objective, budget):
        optimizer.optimize(genome["primers"], 5)

    assert budget.evaluations == 0, (
        f"{method}.optimize now evaluates the shared objective. That is an "
        f"improvement: update this test and the alternatives claim together, "
        f"because the ledger now bounds an alternative search for real"
    )


def test_every_stage_records_a_reason_and_what_it_had_spent(built, genome):
    """A termination without a reason cannot be told from one nobody recorded."""
    result, _budget, _optimizer = _search(built, genome, "dominating-set", allowance=30)

    stages = list(getattr(result, "stage_history", ()) or [])
    assert stages, "the search recorded no stages"
    for stage in stages:
        assert "stage" in stage
        spent = stage.get("search_budget")
        assert isinstance(spent, dict), f"{stage['stage']} recorded no ledger state"
        assert spent["evaluations"] <= 30
        assert "stop_reason" in spent, f"{stage['stage']} recorded no stop reason field"
        assert spent["evaluation_scope"] == "uncached_shared_objective"


def test_an_allowance_of_zero_stops_before_any_evaluation(built, genome):
    _result, budget, _optimizer = _search(built, genome, "dominating-set", allowance=0)

    assert budget.evaluations == 0
    assert budget.describe()["stop_reason"] == "total_evaluation_budget"
