"""Stage 1 chooses the panel, and it was the stage choosing on the wrong thing.

Known Issue 16. `optimize_greedy` takes an `objective`, and supplying it makes
the set-cover select on occupancy-weighted coverage with a background tie-break
instead of on unweighted coverage bins. Commit `59a4ee3` added it under the
heading "The greedy now chooses on the quantity the design is judged on". No
production caller passed it: `hybrid_optimizer.py`, `dominating_set_adapter.py`
and `primer_expansion.py` all omitted it, and the only caller that supplied one
was a test.

It matters exactly where additives matter. Occupancy depends only on the
primer, so an unweighted bin count misranks two candidates by the ratio of
their occupancies, and across the pool the Tm gate admits that ratio is 1.8 on
phi29 at 30 C but 7.8 on equiphi29 at 42 C and 8.3 under DMSO 5% plus betaine
1 M.

**Why it was never wired, measured 2026-09-19 before wiring it.** The scan
recomputed the objective for every candidate at every pick. One
`compute_metrics` call on the Wolbachia design costs 36 ms, so a 2,000-
candidate pool at 12 picks is 14.4 minutes of objective evaluation and the
20,670-candidate inventory is about 2.5 hours -- against 39 s for the whole
`optimize` run as it stands. A correct rule nobody can afford to run is not a
fix.

So the objective arrives with a prescreen, which is the shape this codebase
already uses for the swap repair: rank every candidate by the cheap bin gain,
then score only the leaders with the full objective. The prescreen is not a new
criterion. It is exactly what the greedy used when it had no objective at all.
"""

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_objective import PoolConstraints, PoolObjective
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

    rng = random.Random(20260919)
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

    tmp = tmp_path_factory.mktemp("stage_one_objective")
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


def _factory(prefixes, primers, method="hybrid", **config_overrides):
    settings = dict(
        max_dimer_bp=3,
        max_self_dimer_bp=4,
        extension_reach=3_000,
        allow_dimer_relaxation=False,
        verbose=False,
    )
    settings.update(config_overrides)
    cache = PositionCache([prefixes["target"], prefixes["host"]], primers)
    return OptimizerFactory.create(
        method,
        cache,
        [prefixes["target"]],
        [GENOME],
        [prefixes["host"]],
        [GENOME],
        config=OptimizerConfig(**settings),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )


class TestTheObjectiveReachesStageOne:
    """The PATH, not the attribute. Both ends of this existed before and the
    path between them did not; that is the shape of the Phase 6 defect the
    `attach_search_config` docstring records."""

    def _spy(self, prefixes, primers, method, width=64):
        from neoswga.core import dominating_set_optimizer as ds

        seen = []
        original = ds.DominatingSetOptimizer.optimize_greedy

        def watched(self, *args, **kwargs):
            seen.append(kwargs.get("objective"))
            return original(self, *args, **kwargs)

        ds.DominatingSetOptimizer.optimize_greedy = watched
        try:
            optimizer = _factory(prefixes, primers, method, stage1_objective_width=width)
            optimizer.optimize(primers, target_size=6)
        finally:
            ds.DominatingSetOptimizer.optimize_greedy = original
        return seen

    def test_hybrid_stage_one_receives_one(self, designed):
        prefixes, primers = designed

        seen = self._spy(prefixes, primers, "hybrid")

        assert seen, "Stage 1 never ran"
        assert any(
            o is not None for o in seen
        ), "the greedy that picks the panel still receives objective=None"

    def test_background_aware_stage_one_receives_one(self, designed):
        prefixes, primers = designed

        seen = self._spy(prefixes, primers, "background-aware")

        assert seen and any(o is not None for o in seen)

    def test_dominating_set_receives_one(self, designed):
        prefixes, primers = designed

        seen = self._spy(prefixes, primers, "dominating-set")

        assert seen and any(o is not None for o in seen)

    def test_what_it_receives_can_measure_coverage(self, designed):
        """An objective whose `coverage` returns None is the same as no
        objective: `_objective_gain` treats it as zero gain for every
        candidate, so the scan falls back on its tie-break alone."""
        prefixes, primers = designed

        seen = [o for o in self._spy(prefixes, primers, "hybrid") if o is not None]

        assert seen[0].coverage(primers[:3]) is not None


class TestTheScanIsAffordable:
    """36 ms per evaluation is why a full recomputation was never viable."""

    def _optimizer(self, tmp_path, primers, width):
        prefix = str(tmp_path / "t")
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as db:
            for i, primer in enumerate(primers):
                db.create_dataset(primer, data=np.array([i * 900 + 100], dtype=np.int64))
            db.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))
        return DominatingSetOptimizer(
            PositionCache([prefix], primers),
            fg_prefixes=[prefix],
            fg_seq_lengths=[80_000],
            extension_reach=3_000,
            stage1_objective_width=width,
        )

    def _counting_objective(self, calls):
        def evaluate(panel):
            from types import SimpleNamespace

            calls.append(tuple(sorted(panel)))
            return SimpleNamespace(
                fg_coverage=0.01 * len(panel),
                effective_fg_coverage=0.01 * len(panel),
                selectivity_density=100.0,
                total_bg_sites=len(panel),
            )

        return PoolObjective(evaluate, PoolConstraints())

    def test_the_width_bounds_the_evaluations_per_pick(self, tmp_path):
        primers = [f"{'ACGT' * 3}"[:K] for _ in range(1)] + [
            "".join(c) for c in _distinct_kmers(60)
        ]
        calls = []
        optimizer = self._optimizer(tmp_path, primers, width=8)
        optimizer.optimize_greedy(
            primers, max_primers=4, verbose=False, objective=self._counting_objective(calls)
        )

        distinct = len(set(calls))
        assert distinct <= 8 * 4 + 4 + 1, (
            f"{distinct} distinct panels evaluated; the width was supposed to "
            f"bound this to about 8 per pick"
        )

    def test_an_unbounded_width_scores_everything(self, tmp_path):
        primers = ["".join(c) for c in _distinct_kmers(40)]
        calls = []
        optimizer = self._optimizer(tmp_path, primers, width=None)
        optimizer.optimize_greedy(
            primers, max_primers=2, verbose=False, objective=self._counting_objective(calls)
        )

        assert len(set(calls)) > 40, "None was supposed to restore the full scan"

    def test_the_width_is_configurable_from_the_optimizer_config(self, designed):
        prefixes, primers = designed

        optimizer = _factory(prefixes, primers, stage1_objective_width=16)

        inner = getattr(optimizer, "_hybrid", None)
        target = getattr(inner, "dominating_optimizer", None) if inner else None
        assert target is not None
        assert target.stage1_objective_width == 16


class TestStageOneDoesNotClobberStageTwo:
    """The regression this nearly shipped, and the second time this seam bit.

    Stage 2's swap refinement reads `pool_objective`, which `plan_pool`
    attaches -- that path is the Phase 6 fix recorded in
    `attach_search_config`. Writing Stage 1's objective into the SAME name
    overwrote it, with None on every default run, silently undoing that fix.
    Four tests in `test_the_objective_reaches_the_stage_that_refines.py`
    caught it. Stage 1 uses `stage1_pool_objective` for that reason.
    """

    def test_the_two_names_are_distinct(self, designed):
        from neoswga.core.swap_refinement import attach_search_config

        prefixes, primers = designed
        optimizer = _factory(prefixes, primers, stage1_objective_width=64)

        sentinel = object()
        attach_search_config(optimizer, "pool_objective", sentinel)
        optimizer.optimize(primers, target_size=6)

        assert (
            optimizer.pool_objective is sentinel
        ), "Stage 1 overwrote the objective Stage 2's refinement reads"

    def test_the_inner_optimizer_keeps_it_too(self, designed):
        """The inner object is the one that actually reads it."""
        from neoswga.core.swap_refinement import attach_search_config

        prefixes, primers = designed
        optimizer = _factory(prefixes, primers, stage1_objective_width=64)

        sentinel = object()
        attach_search_config(optimizer, "pool_objective", sentinel)
        optimizer.optimize(primers, target_size=6)

        assert optimizer._hybrid.pool_objective is sentinel


class TestItIsOffByDefault:
    """On measurement, not caution.

    Measured on the Wolbachia design at n=6/12/24, turning it on improves the
    metric it selects on -- effective coverage +0.0073, +0.0139, +0.0583 --
    and costs specificity every time: selectivity density -1.64, -6.50, -7.60,
    with host sites rising from 261 to 456 at n=24. Runtime is 3.5x to 8.4x.

    That is a trade, not an improvement, so the shipped default is unchanged
    and a user who wants coverage asks for it. The same resolution Known Issue
    11 reached for `DEFAULT_REDUNDANCY_THRESHOLD`.
    """

    def test_the_default_width_is_none(self):
        assert OptimizerConfig().stage1_objective_width is None

    def test_so_stage_one_gets_no_objective(self, designed):
        from neoswga.core.swap_refinement import stage1_objective

        prefixes, primers = designed

        assert stage1_objective(_factory(prefixes, primers)) is None

    def test_and_the_delivered_panel_is_unchanged(self, designed):
        """The guarantee that makes shipping this safe: with no width set, the
        greedy runs exactly the scan it ran before."""
        prefixes, primers = designed

        default = _factory(prefixes, primers).optimize(primers, target_size=6)
        explicit_off = _factory(prefixes, primers, stage1_objective_width=None).optimize(
            primers, target_size=6
        )

        assert tuple(default.primers) == tuple(explicit_off.primers)
        assert default.primers, "the baseline delivered nothing to compare"

    def test_a_width_changes_the_panel(self, designed):
        """The other half: if turning it on changed nothing, there would be
        nothing to measure and nothing to offer."""
        prefixes, primers = designed

        off = _factory(prefixes, primers).optimize(primers, target_size=6)
        on = _factory(prefixes, primers, stage1_objective_width=64).optimize(primers, target_size=6)

        assert off.primers and on.primers
        assert tuple(off.primers) != tuple(on.primers) or off.metrics is not None


class TestAnObjectiveThatCannotMeasureIsRefused:
    """The regression this nearly shipped.

    Without reaction conditions there is no temperature at which to evaluate
    occupancy, so `effective_fg_coverage` is None, `_objective_gain` reads zero
    for EVERY candidate, and the greedy has nothing to prefer. It selected an
    empty panel. An objective that cannot measure the accepted quantity is
    worse than no objective at all.
    """

    def test_no_conditions_means_no_stage_one_objective(self, designed):
        from neoswga.core.hybrid_optimizer import HybridBaseOptimizer
        from neoswga.core.swap_refinement import stage1_objective

        prefixes, primers = designed
        optimizer = HybridBaseOptimizer(
            position_cache=PositionCache([prefixes["target"]], primers),
            fg_prefixes=[prefixes["target"]],
            fg_seq_lengths=[GENOME],
            coverage_reach=3_000,
        )

        assert stage1_objective(optimizer) is None

    def test_and_the_panel_is_still_delivered(self, designed):
        from neoswga.core.hybrid_optimizer import HybridBaseOptimizer

        prefixes, primers = designed
        optimizer = HybridBaseOptimizer(
            position_cache=PositionCache([prefixes["target"]], primers),
            fg_prefixes=[prefixes["target"]],
            fg_seq_lengths=[GENOME],
            coverage_reach=3_000,
        )

        result = optimizer.optimize(primers, target_size=5)

        assert result.primers, "the greedy selected nothing at all"

    def test_with_conditions_and_a_width_one_is_built(self, designed):
        from neoswga.core.swap_refinement import stage1_objective

        prefixes, primers = designed
        optimizer = _factory(prefixes, primers, stage1_objective_width=64)

        assert stage1_objective(optimizer) is not None


class TestNothingChangesWithoutAnObjective:
    def test_the_panel_is_what_it_was(self, tmp_path):
        """Every existing caller passes no objective, and must be unaffected."""
        primers = ["".join(c) for c in _distinct_kmers(20)]
        prefix = str(tmp_path / "t")
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as db:
            for i, primer in enumerate(primers):
                db.create_dataset(primer, data=np.array([i * 900 + 100], dtype=np.int64))
            db.create_dataset("#record_starts", data=np.array([0], dtype=np.int64))
        optimizer = DominatingSetOptimizer(
            PositionCache([prefix], primers),
            fg_prefixes=[prefix],
            fg_seq_lengths=[80_000],
            extension_reach=3_000,
        )

        a = optimizer.optimize_greedy(primers, max_primers=5, verbose=False)["primers"]
        b = optimizer.optimize_greedy(primers, max_primers=5, verbose=False)["primers"]

        assert a == b


def _distinct_kmers(n):
    """n distinct k-mers, deterministic and not reverse complements of each
    other, so the HDF5 fixture has one dataset per primer."""
    import random

    rng = random.Random(7)
    out, seen = [], set()
    while len(out) < n:
        candidate = "".join(rng.choice("ACGT") for _ in range(K))
        if candidate in seen or reverse_complement(candidate) in seen:
            continue
        seen.add(candidate)
        out.append(candidate)
    return out
