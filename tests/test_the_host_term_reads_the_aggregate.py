"""`background-aware`'s host term reads `background_aggregate` (Phase 7c).

Known Issue 14 is the reason this file is separate from
`tests/test_per_host_figures.py`, which asserts on the METRICS: *reading the
background is not acting on it*. The issue was found twice in this optimizer --
`_network_refine` carried the host term and `_swap_refine` read the objective,
and neither carried both -- so "the term reads the key" has to be shown on the
stage whose output is delivered, not on a stage that runs and is then
overwritten.

So the order here is deliberate:

1. establish WHICH stage produces the delivered panel, by instrumenting the two
   candidates and comparing their recorded outputs against the result;
2. then show that the charge that stage ranks on moves with the key;
3. then show the delivered panel itself moves.

A test that only did (2) would pass if the deciding stage ignored the term, and
a test that only did (3) would not say what made it move.
"""

from __future__ import annotations

import random

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("networkx")

import numpy as np

from neoswga.core import parameter
from neoswga.core import unified_optimizer as uo
from neoswga.core.background_aggregate import SUM, WORST_CASE
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

FG_LENGTH = 30_000
SMALL_HOST = 30_000
LARGE_HOST = 1_290_000


@pytest.fixture(scope="module", autouse=True)
def _registered():
    uo._ensure_optimizers_registered()


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Twelve candidates over a 30 kb target, a 30 kb host and a 1.29 Mb host.

    The host loads are the point. Six candidates bind only the SMALL host and
    six only the LARGE one, with the large-host binders carrying many more
    sites. Pooled, the large-host binders look like the dirty ones; per host,
    the small-host binders are far denser and are the ones a worst-case design
    should shed. The two readings therefore disagree about which candidate is
    worst, which is what makes a change in the delivered panel attributable.
    """
    rng = random.Random(4242)
    directory = tmp_path_factory.mktemp("host_term_world")
    sequence = list("".join(rng.choice("ACGT") for _ in range(FG_LENGTH)))

    primers: list[str] = []
    while len(primers) < 12:
        candidate = "".join(rng.choice("ACGT") for _ in range(10))
        if candidate not in primers:
            primers.append(candidate)
    for index, primer in enumerate(primers):
        for j in range(2 + (index % 3)):
            start = (index * 2_400 + j * 700) % (FG_LENGTH - 20)
            sequence[start : start + 10] = list(primer)
    target = "".join(sequence)

    fg = str(directory / "fg")
    with h5py.File(f"{fg}_10mer_positions.h5", "w") as handle:
        for primer in primers:
            for key in sorted({primer, reverse_complement(primer)}):
                hits, index = [], target.find(key)
                while index != -1:
                    hits.append(index)
                    index = target.find(key, index + 1)
                handle.create_dataset(key, data=np.array(hits, dtype=np.int64))

    small, large = str(directory / "small"), str(directory / "large")
    with h5py.File(f"{small}_10mer_positions.h5", "w") as low, h5py.File(
        f"{large}_10mer_positions.h5", "w"
    ) as high:
        for index, primer in enumerate(primers):
            # Six on the small host, six on the large one.
            on_small = index % 2 == 0
            small_hits = [i * 900 + 17 for i in range(3)] if on_small else []
            large_hits = [] if on_small else [i * 30_000 + 23 for i in range(30)]
            for key in sorted({primer, reverse_complement(primer)}):
                low.create_dataset(
                    key,
                    data=np.array(small_hits if key == primer else [], dtype=np.int64),
                )
                high.create_dataset(
                    key,
                    data=np.array(large_hits if key == primer else [], dtype=np.int64),
                )

    return {"primers": primers, "fg": fg, "small": small, "large": large}


def _optimizer(world, monkeypatch, aggregate, target=5):
    monkeypatch.setattr(parameter, "background_aggregate", aggregate, raising=False)
    cache = PositionCache(
        [world["fg"], world["small"], world["large"]], world["primers"], on_missing="warn"
    )
    return OptimizerFactory.create(
        name="background-aware",
        position_cache=cache,
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[FG_LENGTH],
        bg_prefixes=[world["small"], world["large"]],
        bg_seq_lengths=[SMALL_HOST, LARGE_HOST],
        config=OptimizerConfig(target_set_size=target, extension_reach=3_000),
        conditions=None,
        bin_size=1_000,
        background_aggregate=aggregate,
    )


# ----------------------------------------------------------------------
# 1. Which stage produces the delivered panel
# ----------------------------------------------------------------------


def _record_stages(optimizer):
    """Wrap both candidate stages so each records what it was given and returned.

    Instrumentation rather than reading the source, because the branch that
    decides is a run-time one: Stage 1 can return at or below `final_count`, in
    which case `optimize` returns early and NEITHER of these stages runs.
    """
    inner = optimizer._hybrid
    seen: dict[str, dict] = {}
    original_prune = inner._prune_background
    original_refine = inner._network_refine

    def prune(primers, *args, **kwargs):
        out = original_prune(primers, *args, **kwargs)
        seen["prune"] = {"in": list(primers), "out": list(out[0])}
        return out

    def refine(primers, *args, **kwargs):
        out = original_refine(primers, *args, **kwargs)
        seen["refine"] = {"in": list(primers), "out": list(out)}
        return out

    inner._prune_background = prune
    inner._network_refine = refine
    return seen


def test_the_delivered_panel_comes_out_of_stage_two(world, monkeypatch):
    """Established, not assumed. This is what Known Issue 14 asks for.

    `_prune_background` (Stage 1.5) runs first and hands its survivors to
    `_network_refine` (Stage 2), whose output is the panel. So a host term that
    acted only in pruning could still be undone by Stage 2, and a claim about
    the term has to name Stage 2.

    Both stages charge the host through `_count_background_sites`, which is why
    one function carries the reduction: the two cannot read the key differently.
    """
    optimizer = _optimizer(world, monkeypatch, SUM)
    seen = _record_stages(optimizer)

    result = optimizer.optimize(world["primers"], target_size=5)

    assert result.primers, "the optimizer selected nothing; the world is wrong"
    assert "prune" in seen, "background pruning did not run; this is not the method under test"
    assert "refine" in seen, "Stage 2 did not run, so Stage 1 returned early"
    assert list(result.primers) == seen["refine"]["out"], (
        "the delivered panel is not Stage 2's output, so a claim about the "
        "Stage 2 host term would be about a stage that does not decide"
    )
    # And pruning feeds it rather than deciding: its output is Stage 2's input.
    assert seen["refine"]["in"] == seen["prune"]["out"]


# ----------------------------------------------------------------------
# 2. The charge that stage ranks on moves with the key
# ----------------------------------------------------------------------


def test_the_host_charge_moves_with_the_key(world, monkeypatch):
    """One function, read by both stages, and it answers differently.

    A candidate with three sites on the 30 kb host charges 3 pooled and
    3/30,000 x 1,320,000 = 132 at worst case; one with thirty sites on the
    1.29 Mb host charges 30 pooled and 30/1,290,000 x 1,320,000 = 31 at worst
    case. So the pooled reading calls the second candidate ten times worse and
    the per-host reading calls the first one four times worse -- the orders are
    opposite, which is what makes the delivered panel's movement attributable
    to the key rather than to noise.
    """
    pooled = _optimizer(world, monkeypatch, SUM)._hybrid
    worst = _optimizer(world, monkeypatch, WORST_CASE)._hybrid
    small_binder, large_binder = world["primers"][0], world["primers"][1]

    assert pooled._count_background_sites([small_binder]) == 3
    assert pooled._count_background_sites([large_binder]) == 30

    assert worst._count_background_sites([small_binder]) == 132
    assert worst._count_background_sites([large_binder]) == 31

    # The ranking is inverted, which is the whole claim.
    assert pooled._count_background_sites([small_binder]) < pooled._count_background_sites(
        [large_binder]
    )
    assert worst._count_background_sites([small_binder]) > worst._count_background_sites(
        [large_binder]
    )


def test_a_single_host_charges_the_same_under_both(world, monkeypatch):
    """Integer-identical, so a one-host `background-aware` run cannot move."""
    monkeypatch.setattr(parameter, "background_aggregate", None, raising=False)
    cache = PositionCache([world["fg"], world["small"]], world["primers"], on_missing="warn")
    built = {}
    for aggregate in (SUM, WORST_CASE):
        built[aggregate] = OptimizerFactory.create(
            name="background-aware",
            position_cache=cache,
            fg_prefixes=[world["fg"]],
            fg_seq_lengths=[FG_LENGTH],
            bg_prefixes=[world["small"]],
            bg_seq_lengths=[SMALL_HOST],
            config=OptimizerConfig(target_set_size=5),
            conditions=None,
            background_aggregate=aggregate,
        )._hybrid

    for primer in world["primers"]:
        assert built[SUM]._count_background_sites([primer]) == built[
            WORST_CASE
        ]._count_background_sites([primer])


def test_the_resolved_value_reaches_the_inner_optimizer(world, monkeypatch):
    """The wrapper resolves it; the object that charges has to hold it.

    Dropping it at this seam is how a compositional background came to be inert
    on every optimizer, which is recorded in `test_host_profile.py`.
    """
    for aggregate in (SUM, WORST_CASE):
        optimizer = _optimizer(world, monkeypatch, aggregate)
        assert optimizer.background_aggregate == aggregate
        assert optimizer._hybrid.background_aggregate == aggregate


def test_a_bare_hybrid_optimizer_defaults_to_sum():
    """Not `worst-case`: the class default must agree with the key's default."""
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    assert HybridOptimizer.background_aggregate == SUM


# ----------------------------------------------------------------------
# 3. The delivered panel moves
# ----------------------------------------------------------------------


def test_the_delivered_panel_moves_with_the_key(world, monkeypatch):
    """The assertion the step is for, on the stage established above.

    Compared against the key ABSENT as well as against `sum`, so "nothing moved
    by default" and "something moved on request" are separate claims.
    """
    absent = _optimizer(world, monkeypatch, None).optimize(world["primers"], target_size=5)
    pooled = _optimizer(world, monkeypatch, SUM).optimize(world["primers"], target_size=5)
    worst = _optimizer(world, monkeypatch, WORST_CASE).optimize(world["primers"], target_size=5)

    assert tuple(absent.primers) == tuple(pooled.primers), "the default moved"
    assert tuple(worst.primers) != tuple(pooled.primers), (
        "the host term did not change the delivered panel, so `worst-case` is "
        "read and not acted on -- Known Issue 14's exact shape"
    )
    assert set(worst.primers) <= set(world["primers"])


def test_the_worst_case_panel_is_cleaner_on_the_worst_host(world, monkeypatch):
    """Moving is not enough; it has to move the right way.

    The quantity the key is about is the worst host's density, so that is what
    is compared -- not the pooled figure, which the panel was not selected on.
    """
    pooled = _optimizer(world, monkeypatch, SUM)
    worst = _optimizer(world, monkeypatch, WORST_CASE)
    pooled_panel = pooled.optimize(world["primers"], target_size=5).primers
    worst_panel = worst.optimize(world["primers"], target_size=5).primers

    # Measured through ONE optimizer, so the comparison is between panels
    # rather than between two readings of the metric.
    scorer = _optimizer(world, monkeypatch, WORST_CASE)
    before = scorer.compute_metrics(list(pooled_panel)).worst_host_selectivity_density
    after = scorer.compute_metrics(list(worst_panel)).worst_host_selectivity_density

    assert before is not None and after is not None
    assert after > before, (
        "the worst-case panel is no better on the worst host than the pooled "
        "one, so the term changed the answer without improving what it is for"
    )
