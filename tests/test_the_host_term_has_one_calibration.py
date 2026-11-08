"""`background-aware`'s host term has ONE calibration, and refuses the other.

Phase 7c first made this optimizer's host term read `background_aggregate`,
charging the worst host's site density carried back onto the pooled length. The
ranking was right and the magnitude was not: that charge is the pooled one
multiplied by `sum(host lengths) / worst host length`, which is 44x in the world
below and 1029x on this project's own three-host panel, while
`background_weight` scales the term against a coverage term and was chosen
against pooled magnitudes. No measurement establishes the weight under the
per-host reading, so the combination is REFUSED rather than left reachable
behind a note. `per_host_metrics.refuse_uncalibrated_host_term` carries the
argument and names what would settle it.

The order here is deliberate, and the first part survives the withdrawal because
it is a fact about this optimizer that the record needs either way.

1. Establish WHICH stage produces the delivered panel, by instrumenting the two
   candidates and comparing their recorded outputs against the result. Known
   Issue 14 is that reading the background is not acting on it, and it was found
   twice in this optimizer, so any claim about a host term has to name the stage
   that decides -- including a claim that the term is pooled.
2. The charge has one reading, and the other raises: at `run_optimization`'s
   gate, at the constructor, and in `host_site_charge` itself.
3. The refusal is narrow. Every other method still honours the key, the per-host
   figures are still computed, and the configured limits still bind per host --
   refusing `worst-case` outright would have thrown away the phase.
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
from neoswga.core.exceptions import InvalidDesignRequest
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
    with (
        h5py.File(f"{small}_10mer_positions.h5", "w") as low,
        h5py.File(f"{large}_10mer_positions.h5", "w") as high,
    ):
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
# 2. The charge has ONE calibration, and the other reading is refused
# ----------------------------------------------------------------------


def test_the_host_charge_is_the_pooled_sum_under_sum(world, monkeypatch):
    """The reading this term is calibrated for, unchanged from `main`."""
    pooled = _optimizer(world, monkeypatch, SUM)._hybrid
    small_binder, large_binder = world["primers"][0], world["primers"][1]

    assert pooled._count_background_sites([small_binder]) == 3
    assert pooled._count_background_sites([large_binder]) == 30
    assert pooled._count_background_sites([small_binder, large_binder]) == 33


def test_worst_case_with_background_aware_and_two_hosts_is_refused(world, monkeypatch):
    """Refused at the gate every command goes through, not deep in a loop.

    `run_optimization` resolves the key and refuses before any optimizer is
    built, because `_dispatch_optimizer` converts a constructor exception into an
    `OptimizationResult.failure` and a `DesignError` must not become one.
    """
    monkeypatch.setattr(parameter, "background_aggregate", WORST_CASE, raising=False)

    with pytest.raises(InvalidDesignRequest) as caught:
        uo.run_optimization(
            method="background-aware",
            candidates=list(world["primers"]),
            cache=PositionCache(
                [world["fg"], world["small"], world["large"]],
                world["primers"],
                on_missing="warn",
            ),
            fg_prefixes=[world["fg"]],
            fg_seq_lengths=[FG_LENGTH],
            bg_prefixes=[world["small"], world["large"]],
            bg_seq_lengths=[SMALL_HOST, LARGE_HOST],
            target_size=5,
            verbose=False,
        )

    message = str(caught.value)
    assert "background-aware" in message
    assert "background_weight" in message
    # The message has to say what would settle it, or it is a dead end rather
    # than a refusal.
    assert "What would settle it" in message
    assert "worst-case" in message and "sum" in message


def test_the_constructor_refuses_too_for_a_caller_that_skipped_the_gate(world, monkeypatch):
    """`BackgroundAwareBaseOptimizer` is public and reachable from a library.

    Belt and braces: the gate above is the one every command crosses, and this
    is the one a direct construction crosses.
    """
    with pytest.raises(InvalidDesignRequest):
        _optimizer(world, monkeypatch, WORST_CASE)


def test_an_ensemble_that_would_run_it_is_refused_as_well(world, monkeypatch):
    """`ensemble` with no explicit list runs the default, which includes it.

    Refusing the named method and letting the ensemble run it would be the
    inert-refusal shape: the combination would still be reached, by the method
    list nobody typed.
    """
    from neoswga.core.per_host_metrics import refuse_uncalibrated_host_term

    for method in ("ensemble", "auto", "all"):
        with pytest.raises(InvalidDesignRequest):
            refuse_uncalibrated_host_term(method, WORST_CASE, 2)

    # An explicit list that leaves it out is fine.
    refuse_uncalibrated_host_term("ensemble", WORST_CASE, 2, ensemble_methods=["hybrid", "network"])


def test_the_default_ensemble_list_still_contains_the_refused_method():
    """The refusal names `background-aware` for an ensemble by assumption.

    If `_dispatch_optimizer`'s default list ever drops it, that assumption
    becomes a refusal of a combination that no longer occurs -- an inert
    refusal, which is Known Issue 8 pointing the other way. Read from the
    dispatch's own source so the two cannot drift silently.
    """
    import inspect

    source = inspect.getsource(uo._dispatch_optimizer)
    head = source[: source.index("ensemble_combine")]

    assert '"background-aware"' in head


def test_nothing_is_refused_with_one_host(world, monkeypatch):
    """One host makes the two readings the same number, so there is no question.

    Both the up-front gate and the constructor let it through, and the charge is
    integer-identical.
    """
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


def test_the_resolved_value_still_reaches_the_inner_optimizer(world, monkeypatch):
    """The wrapper resolves it; the object that charges has to hold it.

    Dropping it at this seam is how a compositional background came to be inert
    on every optimizer, which `test_host_profile.py` records. Asserted at one
    host, where the value is not refused.
    """
    monkeypatch.setattr(parameter, "background_aggregate", None, raising=False)
    cache = PositionCache([world["fg"], world["small"]], world["primers"], on_missing="warn")
    for aggregate in (SUM, WORST_CASE):
        built = OptimizerFactory.create(
            name="background-aware",
            position_cache=cache,
            fg_prefixes=[world["fg"]],
            fg_seq_lengths=[FG_LENGTH],
            bg_prefixes=[world["small"]],
            bg_seq_lengths=[SMALL_HOST],
            config=OptimizerConfig(target_set_size=5),
            conditions=None,
            background_aggregate=aggregate,
        )
        assert built.background_aggregate == aggregate
        assert built._hybrid.background_aggregate == aggregate


def test_a_bare_hybrid_optimizer_defaults_to_sum():
    """Not `worst-case`: the class default must agree with the key's default."""
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    assert HybridOptimizer.background_aggregate == SUM


# ----------------------------------------------------------------------
# 3. The refusal is narrow: every other method still honours the key
# ----------------------------------------------------------------------


def test_another_method_runs_under_worst_case_with_two_hosts(world, monkeypatch):
    """The refusal is about ONE optimizer's weighted host term, nothing else.

    `hybrid` has no host term in its selection -- `background_pruning` is off,
    so `_count_background_sites` is never consulted -- and the per-host figures
    and every configured panel limit still bind per host there. Refusing
    `worst-case` outright would have thrown away the phase.
    """
    monkeypatch.setattr(parameter, "background_aggregate", WORST_CASE, raising=False)
    cache = PositionCache(
        [world["fg"], world["small"], world["large"]], world["primers"], on_missing="warn"
    )
    optimizer = OptimizerFactory.create(
        name="hybrid",
        position_cache=cache,
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[FG_LENGTH],
        bg_prefixes=[world["small"], world["large"]],
        bg_seq_lengths=[SMALL_HOST, LARGE_HOST],
        config=OptimizerConfig(target_set_size=5, extension_reach=3_000),
        conditions=None,
        bin_size=1_000,
        background_aggregate=WORST_CASE,
    )

    assert optimizer.background_aggregate == WORST_CASE
    assert optimizer._hybrid.background_pruning is False

    metrics = optimizer.compute_metrics(list(world["primers"])[:4])
    assert metrics.background_aggregate == WORST_CASE
    assert len(metrics.per_host) == 2
    # And the reduction is there to be read, which is what the limits use.
    assert metrics.worst_host_background_sites is not None
