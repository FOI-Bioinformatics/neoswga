"""Thresholds configured in params.json must reach the optimizer.

`run_optimization` builds an `OptimizerConfig` and hands it to the optimizer
factory. It populated five of the dataclass's fields and left the rest at
their defaults, so `max_dimer_bp`, `max_self_dimer_bp`, `min_tm` and `max_tm`
never arrived however they were configured. Every optimizer read the dataclass
default instead.

This was invisible in unit tests, which construct the config directly and
therefore pass whatever they mean to test. It only shows end to end, and it
shows worst on the clique optimizer, where the consequence is not a slightly
different score but a broken guarantee:

    params.json          max_dimer_bp: 3
    OptimizerConfig      max_dimer_bp: 4   (the dataclass default)

4 is the LOOSER threshold -- fewer pairs count as dimers -- so the
compatibility graph gains edges, and a clique in that graph can contain a pair
that dimerises at the threshold the user asked for. `--optimization-method
clique` exists precisely to guarantee no such pair, and it returned sets
containing them.

The direction matters and is worth stating: a stricter-than-configured
threshold would have been safe, because a clique in a sparser graph is still a
clique in a denser one. Only the loose direction breaks the promise.
"""

import itertools
import random

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("networkx")

import numpy as np

from neoswga.core.dimer import is_dimer_fast
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 40_000
PRIMER_LENGTH = 10


@pytest.fixture(scope="module")
def planted():
    """A pool with enough mutual incompatibility that the threshold matters.

    A pool where almost nothing dimerises cannot distinguish a threshold of 3
    from one of 4, and the test would pass without checking anything.
    """
    rng = random.Random(20260802)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    for i in range(24):
        primer = "".join(rng.choice("ACGT") for _ in range(PRIMER_LENGTH))
        if primer not in primers:
            primers.append(primer)

    for i, primer in enumerate(primers):
        for j in range(3):
            pos = (i * 1_500 + j * 480) % (GENOME_LENGTH - 20)
            seq[pos : pos + PRIMER_LENGTH] = list(primer)

    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture(scope="module")
def prefix(tmp_path_factory, planted):
    path = str(tmp_path_factory.mktemp("cfg") / "target")
    keys = set()
    for primer in planted["primers"]:
        keys.update({primer, reverse_complement(primer)})

    with h5py.File(f"{path}_{PRIMER_LENGTH}mer_positions.h5", "w") as f:
        for key in sorted(keys):
            positions, i = [], planted["seq"].find(key)
            while i != -1:
                positions.append(i)
                i = planted["seq"].find(key, i + 1)
            if positions:
                f.create_dataset(key, data=np.array(positions, dtype=np.int32))
    return path


def run(prefix, planted, **kwargs):
    from neoswga.core.unified_optimizer import run_optimization

    return run_optimization(
        candidates=list(planted["primers"]),
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        target_size=6,
        verbose=False,
        extension_reach=3000,
        **kwargs,
    )


# ----------------------------------------------------------------------
# The guarantee, through the real dispatch path
# ----------------------------------------------------------------------


@pytest.mark.parametrize("max_dimer_bp", [3, 4, 5])
def test_clique_result_is_dimer_free_at_the_configured_threshold(prefix, planted, max_dimer_bp):
    """The property `--optimization-method clique` is chosen for.

    Checked at several thresholds because the failure was direction-dependent:
    the optimizer silently used 4, so a run configured with 3 got a set that
    was compatible at 4 and not at 3.
    """
    result = run(prefix, planted, method="clique", max_dimer_bp=max_dimer_bp)
    assert result.primers, f"clique returned nothing: {result.message}"

    offenders = [
        (a, b)
        for a, b in itertools.combinations(result.primers, 2)
        if is_dimer_fast(a, b, max_dimer_bp)
    ]
    assert not offenders, (
        f"clique set contains {len(offenders)} dimerising pair(s) at "
        f"max_dimer_bp={max_dimer_bp}: {offenders}"
    )


def test_a_stricter_threshold_is_not_silently_ignored(prefix, planted):
    """Tightening the threshold has to be able to change the answer.

    If it cannot, the value is not reaching the compatibility graph -- which is
    exactly what was happening.
    """
    loose = run(prefix, planted, method="clique", max_dimer_bp=6)
    strict = run(prefix, planted, method="clique", max_dimer_bp=2)

    assert loose.primers, "the loose threshold should find a set"

    # A threshold strict enough that no compatible set of the requested size
    # exists is a legitimate outcome, and refusing is the right answer -- the
    # alternative is returning a set that does not meet the constraint, which
    # is the bug this file is about. What must not happen is the strict run
    # quietly producing the same answer as the loose one.
    if strict.primers:
        assert set(strict.primers) != set(loose.primers) or strict.score != loose.score
    else:
        assert "dimer" in strict.message.lower() or "target size" in strict.message.lower()


# ----------------------------------------------------------------------
# The config plumbing itself
# ----------------------------------------------------------------------


CONFIGURABLE = ["max_dimer_bp", "max_self_dimer_bp", "min_tm", "max_tm"]


@pytest.mark.parametrize("field", CONFIGURABLE)
def test_optimizer_config_carries_the_field(field):
    """These are all settable in params.json, so the config must have a home
    for each. `max_self_dimer_bp` had none, and the clique optimizer's
    `getattr(self.config, 'max_self_dimer_bp', ...)` fell back every time.
    """
    import dataclasses

    from neoswga.core.base_optimizer import OptimizerConfig

    names = {f.name for f in dataclasses.fields(OptimizerConfig)}
    assert field in names


@pytest.mark.parametrize("field,value", [("max_dimer_bp", 2), ("min_tm", 12.0), ("max_tm", 55.0)])
def test_run_optimization_forwards_the_field(monkeypatch, prefix, planted, field, value):
    """Capture the config the factory actually receives.

    Asserting on the returned primer set would be indirect and could pass for
    the wrong reason; this checks the value arrives.
    """
    from neoswga.core import unified_optimizer as uo

    seen = {}
    original = uo.OptimizerFactory.create

    def spy(*args, **kwargs):
        config = kwargs.get("config")
        if config is not None:
            seen["config"] = config
        return original(*args, **kwargs)

    monkeypatch.setattr(uo.OptimizerFactory, "create", spy)

    run(prefix, planted, method="hybrid", **{field: value})

    assert "config" in seen, "no config reached OptimizerFactory.create"
    assert getattr(seen["config"], field) == value, (
        f"{field}={value} did not reach the optimizer; " f"got {getattr(seen['config'], field)}"
    )


def test_params_json_values_reach_the_config(monkeypatch, prefix, planted):
    """The CLI route: params.json -> parameter globals -> OptimizerConfig.

    A user setting max_dimer_bp in params.json and passing no flag should get
    that value, not the dataclass default.
    """
    from neoswga.core import parameter
    from neoswga.core import unified_optimizer as uo

    monkeypatch.setattr(parameter, "max_dimer_bp", 2, raising=False)

    seen = {}
    original = uo.OptimizerFactory.create

    def spy(*args, **kwargs):
        if kwargs.get("config") is not None:
            seen["config"] = kwargs["config"]
        return original(*args, **kwargs)

    monkeypatch.setattr(uo.OptimizerFactory, "create", spy)
    run(prefix, planted, method="hybrid")

    assert seen["config"].max_dimer_bp == 2


def test_network_optimizer_receives_the_configured_dimer_threshold(monkeypatch, prefix, planted):
    """The config reaches the adapter and stops there.

    `NetworkBaseOptimizer.__init__` builds the inner `NetworkOptimizer` from an
    explicit argument list that omits `max_dimer_bp`, so the object that
    actually scores pairs fell back to its own default of 4 whatever
    params.json asked for. As in the module docstring, 4 is the looser
    threshold: pairs the user wanted treated as dimers scored as clean.

    Asserting on the inner object rather than the config is the point. The
    config already carried the value; the gap is one constructor further in.
    """
    from neoswga.core import unified_optimizer as uo

    seen = {}
    original = uo.OptimizerFactory.create

    def spy(*args, **kwargs):
        optimizer = original(*args, **kwargs)
        seen["optimizer"] = optimizer
        return optimizer

    monkeypatch.setattr(uo.OptimizerFactory, "create", spy)
    run(prefix, planted, method="network", max_dimer_bp=3)

    optimizer = seen.get("optimizer")
    assert optimizer is not None, "no optimizer reached the spy"
    assert optimizer._network.max_dimer_bp == 3, (
        "configured max_dimer_bp=3 did not reach the inner NetworkOptimizer; "
        f"it is using {optimizer._network.max_dimer_bp}"
    )


# ----------------------------------------------------------------------
# The same question, asked of every method and every scoring parameter
# ----------------------------------------------------------------------
#
# The test above asks it once, for `network` and `max_dimer_bp`, and that is
# how the ninth instance of this defect class was found. The tenth, eleventh
# and twelfth were sitting one constructor further along the same path and went
# unnoticed because nothing asked the same question of `hybrid` -- the DEFAULT
# method -- or of the other four parameters.
#
# What reaches the user is decided by the object that scores primers, not by
# the adapter the factory returns and not by the config the adapter holds. So
# these assert on the innermost engine, and they are parameterised so that a
# newly wired parameter, or a newly registered method, produces a named failing
# test rather than a silent drop.
#
# Measured on this tree before the fix (--application clinical, max_dimer_bp 3,
# --max-extension 25000):
#
#     parameter          requested  network  hybrid  background-aware
#     tm_weight               0.30     0.30    0.00              0.00
#     dimer_penalty           0.45     0.45    0.00              0.00
#     uniformity_weight       0.05     0.05    0.00              0.00
#     max_dimer_bp               3        3       4                 4
#     max_extension          25000    25000   25000             70000
#
# while the run log said "selection weights applied: tm=0.30, uniformity=0.05,
# dimer_penalty=0.45".


#: Methods whose selection is ultimately performed by a `NetworkOptimizer`,
#: and the attribute path from the factory-returned adapter to it.
_SCORING_ENGINE_PATH = {
    "network": ("_network",),
    "hybrid": ("_hybrid", "network_optimizer"),
    "background-aware": ("_hybrid", "network_optimizer"),
}


def scoring_engine(optimizer, method):
    """The object that actually scores primers for `method`.

    Asserting on the adapter, or on the `OptimizerConfig` it holds, is what let
    these defects through: the value arrives there and is then dropped on the
    way in.
    """
    engine = optimizer
    for attribute in _SCORING_ENGINE_PATH[method]:
        engine = getattr(engine, attribute, None)
        assert engine is not None, (
            f"{method}: could not reach the scoring engine, "
            f"{type(optimizer).__name__} has no '{attribute}'. If the "
            f"attribute was renamed, update _SCORING_ENGINE_PATH."
        )
    return engine


def build(prefix, planted, method, **kwargs):
    """Construct `method` through the real dispatch path and return it."""
    from neoswga.core import unified_optimizer as uo

    seen = {}
    original = uo.OptimizerFactory.create

    def spy(*args, **create_kwargs):
        optimizer = original(*args, **create_kwargs)
        seen["optimizer"] = optimizer
        return optimizer

    uo.OptimizerFactory.create = staticmethod(spy)
    try:
        run(prefix, planted, method=method, **kwargs)
    finally:
        uo.OptimizerFactory.create = original

    assert "optimizer" in seen, f"{method}: no optimizer reached the spy"
    return seen["optimizer"]


#: (parameter, configured value, attribute on the scoring engine).
#: Every one of these is settable from params.json or a CLI flag, and every one
#: is read by `NetworkOptimizer` when it scores a candidate.
SCORING_PARAMETERS = [
    ("tm_weight", 0.30, "tm_weight"),
    ("dimer_penalty", 0.45, "dimer_penalty"),
    ("uniformity_weight", 0.05, "uniformity_weight"),
    ("max_dimer_bp", 3, "max_dimer_bp"),
    ("max_extension", 25_000, "max_extension"),
]


@pytest.mark.parametrize("method", sorted(_SCORING_ENGINE_PATH))
@pytest.mark.parametrize("name,value,attribute", SCORING_PARAMETERS, ids=lambda p: str(p))
def test_the_scoring_engine_uses_the_configured_value(
    prefix, planted, method, name, value, attribute
):
    """A parameter the user set has to reach the object that scores primers.

    `dimer_penalty` is the one worth reading twice: it is the ONLY dimer
    consideration in hybrid's stage-2 refinement, so when it does not arrive the
    default method weighs dimers at zero whatever `--application` reports.
    """
    optimizer = build(prefix, planted, method, **{name: value})
    engine = scoring_engine(optimizer, method)

    actual = getattr(engine, attribute, None)
    assert actual == value, (
        f"{method}: {name}={value} did not reach {type(engine).__name__}."
        f"{attribute}; it is using {actual!r}. The value reaches the adapter and "
        f"is dropped when the adapter builds its inner engine."
    )


# ----------------------------------------------------------------------
# Optimizer-specific configuration
# ----------------------------------------------------------------------
#
# Two optimizers declare their own `OptimizerConfig` subclass and then branch on
# `isinstance` to read it. `OptimizerFactory.create` always built the BASE class,
# so both branches were constant-false in production and every subclass field sat
# at its dataclass default, unreachable from params.json or the CLI.
#
# `CliqueOptimizerConfig.max_pool_size` is the one that costs something: it
# defaults to 200, so with the shipped `max_primer` of 500 the optimizer sold on
# a hard dimer-free guarantee silently discarded 300 candidates and ranked the
# rest on raw foreground abundance, with no background term.
#
# `CliqueOptimizerConfig` is instantiated nowhere outside tests/test_clique_optimizer.py,
# which is why nothing noticed.


def test_the_clique_optimizer_receives_its_own_config_class():
    """The factory has to build the subclass the optimizer branches on.

    A BASE `OptimizerConfig` is passed deliberately, because that is what
    production does -- `run_optimization` always builds the base class. Passing
    `None` here would take the `config or CliqueOptimizerConfig()` default in
    `__init__` and the test would pass against the broken factory.
    """
    from neoswga.core import unified_optimizer as uo
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.clique_optimizer import CliqueOptimizerConfig
    from neoswga.core.optimizer_factory import OptimizerFactory

    uo._ensure_optimizers_registered()

    class _Cache:
        def get_positions(self, *args, **kwargs):
            import numpy as np

            return np.array([], dtype=np.int64)

    optimizer = OptimizerFactory.create(
        name="clique",
        position_cache=_Cache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[10_000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=6, max_dimer_bp=3),
        conditions=None,
    )
    assert optimizer.config.max_dimer_bp == 3, (
        "upgrading the config to the optimizer's own subclass must carry the "
        "base fields across, not reset them to defaults"
    )
    assert isinstance(optimizer.clique_config, CliqueOptimizerConfig), (
        "clique fell back to a throwaway default config; its isinstance branch "
        "is constant-false whenever the factory hands it a base OptimizerConfig"
    )


@pytest.mark.parametrize("value", [50, 400])
def test_max_pool_size_configured_in_params_json_reaches_the_clique_optimizer(
    monkeypatch, prefix, planted, value
):
    """A pool cap the user set has to be the cap the optimizer applies.

    Pinned at both a smaller and a larger value than the 200 default, because a
    test that only tightens it would pass against code that ignores the setting
    and truncates to 200 anyway whenever the pool is smaller than that.
    """
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "max_pool_size", value, raising=False)
    optimizer = build(prefix, planted, "clique")

    assert optimizer.clique_config.max_pool_size == value, (
        f"max_pool_size={value} did not reach the clique optimizer; "
        f"it is using {optimizer.clique_config.max_pool_size}"
    )
