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
