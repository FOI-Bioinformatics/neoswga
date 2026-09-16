"""The focused evaluator must agree with the full one, exactly.

Two definitions of the same quantity drift while each stays self-consistent.
`compute_pool_metrics` exists only because `compute_metrics` computes a great
deal the pool search never reads, so the moment the two disagree on a field they
both supply, the cheaper one is wrong.

Panels are drawn at random rather than hand-picked: a fixed panel would pin the
one composition that happened to agree when this was written.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_metrics import compute_pool_metrics
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

GENOME_LENGTH = 60_000
K = 10
FIELDS = (
    "fg_coverage",
    "effective_fg_coverage",
    "selectivity_density",
    "total_bg_sites",
    "max_gap",
)


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """A target and a background, both with real sites for every candidate."""
    rng = random.Random(20260916)
    target = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    host = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers, seen = [], set()
    while len(primers) < 24:
        p = "".join(rng.choice("ACGT") for _ in range(K))
        if p in seen or reverse_complement(p) in seen:
            continue
        seen.add(p)
        primers.append(p)
        for _ in range(rng.randint(2, 9)):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            target[pos : pos + K] = list(p)
        for _ in range(rng.randint(0, 4)):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            host[pos : pos + K] = list(p)

    tmp = tmp_path_factory.mktemp("pool_metrics")
    prefixes = {}
    for name, seq in (("target", "".join(target)), ("host", "".join(host))):
        prefix = str(tmp / name)
        prefixes[name] = prefix
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as f:
            for p in primers:
                for key in {p, reverse_complement(p)}:
                    hits, i = [], seq.find(key)
                    while i != -1:
                        hits.append(i)
                        i = seq.find(key, i + 1)
                    f.create_dataset(key, data=np.array(hits, dtype=np.int64))

    _ensure_optimizers_registered()
    cache = PositionCache([prefixes["target"], prefixes["host"]], primers)
    optimizer = OptimizerFactory.create(
        "hybrid",
        cache,
        [prefixes["target"]],
        [GENOME_LENGTH],
        [prefixes["host"]],
        [GENOME_LENGTH],
        config=OptimizerConfig(extension_reach=3_000, verbose=False),
        conditions=None,
        polymerase="phi29",
    )
    return optimizer, primers


@pytest.mark.parametrize("seed", range(12))
def test_a_random_panel_gets_the_same_five_numbers(built, seed):
    optimizer, primers = built
    rng = random.Random(seed)
    panel = rng.sample(primers, rng.randint(1, len(primers)))

    full = optimizer.compute_metrics(panel)
    light = compute_pool_metrics(optimizer, panel)

    for field in FIELDS:
        expected, actual = getattr(full, field), getattr(light, field)
        if expected is None or actual is None:
            assert expected is actual, f"{field}: {expected!r} against {actual!r}"
        else:
            assert actual == pytest.approx(expected), field


def test_an_empty_panel_agrees_too(built):
    optimizer, _ = built

    full = optimizer.compute_metrics([])
    light = compute_pool_metrics(optimizer, [])

    for field in FIELDS:
        expected, actual = getattr(full, field), getattr(light, field)
        if expected is None or actual is None:
            assert expected is actual, field
        else:
            assert actual == pytest.approx(expected), field


def test_reaching_for_a_sixth_field_is_an_error(built):
    """A default here would be indistinguishable from a measurement."""
    optimizer, primers = built
    light = compute_pool_metrics(optimizer, primers[:3])

    with pytest.raises(AttributeError):
        light.dimer_risk


@pytest.fixture(scope="module")
def with_conditions(built):
    """The same pair under real chemistry, so occupancy weighting is exercised.

    Without reaction conditions `effective_fg_coverage` is `None` on both sides
    and the agreement above says nothing about the occupancy path, which is the
    one the default coverage metric uses.
    """
    from neoswga.core.reaction_conditions import ReactionConditions

    optimizer, primers = built
    _ensure_optimizers_registered()
    return (
        OptimizerFactory.create(
            "hybrid",
            optimizer.cache,
            optimizer.fg_prefixes,
            optimizer.fg_seq_lengths,
            optimizer.bg_prefixes,
            optimizer.bg_seq_lengths,
            config=OptimizerConfig(extension_reach=3_000, verbose=False),
            conditions=ReactionConditions(temp=30.0),
            polymerase="phi29",
        ),
        primers,
    )


@pytest.mark.parametrize("seed", range(8))
def test_occupancy_weighted_coverage_agrees_as_well(with_conditions, seed):
    optimizer, primers = with_conditions
    rng = random.Random(1000 + seed)
    panel = rng.sample(primers, rng.randint(1, len(primers)))

    full = optimizer.compute_metrics(panel)
    light = compute_pool_metrics(optimizer, panel)

    assert full.effective_fg_coverage is not None, "the occupancy path did not run"
    assert light.effective_fg_coverage == pytest.approx(full.effective_fg_coverage)
    for field in FIELDS:
        expected, actual = getattr(full, field), getattr(light, field)
        assert actual == pytest.approx(expected), field
