"""Choosing reaction conditions for specificity, from a measurement.

Two things in the codebase claimed to optimise conditions for specificity and
neither looked at a genome.

`additive_optimizer._evaluate_combination` scored its "specificity" objective
as `effective_binding_rate * (1 / koff_factor)` of a **synthetic** primer built
from a length and a GC fraction -- literally `"GG...CC...AA...TT"`. No primer
set, no background, nothing that could distinguish a specific design from an
indiscriminate one. It measured binding stability, which is not specificity;
a primer that binds everything scores well on it.

The `optimize-conditions` CLI did not use that optimiser at all. It was a
GC-content lookup table (`reaction_conditions.recommend_conditions`) taking
`--fg` with no `--bg`, so it could not have measured specificity even in
principle.

With occupancy-weighted selectivity available, the question becomes answerable:
sweep a temperature x additive grid, compute selectivity against the actual
background for the actual primer set, and report it against the amplification
it costs. Both are needed. Stringency buys discrimination by pushing duplexes
off their transition, which is the same thing as pushing them toward not being
bound at all -- so an optimum exists and a search reporting only selectivity
would walk straight past it to a reaction that amplifies nothing.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.condition_sweep import ConditionPoint, sweep_conditions
from neoswga.core.mismatch_counts import canonical_kmer
from neoswga.core.thermodynamics import reverse_complement

GENOME = 60_000
K = 10


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Foreground and background with a Tm-coherent primer set.

    Tm coherence matters: a set melting 25 C apart has no temperature at which
    all of it is near its transition, so no condition can be optimal for it and
    the sweep would be measuring nothing. See test_occupancy_selectivity.
    """
    rng = random.Random(99)
    directory = tmp_path_factory.mktemp("sweep")

    fg_seq = [rng.choice("ACGT") for _ in range(GENOME)]
    bg_seq = [rng.choice("ACGT") for _ in range(GENOME)]

    primers = []
    for i in range(6):
        bases = ["G", "C"] * 2 + ["A", "T"] * 3
        rng.shuffle(bases)
        primer = "".join(bases) + "".join(rng.choice("AT") for _ in range(K - len(bases)))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(4):
            at = (i * 7_000 + j * 800) % (GENOME - 20)
            fg_seq[at : at + K] = list(primer)
        bg_seq[i * 6_000 : i * 6_000 + K] = list(primer)
        for m in range(3):
            variant = list(primer)
            variant[m] = "ACGT"[("ACGT".index(variant[m]) + 1) % 4]
            at = i * 6_000 + 400 + m * 120
            bg_seq[at : at + K] = variant

    prefixes = {}
    for name, seq in (("fg", "".join(fg_seq)), ("bg", "".join(bg_seq))):
        prefix = str(directory / name)
        prefixes[name] = prefix
        counts = {}
        for i in range(len(seq) - K + 1):
            key = canonical_kmer(seq[i : i + K])
            counts[key] = counts.get(key, 0) + 1
        with open(f"{prefix}_{K}mer_all.txt", "w") as handle:
            for kmer, count in counts.items():
                handle.write(f"{kmer} {count}\n")
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as handle:
            for key in sorted({x for p in primers for x in (p, reverse_complement(p))}):
                hits, at = [], seq.find(key)
                while at != -1:
                    hits.append(at)
                    at = seq.find(key, at + 1)
                if hits:
                    handle.create_dataset(key, data=np.array(hits, dtype=np.int32))

    return {"primers": primers, "fg": prefixes["fg"], "bg": prefixes["bg"]}


def run_sweep(world, **kwargs):
    return sweep_conditions(
        primers=world["primers"],
        fg_prefixes=[world["fg"]],
        bg_prefixes=[world["bg"]],
        polymerase="phi29",
        **kwargs,
    )


# ----------------------------------------------------------------------
# It measures something real
# ----------------------------------------------------------------------


def test_the_sweep_covers_the_grid_it_was_given(world):
    points = run_sweep(world, temperatures=[28.0, 32.0], dmso_percentages=[0.0, 4.0])

    assert len(points) == 4
    assert all(isinstance(p, ConditionPoint) for p in points)
    assert {(p.temp, p.dmso_percent) for p in points} == {
        (28.0, 0.0),
        (28.0, 4.0),
        (32.0, 0.0),
        (32.0, 4.0),
    }


def test_every_point_carries_both_halves_of_the_trade(world):
    """Selectivity alone would send a search to a reaction that amplifies
    nothing, since the way stringency buys discrimination is by pushing
    duplexes toward not being bound."""
    for point in run_sweep(world, temperatures=[30.0], dmso_percentages=[0.0]):
        assert point.selectivity > 0
        assert point.amplification > 0


def test_selectivity_responds_to_the_grid(world):
    """The whole point. Under the old objective this was a property of a
    synthetic primer and could not vary with a real background at all."""
    points = run_sweep(world, temperatures=[26.0, 30.0, 34.0, 38.0], dmso_percentages=[0.0])
    ratios = [p.selectivity for p in sorted(points, key=lambda p: p.temp)]

    assert len(set(ratios)) > 1, "selectivity was flat across the temperature grid"
    for cooler, hotter in zip(ratios, ratios[1:]):
        assert hotter > cooler


def test_an_additive_moves_the_answer(world):
    points = run_sweep(world, temperatures=[30.0], dmso_percentages=[0.0, 5.0])
    by_dmso = {p.dmso_percent: p.selectivity for p in points}

    assert by_dmso[5.0] != by_dmso[0.0]


# ----------------------------------------------------------------------
# The trade-off is visible rather than implied
# ----------------------------------------------------------------------


def test_amplification_falls_as_stringency_rises(world):
    """The cost side. If it did not fall, there would be no trade to report
    and the recommendation would just be "as hot as possible"."""
    points = sorted(run_sweep(world, temperatures=[24.0, 32.0, 39.0]), key=lambda p: p.temp)
    amps = [p.amplification for p in points]

    assert amps[-1] < amps[0]


def test_the_recommendation_is_not_simply_the_hottest_point(world):
    """A search on selectivity alone walks to the end of the grid every time.

    The recommendation weighs the amplification it costs, so it should sit
    inside the range rather than at its edge for a grid that spans from
    saturated to melted.
    """
    from neoswga.core.condition_sweep import recommend_conditions_for_set

    temperatures = [22.0, 26.0, 30.0, 34.0, 38.0, 40.0]
    best = recommend_conditions_for_set(
        primers=world["primers"],
        fg_prefixes=[world["fg"]],
        bg_prefixes=[world["bg"]],
        polymerase="phi29",
        temperatures=temperatures,
        dmso_percentages=[0.0],
    )

    assert best.temp != temperatures[-1], (
        "the recommendation went to the hottest point on the grid, which means "
        "the amplification cost is not being weighed"
    )


def test_a_report_can_be_rendered(world):
    from neoswga.core.condition_sweep import format_sweep_table

    points = run_sweep(world, temperatures=[30.0, 36.0], dmso_percentages=[0.0, 4.0])
    text = format_sweep_table(points)

    assert "selectivity" in text.lower()
    assert "amplification" in text.lower()
    assert text.count("\n") >= len(points)


# ----------------------------------------------------------------------
# One implementation of the weighted load
# ----------------------------------------------------------------------


def test_the_sweep_and_the_optimizer_agree(world):
    """Three call sites now want an occupancy-weighted load: the optimizer's
    metrics, `evaluate-set`, and this sweep. They share one function, because
    three implementations of one number is how they come to disagree.
    """
    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import ReactionConditions

    temp = 30.0
    cache = PositionCache([world["fg"], world["bg"]], world["primers"])
    optimizer = DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[world["bg"]],
        bg_seq_lengths=[GENOME],
        config=OptimizerConfig(target_set_size=4),
        conditions=ReactionConditions(temp=temp, polymerase="phi29"),
    )
    from_metrics = optimizer.compute_metrics(world["primers"]).selectivity_ratio
    from_sweep = run_sweep(world, temperatures=[temp], dmso_percentages=[0.0])[0].selectivity

    assert from_sweep == pytest.approx(from_metrics, rel=1e-6)


def test_a_grid_past_the_enzyme_range_is_trimmed_and_reported(world, caplog):
    """phi29 stops at 40 C, and a sweep is exactly where a wide grid gets
    specified. Skipping is right; skipping silently is not -- an optimum at the
    edge of a quietly shortened grid reads as a real optimum.
    """
    import logging

    with caplog.at_level(logging.WARNING):
        points = run_sweep(world, temperatures=[30.0, 36.0, 55.0], dmso_percentages=[0.0])

    assert len(points) == 2
    assert "operating range" in caplog.text


def test_a_grid_entirely_outside_the_range_raises(world):
    """Returning an empty sweep would leave the caller to infer why."""
    with pytest.raises(ValueError, match="operating range"):
        run_sweep(world, temperatures=[70.0, 80.0], dmso_percentages=[0.0])


def test_a_fully_saturated_grid_says_the_ranking_is_amplification_only():
    """On a small background every primer can have zero near-matches, so every
    point returns MAX_SELECTIVITY and the ordering is decided entirely by
    amplification.

    That is a legitimate answer, but a table of identical selectivity values
    reads as a specificity result unless it says otherwise -- and
    MAX_SELECTIVITY printed as 1000000.00 reads as a measurement rather than a
    stand-in for "unbounded".
    """
    from neoswga.core.base_optimizer import MAX_SELECTIVITY
    from neoswga.core.condition_sweep import ConditionPoint, format_sweep_table

    saturated = [
        ConditionPoint(t, 0.0, 0.0, MAX_SELECTIVITY, amp, 5.0, 0.0)
        for t, amp in ((28.0, 0.4), (34.0, 0.5))
    ]
    text = format_sweep_table(saturated)

    assert "no bg binding" in text
    assert "amplification alone" in text
    assert "1000000" not in text


def test_the_additive_optimizer_can_use_a_real_selectivity_measurement():
    """`optimize_for="specificity"` scored a synthetic primer's binding
    stability -- which an indiscriminate primer also scores well on, since
    binding tightly to everything is still binding tightly. It can now be given
    a real measurement.
    """
    from neoswga.core.additive_optimizer import AdditiveOptimizer

    seen = []

    def fake_selectivity(conditions):
        seen.append(conditions.temp)
        return 42.0

    optimizer = AdditiveOptimizer(polymerase="phi29", selectivity_fn=fake_selectivity)
    score, _effects, _warnings = optimizer._evaluate_combination(
        {"dmso_percent": 0.0},
        primer_length=10,
        template_gc=0.5,
        primer_gc=0.5,
        optimize_for="specificity",
    )

    assert score == 42.0
    assert seen, "the selectivity function was not consulted"


def test_the_proxy_objective_says_it_is_a_proxy(caplog):
    """Without a measurement the objective still runs, but a caller must not
    read its output as specificity."""
    import logging

    from neoswga.core.additive_optimizer import AdditiveOptimizer

    optimizer = AdditiveOptimizer(polymerase="phi29")
    with caplog.at_level(logging.WARNING):
        optimizer._evaluate_combination(
            {"dmso_percent": 0.0},
            primer_length=10,
            template_gc=0.5,
            primer_gc=0.5,
            optimize_for="specificity",
        )

    assert "not foreground" in caplog.text or "synthetic primer" in caplog.text
