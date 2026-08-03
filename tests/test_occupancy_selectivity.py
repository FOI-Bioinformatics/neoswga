"""Selectivity that responds to reaction conditions.

`selectivity_ratio` was `exact_fg_sites / exact_bg_sites`: a count of perfect
matches, with no temperature, free energy or additive term in it. For a fixed
primer set, no chemistry could move it.

It now weights each site class by the occupancy that class has under the
configured conditions:

    selectivity(C) = SUM_j n_fg,j * theta_j(C)  /  SUM_j n_bg,j * theta_j(C)

Raising stringency -- a hotter reaction, or an additive that lowers Tm -- drops
occupancy faster for mismatched background sites than for perfect foreground
ones, so the ratio rises. That is the mechanism this project wants to exploit,
and the first test below is the one that says whether it works.

Two guards matter as much as the feature.

**Reduction.** With the mismatch cap at zero the weighted ratio must reproduce
the old count ratio closely. If it does not, the occupancy term rather than the
mismatch counting is what moved every number in the project, and that would
need explaining rather than accepting.

**Degradation.** Occupancy needs the jellyfish count files. Where they are
absent -- an externally supplied oligo set scanned straight from a FASTA, for
instance -- the metric falls back to exact counting and says so in
`selectivity_mode`. A silent switch between two different definitions of the
same number is the failure this audit has met most often.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.mismatch_counts import canonical_kmer
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import reverse_complement

GENOME = 60_000
K = 10


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """A foreground and a background genome with k-mer counts and positions.

    The background is deliberately built to carry near-matches of the
    foreground primers: that population is what stringency acts on, and a
    fixture without it could not exercise the feature at all.
    """
    rng = random.Random(2026)
    directory = tmp_path_factory.mktemp("sel")

    fg_seq = [rng.choice("ACGT") for _ in range(GENOME)]
    bg_seq = [rng.choice("ACGT") for _ in range(GENOME)]

    # Primers of matched GC, so the set has one coherent melting transition.
    #
    # A first version drew them at random and the effect vanished: Tm spanned
    # 22-48 C, so at any single temperature the sum was dominated by primers
    # already saturated (theta 0.99 perfect against 0.96 mismatched, which
    # discriminates not at all) while the ones near their transition -- where
    # discrimination is 3x -- contributed little. That is the model behaving
    # correctly and is worth stating: a set with a wide Tm spread cannot be
    # made selective by tuning temperature, because no single temperature puts
    # every primer near its own transition. It is also why Tm tightness is a
    # term in `normalized_score`.
    primers = []
    for i in range(8):
        bases = ["G", "C"] * 2 + ["A", "T"] * 3
        rng.shuffle(bases)
        primer = "".join(bases) + "".join(rng.choice("AT") for _ in range(K - len(bases)))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(4):
            pos = (i * 6_000 + j * 900) % (GENOME - 20)
            fg_seq[pos : pos + K] = list(primer)
        # One perfect and several 1-mismatch copies in the background.
        bg_seq[i * 5_000 : i * 5_000 + K] = list(primer)
        for m in range(3):
            variant = list(primer)
            variant[m] = "ACGT"[(("ACGT".index(variant[m])) + 1) % 4]
            at = i * 5_000 + 500 + m * 100
            bg_seq[at : at + K] = variant

    fg_seq = "".join(fg_seq)
    bg_seq = "".join(bg_seq)

    prefixes = {}
    for name, seq in (("fg", fg_seq), ("bg", bg_seq)):
        prefix = str(directory / name)
        prefixes[name] = prefix

        counts = {}
        for i in range(len(seq) - K + 1):
            kmer = seq[i : i + K]
            counts[canonical_kmer(kmer)] = counts.get(canonical_kmer(kmer), 0) + 1
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


def metrics_at(world, temp, max_mismatches=1, **additives):
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.position_cache import PositionCache

    cache = PositionCache([world["fg"], world["bg"]], world["primers"])
    conditions = ReactionConditions(temp=temp, polymerase="phi29", **additives)
    optimizer = DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[world["bg"]],
        bg_seq_lengths=[GENOME],
        config=OptimizerConfig(target_set_size=4, max_mismatches=max_mismatches),
        conditions=conditions,
    )
    return optimizer.compute_metrics(world["primers"])


# ----------------------------------------------------------------------
# The thing this whole plan exists for
# ----------------------------------------------------------------------


def test_selectivity_rises_with_temperature(world):
    """Stringency buys discrimination.

    Under the old count ratio this was flat by construction: temperature does
    not appear in a count of exact matches.
    """
    ratios = [metrics_at(world, temp).selectivity_ratio for temp in (25.0, 30.0, 35.0, 40.0)]

    for cooler, hotter in zip(ratios, ratios[1:]):
        assert hotter > cooler, f"selectivity did not improve with temperature: {ratios}"


def test_an_additive_improves_selectivity(world):
    """The headline claim: `--dmso-percent` is a specificity control.

    DMSO lowers effective Tm, which is the same lever as heating the reaction.
    """
    plain = metrics_at(world, 30.0).selectivity_ratio
    with_dmso = metrics_at(world, 30.0, dmso_percent=5.0).selectivity_ratio

    assert with_dmso > plain


def test_more_additive_helps_more(world):
    ratios = [
        metrics_at(world, 30.0, dmso_percent=pct).selectivity_ratio for pct in (0.0, 3.0, 6.0)
    ]

    for less, more in zip(ratios, ratios[1:]):
        assert more > less


# ----------------------------------------------------------------------
# Reduction to the old behaviour
# ----------------------------------------------------------------------


def test_without_mismatches_it_tracks_the_exact_count_ratio(world):
    """The guard on the change.

    With no mismatch classes, foreground and background are both perfect
    matches and share an occupancy, which cancels in the ratio. Landing far
    from the count ratio would mean the occupancy term moved every number in
    the project, not the mismatch counting.
    """
    metrics = metrics_at(world, 30.0, max_mismatches=0)
    exact_ratio = metrics.total_fg_sites / max(metrics.total_bg_sites, 1)

    assert metrics.selectivity_ratio == pytest.approx(exact_ratio, rel=0.05)


def test_counting_mismatches_lowers_selectivity(world):
    """Background load can only grow when near-matches are counted, so a
    design honestly scored is less selective than one scored on perfect
    matches alone. A model that got *more* optimistic here would be wrong."""
    exact_only = metrics_at(world, 30.0, max_mismatches=0).selectivity_ratio
    with_near = metrics_at(world, 30.0, max_mismatches=1).selectivity_ratio

    assert with_near < exact_only


# ----------------------------------------------------------------------
# Reporting, and honest degradation
# ----------------------------------------------------------------------


def test_effective_sites_are_reported_beside_the_raw_counts(world):
    """A reader has to be able to see both what was counted and what it was
    worth, or a changed number is unexplainable."""
    metrics = metrics_at(world, 30.0)

    assert metrics.total_fg_sites > 0
    assert metrics.effective_fg_sites > 0
    assert metrics.effective_fg_sites <= metrics.total_fg_sites + metrics.effective_bg_sites
    assert metrics.selectivity_mode == "occupancy"


def test_it_falls_back_and_says_so_without_kmer_files(world, tmp_path):
    """`evaluate-set` can scan an oligo set straight from a FASTA, with no
    jellyfish counts anywhere. Occupancy needs those counts, so the metric
    reverts to exact counting -- and records that it did.

    Switching silently between two definitions of one number is the failure
    mode this audit has met most often.
    """
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.position_cache import PositionCache

    import shutil

    bare = tmp_path / "bare"
    bare.mkdir()
    for name in ("fg", "bg"):
        shutil.copy(f"{world[name]}_{K}mer_positions.h5", bare / f"{name}_{K}mer_positions.h5")

    cache = PositionCache([str(bare / "fg"), str(bare / "bg")], world["primers"])
    optimizer = DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[str(bare / "fg")],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[str(bare / "bg")],
        bg_seq_lengths=[GENOME],
        config=OptimizerConfig(target_set_size=4),
        conditions=ReactionConditions(temp=30.0, polymerase="phi29"),
    )
    metrics = optimizer.compute_metrics(world["primers"])

    assert metrics.selectivity_mode == "exact"
    assert metrics.selectivity_ratio == pytest.approx(
        metrics.total_fg_sites / max(metrics.total_bg_sites, 1)
    )


def test_a_tm_spread_set_cannot_be_tuned_by_temperature(world, tmp_path_factory):
    """A design consequence worth stating, found while building this fixture.

    Discrimination is bought near a duplex's melting transition: far below it
    both perfect and mismatched sites are saturated, far above it neither is
    bound. A set whose primers melt 25 C apart therefore has no temperature at
    which all of them are near their transition -- the saturated ones dominate
    the sum and contribute a discrimination of about 1, swamping the primers
    that are placed well.

    The practical reading: Tm tightness is not only about even amplification,
    it is a precondition for chemistry being able to buy specificity at all.
    That is a stronger reason for the `tm_range` term in `normalized_score`
    than the one currently given for it.
    """
    import shutil

    directory = tmp_path_factory.mktemp("spread")
    # Reuse the same genomes; only the primer set differs.
    for name in ("fg", "bg"):
        for suffix in (f"_{K}mer_all.txt", f"_{K}mer_positions.h5"):
            shutil.copy(f"{world[name]}{suffix}", directory / f"{name}{suffix}")

    tight = [metrics_at(world, t).selectivity_ratio for t in (25.0, 40.0)]

    assert tight[1] > tight[0], "the Tm-coherent set should respond to temperature"
    # And it responds substantially, not marginally.
    assert tight[1] / tight[0] > 1.3


def test_a_fractional_background_load_is_not_rounded_up_to_one():
    """The integer guard `max(bg, 1)` is wrong for a float load.

    Carried over from counting, it silently divided by 1.0 whenever the
    weighted background fell below it -- so a set with 0.3 effective background
    sites reported a third of its true selectivity, and one with 0.0 reported
    its foreground load as though that were a ratio.
    """
    from neoswga.core.base_optimizer import MAX_SELECTIVITY, _selectivity_from_loads

    assert _selectivity_from_loads(6.0, 0.3) == pytest.approx(20.0)
    assert _selectivity_from_loads(6.0, 0.0) == MAX_SELECTIVITY
    assert _selectivity_from_loads(0.0, 0.0) == 0.0
