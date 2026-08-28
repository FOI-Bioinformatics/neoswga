"""Two reported numbers that did not mean what they said.

**`mean_tm` was the Wallace rule.** `4*GC + 2*AT` is a 1979 approximation for
14-20mers at 1 M salt. The optimizer reported it for equiphi29 12-mer designs
run at 42 C in a 10 mM Mg buffer, where it reads about 11 C low: a set whose
salt-corrected nearest-neighbour Tm is 45.9 C was published as 34.9 C. That
inverts the one check an equiphi29 design exists to pass -- do the primers melt
above the reaction temperature -- and it disagreed with the Tm the occupancy
model was using at the same moment, which is the real one.

**`selectivity_ratio` is a ratio of counts, not of densities.** It carries no
genome length, so it moves with how much background sequence the user happened
to supply. Swapping a single human chromosome (46.7 Mb) for the whole genome
(3.1 Gb) drops it about 66x with nothing about the primers changed, and the
number that governs enrichment -- sites per base of target against sites per
base of host -- is not recoverable from it without knowing both lengths.

The count ratio is kept: it is what the objective is scored on, and within one
run the background size is fixed, so selection is unaffected. The density ratio
is reported beside it.
"""

import pytest

from neoswga.core.base_optimizer import (
    MAX_SELECTIVITY,
    OptimizerConfig,
    _selectivity_density_from_loads,
)
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import calculate_tm_basic, wallace_tm

# A real 12-mer from the Prevotella/equiphi29 validation set. Wallace puts it at
# 34 C; under the configured buffer it melts at about 47 C.
PRIMER = "AAAAGGGCGTTA"


def _optimizer(conditions=None):
    """Bare optimizer: `_estimate_tm` reads no position data."""
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter

    return DominatingSetAdapter(
        position_cache=None,
        fg_prefixes=["fg"],
        fg_seq_lengths=[1_000_000],
        config=OptimizerConfig(target_set_size=4),
        conditions=conditions,
    )


# ----------------------------------------------------------------------
# mean_tm
# ----------------------------------------------------------------------


def test_estimate_tm_uses_the_configured_reaction_conditions():
    """The reported Tm must be the one the occupancy model uses."""
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29", mg_conc=10.0)
    optimizer = _optimizer(conditions)

    assert optimizer._estimate_tm(PRIMER) == pytest.approx(
        conditions.calculate_effective_tm(PRIMER), abs=1e-6
    )


def test_estimate_tm_is_not_the_wallace_rule():
    """Pinned separately: the two agree only by coincidence, never here."""
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29", mg_conc=10.0)

    reported = _optimizer(conditions)._estimate_tm(PRIMER)

    assert reported > wallace_tm(PRIMER) + 5.0
    # The design question this number exists to answer.
    assert reported > conditions.temp


def test_estimate_tm_without_conditions_falls_back_to_nearest_neighbour():
    """No buffer to correct for is not a reason to return a worse estimate."""
    optimizer = _optimizer(conditions=None)

    assert optimizer._estimate_tm(PRIMER) == pytest.approx(
        calculate_tm_basic(PRIMER), abs=1e-6
    )


# ----------------------------------------------------------------------
# selectivity density
# ----------------------------------------------------------------------


def test_density_ratio_is_invariant_to_background_size():
    """The point of the metric.

    A background twice as long carries about twice the binding load, so the
    count ratio halves while nothing about the primers changed. Sites per base
    is the quantity that survives the substitution.
    """
    small = _selectivity_density_from_loads(1000.0, 3_000_000, 100.0, 46_000_000)
    large = _selectivity_density_from_loads(1000.0, 3_000_000, 6_700.0, 3_100_000_000)

    assert small == pytest.approx(large, rel=0.02)


def test_density_ratio_normalises_by_both_lengths():
    fg_per_base = 1000.0 / 3_000_000
    bg_per_base = 100.0 / 46_000_000

    assert _selectivity_density_from_loads(
        1000.0, 3_000_000, 100.0, 46_000_000
    ) == pytest.approx(fg_per_base / bg_per_base)


def test_density_ratio_reports_max_when_no_background_binds():
    """Same convention as the count ratio: unbounded, not large."""
    assert (
        _selectivity_density_from_loads(1000.0, 3_000_000, 0.0, 46_000_000)
        == MAX_SELECTIVITY
    )


def test_density_ratio_is_zero_without_a_measurable_background():
    """A background of no length is an absent measurement, not a clean one."""
    assert _selectivity_density_from_loads(1000.0, 3_000_000, 5.0, 0) == 0.0


# ----------------------------------------------------------------------
# reported beside the count ratio
# ----------------------------------------------------------------------

GENOME_FG = 40_000
# Deliberately unequal: with equal lengths the two ratios coincide and the
# test could not tell them apart.
GENOME_BG = 120_000
K = 10


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """Foreground and background with k-mer counts and binding positions."""
    import random

    h5py = pytest.importorskip("h5py")
    import numpy as np

    from neoswga.core.mismatch_counts import canonical_kmer
    from neoswga.core.thermodynamics import reverse_complement

    rng = random.Random(7)
    directory = tmp_path_factory.mktemp("density")

    fg_seq = "".join(rng.choice("ACGT") for _ in range(GENOME_FG))
    bg_seq = "".join(rng.choice("ACGT") for _ in range(GENOME_BG))
    primers = sorted({fg_seq[i : i + K] for i in range(0, 4000, 137)})[:6]

    prefixes = {}
    for name, seq in (("fg", fg_seq), ("bg", bg_seq)):
        prefix = str(directory / name)
        prefixes[name] = prefix

        counts = {}
        for i in range(len(seq) - K + 1):
            key = canonical_kmer(seq[i : i + K])
            counts[key] = counts.get(key, 0) + 1
        with open(f"{prefix}_{K}mer_all.txt", "w") as handle:
            for key, value in counts.items():
                handle.write(f"{key} {value}\n")

        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as handle:
            for primer in primers:
                for key in {primer, reverse_complement(primer)}:
                    hits, at = [], seq.find(key)
                    while at != -1:
                        hits.append(at)
                        at = seq.find(key, at + 1)
                    if hits:
                        handle.create_dataset(key, data=np.array(hits, dtype=np.int32))

    return {"primers": primers, "fg": prefixes["fg"], "bg": prefixes["bg"]}


def _metrics(world):
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.position_cache import PositionCache

    cache = PositionCache([world["fg"], world["bg"]], world["primers"])
    optimizer = DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[world["fg"]],
        fg_seq_lengths=[GENOME_FG],
        bg_prefixes=[world["bg"]],
        bg_seq_lengths=[GENOME_BG],
        config=OptimizerConfig(target_set_size=4, max_mismatches=1),
        conditions=ReactionConditions(temp=30.0, polymerase="phi29", mg_conc=10.0),
    )
    return optimizer.compute_metrics(world["primers"])


def test_metrics_report_the_density_ratio(world):
    metrics = _metrics(world)

    expected = _selectivity_density_from_loads(
        metrics.effective_fg_sites, GENOME_FG, metrics.effective_bg_sites, GENOME_BG
    )
    assert metrics.selectivity_density == pytest.approx(expected)


def test_density_and_count_ratios_differ_when_the_genomes_do(world):
    """The whole reason for reporting both."""
    metrics = _metrics(world)

    assert metrics.selectivity_density == pytest.approx(
        metrics.selectivity_ratio * (GENOME_BG / GENOME_FG), rel=1e-6
    )


def test_to_dict_carries_the_density_ratio_and_its_lengths(world):
    """A density is unreadable without the lengths it was taken over."""
    payload = _metrics(world).to_dict()

    assert payload["selectivity_density"] == pytest.approx(
        _metrics(world).selectivity_density
    )
    assert payload["fg_total_length"] == GENOME_FG
    assert payload["bg_total_length"] == GENOME_BG
