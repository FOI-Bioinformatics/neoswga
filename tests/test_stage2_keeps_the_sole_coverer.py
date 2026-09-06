"""Stage-2 refinement must not drop the primer that solely covers a region.

`_network_refine` trades network connectivity against coverage cost, and the
coverage half of that trade is "how many bins would we lose that nothing else
covers". It counted every bin of every primer as uniquely covered, so the term
it computed was the primer's TOTAL bin count, and the ranking -- which reads a
low count as cheap to remove -- preferred discarding the primer with the fewest
bins. That is typically the one primer holding a region on its own.

The cause is a dataclass subtlety rather than the arithmetic.
`CoverageRegion` defines `__hash__` over `(chromosome, start, end)` but inherits
the generated `__eq__`, which compares all four fields including `covered_by`:

    hash equal : True
    eq         : False
    set size   : 2

`_coverage_bins_by_primer` built a SEPARATE `BipartiteGraph` per primer, so the
same genomic bin reached the occupancy counter once per primer, as a distinct
key each time, and every count came back 1. Measured on eight primers from the
plasmid example before the fix: 50 per-primer bins, 50 distinct counter keys,
zero keys with a count above 1.

The comment above the ranking records that this term exists to prevent a 94.0
to 64.9 per cent coverage regression. It was selecting for it.

The fix reuses what `BipartiteGraph` already does correctly for the greedy
path: one graph over all the primers, `primer_to_regions` for a primer's bins
and `region_to_primers` for a bin's occupancy. That is why greedy never had this
bug and stage 2 did.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("networkx")

import numpy as np

from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 300_000
PRIMER_LENGTH = 10

#: Primers 0-5 bind the same six places in the left half, so each is redundant:
#: removing any one of them costs no coverage at all. The last primer is the
#: only thing binding the right end, so removing it costs a whole region.
REDUNDANT_COUNT = 6
SHARED_SITES = [1_000, 6_000, 11_000, 16_000, 21_000, 26_000]
SOLE_SITES = [250_000, 262_000, 274_000]


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """A genome where one primer is load-bearing and the rest are not."""
    rng = random.Random(20260903)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    while len(primers) < REDUNDANT_COUNT + 1:
        candidate = "".join(rng.choice("ACGT") for _ in range(PRIMER_LENGTH))
        if candidate not in primers and reverse_complement(candidate) not in primers:
            primers.append(candidate)

    sole = primers[-1]
    placements = {p: list(SHARED_SITES) for p in primers[:REDUNDANT_COUNT]}
    placements[sole] = list(SOLE_SITES)

    prefix = str(tmp_path_factory.mktemp("stage2") / "target")
    with h5py.File(f"{prefix}_{PRIMER_LENGTH}mer_positions.h5", "w") as handle:
        for primer, sites in placements.items():
            handle.create_dataset(primer, data=np.array(sites, dtype=np.int64))

    return {
        "prefix": prefix,
        "primers": primers,
        "sole": sole,
        "cache": PositionCache([prefix], primers),
    }


def make_optimizer(world):
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    return HybridOptimizer(
        position_cache=world["cache"],
        fg_prefixes=[world["prefix"]],
        fg_seq_lengths=[GENOME_LENGTH],
        bin_size=1_000,
        coverage_reach=3_000,
    )


# ---------------------------------------------------------------------------
# The occupancy counter
# ---------------------------------------------------------------------------


def test_a_bin_covered_by_several_primers_is_counted_once(world):
    """Six primers binding the same sites must share their bins, not own them.

    This is the measurement the defect showed up in, asked directly: if every
    bin comes back with an occupancy of 1, the "uniquely covered" term is really
    a total-bins term and the ranking built on it is inverted.
    """
    optimizer = make_optimizer(world)
    redundant = world["primers"][:REDUNDANT_COUNT]

    bins_by_primer = optimizer._coverage_bins_by_primer(redundant)
    occupancy = optimizer._bin_occupancy(bins_by_primer)

    assert bins_by_primer[redundant[0]], "the first primer covers no bins at all"
    shared = [key for key, count in occupancy.items() if count > 1]
    assert shared, (
        "no bin is recorded as covered by more than one primer, although all "
        f"{REDUNDANT_COUNT} bind identical positions. The occupancy keys are not "
        "being shared, so every bin reads as uniquely covered."
    )


def test_the_sole_coverer_owns_bins_that_no_other_primer_owns(world):
    """The load-bearing primer must be the one with uniquely covered bins."""
    optimizer = make_optimizer(world)
    primers = world["primers"]

    bins_by_primer = optimizer._coverage_bins_by_primer(primers)
    occupancy = optimizer._bin_occupancy(bins_by_primer)

    def unique_bins(primer):
        return sum(
            1 for b in bins_by_primer.get(primer, ()) if occupancy[optimizer._bin_key(b)] == 1
        )

    sole_unique = unique_bins(world["sole"])
    redundant_unique = [unique_bins(p) for p in primers[:REDUNDANT_COUNT]]

    assert sole_unique > 0, "the sole coverer owns no unique bins"
    assert all(count == 0 for count in redundant_unique), (
        "a primer whose every site is shared with five others still reports "
        f"uniquely covered bins: {redundant_unique}"
    )


# ---------------------------------------------------------------------------
# The consequence for a delivered pool
# ---------------------------------------------------------------------------


def test_refinement_keeps_the_primer_that_solely_covers_a_region(world):
    """The whole point of the coverage term in the stage-2 ranking.

    Asking for five of seven forces two removals. Both should come from the six
    interchangeable primers, never from the one holding the right end alone.
    """
    optimizer = make_optimizer(world)
    kept = optimizer._network_refine(list(world["primers"]), target_count=5, verbose=False)

    assert len(kept) == 5, f"expected five primers, got {len(kept)}: {kept}"
    assert world["sole"] in kept, (
        f"stage 2 dropped {world['sole']}, the only primer covering "
        f"{SOLE_SITES[0]}-{SOLE_SITES[-1]}, and kept six interchangeable ones: {kept}"
    )
