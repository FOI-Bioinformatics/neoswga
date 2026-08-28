"""Asking for more primers must not produce a worse design.

The hybrid optimizer runs set cover (Stage 1) and then refines by removing
primers to maximise amplification-network connectivity (Stage 2). Stage 2's
removal criterion was `connectivity + pred_amp / 100` -- coverage appeared
nowhere in it.

That is not a neutral omission. The primer whose removal leaves the best
network score is the one that is most peripheral in the amplification graph,
and a primer is peripheral precisely when its binding sites are far from
everything else -- which is to say, when it is the only thing covering its
region. Stage 2 was therefore biased towards removing exactly the primers
carrying unique coverage. Measured on a 100 kb target, it took a Stage-1 set
covering 94.0% down to 64.9%.

The user-visible consequence was that coverage stopped being monotone in
`--num-primers`. Stage-1 count is sized from the requested count, so a larger
request builds a larger Stage-1 pool and hands Stage 2 more to remove:

    asked  stage1_n  stage1_cov   final_n  final_cov
        8        16       0.910         8      0.687
       10        18       0.940        10      0.649   <- more primers, less coverage

Stage-1 coverage rises monotonically across that table. All of the loss is in
Stage 2.

Two changes are pinned here. The removal criterion now weighs the coverage a
removal costs, so a redundant primer is dropped in preference to a uniquely
covering one. And the result is floored at the Stage-1 greedy prefix of the
same size: greedy set cover selects in a deterministic order, so its first N
primers form a chain whose coverage is non-decreasing in N by construction,
which makes it a monotone baseline the refinement is never allowed to fall
below.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("networkx")

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 100_000
PRIMER_LENGTH = 10
COVERAGE_REACH = 3_000
N_CANDIDATES = 40


@pytest.fixture(scope="module")
def planted():
    """A target large enough that coverage does not saturate.

    The bundled 6 kb plasmid cannot exercise this at all: a 3 kb reach covers
    it with two primers, so every request saturates and the count never binds.
    """
    rng = random.Random(4242)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    for i in range(N_CANDIDATES):
        primer = "".join(rng.choice("ACGT") for _ in range(PRIMER_LENGTH))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(3):
            pos = (i * (GENOME_LENGTH // N_CANDIDATES) + j * 900) % (GENOME_LENGTH - 20)
            seq[pos : pos + PRIMER_LENGTH] = list(primer)

    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture(scope="module")
def optimizer(tmp_path_factory, planted):
    prefix = str(tmp_path_factory.mktemp("stage2") / "target")
    keys = set()
    for primer in planted["primers"]:
        keys.update({primer, reverse_complement(primer)})

    with h5py.File(f"{prefix}_{PRIMER_LENGTH}mer_positions.h5", "w") as f:
        for key in sorted(keys):
            positions, i = [], planted["seq"].find(key)
            while i != -1:
                positions.append(i)
                i = planted["seq"].find(key, i + 1)
            if positions:
                f.create_dataset(key, data=np.array(positions, dtype=np.int32))

    cache = PositionCache([prefix], planted["primers"])
    return HybridOptimizer(
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        bin_size=10_000,
        coverage_reach=COVERAGE_REACH,
    )


REQUESTS = [4, 6, 8, 10, 12]


@pytest.fixture(scope="module")
def results(optimizer, planted):
    return {
        n: optimizer.optimize(planted["primers"], final_count=n, verbose=False) for n in REQUESTS
    }


# ----------------------------------------------------------------------
# The property
# ----------------------------------------------------------------------


def test_coverage_is_monotone_in_the_requested_count(results):
    """More primers must not mean less genome covered.

    This is the behaviour a user reasons about when choosing a set size, and
    it was violated: 10 primers covered less than 8.
    """
    table = [(n, len(results[n].primers), round(results[n].final_coverage, 4)) for n in REQUESTS]

    for (n_small, _, cov_small), (n_large, _, cov_large) in zip(table, table[1:]):
        assert cov_large >= cov_small - 1e-9, (
            f"asking for {n_large} primers covered {cov_large} where "
            f"{n_small} covered {cov_small}; full table: {table}"
        )


def test_the_requested_count_is_still_honoured(results):
    """The fix must not buy monotonicity by returning a different size."""
    for n in REQUESTS:
        assert len(results[n].primers) == n


# ----------------------------------------------------------------------
# The floor that makes it a guarantee rather than a tendency
# ----------------------------------------------------------------------


def test_refinement_never_falls_below_the_stage1_prefix(optimizer, planted, results):
    """Greedy set cover selects in a deterministic order, so its first N
    primers are the same N whatever `stage1_count` was. That prefix is a
    coverage-monotone baseline, and the refined set is never allowed to be
    worse than it.

    Without the floor, "Stage 2 usually keeps coverage" is a tendency that
    holds on the fixtures someone happened to test.
    """
    for n in REQUESTS:
        result = results[n]
        prefix = list(result.stage1_ordered_primers)[:n]
        assert len(prefix) == n, "stage-1 order was not recorded"

        floor = optimizer._calculate_coverage(prefix)
        assert result.final_coverage >= floor - 1e-9, (
            f"refined set for n={n} covers {result.final_coverage:.4f}, below the "
            f"stage-1 greedy prefix at {floor:.4f}"
        )


def test_stage1_order_is_a_prefix_chain(optimizer, planted):
    """The property the floor rests on: the greedy order does not depend on
    how many primers were requested, so prefixes nest."""
    small = optimizer.optimize(planted["primers"], final_count=4, verbose=False)
    large = optimizer.optimize(planted["primers"], final_count=12, verbose=False)

    a = list(small.stage1_ordered_primers)
    b = list(large.stage1_ordered_primers)
    shared = min(len(a), len(b))

    assert a[:shared] == b[:shared], "greedy selection order changed with the request size"


def test_stage1_prefix_coverage_is_itself_monotone(optimizer, planted):
    """Guards the baseline. If the prefix chain were not monotone, flooring
    against it would not deliver monotone results."""
    result = optimizer.optimize(planted["primers"], final_count=12, verbose=False)
    order = list(result.stage1_ordered_primers)

    covs = [optimizer._calculate_coverage(order[:n]) for n in REQUESTS if n <= len(order)]
    for a, b in zip(covs, covs[1:]):
        assert b >= a - 1e-9, f"greedy prefix coverage fell: {covs}"


# ----------------------------------------------------------------------
# Stage 2 must still do its job
# ----------------------------------------------------------------------


def test_refinement_still_improves_network_connectivity(optimizer, planted):
    """Stage 2 exists to pick a connected, amplifiable subset. Weighting
    coverage into the removal must not turn it into a second set-cover pass
    that ignores the network entirely.
    """
    result = optimizer.optimize(planted["primers"], final_count=8, verbose=False)

    refined = optimizer._build_network(list(result.primers)).get_statistics()
    prefix = optimizer._build_network(list(result.stage1_ordered_primers)[:8]).get_statistics()

    assert refined["largest_component"] >= prefix["largest_component"], (
        "refinement no longer improves on the coverage-only prefix; "
        f"refined={refined['largest_component']} prefix={prefix['largest_component']}"
    )


def test_refinement_returns_a_subset_of_stage_one(results):
    for n in REQUESTS:
        assert set(results[n].primers) <= set(results[n].stage1_primers)
