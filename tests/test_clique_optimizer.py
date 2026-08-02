"""The clique optimizer, and the guarantee that justifies it.

The four optimizers this tool ships all *penalise* primer-dimers: a dimerising
pair costs score, and a bad enough pair usually loses. None of them can promise
a set contains no dimerising pair, because a greedy search trades that cost off
against coverage and can always decide the coverage was worth it.

swga 1.0 (Clarke 2017) makes the promise structurally instead. Model primers as
vertices and join two with an edge when they do NOT dimerise; a clique in that
graph is by construction a set of mutually compatible primers, so the guarantee
falls out of the representation rather than out of a threshold. A dimerising
pair consumes primer that should be priming the template and gives the
polymerase a substrate that is not the target, so this is worth having as an
option even though it costs an exhaustive search.

`clique_optimizer.py` was deleted in commit 08044d9 for being unused, not for
being wrong. These tests cover the restoration, and two defects that had to be
fixed on the way back in: an enumeration cap that could not stop a single large
clique from expanding to millions of subsets, and a selection objective that
disagreed with the metric the result was reported and status-flagged on.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")
nx = pytest.importorskip("networkx")

import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.clique_optimizer import (
    CliqueOptimizer,
    CliqueOptimizerConfig,
    build_compatibility_graph,
    enumerate_dimer_free_sets,
)
from neoswga.core.dimer import is_dimer_fast
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 40_000
PRIMER_LENGTH = 10
MAX_DIMER_BP = 3


@pytest.fixture(scope="module")
def planted():
    """A genome with primers planted in it, including deliberate dimer pairs.

    Half the pool is seeded as reverse complements of the other half. Those
    pairs are maximally complementary, so any set containing both members is
    exactly what a dimer-free guarantee has to exclude -- and a greedy
    optimizer weighing a dimer penalty against coverage might not.
    """
    rng = random.Random(31337)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    for i in range(10):
        primer = "".join(rng.choice("ACGT") for _ in range(PRIMER_LENGTH))
        if primer in primers:
            continue
        primers.append(primer)
        # Its reverse complement: a guaranteed dimer partner.
        rc = reverse_complement(primer)
        if rc not in primers:
            primers.append(rc)

    for i, primer in enumerate(primers):
        for j in range(3):
            pos = (i * 1_700 + j * 550) % (GENOME_LENGTH - 20)
            seq[pos : pos + PRIMER_LENGTH] = list(primer)

    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture(scope="module")
def optimizer(tmp_path_factory, planted):
    prefix = str(tmp_path_factory.mktemp("clique") / "target")
    # The pool deliberately contains reverse-complement pairs, so a primer and
    # its partner ask for the same two HDF5 keys. Collect them once.
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
    return CliqueOptimizer(
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        config=CliqueOptimizerConfig(target_set_size=4, max_dimer_bp=MAX_DIMER_BP),
    )


# ----------------------------------------------------------------------
# The guarantee
# ----------------------------------------------------------------------


def test_the_returned_set_contains_no_dimerising_pair(optimizer, planted):
    """The whole reason this optimizer exists.

    The pool is seeded with reverse-complement pairs, so a set that ignored
    compatibility would very likely contain one.
    """
    result = optimizer.optimize(planted["primers"], target_size=4)
    assert result.primers, f"optimizer returned nothing: {result.message}"

    chosen = list(result.primers)
    for i in range(len(chosen)):
        for j in range(i + 1, len(chosen)):
            assert not is_dimer_fast(
                chosen[i], chosen[j], MAX_DIMER_BP
            ), f"{chosen[i]} and {chosen[j]} dimerise, which a clique cannot contain"


def test_reverse_complement_pairs_are_never_both_selected(optimizer, planted):
    """The sharpest case: a primer and its exact reverse complement."""
    chosen = set(optimizer.optimize(planted["primers"], target_size=4).primers)

    for primer in chosen:
        assert reverse_complement(primer) not in chosen


def test_compatibility_graph_never_joins_a_dimerising_pair(planted):
    """The guarantee is a property of the graph, so test it there too --
    everything downstream inherits correctness from this one invariant."""
    G = build_compatibility_graph(planted["primers"], max_dimer_bp=MAX_DIMER_BP)

    for a, b in G.edges():
        assert not is_dimer_fast(a, b, MAX_DIMER_BP)


def test_every_enumerated_set_is_a_clique(planted):
    G = build_compatibility_graph(planted["primers"], max_dimer_bp=MAX_DIMER_BP)

    for primer_set in enumerate_dimer_free_sets(G, target_size=3, max_cliques=50):
        for i in range(len(primer_set)):
            for j in range(i + 1, len(primer_set)):
                assert G.has_edge(primer_set[i], primer_set[j])


# ----------------------------------------------------------------------
# The enumeration cap has to actually cap
# ----------------------------------------------------------------------


def test_enumeration_respects_the_cap_inside_a_single_large_clique():
    """The defect found on restoration.

    The cap was tested only after each clique finished expanding, so one large
    clique ran `combinations` to exhaustion first. A complete graph on 40
    primers has a single clique of 40; at target size 6 that is 3.8 million
    subsets, every one appended to a list and a set before anything checked the
    limit. A cap meant to bound the work instead bounded only how many cliques
    got to blow up.

    A complete graph is the worst case and the right test: there is exactly one
    maximal clique, so nothing but the inner check can stop the expansion.
    """
    G = nx.complete_graph(["P%02d" % i for i in range(40)])

    sets = enumerate_dimer_free_sets(G, target_size=6, max_cliques=25)

    assert len(sets) == 25


def test_enumeration_reports_when_it_was_not_exhaustive(caplog):
    """Silent truncation reads as "searched everything", which for a method
    whose selling point is exhaustiveness is the wrong thing to imply."""
    import logging

    G = nx.complete_graph(["P%02d" % i for i in range(20)])

    with caplog.at_level(logging.INFO):
        enumerate_dimer_free_sets(G, target_size=4, max_cliques=10)

    assert "NOT exhaustive" in caplog.text


# ----------------------------------------------------------------------
# Selection agrees with what gets reported
# ----------------------------------------------------------------------


def test_result_is_scored_on_the_comparable_metric(optimizer, planted):
    """`normalized_score` is the only value comparable across optimizers, and
    it is what an ensemble run picks a winner by.

    Selection used a raw fg/bg site ratio while the result's PARTIAL status was
    decided on `fg_coverage`, so the set chosen was not the best set by the
    measure it was then judged on.
    """
    result = optimizer.optimize(planted["primers"], target_size=4)

    assert result.score == pytest.approx(result.metrics.normalized_score())
    assert 0.0 <= result.score <= 1.0


def test_reported_metrics_describe_the_returned_set(optimizer, planted):
    """A mismatch here means the metrics belong to some other candidate set."""
    result = optimizer.optimize(planted["primers"], target_size=4)

    recomputed = optimizer.compute_metrics(list(result.primers))
    assert result.metrics.fg_coverage == pytest.approx(recomputed.fg_coverage)


# ----------------------------------------------------------------------
# Contracts shared with the other optimizers
# ----------------------------------------------------------------------


def test_it_selects_only_from_the_candidates(optimizer, planted):
    result = optimizer.optimize(planted["primers"], target_size=4)
    assert set(result.primers) <= set(planted["primers"])


def test_it_returns_the_requested_size(optimizer, planted):
    result = optimizer.optimize(planted["primers"], target_size=3)
    assert len(result.primers) == 3


def test_fixed_primers_are_kept(optimizer, planted):
    """The iterative-design workflow depends on validated oligos surviving."""
    graph = build_compatibility_graph(planted["primers"], max_dimer_bp=MAX_DIMER_BP)
    compatible_pair = next(iter(graph.edges()))

    result = optimizer.optimize(
        planted["primers"], target_size=4, fixed_primers=list(compatible_pair)
    )

    assert set(compatible_pair) <= set(result.primers)


def test_mutually_incompatible_fixed_primers_are_refused(optimizer, planted):
    """Silently dropping one would hand back a set the user did not ask for."""
    primer = planted["primers"][0]
    result = optimizer.optimize(
        planted["primers"], target_size=4, fixed_primers=[primer, reverse_complement(primer)]
    )

    assert not result.primers
    assert "dimer" in result.message.lower()


def test_an_impossible_target_size_fails_rather_than_inventing_a_set(optimizer, planted):
    result = optimizer.optimize(planted["primers"][:2], target_size=8)
    assert not result.primers


def test_registered_under_its_method_name():
    """`--optimization-method clique` has to resolve."""
    import neoswga.core.clique_optimizer  # noqa: F401  (registration side effect)
    from neoswga.core.optimizer_factory import OptimizerFactory

    assert "clique" in [n.lower() for n in OptimizerFactory.list_optimizers()]
