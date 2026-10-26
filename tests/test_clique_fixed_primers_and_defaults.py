"""Two small clique defects: an off-by-one on fixed primers, and split defaults.

`optimize` reduced its search target by the number of fixed primers with
`max(1, target - len(fixed))`. The floor of 1 meant that fixing as many primers
as were asked for still requested one more from the search, so `--num-primers 4`
with four fixed primers returned five. Zero is the correct floor: fixing the
whole set is a complete answer.

Separately, `build_compatibility_graph`'s own defaults were `max_dimer_bp=3,
max_self_dimer_bp=4` while `OptimizerConfig`'s are 4 and 5. A direct caller of
the helper therefore got a STRICTER graph than the optimizer built from the same
unset configuration -- two answers to "what counts as a dimer here", declared a
few lines apart in the same file.
"""

import dataclasses
import inspect

import pytest

pytest.importorskip("networkx")

import itertools

import numpy as np

from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.clique_optimizer import CliqueOptimizer, build_compatibility_graph
from neoswga.core.dimer import is_dimer_fast


class _StubCache:
    """Positions are irrelevant to the dimer guarantee, so stub them.

    `_score_set` is overridden in each test to steer selection, so the cache is
    never consulted for coverage; this keeps the guarantee tests free of a
    genome.
    """

    def __init__(self, *args, **kwargs):
        pass

    def get_positions(self, *args, **kwargs):
        return np.array([], dtype=np.int64)


def _clique():
    return CliqueOptimizer(position_cache=_StubCache(), fg_prefixes=["x"], fg_seq_lengths=[1000])


def test_a_candidate_that_dimers_with_a_fixed_primer_is_never_delivered():
    """The clique guarantee must hold between fixed primers and the clique.

    Fixed primers are removed from the candidate pool before the graph is
    built, so the old `if fp in G` filter was always False and never excluded a
    candidate that heterodimers with a fixed primer. A clique drawn from the
    unfiltered pool could then be returned as `fixed + clique`, a dimerising
    pair the whole method exists to forbid (proved at the filter step on
    2026-10-26). The filter now screens candidates directly against each fixed
    primer.
    """
    fixed = "AAAAGGGGAAAA"  # does not self-dimer
    partner = "TTTTCCCCTTTT"  # heterodimers with `fixed`
    benign = ["ACACACATACAC", "ATACACAACAAT", "AATACAATACAA", "ACAATACAAATA", "TACAAATAACAA"]
    assert not is_dimer_fast(fixed, fixed, 3)
    assert is_dimer_fast(fixed, partner, 3)

    opt = _clique()
    # Make the dimerising partner the most attractive candidate, so only the
    # filter can keep it out of the delivered set.
    opt._score_set = lambda s: 1.0 if partner in s else 0.0
    result = opt.optimize(candidates=benign + [partner], num_primers=4, fixed_primers=[fixed])

    delivered = list(result.primers)
    pairs = list(itertools.combinations(delivered, 2))
    dimers = [(a, b) for a, b in pairs if is_dimer_fast(a, b, 3)]
    assert not dimers, f"the delivered set contains a dimerising pair: {dimers}"
    assert partner not in delivered, "the candidate that dimers with the fixed primer was delivered"


def test_a_self_dimering_fixed_primer_is_refused():
    """A fixed primer that self-dimers cannot be in a dimer-free set.

    The method held candidates to this (it drops self-dimering nodes) and
    fixed-fixed pairs to it (the mutual check), but let a self-dimering fixed
    primer through with a warning and delivered it anyway. It now refuses, the
    same answer the mutual-dimer check already gives.
    """
    self_dimer = "GGGGCCCCGGGG"
    assert is_dimer_fast(self_dimer, self_dimer, 4)
    benign = ["ACACACATACAC", "ATACACAACAAT", "AATACAATACAA", "ACAATACAAATA"]

    result = _clique().optimize(candidates=benign, num_primers=4, fixed_primers=[self_dimer])
    assert result.status == OptimizationStatus.ERROR
    assert "self-dimer" in (result.message or "").lower()


@pytest.mark.parametrize("field", ["max_dimer_bp", "max_self_dimer_bp"])
def test_the_helper_default_matches_the_config_default(field):
    """One answer to "what counts as a dimer", not two."""
    helper = inspect.signature(build_compatibility_graph).parameters[field].default
    config = {f.name: f.default for f in dataclasses.fields(OptimizerConfig)}[field]

    assert helper == config, (
        f"build_compatibility_graph defaults {field}={helper} while "
        f"OptimizerConfig defaults it to {config}; a direct caller and the "
        "optimizer disagree about the threshold"
    )


def _target_after_fixing(target, n_fixed):
    """The production expression, reached through the optimizer.

    Kept as a direct check of the arithmetic because reaching the branch
    end-to-end needs a genome, a cache and a pool that admits a clique, none of
    which this off-by-one depends on.
    """
    from neoswga.core import clique_optimizer

    source = inspect.getsource(clique_optimizer.CliqueOptimizer.optimize)
    assert "max(0, target - len(fixed))" in source, (
        "the fixed-primer reduction is not floored at 0; `max(1, ...)` asks the "
        "search for one more primer than was requested once the fixed set fills "
        "the target"
    )
    return max(0, target - n_fixed)


@pytest.mark.parametrize(
    "target,n_fixed,expected",
    [(6, 0, 6), (6, 2, 4), (4, 4, 0), (4, 5, 0)],
)
def test_fixing_the_whole_set_asks_the_search_for_nothing_more(target, n_fixed, expected):
    assert _target_after_fixing(target, n_fixed) == expected


# ----------------------------------------------------------------------
# The free-energy floor (`max_dimer_dg`), which clique read nowhere
# ----------------------------------------------------------------------
#
# `clique` is documented as "the only method that guarantees no dimerising
# pair", and until 2026-10-26 it read `max_dimer_dg` nowhere at all: zero
# occurrences in the module, and no such field on `CliqueOptimizerConfig`. A
# user who configured a free-energy floor and chose this method got base-pair
# screening alone, silently -- Known Issue 8, in the shape where the option is
# declared on the base config and one method never looks at it.
#
# These fixtures are chosen so the FLOOR IS THE ONLY DISCRIMINATOR: across the
# whole pool no pair exceeds `max_dimer_bp=3` and none self-dimers at 4, and
# exactly one pair -- `_DG_PARTNER` with `_DG_FIXED` -- falls below a -1.0
# kcal/mol floor at 37 C. So anything the floor excludes, only the floor could
# have excluded, and the no-floor control must deliver it.

_DG_FIXED = "GGATCACAGTCT"
_DG_PARTNER = "TTCCAGATAGCG"
_DG_BENIGN = ["ACATTTTTTCGT", "CTCACAATAGAA", "TTTGAATACTCA"]
_DG_FLOOR = -1.0


def _clique_with(max_dimer_dg):
    """A clique optimizer with the floor set or unset, and no conditions.

    With no `ReactionConditions` attached the screen temperature is the 37 C
    default, which is what the fixtures above were chosen at.
    """
    return CliqueOptimizer(
        position_cache=_StubCache(),
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        config=OptimizerConfig(target_set_size=4, max_dimer_bp=3, max_dimer_dg=max_dimer_dg),
    )


def test_the_fixture_isolates_the_free_energy_floor():
    """Guard the fixtures themselves: only the floor separates these primers.

    If a later edit to these sequences made any pair fail the base-pair screen,
    the two tests below would still pass while proving nothing about the floor.
    """
    from neoswga.core.dimer import is_dimer_thermodynamic

    pool = [*_DG_BENIGN, _DG_PARTNER, _DG_FIXED]
    assert not [p for p in pool if is_dimer_fast(p, p, 4)], "a fixture self-dimers"
    bp_pairs = [(a, b) for a, b in itertools.combinations(pool, 2) if is_dimer_fast(a, b, 3)]
    assert not bp_pairs, f"the base-pair screen already separates {bp_pairs}"
    dg_pairs = [
        (a, b)
        for a, b in itertools.combinations(pool, 2)
        if is_dimer_thermodynamic(a, b, delta_g_threshold=_DG_FLOOR, temperature=37.0)
    ]
    assert dg_pairs == [
        (_DG_PARTNER, _DG_FIXED)
    ], f"expected exactly the partner/fixed pair below the floor, got {dg_pairs}"


def test_a_candidate_that_only_fails_the_free_energy_floor_is_never_delivered():
    """The configured floor has to exclude it, since nothing else can."""
    opt = _clique_with(_DG_FLOOR)
    # Make it the most attractive candidate, so only the screen keeps it out.
    opt._score_set = lambda s: 1.0 if _DG_PARTNER in s else 0.0
    result = opt.optimize(
        candidates=[*_DG_BENIGN, _DG_PARTNER], num_primers=4, fixed_primers=[_DG_FIXED]
    )

    assert _DG_PARTNER not in list(result.primers), (
        "the candidate that falls below the configured max_dimer_dg against the "
        "fixed primer was delivered; clique applied no free-energy floor"
    )


def test_without_a_floor_the_same_candidate_is_delivered():
    """The control, and the proof that the default path is unchanged.

    `max_dimer_dg` ships unset. With no floor this pair passes every screen
    clique applies, so it must still be delivered -- otherwise the test above
    would pass for some unrelated reason and the fix would be changing the
    default behaviour rather than honouring an opted-in setting.
    """
    opt = _clique_with(None)
    opt._score_set = lambda s: 1.0 if _DG_PARTNER in s else 0.0
    result = opt.optimize(
        candidates=[*_DG_BENIGN, _DG_PARTNER], num_primers=4, fixed_primers=[_DG_FIXED]
    )

    assert _DG_PARTNER in list(result.primers), (
        "with no floor configured this candidate clears every screen and should "
        "still be selected; the floor must be opt-in"
    )


def test_fixed_primers_that_only_fail_the_free_energy_floor_are_refused():
    """Fixed primers are held to the same screen as the candidates.

    A fixed pair below the configured floor cannot belong to a dimer-free set
    any more than one exceeding `max_dimer_bp` can, so the method refuses
    rather than deliver it.
    """
    result = _clique_with(_DG_FLOOR).optimize(
        candidates=list(_DG_BENIGN), num_primers=4, fixed_primers=[_DG_FIXED, _DG_PARTNER]
    )
    assert result.status == OptimizationStatus.ERROR
    assert "dimer" in (result.message or "").lower()


def test_the_floor_only_ever_removes_edges():
    """It can tighten the graph and never loosen it, at any temperature.

    `max_dimer_dg` is documented as only able to tighten the screen. The floor
    is consulted exclusively for pairs the base-pair matrix already cleared, so
    the edge count must be monotone in it.
    """
    import random

    random.seed(3)
    primers = list(
        dict.fromkeys("".join(random.choice("ACGT") for _ in range(14)) for _ in range(60))
    )

    def edges(max_dimer_dg, temp):
        graph = build_compatibility_graph(
            primers,
            max_dimer_bp=7,
            max_self_dimer_bp=8,
            check_self_dimers=True,
            max_dimer_dg=max_dimer_dg,
            temp=temp,
        )
        return graph.number_of_edges()

    unscreened = edges(None, 30.0)
    for temp in (30.0, 65.0):
        for floor in (-3.0, -6.0):
            assert (
                edges(floor, temp) <= unscreened
            ), f"floor {floor} at {temp} C added edges; it may only tighten"
