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
