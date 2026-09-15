"""Compatibility against the selected panel, not against the whole pool.

Task 4 of the condition-aware pool design plan.

`dimer_matrix.build` allocates an n-by-n boolean array over every candidate. At
the shortlist sizes the pipeline used to deliver that was unremarkable, but
`candidate_retention="all_qc"` hands the search an order of magnitude more:

    2,000 candidates ->     4 MB
   20,670 candidates ->   427 MB
   50,000 candidates -> 2,500 MB

The Wolbachia design has 20,670 hard-QC survivors, so retaining them all makes
the matrix the largest single allocation in the run, and it grows quadratically
from there.

Almost none of it is read. The greedy screens each candidate against the
SELECTED panel, which is tens of primers, so it touches a thin slice of a
quadratic structure. Computing those pairs on demand and caching them costs a
bounded amount of memory and answers exactly the same question.

This is the interface `_would_dimerise` already uses, so it is a substitution
rather than a new code path.
"""

import pytest

from neoswga.core.dimer_matrix import build
from neoswga.core.lazy_dimer import LazyDimerCompatibility

# Complementary runs long enough to dimerise at the default threshold.
A = "AAAAAAAAAAAA"
T = "TTTTTTTTTTTT"
C = "CCCCCCCCCCCC"
G = "GGGGGGGGGGGG"
MIXED = "ACGTACGTACGT"

POOL = [A, T, C, G, MIXED]


@pytest.mark.parametrize("max_dimer_bp", [3, 4, 5])
def test_it_agrees_with_the_materialised_matrix(max_dimer_bp):
    """The substitution has to answer identically, or it is a different screen."""
    matrix = build(POOL, max_dimer_bp)
    lazy = LazyDimerCompatibility(max_dimer_bp)

    for candidate in POOL:
        for size in range(len(POOL)):
            selected = [p for p in POOL if p != candidate][:size]
            assert lazy.dimerises(candidate, selected) == matrix.dimerises(
                candidate, selected
            ), f"disagreement for {candidate} against {selected}"


def test_an_empty_panel_has_nothing_to_conflict_with():
    assert LazyDimerCompatibility(3).dimerises(A, []) is False


def test_a_self_complementary_pair_is_detected():
    assert LazyDimerCompatibility(3).dimerises(A, [T]) is True


def test_the_cache_is_bounded():
    """A long search must not accumulate a pair for everything it ever compared."""
    lazy = LazyDimerCompatibility(3, cache_size=8)

    for i in range(200):
        lazy.dimerises(f"{'ACGT' * 3}", [f"AAAC{i:04d}AAAA"[:12]])

    assert len(lazy._cache) <= 8


def test_a_pair_is_computed_once():
    lazy = LazyDimerCompatibility(3)
    lazy.dimerises(A, [T])
    before = lazy.computations
    lazy.dimerises(A, [T])

    assert lazy.computations == before


def test_the_pair_key_does_not_depend_on_argument_order():
    """A dimerises with B exactly when B dimerises with A."""
    lazy = LazyDimerCompatibility(3)
    lazy.dimerises(A, [T])
    before = lazy.computations
    lazy.dimerises(T, [A])

    assert lazy.computations == before, "the same pair was computed twice"


def test_it_allocates_nothing_quadratic_in_the_pool():
    """The point of the substitution.

    Constructing it must not depend on pool size at all: it is handed a
    threshold, not a pool.
    """
    import inspect

    signature = inspect.signature(LazyDimerCompatibility.__init__)
    assert "primers" not in signature.parameters
    assert "candidates" not in signature.parameters


def test_a_large_pool_gets_the_lazy_screen(monkeypatch):
    """The substitution has to happen, or it is a module nothing calls."""
    from neoswga.core.dominating_set_optimizer import (
        LAZY_DIMER_POOL_THRESHOLD,
        DominatingSetOptimizer,
    )

    optimizer = DominatingSetOptimizer.__new__(DominatingSetOptimizer)
    optimizer.max_dimer_bp = 3
    optimizer.relax_dimer_constraint_when_stuck = False

    def _oligo(index):
        """Distinct 12-mers: base-4 digits, so no two indices collide.

        A first attempt translated decimal digits to bases, which mapped 0, 4
        and 8 to the same letter; the pool deduplicated below the threshold and
        the test passed for the wrong reason.
        """
        letters = "ACGT"
        out = []
        for _ in range(12):
            out.append(letters[index % 4])
            index //= 4
        return "".join(out)

    small = [_oligo(i) for i in range(10)]
    assert not isinstance(
        optimizer._build_dimer_matrix_for_greedy(small), LazyDimerCompatibility
    ), "a small pool should keep the materialised matrix"

    big = [_oligo(i) for i in range(LAZY_DIMER_POOL_THRESHOLD + 1)]
    assert len(set(big)) == len(big), "the fixture pool must be genuinely distinct"
    assert isinstance(
        optimizer._build_dimer_matrix_for_greedy(big), LazyDimerCompatibility
    ), "a pool past the threshold still allocates a quadratic matrix"


def test_the_threshold_is_below_a_real_retained_pool():
    """20,670 hard-QC candidates is the measured Wolbachia figure."""
    from neoswga.core.dominating_set_optimizer import LAZY_DIMER_POOL_THRESHOLD

    assert LAZY_DIMER_POOL_THRESHOLD < 20_670
