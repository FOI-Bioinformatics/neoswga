"""The vectorised dimer relation agrees with is_dimer_fast on every pair.

Finding A1. is_dimer_fast costs 105 us per pair, which is 7.9 minutes for the
4.5 million pairs of a 3000-candidate pool, so it cannot be used as a pre-screen
for the thermodynamic check. The matrix form is an exact reformulation, not an
approximation: a common substring of length t between seq_1 and revcomp(seq_2)
exists exactly when some t-mer of seq_1 is the reverse complement of some t-mer
of seq_2.
"""

import itertools
import random

import pytest

from neoswga.core import dimer_matrix
from neoswga.core.dimer import is_dimer_fast


def _pool(n, k, seed):
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(k)) for _ in range(n)]


@pytest.mark.parametrize("max_dimer_bp", [3, 4, 5, 6])
@pytest.mark.parametrize("k", [8, 12, 16])
def test_matrix_agrees_with_is_dimer_fast(max_dimer_bp, k):
    primers = _pool(40, k, seed=k * 100 + max_dimer_bp)
    matrix = dimer_matrix.build(primers, max_dimer_bp)
    for i, j in itertools.combinations(range(len(primers)), 2):
        expected = is_dimer_fast(primers[i], primers[j], max_dimer_bp)
        assert bool(matrix.pairs[i, j]) == expected, (primers[i], primers[j], max_dimer_bp)


def test_matrix_is_symmetric():
    primers = _pool(30, 12, seed=7)
    matrix = dimer_matrix.build(primers, 3)
    assert (matrix.pairs == matrix.pairs.T).all()


def test_dimerises_against_a_selected_set():
    primers = ["AAAATTTTAAAA", "TTTTAAAATTTT", "CGCGCGCGCGCG"]
    matrix = dimer_matrix.build(primers, 3)
    assert matrix.dimerises("AAAATTTTAAAA", ["TTTTAAAATTTT"]) is True
    assert matrix.dimerises("CGCGCGCGCGCG", []) is False


def test_flagged_pairs_yields_upper_triangle_only():
    primers = _pool(20, 12, seed=11)
    matrix = dimer_matrix.build(primers, 3)
    pairs = list(matrix.flagged_pairs())
    assert all(i < j for i, j in pairs)
    assert len(pairs) == int(matrix.pairs.sum()) // 2


def test_lowercase_input_is_normalised():
    """One lowercase primer defeated every dimer check once already; see
    docs/AUDIT_optimize_step_2026-09-03.md finding 20."""
    matrix = dimer_matrix.build(["aaaattttaaaa", "TTTTAAAATTTT"], 3)
    assert bool(matrix.pairs[0, 1]) is True


def test_a_primer_shorter_than_t_never_dimerises():
    matrix = dimer_matrix.build(["ACG", "CGT"], 3)
    assert bool(matrix.pairs[0, 1]) is False


def test_a_threshold_no_primer_can_reach_short_circuits():
    """params.schema.json allows max_dimer_bp up to 15, which is 4**16 codes.

    No 12-mer contains a run of 16 bases, so the answer is all-False and the
    code space must not be allocated to say so.
    """
    matrix = dimer_matrix.build(["AAACCCGGGTTT", "ACCCGGGTTTAA"], 15)
    assert matrix.pairs.shape == (2, 2)
    assert not matrix.pairs.any()
