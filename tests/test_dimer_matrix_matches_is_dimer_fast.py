"""The vectorised dimer relation agrees with is_dimer_fast on every pair.

Finding A1. Measured on the real 449-candidate E. coli pool
(runs/gc_tiers/mid_ecoli/step3_df.csv) on an unloaded machine, is_dimer_fast costs
12.9 us per pair against secondary_structure.check_heterodimer's 627.1 us per pair
-- 49x cheaper. Neither is affordable 4.5 million times over (the pair count of a
3000-candidate pool): about 1.0 minute pairwise for the substring test, about 47
minutes for the thermodynamic one. It is the thermodynamic screen's cost this
module lets an optimizer avoid paying by pre-screening with the substring test
first. The matrix form is an exact reformulation for primers over ACGT, not an
approximation: a common substring of length t between seq_1 and revcomp(seq_2)
exists exactly when some t-mer of seq_1 is the reverse complement of some t-mer
of seq_2. Outside ACGT it is conservative rather than exact; see
test_ambiguity_codes_diverge_from_the_oracle_by_missing_a_pair below.
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


def test_ambiguity_codes_diverge_from_the_oracle_by_missing_a_pair():
    """The reformulation is exact only for primers over ACGT.

    is_dimer_fast matches "N" against "N" by plain character equality, so an
    aligned pair of ambiguity codes can extend a run past the threshold. This
    module's t-mer codes skip any t-mer containing a character outside ACGT,
    so it contributes nothing there. The two answers disagree on this pair and
    only in this direction: the matrix can miss a pair the oracle flags, never
    invent one the oracle does not. This is deliberate -- see the module
    docstring -- and pinned here so the divergence cannot change silently.
    """
    assert is_dimer_fast("ACGNT", "ANCGT", 3) is True
    matrix = dimer_matrix.build(["ACGNT", "ANCGT"], 3)
    assert bool(matrix.pairs[0, 1]) is False


def test_a_threshold_too_loose_for_the_code_space_raises():
    """params.schema.json allows max_dimer_bp up to 15 and max_k up to 30; the
    EquiPhi29 and Bst presets alone reach 15-25 base primers, so a combination
    that needs more t-mer codes than MAX_CODES allocates is realistic, not
    hypothetical. The error names both the parameter and the fallback so a
    caller (Tasks 2 and 5) can catch it and route to is_dimer_fast pairwise.
    """
    with pytest.raises(ValueError, match="max_dimer_bp") as excinfo:
        dimer_matrix.build(["A" * 20, "T" * 20], 10)
    assert "is_dimer_fast" in str(excinfo.value)


def test_a_pair_sharing_exactly_256_codes_is_still_flagged():
    """The matrix product accumulates in bool, not uint8.

    numpy accumulates an integer matmul in the input dtype, so the uint8 form
    this replaced summed a pair sharing exactly 256 (or 512, ...) distinct
    t-mer codes to 0 and reported it as not dimerising -- a wrong answer with
    no error and no warning. Constructed here rather than searched for: the
    concatenation of all 256 4-mers contains every 4-mer, so two orderings of
    that concatenation share all 256 codes exactly.

    Not reachable with a legal primer -- the shared count is bounded by the
    number of distinct t-mers in a primer and params.schema.json caps max_k at
    30 -- so this pins the arithmetic rather than a live design defect. The
    boolean form is also the faster of the two at the code-space ceiling; see
    the comment in dimer_matrix.build.
    """
    kmers = ["".join(p) for p in itertools.product("ACGT", repeat=4)]
    first = "".join(kmers)
    second = "".join(reversed(kmers))
    assert first != second

    matrix = dimer_matrix.build([first, second], 3)
    assert bool(matrix.pairs[0, 1]) is True
    assert is_dimer_fast(first, second, 3) is True


# ---------------------------------------------------------------------------
# A threshold the representation cannot hold must be refused, not ignored
# ---------------------------------------------------------------------------


def test_the_schema_refuses_a_threshold_the_matrix_cannot_represent():
    """`max_dimer_bp` 8 and above silently disabled the screen for a whole run.

    `MAX_CODES = 4**8` caps the representation at 7. Above it `build` raises,
    `DominatingSetOptimizer._build_dimer_matrix_for_greedy` catches the
    ValueError, logs `Dimer-aware selection disabled` at WARNING and returns
    None -- so the guard is off for the rest of the run, with one line among
    many to say so. Plan 2's threshold sweep confirmed this in all three pools
    at 8 and in no run at 3 through 7.

    The schema permitted up to 15, so a user could configure a value that
    removes the screen. It is now refused at validation time, where the message
    reaches the person who wrote it.
    """
    from neoswga.core.schema import load_schema

    prop = load_schema()["properties"]["max_dimer_bp"]
    assert prop["maximum"] == 7, "the schema admits a threshold the guard cannot enforce"
    assert prop["minimum"] == 1


def test_the_representation_ceiling_and_the_schema_agree():
    """Pin the two together, so raising MAX_CODES without the schema, or the
    schema without MAX_CODES, fails here rather than in a user's run.

    `build` uses `t = max_dimer_bp + 1` codes, so the highest representable
    threshold is one below log4(MAX_CODES), not equal to it. Derived here
    rather than written as a literal, and confirmed against `build` itself
    below so the derivation cannot drift from the behaviour.
    """
    import math

    from neoswga.core import dimer_matrix
    from neoswga.core.schema import load_schema

    highest = int(math.log(dimer_matrix.MAX_CODES, 4)) - 1
    assert load_schema()["properties"]["max_dimer_bp"]["maximum"] == highest

    primers = ["ACGTACGTACGT", "TTTTTTTTTTTT"]
    dimer_matrix.build(primers, highest)
    with pytest.raises(ValueError):
        dimer_matrix.build(primers, highest + 1)


def test_a_validator_run_rejects_the_disabling_value(tmp_path):
    """End to end through the validator a user actually runs."""
    from neoswga.core.param_validator import ParamValidator, ValidationLevel

    messages = ParamValidator().validate_params(
        {
            "data_dir": str(tmp_path),
            "fg_genomes": ["fg.fna"],
            "fg_prefixes": ["fg"],
            "max_dimer_bp": 8,
        }
    )
    errors = [m for m in messages if m.level == ValidationLevel.ERROR]
    assert any("max_dimer_bp" in (m.parameter or "") for m in errors), [
        (m.level, m.parameter, m.message) for m in messages
    ]
