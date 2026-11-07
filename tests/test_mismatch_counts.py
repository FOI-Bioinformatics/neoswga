"""Counting background sites by how well they match.

Selectivity counts only exact matches, so the population that stringency
actually discriminates against -- sites carrying a mismatch or two -- is
invisible to the model. Occupancy has nothing to weight and additives have
nothing to act on.

Counting those sites turns out to be nearly free. The jellyfish `*_all.txt`
files already hold a count for every k-mer in the genome, so the 1-mismatch
background load of a primer is the sum over its `3k` neighbours: 30 lookups for
a 10-mer, not a rescan.

Two details decide whether the numbers mean anything.

**Canonicalisation.** Jellyfish runs with `-C`, so the file holds only the
lexicographically smaller of each k-mer and its reverse complement -- never
both. A primer and its reverse complement are the same physical site on
double-stranded DNA and must count the same. About half of the `3k` neighbours
of any primer are non-canonical, so looking them up as written would drop half
the background load and make selectivity look roughly twice as good as it is.
The existing `filter._get_rate_for_one_file` has exactly this shape (it tests
`parts[0] in primer_set`); it is latent there only because the default
candidate list is itself read out of those files and so is already canonical.

**Double counting.** Variants are summed over distinct canonical forms, because
two different 1-mismatch variants can share one canonical form when a primer is
near-palindromic. That handles collisions WITHIN a class. Collisions ACROSS
classes are the other half, fixed 2026-11-01 and covered at the bottom of this
file: the same stored k-mer can be a 1-mismatch neighbour of the primer read
one way and a 2-mismatch neighbour read the other, because the file holds only
one of each k-mer and its reverse complement. One physical site then landed in
two classes and `occupancy.weighted_site_load` weighted it at two different
occupancies.
"""

import pytest

from neoswga.core.mismatch_counts import (
    canonical_kmer,
    mismatch_class_counts,
    one_mismatch_variants,
)
from neoswga.core.thermodynamics import reverse_complement


@pytest.fixture
def genome(tmp_path):
    """A k-mer count file in the layout jellyfish -C writes: canonical only."""
    counts = {
        "AAAACCCCGG": 7,  # the primer under test, canonical
        "AAAACCCCGT": 3,  # one mismatch at the last base
        "TAAACCCCGG": 5,  # one mismatch at the first base
        "GGGGTTTTAA": 11,  # unrelated
    }
    # Store each under its canonical form, as jellyfish would.
    prefix = tmp_path / "g"
    with open(f"{prefix}_10mer_all.txt", "w") as fh:
        for kmer, count in counts.items():
            fh.write(f"{canonical_kmer(kmer)} {count}\n")
    return {"prefix": str(prefix), "counts": counts}


# ----------------------------------------------------------------------
# Canonicalisation
# ----------------------------------------------------------------------


def test_canonical_form_is_the_lexicographic_minimum():
    """Jellyfish's `-C` convention, which the count files follow."""
    assert canonical_kmer("GTTACATCGA") == "GTTACATCGA"
    assert canonical_kmer("TCGATGTAAC") == "GTTACATCGA"


def test_a_primer_and_its_reverse_complement_are_one_site(genome):
    """They are the same duplex on double-stranded DNA.

    Counting them differently is the bug this module exists to avoid: roughly
    half of every primer's mismatch neighbours are non-canonical, so a naive
    lookup silently halves the background load.
    """
    primer = "AAAACCCCGG"
    forward = mismatch_class_counts(primer, [genome["prefix"]], max_mismatches=1)
    reverse = mismatch_class_counts(
        reverse_complement(primer), [genome["prefix"]], max_mismatches=1
    )

    assert forward == reverse


# ----------------------------------------------------------------------
# Variant generation
# ----------------------------------------------------------------------


def test_there_are_three_variants_per_position():
    variants = one_mismatch_variants("ACGT")

    assert len(variants) == 3 * 4
    assert "ACGT" not in variants, "the sequence itself is not a mismatch of itself"


def test_every_variant_differs_at_exactly_one_position():
    seq = "AAAACCCCGG"
    for variant in one_mismatch_variants(seq):
        differences = sum(1 for a, b in zip(seq, variant) if a != b)
        assert differences == 1


# ----------------------------------------------------------------------
# The counts
# ----------------------------------------------------------------------


def test_exact_class_matches_the_file(genome):
    counts = mismatch_class_counts("AAAACCCCGG", [genome["prefix"]], max_mismatches=0)

    assert counts[0] == 7


def test_one_mismatch_class_sums_the_neighbours(genome):
    """Only the two planted neighbours are present in the file; every other
    neighbour is absent and contributes nothing."""
    counts = mismatch_class_counts("AAAACCCCGG", [genome["prefix"]], max_mismatches=1)

    assert counts[0] == 7
    assert counts[1] == 3 + 5


def test_the_primer_itself_is_not_counted_as_its_own_neighbour(genome):
    """Class 1 must exclude the exact match, or the two classes overlap and the
    same site is weighted twice at two different occupancies."""
    counts = mismatch_class_counts("AAAACCCCGG", [genome["prefix"]], max_mismatches=1)

    assert counts[1] == 8
    assert counts[0] + counts[1] == 15


def test_a_primer_absent_from_the_genome_counts_zero(genome):
    counts = mismatch_class_counts("TTTTTTTTTT", [genome["prefix"]], max_mismatches=1)

    assert counts[0] == 0


def test_the_cap_is_respected(genome):
    counts = mismatch_class_counts("AAAACCCCGG", [genome["prefix"]], max_mismatches=0)

    assert set(counts) == {0}, "asking for no mismatches must not report class 1"


def test_counts_sum_across_prefixes(genome, tmp_path):
    """Multiple background genomes contribute additively -- a primer binding
    two hosts carries the load of both."""
    second = tmp_path / "h"
    with open(f"{second}_10mer_all.txt", "w") as fh:
        fh.write(f"{canonical_kmer('AAAACCCCGG')} 4\n")

    counts = mismatch_class_counts("AAAACCCCGG", [genome["prefix"], str(second)], max_mismatches=0)

    assert counts[0] == 7 + 4


def test_a_missing_count_file_is_not_silently_zero(tmp_path):
    """Zero from an absent file is indistinguishable from zero background
    binding, and the second is a strong claim. Absent inputs raise."""
    with pytest.raises(FileNotFoundError):
        mismatch_class_counts("AAAACCCCGG", [str(tmp_path / "nonexistent")], max_mismatches=0)


# ----------------------------------------------------------------------
# The property that makes this worth doing
# ----------------------------------------------------------------------


def test_mismatched_sites_outnumber_exact_ones_at_genome_scale():
    """Why exact-match counting understates background so badly.

    A primer has `3k` one-mismatch neighbours, and in a genome large enough for
    each to occur, the near-match population approaches that multiple of the
    exact one. Measured on random sequence:

        50 kb    2.5x
        500 kb  15.4x
        5 Mb    27.0x     (approaching the 30 neighbours of a 10-mer)

    So at bacterial scale a model counting only perfect matches is looking at
    roughly one twenty-eighth of the background a primer actually binds -- not
    an imprecise estimate of the load, but the wrong population.

    50 kb here because it is fast and the effect is already unambiguous.
    """
    import os
    import random
    import tempfile
    from collections import Counter

    from neoswga.core.mismatch_counts import canonical_kmer

    size = 50_000
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(size))

    counts = Counter(canonical_kmer(seq[i : i + 10]) for i in range(size - 9))
    directory = tempfile.mkdtemp()
    prefix = os.path.join(directory, "g")
    with open(f"{prefix}_10mer_all.txt", "w") as handle:
        for kmer, count in counts.items():
            handle.write(f"{kmer} {count}\n")

    probes = [seq[i : i + 10] for i in range(0, size - 9, 250)][:200]
    rows = [mismatch_class_counts(p, [prefix], max_mismatches=1) for p in probes]

    exact = sum(r[0] for r in rows)
    near = sum(r[1] for r in rows)

    assert exact > 0, "the probes should match the genome they came from"
    assert near > exact, (
        f"near-match load {near} did not exceed exact load {exact}; exact "
        f"counting would not be understating background at all"
    )


# ----------------------------------------------------------------------
# Disjointness ACROSS classes (2026-11-01)
# ----------------------------------------------------------------------
#
# A site is one physical duplex and must be weighted at one occupancy, so it
# belongs to the class of its BEST match. Canonicalisation is what makes that
# non-obvious: a count file holds only one of each k-mer and its reverse
# complement, so the SAME stored k-mer can be a 1-mismatch neighbour of the
# primer read one way and a 2-mismatch neighbour read the other.
#
# `ACTGTATAGT` is such a primer, verified rather than assumed. Its reverse
# complement is `ACTATACAGT`, and the canonical form `ACTATAAAGT` is reached
# both as the canonicalisation of the distance-1 variant `ACTTTATAGT` and
# directly as the distance-2 variant `ACTATAAAGT`.

COLLIDING_PRIMER = "ACTGTATAGT"
#: The one stored k-mer both classes can reach, and the distinctive count it
#: carries so a double count is visible in the total rather than inferred.
SHARED_CANONICAL = "ACTATAAAGT"
SHARED_COUNT = 1000


@pytest.fixture
def colliding_genome(tmp_path):
    """A count file holding the shared form and nothing that overlaps it."""
    prefix = tmp_path / "c"
    with open(f"{prefix}_10mer_all.txt", "w") as fh:
        fh.write(f"{canonical_kmer(COLLIDING_PRIMER)} 2\n")
        fh.write(f"{canonical_kmer(SHARED_CANONICAL)} {SHARED_COUNT}\n")
    return str(prefix)


def test_the_shared_form_is_reachable_from_both_classes(colliding_genome):
    """The premise of the next three tests, asserted rather than trusted."""
    from neoswga.core.mismatch_counts import _variants_at_distance

    at_one = {canonical_kmer(v) for v in _variants_at_distance(COLLIDING_PRIMER, 1)}
    at_two = {canonical_kmer(v) for v in _variants_at_distance(COLLIDING_PRIMER, 2)}

    assert SHARED_CANONICAL in at_one
    assert SHARED_CANONICAL in at_two


def test_a_site_reachable_at_two_distances_is_counted_once(colliding_genome):
    """Before the fix its count appeared in class 1 AND class 2."""
    counts = mismatch_class_counts(COLLIDING_PRIMER, [colliding_genome], max_mismatches=2)

    assert counts[1] + counts[2] == SHARED_COUNT, (
        f"the shared site contributes {counts[1] + counts[2]} across classes "
        f"1 and 2, but the file holds {SHARED_COUNT} of it"
    )


def test_it_counts_in_the_nearer_class(colliding_genome):
    """Best match wins: a 1-mismatch duplex is not a 2-mismatch duplex, and
    the occupancy it is weighted at differs by a whole mismatch penalty."""
    counts = mismatch_class_counts(COLLIDING_PRIMER, [colliding_genome], max_mismatches=2)

    assert counts[1] == SHARED_COUNT
    assert counts[2] == 0


def test_the_exact_match_still_wins_over_both(colliding_genome):
    """Class 0 claims first, which is what the old in-loop filter did for the
    primer's own canonical form and what the ordering now does generally."""
    counts = mismatch_class_counts(COLLIDING_PRIMER, [colliding_genome], max_mismatches=2)

    assert counts[0] == 2


def test_the_default_depth_is_unchanged(colliding_genome):
    """At `max_mismatches=1` classes 0 and 1 were already disjoint, so no
    delivered panel moves: the fix is reachable only above the default."""
    counts = mismatch_class_counts(COLLIDING_PRIMER, [colliding_genome], max_mismatches=1)

    assert counts == {0: 2, 1: SHARED_COUNT}


def test_no_class_shares_a_canonical_form_with_another():
    """Over the primers where the collision actually occurs, not over random
    ones: a property test that never reaches the case proves nothing."""
    from neoswga.core.mismatch_counts import _variants_at_distance

    def classes(primer):
        claimed, out = set(), {}
        for distance in range(0, 3):
            if distance == 0:
                forms = {canonical_kmer(primer)}
            else:
                forms = {canonical_kmer(v) for v in _variants_at_distance(primer, distance)}
            out[distance] = forms - claimed
            claimed |= forms
        return out

    colliding = []
    for primer in ("ACTGTATAGT", "TATGCGTATA", "ATCTGCAAAT"):
        at_one = {canonical_kmer(v) for v in _variants_at_distance(primer, 1)}
        at_two = {canonical_kmer(v) for v in _variants_at_distance(primer, 2)}
        if at_one & at_two:
            colliding.append(primer)

    assert colliding, "the fixtures should be primers that actually collide"
    for primer in colliding:
        by_class = classes(primer)
        for left in by_class:
            for right in by_class:
                if left < right:
                    assert not (
                        by_class[left] & by_class[right]
                    ), f"{primer}: classes {left} and {right} share a form"
