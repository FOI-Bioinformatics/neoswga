"""Where the mismatch sits, when the model is asked to care, and not otherwise.

Phase 4b. The shipped model subtracts one number per mismatch
(`occupancy.mismatch_tm`, registry status `assumed`). This file pins the
alternative, which is reachable through `mismatch_model` and off by default.

Six things are checked, and the fourth is the one that matters most:

1. **The oracle.** One short duplex with one internal mismatch equals a value
   computed by hand from the published table, with the arithmetic written out
   below so a reader can check it without running anything. A perfect duplex
   equals what the existing perfect-match walk returns, exactly.
1b. **The table's domain is internal positions.** A mismatch at either terminus
   does not use it, and the clamp that used to substitute zero penalty there is
   now a guard that measurably never binds.
2. **The reduction.** With every mismatch weighted alike, the per-neighbour
   load equals the uniform load. The two then differ only in the order the sum
   is taken, so the tolerance is floating-point and is stated as such.
3. **Canonicalisation is not re-derived.** A near-palindromic primer has
   neighbours that collapse onto one canonical form, and that form is one site
   group counted once.
4. **The default path does not move.** With the key unset, or set to
   `uniform`, every number and every mode string is what it was.
5. **Identity and provenance.** The model enters `request_hash` and the run
   manifest, every consumer of a weighted load records which model produced it,
   and the extrapolation's size is stated wherever the model runs.
"""

import json
import math
import zlib

import pytest

from neoswga.core import mismatch_model as mm
from neoswga.core import occupancy
from neoswga.core.mismatch_counts import canonical_kmer, mismatch_class_counts
from neoswga.core.mismatch_sites import neighbour_sites
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import (
    DELTA_G_MISMATCH,
    NN_INIT_CORRECTIONS,
    calculate_enthalpy_entropy,
    compute_free_energy_for_two_strings,
)
from neoswga.core.utility import complement

PRIMER = "ACCACAGATAGC"


@pytest.fixture
def conditions():
    return ReactionConditions(temp=30.0, na_conc=50.0, mg_conc=10.0)


@pytest.fixture
def table(tmp_path):
    """A canonical-only count table, in the layout `jellyfish -C` writes.

    Holds the primer, two one-mismatch neighbours and an unrelated k-mer, so a
    load over it is a small number a reader can follow.
    """
    counts = {
        PRIMER: 9,
        "ACCACAGATAGA": 4,  # one mismatch on the 3'-terminal base
        "ACCACAGATTGC": 6,  # one mismatch three bases in from the 3' end
        "GGGGTTTTAACC": 13,  # unrelated
    }
    prefix = tmp_path / "fg"
    with open(f"{prefix}_12mer_all.txt", "w") as handle:
        for kmer, count in counts.items():
            handle.write(f"{canonical_kmer(kmer)} {count}\n")
    return str(prefix)


@pytest.fixture(scope="module")
def counting_table(tmp_path_factory):
    """A table holding every neighbour of every primer the counting test uses.

    Built so the class sums are large numbers rather than a handful: a table
    where almost every lookup misses cannot distinguish a correct grouping from
    one that looks its neighbours up in the wrong form. Counts are deterministic
    and distinct per canonical form, so a form counted twice shows up as a sum
    that is wrong by a specific amount rather than by chance.

    One line per CANONICAL form, as `jellyfish -C` writes it, which means
    aggregating before writing: two neighbours sharing a canonical form are one
    row with one count, and that is exactly the case under test.
    """
    from neoswga.core.mismatch_counts import _variants_at_distance

    primers = (PRIMER, *MOSTLY_NON_CANONICAL)
    rows: dict[str, int] = {}
    for primer in primers:
        for kmer in (primer, *_variants_at_distance(primer, 1)):
            canonical = canonical_kmer(kmer)
            if canonical not in rows:
                # Deterministic and never zero, so every group contributes to
                # the sum under test. `hash()` is salted per process for
                # strings, so it is not usable for a count a reader may want to
                # reproduce; crc32 is.
                rows[canonical] = 1 + zlib.crc32(canonical.encode()) % 97

    prefix = tmp_path_factory.mktemp("counting") / "cnt"
    with open(f"{prefix}_12mer_all.txt", "w") as handle:
        for canonical, count in sorted(rows.items()):
            handle.write(f"{canonical} {count}\n")
    return str(prefix)


# ----------------------------------------------------------------------
# 1. The oracle
# ----------------------------------------------------------------------


def test_one_internal_mismatch_equals_the_hand_computed_table_value():
    """Arithmetic written out, so the test checks the table and not itself.

    Primer 5'-G C A T G-3' against the genomic site 5'-G C G T G-3'. The primer
    anneals to that site's complement, read in the same index order, so the
    duplex is

        x (primer,   5'->3'):  G  C  A  T  G
        y (template, 3'->5'):  C  G  C  A  C

    Position 3 is A:C, a mismatch; the other four are Watson-Crick.

    Terminal corrections, applied only at a Watson-Crick end:
        G/C at the 5' end   +0.98
        G/C at the 3' end   +0.98

    Nearest-neighbour doublets, as `x[i-1]x[i]/y[i-1]y[i]`:
        GC/CG   -2.24   tabulated as written
        CA/GC   +0.75   tabulated as written
        AT/CA     --    absent as written; the same stack rotated 180 degrees
                        is AC/TA, +0.77
        TG/AC   -1.45   tabulated as written

    No palindrome correction: complement(x) is CGTAC and reverse(y) is CACGC.

    Total: 0.98 + 0.98 - 2.24 + 0.75 + 0.77 - 1.45 = -0.21 kcal/mol.
    """
    expected = (
        NN_INIT_CORRECTIONS["G/C"]
        + NN_INIT_CORRECTIONS["G/C"]
        + DELTA_G_MISMATCH["GC/CG"]
        + DELTA_G_MISMATCH["CA/GC"]
        + DELTA_G_MISMATCH["AC/TA"]
        + DELTA_G_MISMATCH["TG/AC"]
    )

    assert expected == pytest.approx(-0.21, abs=1e-12)
    assert mm.duplex_delta_g("GCATG", "GCGTG") == pytest.approx(expected, abs=1e-12)


def test_the_rotation_is_what_covers_the_three_prime_side_doublet():
    """Without it that doublet falls through to the flat penalty.

    This is the reason the new module walks the duplex itself: the shipped walk
    has no rotation, so it charges `penalty` for one doublet of every mismatch.
    """
    assert "AT/CA" not in DELTA_G_MISMATCH
    assert "AC/TA" in DELTA_G_MISMATCH

    with_rotation = mm.duplex_delta_g("GCATG", "GCGTG")
    without_rotation = with_rotation - DELTA_G_MISMATCH["AC/TA"] + 4.0

    assert without_rotation > with_rotation
    assert compute_free_energy_for_two_strings("GCATG", "CGCAC") == pytest.approx(
        without_rotation, abs=1e-12
    )


@pytest.mark.parametrize("primer", ["GCATG", PRIMER, "AAAACCCCGG", "ACGTACGT"])
def test_a_perfect_duplex_equals_the_existing_perfect_match_path(primer):
    """Exactly, not approximately. A perfect duplex reaches neither the
    rotation nor the early stop, so the two walks must agree bit for bit."""
    assert mm.duplex_delta_g(primer, primer) == compute_free_energy_for_two_strings(
        primer, complement(primer)
    )


def test_a_mismatch_never_stabilises_a_site_above_its_perfect_match(conditions):
    """The guard against the model turning a rejection into a pass.

    At a terminal position the nearest-neighbour arithmetic can net out
    positive, because a mismatched end also loses that end's initiation
    correction. An unclamped shift there would raise a mismatched site's
    occupancy above the perfect-match site's.
    """
    dh, _ = calculate_enthalpy_entropy(PRIMER)
    tm = conditions.calculate_effective_tm(PRIMER)

    shifts = [
        mm.mismatch_delta_tm(PRIMER, reading, dh_kcal=dh, tm=tm)
        for group in neighbour_sites(PRIMER, [], 1, tables=[])
        for reading in group.readings
        if reading.pairs
    ]

    assert shifts
    assert max(shifts) <= 0.0
    assert min(shifts) < -1.0, "nothing was destabilised at all; the walk is not running"


def test_the_duplex_term_alone_does_not_single_out_the_three_prime_base(conditions):
    """Why the 3'-terminal rule has to be a separate term.

    The duplex term cannot express "a mismatched 3' terminus is largely a
    non-site": that position is outside the internal table's domain and takes
    the uniform penalty, which is milder than the model charges most interior
    positions. So the extension rule is additional rather than redundant.
    """
    dh, _ = calculate_enthalpy_entropy(PRIMER)
    tm = conditions.calculate_effective_tm(PRIMER)

    by_offset: dict[int, list[float]] = {}
    for group in neighbour_sites(PRIMER, [], 1, tables=[]):
        for reading in group.readings:
            if not reading.pairs:
                continue
            offset = reading.offsets_from_three_prime[0]
            by_offset.setdefault(offset, []).append(
                mm.mismatch_delta_tm(PRIMER, reading, dh_kcal=dh, tm=tm, uniform_penalty=4.0)
            )

    mean = {off: sum(v) / len(v) for off, v in by_offset.items()}
    interior = [value for off, value in mean.items() if off not in (0, len(PRIMER) - 1)]

    assert mean[0] > min(interior), (
        "the duplex term already penalises the 3'-terminal base hardest, which "
        "would make the separate extension rule redundant rather than additional"
    )


# ----------------------------------------------------------------------
# 1b. The table's domain is internal positions
# ----------------------------------------------------------------------


TERMINAL_SPREAD = ("GCATG", "ACCACAGATAGC", "ACGTACGTACGT", "AAAAAACCCCCC", "GCGCGCATATATCGCG")


@pytest.mark.parametrize("primer", TERMINAL_SPREAD)
def test_a_terminal_mismatch_takes_the_uniform_penalty_not_the_table(primer, conditions):
    """`DELTA_G_MISMATCH` is an INTERNAL single-mismatch set.

    A terminal position has one flanking stack instead of two, and
    dangling-end / terminal-mismatch parameters are the applicable set, which
    this repository does not carry. Applying the internal table there measured
    as a NET STABILISING result for one neighbour of a 12-mer -- a mismatched
    host site weighted as well as a perfect one. So the model declines to use
    the table at a terminus and charges the uniform penalty instead, exactly as
    it already declines for a tandem-mismatch doublet.
    """
    dh, _ = calculate_enthalpy_entropy(primer)
    tm = conditions.calculate_effective_tm(primer)
    penalty = 4.0

    terminal, interior = [], []
    for group in neighbour_sites(primer, [], 1, tables=[]):
        for reading in group.readings:
            if not reading.pairs:
                continue
            shift = mm.mismatch_delta_tm(
                primer, reading, dh_kcal=dh, tm=tm, uniform_penalty=penalty
            )
            (terminal if reading.is_at_a_duplex_terminus else interior).append(shift)

    assert len(terminal) == 6, "three neighbours at each of the two termini"
    for shift in terminal:
        assert shift == pytest.approx(-penalty, abs=1e-12)

    assert interior, "a primer with no interior position would not test anything"


def test_the_terminal_fallback_follows_the_runs_resolved_penalty(conditions):
    """Not the hardcoded 4.0. The uniform path resolves `mismatch_penalty`
    once per run and the fallback must be the same number, or the two halves of
    one model would charge a terminal mismatch differently."""
    dh, _ = calculate_enthalpy_entropy(PRIMER)
    tm = conditions.calculate_effective_tm(PRIMER)
    reading = next(
        r
        for group in neighbour_sites(PRIMER, [], 1, tables=[])
        for r in group.readings
        if r.is_at_a_duplex_terminus
    )

    assert mm.mismatch_delta_tm(
        PRIMER, reading, dh_kcal=dh, tm=tm, uniform_penalty=7.5
    ) == pytest.approx(-7.5)


@pytest.mark.parametrize("primer", TERMINAL_SPREAD)
def test_the_clamp_never_binds_at_an_internal_position(primer, conditions):
    """The clamp is a guard, not a model, and this is what makes that measured.

    Before the terminal positions were taken out of the table's domain, the
    clamp substituted EXACTLY ZERO penalty -- the most favourable value
    available -- for 1 of 36 neighbours of a 12-mer, all at the 5' terminus.
    If it ever binds at an interior position that is a finding about the
    nearest-neighbour arithmetic and not something to clamp away, so this test
    reports where.
    """
    dh, _ = calculate_enthalpy_entropy(primer)
    tm = conditions.calculate_effective_tm(primer)
    tm_k = tm + 273.15
    perfect = mm.duplex_delta_g(primer, primer)

    bound = []
    for group in neighbour_sites(primer, [], 1, tables=[]):
        for reading in group.readings:
            if not reading.pairs or reading.is_at_a_duplex_terminus:
                continue
            raw = tm_k * (mm.duplex_delta_g(primer, reading.neighbour) - perfect) / dh
            if raw > 0:
                bound.append((reading.neighbour, reading.offsets_from_three_prime, raw))

    assert not bound, (
        f"the clamp binds at an INTERNAL position for {primer}: {bound}. That is "
        f"a finding about the nearest-neighbour arithmetic at an interior "
        f"position, not something the clamp should be hiding"
    )


@pytest.mark.parametrize("primer", TERMINAL_SPREAD)
@pytest.mark.parametrize("model", [mm.POSITION_DEPENDENT, mm.POSITION_DEPENDENT_THREE_PRIME])
def test_no_mismatched_neighbour_is_weighted_at_or_above_its_perfect_match(
    primer, model, conditions
):
    """The guard, stated as the thing it protects rather than as the clamp.

    A one-mismatch site weighted at the perfect-match occupancy is a mismatch
    turning a rejection into a pass, whatever arithmetic produced it.
    """
    from neoswga.core.occupancy import site_occupancy

    dh, _ = calculate_enthalpy_entropy(primer)
    tm = conditions.calculate_effective_tm(primer)
    matched = site_occupancy(dh, tm, conditions.temp)
    three_prime = model == mm.POSITION_DEPENDENT_THREE_PRIME

    for group in neighbour_sites(primer, [], 1, tables=[]):
        for reading in group.readings:
            if not reading.pairs:
                continue
            shift = mm.mismatch_delta_tm(primer, reading, dh_kcal=dh, tm=tm, uniform_penalty=4.0)
            weight = site_occupancy(dh, tm + shift, conditions.temp)
            if three_prime:
                weight *= mm._extension_factor(reading, mm.THREE_PRIME_WINDOW)
            assert weight < matched, (
                f"{reading.neighbour} at offsets {reading.offsets_from_three_prime} "
                f"is bound {weight} against the perfect match's {matched}"
            )


# ----------------------------------------------------------------------
# 2. The reduction
# ----------------------------------------------------------------------


def test_weighting_every_position_alike_reproduces_the_uniform_load(table, conditions):
    """The two paths then differ only in the order the sum is taken.

    The uniform path sums the counts in a mismatch class and multiplies the
    total by one occupancy; the per-neighbour path multiplies each group's
    count by the same occupancy and sums. Those are the same arithmetic in a
    different order, so the tolerance is floating-point reassociation and
    nothing else: 1e-9 relative, which is some six orders of magnitude tighter
    than any difference the model itself produces (the shifts above move Tm by
    6 to 19 C against the uniform 4). A failure at this tolerance means the two
    paths disagree about the counting, not about the position.
    """
    uniform = occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0, mm.UNIFORM)
    flattened = mm.position_dependent_site_load(
        [PRIMER], [table], conditions, 1, 4.0, uniform_delta_tm=True
    )

    assert uniform > 0
    assert flattened == pytest.approx(uniform, rel=1e-9)


def test_the_position_dependent_load_differs_from_the_uniform_one(table, conditions):
    """If it did not, the key would be inert whatever it was wired to."""
    uniform = occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0, mm.UNIFORM)
    modelled = occupancy.weighted_site_load(
        [PRIMER], [table], conditions, 1, 4.0, mm.POSITION_DEPENDENT
    )

    assert modelled != pytest.approx(uniform, rel=1e-6)
    assert modelled < uniform, (
        "every per-neighbour shift measured on this primer is larger than the "
        "uniform 4 C, so the modelled load must be the smaller number"
    )


def test_the_three_prime_rule_only_moves_terminal_mismatches(table, conditions):
    """Separately switchable, and it must not touch the duplex term."""
    duplex_only = occupancy.weighted_site_load(
        [PRIMER], [table], conditions, 1, 4.0, mm.POSITION_DEPENDENT
    )
    with_rule = occupancy.weighted_site_load(
        [PRIMER], [table], conditions, 1, 4.0, mm.POSITION_DEPENDENT_THREE_PRIME
    )

    assert with_rule < duplex_only, (
        "the table holds a site with a 3'-terminal mismatch, so the rule has " "something to act on"
    )
    # Exact-match sites carry no mismatch, so the rule can never remove all of
    # the load.
    assert with_rule > 0


def test_the_exact_match_class_is_untouched_by_either_term(table, conditions):
    """At depth 0 there is no mismatch, so all three models agree exactly."""
    loads = {
        model: occupancy.weighted_site_load([PRIMER], [table], conditions, 0, 4.0, model)
        for model in mm.MISMATCH_MODELS
    }

    assert len(set(loads.values())) == 1
    assert loads[mm.UNIFORM] > 0


# ----------------------------------------------------------------------
# 3. Canonicalisation is not re-derived
# ----------------------------------------------------------------------


def test_a_palindromic_neighbour_is_one_site_group_not_two():
    """Two neighbours that are each other's reverse complement are one count.

    `canonical_kmer` is the single definition of the form a table stores, and
    re-deriving it is how a near-palindromic primer's background load comes to
    be counted twice. The group carries both readings and the count once, so
    summing `count` over the groups cannot double count whatever a scoring
    function does with the readings.
    """
    primer = "ACGTACGTACGT"
    groups = neighbour_sites(primer, [], 1, tables=[])

    canonicals = [group.canonical for group in groups]
    assert len(canonicals) == len(set(canonicals)), "a canonical form was returned twice"

    doubled = [group for group in groups if len(group.readings) > 1]
    assert doubled, (
        "this primer was chosen because it is palindromic enough to collapse "
        "neighbours onto one canonical form; it no longer does, so the test "
        "no longer checks what it claims"
    )
    for group in doubled:
        readings = {reading.neighbour for reading in group.readings}
        assert len(readings) == len(group.readings)
        assert len({canonical_kmer(n) for n in readings}) == 1


# Primers whose neighbours are mostly NON-canonical, so looking one up as
# written misses. For `ACCACAGATAGC` only 1 of 36 is non-canonical, which left
# the counting test below passing under the very mutation it is written
# about -- grouping by the raw variant instead of its canonical form -- because
# a raw lookup hits anyway 35 times out of 36. For these, 18 of 36 are
# non-canonical and collapse to 19 groups carrying two readings each.
MOSTLY_NON_CANONICAL = ("ACGTACGTACGT", "ACGCGCGCGCGT", "ACGTTGCAACGT")


@pytest.mark.parametrize("primer", (PRIMER, *MOSTLY_NON_CANONICAL))
def test_the_group_counts_match_the_uniform_class_counts(primer, counting_table):
    """The same counting, grouped differently. Checked per class, because that
    is the level at which the uniform path deduplicates.

    **What gives this test its power is the dense table, not the primer.** The
    first version used a four-row fixture, so almost every lookup missed under
    the correct grouping as well as under a wrong one, and the invariant held
    under the very mutation it is written about -- grouping by the raw variant
    instead of its canonical form. Against `counting_table`, which holds a
    count for every neighbour of every primer here, the same mutation moves
    `ACCACAGATAGC`'s distance-1 total from 1821 to 1796 and the test fails.

    Measured, and worth knowing before adding more primers: a perfectly
    self-complementary primer does NOT catch that mutation, whatever the table.
    Its 36 neighbours form 18 reverse-complement pairs with one canonical
    member each, so summing the 18 raw forms that happen to be canonical gives
    exactly the same total as summing the 18 distinct canonical forms. They are
    parametrised anyway because they exercise the two-reading grouping, which
    `test_a_palindromic_neighbour_is_one_site_group_not_two` is the test for.
    """
    classes = mismatch_class_counts(primer, [counting_table], 1)
    groups = neighbour_sites(primer, [counting_table], 1)

    for distance, total in classes.items():
        regrouped = sum(g.count for g in groups if g.distance == distance)
        assert regrouped == total, distance


@pytest.mark.parametrize("primer", MOSTLY_NON_CANONICAL)
def test_half_of_these_primers_neighbours_are_non_canonical(primer):
    """The premise the test above rests on, asserted so it cannot rot.

    `mismatch_sites`'s docstring says "roughly half of any primer's 3k
    neighbours are non-canonical". That is true only for a self-complementary
    primer; for an ordinary one it is 1 in 36. If these primers stop being
    self-complementary the counting test silently loses its power again.
    """
    from neoswga.core.thermodynamics import reverse_complement

    variants = {
        reading.neighbour
        for group in neighbour_sites(primer, [], 1, tables=[])
        for reading in group.readings
        if reading.pairs
    }
    non_canonical = [v for v in variants if canonical_kmer(v) != v]

    assert primer == reverse_complement(primer), "this primer is not self-complementary"
    assert len(non_canonical) >= len(variants) // 3, (
        f"only {len(non_canonical)} of {len(variants)} neighbours are "
        f"non-canonical; this primer no longer exercises the canonicalisation"
    )


def test_the_primers_own_canonical_form_is_not_a_mismatch_group():
    """It would be weighted at two different occupancies if it were."""
    primer = "ACGTACGTACGT"
    own = canonical_kmer(primer)
    groups = neighbour_sites(primer, [], 1, tables=[])

    assert [g.canonical for g in groups].count(own) == 1
    assert next(g for g in groups if g.canonical == own).distance == 0


def test_a_missing_table_makes_the_load_unavailable_rather_than_zero(tmp_path, conditions):
    """A neighbour whose count cannot be looked up makes the load unknown.

    Zero sites and an unreadable table are different claims, and the second is
    a strong one to make from a missing input. Both models must refuse the same
    way, or switching the key would change whether a run fails.
    """
    absent = str(tmp_path / "nothing")

    for model in mm.MISMATCH_MODELS:
        with pytest.raises((FileNotFoundError, OSError)):
            occupancy.weighted_site_load([PRIMER], [absent], conditions, 1, 4.0, model)


# ----------------------------------------------------------------------
# 4. The default path does not move
# ----------------------------------------------------------------------


def test_an_unset_key_gives_the_uniform_load_and_the_uniform_mode(table, conditions, monkeypatch):
    """Proved by comparison, not by reading the dispatch.

    The load with the key absent must equal the load with the model named
    explicitly, and the reported mode must be the string every existing output
    already carries.
    """
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "mismatch_model", None, raising=False)

    assert mm.resolve_mismatch_model() == mm.UNIFORM
    assert mm.site_load_mode(None) == "occupancy"
    assert mm.site_load_mode(mm.UNIFORM) == "occupancy"

    unset = occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0)
    named = occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0, mm.UNIFORM)

    assert unset == named


def test_the_key_set_to_uniform_changes_nothing(table, conditions, monkeypatch):
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "mismatch_model", None, raising=False)
    before = occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0)

    monkeypatch.setattr(parameter, "mismatch_model", "uniform", raising=False)
    after = occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0)

    assert after == before


def test_a_position_dependent_load_is_never_reported_as_the_uniform_one():
    """The `mode` return exists for this."""
    assert mm.site_load_mode(mm.POSITION_DEPENDENT) == "occupancy-position-dependent"
    assert (
        mm.site_load_mode(mm.POSITION_DEPENDENT_THREE_PRIME)
        == "occupancy-position-dependent-3prime"
    )
    assert mm.site_load_mode(mm.POSITION_DEPENDENT) != "occupancy"


def test_an_unrecognised_model_is_refused_rather_than_falling_back():
    """Running the old model for someone who asked for the new one is the
    inert-option defect with an extra step."""
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest):
        mm.resolve_mismatch_model("position-aware")


def test_the_assessment_evidence_follows_the_model_that_ran():
    """The hardcoded line described the uniform model on every assessment."""
    from neoswga.core.panel_evaluation import _mismatch_evidence

    uniform = _mismatch_evidence("occupancy")
    assert set(uniform) == {"mismatch_penalty"}
    assert "uniform" in uniform["mismatch_penalty"]

    modelled = _mismatch_evidence("occupancy-position-dependent")
    assert set(modelled) == {"mismatch_model"}
    assert "estimated" in modelled["mismatch_model"]
    assert "37 C" in modelled["mismatch_model"]
    assert "3'-terminal" not in modelled["mismatch_model"]

    with_rule = _mismatch_evidence("occupancy-position-dependent-3prime")
    assert "3'-terminal" in with_rule["mismatch_model"]

    # An exact-count fallback is still the uniform record: no model ran.
    assert set(_mismatch_evidence("exact")) == {"mismatch_penalty"}
    assert set(_mismatch_evidence(None)) == {"mismatch_penalty"}


def test_both_new_quantities_have_an_evidence_record():
    """A model with no record is a number with no evidence behind it that looks
    exactly like one with evidence."""
    from neoswga.core.model_evidence import load_evidence

    evidence = load_evidence()

    duplex = evidence["mismatch_duplex_delta_g"]
    assert duplex.status == "estimated", (
        "the table is measured; THIS USE of it is an extrapolation, and the "
        "status describes the claim rather than the table"
    )
    assert "37 C" in duplex.temperature_domain
    assert "1 M NaCl" in duplex.buffer_domain or "1 M Na" in duplex.buffer_domain

    terminal = evidence["three_prime_mismatch_extension"]
    assert terminal.status == "assumed"
    assert "would settle it" in terminal.notes or "What would settle it" in terminal.notes

    # The uniform record stays, for the path that still ships.
    assert evidence["mismatch_penalty"].status == "assumed"


def test_the_modelled_background_path_asks_for_the_uniform_model():
    """A compositional profile has no neighbours to weight, so a
    position-dependent foreground over a uniform modelled background would be a
    ratio between two different models. Pinned at the source, because the two
    halves of that ratio are computed in different modules."""
    import ast
    import inspect

    from neoswga.core import base_optimizer

    source = inspect.getsource(base_optimizer.BaseOptimizer._modelled_site_load)
    tree = ast.parse(source.strip())
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "weighted_site_load"
    ]

    assert calls
    for call in calls:
        names = [getattr(arg, "id", None) for arg in call.args]
        assert "UNIFORM" in names, (
            "the modelled-background foreground load must request the uniform "
            "model explicitly, not inherit the configured one"
        )


def test_the_schema_offers_exactly_the_models_the_code_accepts():
    """A documented value the code refuses, or a value the code accepts and the
    schema rejects, are the same defect in opposite directions."""
    import pathlib

    schema = json.loads(
        (
            pathlib.Path(__file__).resolve().parent.parent
            / "neoswga"
            / "core"
            / "schema"
            / "params.schema.json"
        ).read_text()
    )

    assert tuple(schema["properties"]["mismatch_model"]["enum"]) == mm.MISMATCH_MODELS
    assert schema["properties"]["mismatch_model"]["default"] == mm.UNIFORM


def test_the_three_prime_factor_is_a_fraction():
    """A factor above 1 would make a terminal mismatch an advantage."""
    assert 0.0 <= mm.THREE_PRIME_EXTENSION_FACTOR < 1.0
    assert mm.THREE_PRIME_WINDOW >= 1
    assert not math.isnan(mm.THREE_PRIME_EXTENSION_FACTOR)


# ----------------------------------------------------------------------
# 5. Identity and provenance
# ----------------------------------------------------------------------


MINIMAL_PARAMS = {
    "schema_version": 2,
    "fg_genomes": ["fg.fna"],
    "fg_prefixes": ["fg"],
    "bg_genomes": [],
    "bg_prefixes": [],
    "data_dir": "./",
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "min_fg_freq": 1e-5,
    "max_bg_freq": 5e-6,
    "min_tm": 15,
    "max_tm": 45,
    "max_gini": 0.6,
    "max_primer": 500,
    "cpus": 1,
    "iterations": 8,
    "max_sets": 5,
    "fg_seq_lengths": [1_000_000],
}


def _request(**overrides):
    from neoswga.core.design_request import resolve_design_request

    return resolve_design_request({**MINIMAL_PARAMS, **overrides})


def test_the_model_is_part_of_the_requests_identity():
    """`request_hash` is documented as covering every setting that affects a
    design. The model changes how every background site is weighted, so two
    designs that deliver different panels recorded one identity without it."""
    uniform = _request().request_hash
    duplex = _request(mismatch_model="position-dependent").request_hash
    with_rule = _request(mismatch_model="position-dependent-3prime").request_hash

    assert len({uniform, duplex, with_rule}) == 3


def test_an_absent_key_and_an_explicit_uniform_are_one_identity():
    """They deliver byte-identical designs, and the hash is a design's identity.

    This follows `candidate_retention`, which resolves absent to `all_qc`,
    rather than `optimization_method`, which keeps None so a CLI flag can beat
    the file. There is no flag for this key, so there is no precedence to keep.
    """
    assert _request().request_hash == _request(mismatch_model="uniform").request_hash
    assert _request().mismatch_model == mm.UNIFORM


def test_an_unrecognised_model_is_refused_before_any_search():
    """`resolve_design_request` is where everything knowable from the request
    alone is refused, and this is knowable there."""
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest) as caught:
        _request(mismatch_model="position-aware")

    assert "mismatch_model" in str(caught.value)


def test_the_model_is_recorded_in_the_run_manifest(monkeypatch):
    """A manifest that does not name the model describes numbers it cannot
    account for. Recorded with the reaction, because it decides how every
    background site is weighted."""
    from neoswga.cli._common import _effective_conditions
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "mismatch_model", None, raising=False)
    assert _effective_conditions(parameter)["mismatch_model"] == mm.UNIFORM

    monkeypatch.setattr(parameter, "mismatch_model", "position-dependent", raising=False)
    assert _effective_conditions(parameter)["mismatch_model"] == "position-dependent"


def test_a_change_of_model_between_steps_is_condition_drift(tmp_path, monkeypatch):
    """`filter` under one model and `optimize` under another selected a pool
    under one weighting and scored it under a different one. That is exactly
    what the drift check exists for, and it said nothing before."""
    from neoswga.cli._common import _effective_conditions, warn_on_condition_drift
    from neoswga.core import parameter
    from neoswga.core.run_manifest import write_manifest

    monkeypatch.setattr(parameter, "data_dir", str(tmp_path), raising=False)
    monkeypatch.setattr(parameter, "mismatch_model", "position-dependent", raising=False)
    write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions=_effective_conditions(parameter),
    )

    monkeypatch.setattr(parameter, "mismatch_model", None, raising=False)

    assert "mismatch_model" in warn_on_condition_drift(parameter, reference_step="filter")


def test_the_extrapolations_size_is_stated_under_the_model_and_not_otherwise():
    """No threshold refuses a design on this distance, so stating it is the
    whole of what can be done about it. An assessment that does not state it
    reads as if the distance were zero."""
    assert mm.extrapolation_notice(mm.UNIFORM, 63.0) == ""
    assert mm.extrapolation_notice(None, 63.0) == ""

    hot = mm.extrapolation_notice(mm.POSITION_DEPENDENT, 63.0)
    assert "26.0 C above" in hot and "37 C" in hot
    assert mm.TABLE_REFERENCE_BUFFER in hot

    cool = mm.extrapolation_notice(mm.POSITION_DEPENDENT, 30.0)
    assert "7.0 C below" in cool

    unknown = mm.extrapolation_notice(mm.POSITION_DEPENDENT, None)
    assert "not recorded" in unknown


def test_the_notice_reaches_the_assessment_evidence():
    from neoswga.core.panel_evaluation import _mismatch_evidence

    hot = _mismatch_evidence("occupancy-position-dependent", 63.0)["mismatch_model"]
    assert "26.0 C above" in hot
    assert "terminus" in hot, "the domain limit must be stated where the model is described"

    assert "C above" not in _mismatch_evidence("occupancy", 63.0)["mismatch_penalty"]


# ----------------------------------------------------------------------
# 6. Every consumer of a weighted load says which model produced it
# ----------------------------------------------------------------------
#
# The schema says a non-uniform load "is reported with its own
# selectivity_mode so it cannot be read as the uniform one". That was true of
# the optimizer and of nothing else: five consumers picked the key up from the
# module global and recorded nothing, and `report/` read no mode at all. One
# test per consumer, and each asserts the uniform case as well, because "an
# added field and nothing else moves" is the other half of the promise.


@pytest.mark.parametrize(
    "model,expected",
    [
        (None, "occupancy"),
        ("uniform", "occupancy"),
        ("position-dependent", "occupancy-position-dependent"),
        ("position-dependent-3prime", "occupancy-position-dependent-3prime"),
    ],
)
def test_evaluate_set_reports_the_model(model, expected, monkeypatch):
    from neoswga.cli.evaluate import _load_mode
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "mismatch_model", model, raising=False)
    assert _load_mode() == expected


def test_the_evaluate_set_result_carries_the_mode_field():
    """Source-pinned: the figure and its mode are built in one dict literal,
    so a reader cannot get one without the other."""
    import inspect

    from neoswga.cli import evaluate

    source = inspect.getsource(evaluate)
    index = source.index('"occupancy_selectivity_ratio"')
    assert '"occupancy_selectivity_mode"' in source[index : index + 800]


@pytest.mark.parametrize(
    "model,expected",
    [(None, "occupancy"), ("position-dependent", "occupancy-position-dependent")],
)
def test_reference_panel_evaluation_publishes_the_model(model, expected, monkeypatch, table):
    """`weighted_site_load` is published beside measured site counts under one
    field name, so a modelled load that said nothing read as measured."""
    from neoswga.core import parameter
    from neoswga.core import reference_panel_evaluation as rpe

    monkeypatch.setattr(parameter, "mismatch_model", model, raising=False)

    record = rpe.ReferenceRecord(
        prefix=table,
        role="target",
        length=1_000_000,
        status="measured",
        source="table",
        sites=rpe.Measurement("sites", 1.0, "sites"),
        site_density=rpe.Measurement("sites_per_mb", 1.0, "sites/Mb"),
        coverage=rpe.Measurement("coverage", 1.0, "fraction"),
        mean_gap=rpe.Measurement("mean_gap_bp", 1.0, "bp"),
        max_gap=rpe.Measurement("max_gap_bp", 1.0, "bp"),
        gap_gini=rpe.Measurement("gap_gini", 0.1, "dimensionless"),
        weighted_load=rpe.Measurement("weighted_site_load", None, "sites", unavailable="pending"),
    )
    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    updated = rpe._with_weighted_load(record, [PRIMER], conditions, 1, rpe.PanelSources())

    published = updated.as_dict()
    assert published["weighted_site_load_mode"] == expected
    assert published["weighted_site_load"]["value"] is not None

    basis = published["weighted_site_load"]["basis"]
    if model is None:
        assert basis == "<= 1 mismatch(es) at 30.0 C", "the uniform basis string moved"
    else:
        assert "position-dependent" in basis and "C below" in basis


@pytest.mark.parametrize(
    "model,expected",
    [(None, "occupancy"), ("position-dependent", "occupancy-position-dependent")],
)
def test_the_condition_sweep_records_the_model(model, expected, monkeypatch, table):
    """A sweep is read as "this reaction buys this much discrimination", and
    how much it reports depends on the model."""
    from neoswga.core import parameter
    from neoswga.core.condition_sweep import sweep_conditions

    monkeypatch.setattr(parameter, "mismatch_model", model, raising=False)
    result = sweep_conditions([PRIMER], [table], [table], polymerase="phi29", temperatures=(30.0,))
    points = result.points if hasattr(result, "points") else result

    assert points
    for point in points:
        assert point.selectivity_mode == expected


@pytest.mark.parametrize(
    "model,expected",
    [(None, "occupancy"), ("position-dependent", "occupancy-position-dependent")],
)
def test_the_filter_records_which_model_ordered_the_pool(model, expected, monkeypatch):
    """The ranking's output is the ORDER of step2_df.csv, so the file that
    describes the funnel is where the model belongs."""
    import inspect

    from neoswga.core import parameter, pipeline
    from neoswga.core.mismatch_model import resolve_mismatch_model, site_load_mode

    monkeypatch.setattr(parameter, "mismatch_model", model, raising=False)
    assert site_load_mode(resolve_mismatch_model()) == expected

    source = inspect.getsource(pipeline._write_filter_stats)
    assert (
        '"occupancy_ranking_mode"' in source
    ), "filter_stats.json must name the model that ordered the candidates"
    assert "_write_filter_stats(" in inspect.getsource(pipeline.step2), (
        "the helper must still be called from step2, or filter_stats.json is " "not written at all"
    )


def test_every_key_the_filter_writes_is_a_field_the_report_reads():
    """`FilteringStats(**stats)` is built straight from `filter_stats.json`.

    So a key the filter writes and the dataclass does not declare is a
    TypeError in `report`, not an ignored field. Adding
    `occupancy_ranking_mode` to the file broke the funnel test exactly that
    way; this pins the coupling so the next key does not.
    """
    import dataclasses

    from neoswga.core.report.metrics import FilteringStats

    declared = {f.name for f in dataclasses.fields(FilteringStats)}

    assert "occupancy_ranking_mode" in declared
    # Not a stage count, so it must never appear in the rendered funnel.
    stats = FilteringStats(total_kmers=10, final_candidates=2)
    assert "occupancy_ranking_mode" not in {label for label, _ in stats.as_funnel()}


def test_improve_set_keys_its_load_cache_on_the_model():
    """A cache keyed on everything but the model answers with the other
    model's number the moment the key changes in one process."""
    import inspect

    from neoswga.core.set_improvement import MemoisedSources

    source = inspect.getsource(MemoisedSources.weighted_load)
    assert "resolve_mismatch_model" in source


@pytest.mark.parametrize(
    "mode,shown",
    [("occupancy", False), ("exact", False), ("occupancy-position-dependent", True)],
)
def test_the_report_reads_and_shows_the_mode(mode, shown):
    """`grep selectivity_mode neoswga/core/report/` used to return nothing, so
    a figure from an extrapolated table rendered identically to a uniform one."""
    from neoswga.core.report.metrics import SpecificityMetrics

    assert SpecificityMetrics().selectivity_mode == "occupancy"

    import inspect

    from neoswga.core.report import metrics as report_metrics
    from neoswga.core.report import technical_report

    assert "selectivity_mode" in inspect.getsource(report_metrics)

    narrative = inspect.getsource(technical_report.render_technical_report)
    assert "selectivity_mode" in narrative
    assert (
        "occupancy-position-dependent" in narrative
    ), "the report must name the model it is warning about"
    del mode, shown


def test_the_notice_is_logged_once_per_run(caplog, table, conditions):
    """Warning level, because a number produced outside its source's domain
    must not pass silently; once, because a search evaluates thousands of
    panels."""
    import logging

    mm._LOGGED.clear()
    with caplog.at_level(logging.WARNING, logger="neoswga.core.mismatch_model"):
        for _ in range(3):
            occupancy.weighted_site_load(
                [PRIMER], [table], conditions, 1, 4.0, mm.POSITION_DEPENDENT
            )
    warnings = [r for r in caplog.records if "mismatch_model=" in r.message]
    assert len(warnings) == 1, [r.message for r in warnings]

    mm._LOGGED.clear()
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="neoswga.core.mismatch_model"):
        occupancy.weighted_site_load([PRIMER], [table], conditions, 1, 4.0, mm.UNIFORM)
    assert not [r for r in caplog.records if "mismatch_model=" in r.message]
