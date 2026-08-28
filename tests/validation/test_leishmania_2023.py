"""Set-level benchmark: Dwivedi-Yu et al. (2023) Leishmania braziliensis.

Data: PLOS Pathogens 19(3):e1011230, S1 Table (both sheets) for the sets, Results
for the outcomes. Four sets of ten 8-mers drawn from a common pool of 23, two of
which worked and two of which did not.

Three things make this one worth carrying.

**It is the closest thing in the suite to a controlled A/B.** PS1 and PS2 are
*entirely disjoint* - twenty primers, no overlap - and PS3/PS4 are hybrids of the
two. PS1 won and PS2 lost; PS4, which takes six of its ten from PS1, won; PS3,
which takes six from PS2, lost. Membership in the PS1 pool tracks the outcome
perfectly, which is a much cleaner signal than the usual "these four of twelve
worked".

**Sequence composition does not explain it, and runs backwards.** The winners are
not distinguished by mean Tm, Tm spread, GC content or 3' clamp. PS4 won with the
*widest* Tm spread and the *lowest* GC of the four. Since the composite score
carries a Tm-uniformity term, this is a direct counter-example to the premise
that uniform Tm within a set is what makes it work, and it is recorded as such.

**It is the only multi-background set here** (human plus two skin commensals) and
the only 8-mer set against a GC-rich target, so it also guards against the suite
quietly becoming a collection of AT-rich bacterial designs.

What this dataset cannot test is binding geometry, which is where the other
benchmarks locate the signal - that would need the 32 Mb *L. braziliensis*
assembly, which is too large to carry as a fixture. The honest reading is
therefore narrow: whatever separated these sets, it was not primer composition.
"""

import json
import statistics
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "dwivedi_yu_2023_leishmania.json"


@pytest.fixture(scope="module")
def doc():
    return json.loads(DATA.read_text())


@pytest.fixture(scope="module")
def sets(doc):
    return doc["sets"]


@pytest.fixture(scope="module")
def winners(sets):
    return {n for n, v in sets.items() if v["outcome"]["effective"]}


# ----------------------------------------------------------------------
# Fixture integrity
# ----------------------------------------------------------------------


def test_fixture_is_intact(doc, sets):
    assert set(sets) == {"PS1", "PS2", "PS3", "PS4"}
    pool = set(doc["candidate_pool"])
    assert len(pool) == 23
    for name, s in sets.items():
        assert len(s["primers"]) == 10, name
        assert all(set(p) <= set("ACGT") for p in s["primers"]), name
        assert all(len(p) == 8 for p in s["primers"]), name
        assert set(s["primers"]) <= pool, name


def test_two_sets_worked_and_two_did_not(winners):
    assert winners == {"PS1", "PS4"}


# ----------------------------------------------------------------------
# The near-controlled contrast
# ----------------------------------------------------------------------


def test_the_two_parent_sets_are_disjoint(sets):
    """PS1 and PS2 share no primers, which is what makes the contrast clean."""
    assert not set(sets["PS1"]["primers"]) & set(sets["PS2"]["primers"])


def test_outcome_tracks_membership_in_the_winning_pool(sets, winners):
    """Every set that drew mostly from PS1 worked; every set that drew mostly
    from PS2 did not. With PS1 and PS2 disjoint this is a real split rather than
    a restatement of the labels.
    """
    ps1 = set(sets["PS1"]["primers"])
    ps2 = set(sets["PS2"]["primers"])

    for name, s in sets.items():
        primers = set(s["primers"])
        from_ps1, from_ps2 = len(primers & ps1), len(primers & ps2)
        if name in winners:
            assert from_ps1 >= 6 and from_ps2 <= 2, name
        else:
            assert from_ps2 >= 6 and from_ps1 <= 2, name


# ----------------------------------------------------------------------
# Composition does not separate them
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def composition(sets):
    from neoswga.core import thermodynamics as thermo
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    out = {}
    for name, s in sets.items():
        tms = [conditions.calculate_effective_tm(p) for p in s["primers"]]
        out[name] = {
            "tm_mean": statistics.mean(tms),
            "tm_spread": statistics.pstdev(tms),
            "gc_mean": statistics.mean(thermo.gc_content(p) for p in s["primers"]),
        }
    return out


@pytest.mark.parametrize("metric", ["tm_mean", "tm_spread", "gc_mean"])
def test_no_composition_metric_separates_the_winners(composition, winners, metric):
    """None of these draws a line between the two that worked and the two that did not."""
    w = [composition[n][metric] for n in winners]
    o = [composition[n][metric] for n in composition if n not in winners]
    assert not (max(w) < min(o) or min(w) > max(o)), (
        f"{metric} unexpectedly separates the Leishmania sets; if this is real "
        "the conclusion that composition does not explain the outcome needs revisiting"
    )


def test_the_winner_has_the_worst_tm_uniformity(composition, winners):
    """A direct counter-example to weighting Tm uniformity.

    PS4 worked while having the widest Tm spread of all four sets. The composite
    score rewards uniform Tm within a set; on this dataset that term points the
    wrong way. One counter-example is not grounds for removing it, but it is
    grounds for not treating it as established.
    """
    widest = max(composition, key=lambda n: composition[n]["tm_spread"])
    assert widest in winners


# ----------------------------------------------------------------------
# Coverage of the suite
# ----------------------------------------------------------------------


def test_this_is_the_only_multi_background_design(doc):
    """Guards a gap the other benchmarks share: all of them use one background."""
    assert len(doc["source"]["backgrounds"]) > 1

    here = Path(__file__).parent / "data"
    others = [json.loads(p.read_text()) for p in here.glob("*.json") if p.name != DATA.name]
    assert all(
        len(o.get("source", {}).get("backgrounds", [])) <= 1 for o in others
    ), "another multi-background dataset has been added; this test can be relaxed"


def test_target_is_gc_rich_unlike_the_other_benchmarks(sets):
    """The other set-level benchmarks target AT-rich genomes.

    Borrelia primers are near-zero GC; these average around 60 %, so the suite
    covers both ends rather than only the AT-rich case that SWGA is easiest on.
    """
    from neoswga.core import thermodynamics as thermo

    gc = [thermo.gc_content(p) for s in sets.values() for p in s["primers"]]
    assert statistics.mean(gc) > 0.5
