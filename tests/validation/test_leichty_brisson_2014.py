"""The original SWGA paper: Leichty and Brisson (2014) Genetics 198:473.

Three primer sets with measured outcomes: `B31-BL21` (20 primers, *Borrelia
burgdorferi* against *E. coli*), and `SR1` / `SR2` (10 and 2 primers,
*Wolbachia pipientis* against *Drosophila melanogaster*).

Two things make this worth carrying alongside the later benchmarks.

**It is a nested-set experiment.** SR2 is literally the two least-permissive
primers of SR1 - same design, one a subset of the other. SR2 was *more*
selective: with restriction digest it reached 1264-fold target amplification
against 9-fold background (~140x) and up to 70 % of reads on target, against
2.4-8.8 % unamplified. Fewer, more selective primers beat more primers on
enrichment.

**It documents the other half of a trade-off.** Clarke et al. (2017) later
benchmarked this same SR2 set as "Leichty and Brisson 2014" and found it gave
*uneven* coverage and did not improve sequencing efficiency - despite having the
best binding density of the five sets they compared. Both results are correct
and they are measuring different things: SR2 maximises enrichment and sacrifices
evenness.

That trade-off is invisible from either paper alone, and it is the clearest
evidence in the whole benchmark suite that a single composite score cannot serve
both goals. Which end of it you want depends on whether you need target DNA or
even genome coverage.

Note evenness was never a design criterion here - the selection script scored
target/background k-mer frequency ratio with a Tm cutoff, predating the swga
toolkit's Gini index entirely.
"""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "leichty_brisson_2014.json"
CLARKE_WOL = Path(__file__).parent / "data" / "clarke_2017_wolbachia.json"


@pytest.fixture(scope="module")
def sets():
    return json.loads(DATA.read_text())["sets"]


def test_fixture_is_intact(sets):
    assert set(sets) == {"B31-BL21", "SR1", "SR2"}
    for name, s in sets.items():
        assert s["primers"], name
        assert all(set(p) <= set("ACGT") for p in s["primers"]), name


def test_set_sizes_match_the_paper(sets):
    assert len(sets["B31-BL21"]["primers"]) == 20
    assert len(sets["SR1"]["primers"]) == 10
    assert len(sets["SR2"]["primers"]) == 2


def test_sr2_is_a_strict_subset_of_sr1(sets):
    """The nested design is what makes the size comparison interpretable."""
    assert set(sets["SR2"]["primers"]) < set(sets["SR1"]["primers"])


def test_the_smaller_nested_set_was_more_selective(sets):
    """Fewer, more selective primers beat more primers on enrichment."""
    assert sets["SR2"]["outcome"]["more_selective_than"] == "SR1"
    assert sets["SR2"]["outcome"]["selectivity"] == pytest.approx(140.0, rel=0.1)
    assert sets["SR2"]["outcome"]["fold_target"] > sets["SR2"]["outcome"]["fold_background"]


def test_borrelia_primers_are_extremely_at_rich(sets):
    """A sanity check on the fixture, and on why this target is tractable.

    B. burgdorferi is very AT-rich relative to E. coli, so the selective motifs
    are near-homopolymeric runs of A and T. Their melting temperatures are all
    around 16-22 C, which only works because phi29 runs at 30 C.
    """
    gc = [(p.count("G") + p.count("C")) / len(p) for p in sets["B31-BL21"]["primers"]]
    assert max(gc) < 0.15, "expected near-zero GC primers for the Borrelia set"


# ----------------------------------------------------------------------
# Cross-referencing the same set between two papers
# ----------------------------------------------------------------------


def test_sr2_is_the_set_clarke_benchmarked(sets):
    """Clarke's 'Leichty and Brisson 2014' row is this set."""
    clarke = json.loads(CLARKE_WOL.read_text())["sets"]
    assert set(clarke["Leichty and Brisson 2014"]["primers"]) == set(sets["SR2"]["primers"])


def test_the_same_set_is_selective_but_uneven(sets):
    """The trade-off, visible only by combining the two papers.

    Leichty measured SR2 as highly selective (~140x). Clarke measured the same
    set as having the worst evenness of the five Wolbachia sets they compared,
    failing to improve sequencing efficiency despite competitive density.

    Neither is wrong. Enrichment and even coverage are different objectives, and
    this set trades the second for the first.
    """
    clarke = json.loads(CLARKE_WOL.read_text())["sets"]
    leichty_in_clarke = clarke["Leichty and Brisson 2014"]

    # Selective, per Leichty.
    assert sets["SR2"]["outcome"]["selectivity"] > 100

    # Yet the worst evenness of the sets Clarke compared, and not effective there.
    assert leichty_in_clarke["published"]["fg_gini"] == max(
        s["published"]["fg_gini"] for s in clarke.values()
    )
    assert not leichty_in_clarke["outcome"]["effective"]


def test_evenness_was_not_a_design_criterion_here():
    """Context for why this set is uneven: nothing selected for evenness."""
    doc = json.loads(DATA.read_text())
    assert "evenness was never a design criterion" in doc["source"]["notes"]
