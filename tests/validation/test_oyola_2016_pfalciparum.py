"""Set-level benchmark: Oyola et al. (2016) Plasmodium falciparum from blood spots.

Data: *Malaria Journal* 15:597, Table S1 (Additional file 1) for the sequences,
Results for the outcomes. One published set, `Probe_10`, used on 156 patients:
>18-fold enrichment, 73.2 % of reads on target against <1 % for standard WGA.

Two reasons this dataset earns its place, and they are unrelated.

**It replicates the nested-set result independently.** Three nested pools were
compared - `Probe_10` (first 10 primers), `Probe_20` (first 20), `Probe_28` (all
28) - and the smallest won, on yield and cost. Leichty and Brisson (2014) found
the same thing with SR2 beating the SR1 it was drawn from. Two independent groups,
different organisms, same direction: adding primers to a working pool made it
worse. Only `Probe_10`'s sequences were published, so the losing pools cannot be
scored here; the structure is recorded, not measured.

**It is a set this tool cannot design.** *P. falciparum* is about 19 % GC and the
working primers are pure A/T - one is twelve consecutive adenines. Our GC clamp
requires at least one G or C in the final five bases, which no zero-GC primer can
satisfy, so all ten are rejected. The adaptive GC filter does not rescue them
either: it floors its lower bound at 0.20, which for a 0.19-GC genome excludes the
genome's own average composition.

That is the most consequential thing to come out of the validation suite. The GC
clamp is standard PCR primer-design doctrine, and it is imported into a regime it
does not fit: SWGA against an AT-rich target at 30 C, where AT-rich primers are
the entire point. The tests below pin the current behaviour rather than assert the
desired behaviour, so that a deliberate change to the filters shows up here as a
failure to be reviewed rather than passing silently.
"""

import json
import statistics
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "oyola_2016_pfalciparum.json"
LEICHTY = Path(__file__).parent / "data" / "leichty_brisson_2014.json"

PFALCIPARUM_GC = 0.19


@pytest.fixture(scope="module")
def doc():
    return json.loads(DATA.read_text())


@pytest.fixture(scope="module")
def probe10(doc):
    return doc["sets"]["Probe_10"]["primers"]


# ----------------------------------------------------------------------
# Fixture integrity
# ----------------------------------------------------------------------


def test_fixture_is_intact(doc, probe10):
    assert len(probe10) == 10
    for p in probe10:
        assert set(p) <= set("ACGT"), p
        assert 8 <= len(p) <= 12, p


def test_phosphorothioate_bonds_were_stripped(doc, probe10):
    """The source table writes one 3' PTO bond as '*'; sequences here are bare.

    Worth pinning because the modification is real chemistry, not notation - it
    protects the 3' end from phi29's exonuclease. `neoswga.core.export` adds it
    back on ordering sheets (two bonds by default, where this paper used one).
    """
    assert "phosphorothioate" in doc["source"]["modification"]
    assert all("*" not in p for p in probe10)


# ----------------------------------------------------------------------
# The nested-pool result, replicated
# ----------------------------------------------------------------------


def test_the_smallest_of_three_nested_pools_was_chosen(doc):
    outcome = doc["sets"]["Probe_10"]["outcome"]
    assert outcome["effective"]
    assert set(outcome["chosen_over"]) == {"Probe_20", "Probe_28"}


def test_two_independent_papers_agree_that_smaller_nested_sets_won():
    """Oyola's Probe_10 over Probe_20/28, and Leichty's SR2 over SR1.

    Different organisms, different groups, different decades of the method. Both
    found that the smaller of two nested pools performed better, which is the
    opposite of what "more primers means more coverage" would predict. Neither is
    conclusive alone; together they are worth not designing against.
    """
    oyola = json.loads(DATA.read_text())["sets"]["Probe_10"]["outcome"]
    leichty = json.loads(LEICHTY.read_text())["sets"]["SR2"]["outcome"]

    assert oyola["chosen_over"], "Oyola: smaller pool preferred over larger nested pools"
    assert leichty["more_selective_than"] == "SR1"


# ----------------------------------------------------------------------
# The set this tool cannot design
# ----------------------------------------------------------------------


def test_the_working_primers_are_pure_at(probe10):
    from neoswga.core import thermodynamics as thermo

    assert all(thermo.gc_content(p) == 0.0 for p in probe10)
    assert any(p == "A" * len(p) for p in probe10), "expected a pure homopolymer"


def test_gc_clamp_rejects_the_entire_published_set(probe10):
    """Current behaviour, pinned deliberately.

    The clamp requires 1-3 G/C in the last five bases. A zero-GC primer cannot
    satisfy it at any threshold, so this is structural rather than a matter of
    tuning. If the clamp is ever made conditional on target GC, this test should
    fail and be updated as part of that review.
    """
    from neoswga.core.adaptive_filters import GCClampFilter

    clamp = GCClampFilter()
    assert sum(clamp.passes(p) for p in probe10) == 0


def test_adaptive_gc_filter_excludes_the_genome_it_adapts_to(probe10):
    """The lower bound floors at 0.20, above *P. falciparum*'s own 0.19 GC.

    An adaptive filter centred on the target should admit sequences at the
    target's composition. Here it cannot, so the adaptive path does not rescue
    this set either.
    """
    from neoswga.core.adaptive_filters import AdaptiveGCFilter

    adaptive = AdaptiveGCFilter(PFALCIPARUM_GC)
    assert adaptive.gc_min > PFALCIPARUM_GC
    assert sum(adaptive.passes(p) for p in probe10) == 0


def test_only_the_at_rich_targets_are_affected():
    """Scopes the finding: this is not a general filter problem.

    Every other published set in the suite passes the GC clamp almost entirely.
    The two that fail outright are both against very AT-rich targets - this one
    and Leichty's *Borrelia* set - which is what makes the clamp's unconditional
    application the thing to question rather than the clamp itself.
    """
    from neoswga.core.adaptive_filters import GCClampFilter

    clamp = GCClampFilter()
    here = Path(__file__).parent / "data"

    rejected, retained = [], []
    for path in sorted(here.glob("*.json")):
        doc = json.loads(path.read_text())
        for name, s in doc.get("sets", {}).items():
            primers = s["primers"]
            passing = sum(clamp.passes(p) for p in primers)
            label = f"{path.stem}:{name}"
            (rejected if passing == 0 else retained).append(label)

    assert set(rejected) == {
        "leichty_brisson_2014:B31-BL21",
        "oyola_2016_pfalciparum:Probe_10",
    }, f"the set of fully-rejected published sets changed: {sorted(rejected)}"

    # And the rest are not marginal - they pass overwhelmingly.
    assert len(retained) > 25


def test_rejected_sets_are_the_at_rich_ones():
    from neoswga.core import thermodynamics as thermo

    oyola = json.loads(DATA.read_text())["sets"]["Probe_10"]["primers"]
    borrelia = json.loads(LEICHTY.read_text())["sets"]["B31-BL21"]["primers"]

    for primers in (oyola, borrelia):
        assert statistics.mean(thermo.gc_content(p) for p in primers) < 0.15
