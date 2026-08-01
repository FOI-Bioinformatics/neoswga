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

**It was a set this tool could not design, and it drove a fix.** *P. falciparum*
is about 19 % GC and the working primers are pure A/T - one is twelve consecutive
adenines. Every GC rule rejected them, in two different ways:

  - the GC *clamp* in `adaptive_filters` required at least one G or C in the
    final five bases, which no zero-GC primer can satisfy at any threshold;
  - the adaptive GC *window* floored its lower bound at a positive value (0.20 in
    `adaptive_filters`, `max(0.15, genome_gc - tol)` in `parameter`), so for a
    0.19-GC genome the "adapted" window excluded the target's own composition.

The second was the binding constraint in the production path, where the clamp in
`filter.filter_extra` was already conditional on target GC. Both sides now branch
on the same shared thresholds, and the published set passes both GC rules.

The GC clamp is standard PCR primer-design practice imported into a regime it does
not fit: SWGA against an AT-rich target at 30 C, where AT-rich primers are the
entire point. The fix is deliberately narrow - it moves only compositionally
extreme targets, and tests below pin that GC-normal behaviour is unchanged.

The fix is also **necessary but not sufficient**: `filter_extra` still keeps only
3 of the 10, the rest held back by the homopolymer and self-dimer rules, which
this change did not touch. That boundary is pinned rather than moved, since
relaxing either is a separate question with its own evidence.
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


def test_gc_clamp_is_waived_on_at_rich_targets(probe10):
    """A zero-GC primer cannot satisfy a GC clamp at any threshold.

    The clamp is PCR primer-design practice and does not transfer to SWGA
    against an AT-rich genome. Told the target GC, it now waives the minimum
    and admits the published set.
    """
    from neoswga.core.adaptive_filters import GCClampFilter

    informed = GCClampFilter(genome_gc=PFALCIPARUM_GC)
    assert sum(informed.passes(p) for p in probe10) == len(probe10)


def test_gc_clamp_is_unchanged_when_the_target_is_not_at_rich(probe10):
    """The waiver is narrow: it applies only to compositionally extreme targets."""
    from neoswga.core.adaptive_filters import GCClampFilter

    for gc in (0.40, 0.50, 0.65):
        assert GCClampFilter(genome_gc=gc).min_gc == 1
        assert sum(GCClampFilter(genome_gc=gc).passes(p) for p in probe10) == 0


def test_adaptive_gc_filter_admits_the_genome_it_adapts_to(probe10):
    """The old lower bound floored at 0.20, above *P. falciparum*'s own 0.19 GC.

    An adaptive filter that cannot admit sequences at the target's composition
    is not adapting. On extreme targets the bound on the matching side is now
    released rather than centred on the genome mean.
    """
    from neoswga.core.adaptive_filters import AdaptiveGCFilter

    adaptive = AdaptiveGCFilter(PFALCIPARUM_GC)
    assert adaptive.gc_min <= PFALCIPARUM_GC
    assert sum(adaptive.passes(p) for p in probe10) == len(probe10)


def test_adaptive_gc_window_is_unchanged_for_normal_targets():
    """Regression guard: only extreme targets move."""
    from neoswga.core.adaptive_filters import AdaptiveGCFilter

    normal = AdaptiveGCFilter(0.50)
    assert (normal.gc_min, normal.gc_max) == (0.35, 0.65)


def test_the_gc_rules_no_longer_reject_any_published_set():
    """Scopes the fix: it rescues the two AT-rich sets and disturbs nothing else.

    Before the change, `oyola_2016:Probe_10` and `leichty_brisson_2014:B31-BL21`
    were rejected outright by the GC clamp while every other published set passed
    almost completely. Given each set's own target GC, none is now rejected.
    """
    from neoswga.core import thermodynamics as thermo
    from neoswga.core.adaptive_filters import GCClampFilter

    here = Path(__file__).parent / "data"
    rejected, checked = [], 0
    for path in sorted(here.glob("*.json")):
        doc = json.loads(path.read_text())
        for name, s in doc.get("sets", {}).items():
            primers = s["primers"]
            # Stand in for the target's composition with the set's own mean GC;
            # SWGA primers track their target closely enough for this purpose.
            set_gc = statistics.mean(thermo.gc_content(p) for p in primers)
            clamp = GCClampFilter(genome_gc=set_gc)
            checked += 1
            if not any(clamp.passes(p) for p in primers):
                rejected.append(f"{path.stem}:{name}")

    assert rejected == [], f"published sets still rejected outright: {rejected}"
    assert checked > 25


def test_other_filter_rules_still_reject_most_of_this_set(probe10):
    """The GC fix is necessary but not sufficient, and that is worth pinning.

    With the GC window and clamp corrected, `filter_extra` keeps 3 of the 10.
    The rest are held back by rules untouched by this change:

      - the homopolymer rule rejects the twelve-adenine primer;
      - the self-dimer rule rejects the pure (AT)n primers, which are genuinely
        self-complementary -- though at 30 C a 10 bp AT duplex melts well below
        the reaction temperature, so counting paired bases without a
        thermodynamic threshold over-rejects in this regime.

    Whether to relax either for AT-rich targets is a separate question with its
    own evidence, deliberately not decided here. This test records where the
    boundary currently sits so that answering it later is a visible change.
    """
    from neoswga.core import dimer
    from neoswga.core import filter as filter_module
    from neoswga.core.filter import _has_homopolymer_run

    homopolymer = [p for p in probe10 if _has_homopolymer_run(p)]
    self_dimer = [p for p in probe10 if dimer.is_dimer_fast(p, p, 4)]

    assert homopolymer == ["A" * 12]
    assert len(self_dimer) == 6
    assert set(homopolymer) & set(self_dimer) == set()

    # 10 primers, minus 1 homopolymer, minus 6 self-dimering, leaves 3.
    assert len(probe10) - len(homopolymer) - len(self_dimer) == 3
    assert filter_module is not None


def test_rejected_sets_are_the_at_rich_ones():
    from neoswga.core import thermodynamics as thermo

    oyola = json.loads(DATA.read_text())["sets"]["Probe_10"]["primers"]
    borrelia = json.loads(LEICHTY.read_text())["sets"]["B31-BL21"]["primers"]

    for primers in (oyola, borrelia):
        assert statistics.mean(thermo.gc_content(p) for p in primers) < 0.15


def test_production_gc_window_reaches_zero_on_at_rich_targets():
    """The rule that actually governs `neoswga filter`.

    `parameter.get_params` derives the primer GC window from the foreground
    genome. Its old lower bound of `max(0.15, genome_gc - tolerance)` was the
    binding constraint on this set in the production path -- the GC clamp there
    was already conditional. Both sides are now consistent.
    """
    from neoswga.core.parameter import EXTREME_AT_GENOME_GC, EXTREME_GC_GENOME_GC

    assert EXTREME_AT_GENOME_GC == 0.30
    assert EXTREME_GC_GENOME_GC == 0.70
    assert PFALCIPARUM_GC < EXTREME_AT_GENOME_GC


def test_filter_and_adaptive_paths_use_the_same_thresholds():
    """Two GC-adaptive implementations exist; they must not drift apart.

    `filter.filter_extra` and `adaptive_filters` both branch on target GC. They
    now read the same constants rather than each carrying its own literal, which
    is how the two ended up disagreeing in the first place.
    """
    import inspect

    from neoswga.core import adaptive_filters
    from neoswga.core import filter as filter_module

    for module in (filter_module, adaptive_filters):
        src = inspect.getsource(module)
        assert "EXTREME_AT_GENOME_GC" in src, module.__name__
