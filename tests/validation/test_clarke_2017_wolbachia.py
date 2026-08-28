"""Third set-level benchmark: Clarke et al. (2017) Wolbachia sets.

Data: Clarke et al. (2017) Bioinformatics 33(14):2071-2077, Table 1 /
Supplementary Table 1, with outcomes quoted from Results section 3.1. Five sets
against *W. pipientis* in a *D. melanogaster* background, including the original
hand-designed set from Leichty and Brisson (2014).

This benchmark is unusually clean because the authors built it as a controlled
contrast: within each of two primer melting-temperature ranges they took the
*most selective* set and the *most even* set. That directly pits the two design
criteria against each other.

The result favours evenness. `TmL/Even` achieved 10x coverage over 60-75% of the
genome against 2.8% unamplified, and it did so with **worse** mean spacing
(6853 bp) than either `TmL/Selective` (5327 bp) or the Leichty set (5305 bp) --
both of which "failed to improve sequencing efficiency due an unevenness of
coverage".

Combined with the other two set-level benchmarks:

  Prevotella (Dwivedi-Yu)  evenness / worst gap predict; density does not
  Wolbachia  (Clarke)      evenness predicts; density runs backwards
  M. tuberculosis (Clarke) density predicts; evenness runs backwards

The odd one out is *M. tuberculosis*, and the paper says why: there the primer
pool was pre-filtered for even binding, which removed most of the variance in
evenness and left density as the limiting property. Two of three datasets favour
evenness; all three are consistent with "whichever property is limiting is the
one that predicts".
"""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "clarke_2017_wolbachia.json"


@pytest.fixture(scope="module")
def sets():
    return json.loads(DATA.read_text())["sets"]


def test_fixture_is_intact(sets):
    assert len(sets) == 5
    for name, s in sets.items():
        assert len(s["primers"]) == s["published"]["size"], name
        assert all(set(p) <= set("ACGT") for p in s["primers"]), name


def test_exactly_one_set_was_effective(sets):
    effective = {n for n, v in sets.items() if v["outcome"]["effective"]}
    assert effective == {"TmL/Even"}


def test_the_winner_was_the_one_selected_for_evenness(sets):
    """The design directly contrasts evenness against selectivity."""
    winner = next(n for n, v in sets.items() if v["outcome"]["effective"])
    assert sets[winner]["outcome"]["selected_for"] == "evenness"
    assert sets[winner]["outcome"]["tm_range"] == "low"


def test_the_winner_has_worse_density_than_sets_it_beat(sets):
    """Density is not what separated them - this is the point of the benchmark.

    TmL/Even won with 6853 bp mean spacing while TmL/Selective (5327) and the
    hand-designed Leichty set (5305) were denser and did not.
    """
    winner = sets["TmL/Even"]["published"]["fg_mean_distance"]
    for beaten in ("TmL/Selective", "Leichty and Brisson 2014"):
        assert sets[beaten]["published"]["fg_mean_distance"] < winner, beaten


def test_evenness_separates_the_two_low_tm_sets(sets):
    """Holding Tm range fixed, Gini is what distinguishes them."""
    even = sets["TmL/Even"]["published"]["fg_gini"]
    selective = sets["TmL/Selective"]["published"]["fg_gini"]
    assert even < selective
    assert sets["TmL/Even"]["outcome"]["effective"]
    assert not sets["TmL/Selective"]["outcome"]["effective"]


def test_high_tm_sets_failed_regardless_of_selection_criterion(sets):
    """A second, independent failure mode: too few binding sites.

    Both TmH sets failed, whether chosen for evenness or selectivity, and both
    have mean spacing roughly twice the TmL sets and max gaps over 110 kb. Tm
    range set the number of usable primers before either criterion could apply.
    """
    for name in ("TmH/Even", "TmH/Selective"):
        assert sets[name]["outcome"]["tier"] == "ineffective"
        assert sets[name]["published"]["fg_max_distance"] > 100_000
        assert sets[name]["published"]["fg_mean_distance"] > 10_000


def test_hand_designed_set_underperformed_the_tool(sets):
    """Leichty and Brisson (2014), included by the authors as a comparison."""
    leichty = sets["Leichty and Brisson 2014"]
    assert not leichty["outcome"]["effective"]
    # It has the worst evenness of any set here despite competitive density.
    assert leichty["published"]["fg_gini"] == max(s["published"]["fg_gini"] for s in sets.values())


# ----------------------------------------------------------------------
# Cross-dataset synthesis
# ----------------------------------------------------------------------


def test_two_of_three_benchmarks_favour_evenness_over_density(sets):
    """The combined picture across all three set-level datasets.

    Guards the synthesis in docs/validation/: evenness wins where it varies,
    density wins where evenness was already controlled. If a future dataset
    breaks this, the guidance in the swga-chemistry skill needs revisiting.
    """
    here = Path(__file__).parent / "data"
    prev = json.loads((here / "dwivedi_yu_2023_prevotella.json").read_text())["sets"]
    mtb = json.loads((here / "clarke_2017_mtb.json").read_text())["sets"]

    # Wolbachia: winner is less dense than sets it beat.
    wol_density_backwards = sets["TmL/Even"]["published"]["fg_mean_distance"] > min(
        sets[n]["published"]["fg_mean_distance"]
        for n in sets
        if not sets[n]["outcome"]["effective"]
    )

    # Prevotella: winners are not the densest either.
    prev_winners = {"Prev03", "Prev06"}
    prev_density_backwards = min(
        prev[n]["published"]["fg_mean_distance"] for n in prev_winners
    ) > min(prev[n]["published"]["fg_mean_distance"] for n in prev if n not in prev_winners)

    # M. tuberculosis: winners ARE the densest.
    mtb_winners = {n for n, v in mtb.items() if v["outcome"]["most_effective"]}
    mtb_density_forward = max(mtb[n]["published"]["fg_mean_distance"] for n in mtb_winners) < min(
        mtb[n]["published"]["fg_mean_distance"] for n in mtb if n not in mtb_winners
    )

    assert wol_density_backwards and prev_density_backwards
    assert mtb_density_forward, (
        "M. tuberculosis is the dataset where density does predict, because its "
        "primer pool was pre-filtered for even binding"
    )
