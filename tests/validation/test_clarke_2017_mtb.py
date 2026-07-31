"""Second set-level benchmark: Clarke et al. (2017) M. tuberculosis primer sets.

Data: Clarke EL, Sundararaman SA, Seifert SN, Bushman FD, Hahn BH, Brisson D
(2017) Bioinformatics 33(14):2071-2077 -- Table 2 (set characteristics, with the
four most effective sets footnoted) and Supplementary Table 2 (primer
sequences). Twelve sets against *M. tuberculosis* in a human background,
including two chosen deliberately as negative controls.

**Why this matters: it disagrees with the Prevotella benchmark.**

On the Dwivedi-Yu *Prevotella* data, mean binding distance did not predict
success (rho -0.09) while the worst coverage hole and Gini did. Here the
opposite holds: the four most effective sets are the four *densest*, and Gini
runs the wrong way -- the winners have slightly *higher* Gini than the sets they
beat.

The paper explains it. Clarke pre-filtered the primer pool for even binding
before assembling sets, "constrain[ing] the range of set binding evenness ... by
removing primers that cluster on repeat regions". With evenness already
controlled, it had little residual variance left to explain the outcome, and
density dominated. Dwivedi-Yu did not control evenness, so it varied across
their sets and coverage holes dominated instead.

So neither "density predicts success" nor "evenness predicts success" is
universally true. What predicts is whichever is currently *limiting*. That is a
more useful rule than either single-dataset conclusion, and it is only visible
because there are now two benchmarks that disagree.

A note on `fg_bg_ratio`: it is defined identically in both papers
(fg_mean_distance / bg_mean_distance), yet it separates winners in *opposite
directions* -- lower is better here, higher is better in Prevotella. A metric
that flips direction between experiments is not a predictor, and any analysis
that supplies the expected direction per dataset is circular.
"""

import json
import statistics
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "clarke_2017_mtb.json"
PREV = Path(__file__).parent / "data" / "dwivedi_yu_2023_prevotella.json"


@pytest.fixture(scope="module")
def sets():
    return json.loads(DATA.read_text())["sets"]


@pytest.fixture(scope="module")
def winners(sets):
    return {n for n, v in sets.items() if v["outcome"]["most_effective"]}


def _separates(sets, winners, metric, lower_is_better):
    w = [sets[n]["published"][metric] for n in winners]
    o = [sets[n]["published"][metric] for n in sets if n not in winners]
    return (max(w) < min(o)) if lower_is_better else (min(w) > max(o))


# ----------------------------------------------------------------------
# Fixture integrity
# ----------------------------------------------------------------------


def test_fixture_is_intact(sets):
    assert len(sets) == 12
    for name, s in sets.items():
        assert len(s["primers"]) == s["published"]["size"], name
        assert all(set(p) <= set("ACGT") for p in s["primers"]), name


def test_the_four_footnoted_sets_are_the_winners(sets, winners):
    """Table 2 footnote: 'The sets that most effectively amplified Mycobacterium'."""
    assert winners == {"Mtb4", "Mtb6", "Mtb8", "Mtb9"}


def test_negative_controls_are_labelled(sets):
    """MtbUneven and MtbSparse were chosen to be bad, and should look bad."""
    controls = {n for n, v in sets.items() if v["outcome"]["negative_control"]}
    assert controls == {"MtbUneven", "MtbSparse"}
    for n in controls:
        assert not sets[n]["outcome"]["most_effective"]

    # MtbUneven was picked for the highest Gini, MtbSparse for the sparsest binding.
    assert sets["MtbUneven"]["published"]["fg_gini"] == max(
        s["published"]["fg_gini"] for s in sets.values()
    )
    assert sets["MtbSparse"]["published"]["fg_mean_distance"] == max(
        s["published"]["fg_mean_distance"] for s in sets.values()
    )


# ----------------------------------------------------------------------
# What predicts here
# ----------------------------------------------------------------------


def test_binding_density_separates_the_winners(sets, winners):
    """The four best sets are the four densest, with no overlap."""
    assert _separates(sets, winners, "fg_mean_distance", lower_is_better=True)


def test_evenness_does_not_separate_the_winners(sets, winners):
    """Gini runs the wrong way here, because the pool was pre-filtered for it.

    This is the mirror image of the Prevotella result and the reason neither
    dataset alone supports a general rule.
    """
    assert not _separates(sets, winners, "fg_gini", lower_is_better=True)

    w = statistics.median(sets[n]["published"]["fg_gini"] for n in winners)
    o = statistics.median(sets[n]["published"]["fg_gini"] for n in sets if n not in winners)
    assert w > o, (
        "the winners are expected to have slightly HIGHER (worse) Gini here; if "
        "that has changed, re-examine the combined synthesis"
    )


def test_worst_gap_does_not_separate_the_winners(sets, winners):
    """max_gap was the strongest predictor on Prevotella; it is not here."""
    assert not _separates(sets, winners, "fg_max_distance", lower_is_better=True)


# ----------------------------------------------------------------------
# The two benchmarks disagree, and that is the finding
# ----------------------------------------------------------------------


def test_no_single_metric_predicts_across_both_datasets(sets, winners):
    """The headline result of having two benchmarks rather than one.

    Density separates on Clarke and not on Dwivedi-Yu; the worst gap and Gini
    separate on Dwivedi-Yu and not on Clarke. Any guidance that names one metric
    as universally dominant is contradicted by one of the two.
    """
    prev = json.loads(PREV.read_text())["sets"]
    prev_winners = {"Prev03", "Prev06"}

    density_here = _separates(sets, winners, "fg_mean_distance", True)
    density_there = _separates(prev, prev_winners, "fg_mean_distance", True)
    holes_here = _separates(sets, winners, "fg_max_distance", True)
    holes_there = _separates(prev, prev_winners, "fg_max_distance", True)

    assert density_here and not density_there
    assert holes_there and not holes_here


def test_selectivity_ratio_flips_direction_between_datasets(sets, winners):
    """fg_bg_ratio is defined identically in both papers yet points both ways.

    Lower is better here; higher is better on Prevotella. Reporting that it
    "separates" in both is only possible by supplying the expected direction per
    dataset, which is circular. Recorded so that trap is not walked into again.
    """
    prev = json.loads(PREV.read_text())["sets"]
    prev_winners = {"Prev03", "Prev06"}

    def direction(s, w):
        wm = statistics.median(s[n]["published"]["fg_bg_ratio"] for n in w)
        om = statistics.median(s[n]["published"]["fg_bg_ratio"] for n in s if n not in w)
        return "lower" if wm < om else "higher"

    assert direction(sets, winners) == "lower"
    assert direction(prev, prev_winners) == "higher"


def test_evenness_was_constrained_by_design_here(sets):
    """The mechanism behind the disagreement, made checkable.

    Clarke pre-filtered the primer pool for even binding, so Gini varies far
    less across their sets than across Dwivedi-Yu's. With little variance left,
    it cannot explain the outcome - which is why density dominates here.
    """
    prev = json.loads(PREV.read_text())["sets"]
    real = [
        v["published"]["fg_gini"] for n, v in sets.items() if not v["outcome"]["negative_control"]
    ]
    theirs = [v["published"]["fg_gini"] for v in prev.values()]

    assert max(real) - min(real) < max(theirs) - min(theirs), (
        "Clarke's Gini spread should be narrower than Dwivedi-Yu's; if it is not, "
        "the explanation for why the two benchmarks disagree needs revisiting"
    )
