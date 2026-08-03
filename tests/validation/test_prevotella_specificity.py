"""The first validation against a measured specificity outcome, and what it says.

`dwivedi_yu_2023_prevotella.json` is the only fixture carrying a continuous
measured outcome: six sets with percent-of-reads-on-target spanning 1.2% to
72.2%, against human background. It is the natural test of whether predicted
selectivity means anything.

It needs genomes, which are not committed. `scripts/fetch_reference_genomes.py`
fetches Prevotella (3 Mb) and human chr21 (46 Mb, a subset standing in for the
3.1 Gb background); these tests skip without them.

The result is negative, and the reason matters more than the result.

    occupancy selectivity vs measured % on target   rho = -0.086  p = 0.872

The model does not rank these sets. Neither does anything else built on
background binding, including the metric the authors published:

    fg_bg_ratio        rho = +0.600   p = 0.208
    bg_mean_distance   rho = -0.086   p = 0.872
    fg_max_distance    rho = -0.829   p = 0.042     <- the one that predicts

What separates these six experiments is foreground gap structure, not
specificity. And there is a structural reason no specificity metric could have
worked here: five of the six sets are made entirely of 8-mers, and every 8-mer
occurs about 97 times in Prevotella and 1219 times in chr21. All 32,896
canonical 8-mers are present in both genomes. At that length a primer's
background load is set by how much sequence it is counted against, not by its
sequence -- so every set looks alike, and the ratios cluster in a band too
narrow to rank.

The consequence for design is sharper than the failed correlation. **Additives
cannot buy specificity for a set of 7-8mers**, because those primers bind
everywhere regardless of stringency; there is no mismatched population to
discriminate against when the perfect matches already number in the thousands.
The lever needs primers long enough for exact matches to be rare -- at k=11 the
mean count falls to 2.5 and 21.3 -- which is the equiphi29 regime, not phi29's
6-12mers.

So this benchmark does not validate the occupancy model and cannot: its outcome
variance lives in a foreground property, and its primers are too short for
sequence-specific background discrimination to exist. A dataset that could
validate it would need longer primers and measured on-target yield. Recording
that is more useful than recording a correlation of -0.086.
"""

import json
import os
from pathlib import Path

import pytest

pytest.importorskip("scipy")

from scipy.stats import spearmanr

DATA = Path(__file__).parent / "data" / "dwivedi_yu_2023_prevotella.json"
GENOMES = Path(__file__).parent / "genomes"

FG_PREFIX = GENOMES / "prevotella"
BG_PREFIX = GENOMES / "human_chr21"

pytestmark = pytest.mark.skipif(
    not (GENOMES / "prevotella_8mer_all.txt").exists(),
    reason="reference genomes absent; run scripts/fetch_reference_genomes.py "
    "prevotella human_chr21 and neoswga count-kmers",
)


@pytest.fixture(scope="module")
def sets():
    return sorted(json.loads(DATA.read_text())["sets"].items())


@pytest.fixture(scope="module")
def predicted(sets):
    """Occupancy-weighted selectivity per set, at the published conditions."""
    from neoswga.core.base_optimizer import _selectivity_from_loads
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    values = []
    for _name, entry in sets:
        fg = weighted_site_load(entry["primers"], [str(FG_PREFIX)], conditions, 1)
        bg = weighted_site_load(entry["primers"], [str(BG_PREFIX)], conditions, 1)
        values.append(_selectivity_from_loads(fg, bg))
    return values


# ----------------------------------------------------------------------
# Why no specificity metric can work on this benchmark
# ----------------------------------------------------------------------


def test_every_short_kmer_is_present_in_both_genomes():
    """The structural reason, and the finding worth keeping.

    All 32,896 canonical 8-mers occur in both a 3 Mb bacterium and a 46 Mb
    human chromosome. A primer that short has no sequence-specific background
    load -- only a size-determined one -- so a metric built on binding counts
    cannot separate two sets of 8-mers however it weights them.
    """
    from neoswga.core.mismatch_counts import load_kmer_counts

    for k, expected_canonical in ((7, 8_192), (8, 32_896)):
        background = load_kmer_counts(str(BG_PREFIX), k)
        assert len(background) == expected_canonical, (
            f"expected every canonical {k}-mer present in chr21; "
            f"found {len(background)} of {expected_canonical}"
        )


def test_short_primers_bind_the_background_thousands_of_times():
    """Which is why stringency has nothing to work with: there is no rare
    perfect-match population to protect against a common mismatched one."""
    from neoswga.core.mismatch_counts import load_kmer_counts

    counts = load_kmer_counts(str(BG_PREFIX), 8)
    mean_occurrences = sum(counts.values()) / len(counts)

    assert mean_occurrences > 500, (
        f"8-mers occur {mean_occurrences:.0f} times on average in chr21; the "
        f"claim that short primers saturate a large background rests on this"
    )


def test_longer_primers_are_rare_enough_for_sequence_to_matter():
    """Where the lever does apply. At k=11 a primer occurs a couple of dozen
    times rather than a couple of thousand, so which sequence it is starts to
    determine its background load."""
    from neoswga.core.mismatch_counts import load_kmer_counts

    counts = load_kmer_counts(str(BG_PREFIX), 11)
    mean_occurrences = sum(counts.values()) / len(counts)

    assert mean_occurrences < 50


# ----------------------------------------------------------------------
# The negative result itself, recorded rather than hidden
# ----------------------------------------------------------------------


def test_predicted_selectivity_does_not_rank_these_sets(sets, predicted):
    """Stated as a fact about this benchmark, not asserted as desirable.

    If a future change makes the model rank them, this test fails and the
    documented conclusion needs revisiting -- which is the point of pinning a
    negative result rather than omitting it.
    """
    measured = [entry["outcome"]["pct_reads_on_target"] for _name, entry in sets]
    rho = spearmanr(predicted, measured)

    assert abs(rho.statistic) < 0.5, (
        f"predicted selectivity now correlates with outcome (rho={rho.statistic:+.3f}); "
        f"the documented conclusion that this benchmark cannot validate a "
        f"specificity model needs revisiting"
    )


def test_the_published_background_metric_does_no_better(sets):
    """The check that stops this being read as a defect in our model.

    The authors' own foreground/background ratio does not significantly predict
    their own outcomes either. Whatever separated these experiments, it was not
    background binding.
    """
    measured = [entry["outcome"]["pct_reads_on_target"] for _name, entry in sets]
    published = [entry["published"]["fg_bg_ratio"] for _name, entry in sets]

    rho = spearmanr(published, measured)
    assert rho.pvalue > 0.05


def test_a_foreground_gap_metric_is_what_predicts_here(sets):
    """Consistent with docs/validation/published_primer_sets.md: whichever
    property is currently limiting is the one that predicts. In these six
    experiments that was coverage evenness, not specificity."""
    measured = [entry["outcome"]["pct_reads_on_target"] for _name, entry in sets]
    max_gap = [entry["published"]["fg_max_distance"] for _name, entry in sets]

    rho = spearmanr(max_gap, measured)
    assert rho.statistic < 0, "larger maximum gap should mean worse outcome"
    assert rho.pvalue < 0.05


def test_the_predicted_values_are_too_close_together_to_rank(predicted):
    """The mechanism behind the null result, made checkable.

    The six predictions span a narrow band because the sets are made of
    primers too short to differ in background load. A metric whose spread is
    this small is not ranking anything, whatever its correlation happens to be.
    """
    spread = max(predicted) / min(predicted)

    assert spread < 2.0, (
        f"predicted selectivity spans {spread:.1f}x across sets whose measured "
        f"outcomes differ 60-fold"
    )


# ----------------------------------------------------------------------
# Where the lever does work
# ----------------------------------------------------------------------

needs_12mers = pytest.mark.skipif(
    not (GENOMES / "prevotella_12mer_all.txt").exists(),
    reason="12-mer counts absent; run neoswga count-kmers with min_k=max_k=12",
)

needs_length_series = pytest.mark.skipif(
    not all((GENOMES / f"prevotella_{k}mer_all.txt").exists() for k in (10, 11, 12)),
    reason="10/11/12-mer counts absent; run neoswga count-kmers at each length",
)

# Sets the pipeline itself produced from Prevotella against human chr21,
# equiphi29 at 42 C, eight primers each, candidate pool 2000. Fixed here rather
# than regenerated so the tests do not depend on a 2-minute filter run per
# length. tests/validation/genomes/equiphi29_{k}mer.json reproduces them.
PIPELINE_SETS = {
    10: [
        "AAGGGCGTTA",
        "ACTCGTCAAC",
        "AAGACGCCTA",
        "ATGCTGACGA",
        "CGAGTGCAAA",
        "AGCGTATTGC",
        "CAACTGACGA",
        "ACCGATACCA",
    ],
    11: [
        "ACGTAAGAACC",
        "GGTACGATGGA",
        "ATTAGGTGCGA",
        "AACATCTTCGA",
        "AAGCGTATGAC",
        "ATCCTCGTTGA",
        "GATATGACGGA",
        "CACGAATCTGA",
    ],
    12: [
        "CATCAACGATAA",
        "CCCAGAACGATA",
        "GCGTATCAACCA",
        "TATGCGCTGCAA",
        "GTCAAAGTCGAA",
        "ACGAGTTACCAT",
        "CAACACGCTTAC",
        "CGTCAGTTTATC",
    ],
}

PIPELINE_12MERS = PIPELINE_SETS[12]


def _selectivity_curve(primers, temperatures):
    from neoswga.core.base_optimizer import _selectivity_from_loads
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    values = []
    for temp in temperatures:
        conditions = ReactionConditions(temp=temp, polymerase="equiphi29")
        fg = weighted_site_load(primers, [str(FG_PREFIX)], conditions, 1)
        bg = weighted_site_load(primers, [str(BG_PREFIX)], conditions, 1)
        values.append(_selectivity_from_loads(fg, bg))
    return values


def test_short_primers_do_not_respond_to_stringency_at_all(sets):
    """The mechanism behind the null result above, shown rather than argued.

    For a set of 7-8mers, foreground and background carry the same mixture of
    exact and mismatched sites -- both saturated -- so the occupancy terms
    cancel in the ratio and no temperature changes it. Measured across 38-50 C
    the published sets are flat to two decimal places:

        Prev06 (7,8-mers)   0.41  0.41  0.41  0.41
        Prev03 (8-mers)     0.30  0.30  0.30  0.30
    """
    for name, entry in sets:
        if max(len(p) for p in entry["primers"]) > 9:
            continue
        curve = _selectivity_curve(entry["primers"], [38.0, 42.0, 46.0, 50.0])
        assert (
            curve[-1] / curve[0] < 1.05
        ), f"{name} responded to stringency: {[round(v, 3) for v in curve]}"


@needs_length_series
@pytest.mark.parametrize("k", sorted(PIPELINE_SETS))
def test_longer_primers_do_respond_to_stringency(k):
    """The design consequence, on the same target and background.

    At 10 bases and above the foreground is dominated by exact matches and the
    background by mismatched ones, so the two no longer share a composition and
    stringency separates them. Measured across 38-50 C:

        k=10   0.53  0.61  0.68  0.71    1.36x
        k=11   0.51  0.59  0.67  0.70    1.36x
        k=12   0.59  0.66  0.78  0.88    1.51x

    against 1.00x for the published 7-8mer sets. This is what makes the
    additive lever a real control in the equiphi29 regime and not in phi29's
    6-12mers.
    """
    curve = _selectivity_curve(PIPELINE_SETS[k], [38.0, 42.0, 46.0, 50.0])

    assert curve == sorted(curve), f"k={k} not monotone in temperature: {curve}"
    assert curve[-1] / curve[0] > 1.25, f"k={k} got {[round(v, 3) for v in curve]}"


@needs_length_series
def test_the_two_regimes_are_separated_by_a_clear_margin():
    """The trend itself, stated as the thing that is actually robust.

    The fold response is set-dependent -- a six-primer 12-mer set measured
    2.05x where the eight-primer one measures 1.51x -- so no particular fold
    value is worth pinning. What is stable is the gap: every set of 10-12mers
    responds to stringency and every set of 7-8mers is flat to two decimal
    places, with nothing in between. A change that narrowed this gap would mean
    the length threshold had moved, which is the design-relevant fact.
    """
    published = json.loads(DATA.read_text())["sets"]
    short_folds = []
    for entry in published.values():
        if max(len(p) for p in entry["primers"]) > 9:
            continue
        curve = _selectivity_curve(entry["primers"], [38.0, 50.0])
        short_folds.append(curve[-1] / curve[0])

    long_folds = [
        _selectivity_curve(primers, [38.0, 50.0])[-1] / _selectivity_curve(primers, [38.0, 50.0])[0]
        for primers in PIPELINE_SETS.values()
    ]

    assert short_folds, "no short-primer sets in the fixture"
    assert min(long_folds) > max(short_folds) + 0.2, (
        f"the regimes have merged: 7-8mers fold {[round(f, 3) for f in short_folds]}, "
        f"10-12mers fold {[round(f, 3) for f in long_folds]}"
    )


@needs_12mers
def test_exact_counting_calls_this_set_perfectly_specific():
    """Why the occupancy model is worth its cost, on a real design.

    The pipeline's 12-mer set has zero exact matches in the background, so a
    count-based selectivity reports no background binding whatsoever. Weighted
    by occupancy the same set measures below 1 -- it binds the human background
    more than the target, entirely through near-matches a count cannot see.
    """
    from neoswga.core.mismatch_counts import mismatch_class_counts
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    exact_bg = sum(mismatch_class_counts(p, [str(BG_PREFIX)], 0)[0] for p in PIPELINE_12MERS)
    assert exact_bg == 0, "the premise is that exact counting sees no background here"

    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")
    weighted_bg = weighted_site_load(PIPELINE_12MERS, [str(BG_PREFIX)], conditions, 1)
    assert weighted_bg > 10, f"expected substantial near-match load; got {weighted_bg:.2f}"


@needs_length_series
def test_the_two_models_disagree_in_direction_where_both_have_data():
    """The sharper version of the case above, with no zero to hide behind.

    At k=11 and 12 the background has no exact matches at all, so the count
    model reports perfect specificity from an empty denominator -- visibly
    uninformative. At k=10 it has data: 442 foreground and 51 background sites,
    a count selectivity near 9, which reads as a good design. Occupancy
    weighting puts the same set below 1, on the other side of the decision
    boundary, because 51 exact background sites are accompanied by a much
    larger near-match population that a count does not see.

    Both models are confident and they disagree about whether the set is
    selective at all. That is the case for carrying the occupancy term.
    """
    from neoswga.core.base_optimizer import _selectivity_from_loads
    from neoswga.core.mismatch_counts import mismatch_class_counts
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    primers = PIPELINE_SETS[10]
    exact_fg = sum(mismatch_class_counts(p, [str(FG_PREFIX)], 0)[0] for p in primers)
    exact_bg = sum(mismatch_class_counts(p, [str(BG_PREFIX)], 0)[0] for p in primers)
    assert exact_bg > 0, "the premise is that the count model has data here"

    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")
    occupancy = _selectivity_from_loads(
        weighted_site_load(primers, [str(FG_PREFIX)], conditions, 1),
        weighted_site_load(primers, [str(BG_PREFIX)], conditions, 1),
    )

    assert (
        exact_fg / exact_bg > 2.0
    ), f"count model should look favourable; got {exact_fg}/{exact_bg}"
    assert occupancy < 1.0, f"occupancy should disagree; got {occupancy:.3f}"
