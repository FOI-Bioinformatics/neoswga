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
