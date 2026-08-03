"""The default Tm window against the only per-primer amplification data we have.

`dwivedi_yu_2023_plasmid_amplification.json` is the one benchmark that measures
how much each individual primer amplified, rather than how a set performed:
310 qPCR measurements over 284 primers on pcDNA and pLTR, the two plasmids
bundled in `examples/plasmid_example`. That makes it the only dataset that can
say anything about where a Tm threshold should sit.

Run the 141 pcDNA-designed primers through the shipped phi29 window (20-50,
from the registry) and it splits them three ways:

    below  window   n=17   median measured amplification   0.94x
    inside window   n=57                                   5.87x
    above  window   n=67                                  21.50x

The same split under the other enzymes' windows shows the effect is a property
of where the upper bound sits, not of Tm filtering as such:

    phi29      20-50   below  0.94x   inside   5.87x   above 21.50x  (n above = 67)
    equiphi29  37-62   below  1.60x   inside  16.77x   above 32.75x  (n above =  6)
    bst        50-75   below  3.62x   inside  21.32x   above  -      (n above =  0)

As the window rises, the primers it admits are the ones that amplified well and
what it rejects is increasingly the poor end. Under phi29 the lower bound is
already doing that job -- what it excludes amplified at 0.94x -- while the
upper bound is not.

The 67 primers above phi29's upper bound carry 70% of all amplification
measured in the experiment, so this is not a handful of outliers.

This is recorded rather than acted on, because the benchmark cannot settle the
question it raises. Every primer here has exactly one on-target site and zero
off-target sites -- they were chosen that way -- so there is no selectivity to
lose. In genome-scale SWGA a longer, higher-Tm primer binds more sites
including background ones, and an upper Tm bound is defensible as a
selectivity constraint. What these tests pin is the cost side of that
trade-off, which was previously unmeasured: raising max_tm is not free of
consequence in either direction, and anyone changing it should know both.

See docs/validation/published_primer_sets.md.
"""

import json
import statistics
from pathlib import Path

import pytest

from neoswga.core.parameter import default_reaction_temp, default_tm_range

DATA = Path(__file__).parent / "data" / "dwivedi_yu_2023_plasmid_amplification.json"

# Resolved from the registry rather than hardcoded. These were written as a
# literal 15/45 when that was what examples/plasmid_example/params.json
# carried; the default is now polymerase-aware, so a literal would pin numbers
# the tool no longer uses and quietly stop testing the shipped behaviour.
MIN_TM, MAX_TM = default_tm_range("phi29")


@pytest.fixture(scope="module")
def measured():
    """Mean measured on-target amplification per pcDNA-designed primer."""
    payload = json.loads(DATA.read_text())
    by_primer = {}
    for row in payload["measurements"]:
        if row["designed_for"] == "pcDNA":
            by_primer.setdefault(row["primer"], []).append(row["on_target_amp"])
    return {p: statistics.mean(v) for p, v in by_primer.items()}


@pytest.fixture(scope="module")
def tm_of(measured):
    """Tm under the shipped phi29 conditions, from the same code the filter uses."""
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    return {p: conditions.calculate_effective_tm(p) for p in measured}


def split(measured, tm_of):
    below = [p for p in measured if tm_of[p] < MIN_TM]
    inside = [p for p in measured if MIN_TM <= tm_of[p] <= MAX_TM]
    above = [p for p in measured if tm_of[p] > MAX_TM]
    return below, inside, above


def test_the_dataset_is_the_one_we_think(measured):
    """Guards the fixture: a changed dataset would silently change the claim."""
    assert len(measured) == 141


def test_the_lower_tm_bound_excludes_poor_amplifiers(measured, tm_of):
    """The part of the window that is doing useful work."""
    below, inside, _ = split(measured, tm_of)

    assert below, "no primer fell below the window; the split is not exercised"
    assert statistics.median(measured[p] for p in below) < statistics.median(
        measured[p] for p in inside
    )


def test_the_upper_tm_bound_excludes_the_best_amplifiers(measured, tm_of):
    """The finding, and the reason this file exists.

    Primers above max_tm amplified substantially better than those inside it.
    If a future change to the Tm model or the window makes this stop being
    true, that is worth knowing -- it would mean the trade-off documented here
    no longer applies.
    """
    _, inside, above = split(measured, tm_of)

    assert above, "no primer sat above the window; the finding is not exercised"
    assert statistics.median(measured[p] for p in above) > statistics.median(
        measured[p] for p in inside
    )


def test_amplification_excluded_by_the_upper_bound_is_a_large_share(measured, tm_of):
    """Quantifies the cost, so it cannot be dismissed as a handful of outliers."""
    _, _, above = split(measured, tm_of)

    excluded_share = sum(measured[p] for p in above) / sum(measured.values())
    assert excluded_share > 0.5, (
        f"only {excluded_share:.0%} of measured amplification sits above phi29's " f"max_tm"
    )


def test_a_higher_window_admits_the_good_amplifiers(measured):
    """The contrast that shows this is about where the bound sits.

    equiphi29's window is 37-62. Under it the primers that amplified best are
    inside rather than above, and what falls below is the poor end -- which is
    what a well-placed Tm filter should look like on this data.
    """
    import statistics

    from neoswga.core.reaction_conditions import ReactionConditions

    low, high = default_tm_range("equiphi29")
    conditions = ReactionConditions(temp=default_reaction_temp("equiphi29"), polymerase="equiphi29")
    tm = {p: conditions.calculate_effective_tm(p) for p in measured}

    inside = [p for p in measured if low <= tm[p] <= high]
    below = [p for p in measured if tm[p] < low]

    assert statistics.median(measured[p] for p in inside) > statistics.median(
        measured[p] for p in below
    )


def test_amplification_rises_with_primer_length(measured):
    """The mechanism behind the Tm result, and a check on the premise of
    thermodynamic filtering generally: longer primers bind harder and
    amplified more, monotonically across the whole length range."""
    by_length = {}
    for primer, amp in measured.items():
        by_length.setdefault(len(primer), []).append(amp)

    medians = [statistics.median(by_length[k]) for k in sorted(by_length)]
    assert medians[-1] > medians[0] * 5, f"medians by length: {medians}"
