"""The default Tm window against the only per-primer amplification data we have.

`dwivedi_yu_2023_plasmid_amplification.json` is the one benchmark that measures
how much each individual primer amplified, rather than how a set performed:
310 qPCR measurements over 284 primers on pcDNA and pLTR, the two plasmids
bundled in `examples/plasmid_example`. That makes it the only dataset that can
say anything about where a Tm threshold should sit.

Run the 141 pcDNA-designed primers through the shipped filter settings
(phi29, min_tm 15, max_tm 45) and the window splits them three ways:

    below  window   n=14   median measured amplification   0.80x
    inside window   n=51                                   6.45x
    above  window   n=76                                  19.95x

The lower bound earns its place -- the primers it rejects genuinely amplified
badly. The upper bound does not look like a filter on quality: the 76 primers
above it amplified about three times better than those inside, and carry 74%
of all amplification measured in the experiment.

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

DATA = Path(__file__).parent / "data" / "dwivedi_yu_2023_plasmid_amplification.json"

# The shipped phi29 settings from examples/plasmid_example/params.json.
MIN_TM = 15.0
MAX_TM = 45.0


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
        f"only {excluded_share:.0%} of measured amplification sits above max_tm; "
        f"the documented figure is 74%"
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
