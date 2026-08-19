"""The per-primer coverage reach, fitted to a measured outcome.

Coverage is meaningless without the reach it is computed at -- the same set
scores 0.418 at 3 kb, 0.836 at 10 kb and ~1.0 at 70 kb -- and neither of the two
reaches this tool uses had ever been measured. `calibrate-reach --bam` exists
for it and has never been run; there is no BAM in the repository.

Clarke et al. (2017) is the substitute. Their `TmL/Even` Wolbachia set is
published with its sequences and reached 10x depth over 60-75% of the genome
against 2.8% unamplified, so the question "what reach reproduces that?" has an
answer. It is 3.0-4.5 kb on the sites found here, and 4.1-6.2 kb under the site
count the paper's own mean-gap statistic implies.

Full write-up, including why three of five sets validate the site-finding to
within one site and two do not: docs/validation/reach_calibration.md

Skips without the genome; fetch with
`python scripts/fetch_reference_genomes.py wolbachia`.
"""

import json
import re
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).parent / "data" / "clarke_2017_wolbachia.json"
GENOME = Path(__file__).parent / "genomes" / "wolbachia.fna"

pytestmark = pytest.mark.skipif(
    not GENOME.exists(),
    reason="run scripts/fetch_reference_genomes.py wolbachia",
)

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(seq):
    return seq.translate(_COMP)[::-1]


@pytest.fixture(scope="module")
def genome():
    return "".join(
        line.strip().upper() for line in GENOME.read_text().splitlines() if not line.startswith(">")
    )


@pytest.fixture(scope="module")
def sets():
    return json.loads(DATA.read_text())["sets"]


def _sites(genome, primers):
    """Every position a primer can anneal: the primer or its reverse complement
    on the plus strand. Circular, so a match spanning the origin is found."""
    found = set()
    for primer in primers:
        padded = genome + genome[: len(primer) - 1]
        for pattern in {primer, _rc(primer)}:
            found.update(
                m.start() for m in re.finditer(f"(?={pattern})", padded) if m.start() < len(genome)
            )
    return sorted(found)


def _coverage(sites, genome_len, reach):
    occupied = np.zeros(genome_len, dtype=bool)
    for pos in sites:
        occupied[max(0, pos - reach) : min(genome_len, pos + reach)] = True
    return occupied.sum() / genome_len


# ----------------------------------------------------------------------
# The site-finding is checked against the paper before anything is fitted
# ----------------------------------------------------------------------


@pytest.mark.parametrize("name", ["TmH/Even", "TmH/Selective", "Leichty and Brisson 2014"])
def test_site_counts_match_the_published_statistics(genome, sets, name):
    """Three of the five sets reproduce the paper's own mean binding distance to
    within one site, which is what makes the fit below credible."""
    entry = sets[name]
    found = len(_sites(genome, entry["primers"]))
    implied = round(len(genome) / entry["published"]["fg_mean_distance"])

    assert abs(found - implied) <= 1, f"{found} found against {implied} implied"


@pytest.mark.parametrize("name", ["TmL/Even", "TmL/Selective"])
def test_the_two_multi_primer_sets_do_not_match_and_that_is_recorded(genome, sets, name):
    """The unexplained half, pinned so it is not mistaken for agreement.

    Both are the multi-primer TmL sets and both carry roughly 1.8x the sites
    their published mean gap implies. No convention tried reconciles them, there
    are no reverse-complement pairs within either, and every primer is found
    with a plausible count. It is why the calibrated range is 3.0-6.2 kb rather
    than 3.0-4.5.
    """
    entry = sets[name]
    found = len(_sites(genome, entry["primers"]))
    implied = round(len(genome) / entry["published"]["fg_mean_distance"])

    assert 1.6 <= found / implied <= 2.0, f"the discrepancy has changed: {found} vs {implied}"


# ----------------------------------------------------------------------
# The fit
# ----------------------------------------------------------------------


def test_the_effective_set_reaches_its_measured_breadth_at_the_default_reach(genome, sets):
    """`TmL/Even` achieved 10x over 60-75%. At the shipped 3 kb selection reach
    the model puts it at the bottom of that band, and it leaves the band by
    4.5 kb -- so the default is corroborated rather than merely assumed."""
    sites = _sites(genome, sets["TmL/Even"]["primers"])

    assert _coverage(sites, len(genome), 3000) == pytest.approx(0.62, abs=0.03)
    assert 0.60 <= _coverage(sites, len(genome), 3000) <= 0.75
    assert _coverage(sites, len(genome), 4500) >= 0.72


def test_the_reporting_reach_is_outside_the_calibrated_band(genome, sets):
    """`headline_coverage` reports at a 10 kb measured product length. That is a
    defensible quantity in its own right, but this outcome does not support it:
    the same set reads 0.96 there against 0.62 at the calibrated reach."""
    sites = _sites(genome, sets["TmL/Even"]["primers"])

    assert _coverage(sites, len(genome), 10000) > 0.90


@pytest.mark.parametrize("name", ["TmH/Even", "TmH/Selective"])
def test_the_ineffective_sets_stay_well_below_at_the_calibrated_reach(genome, sets, name):
    """The harder test: an absolute fit could be coincidence, an ordering is
    less likely to be. Both ineffective sets need 6.8-8.3 kb just to reach 60%.
    """
    ineffective = _coverage(_sites(genome, sets[name]["primers"]), len(genome), 3000)
    effective = _coverage(_sites(genome, sets["TmL/Even"]["primers"]), len(genome), 3000)

    assert ineffective < 0.45
    assert effective - ineffective > 0.15


def test_coverage_does_not_separate_the_two_low_tm_sets(genome, sets):
    """And the limit of what this metric can do, recorded so it is not oversold.

    `TmL/Even` and `TmL/Selective` were built as a controlled contrast between
    evenness and selectivity, and evenness separated them in the lab. Coverage
    within reach is not an evenness metric and puts them 0.02 apart; `max_gap`
    and `gap_gini` are what carry that signal.
    """
    even = _coverage(_sites(genome, sets["TmL/Even"]["primers"]), len(genome), 3000)
    selective = _coverage(_sites(genome, sets["TmL/Selective"]["primers"]), len(genome), 3000)

    assert abs(even - selective) < 0.05
    assert sets["TmL/Even"]["published"]["fg_gini"] < sets["TmL/Selective"]["published"]["fg_gini"]
