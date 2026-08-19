"""One set of grading thresholds, calibrated against the published benchmarks.

`report/quality.py` and `results_interpreter.py` each carried their own table
and disagreed. Enrichment "excellent" was 500x in one and 200x in the other, so
`neoswga report` and `neoswga interpret` could hand a user different verdicts on
the same result. Uniformity disagreed more sharply still: one graded 1 - Gini
against a 0.85 bar, the other Gini against a 0.30 bar.

Worse than the disagreement, both were miscalibrated. Grading every set in the
validation suite -- the eleven with published Gini and a measured outcome:

    set                   gini   outcome            old report   old interpret
    Prev03               0.533   96x, success       poor         acceptable
    TmL/Even             0.537   effective          poor         acceptable
    Prev06               0.552   120x, success      poor         acceptable
    Prev02               0.566   2.1x, failure      poor         acceptable
    Prev04               0.623   13x, success       critical     poor
    Leichty 2014         0.712   partial            critical     poor

No design that has ever been published reached "good". `Prev06` is the best set
in the entire suite -- 120-fold enrichment, 72% of reads on target -- and the
report called its uniformity poor. A scale on which nothing achievable scores
well is not a scale; it is a constant.

The thresholds below are read off that distribution instead. Observed Gini
spans 0.533 to 0.712 across every published set, and `BENCHMARK_GINI_GOOD`
(0.56, in `coverage.py`) is the cut already drawn between the two Prevotella
winners and the four that failed.

**What this grade is not.** The Gini grade is descriptive, not predictive: it
does not separate winners from losers on the available data. `TmL/Even` (effective) and
`TmH/Even` (ineffective) have identical Gini of 0.537 and opposite outcomes;
`Prev04` succeeded at 0.623 while `Prev02` failed at 0.566; and on Clarke's
M. tuberculosis sets evenness ran backwards outright. This grade says where a
design sits in the published distribution, nothing more. The position the
project takes elsewhere -- surface `mean_gap`, `max_gap` and `gap_gini` with
their benchmark context rather than compressing them into one verdict -- is the
one to trust; see docs/validation/published_primer_sets.md.
"""

from typing import Dict

# ---------------------------------------------------------------------------
# Coverage: fraction of the target within reach of a selected binding site.
#
# Read at the reach the figure was computed at -- 95% at a 3 kb per-primer
# reach and 95% at 10 kb are different designs. See
# docs/validation/reach_calibration.md.
# ---------------------------------------------------------------------------
COVERAGE = {
    "excellent": 0.95,
    "good": 0.85,
    "acceptable": 0.70,
    "poor": 0.50,
}

# ---------------------------------------------------------------------------
# Binding-site evenness, as Gini (LOWER is better).
#
# Calibrated to the published distribution: 0.53 is the best observed anywhere,
# 0.56 is the cut between the Prevotella winners and the sets that failed, and
# 0.71 is the worst published set. Descriptive, not predictive -- see the module
# docstring.
# ---------------------------------------------------------------------------
GINI = {
    "excellent": 0.53,
    "good": 0.56,
    "acceptable": 0.62,
    "poor": 0.70,
}

# ---------------------------------------------------------------------------
# Enrichment.
#
# CAVEAT, and it is a real one: the two consumers do not grade the same
# quantity. `report/metrics.py` computes a site-DENSITY ratio (target sites per
# Mb over background sites per Mb), while `results_interpreter` estimates a
# FOLD enrichment from amplification factors. They were nonetheless graded with
# near-identical tables, which is how a 500-versus-200 disagreement at
# "excellent" went unnoticed.
#
# The bar is set from measured outcomes rather than from either old table. The
# best fold enrichment in the suite is 120x (Prev06); 96x also succeeded, 13x
# worked modestly, and everything at 2x failed. A bar at 500x, or at 200x, sat
# above every result ever published.
#
# Reconciling the two quantities is open work; until then this is one scale
# applied to two things, and the report labels which it is showing.
# ---------------------------------------------------------------------------
ENRICHMENT = {
    "excellent": 100.0,
    "good": 50.0,
    "acceptable": 10.0,
    "poor": 2.0,
}

# ---------------------------------------------------------------------------
# Within-set Tm spread (degrees C, lower is better) and dimer risk (0-1).
# Neither is calibrated against an outcome; both are operating points.
# On Dwivedi-Yu's Leishmania sets the Tm-uniformity term ran backwards -- PS4
# won with the widest spread of the four -- so treat this one as a preference,
# not a finding.
# ---------------------------------------------------------------------------
TM_RANGE = {
    "excellent": 3.0,
    "good": 5.0,
    "acceptable": 8.0,
    "poor": 12.0,
}

DIMER_RISK = {
    "excellent": 0.1,
    "good": 0.2,
    "acceptable": 0.35,
    "poor": 0.5,
}

RATING_ORDER = ("excellent", "good", "acceptable", "poor")


def as_uniformity_score(thresholds: Dict[str, float]) -> Dict[str, float]:
    """Gini thresholds expressed as a `1 - Gini` score, higher-is-better.

    `report/quality.py` reports uniformity that way round. Deriving it here
    rather than restating it is what stops the two drifting apart again.
    """
    return {label: 1.0 - value for label, value in thresholds.items()}
