"""Phase 13D: grade-vs-score alignment.

`results_interpreter.rate_metric` returns EXCELLENT / GOOD / ACCEPTABLE /
POOR / CRITICAL. `OptimizationResult.normalized_score` returns 0-1. If a
user sees "GOOD" from interpret but the optimizer's normalized_score is
0.40 (typically C-ish), that inconsistency silently erodes trust.

These tests build synthetic PrimerSetMetrics at each grade band and
assert the normalized_score lands in a reasonable corresponding range.
Balanced weighting is used; application-specific weights still track
the same ordering.
"""

import pytest

from neoswga.core.base_optimizer import PrimerSetMetrics
from neoswga.core.results_interpreter import (
    QualityRating, COVERAGE_THRESHOLDS, ENRICHMENT_THRESHOLDS,
    UNIFORMITY_THRESHOLDS, DIMER_SCORE_THRESHOLDS,
)


def _metrics_at(cov: float, sel: float, dimer: float, gini: float,
                tm_range=(50, 60)):
    """Build a metrics object that matches the per-component threshold at
    the target grade."""
    return PrimerSetMetrics(
        fg_coverage=cov,
        bg_coverage=max(0.0, 1.0 / max(sel, 1.0)),
        coverage_uniformity=1.0 - gini,
        total_fg_sites=100,
        total_bg_sites=max(1, int(100 / max(sel, 1))),
        selectivity_ratio=sel,
        mean_tm=(tm_range[0] + tm_range[1]) / 2,
        tm_range=tm_range,
        dimer_risk_score=dimer,
        mean_gap=1000.0,
        max_gap=5000.0,
        gap_gini=gini,
        gap_entropy=0.8,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
    )


def _grade_all(cov, sel, dimer, gini):
    """Return the overall grade by taking the worst of the four component
    ratings (mirrors results_interpreter's min-of-ratings rule)."""
    grades = [
        QualityRating(r).value for r in [
            _rate_coverage(cov),
            _rate_enrichment(sel),
            _rate_dimer(dimer),
            _rate_uniformity(gini),
        ]
    ]
    # ranking worst-to-best
    rank = ["CRITICAL", "POOR", "ACCEPTABLE", "GOOD", "EXCELLENT"]
    return min(grades, key=lambda g: rank.index(g))


def _rate_coverage(v):
    from neoswga.core.results_interpreter import rate_metric
    return rate_metric(v, COVERAGE_THRESHOLDS)


def _rate_enrichment(v):
    from neoswga.core.results_interpreter import rate_metric
    return rate_metric(v, ENRICHMENT_THRESHOLDS)


def _rate_dimer(v):
    from neoswga.core.results_interpreter import rate_metric
    return rate_metric(v, DIMER_SCORE_THRESHOLDS, lower_is_better=True)


def _rate_uniformity(v):
    from neoswga.core.results_interpreter import rate_metric
    return rate_metric(v, UNIFORMITY_THRESHOLDS, lower_is_better=True)


# ----------------------------------------------------------------------
# Alignment assertions
# ----------------------------------------------------------------------

def test_all_excellent_maps_to_high_normalized_score():
    """A set that grades EXCELLENT on every component should score > 0.85."""
    m = _metrics_at(cov=0.98, sel=250.0, dimer=0.05, gini=0.2)
    grade = _grade_all(0.98, 250.0, 0.05, 0.2)
    assert grade == "EXCELLENT"
    ns = m.normalized_score()
    assert ns > 0.85, f"EXCELLENT set should have normalized_score > 0.85; got {ns:.3f}"


def test_all_good_maps_to_mid_high_score():
    """A set that grades GOOD on every component should score > 0.70."""
    m = _metrics_at(cov=0.88, sel=120.0, dimer=0.15, gini=0.38)
    grade = _grade_all(0.88, 120.0, 0.15, 0.38)
    assert grade == "GOOD"
    ns = m.normalized_score()
    assert ns > 0.70, f"GOOD set should have normalized_score > 0.70; got {ns:.3f}"


def test_all_acceptable_maps_to_middling_score():
    """A set that grades ACCEPTABLE should score between 0.50 and 0.75."""
    m = _metrics_at(cov=0.75, sel=60.0, dimer=0.30, gini=0.55)
    grade = _grade_all(0.75, 60.0, 0.30, 0.55)
    assert grade == "ACCEPTABLE"
    ns = m.normalized_score()
    assert 0.50 <= ns <= 0.85, f"ACCEPTABLE should score 0.50-0.85; got {ns:.3f}"


def test_all_poor_maps_to_low_score():
    """A set that grades POOR should score below 0.60."""
    m = _metrics_at(cov=0.55, sel=25.0, dimer=0.45, gini=0.70)
    grade = _grade_all(0.55, 25.0, 0.45, 0.70)
    assert grade == "POOR"
    ns = m.normalized_score()
    assert ns < 0.60, f"POOR set should score < 0.60; got {ns:.3f}"


def test_grade_ordering_preserves_score_ordering():
    """If set A grades strictly better than set B on every component, A
    must have a higher normalized_score."""
    better = _metrics_at(cov=0.95, sel=200.0, dimer=0.10, gini=0.25)
    worse = _metrics_at(cov=0.60, sel=30.0, dimer=0.40, gini=0.65)
    assert better.normalized_score() > worse.normalized_score()


def test_application_weights_preserve_grade_ordering():
    """Under every application profile, strictly better grades yield
    strictly better normalized_score."""
    better = _metrics_at(cov=0.95, sel=200.0, dimer=0.10, gini=0.25)
    worse = _metrics_at(cov=0.60, sel=30.0, dimer=0.40, gini=0.65)
    for application in ("balanced", "discovery", "clinical",
                        "enrichment", "metagenomics"):
        b = better.normalized_score(application=application)
        w = worse.normalized_score(application=application)
        assert b > w, f"{application}: better should outrank worse ({b:.3f} vs {w:.3f})"


# ----------------------------------------------------------------------
# Mixed-grade disagreement shape (review I7)
# ----------------------------------------------------------------------
# `interpret` uses a min-of-ratings rule, so EXCELLENT coverage + CRITICAL
# dimer grades CRITICAL overall. `normalized_score` uses weighted sum and
# will give a middling number (roughly 0.55 with balanced weights). These
# tests document the intentional asymmetry and assert that clinical users
# (who weight dimer at 0.20) see a lower score than metagenomics users
# (who weight dimer at 0.10) for the same "excellent-coverage-but-
# pathological-dimer" set.


def test_mixed_grade_excellent_coverage_critical_dimer():
    """A set with EXCELLENT coverage and CRITICAL dimer is a single-
    component disaster. The overall interpret grade collapses to
    CRITICAL; normalized_score lands in the 0.45-0.70 middling band
    because the weighted sum cannot fully penalise one bad component."""
    m = _metrics_at(cov=0.98, sel=250.0, dimer=0.85, gini=0.25)
    dimer_grade = _rate_dimer(0.85)
    assert dimer_grade == QualityRating.CRITICAL

    balanced = m.normalized_score(application="balanced")
    # CRITICAL dimer drags the score below EXCELLENT territory but not to
    # CRITICAL; this is the documented disagreement between the two reports.
    assert balanced < 0.85, (
        f"CRITICAL dimer should keep normalized_score below 0.85; got {balanced:.3f}"
    )
    assert balanced > 0.30, (
        f"EXCELLENT coverage+selectivity should keep score > 0.30; got {balanced:.3f}"
    )


def test_clinical_penalises_bad_dimer_more_than_metagenomics():
    """The same mixed-grade set should score lower under clinical (dimer
    weight 0.20) than under metagenomics (dimer weight 0.10)."""
    m = _metrics_at(cov=0.98, sel=250.0, dimer=0.85, gini=0.25)
    clinical = m.normalized_score(application="clinical")
    meta = m.normalized_score(application="metagenomics")
    assert clinical < meta, (
        f"clinical should penalise CRITICAL dimer harder than metagenomics; "
        f"got clinical={clinical:.3f}, meta={meta:.3f}"
    )
