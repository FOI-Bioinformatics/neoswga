"""Phase 14C: --application tunes normalized_score weights.

A primer set that scores well under the 'balanced' default may rank
lower for 'clinical' use (which penalises background hits harder) and
higher for 'metagenomics' (which favours coverage).
"""

import pytest

from neoswga.core.base_optimizer import PrimerSetMetrics


def _metrics(fg_cov=0.9, selectivity=50.0, dimer=0.1, gini=0.3, tm_range=(50, 60)):
    return PrimerSetMetrics(
        fg_coverage=fg_cov,
        bg_coverage=0.01,
        coverage_uniformity=1.0 - gini,
        total_fg_sites=100,
        total_bg_sites=max(1, int(100 / selectivity)) if selectivity else 100,
        selectivity_ratio=selectivity,
        mean_tm=(tm_range[0] + tm_range[1]) / 2,
        tm_range=tm_range,
        dimer_risk_score=dimer,
        mean_gap=1000.0,
        max_gap=5000.0,
        gap_gini=gini,
        gap_entropy=0.5,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
    )


def test_balanced_weights_default():
    m = _metrics()
    # balanced is the default behavior
    a = m.normalized_score()
    b = m.normalized_score(application="balanced")
    assert abs(a - b) < 1e-9


def test_clinical_penalises_low_selectivity():
    """Two sets with the same coverage but different selectivity:
    clinical must prefer the high-selectivity one more strongly than
    metagenomics does."""
    high_sel = _metrics(fg_cov=0.90, selectivity=100.0, dimer=0.1)
    low_sel = _metrics(fg_cov=0.90, selectivity=10.0, dimer=0.1)

    clinical_delta = (
        high_sel.normalized_score(application="clinical")
        - low_sel.normalized_score(application="clinical")
    )
    meta_delta = (
        high_sel.normalized_score(application="metagenomics")
        - low_sel.normalized_score(application="metagenomics")
    )

    assert clinical_delta > meta_delta, (
        f"clinical should favour selectivity more than metagenomics; "
        f"got clinical_delta={clinical_delta:.3f}, meta_delta={meta_delta:.3f}"
    )


def test_metagenomics_favours_coverage():
    """Set A = high coverage but weak selectivity.
    Set B = low coverage but strong selectivity.
    Metagenomics should rank A higher; clinical should rank B higher."""
    A = _metrics(fg_cov=0.95, selectivity=20.0, dimer=0.1)
    B = _metrics(fg_cov=0.60, selectivity=100.0, dimer=0.1)

    meta_A = A.normalized_score(application="metagenomics")
    meta_B = B.normalized_score(application="metagenomics")
    clinical_A = A.normalized_score(application="clinical")
    clinical_B = B.normalized_score(application="clinical")

    assert meta_A > meta_B, f"meta should prefer high-coverage A; got {meta_A} vs {meta_B}"
    assert clinical_B > clinical_A, (
        f"clinical should prefer high-selectivity B; got clinical_A={clinical_A}, clinical_B={clinical_B}"
    )


def test_discovery_weights_emphasize_coverage():
    """Discovery weights favour coverage most (0.50)."""
    m = _metrics(fg_cov=0.95, selectivity=30.0, dimer=0.1)
    disc = m.normalized_score(application="discovery")
    bal = m.normalized_score(application="balanced")
    # Weights differ but both should be in [0, 1]
    assert 0.0 <= disc <= 1.0
    assert 0.0 <= bal <= 1.0
    # With 95% coverage and moderate everything else, discovery should
    # score at least as high as balanced (cov_w bumped from 0.35 to 0.50).
    assert disc >= bal - 1e-6


def test_unknown_application_falls_back_to_defaults():
    m = _metrics()
    a = m.normalized_score()
    b = m.normalized_score(application="does-not-exist")
    assert abs(a - b) < 1e-9
