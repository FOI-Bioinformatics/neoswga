"""Coverage + regression tests for neoswga.core.efficiency_predictor.

Includes a guard against the bug where _calculate_ml_score silently returned a
constant 0.5 (it called a nonexistent compute_primer_features() and the
classifier-only predict_proba() on a RandomForestRegressor).
"""

import numpy as np
import pytest

from neoswga.core.efficiency_predictor import (
    EfficiencyPrediction,
    EfficiencyPredictor,
    Recommendation,
)


class _StubCache:
    """Position cache returning a few foreground sites, none in background."""

    def __init__(self, positions=None):
        self._positions = positions or {}

    def get_positions(self, prefix, primer, strand):
        return np.asarray(self._positions.get((prefix, strand), []), dtype=np.int64)


@pytest.fixture
def predictor():
    return EfficiencyPredictor(
        position_cache=_StubCache(
            {
                ("fg", "forward"): [100, 4000, 8000],
                ("fg", "reverse"): [2000, 6000],
            }
        ),
        fg_prefixes=["fg"],
        fg_seq_lengths=[10000],
        bg_prefixes=[],
        bg_seq_lengths=[],
    )


def test_ml_score_is_real_not_constant(predictor):
    """Regression guard: the RF ML score must reflect the model, not 0.5."""
    good = ["ATCGATCGATCG", "GCTAGCTAGCTA", "ACGTACGTACGT"]
    bad = ["AAAAAAAAAAAA", "TTTTTTTTTTTT", "GGGGGGGGGGGG"]
    s_good, e_good = predictor._calculate_ml_score(good, verbose=False)
    s_bad, e_bad = predictor._calculate_ml_score(bad, verbose=False)

    assert 0.0 <= s_good <= 1.0 and 0.0 <= s_bad <= 1.0
    # Not the old constant-0.5 fallback for both inputs.
    assert not (s_good == 0.5 and s_bad == 0.5)
    # The model discriminates balanced primers from homopolymer junk.
    assert s_good > s_bad
    assert e_good > e_bad


def test_predict_empty_returns_do_not_synthesize(predictor):
    pred = predictor.predict([], verbose=False)
    assert isinstance(pred, EfficiencyPrediction)
    assert pred.recommendation == Recommendation.DO_NOT_SYNTHESIZE
    assert pred.confidence_score == 0.0


def test_predict_end_to_end_returns_valid_prediction(predictor):
    pred = predictor.predict(["ATCGATCGATCG", "GCTAGCTAGCTA"], verbose=False)
    assert isinstance(pred, EfficiencyPrediction)
    assert 0.0 <= pred.confidence_score <= 1.0
    assert isinstance(pred.recommendation, Recommendation)
    assert 0.0 <= pred.ml_score <= 1.0
    assert pred.predicted_enrichment >= 0.0
    assert isinstance(pred.limiting_factors, list)


def test_prediction_to_dict_and_str(predictor):
    pred = predictor.predict(["ATCGATCGATCG", "GCTAGCTAGCTA"], verbose=False)
    d = pred.to_dict()
    assert isinstance(d, dict)
    assert "confidence_score" in d and "recommendation" in d
    assert isinstance(str(pred), str) and str(pred)


def test_recommendation_tracks_confidence(predictor):
    """A strong set should not be rated worse than a weak one."""
    strong = predictor.predict(["ATCGATCGATCG", "GCTAGCTAGCTA", "ACGTACGTACGT"], verbose=False)
    weak = predictor.predict(["AAAAAAAAAAAA", "TTTTTTTTTTTT"], verbose=False)
    assert strong.confidence_score >= weak.confidence_score
