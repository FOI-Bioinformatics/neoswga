"""Coverage for neoswga.core.pipeline_qa_integration scoring helpers."""

import pandas as pd

from neoswga.core.pipeline_qa_integration import QAFilterConfig, combine_rf_qa_scores


def test_qa_filter_config_defaults():
    cfg = QAFilterConfig()
    assert 0.0 <= cfg.qa_weight <= 1.0
    assert abs((cfg.qa_weight + cfg.rf_weight) - 1.0) < 1e-9
    assert cfg.enable_three_prime is True


def test_combine_rf_qa_scores_adds_composite():
    df = pd.DataFrame({"primer": ["ATCGATCG", "GCTAGCTA", "TTTTAAAA"], "score": [12.0, 6.0, 0.0]})
    qa_scores = {"ATCGATCG": 0.9, "GCTAGCTA": 0.5}  # third primer absent -> default 0.5
    out = combine_rf_qa_scores(df, qa_scores, rf_weight=0.7, qa_weight=0.3, verbose=False)

    assert "qa_score" in out.columns
    assert "composite_score" in out.columns
    # absent primer falls back to the 0.5 default QA score
    assert out.loc[out["primer"] == "TTTTAAAA", "qa_score"].iloc[0] == 0.5
    # RF normalised to [0,1]: best RF score (12.0) maps to 1.0
    best = out.loc[out["primer"] == "ATCGATCG"].iloc[0]
    assert best["rf_score_norm"] == 1.0
    # composite is the weighted blend
    expected = 0.7 * 1.0 + 0.3 * 0.9
    assert abs(best["composite_score"] - expected) < 1e-9


def test_combine_handles_flat_rf_scores():
    # All-equal RF scores must not divide by zero (rf range == 0).
    df = pd.DataFrame({"primer": ["ATCGATCG", "GCTAGCTA"], "score": [5.0, 5.0]})
    out = combine_rf_qa_scores(df, {"ATCGATCG": 0.8}, verbose=False)
    assert out["composite_score"].notna().all()
    assert (out["rf_score_norm"] == 1.0).all()
