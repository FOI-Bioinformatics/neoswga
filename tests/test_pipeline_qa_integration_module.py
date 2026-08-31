"""Coverage for neoswga.core.pipeline_qa_integration scoring helpers."""

import pandas as pd
import pytest

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


# ----------------------------------------------------------------------
# Step wrappers used by `--enable-qa`
# ----------------------------------------------------------------------


def _step2_frame():
    return pd.DataFrame(
        {
            "primer": [
                "ATCGATCGATCG",
                "GCTAGCTAGCTA",
                "AAATTTAAATTT",
                "ACGTACGTACGT",
            ],
            "fg_freq": [1e-4] * 4,
            "bg_freq": [1e-7] * 4,
            "gini": [0.3] * 4,
        }
    )


def test_apply_qa_to_step2_output_rewrites_the_csv(tmp_path):
    from neoswga.core.pipeline_qa_integration import apply_qa_to_step2_output

    _step2_frame().to_csv(tmp_path / "step2_df.csv")

    result = apply_qa_to_step2_output(str(tmp_path), verbose=False)

    written = pd.read_csv(tmp_path / "step2_df.csv", index_col=0)
    assert "qa_score" in written.columns
    assert len(written) == result.primers_out
    assert set(written["primer"]) <= set(_step2_frame()["primer"])
    assert (tmp_path / "qa_report.txt").is_file()


def test_apply_qa_to_step2_output_corrects_the_filter_funnel(tmp_path):
    """The funnel is written before QA runs, so its last stage goes stale."""
    import json

    from neoswga.core.pipeline_qa_integration import apply_qa_to_step2_output

    _step2_frame().to_csv(tmp_path / "step2_df.csv")
    (tmp_path / "filter_stats.json").write_text(
        json.dumps({"total_kmers": 100, "after_gini": 4, "final_candidates": 4})
    )

    result = apply_qa_to_step2_output(str(tmp_path), verbose=False)

    stats = json.loads((tmp_path / "filter_stats.json").read_text())
    assert stats["final_candidates"] == result.primers_out
    assert stats["after_qa"] == result.primers_out
    assert stats["total_kmers"] == 100, "the earlier stages were rewritten"


def test_apply_qa_to_step2_output_without_a_funnel_file(tmp_path):
    from neoswga.core.pipeline_qa_integration import apply_qa_to_step2_output

    _step2_frame().to_csv(tmp_path / "step2_df.csv")
    apply_qa_to_step2_output(str(tmp_path), verbose=False)
    assert not (tmp_path / "filter_stats.json").exists()


def test_apply_qa_to_step2_output_keeps_the_pool_when_qa_rejects_all(tmp_path):
    """An empty pool is a failed run, not a silently emptied CSV."""
    from neoswga.core.pipeline_qa_integration import QAFilterConfig, apply_qa_to_step2_output

    _step2_frame().to_csv(tmp_path / "step2_df.csv")
    config = QAFilterConfig(max_hub_degree=-1)  # every primer counts as a hub

    with pytest.raises(ValueError, match="QA"):
        apply_qa_to_step2_output(str(tmp_path), config=config, verbose=False)

    kept = pd.read_csv(tmp_path / "step2_df.csv", index_col=0)
    assert len(kept) == 4, "the pre-QA candidates were overwritten"


def test_apply_qa_to_step3_output_uses_step2_qa_scores(tmp_path):
    from neoswga.core.pipeline_qa_integration import apply_qa_to_step3_output

    step2 = _step2_frame()
    step2["qa_score"] = [0.9, 0.8, 0.7, 0.6]
    step2.to_csv(tmp_path / "step2_df.csv")

    step3 = step2[["primer"]].copy()
    step3["on.target.pred"] = [12.0, 6.0, 3.0, 0.0]
    step3.set_index("primer").to_csv(tmp_path / "step3_df.csv")

    out = apply_qa_to_step3_output(str(tmp_path), verbose=False)

    written = pd.read_csv(tmp_path / "step3_df.csv")
    assert "composite_score" in written.columns
    assert written["composite_score"].is_monotonic_decreasing
    assert out.loc[out["primer"] == "ATCGATCGATCG", "qa_score"].iloc[0] == 0.9


def test_apply_qa_to_step3_output_scores_primers_itself_when_step2_has_none(tmp_path):
    """`score --enable-qa` after a plain `filter` must still produce QA scores."""
    from neoswga.core.pipeline_qa_integration import (
        QAAwareOptimizer,
        QAFilterConfig,
        apply_qa_to_step3_output,
    )
    from neoswga.core.reaction_conditions import get_standard_conditions

    step2 = _step2_frame()
    step2.to_csv(tmp_path / "step2_df.csv")
    step3 = step2[["primer"]].copy()
    step3["on.target.pred"] = [12.0, 6.0, 3.0, 0.0]
    step3.set_index("primer").to_csv(tmp_path / "step3_df.csv")

    # Conditions are passed explicitly: the helper otherwise resolves them from
    # the parameter globals, which another test in the session may have set.
    conditions = get_standard_conditions()
    out = apply_qa_to_step3_output(str(tmp_path), conditions=conditions, verbose=False)

    expected = dict(
        QAAwareOptimizer(QAFilterConfig(), conditions).rank_by_quality(
            step2["primer"].tolist(), verbose=False
        )
    )
    for primer, score in expected.items():
        assert out.loc[out["primer"] == primer, "qa_score"].iloc[0] == pytest.approx(score)


def test_combine_rf_qa_scores_reads_the_pipeline_score_column():
    """step3_df.csv names the RF score `on.target.pred`, not `score`."""
    df = pd.DataFrame({"primer": ["ATCGATCG", "GCTAGCTA"], "on.target.pred": [12.0, 6.0]})
    out = combine_rf_qa_scores(df, {"ATCGATCG": 0.9}, verbose=False)
    assert out.loc[out["primer"] == "ATCGATCG", "rf_score_norm"].iloc[0] == 1.0


def test_qa_prefiltered_candidates_is_none_when_the_flag_is_off(tmp_path):
    from types import SimpleNamespace

    from neoswga.core.pipeline_qa_integration import qa_prefiltered_candidates

    module = SimpleNamespace(enable_qa=False, data_dir=str(tmp_path))
    assert qa_prefiltered_candidates(module, verbose=False) is None


def test_qa_prefiltered_candidates_returns_a_subset_of_the_step3_pool(tmp_path):
    from types import SimpleNamespace

    from neoswga.core.pipeline_qa_integration import qa_prefiltered_candidates

    primers = _step2_frame()["primer"].tolist()
    pd.DataFrame({"primer": primers}).set_index("primer").to_csv(tmp_path / "step3_df.csv")

    module = SimpleNamespace(enable_qa=True, data_dir=str(tmp_path))
    candidates = qa_prefiltered_candidates(module, verbose=False)

    assert candidates is not None
    assert set(candidates) <= set(primers)
    assert candidates, "an empty pool must fall back to None, not an empty list"
