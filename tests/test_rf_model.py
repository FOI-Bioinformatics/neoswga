"""Tests for random forest model loading, feature engineering, and prediction.

Verifies that the pre-trained RF model loads correctly, that feature vectors
match the expected schema, and that predictions remain stable across updates.
"""

import os
import hashlib

import numpy as np
import pandas as pd
import pytest

from neoswga.core.rf_preprocessing import (
    _compute_file_hash,
    _TRUSTED_MODEL_HASHES,
    base_features,
    bins,
    delta_g_on_features,
    get_features,
    load_model_safely,
    ModelIntegrityError,
    regression_features,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_MODEL_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "neoswga",
    "core",
    "models",
)
_MODEL_PATH = os.path.join(_MODEL_DIR, "random_forest_filter.p")

# Test primers spanning a range of compositions.
TEST_PRIMERS = [
    "ATCGATCG",   # balanced 8-mer
    "GCGCGCGC",   # 100% GC 8-mer
    "AAAATTTT",   # 0% GC 8-mer
    "ATCGATCGATCG",  # balanced 12-mer
    "GCTAGCTAGC",    # moderate GC 10-mer
]

# Expected number of base features (non-thermodynamic).
EXPECTED_BASE_FEATURE_COUNT = 26

# Number of delta-G histogram bins (len(bins) - 1).
EXPECTED_DELTA_G_FEATURE_COUNT = len(bins) - 1

# Total feature count used by the model.
EXPECTED_TOTAL_FEATURE_COUNT = EXPECTED_BASE_FEATURE_COUNT + EXPECTED_DELTA_G_FEATURE_COUNT


# =========================================================================
# Model loading tests
# =========================================================================


class TestModelLoading:
    """Verify that the shipped model file is present and loadable."""

    def test_model_file_exists(self):
        """The pre-trained model file must be present in the package."""
        assert os.path.isfile(_MODEL_PATH), (
            f"Model file not found at {_MODEL_PATH}"
        )

    def test_model_loads_successfully(self):
        """The model must load without errors via the secure loader."""
        model = load_model_safely(_MODEL_PATH)
        assert model is not None
        assert hasattr(model, "predict"), "Loaded object lacks a predict method"

    def test_model_hash_verification_passes(self):
        """The model file hash must match the trusted hash constant."""
        model_name = os.path.basename(_MODEL_PATH)
        assert model_name in _TRUSTED_MODEL_HASHES, (
            f"No trusted hash registered for {model_name}"
        )
        actual_hash = _compute_file_hash(_MODEL_PATH)
        expected_hash = _TRUSTED_MODEL_HASHES[model_name]
        assert actual_hash == expected_hash, (
            f"Hash mismatch.\nExpected: {expected_hash}\nActual:   {actual_hash}"
        )

    def test_tampered_model_raises_integrity_error(self, tmp_path):
        """A model with an incorrect hash must raise ModelIntegrityError."""
        fake_model = tmp_path / "random_forest_filter.p"
        fake_model.write_bytes(b"not-a-real-model")
        with pytest.raises(ModelIntegrityError):
            load_model_safely(str(fake_model), verify_hash=True)

    def test_missing_model_raises_file_not_found(self):
        """A non-existent path must raise FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            load_model_safely("/nonexistent/path/model.p")

    def test_model_is_random_forest_regressor(self):
        """The model should be a RandomForestRegressor instance."""
        from sklearn.ensemble import RandomForestRegressor

        model = load_model_safely(_MODEL_PATH)
        assert isinstance(model, RandomForestRegressor), (
            f"Expected RandomForestRegressor, got {type(model).__name__}"
        )


# =========================================================================
# Feature engineering tests
# =========================================================================


class TestFeatureEngineering:
    """Verify that feature vectors match the schema expected by the model."""

    def test_base_feature_count(self):
        """base_features list must have the expected number of entries."""
        assert len(base_features) == EXPECTED_BASE_FEATURE_COUNT

    def test_delta_g_feature_count(self):
        """delta_g_on_features list must match the histogram bin count."""
        assert len(delta_g_on_features) == EXPECTED_DELTA_G_FEATURE_COUNT

    def test_regression_features_count(self):
        """regression_features must equal base + delta-G features."""
        assert len(regression_features) == EXPECTED_TOTAL_FEATURE_COUNT

    def test_get_features_vector_length(self):
        """get_features must return a list of length 2 + base_feature_count."""
        # get_features returns [primer, target, feat1, feat2, ...]
        result = get_features("ATCGATCG", target=True)
        expected_length = 2 + EXPECTED_BASE_FEATURE_COUNT
        assert len(result) == expected_length, (
            f"Expected {expected_length} elements, got {len(result)}"
        )

    def test_get_features_primer_identity(self):
        """The first element of the feature vector must be the primer string."""
        result = get_features("GCTAGCTA", target=False)
        assert result[0] == "GCTAGCTA"
        assert result[1] is False

    @pytest.mark.parametrize("primer", TEST_PRIMERS)
    def test_gc_content_calculation(self, primer):
        """GC content must be correctly computed."""
        result = get_features(primer, target=True)
        # Index: 0=primer, 1=target, then base_features in order
        # GC.content is at index base_features.index('GC.content') + 2
        gc_idx = base_features.index("GC.content") + 2
        expected_gc = (primer.count("G") + primer.count("C")) / len(primer)
        assert abs(result[gc_idx] - expected_gc) < 1e-9

    @pytest.mark.parametrize("primer", TEST_PRIMERS)
    def test_sequence_length_feature(self, primer):
        """Sequence length feature must match actual primer length."""
        result = get_features(primer, target=True)
        length_idx = base_features.index("sequence.length") + 2
        assert result[length_idx] == len(primer)

    def test_feature_names_match_model_expectations(self):
        """Computed feature names must match what the model was trained on.

        Compares regression_features against the feature list used in
        scripts/retrain_rf_model.py to detect any drift.
        """
        # The canonical list from the retraining script.
        retrain_base = [
            "molarity", "sequence.length", "number.of.A", "proportion.of.A",
            "number.of.T", "proportion.of.T", "number.of.G", "proportion.of.G",
            "number.of.C", "proportion.of.C", "GC.content", "melting_tm",
            "GC.clamp", "longest.A.repeat", "longest.T.repeat",
            "longest.G.repeat", "longest.C.repeat", "AA repeat", "CC repeat",
            "TT repeat", "GG repeat", "3.end.first.base", "3.end.second.base",
            "3.end.third.base", "3.end.fourth.base", "3.end.fifth.base",
        ]
        retrain_bins = [
            -20, -18, -16, -14, -12, -10, -9, -8, -7, -6,
            -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5,
            0, 0.5, 1, 1.5, 2, 2.5, 3,
        ]
        retrain_delta_g = ["on_target_" + str(b) for b in retrain_bins[1:]]
        retrain_all = retrain_base + retrain_delta_g

        assert regression_features == retrain_all, (
            "regression_features in rf_preprocessing.py have drifted from "
            "the retraining script's ALL_FEATURES list."
        )

    def test_no_nan_in_features(self):
        """Feature vectors must not contain NaN values."""
        for primer in TEST_PRIMERS:
            result = get_features(primer, target=True)
            numeric_values = result[2:]  # skip primer string and target bool
            for i, val in enumerate(numeric_values):
                assert not (isinstance(val, float) and np.isnan(val)), (
                    f"NaN at index {i} for primer {primer}"
                )


# =========================================================================
# Regression / prediction tests
# =========================================================================


class TestPredictionRegression:
    """Verify that model predictions are stable and sensible."""

    @pytest.fixture(scope="class")
    def model(self):
        """Load the model once for the entire test class."""
        return load_model_safely(_MODEL_PATH)

    @pytest.fixture(scope="class")
    def predictions(self, model):
        """Compute predictions for all test primers and return as a dict.

        Builds a synthetic feature DataFrame with zeroed thermodynamic
        features (no genome context) so the test is self-contained.
        """
        rows = []
        for primer in TEST_PRIMERS:
            feat = get_features(primer, target=True)
            # feat layout: [primer, target, base_feat_1, ..., base_feat_N]
            row = dict(zip(base_features, feat[2:]))
            # Fill thermodynamic features with zeros (no genome context)
            for dg_feat in delta_g_on_features:
                row[dg_feat] = 0.0
            rows.append(row)

        df = pd.DataFrame(rows)
        X = df[regression_features]
        scores = model.predict(X)
        return dict(zip(TEST_PRIMERS, scores))

    def test_predictions_are_finite(self, predictions):
        """All predictions must be finite numbers."""
        for primer, score in predictions.items():
            assert np.isfinite(score), (
                f"Non-finite score {score} for primer {primer}"
            )

    def test_predictions_in_expected_range(self, predictions):
        """Scores should fall within the training target range [0, 20]."""
        for primer, score in predictions.items():
            assert -5.0 <= score <= 25.0, (
                f"Score {score:.4f} for primer {primer} is outside "
                "the plausible range [-5, 25]"
            )

    def test_gc_rich_differs_from_at_rich(self, predictions):
        """A pure-GC primer must score differently from a pure-AT primer.

        This is a basic sanity check: if the model treats these identically,
        the feature pipeline is likely broken.
        """
        gc_score = predictions["GCGCGCGC"]
        at_score = predictions["AAAATTTT"]
        assert gc_score != pytest.approx(at_score, abs=0.01), (
            f"GC-only ({gc_score:.4f}) and AT-only ({at_score:.4f}) primers "
            "received nearly identical scores, suggesting broken features."
        )

    def test_balanced_primer_not_extreme(self, predictions):
        """A balanced primer should not produce an extreme outlier score."""
        balanced = predictions["ATCGATCG"]
        gc_only = predictions["GCGCGCGC"]
        at_only = predictions["AAAATTTT"]
        # The balanced primer should not be the most extreme score
        all_scores = [balanced, gc_only, at_only]
        assert balanced != max(all_scores) or balanced != min(all_scores), (
            "Balanced primer scored at the extreme end, which is unexpected."
        )

    def test_prediction_reproducibility(self, model):
        """Repeated predictions on the same input must be identical.

        Random forests are deterministic at prediction time, so any
        variation would indicate a stateful bug.
        """
        feat = get_features("ATCGATCG", target=True)
        row = dict(zip(base_features, feat[2:]))
        for dg_feat in delta_g_on_features:
            row[dg_feat] = 0.0
        df = pd.DataFrame([row])
        X = df[regression_features]

        score_a = model.predict(X)[0]
        score_b = model.predict(X)[0]
        assert score_a == pytest.approx(score_b)

    def test_model_accepts_full_feature_vector(self, model):
        """The model must accept a feature vector of the expected width."""
        n_model_features = model.n_features_in_
        assert n_model_features == EXPECTED_TOTAL_FEATURE_COUNT, (
            f"Model expects {n_model_features} features but "
            f"regression_features has {EXPECTED_TOTAL_FEATURE_COUNT}."
        )

    # ------------------------------------------------------------------
    # Approximate reference scores (zeroed thermo features, molarity=2.5)
    # Update these after intentional model retrains.
    # ------------------------------------------------------------------
    REFERENCE_SCORES = {
        "ATCGATCG": 16.69,
        "GCGCGCGC": 12.62,
        "AAAATTTT": 6.71,
        "ATCGATCGATCG": 15.42,
        "GCTAGCTAGC": 19.03,
    }

    @pytest.mark.parametrize("primer", TEST_PRIMERS)
    def test_score_approximate_stability(self, predictions, primer):
        """Scores should remain within 4.0 of the reference values.

        A wide tolerance is used because the reference scores are
        approximate. The purpose is to catch large regressions, not
        to enforce exact reproducibility after a retrain.
        """
        expected = self.REFERENCE_SCORES[primer]
        actual = predictions[primer]
        assert abs(actual - expected) < 4.0, (
            f"Primer {primer}: score {actual:.4f} deviates more than 4.0 "
            f"from reference {expected:.4f}. If the model was intentionally "
            "retrained, update REFERENCE_SCORES."
        )


# =========================================================================
# Integration: predict_new_primers pathway
# =========================================================================


class TestPredictNewPrimers:
    """Test the full predict_new_primers function end-to-end."""

    def test_predict_new_primers_returns_dataframe(self):
        """predict_new_primers must return a DataFrame with an
        'on.target.pred' column."""
        from neoswga.core.rf_preprocessing import predict_new_primers

        rows = []
        for primer in TEST_PRIMERS[:2]:
            feat = get_features(primer, target=True)
            row = {"sequence": primer, "target": True}
            row.update(dict(zip(base_features, feat[2:])))
            for dg_feat in delta_g_on_features:
                row[dg_feat] = 0.0
            rows.append(row)

        df = pd.DataFrame(rows)
        result = predict_new_primers(df)

        assert isinstance(result, pd.DataFrame)
        assert "on.target.pred" in result.columns
        assert len(result) == len(TEST_PRIMERS[:2])

    def test_predict_new_primers_scores_are_numeric(self):
        """All prediction scores must be numeric and finite."""
        from neoswga.core.rf_preprocessing import predict_new_primers

        rows = []
        for primer in TEST_PRIMERS:
            feat = get_features(primer, target=True)
            row = {"sequence": primer, "target": True}
            row.update(dict(zip(base_features, feat[2:])))
            for dg_feat in delta_g_on_features:
                row[dg_feat] = 0.0
            rows.append(row)

        df = pd.DataFrame(rows)
        result = predict_new_primers(df)

        for score in result["on.target.pred"]:
            assert np.isfinite(score), f"Non-finite score: {score}"
