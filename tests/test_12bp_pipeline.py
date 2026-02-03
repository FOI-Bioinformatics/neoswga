"""
Pytest tests for 12bp pipeline validation.

Tests the output quality of the Francisella 12bp SWGA primer design pipeline.
Run after executing the pipeline with:
    cd tests/integration/francisella_e2e/12bp
    neoswga filter -j params.json
    neoswga score -j params.json
    neoswga optimize -j params.json
"""

import pytest
import pandas as pd
import os
from pathlib import Path

# Test directory relative to repo root
TEST_DIR = Path(__file__).parent / "integration" / "francisella_e2e" / "12bp"


class Test12bpPipelineResults:
    """Tests for 12bp pipeline output quality."""

    @pytest.fixture(scope="class")
    def step2_df(self):
        """Load filter step results."""
        path = TEST_DIR / "step2_df.csv"
        if not path.exists():
            pytest.skip(f"step2_df.csv not found at {path}")
        return pd.read_csv(path)

    @pytest.fixture(scope="class")
    def step3_df(self):
        """Load score step results."""
        path = TEST_DIR / "step3_df.csv"
        if not path.exists():
            pytest.skip(f"step3_df.csv not found at {path}")
        return pd.read_csv(path)

    @pytest.fixture(scope="class")
    def step4_df(self):
        """Load optimize step results."""
        path = TEST_DIR / "step4_improved_df.csv"
        if not path.exists():
            pytest.skip(f"step4_improved_df.csv not found at {path}")
        return pd.read_csv(path)

    # ============================================================
    # Primer Set Tests
    # ============================================================

    def test_primer_count(self, step4_df):
        """Pipeline should select at least 6 primers."""
        assert len(step4_df) >= 6, f"Expected at least 6 primers, got {len(step4_df)}"

    def test_primer_length(self, step4_df):
        """All primers should be exactly 12bp."""
        for primer in step4_df["primer"]:
            assert len(primer) == 12, f"Primer {primer} is {len(primer)}bp, expected 12bp"

    def test_primer_valid_bases(self, step4_df):
        """All primers should only contain A, T, G, C."""
        valid_bases = {"A", "T", "G", "C"}
        for primer in step4_df["primer"]:
            invalid = set(primer) - valid_bases
            assert not invalid, f"Primer {primer} contains invalid bases: {invalid}"

    def test_primers_unique(self, step4_df):
        """All selected primers should be unique."""
        primers = step4_df["primer"].tolist()
        assert len(primers) == len(set(primers)), "Duplicate primers detected"

    # ============================================================
    # Coverage Tests
    # ============================================================

    def test_coverage_threshold(self, step4_df):
        """Coverage should be at least 70%."""
        if "coverage" not in step4_df.columns:
            pytest.skip("coverage column not present")
        coverage = step4_df["coverage"].iloc[0]
        assert coverage >= 0.70, f"Coverage {coverage:.1%} is below 70% threshold"

    def test_no_background_binding(self, step4_df):
        """Background coverage should be less than 1%."""
        if "bg_coverage" not in step4_df.columns:
            pytest.skip("bg_coverage column not present")
        bg_cov = step4_df["bg_coverage"].iloc[0]
        assert bg_cov < 0.01, f"Background coverage {bg_cov:.1%} exceeds 1%"

    # ============================================================
    # Quality Score Tests
    # ============================================================

    def test_amplification_score(self, step4_df):
        """Amplification score should be at least 1000x."""
        if "score" not in step4_df.columns:
            pytest.skip("score column not present")
        min_score = step4_df["score"].min()
        assert min_score >= 1000, f"Min amplification {min_score:.0f}x is below 1000x"

    def test_selectivity(self, step4_df):
        """Selectivity should be positive (foreground > background)."""
        if "selectivity" not in step4_df.columns:
            pytest.skip("selectivity column not present")
        selectivity = step4_df["selectivity"].iloc[0]
        assert selectivity > 1, f"Selectivity {selectivity:.1f}x indicates poor enrichment"

    # ============================================================
    # Filter Quality Tests
    # ============================================================

    def test_gini_threshold(self, step2_df):
        """Gini index should not exceed 0.6 for even binding."""
        max_gini = step2_df["gini"].max()
        assert max_gini <= 0.6, f"Max Gini {max_gini:.2f} exceeds 0.6 threshold"

    def test_filter_reduces_candidates(self, step2_df):
        """Filter step should reduce candidates to target count."""
        assert len(step2_df) <= 500, f"Filter output {len(step2_df)} exceeds max_primer=500"

    # ============================================================
    # Scoring Quality Tests
    # ============================================================

    def test_amplification_prediction(self, step3_df):
        """All scored primers should have amplification prediction."""
        if "on.target.pred" in step3_df.columns:
            col = "on.target.pred"
        elif "amp_pred" in step3_df.columns:
            col = "amp_pred"
        else:
            pytest.skip("No amplification prediction column found")

        # Check for valid predictions
        assert step3_df[col].notna().all(), "Some primers missing amplification prediction"

    def test_scoring_filters_poor_candidates(self, step2_df, step3_df):
        """Scoring step may filter out poorly predicted primers."""
        # This is informational - scoring can reduce candidates
        n_before = len(step2_df)
        n_after = len(step3_df)
        assert n_after > 0, "Scoring eliminated all candidates"


class Test12bpPipelineParams:
    """Tests for pipeline parameters and configuration."""

    @pytest.fixture(scope="class")
    def params_path(self):
        """Get params.json path."""
        path = TEST_DIR / "params.json"
        if not path.exists():
            pytest.skip(f"params.json not found at {path}")
        return path

    def test_params_file_exists(self, params_path):
        """params.json should exist."""
        assert params_path.exists()

    def test_params_valid_json(self, params_path):
        """params.json should be valid JSON."""
        import json
        with open(params_path) as f:
            params = json.load(f)
        assert isinstance(params, dict)

    def test_bloom_filter_enabled(self, params_path):
        """Bloom filter should be enabled for performance."""
        import json
        with open(params_path) as f:
            params = json.load(f)

        # Either use_bloom_filter or bg_bloom should be set
        has_bloom = params.get("use_bloom_filter", False) or params.get("bg_bloom")
        assert has_bloom, "Bloom filter not configured (slow performance)"

    def test_kmer_length_correct(self, params_path):
        """K-mer length should be 12bp."""
        import json
        with open(params_path) as f:
            params = json.load(f)

        min_k = params.get("min_k", 12)
        max_k = params.get("max_k", 12)
        assert min_k == 12, f"min_k should be 12, got {min_k}"
        assert max_k == 12, f"max_k should be 12, got {max_k}"
