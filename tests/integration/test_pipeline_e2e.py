"""End-to-end integration test using plasmid example data.

Tests the filter -> score -> optimize pipeline with real genome data.
Requires pre-generated k-mer count files (no jellyfish dependency).
"""
import json
import os
import shutil
import tempfile
import pytest
import pandas as pd


EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'plasmid_example')


@pytest.fixture
def pipeline_workdir():
    """Create a temp copy of the plasmid example for testing."""
    if not os.path.isdir(EXAMPLE_DIR):
        pytest.skip("Plasmid example data not available")

    tmpdir = tempfile.mkdtemp(prefix='neoswga_integration_')
    for fname in os.listdir(EXAMPLE_DIR):
        src = os.path.join(EXAMPLE_DIR, fname)
        if os.path.isfile(src):
            shutil.copy2(src, tmpdir)

    # Ensure optional params that the example params.json omits are present,
    # so the pipeline uses sensible defaults rather than receiving None.
    params_path = os.path.join(tmpdir, 'params.json')
    with open(params_path) as fh:
        params_data = json.load(fh)
    params_data.setdefault('min_amp_pred', 0.0)  # Accept all primers for small plasmid
    params_data.setdefault('schema_version', 1)
    with open(params_path, 'w') as fh:
        json.dump(params_data, fh, indent=2)

    original_cwd = os.getcwd()
    os.chdir(tmpdir)

    yield tmpdir

    os.chdir(original_cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.mark.integration
@pytest.mark.slow
class TestPipelineE2E:
    """End-to-end pipeline integration tests."""

    def test_step2_filter_produces_output(self, pipeline_workdir, reset_pipeline_state):
        """Step 2 (filter) should produce a non-empty step2_df.csv."""
        step2_output = os.path.join(pipeline_workdir, 'step2_df.csv')
        if os.path.exists(step2_output):
            os.remove(step2_output)

        params_file = os.path.join(pipeline_workdir, 'params.json')
        reset_pipeline_state(params_file)

        from neoswga.core.pipeline import step2
        df = step2()

        assert os.path.exists(step2_output), "step2_df.csv was not created"
        assert len(df) > 0, "Filter produced no candidate primers"
        assert 'primer' in df.columns, f"Expected 'primer' column, got: {list(df.columns)}"

    def test_step2_primers_are_valid_dna(self, pipeline_workdir):
        """All primers in step2_df.csv should be valid DNA sequences."""
        step2_csv = os.path.join(pipeline_workdir, 'step2_df.csv')
        assert os.path.exists(step2_csv), "step2_df.csv required — run filter step first"

        df = pd.read_csv(step2_csv)
        valid_bases = set('ATCG')
        for primer in df['primer'].head(20):
            assert set(primer.upper()).issubset(valid_bases), \
                f"Invalid DNA sequence in step2 output: {primer}"

    def test_step3_score_produces_output(self, pipeline_workdir, reset_pipeline_state):
        """Step 3 (score) should produce a non-empty step3_df.csv."""
        step3_output = os.path.join(pipeline_workdir, 'step3_df.csv')
        if os.path.exists(step3_output):
            os.remove(step3_output)

        step2_csv = os.path.join(pipeline_workdir, 'step2_df.csv')
        assert os.path.exists(step2_csv), "step2_df.csv required for step3"

        params_file = os.path.join(pipeline_workdir, 'params.json')
        reset_pipeline_state(params_file)

        from neoswga.core.pipeline import step3
        df = step3()

        assert os.path.exists(step3_output), "step3_df.csv was not created"
        assert len(df) > 0, "Score step produced no primers"
        # step3 returns a DataFrame with primer sequences as the index
        assert 'on.target.pred' in df.columns, \
            f"Expected 'on.target.pred' column in step3 output, got: {list(df.columns)}"

    def test_step3_has_amp_pred_scores(self, pipeline_workdir):
        """step3_df.csv should include amplification prediction scores."""
        step3_csv = os.path.join(pipeline_workdir, 'step3_df.csv')
        assert os.path.exists(step3_csv), "step3_df.csv required"

        df = pd.read_csv(step3_csv, index_col=0)
        assert 'on.target.pred' in df.columns, \
            f"Expected 'on.target.pred' column in step3 output, got: {list(df.columns)}"
        assert df['on.target.pred'].notna().any(), "All amp_pred scores are NaN"

    def test_step4_optimize_produces_output(self, pipeline_workdir, reset_pipeline_state):
        """Step 4 (optimize) should produce a non-empty step4_improved_df.csv."""
        step4_output = os.path.join(pipeline_workdir, 'step4_improved_df.csv')
        if os.path.exists(step4_output):
            os.remove(step4_output)

        step3_csv = os.path.join(pipeline_workdir, 'step3_df.csv')
        assert os.path.exists(step3_csv), "step3_df.csv required for step4"

        params_file = os.path.join(pipeline_workdir, 'params.json')
        reset_pipeline_state(params_file)

        from neoswga.core.unified_optimizer import optimize_step4
        primer_sets, scores, _cache = optimize_step4(
            optimization_method='dominating-set',
            verbose=False,
        )

        assert os.path.exists(step4_output), "step4_improved_df.csv was not created"
        assert len(primer_sets) > 0, "Optimize step returned no primer sets"
        assert len(primer_sets[0]) > 0, "Best primer set is empty"

    def test_step4_primers_are_valid_dna(self, pipeline_workdir):
        """Optimized primers in step4_improved_df.csv should be valid DNA."""
        step4_csv = os.path.join(pipeline_workdir, 'step4_improved_df.csv')
        assert os.path.exists(step4_csv), "step4_improved_df.csv required"

        df = pd.read_csv(step4_csv)
        assert 'primer' in df.columns, f"Expected 'primer' column, got: {list(df.columns)}"

        valid_bases = set('ATCG')
        for primer in df['primer']:
            assert set(primer.upper()).issubset(valid_bases), \
                f"Invalid DNA sequence in step4 output: {primer}"

    def test_full_pipeline_runs_sequentially(self, pipeline_workdir, reset_pipeline_state):
        """Run steps 2 -> 3 -> 4 end-to-end from scratch."""
        for fname in ['step2_df.csv', 'step3_df.csv', 'step4_improved_df.csv']:
            path = os.path.join(pipeline_workdir, fname)
            if os.path.exists(path):
                os.remove(path)

        params_file = os.path.join(pipeline_workdir, 'params.json')
        reset_pipeline_state(params_file)

        from neoswga.core.pipeline import step2, step3
        from neoswga.core.unified_optimizer import optimize_step4

        df2 = step2()
        assert len(df2) > 0, "Step 2 produced no primers"

        reset_pipeline_state(params_file)
        df3 = step3()
        assert len(df3) > 0, "Step 3 produced no primers"

        reset_pipeline_state(params_file)
        primer_sets, scores, _cache = optimize_step4(
            optimization_method='dominating-set',
            verbose=False,
        )
        assert len(primer_sets) > 0 and len(primer_sets[0]) > 0, \
            "Step 4 produced no primers"

        step4_csv = os.path.join(pipeline_workdir, 'step4_improved_df.csv')
        df4 = pd.read_csv(step4_csv)
        assert len(df4) > 0, "step4_improved_df.csv is empty"
        assert 'primer' in df4.columns
        assert 'score' in df4.columns
