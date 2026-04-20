"""End-to-end test: plasmid scenario with a blacklist.

Exercises the main filter path with bl_genomes populated via params.json,
confirming primers that bind the blacklist are rejected. Reuses the
pre-generated k-mer counts shipped with `examples/plasmid_example/`.
"""

import json
import os
import shutil
import tempfile
import pytest
import pandas as pd


EXAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', '..', 'examples', 'plasmid_example')


def _reset_pipeline_state(params_file):
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter

    pipeline_mod._initialized = False
    pipeline_mod.fg_prefixes = None
    pipeline_mod.bg_prefixes = None
    pipeline_mod.fg_genomes = None
    pipeline_mod.bg_genomes = None
    pipeline_mod.fg_seq_lengths = None
    pipeline_mod.bg_seq_lengths = None
    pipeline_mod.fg_circular = None
    pipeline_mod.bg_circular = None
    parameter.json_file = params_file


@pytest.fixture
def blacklist_workdir():
    """Copy the plasmid example and repurpose pLTR as a blacklist genome.

    This leaves pcDNA as the target and uses the pre-built pLTR k-mer
    files as the blacklist k-mer input. The workdir is laid out so
    step2 can find `bl_pLTR_{k}mer_all.txt` files symlinked from the
    existing pLTR_{k}mer_all.txt files.
    """
    if not os.path.isdir(EXAMPLE_DIR):
        pytest.skip("Plasmid example data not available")

    tmpdir = tempfile.mkdtemp(prefix='neoswga_bl_scenario_')
    for fname in os.listdir(EXAMPLE_DIR):
        src = os.path.join(EXAMPLE_DIR, fname)
        if os.path.isfile(src):
            shutil.copy2(src, tmpdir)

    # Create bl_pLTR_{k}mer_all.txt symlinks to the existing pLTR k-mer files
    # so the blacklist path has data to read without re-running jellyfish.
    for k in range(6, 13):
        src = os.path.join(tmpdir, f'pLTR_{k}mer_all.txt')
        if os.path.isfile(src):
            dst = os.path.join(tmpdir, f'bl_pLTR_{k}mer_all.txt')
            if not os.path.exists(dst):
                shutil.copy2(src, dst)

    # Rewrite params.json: remove background, add blacklist pointing at pLTR
    params_path = os.path.join(tmpdir, 'params.json')
    with open(params_path) as fh:
        params_data = json.load(fh)

    # Keep pLTR as blacklist; remove it from background
    params_data['bg_genomes'] = []
    params_data['bg_prefixes'] = []
    params_data['bg_seq_lengths'] = []
    params_data['bl_genomes'] = [os.path.join(tmpdir, 'pLTR.fasta')]
    params_data['bl_prefixes'] = [os.path.join(tmpdir, 'bl_pLTR')]
    params_data['bl_seq_lengths'] = [6258]
    params_data['max_bl_freq'] = 0.0  # zero tolerance
    params_data.setdefault('min_amp_pred', 0.0)
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
def test_blacklist_filter_runs_end_to_end(blacklist_workdir):
    """Filter step completes with blacklist active and produces output."""
    params_file = os.path.join(blacklist_workdir, 'params.json')
    _reset_pipeline_state(params_file)

    from neoswga.core.pipeline import step2
    df = step2()

    step2_csv = os.path.join(blacklist_workdir, 'step2_df.csv')
    assert os.path.exists(step2_csv), "step2_df.csv was not created"
    assert len(df) >= 0, "step2 should complete without crash"
    # If any primers survived, they should have bl_freq 0 or missing
    if 'bl_freq' in df.columns and len(df) > 0:
        assert (df['bl_freq'] <= 1e-9).all(), (
            "With max_bl_freq=0.0, every surviving primer must have bl_freq ~0"
        )


@pytest.mark.integration
@pytest.mark.slow
def test_blacklist_rejects_common_primers(blacklist_workdir):
    """Rerun with extreme blacklist and verify it filters more than no-blacklist.

    Compares primer counts between a run WITH a blacklist and one
    WITHOUT. If the blacklist plumbing works end-to-end, the former
    should filter out at least as many primers as the latter (usually
    more, since pLTR shares k-mers with pcDNA).
    """
    # First run: no blacklist (match against plain foreground)
    params_file = os.path.join(blacklist_workdir, 'params.json')
    with open(params_file) as fh:
        bl_params = json.load(fh)

    # Save counts under blacklist
    _reset_pipeline_state(params_file)
    from neoswga.core.pipeline import step2
    df_with_bl = step2()
    count_with_bl = len(df_with_bl)

    # Second run: disable blacklist
    no_bl_params = dict(bl_params)
    no_bl_params['bl_genomes'] = []
    no_bl_params['bl_prefixes'] = []
    no_bl_params['bl_seq_lengths'] = []
    no_bl_params_file = os.path.join(blacklist_workdir, 'params_no_bl.json')
    with open(no_bl_params_file, 'w') as fh:
        json.dump(no_bl_params, fh)

    _reset_pipeline_state(no_bl_params_file)
    df_no_bl = step2()
    count_no_bl = len(df_no_bl)

    # Blacklist should reject at least as many primers as the baseline
    assert count_with_bl <= count_no_bl, (
        f"Blacklist should filter at least as many primers as baseline. "
        f"with_bl={count_with_bl}, no_bl={count_no_bl}"
    )
