"""Phase 17E — plasmid golden-snapshot regression test.

Nightly CI previously only checked that `step4_improved_df.csv` exists
after running the full pipeline on examples/plasmid_example/. A silent
scoring or coverage drift would still leave the file in place, so
regressions could slip in unnoticed. This test pins the expected
outputs for seed=42 on the plasmid example so any behaviour change in
the pipeline — optimizer swap, coverage math error, scoring model
drift — trips the snapshot.

The test is marked @pytest.mark.integration so it runs in nightly but
stays out of the fast-feedback unit test loop.
"""

import json
import os
import shutil
import tempfile
from pathlib import Path

import pandas as pd
import pytest


EXAMPLE_DIR = os.path.join(
    os.path.dirname(__file__), '..', '..', 'examples', 'plasmid_example'
)


def _reset_pipeline_state(params_file):
    """Reset global pipeline state and set json_file for a fresh run."""
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
def plasmid_workdir():
    """Copy the plasmid example into a tmpdir so the pipeline runs in
    isolation and cleanup is automatic."""
    if not os.path.isdir(EXAMPLE_DIR):
        pytest.skip("Plasmid example data not available")

    tmpdir = tempfile.mkdtemp(prefix='neoswga_golden_')
    for fname in os.listdir(EXAMPLE_DIR):
        src = os.path.join(EXAMPLE_DIR, fname)
        if os.path.isfile(src):
            shutil.copy2(src, tmpdir)

    params_path = os.path.join(tmpdir, 'params.json')
    with open(params_path) as fh:
        params_data = json.load(fh)
    # The example params.json omits min_amp_pred; the small plasmid needs
    # it permissive so the RF score doesn't truncate to zero.
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
def test_plasmid_golden_snapshot(plasmid_workdir):
    """Full pipeline (filter -> score -> optimize) on the plasmid example
    with seed=42 must produce a primer set whose key metrics sit inside
    the locked-in tolerance bounds.

    The assertions are on metrics rather than exact primer sequences
    because the tiling / dominating-set optimizers can pick permutation-
    equivalent sets that differ by a single primer when coverage ties.
    Metric values are stable across optimizers that satisfy the same
    objective, so asserting on them catches real drift without making
    the test brittle to a legitimate tie-break change.

    Failing modes this catches:
    - Optimizer default silently swapped (would change the chosen set size
      and coverage pattern).
    - Coverage math error (would push fg_coverage away from 1.0).
    - Scoring model drift (would change the primer ranking that feeds
      into optimization).
    - Extension-reach default regressed to processivity (would inflate
      coverage on a 6 kb plasmid).
    """
    params_file = os.path.join(plasmid_workdir, 'params.json')
    for fname in ['step2_df.csv', 'step3_df.csv',
                  'step4_improved_df.csv',
                  'step4_improved_df_summary.json']:
        path = os.path.join(plasmid_workdir, fname)
        if os.path.exists(path):
            os.remove(path)

    _reset_pipeline_state(params_file)

    from neoswga.core.pipeline import step2, step3
    from neoswga.core.unified_optimizer import optimize_step4

    df2 = step2()
    assert len(df2) > 0, "Step 2 produced no candidates — filter regression"

    _reset_pipeline_state(params_file)
    df3 = step3()
    assert len(df3) > 0, "Step 3 produced no scored primers"

    _reset_pipeline_state(params_file)
    primer_sets, scores, _cache = optimize_step4(
        optimization_method='dominating-set',
        verbose=False,
        seed=42,
    )
    assert len(primer_sets) > 0, "Step 4 returned no primer sets"
    best_set = primer_sets[0]
    assert len(best_set) > 0, "Best primer set is empty"

    # Snapshot 1: primer count within sensible bounds for this plasmid.
    # Dominating-set with phi29 3 kb per-primer reach can legitimately
    # pick a single primer for a 6.2 kb plasmid — a 6 kb window covers
    # essentially the whole circular target. Tiling / hybrid on the same
    # data historically pick 2-3 primers. Accept 1-6; lower bound 1
    # guards against "optimizer returned empty set with status=SUCCESS".
    assert 1 <= len(best_set) <= 6, (
        f"Expected 1-6 primers for plasmid; got {len(best_set)}: {best_set}"
    )

    # Snapshot 2: read back summary JSON and check metrics.
    summary_path = os.path.join(plasmid_workdir, 'step4_improved_df_summary.json')
    assert os.path.exists(summary_path), (
        "step4_improved_df_summary.json missing — optimizer did not "
        "write the summary"
    )
    with open(summary_path) as fh:
        summary = json.load(fh)

    metrics = summary.get('metrics', {})
    fg_cov = metrics.get('fg_coverage', 0.0)
    selectivity = metrics.get('selectivity_ratio', 0.0)
    bg_cov = metrics.get('bg_coverage', 0.0)
    dimer_risk = metrics.get('dimer_risk_score', 0.0)

    # Snapshot 3: fg_coverage must be high. The plasmid is 6.2 kb and
    # phi29 reach is 3 kb, so a 2-3 primer set trivially covers it.
    # Saturation-bounded, but we still assert the coverage number is
    # non-degenerate (the Phase 17B warning fires elsewhere, not here).
    assert fg_cov >= 0.95, (
        f"Plasmid fg_coverage should be >=0.95, got {fg_cov}. "
        f"Regression in coverage computation or optimizer dispatch."
    )

    # Snapshot 4: selectivity. For the pcDNA vs pLTR design, selected
    # primers are fg-unique (bg_count=0 for the chosen set historically),
    # so selectivity should be > 1 and bg_coverage low.
    assert selectivity >= 1.0, (
        f"Selectivity should be >=1.0 on pcDNA vs pLTR; got {selectivity}. "
        f"Regression in scoring / filter thresholds."
    )
    assert bg_cov <= 0.20, (
        f"bg_coverage should stay low on pcDNA vs pLTR; got {bg_cov}. "
        f"Possibly a regression in background filter or selectivity weight."
    )

    # Snapshot 5: dimer risk stays low on such a tiny set.
    assert dimer_risk <= 0.10, (
        f"dimer_risk_score should be small on 2-3 primer plasmid set; "
        f"got {dimer_risk}. Regression in dimer scoring."
    )

    # Snapshot 6: step4_improved_df.csv schema is stable.
    step4_csv = os.path.join(plasmid_workdir, 'step4_improved_df.csv')
    df4 = pd.read_csv(step4_csv)
    for col in ('primer', 'score', 'coverage'):
        assert col in df4.columns, (
            f"step4_improved_df.csv missing column '{col}'; got {list(df4.columns)}. "
            f"Schema regression."
        )


@pytest.mark.integration
@pytest.mark.slow
def test_plasmid_snapshot_is_reproducible(plasmid_workdir):
    """Running the same seed twice must produce the same primer set.
    Guards against non-threaded RNG that would cause drift between
    nightly runs even without code changes."""
    params_file = os.path.join(plasmid_workdir, 'params.json')

    def run_once():
        for fname in ['step2_df.csv', 'step3_df.csv',
                      'step4_improved_df.csv',
                      'step4_improved_df_summary.json']:
            path = os.path.join(plasmid_workdir, fname)
            if os.path.exists(path):
                os.remove(path)
        _reset_pipeline_state(params_file)
        from neoswga.core.pipeline import step2, step3
        from neoswga.core.unified_optimizer import optimize_step4
        step2()
        _reset_pipeline_state(params_file)
        step3()
        _reset_pipeline_state(params_file)
        primer_sets, _scores, _cache = optimize_step4(
            optimization_method='dominating-set',
            verbose=False,
            seed=42,
        )
        return tuple(sorted(primer_sets[0]))

    first = run_once()
    second = run_once()
    assert first == second, (
        f"Seed=42 produced different primer sets across two runs:\n"
        f"  first:  {first}\n  second: {second}\n"
        f"RNG is not properly threaded through all optimizers."
    )
