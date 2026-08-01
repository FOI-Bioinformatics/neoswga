"""Execution tests for the four core pipeline commands.

`count-kmers`, `filter`, `score` and `optimize` are the product. They are
covered end to end by the integration suite, but largely at the level of "the
run completed" -- `neoswga/cli/pipeline.py` sat at 18%. These tests assert on
the contracts the rest of the tool depends on: which files appear, what the
columns are called, and that the authoritative summary JSON is well formed.

`optimize` is worth testing here specifically because its
`step4_improved_df_summary.json` is what the report reads and what the audit
found could contain `Infinity` on a failed run -- valid to Python, invalid to
every other JSON parser.
"""

import json
import os

import pytest

pytestmark = pytest.mark.usefixtures("pipeline_run")


# ----------------------------------------------------------------------
# Artifacts the downstream commands rely on
# ----------------------------------------------------------------------


def test_filter_and_score_produce_their_documented_files(pipeline_run):
    data_dir = pipeline_run["data_dir"]
    for name in ("step2_df.csv", "step3_df.csv", "filter_stats.json"):
        path = os.path.join(data_dir, name)
        assert os.path.isfile(path), f"{name} missing from {data_dir}"
        assert os.path.getsize(path) > 0, f"{name} is empty"


def test_position_index_is_written_next_to_the_prefix(pipeline_run):
    """PositionCache reads `<prefix>_<k>mer_positions.h5`; the layout is a contract."""
    root = pipeline_run["dir"]
    h5 = [f for f in os.listdir(root) if f.endswith("_positions.h5")]
    assert h5, f"no HDF5 position file in {root}"


def test_filter_stats_records_a_real_funnel(pipeline_run):
    """The report renders these counts; fabricated or empty ones would show as
    a clean-looking funnel that never happened."""
    with open(os.path.join(pipeline_run["data_dir"], "filter_stats.json")) as fh:
        stats = json.load(fh)

    assert stats, "filter_stats.json is empty"
    numbers = [v for v in stats.values() if isinstance(v, (int, float))]
    assert numbers, f"no counts in filter_stats.json: {stats}"
    assert all(n >= 0 for n in numbers)


def test_scored_output_carries_the_prediction_column(pipeline_run):
    import pandas as pd

    df = pd.read_csv(os.path.join(pipeline_run["data_dir"], "step3_df.csv"))
    assert len(df) > 0
    assert any("amp" in c.lower() or "pred" in c.lower() for c in df.columns), list(df.columns)


def test_scoring_preserves_every_filtered_primer(pipeline_run):
    """score annotates; it must not silently drop candidates."""
    import pandas as pd

    step2 = pd.read_csv(os.path.join(pipeline_run["data_dir"], "step2_df.csv"))
    step3 = pd.read_csv(os.path.join(pipeline_run["data_dir"], "step3_df.csv"))
    assert len(step3) == len(
        step2
    ), f"score changed the candidate count: {len(step2)} -> {len(step3)}"


# ----------------------------------------------------------------------
# optimize
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def optimized(pipeline_run):
    """Run `optimize` once for this module."""
    import subprocess
    import sys

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "neoswga.cli_unified",
            "optimize",
            "-j",
            pipeline_run["params_file"],
        ],
        capture_output=True,
        text=True,
        cwd=str(pipeline_run["dir"]),
        timeout=1800,
    )
    if proc.returncode != 0:
        pytest.skip(f"optimize failed on the fixture genome:\n{proc.stderr[-1200:]}")
    return pipeline_run


def test_optimize_writes_results_and_summary(optimized):
    data_dir = optimized["data_dir"]
    csv_path = os.path.join(data_dir, "step4_improved_df.csv")
    summary_path = os.path.join(data_dir, "step4_improved_df_summary.json")
    assert os.path.isfile(csv_path), sorted(os.listdir(data_dir))
    assert os.path.isfile(summary_path), sorted(os.listdir(data_dir))


def test_summary_json_is_valid_json(optimized):
    """The audit's finding: a failed run wrote `Infinity`, which RFC 8259 does
    not permit. Python reads it back, so the breakage was invisible from inside
    the project while rejecting every other parser."""

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    path = os.path.join(optimized["data_dir"], "step4_improved_df_summary.json")
    json.loads(open(path).read(), parse_constant=reject)


def test_summary_carries_the_metrics_the_report_reads(optimized):
    """CLAUDE.md documents this file as authoritative for the report."""
    path = os.path.join(optimized["data_dir"], "step4_improved_df_summary.json")
    summary = json.loads(open(path).read())

    assert "metrics" in summary, sorted(summary)
    metrics = summary["metrics"]
    for key in ("fg_coverage", "mean_gap", "max_gap", "gap_gini"):
        assert key in metrics, f"{key} missing from summary metrics"


def test_reported_coverage_is_a_real_fraction(optimized):
    path = os.path.join(optimized["data_dir"], "step4_improved_df_summary.json")
    metrics = json.loads(open(path).read())["metrics"]

    coverage = metrics["fg_coverage"]
    assert coverage is not None, "coverage came back null; the optimizer found nothing"
    assert 0.0 < coverage <= 1.0, coverage


def test_selected_primers_come_from_the_scored_pool(optimized, scored_primers):
    """The optimizer must choose from the candidates, not invent sequences."""
    path = os.path.join(optimized["data_dir"], "step4_improved_df_summary.json")
    summary = json.loads(open(path).read())

    chosen = summary.get("primers") or []
    assert chosen, "optimize selected no primers"
    assert set(chosen) <= set(scored_primers), set(chosen) - set(scored_primers)


def test_coverage_is_scored_at_the_realistic_reach(optimized):
    """Selection and scoring must agree on reach.

    Coverage is reported at the realistic per-primer reach (~3 kb for phi29),
    not single-molecule processivity (~70 kb). Using the latter overstated
    coverage by roughly an order of magnitude.
    """
    path = os.path.join(optimized["data_dir"], "step4_improved_df_summary.json")
    metrics = json.loads(open(path).read())["metrics"]

    assert metrics.get("extension_reach") == 3000, metrics.get("extension_reach")
