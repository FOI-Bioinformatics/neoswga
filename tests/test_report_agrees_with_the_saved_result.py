"""The rendered report and the saved result describe the same panel.

Task 5 of the 2026-09-21 valid-design plan: "require saved JSON and rendered
report to agree for the exact exported panel".

The report computes its own estimates from the results CSV and then overrides
them with the optimizer summary where the summary has an opinion. That is the
right order, and it leaves two ways for the two to disagree.

A key the summary does not carry keeps the report's own estimate, which is a
different quantity computed from different inputs. `from_optimizer` exists to
say which the reader is looking at, so what matters is that the flag is honest.

And a key the summary does not carry can take a FAVOURABLE literal default.
`mean_gap`, `max_gap`, `gap_gini` and `gap_entropy` are read with
`.get(key, 0.0)`, and zero is the best possible value for each: a `max_gap` of
0.0 says the panel has no coverage hole at all. An absent measurement must not
read as the best one.
"""

import json

import pytest

from neoswga.core.report.metrics import collect_pipeline_metrics

PRIMERS = ["AAAACCCCGGGG", "TTTTGGGGCCCC", "ACGTACGTACGT"]

SUMMARY_METRICS = {
    "fg_coverage": 0.6211,
    "bg_coverage": 0.0402,
    "total_fg_sites": 240,
    "total_bg_sites": 31,
    "selectivity_ratio": 7.74,
    "selectivity_density": 21.3,
    "mean_gap": 1800.0,
    "max_gap": 9100.0,
    "gap_gini": 0.41,
    "gap_entropy": 3.2,
    "coverage_uniformity": 0.37,
    "extension_reach": 3000,
}


def build_results_dir(tmp_path, metrics=None, drop=()):
    """A results directory as the pipeline leaves one."""
    import csv

    directory = tmp_path / "results"
    directory.mkdir()

    with open(directory / "step4_improved_df.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["primer", "set_index", "fg_count", "bg_count"])
        for primer in PRIMERS:
            writer.writerow([primer, 0, 80, 10])

    payload = dict(SUMMARY_METRICS if metrics is None else metrics)
    for key in drop:
        payload.pop(key, None)
    (directory / "step4_improved_df_summary.json").write_text(
        json.dumps({"primers": PRIMERS, "metrics": payload})
    )

    (directory / "params.json").write_text(
        json.dumps(
            {
                "fg_genomes": ["target.fna"],
                "fg_prefixes": ["target"],
                "fg_seq_lengths": [100000],
                "bg_prefixes": [],
                "polymerase": "phi29",
                "reaction_temp": 30.0,
                "min_k": 12,
                "max_k": 12,
            }
        )
    )
    return directory


# ---------------------------------------------------------------------------
# Agreement
# ---------------------------------------------------------------------------


def test_the_report_shows_the_panel_the_result_holds(tmp_path):
    metrics = collect_pipeline_metrics(str(build_results_dir(tmp_path)))

    assert [primer.sequence for primer in metrics.primers] == PRIMERS


@pytest.mark.parametrize(
    "summary_key,report_path",
    [
        ("fg_coverage", "coverage.overall_coverage"),
        ("mean_gap", "coverage.mean_gap"),
        ("max_gap", "coverage.max_gap"),
        ("gap_gini", "coverage.gap_gini"),
        ("gap_entropy", "coverage.gap_entropy"),
        ("selectivity_ratio", "specificity.selectivity_ratio"),
        ("selectivity_density", "specificity.enrichment_ratio"),
        ("bg_coverage", "specificity.bg_coverage"),
        ("total_fg_sites", "specificity.target_sites"),
        ("total_bg_sites", "specificity.background_sites"),
    ],
)
def test_every_rendered_quantity_equals_the_saved_one(tmp_path, summary_key, report_path):
    """One number, two readers. A rendered figure and a stored one must match."""
    metrics = collect_pipeline_metrics(str(build_results_dir(tmp_path)))

    section, attribute = report_path.split(".")
    rendered = getattr(getattr(metrics, section), attribute)

    assert rendered == pytest.approx(SUMMARY_METRICS[summary_key]), (
        f"{report_path} renders {rendered} for a saved {summary_key} of "
        f"{SUMMARY_METRICS[summary_key]}"
    )


def test_the_report_says_when_a_figure_came_from_the_optimizer(tmp_path):
    """`from_optimizer` is how a reader tells a measurement from an estimate."""
    metrics = collect_pipeline_metrics(str(build_results_dir(tmp_path)))

    assert metrics.coverage.from_optimizer is True
    assert metrics.specificity.from_optimizer is True


def test_the_reach_travels_with_the_coverage_figure(tmp_path):
    """Coverage is uninterpretable without it: 41.3 percent at 1 kb, 93.5 at 5."""
    metrics = collect_pipeline_metrics(str(build_results_dir(tmp_path)))

    assert metrics.coverage.extension_reach == 3000


# ---------------------------------------------------------------------------
# An absent measurement must not read as the best one
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("key", ["max_gap", "mean_gap", "gap_gini", "gap_entropy"])
def test_a_missing_gap_measurement_does_not_render_as_the_best_value(tmp_path, key):
    """Zero is the FAVOURABLE value for each of these.

    A `max_gap` of 0.0 says the panel leaves no hole anywhere, which is the
    best result available and the one a reader would most like to see. Reaching
    it because the summary did not carry the key is the silent-zero pattern
    this project has recorded four separate routes into.

    A summary written by an older version, or by a path that did not compute
    gaps, is exactly the case. The remedy is that absent renders as absent.
    """
    directory = build_results_dir(tmp_path, drop=[key])
    metrics = collect_pipeline_metrics(str(directory))

    rendered = getattr(metrics.coverage, key)

    assert rendered != 0.0, (
        f"{key} is absent from the summary and renders as 0.0, which is the "
        "best possible value for it"
    )


def test_a_summary_with_no_metrics_block_leaves_the_estimate_flagged(tmp_path):
    """No measurements at all is a legitimate state, and must be visible.

    The report still has its own estimates from the CSV. What it must not do is
    present them as though the optimizer had measured them.
    """
    directory = build_results_dir(tmp_path)
    (directory / "step4_improved_df_summary.json").write_text(json.dumps({"primers": PRIMERS}))

    metrics = collect_pipeline_metrics(str(directory))

    assert metrics.coverage.from_optimizer is False
    assert metrics.specificity.from_optimizer is False


def test_a_null_measurement_stays_null_rather_than_becoming_zero(tmp_path):
    """`json_safe` writes null for a non-finite quantity, deliberately.

    A failed optimization carries `max_gap` of infinity, which is honest: no
    coverage really is unboundedly bad. JSON cannot carry it, so null is
    written. Null must not then be read back as zero, which would turn the
    worst possible result into the best.
    """
    payload = dict(SUMMARY_METRICS)
    payload["max_gap"] = None
    directory = build_results_dir(tmp_path, metrics=payload)

    metrics = collect_pipeline_metrics(str(directory))

    assert metrics.coverage.max_gap != 0.0
