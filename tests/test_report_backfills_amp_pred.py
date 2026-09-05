"""The report's per-primer quality figure must not be the same for everyone.

`report/metrics.py` builds each `PrimerMetrics` from a row of
`step4_improved_df.csv`, and reads `amp_pred` with

    row.get("amp_pred", row.get("on.target.pred", 0))

That file carries seventeen columns and neither of those. So `amp_pred` came
back 0.0 for every primer in every run -- verified on a real 160-primer panel,
where the set of distinct values was exactly `{0.0}`.

Two user-visible surfaces read it:

  - `executive_summary.py:503` computes
    `quality_score = primer.amp_pred if primer.amp_pred > 0 else 0.5`, so the
    star rating is `int(0.5 * 5 + 0.5)` = 3 for every primer in every report.
    A per-primer column that never varies is worse than no column.
  - `technical_report.py:186` weights `amp_score` at 0.2 of the contribution
    score, permanently contributing zero.

The module already solves this for other columns: `_backfill_primer_tm` and
`_backfill_primer_sites` read from `params.json` and `step2_df.csv` when the
step-4 row does not carry a measurement. `amp_pred` lives in `step3_df.csv` and
was simply never backfilled.

This is NOT an argument that `amp_pred` should influence primer SELECTION. It
should not, and `tests/test_step3_order_is_deterministic.py` records the
evidence. Reporting a number the pipeline computed is a different question from
selecting on it.
"""

import pytest


def _write(tmp_path, name, header, rows):
    (tmp_path / name).write_text(header + "\n" + "\n".join(rows) + "\n")


@pytest.fixture
def results_dir(tmp_path):
    """A step-4 output with no score column, and a step-3 output with one."""
    _write(
        tmp_path,
        "step4_improved_df.csv",
        "primer,set_index,score,coverage,gini,mean_tm,optimizer",
        [
            "AAACCCGGGT,0,1.0,0.9,0.4,32.1,hybrid",
            "ACCCGGGTTT,0,1.0,0.9,0.4,31.8,hybrid",
            "AAGGCCTTAC,0,1.0,0.9,0.4,33.2,hybrid",
        ],
    )
    _write(
        tmp_path,
        "step3_df.csv",
        "primer,gini,on.target.pred",
        [
            "AAACCCGGGT,0.4,19.0",
            "ACCCGGGTTT,0.4,11.0",
            "AAGGCCTTAC,0.4,15.0",
        ],
    )
    return tmp_path


def test_the_step4_file_really_lacks_the_column():
    """Guard the guard, against the shipped output rather than the fixture."""
    import pathlib

    real = pathlib.Path("runs/gc_tiers/mid_ecoli/step4_improved_df.csv")
    if not real.exists():
        pytest.skip("no archived run available")
    header = real.read_text().splitlines()[0]
    assert "amp_pred" not in header and "on.target.pred" not in header


def test_amp_pred_is_backfilled_from_step3(results_dir):
    from neoswga.core.report.metrics import collect_pipeline_metrics

    metrics = collect_pipeline_metrics(str(results_dir))
    values = {p.sequence: p.amp_pred for p in metrics.primers if p.sequence}

    assert values, "no primers were collected"
    assert len(set(values.values())) > 1, (
        f"every primer got the same amp_pred ({values}); the report's quality "
        "column cannot distinguish anything"
    )


def test_the_ordering_of_the_backfilled_values_is_preserved(results_dir):
    """19.0 > 15.0 > 11.0 must survive whatever normalisation is applied."""
    from neoswga.core.report.metrics import collect_pipeline_metrics

    metrics = collect_pipeline_metrics(str(results_dir))
    values = {p.sequence: p.amp_pred for p in metrics.primers if p.sequence}

    assert values["AAACCCGGGT"] > values["AAGGCCTTAC"] > values["ACCCGGGTTT"]


def test_a_missing_step3_file_is_not_an_error(tmp_path):
    """A report over a directory holding only step-4 output must still build."""
    _write(
        tmp_path,
        "step4_improved_df.csv",
        "primer,set_index,score,coverage,gini,mean_tm,optimizer",
        ["AAACCCGGGT,0,1.0,0.9,0.4,32.1,hybrid"],
    )
    from neoswga.core.report.metrics import collect_pipeline_metrics

    metrics = collect_pipeline_metrics(str(tmp_path))
    assert metrics.primers


def test_the_star_rating_varies_between_primers(results_dir):
    """The surface a reader actually sees.

    Every report ever produced showed three stars against every primer.
    """
    from neoswga.core.report.executive_summary import _format_primer_row
    from neoswga.core.report.metrics import collect_pipeline_metrics

    metrics = collect_pipeline_metrics(str(results_dir))
    rendered = [_format_primer_row(i, p) for i, p in enumerate(metrics.primers)]
    star_counts = {row.count("★") for row in rendered}

    assert len(star_counts) > 1, (
        f"every primer rendered with the same star rating ({star_counts}); the "
        "column is decoration rather than information"
    )
