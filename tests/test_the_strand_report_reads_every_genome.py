"""Ten strand figures reached disk and two reached the report, one of them wrong.

Phase 7.1 of `docs/reports/PLAN_2026-10-28_pool_design_trustworthiness.md`;
measured in `docs/validation/2026-10-29-the-strand-report.md`.

`base_optimizer._compute_metrics` writes
`metrics.strand_stats` into `step4_improved_df_summary.json`: five figures for
every foreground genome and five for the host. The quality report rendered the
two headline scalars of the FIRST target and nothing else, so the host's strand
structure and three of five figures per genome were computed, serialised, and
read by nobody -- Known Issue 8's shape inside the report.

Worse, one of the two it did render was wrong in the favourable direction.
Alternation is the fraction of ADJACENT site pairs on opposite strands, so a
panel with fewer than two sites on a genome is 0/0 and
`compute_strand_alternation_stats` writes `null` for it on purpose
(`tests/test_strand_scores_say_when_they_are_unmeasurable.py` pins that). The
report converted it with a literal 0.0 default and rendered "0.00": "could not
be measured" shown as "measured: no alternation at all". Known Issues 5, 6
and 13.

Nothing here introduces a threshold. Known Issue 18 keeps the strand scalars
unconstrainable and no cut on them has a reference, so the section flags one
FACT (a genome with sites on a single strand) and one MECHANISM (the widest
stretch with no convergent pair, in units of the polymerase's reach). Both are
reports, not limits.

The gap this closes is the one
`docs/validation/2026-09-22-swga-comparative-review.md` records against this
side: COATswga ranks on the mean of its two strand coverages, so 0.99/0.30
scores 0.645 and beats a balanced 0.60/0.60 with nothing flagging the
imbalance, and NeoSWGA had the data and did not report it either.
"""

import json
import re

import pytest

from neoswga.core.report.metrics import collect_pipeline_metrics
from neoswga.core.report.technical_report import _render_strand_balance

CSV = (
    "primer,set_index,score,normalized_score,coverage,bg_coverage,selectivity,"
    "mean_gap,max_gap,coverage_uniformity,gap_gini,gap_entropy,dimer_risk_score,"
    "strand_alternation,strand_coverage_ratio,mean_tm,optimizer\n"
    "AAAACCCCGGGG,0,1.0,1.0,0.5,0.01,10.0,500,4000,0.4,0.4,2.0,0.1,,0.04,31.0,hybrid\n"
    "TTTTGGGGCCCC,0,1.0,1.0,0.5,0.01,10.0,500,4000,0.4,0.4,2.0,0.1,,0.04,31.0,hybrid\n"
)

TARGET_STATS = {
    "strand_alternation_gap_mean": 61_000.0,
    "strand_alternation_gap_max": 412_000.0,
    # 0/0 on this panel: absent, not zero.
    "strand_alternation_score": None,
    "strand_coverage_ratio": 0.04,
    "longest_same_strand_run": 25,
}

HOST_STATS = {
    "strand_alternation_gap_mean": 9_100.0,
    "strand_alternation_gap_max": 48_000.0,
    "strand_alternation_score": 0.71,
    "strand_coverage_ratio": 0.95,
    "longest_same_strand_run": 3,
}


def _params(**overrides):
    params = {
        "fg_prefixes": ["target"],
        "bg_prefixes": ["host"],
        "fg_seq_lengths": [1_270_000],
        "bg_seq_lengths": [4_000_000],
        "fg_genomes": ["target.fna"],
        "bg_genomes": ["host.fna"],
    }
    params.update(overrides)
    return params


def _summary(strand_stats=None, *, reach=3000, **metric_overrides):
    metrics = {
        "fg_coverage": 0.5,
        "bg_coverage": 0.01,
        "total_fg_sites": 26,
        "total_bg_sites": 40,
        "selectivity_density": 4.2,
        "strand_alternation_score": None,
        "strand_coverage_ratio": 0.04,
    }
    if reach is not None:
        metrics["extension_reach"] = reach
    if strand_stats is not None:
        metrics["strand_stats"] = strand_stats
    metrics.update(metric_overrides)
    return {
        "primers": ["AAAACCCCGGGG", "TTTTGGGGCCCC"],
        "score": 1.0,
        "status": "success",
        "optimizer_name": "hybrid",
        "num_primers": 2,
        "metrics": metrics,
    }


def _rendered(tmp_path, summary, params=None):
    """The strand section a report would carry for this run, tags stripped."""
    root = tmp_path / "work"
    root.mkdir()
    (root / "params.json").write_text(json.dumps(params or _params()), encoding="utf-8")
    (root / "step4_improved_df.csv").write_text(CSV, encoding="utf-8")
    (root / "step4_improved_df_summary.json").write_text(json.dumps(summary), encoding="utf-8")

    metrics = collect_pipeline_metrics(str(root))
    html = _render_strand_balance(
        metrics.uniformity,
        metrics.coverage.extension_reach if metrics.coverage else None,
    )
    return metrics, html, re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html)).strip()


def _both_genomes():
    return {"target": dict(TARGET_STATS), "host": dict(HOST_STATS)}


class TestAnUnmeasurableScoreSaysSo:
    def test_a_null_alternation_score_does_not_become_zero(self, tmp_path):
        metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        assert metrics.uniformity.strand_alternation_score is None
        assert "not measurable" in text

    def test_the_alternation_row_carries_no_number(self, tmp_path):
        """The defect rendered `Strand alternation score 0.00`."""
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        row = text.split("Strand alternation score", 1)[1]
        assert not row.lstrip().startswith("0.00")

    def test_a_measured_score_still_renders_as_a_number(self, tmp_path):
        """A real zero must stay distinguishable from an absent one."""
        _metrics, _html, text = _rendered(
            tmp_path,
            _summary(_both_genomes(), strand_alternation_score=0.0),
        )

        row = text.split("Strand alternation score", 1)[1]
        assert row.lstrip().startswith("0.00")


class TestEveryGenomeAndEveryFigure:
    def test_the_host_reaches_the_report(self, tmp_path):
        metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        assert set(metrics.uniformity.strand_stats) == {"target", "host"}
        assert "host" in text

    @pytest.mark.parametrize(
        "figure",
        ["48,000", "9,100", "0.95", "0.71", "412,000", "61,000", "0.04", "25"],
    )
    def test_each_of_the_ten_figures_is_rendered(self, tmp_path, figure):
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        assert figure in text

    def test_a_prefix_nobody_measured_gets_no_row(self, tmp_path):
        """Absent is not zero: an unmeasured genome must not render as one."""
        params = _params(fg_prefixes=["target", "second_target"])
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()), params)

        assert "second_target" not in text

    def test_a_run_with_no_strand_stats_still_renders_the_headline(self, tmp_path):
        """Directories written before `strand_stats` existed are supported."""
        _metrics, _html, text = _rendered(tmp_path, _summary(None))

        assert "Strand alternation score" in text
        assert "Genome" not in text


class TestRolesComeFromParams:
    def test_the_target_and_the_host_are_labelled(self, tmp_path):
        metrics, _html, _text = _rendered(tmp_path, _summary(_both_genomes()))

        assert metrics.uniformity.strand_roles == {"target": "target", "host": "host"}

    def test_a_prefix_params_does_not_place_is_unknown_not_guessed(self, tmp_path):
        stats = {"orphan": dict(TARGET_STATS)}
        _metrics, _html, text = _rendered(tmp_path, _summary(stats))

        assert "unknown" in text
        assert "orphan" in text


class TestWhatIsFlagged:
    def _starved(self, role_prefix="target"):
        stats = _both_genomes()
        stats[role_prefix]["strand_coverage_ratio"] = 0.0
        # No opposite-strand pair anywhere: the gaps hold the genome length.
        stats[role_prefix]["strand_alternation_gap_max"] = 1_270_000.0
        stats[role_prefix]["strand_alternation_gap_mean"] = 1_270_000.0
        return stats

    def test_a_starved_target_strand_is_flagged(self, tmp_path):
        _metrics, _html, text = _rendered(tmp_path, _summary(self._starved("target")))

        assert "target (target) has binding sites on one strand only" in text
        assert "This is the imbalance a mean of the two strand coverages hides" in text

    def test_a_starved_host_strand_reads_the_other_way(self, tmp_path):
        _metrics, _html, text = _rendered(tmp_path, _summary(self._starved("host")))

        # Not the section prose, which names both readings: the host's own note.
        assert "On a host that is the favourable direction" in text
        assert "a mean of the two strand coverages hides" not in text

    def test_a_balanced_panel_is_not_flagged_for_imbalance(self, tmp_path):
        stats = {"target": dict(HOST_STATS)}
        _metrics, _html, text = _rendered(tmp_path, _summary(stats))

        assert "one strand only" not in text

    def test_the_starved_note_replaces_the_gap_note(self, tmp_path):
        """The gap IS the genome length there; printing it as a multiple of the
        reach would dress a conservative encoding up as a measurement."""
        _metrics, _html, text = _rendered(tmp_path, _summary(self._starved("target")))

        assert "423" not in text

    def test_a_gap_wider_than_the_reach_is_reported_against_the_reach(self, tmp_path):
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        assert "137.3x the 3,000 bp reach" in text
        assert "No reference" in text

    def test_a_gap_inside_the_reach_is_not_flagged(self, tmp_path):
        stats = {"target": dict(TARGET_STATS, strand_alternation_gap_max=1_500.0)}
        _metrics, _html, text = _rendered(tmp_path, _summary(stats))

        assert "the 3,000 bp reach" not in text

    def test_no_reach_means_no_multiple_rather_than_a_fabricated_one(self, tmp_path):
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes(), reach=None))

        assert "reach" not in text.split("Genome", 1)[-1]


class TestTheSectionDoesNotOverclaim:
    def test_the_balance_is_named_a_site_count_ratio(self, tmp_path):
        """It is min/max over site COUNTS. COATswga's quantity, which the
        comparative review contrasts it with, is per-strand COVERAGE."""
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        assert "COUNTS" in text
        assert "not per-strand coverage" in text

    def test_the_figures_are_marked_unconstrainable(self, tmp_path):
        _metrics, _html, text = _rendered(tmp_path, _summary(_both_genomes()))

        assert "Known Issue 18" in text
        assert "none has a published threshold" in text

    def test_nothing_renders_when_the_optimizer_measured_nothing(self):
        from neoswga.core.report.metrics import UniformityMetrics

        assert _render_strand_balance(None, 3000) == ""
        assert _render_strand_balance(UniformityMetrics(), 3000) == ""
