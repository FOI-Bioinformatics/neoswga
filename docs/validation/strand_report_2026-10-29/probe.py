"""What the strand section of the quality report says, and what it drops.

Phase 7.1 of `docs/reports/PLAN_2026-10-28_pool_design_trustworthiness.md`.

Builds a results directory whose `step4_improved_df_summary.json` carries the
full `strand_stats` block the optimizer writes -- five figures for the target
and five for the host -- plus the two headline scalars, one of which is `null`
because the panel has fewer than two sites on that genome and alternation is
then 0/0 rather than zero.

Run from the repository root:

    PYTHONPATH=. python docs/validation/strand_report_2026-10-29/probe.py
"""

from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path

CSV = (
    "primer,set_index,score,normalized_score,coverage,bg_coverage,selectivity,"
    "mean_gap,max_gap,coverage_uniformity,gap_gini,gap_entropy,dimer_risk_score,"
    "strand_alternation,strand_coverage_ratio,mean_tm,optimizer\n"
    "AAAACCCCGGGG,0,1.0,1.0,0.5,0.01,10.0,500,4000,0.4,0.4,2.0,0.1,,0.04,31.0,hybrid\n"
    "TTTTGGGGCCCC,0,1.0,1.0,0.5,0.01,10.0,500,4000,0.4,0.4,2.0,0.1,,0.04,31.0,hybrid\n"
)

PARAMS = {
    "fg_prefixes": ["target"],
    "bg_prefixes": ["host"],
    "fg_seq_lengths": [1_270_000],
    "bg_seq_lengths": [4_000_000],
    "fg_genomes": ["target.fna"],
    "bg_genomes": ["host.fna"],
    "polymerase": "equiphi29",
    "reaction_temp": 43.0,
}

# The target has 25 forward sites and 1 reverse site: a starved strand, which
# is the COATswga failure mode (ranking on the mean of the two strand
# coverages hides it). The host is balanced, which on a HOST is the bad news.
SUMMARY = {
    "primers": ["AAAACCCCGGGG", "TTTTGGGGCCCC"],
    "score": 1.0,
    "status": "success",
    "optimizer_name": "hybrid",
    "num_primers": 2,
    "metrics": {
        "fg_coverage": 0.5,
        "bg_coverage": 0.01,
        "total_fg_sites": 26,
        "total_bg_sites": 40,
        "selectivity_density": 4.2,
        "extension_reach": 3000,
        # Alternation is 0/0 on this panel: absent, not zero.
        "strand_alternation_score": None,
        "strand_coverage_ratio": 0.04,
        "strand_stats": {
            "target": {
                "strand_alternation_gap_mean": 61_000.0,
                "strand_alternation_gap_max": 412_000.0,
                "strand_alternation_score": None,
                "strand_coverage_ratio": 0.04,
                "longest_same_strand_run": 25,
            },
            "host": {
                "strand_alternation_gap_mean": 9_100.0,
                "strand_alternation_gap_max": 48_000.0,
                "strand_alternation_score": 0.71,
                "strand_coverage_ratio": 0.95,
                "longest_same_strand_run": 3,
            },
        },
    },
}


def build(root: Path, summary: dict) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    (root / "params.json").write_text(json.dumps(PARAMS, indent=2), encoding="utf-8")
    (root / "step4_improved_df.csv").write_text(CSV, encoding="utf-8")
    (root / "step4_improved_df_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return root


def section(html: str) -> str:
    """The Strand Balance section of the rendered report, tags stripped."""
    match = re.search(r'<div class="section"><h2>Strand Balance.*?</div>', html, re.S)
    if match is None:
        return "(no Strand Balance section)"
    text = re.sub(r"<[^>]+>", " ", match.group(0))
    return re.sub(r"\s+", " ", text).strip()


def _collect(summary: dict):
    """What `collect_pipeline_metrics` makes of this summary on disk."""
    from neoswga.core.report.metrics import collect_pipeline_metrics

    with tempfile.TemporaryDirectory() as tmp:
        root = build(Path(tmp) / "work", summary)
        return collect_pipeline_metrics(str(root))


def _starved() -> dict:
    """The same panel with the reverse strand EMPTY on the target.

    `strand_coverage_ratio` is min/max over the site counts, so 0.0 with any
    site at all means the minority strand carries none. The gap figures then
    hold the genome length, which is the conservative encoding of "no
    convergent pair anywhere".
    """
    summary = json.loads(json.dumps(SUMMARY))
    target = summary["metrics"]["strand_stats"]["target"]
    target["strand_coverage_ratio"] = 0.0
    target["strand_alternation_gap_mean"] = float(PARAMS["fg_seq_lengths"][0])
    target["strand_alternation_gap_max"] = float(PARAMS["fg_seq_lengths"][0])
    target["longest_same_strand_run"] = 26
    summary["metrics"]["strand_coverage_ratio"] = 0.0
    return summary


def main() -> int:
    from neoswga.core.report.technical_report import _render_strand_balance

    metrics = _collect(SUMMARY)
    uniformity = metrics.uniformity
    reach = metrics.coverage.extension_reach if metrics.coverage else None

    print("== what reached the report object ==")
    print(f"strand_alternation_score : {uniformity.strand_alternation_score!r}")
    print(f"strand_coverage_ratio    : {uniformity.strand_coverage_ratio!r}")
    print(f"strand_stats             : {getattr(uniformity, 'strand_stats', '(no field)')!r}")
    print(f"strand_roles             : {getattr(uniformity, 'strand_roles', '(no field)')!r}")
    print(f"reach passed to renderer : {reach!r}")

    print()
    print("== 25 forward sites to 1 reverse on the target ==")
    print(section(_render_strand_balance(uniformity, reach)))

    starved = _collect(_starved())
    print()
    print("== the reverse strand empty on the target ==")
    print(
        section(
            _render_strand_balance(
                starved.uniformity,
                starved.coverage.extension_reach if starved.coverage else None,
            )
        )
    )

    print()
    print("== what the summary on disk carried ==")
    for prefix, stats in SUMMARY["metrics"]["strand_stats"].items():
        print(f"{prefix}: {stats}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
