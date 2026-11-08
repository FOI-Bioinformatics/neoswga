#!/usr/bin/env python3
"""Do the background-measured panel limits bind differently under `worst-case`?

`docs/validation/2026-11-08-background-aggregate-measured.md` compared
`background_aggregate: sum` against `worst-case` on the D3 design and found
little difference -- with **no panel limits configured**. Its own "Not
measured" list opens with the reason that matters:

> Panel limits under `worst-case` on real genomes. Covered by synthetic tests
> only (`tests/test_per_host_figures.py`).

Step 4's whole difference between the two readings is how
`min_selectivity_density`, `max_background_sites` and `max_host_coverage` are
read -- against the pooled figure under `sum`, against the per-host reduction
under `worst-case`. An unlimited run cannot see it. This script runs the arm
that can.

Two directories, one per reading, each built from `count-kmers` onward (the
step-2 gate differs, so the shortlists differ and one directory cannot serve
both). Per reading: each panel size unlimited, then the same size with one
limit set.

**The thresholds are not recommendations.** They are chosen from the unlimited
arms so that the pooled figure satisfies the limit and the worst-host figure
does not, which is the only configuration in which the two readings can be
told apart. `--limit` takes them explicitly for that reason; no default
threshold is offered, because no reference for these limits exists
(CLAUDE.md, "Panel limits").

    PYTHONPATH=. python docs/validation/aggregate_limits_2026-11-08/measure.py \
        --sum-work <dir> --worst-work <dir> --out results.json \
        --sizes 6,12,24 --limit min_selectivity_density=<v>
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

STEP4 = ("step4_improved_df.csv", "step4_improved_df_summary.json", "design_failure.json")

#: The three limits whose reading `background_aggregate` changes. Any other
#: limit key is rejected: this script exists to compare the pooled figure with
#: the per-host reduction, and a limit that has no per-host reduction would
#: read identically under both and make the comparison look null.
BACKGROUND_LIMITS = (
    "min_selectivity_density",
    "max_background_sites",
    "max_host_coverage",
)


def _clear_step4(work: Path) -> None:
    for name in STEP4:
        (work / name).unlink(missing_ok=True)


def _delivered(work: Path) -> list[str]:
    """Set 0's primers, in file order."""
    path = work / "step4_improved_df.csv"
    if not path.exists():
        return []
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return []
    return [
        str(row["primer"]).strip().upper()
        for row in rows
        if str(row.get("set_index", "0")) in ("0", "0.0")
    ]


def _run(work: Path, size: int, limits: dict[str, float] | None, label: str) -> dict:
    params = json.loads((work / "params.json").read_text())
    params["num_primers"] = size
    params["target_set_size"] = size
    for key in BACKGROUND_LIMITS:
        params.pop(key, None)
    for key, value in (limits or {}).items():
        params[key] = value
    name = f"params_{label}_{size}.json"
    (work / name).write_text(json.dumps(params, indent=2))
    _clear_step4(work)

    started = time.time()
    proc = subprocess.run(
        ["neoswga", "optimize", "-j", name], cwd=work, capture_output=True, text=True
    )
    elapsed = round(time.time() - started, 1)

    record = {
        "reading": params.get("background_aggregate"),
        "size": size,
        "limits": dict(limits or {}),
        "label": label,
        "seconds": elapsed,
        "returncode": proc.returncode,
    }
    summary_path = work / "step4_improved_df_summary.json"
    if not summary_path.exists():
        # A refused design writes no summary. That is an outcome, not a gap:
        # record the reason rather than a missing figure.
        record["failed"] = True
        record["stderr_tail"] = proc.stderr[-1500:]
        failure = work / "design_failure.json"
        if failure.exists():
            record["design_failure"] = json.loads(failure.read_text())
        return record

    summary = json.loads(summary_path.read_text())
    metrics = summary.get("metrics") or {}
    if "fg_coverage" not in metrics:
        raise SystemExit(f"no fg_coverage in the summary for {label} n={size}")
    # The acceptance record is its own file, not a key of the summary:
    # `validation_record.panel_assessment` writes
    # step4_improved_df_validation.json, and the report reads each figure back
    # from there in preference to the summary.
    assessment = {}
    validation_path = work / "step4_improved_df_validation.json"
    if validation_path.exists():
        report = json.loads(validation_path.read_text())
        found = report.get("panel_assessment") or report.get("assessment") or report
        if isinstance(found, dict):
            assessment = found
    record["validation_report_present"] = validation_path.exists()
    record["validator_issues"] = [
        issue.get("code")
        for issue in (json.loads(validation_path.read_text()).get("issues") or [])
        if isinstance(issue, dict)
    ] if validation_path.exists() else None
    record.update(
        {
            "primers": _delivered(work),
            "delivered": len(_delivered(work)),
            "fg_coverage": metrics.get("fg_coverage"),
            "worst_target_coverage": metrics.get("worst_target_coverage"),
            # Pooled, under both readings, by contract.
            "selectivity_density": metrics.get("selectivity_density"),
            "total_bg_sites": metrics.get("total_bg_sites"),
            "bg_coverage": metrics.get("bg_coverage"),
            # The reductions over hosts.
            "worst_host_selectivity_density": metrics.get("worst_host_selectivity_density"),
            "worst_host_background_sites": metrics.get("worst_host_background_sites"),
            "worst_host_bg_coverage": metrics.get("worst_host_bg_coverage"),
            "per_host": metrics.get("per_host"),
            "background_aggregate_recorded": metrics.get("background_aggregate")
            or summary.get("background_aggregate"),
            "violations": assessment.get("violations"),
            "blocking_violations": assessment.get("blocking_violations"),
            "advisory_violations": assessment.get("advisory_violations"),
            "qualified": assessment.get("qualified"),
            "panel_regime": summary.get("panel_regime"),
            "assessment_metrics": assessment.get("metrics"),
            "message": summary.get("message"),
            "stdout_tail": proc.stdout[-1200:],
        }
    )
    return record


def _parse_limits(values: list[str]) -> dict[str, float]:
    limits: dict[str, float] = {}
    for item in values:
        key, _, raw = item.partition("=")
        key = key.strip()
        if key not in BACKGROUND_LIMITS:
            raise SystemExit(f"{key} is not one of {BACKGROUND_LIMITS}")
        limits[key] = float(raw)
    return limits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sum-work", required=True)
    ap.add_argument("--worst-work", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sizes", default="6,12,24")
    ap.add_argument(
        "--limit",
        action="append",
        default=[],
        help="key=value, repeatable; only the three background-measured limits",
    )
    ap.add_argument(
        "--unlimited-only",
        action="store_true",
        help="Run only the unlimited arms, to read the figures a threshold is chosen from",
    )
    args = ap.parse_args()

    sizes = [int(s) for s in args.sizes.split(",") if s.strip()]
    limits = _parse_limits(args.limit)
    if not limits and not args.unlimited_only:
        raise SystemExit("pass --limit, or --unlimited-only to read the figures first")

    works = {"sum": Path(args.sum_work).resolve(), "worst-case": Path(args.worst_work).resolve()}
    runs = []
    for size in sizes:
        for reading, work in works.items():
            for label, applied in (("unlimited", None), ("limited", limits)):
                if label == "limited" and args.unlimited_only:
                    continue
                record = _run(work, size, applied, f"{reading}_{label}")
                runs.append(record)
                print(
                    f"  {reading:>10} {label:>9} n={size:<3} "
                    f"delivered={record.get('delivered')} "
                    f"fg={record.get('fg_coverage')} "
                    f"pooled_density={record.get('selectivity_density')} "
                    f"worst_host_density={record.get('worst_host_selectivity_density')} "
                    f"blocking={record.get('blocking_violations')} "
                    f"({record['seconds']}s)",
                    flush=True,
                )

    Path(args.out).write_text(
        json.dumps(
            {
                "design": {
                    "target": "M. tuberculosis H37Rv, 4,411,532 bp, GC 65.6%, 1 record",
                    "hosts": [
                        "S. aureus NCTC 8325, 2,821,361 bp",
                        "human chr1 (hg38), 248,956,422 bp",
                    ],
                    "host_length_ratio": 88.2,
                    "k": 12,
                    "polymerase": "phi29",
                    "reaction_temp": 30.0,
                    "max_primer": 2000,
                },
                "limits": limits,
                "sizes": sizes,
                "runs": runs,
            },
            indent=2,
        )
    )
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
