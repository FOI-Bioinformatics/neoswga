#!/usr/bin/env python3
"""Does the delivered panel move under the directional geometry, on a second pool?

`docs/validation/2026-09-23-does-the-panel-move.md` measured this on ONE pool
-- Prevotella against human chr21 -- and recorded two limitations this script
exists to address: one pool, and a background axis that was never exercised
because that candidate pool has zero exact matches anywhere in chr21 by
construction, so `bg_coverage` read 0.00000000 in every run.

This runs the same comparison on a second real pair with a background that does
answer: *E. coli* K-12 MG1655 (4,639,675 bp) against human chr21
(46,709,983 bp), 12-mers, equiphi29 at 43 C.

One `neoswga optimize` per (geometry, reach, panel size), against one prepared
`step3_df.csv`. Step 4 artifacts are removed between arms: a step4 file left by
a previous arm is a confound measured earlier in this project's history.

    PYTHONPATH=. python docs/validation/geometry_second_pool_2026-11-02/measure.py \
        --work <dir with params.json and a prepared pool> --out results.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

# The arms. The pairing follows the first pool's: both geometries at the
# configuration's own reach, then directional at the 1.52x refitted reach from
# `2026-09-23-reach-refit-directional.md`. 4000 is equiphi29's realistic reach.
# The fourth arm is not in the first pool's design. symmetric@6080 separates
# the two things the refitted comparison changes at once: without it, a panel
# that moves between symmetric@4000 and directional@6080 cannot be attributed
# to the geometry rather than to the wider window.
ARMS = (
    ("symmetric", 4000),
    ("directional", 4000),
    ("directional", 6080),
    ("symmetric", 6080),
)
SIZES = (4, 8, 16, 24)

STEP4 = ("step4_improved_df.csv", "step4_improved_df_summary.json", "design_failure.json")


def _clear_step4(work: Path) -> None:
    for name in STEP4:
        (work / name).unlink(missing_ok=True)


def _write_params(work: Path, geometry: str, reach: int, size: int) -> Path:
    params = json.loads((work / "params.json").read_text())
    params["coverage_geometry"] = geometry
    params["coverage_reach"] = reach
    params["num_primers"] = size
    params["target_set_size"] = size
    out = work / f"params_{geometry}_{reach}_{size}.json"
    out.write_text(json.dumps(params, indent=2))
    return out


def _delivered(work: Path) -> list[str]:
    """Set 0's primers, in the order the file lists them."""
    import csv

    path = work / "step4_improved_df.csv"
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return []
    key = "set_index" if "set_index" in rows[0] else None
    primer_key = next(k for k in ("primer", "primers", "sequence") if k in rows[0])
    out = []
    for row in rows:
        if key is not None and str(row[key]) not in ("0", "0.0"):
            continue
        value = row[primer_key]
        # One row may hold the whole set as a delimited string.
        out.extend(part for part in str(value).replace(";", ",").split(",") if part.strip())
    return [p.strip().upper() for p in out]


def _run(work: Path, geometry: str, reach: int, size: int) -> dict:
    params = _write_params(work, geometry, reach, size)
    _clear_step4(work)
    started = time.time()
    proc = subprocess.run(
        ["neoswga", "optimize", "-j", params.name],
        cwd=work,
        capture_output=True,
        text=True,
    )
    elapsed = time.time() - started
    if proc.returncode != 0:
        return {
            "geometry": geometry,
            "reach": reach,
            "size": size,
            "failed": True,
            "returncode": proc.returncode,
            "stderr_tail": proc.stderr[-2000:],
            "seconds": round(elapsed, 1),
        }
    summary = json.loads((work / "step4_improved_df_summary.json").read_text())
    primers = _delivered(work)
    # The optimizer metrics live under `metrics`; `coverage_geometries` is
    # written beside them. Reading the top level gives None for every figure,
    # which is why this asserts the keys it found rather than defaulting them.
    metrics = summary.get("metrics") or {}
    geometries = summary.get("coverage_geometries") or metrics.get("coverage_geometries") or {}
    if "fg_coverage" not in metrics:
        raise SystemExit(f"no fg_coverage in the summary for {geometry}@{reach} n={size}")
    return {
        "geometry": geometry,
        "reach": reach,
        "size": size,
        "requested": size,
        "delivered": len(primers),
        "primers": primers,
        "fg_coverage": metrics.get("fg_coverage"),
        "effective_fg_coverage": metrics.get("effective_fg_coverage"),
        "bg_coverage": metrics.get("bg_coverage"),
        "selectivity_density": metrics.get("selectivity_density"),
        "selectivity_ratio": metrics.get("selectivity_ratio"),
        "total_bg_sites": metrics.get("total_bg_sites"),
        "total_fg_sites": metrics.get("total_fg_sites"),
        "mean_gap": metrics.get("mean_gap"),
        "max_gap": metrics.get("max_gap"),
        "extension_reach": metrics.get("extension_reach"),
        "score": summary.get("score"),
        "geometries": {
            "symmetric": geometries.get("symmetric"),
            "directional": geometries.get("directional"),
            "inflation": geometries.get("inflation"),
            "reported_geometry": geometries.get("reported_geometry"),
        },
        "seconds": round(elapsed, 1),
    }


def _jaccard(left: list[str], right: list[str]) -> float | None:
    a, b = set(left), set(right)
    if not a and not b:
        return None
    return len(a & b) / len(a | b)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True, help="Directory holding params.json and the pool")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sizes", default=",".join(str(s) for s in SIZES))
    args = ap.parse_args()

    work = Path(args.work).resolve()
    sizes = [int(s) for s in args.sizes.split(",") if s.strip()]

    runs = []
    for size in sizes:
        for geometry, reach in ARMS:
            record = _run(work, geometry, reach, size)
            runs.append(record)
            print(
                f"  {geometry:>11} @{reach} n={size:>2} "
                f"delivered={record.get('delivered')} "
                f"fg={record.get('fg_coverage')} bg={record.get('bg_coverage')} "
                f"({record['seconds']}s)",
                flush=True,
            )

    # Jaccard against the symmetric arm at the configuration's own reach.
    by_key = {(r["geometry"], r["reach"], r["size"]): r for r in runs}
    for record in runs:
        baseline = by_key.get(("symmetric", 4000, record["size"]))
        if baseline is None or record.get("failed") or baseline.get("failed"):
            record["jaccard_against_symmetric_4000"] = None
            continue
        record["jaccard_against_symmetric_4000"] = _jaccard(
            baseline["primers"], record["primers"]
        )

    payload = {
        "pair": {
            "foreground": "E. coli K-12 MG1655, 4,639,675 bp, 1 record, circular",
            "background": "human chr21 (hg38), 46,709,983 bp, 18.6% N",
            "k": 12,
            "polymerase": "equiphi29",
            "reaction_temp": 43.0,
            "candidates": 1000,
        },
        "arms": [f"{g}@{r}" for g, r in ARMS],
        "sizes": sizes,
        "runs": runs,
    }
    Path(args.out).write_text(json.dumps(payload, indent=2))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
