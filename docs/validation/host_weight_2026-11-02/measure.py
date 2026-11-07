"""Does a weighted host axis beat the lexicographic tie-break? (Known Issue 14)

Phase 2.2 of `docs/reports/PLAN_2026-11-02_known_issues.md`. Run from the
repository root against a prepared design directory with a BACKGROUND:

    PYTHONPATH=. python docs/validation/host_weight_2026-11-02/measure.py \
        --params <work>/params.json --out results.json

The question `2026-10-29-stage2-background-weight-does-not-run.md` could not
ask. That record measured the key's DIRECTION against a stub objective and said
plainly that no production objective carried the host-coverage fraction, so the
key never ran on a real design. `PoolMetrics` carries the field since
2026-11-02 (+48% per evaluation, and cheaper than the full metrics the search
already falls back to for a host limit), so the key runs now and the open
question is whether it is worth switching on.

`stage2_background_weight` is an axis of the SWAP score, so every arm sets
`refinement_method="swap"`. Only the weight varies. Stage 2 has been
deterministic since 2026-10-26, so one configuration gives one panel and a
single run per weight is the whole measurement rather than a sample.

Reported per arm: what the panel covers of the target, what it covers of the
HOST (the quantity being weighted), the background site count and the
selectivity density a design is accepted on.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

WEIGHTS = (0.0, 0.25, 0.5, 1.0)


def _run(params_path, overrides, work):
    """One `optimize` run with `overrides` merged into params.json."""
    config = json.loads(pathlib.Path(params_path).read_text())
    config.update(overrides)
    tmp = pathlib.Path(work) / "params_arm.json"
    tmp.write_text(json.dumps(config, indent=2))
    proc = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "optimize", "-j", str(tmp)],
        capture_output=True,
        text=True,
        timeout=7200,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"optimize failed: {proc.stderr[-2500:]}")
    summary = json.loads((pathlib.Path(work) / "step4_improved_df_summary.json").read_text())
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--size", type=int, default=12)
    args = parser.parse_args()

    params_path = pathlib.Path(args.params)
    work = json.loads(params_path.read_text())["data_dir"]

    rows = []
    for weight in WEIGHTS:
        overrides = {
            "refinement_method": "swap",
            "num_primers": args.size,
            "target_set_size": args.size,
        }
        # Unset rather than 0.0 for the control, so the control is the SHIPPED
        # default path and not a weighted axis that happens to weigh nothing.
        if weight:
            overrides["stage2_background_weight"] = weight
        summary = _run(params_path, overrides, work)
        m = summary["metrics"]
        rows.append(
            {
                "weight": weight or None,
                "primers": len(summary.get("primers") or []),
                "fg_coverage": m.get("fg_coverage"),
                "effective_fg_coverage": m.get("effective_fg_coverage"),
                "bg_coverage": m.get("bg_coverage"),
                "total_bg_sites": m.get("total_bg_sites"),
                "selectivity_density": m.get("selectivity_density"),
                "panel": sorted(summary.get("primers") or []),
            }
        )

    control = rows[0]
    for row in rows:
        row["panel_differs_from_control"] = row["panel"] != control["panel"]

    out = {"size": args.size, "arms": rows}
    pathlib.Path(args.out).write_text(json.dumps(out, indent=2))
    print(json.dumps({k: v for k, v in out.items() if k != "arms"}, indent=2))
    for row in rows:
        print(
            f"  weight={str(row['weight']):>5}  n={row['primers']:>2}  "
            f"fg={row['fg_coverage']:.4f}  bg={row['bg_coverage']:.6f}  "
            f"bg_sites={row['total_bg_sites']:>3}  "
            f"sel_dens={row['selectivity_density']:.3f}  "
            f"moved={row['panel_differs_from_control']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
