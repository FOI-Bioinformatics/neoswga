"""Do the two Stage-2 selection terms buy anything in a reaction that can discriminate?

Run from the repo root, against a prepared design directory:

    PYTHONPATH=. python docs/validation/stage2_terms_2026-10-29/measure.py \
        --params <work>/params.json --out results.json

Both terms ship off, each with "unmeasured against a reaction" as the stated
reason. Until 2026-10-26 that was not merely undone but impossible: the
configured polymerase never reached the optimizer, so every run was phi29 at
30 C, and Known Issue 17 records that the phi29 pool cannot discriminate. This
measures them where the reaction can.

Two arms, because the terms act in different Stage 2s:

- `stage2_tm_weight` is an axis of hybrid's removal ranking, which is the
  `network` refinement -- the default. Measured there.
- `stage2_background_weight` is an axis of the SWAP score, so the arm that
  tests it sets `refinement_method="swap"`. Measured against the same pool.

Only the weight varies within an arm. Stage 2 has been deterministic since
2026-10-26 (`_CONNECTIVITY_SEED`), without which re-running one configuration
would give two or three different panels and none of this would mean anything.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import statistics
import subprocess
import sys

WEIGHTS = (0.0, 0.25, 0.5, 1.0)


def _run(params_path, overrides, work):
    """Run `optimize` once with `overrides` merged into params.json."""
    config = json.loads(pathlib.Path(params_path).read_text())
    config.update(overrides)
    tmp = pathlib.Path(work) / "params_arm.json"
    tmp.write_text(json.dumps(config, indent=2))
    proc = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "optimize", "-j", str(tmp)],
        capture_output=True,
        text=True,
        timeout=3600,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"optimize failed: {proc.stderr[-2000:]}")
    summary = json.loads((pathlib.Path(work) / "step4_improved_df_summary.json").read_text())
    import pandas as pd

    frame = pd.read_csv(pathlib.Path(work) / "step4_improved_df.csv")
    panel = frame[frame["set_index"] == 0]["primer"].astype(str).tolist()
    return summary, panel


def _tm_deviation(panel, conditions, target_tm):
    """Mean |effective Tm - target| over the panel, the quantity the Tm term ranks on.

    `stage2_tm_weight` prefers dropping the primer whose effective Tm sits
    furthest from `reaction_temp + 5`, so this is the figure it should improve
    if it does anything.
    """
    deviations = []
    for primer in panel:
        tm = conditions.calculate_effective_tm(primer)
        if tm is not None:
            deviations.append(abs(float(tm) - target_tm))
    return statistics.fmean(deviations) if deviations else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    base = json.loads(pathlib.Path(args.params).read_text())
    work = base["data_dir"]
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(
        polymerase=base.get("polymerase", "phi29"),
        temp=float(base.get("reaction_temp", 30.0)),
    )
    target_tm = float(base.get("reaction_temp", 30.0)) + 5.0

    arms = {
        # The Tm term rides the network refinement, which is the default.
        "stage2_tm_weight": {"refinement_method": "network"},
        # The host term is an axis of the swap score only. It REFUSES on the
        # production objective, which carries no host-coverage fraction, so it
        # is not swept here; the refusal is pinned in
        # tests/test_swap_refinement.py and recorded in
        # docs/validation/2026-10-29-stage2-background-weight-does-not-run.md.
    }

    results = {"target_tm": target_tm, "arms": {}}
    for key, extra in arms.items():
        rows = []
        print(f"\n=== {key}  ({extra['refinement_method']} refinement) ===")
        header = (
            f"{'weight':>7} {'fg_cov':>8} {'eff_cov':>8} {'bg_cov':>8} "
            f"{'bg_sites':>9} {'sel_dens':>9} {'tm_dev':>7} {'n':>3}"
        )
        print(header)
        for weight in WEIGHTS:
            overrides = dict(extra)
            overrides[key] = None if weight == 0.0 else weight
            summary, panel = _run(args.params, overrides, work)
            metrics = summary.get("metrics", {})
            dev = _tm_deviation(panel, conditions, target_tm)
            row = {
                "weight": weight,
                "fg_coverage": metrics.get("fg_coverage"),
                "effective_fg_coverage": metrics.get("effective_fg_coverage"),
                "bg_coverage": metrics.get("bg_coverage"),
                "total_bg_sites": metrics.get("total_bg_sites"),
                "selectivity_density": metrics.get("selectivity_density"),
                "mean_tm_deviation": dev,
                "panel_size": len(panel),
                "panel": panel,
            }
            rows.append(row)
            print(
                f"{weight:>7} {row['fg_coverage']:>8.4f} {row['effective_fg_coverage']:>8.4f} "
                f"{row['bg_coverage']:>8.4f} {row['total_bg_sites']:>9} "
                f"{row['selectivity_density']:>9.4f} {dev:>7.2f} {len(panel):>3}"
            )
        results["arms"][key] = rows

    pathlib.Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
