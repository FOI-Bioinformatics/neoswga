"""Does mixing oligo lengths buy anything? A second pair. (Known Issue 20)

Phase 4 of `docs/reports/PLAN_2026-11-02_known_issues.md`. Run from the
repository root against a directory whose k-mer TABLES are already built for
every k in `--lengths`:

    PYTHONPATH=. python docs/validation/mixed_length_second_pair_2026-11-02/measure.py \
        --params <work>/params.json --out results.json

Known Issue 20 says mixed lengths "on the one pair measured bought nothing":
on Prevotella against human chr21 a mixed k 10-12 design landed BETWEEN the
single-length designs on both axes at panel sizes 8 and 20, and the entry
concludes "choosing k matters far more than choosing whether to mix". One pair
is what this project's own convention calls "no benefit demonstrated". This
adds a second pair so the claim either holds at a larger size or does not.

Four arms per panel size -- each single length in `--lengths`, then all of them
mixed -- and the question is whether the mixed arm lands between the singles on
coverage and on specificity.

**Each arm resets the artifacts a previous arm wrote** (step2, step3, step4,
the candidate inventory, the funnel and the manifest) and keeps only the k-mer
tables, which do not depend on the arm. Reusing a directory without that reset
is a measured confound: arms run in sequence otherwise disagree with the same
arm run in a fresh directory, because the inventory carries across.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys

RESET = (
    "step2_df.csv",
    "step3_df.csv",
    "step4_improved_df.csv",
    "step4_improved_df_summary.json",
    "step4_improved_df_audit.json",
    "step4_improved_df_validation.json",
    "candidate_inventory.sqlite",
    "filter_stats.json",
    "run_manifest.json",
)


def _reset(data_dir: pathlib.Path) -> None:
    for name in RESET:
        (data_dir / name).unlink(missing_ok=True)


def _step(step: str, params: pathlib.Path) -> None:
    proc = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", step, "-j", str(params)],
        capture_output=True,
        text=True,
        timeout=7200,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{step} failed: {proc.stderr[-2000:]}")


def _arm(base_params: dict, data_dir: pathlib.Path, min_k, max_k, size):
    _reset(data_dir)
    config = dict(base_params)
    config.update({"min_k": min_k, "max_k": max_k, "num_primers": size, "target_set_size": size})
    path = data_dir / "params_mixed_arm.json"
    path.write_text(json.dumps(config, indent=2))
    for step in ("filter", "prepare-candidates", "optimize"):
        _step(step, path)
    summary = json.loads((data_dir / "step4_improved_df_summary.json").read_text())
    funnel = json.loads((data_dir / "filter_stats.json").read_text())
    m = summary["metrics"]
    primers = list(summary.get("primers") or [])
    return {
        "lengths": f"{min_k}-{max_k}" if min_k != max_k else str(min_k),
        "size": size,
        "delivered": len(primers),
        "oligo_lengths_delivered": sorted({len(p) for p in primers}),
        "candidates": funnel.get("final_candidates"),
        "fg_coverage": m.get("fg_coverage"),
        "bg_coverage": m.get("bg_coverage"),
        "total_bg_sites": m.get("total_bg_sites"),
        "selectivity_density": m.get("selectivity_density"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--lengths", default="10,11,12")
    parser.add_argument("--sizes", default="8,20")
    args = parser.parse_args()

    base = json.loads(pathlib.Path(args.params).read_text())
    data_dir = pathlib.Path(base["data_dir"])
    lengths = [int(x) for x in args.lengths.split(",")]
    sizes = [int(x) for x in args.sizes.split(",")]

    rows = []
    for size in sizes:
        for k in lengths:
            rows.append(_arm(base, data_dir, k, k, size))
        rows.append(_arm(base, data_dir, min(lengths), max(lengths), size))

    out = {"lengths": lengths, "sizes": sizes, "arms": rows}
    # Does the mixed arm land between the singles, as the first pair found?
    for size in sizes:
        at = [r for r in rows if r["size"] == size]
        mixed = next(r for r in at if "-" in r["lengths"])
        singles = [r for r in at if "-" not in r["lengths"]]
        out[f"between_at_n{size}"] = {
            "coverage": (
                min(s["fg_coverage"] for s in singles)
                <= mixed["fg_coverage"]
                <= max(s["fg_coverage"] for s in singles)
            ),
            "specificity": (
                min(s["selectivity_density"] for s in singles)
                <= mixed["selectivity_density"]
                <= max(s["selectivity_density"] for s in singles)
            ),
        }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=2))

    print(f"{'arm':>7} {'n':>3} {'cand':>6} {'fg_cov':>8} {'bg_cov':>9} {'bg':>4} {'sel_dens':>9}")
    for r in rows:
        print(
            f"{r['lengths']:>7} {r['size']:>3} {str(r['candidates']):>6} "
            f"{r['fg_coverage']:>8.4f} {r['bg_coverage']:>9.6f} "
            f"{str(r['total_bg_sites']):>4} {r['selectivity_density']:>9.3f}"
        )
    for size in sizes:
        b = out[f"between_at_n{size}"]
        print(
            f"  n={size}: mixed between singles -- coverage {b['coverage']}, "
            f"specificity {b['specificity']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
