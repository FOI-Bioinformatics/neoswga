"""Does a larger search budget change the delivered Wolbachia panel?

The second half of plan step 249: compare runs at equal documented search
budgets, and record runs allowed to spend more time.

The retention benchmark stops at `filter`, so search budgets never entered it.
This regenerates a design directory, runs `score`, then runs `plan-pool` at
several swap-evaluation budgets with everything else held fixed, and compares
the delivered panels.

`post_gini` retention is used rather than `all_qc`: the shortlist `plan-pool`
draws from is the same 2,000 candidates in either mode, and the background index
is 18.9 MB against 443 MB, so this measures the search and not the indexing.
"""

import json
import pathlib
import re
import shutil
import subprocess
import sys
import time

EX = pathlib.Path("examples/wolbachia_pool_design")
OUT = EX / (sys.argv[1] if len(sys.argv) > 1 else "search_budget_2026-09-16")
BUDGETS = [
    int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else ("1000", "10000", "100000"))
]
SIZES = (8, 16)


def run(cmd, cwd, log):
    started = time.time()
    with open(log, "w") as handle:
        result = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified", *cmd],
            cwd=cwd,
            stdout=handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return time.time() - started, result.returncode


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    design = OUT / "design"
    if not (design / "step3_df.csv").exists():
        if design.exists():
            shutil.rmtree(design)
        design.mkdir(parents=True)
        base = json.loads((EX / "params.json").read_text())
        params = dict(base)
        params["data_dir"] = str(design.resolve())
        params["candidate_retention"] = "post_gini"
        for key in ("fg_prefixes", "bg_prefixes"):
            params[key] = [str(design.resolve() / pathlib.Path(p).name) for p in base[key]]
        (design / "params.json").write_text(json.dumps(params, indent=2))
        for step in ("count-kmers", "filter", "score"):
            seconds, code = run([step, "-j", "params.json"], design, design / f"{step}.log")
            assert code == 0, f"{step} failed\n{(design / f'{step}.log').read_text()[-2000:]}"
            print(f"{step}: {seconds:.0f}s", flush=True)

    results = {}
    for budget in BUDGETS:
        target = OUT / f"budget_{budget}"
        if target.exists():
            shutil.rmtree(target)
        seconds, code = run(
            [
                "plan-pool",
                "-j",
                "params.json",
                "--primer-length",
                "12",
                "--min-size",
                str(SIZES[0]),
                "--max-size",
                str(SIZES[1]),
                "--coverage-targets",
                "0.5",
                "--min-selectivity-density",
                "1.0",
                "--swap-max-evaluations",
                str(budget),
                "-o",
                str(target.resolve()),
            ],
            design,
            OUT / f"plan_{budget}.log",
        )
        assert (
            code == 0
        ), f"budget {budget} failed\n{(OUT / f'plan_{budget}.log').read_text()[-2500:]}"
        plan = json.loads((target / "pool_plan.json").read_text())
        results[budget] = {
            "seconds": round(seconds, 1),
            "rows": [
                {
                    "requested_size": r["requested_size"],
                    "size": r["size"],
                    "primers": sorted(r["primers"]),
                    "eligible": r["eligible"],
                    "coverage": r["coverage"],
                    "selectivity_density": r["selectivity_density"],
                    "background_sites": r["background_sites"],
                    "repair": r.get("repair"),
                }
                for r in plan["rows"]
            ],
        }
        print(f"budget {budget}: {seconds:.0f}s", flush=True)
        shutil.rmtree(target)

    (OUT / "results.json").write_text(json.dumps(results, indent=2))

    print("\n=== delivered panels by budget ===")
    reference = results[BUDGETS[0]]["rows"]
    for budget in BUDGETS:
        print(f"\nbudget {budget}  ({results[budget]['seconds']}s)")
        for row, base in zip(results[budget]["rows"], reference):
            same = row["primers"] == base["primers"]
            rep = (row.get("repair") or {}).get("method")
            print(
                f"  n={row['requested_size']:>3} size={row['size']:>3} "
                f"cov={row['coverage']!s:>8.8} dens={row['selectivity_density']!s:>8.8} "
                f"bg={row['background_sites']} repair={rep} "
                f"{'same as smallest budget' if same else 'DIFFERENT PANEL'}"
            )


main()
