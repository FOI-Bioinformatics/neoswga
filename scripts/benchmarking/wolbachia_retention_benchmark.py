"""Compare candidate_retention modes on the real Wolbachia/Drosophila pair.

The `legacy` mode was removed on 2026-09-16 after this comparison was made, so
re-running this now measures `post_gini` against `all_qc`. The published legacy
figures are in docs/validation/wolbachia_retention_benchmark_2026-09-16.md.
"""

import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time

EX = pathlib.Path("examples/wolbachia_pool_design")
OUT = EX / "benchmark_2026-09-16"
OUT.mkdir(exist_ok=True)
base = json.loads((EX / "params.json").read_text())


def run(cmd, cwd, log):
    t0 = time.time()
    with open(log, "w") as fh:
        r = subprocess.run(
            ["/usr/bin/time", "-l", sys.executable, "-m", "neoswga.cli_unified", *cmd],
            cwd=cwd,
            stdout=fh,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return time.time() - t0, r.returncode


def peak_mb(log):
    m = re.search(r"(\d+)\s+maximum resident set size", pathlib.Path(log).read_text())
    return int(m.group(1)) / 1e6 if m else None


def counts_from(log):
    text = pathlib.Path(log).read_text()
    out = {}
    m = re.search(r"Design counts: (.+)", text)
    if m:
        for part in m.group(1).split():
            k, _, v = part.partition("=")
            out[k] = int(v)
    m = re.search(r"Background index: (\d+) of (\d+)", text)
    if m:
        out["bg_indexed"] = int(m.group(1))
    return out


results = {}
for mode in ("post_gini", "all_qc"):
    d = OUT / mode
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    p = dict(base)
    p["data_dir"] = str(d.resolve())
    p["candidate_retention"] = mode
    for key in ("fg_prefixes", "bg_prefixes"):
        p[key] = [str(d.resolve() / pathlib.Path(x).name) for x in base[key]]
    (d / "params.json").write_text(json.dumps(p, indent=2))

    t_count, rc1 = run(["count-kmers", "-j", "params.json"], d, d / "count.log")
    t_filter, rc2 = run(["filter", "-j", "params.json"], d, d / "filter.log")
    assert rc1 == 0 and rc2 == 0, (
        f"{mode}: rc={rc1},{rc2}\n" + (d / "filter.log").read_text()[-1500:]
    )

    stats = json.loads((d / "filter_stats.json").read_text())
    h5 = {f.name: f.stat().st_size for f in d.glob("*positions.h5")}
    results[mode] = {
        "count_seconds": round(t_count, 1),
        "filter_seconds": round(t_filter, 1),
        "count_peak_mb": peak_mb(d / "count.log"),
        "filter_peak_mb": peak_mb(d / "filter.log"),
        "funnel": stats,
        "design_counts": counts_from(d / "filter.log"),
        "step2_rows": len((d / "step2_df.csv").read_text().strip().split("\n")) - 1,
        "position_h5_bytes": h5,
        "inventory_bytes": (d / "candidate_inventory.sqlite").stat().st_size,
    }
    print(f"{mode} done: filter {t_filter:.0f}s", flush=True)

(OUT / "results.json").write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
