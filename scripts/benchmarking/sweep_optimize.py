"""Measure how `neoswga optimize` scales in set size, per method.

Drives the CLI, because that is what users run: a library-level harness would
miss the dispatch, the position-cache build and the metric computation that
dominate a real invocation.

Records wall time, peak RSS (via resource.getrusage on the child) and the
coverage the run reports, so the scaling curve and the quality it buys are read
off the same invocation.
"""

from __future__ import annotations

import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path


def run_one(workdir: Path, method: str, size: int, extra: list[str], tag: str) -> dict:
    out = workdir / f"run_{tag}"
    out.mkdir(exist_ok=True)

    # The set size goes on the CLI, not into params.json. `param_validator`
    # caps num_primers/target_set_size at 50 (param_validator.py:95-96) while
    # params.schema.json allows 200, and the validator gate runs first -- so a
    # params.json asking for 64 is rejected outright while `--num-primers 64`
    # runs fine. Using the flag keeps the sweep on the path that works.
    params = json.loads((workdir / "params.json").read_text())
    params["data_dir"] = "./"
    params["num_primers"] = min(size, 8)
    params["target_set_size"] = min(size, 8)
    (workdir / f"params_{tag}.json").write_text(json.dumps(params, indent=1))

    cmd = [
        "neoswga", "optimize",
        "-j", f"params_{tag}.json",
        "--optimization-method", method,
        "--num-primers", str(size),
        "--seed", "42",
        *extra,
    ]

    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    proc = subprocess.run(
        cmd, cwd=workdir, capture_output=True, text=True, env={**os.environ}
    )
    elapsed = time.time() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)

    # maxrss is bytes on macOS, kilobytes on Linux.
    scale = 1 << 20 if sys.platform == "darwin" else 1 << 10
    peak_mb = max(after.ru_maxrss - before.ru_maxrss, after.ru_maxrss) / scale

    summary_path = workdir / "step4_improved_df_summary.json"
    summary = {}
    if summary_path.exists():
        try:
            summary = json.loads(summary_path.read_text())
        except Exception:
            pass

    metrics = summary.get("metrics", {}) or {}
    rec = {
        "tag": tag,
        "method": method,
        "requested_size": size,
        "extra": " ".join(extra),
        "returncode": proc.returncode,
        "seconds": round(elapsed, 3),
        "peak_rss_mb": round(peak_mb, 1),
        "n_primers": summary.get("num_primers"),
        "score": summary.get("score"),
        "fg_coverage": metrics.get("fg_coverage"),
        "effective_fg_coverage": metrics.get("effective_fg_coverage"),
        "selectivity_ratio": metrics.get("selectivity_ratio"),
        "selectivity_density": metrics.get("selectivity_density"),
        "gap_gini": metrics.get("gap_gini"),
        "max_gap": metrics.get("max_gap"),
        "primers": summary.get("primers"),
        "stderr_tail": proc.stderr[-400:] if proc.returncode else "",
    }
    if summary_path.exists():
        summary_path.rename(out / "step4_improved_df_summary.json")
    return rec


def main() -> None:
    workdir = Path(sys.argv[1]).resolve()
    results_path = workdir / "sweep_results.jsonl"

    sizes = [int(s) for s in sys.argv[2].split(",")]
    methods = sys.argv[3].split(",")

    with results_path.open("a") as fh:
        for method in methods:
            for size in sizes:
                tag = f"{method}_S{size}"
                print(f"[run] {tag}", flush=True)
                rec = run_one(workdir, method, size, [], tag)
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                print(
                    f"[done] {tag} {rec['seconds']}s rss={rec['peak_rss_mb']}MB "
                    f"n={rec['n_primers']} cov={rec['fg_coverage']}",
                    flush=True,
                )


if __name__ == "__main__":
    main()
