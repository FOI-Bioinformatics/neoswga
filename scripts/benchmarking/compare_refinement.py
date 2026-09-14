"""Compare refinement on saved pools in fresh, sequential worker processes.

Measures the library optimizer, cache loading and independent panel evaluation;
this is not a CLI or upstream candidate-filtering benchmark. Inputs are read
only. Each worker reports its own peak RSS, avoiding cumulative child maxima.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def worker(case):
    import logging
    import random
    import resource
    from types import SimpleNamespace

    import numpy as np
    import pandas as pd

    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.dimer import is_dimer_fast
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.hybrid_optimizer import HybridOptimizer
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import build_reaction_conditions

    logging.basicConfig(level=logging.INFO)
    random.seed(42)
    np.random.seed(42)
    started = time.monotonic()
    params_path = Path(case["params"])
    params = json.loads(params_path.read_text())

    def resolve(p):
        return str((params_path.parent / p).resolve())

    fg = [resolve(p) for p in params["fg_prefixes"]]
    bg = [resolve(p) for p in params.get("bg_prefixes", [])]
    candidates = list(dict.fromkeys(pd.read_csv(case["pool"])["primer"].tolist()))
    cache = PositionCache(fg + bg, candidates)
    cache_seconds = time.monotonic() - started
    conditions = build_reaction_conditions(SimpleNamespace(**params))
    optimizer = HybridOptimizer(
        cache,
        fg,
        params["fg_seq_lengths"],
        bg,
        params.get("bg_seq_lengths", []),
        coverage_reach=case["reach"],
        max_dimer_bp=params["max_dimer_bp"],
        polymerase=params["polymerase"],
        conditions=conditions,
        background_pruning=case["background_aware"],
        refinement_method=case["mode"],
        swap_max_evaluations=case["evaluations"],
        swap_max_seconds=case["search_seconds"],
        allow_dimer_relaxation=False,
    )
    begin = time.monotonic()
    result = optimizer.optimize(
        candidates, final_count=case["size"], verbose=False, apply_polymerase_multiplier=False
    )
    optimize_seconds = time.monotonic() - begin
    evaluator = DominatingSetAdapter(
        cache,
        fg,
        params["fg_seq_lengths"],
        bg,
        params.get("bg_seq_lengths", []),
        config=OptimizerConfig(
            extension_reach=case["reach"],
            fg_circular=params.get("fg_circular", False),
            max_dimer_bp=params["max_dimer_bp"],
            verbose=False,
        ),
        conditions=conditions,
    )
    metrics = evaluator.compute_metrics(result.primers)
    violations = sum(
        is_dimer_fast(a, b, params["max_dimer_bp"])
        for i, a in enumerate(result.primers)
        for b in result.primers[:i]
    )
    # Same binned base objective used by swap selection and exact benchmarks.
    regions = optimizer._coverage_bins_by_primer(result.primers)
    covered = {optimizer._bin_key(r): r.end - r.start for owned in regions.values() for r in owned}
    return {
        **case,
        "status": "completed",
        "candidates": len(candidates),
        "primers": result.primers,
        "delivered": len(result.primers),
        "violating_pairs": int(violations),
        "binned_covered_bases": sum(covered.values()),
        "fg_coverage": metrics.fg_coverage,
        "effective_fg_coverage": metrics.effective_fg_coverage,
        "background_sites": metrics.total_bg_sites,
        "selectivity_density": metrics.selectivity_density,
        "max_gap": metrics.max_gap,
        "cache_seconds": cache_seconds,
        "optimize_seconds": optimize_seconds,
        "stage2_seconds": result.runtime_stage2,
        "worker_seconds": time.monotonic() - started,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        / (2**20 if sys.platform == "darwin" else 1024),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--params", type=Path)
    parser.add_argument("--pools", nargs="+", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--modes", nargs="+", choices=["network", "swap"], default=["network", "swap"]
    )
    parser.add_argument("--sizes", nargs="+", type=int, default=[6, 12])
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=90)
    parser.add_argument("--evaluations", type=int, default=10000)
    parser.add_argument("--search-seconds", type=float, default=10)
    parser.add_argument("--background-aware", action="store_true")
    args = parser.parse_args()
    if args.worker:
        case = json.loads(args.worker.read_text())
        result = worker(case)
        args.worker.with_suffix(".result.json").write_text(
            json.dumps(result, indent=2, allow_nan=False)
        )
        return
    if not args.params or not args.pools or not args.output:
        parser.error("--params, --pools and --output are required")
    if args.repeats < 1 or any(size < 1 for size in args.sizes):
        parser.error("Repeats and panel sizes must be positive")
    if args.timeout <= 0 or args.evaluations < 0 or args.search_seconds < 0:
        parser.error("Timeout must be positive and search budgets non-negative")
    if args.output.exists():
        parser.error("Output directory already exists; choose a new one to preserve evidence")
    args.output.mkdir(parents=True)
    manifest = {
        "params": str(args.params.resolve()),
        "params_sha256": sha256(args.params),
        "pools": {str(p.resolve()): sha256(p) for p in args.pools},
        "python": sys.version,
        "platform": sys.platform,
        "command": sys.argv,
        "benchmark_sha256": sha256(__file__),
        "git_head": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_sha256": {
            str(p.relative_to(ROOT)): sha256(p) for p in (ROOT / "neoswga/core").glob("*.py")
        },
    }
    params = json.loads(args.params.read_text())
    indexes = {}
    for prefix in params["fg_prefixes"] + params.get("bg_prefixes", []):
        absolute = (args.params.parent / prefix).resolve()
        for index_path in absolute.parent.glob(absolute.name + "_*mer_positions.h5"):
            indexes[str(index_path)] = sha256(index_path)
    manifest["position_index_sha256"] = indexes
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2))
    env = {
        **os.environ,
        "PYTHONHASHSEED": "42",
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
    }
    index = 0
    with (args.output / "results.jsonl").open("w") as results:
        for pool in args.pools:
            for size in args.sizes:
                for repeat in range(args.repeats):
                    # Alternate order to reduce a systematic warm-cache advantage.
                    modes = args.modes if repeat % 2 == 0 else list(reversed(args.modes))
                    for mode in modes:
                        case = dict(
                            params=str(args.params.resolve()),
                            pool=str(pool.resolve()),
                            size=size,
                            repeat=repeat,
                            mode=mode,
                            background_aware=args.background_aware,
                            reach=3000,
                            evaluations=args.evaluations,
                            search_seconds=args.search_seconds,
                        )
                        path = args.output / f"case_{index:03d}.json"
                        path.write_text(json.dumps(case))
                        start = time.monotonic()
                        with path.with_suffix(".log").open("w") as log:
                            try:
                                process = subprocess.run(
                                    [
                                        sys.executable,
                                        str(Path(__file__).resolve()),
                                        "--worker",
                                        str(path.resolve()),
                                    ],
                                    cwd=ROOT,
                                    stdout=log,
                                    stderr=log,
                                    env=env,
                                    timeout=args.timeout,
                                )
                                record = (
                                    json.loads(path.with_suffix(".result.json").read_text())
                                    if process.returncode == 0
                                    else {
                                        **case,
                                        "status": "error",
                                        "returncode": process.returncode,
                                    }
                                )
                            except subprocess.TimeoutExpired:
                                record = {**case, "status": "timeout"}
                        record["process_seconds"] = time.monotonic() - start
                        results.write(json.dumps(record, allow_nan=False) + "\n")
                        results.flush()
                        print(
                            json.dumps(
                                {
                                    k: record.get(k)
                                    for k in [
                                        "pool",
                                        "size",
                                        "repeat",
                                        "mode",
                                        "status",
                                        "delivered",
                                        "fg_coverage",
                                        "process_seconds",
                                    ]
                                }
                            ),
                            flush=True,
                        )
                        index += 1


if __name__ == "__main__":
    main()
