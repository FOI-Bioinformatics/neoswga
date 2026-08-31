"""Optimality gap of the greedy set cover on tier M, against an exact bound.

Compares, on one candidate pool and one bin set:
  - what `DominatingSetOptimizer.optimize_greedy` selects
  - the exact ILP optimum for the same fixed budget
  - the LP relaxation (a valid upper bound even when the ILP is stopped early)
  - random sets of the same size, as a floor

Everything is measured in covered bases over the same BipartiteGraph bins, so
the three numbers are directly comparable.
"""

import json
import logging
import pathlib
import random
import sys
import time

import pandas as pd

logging.basicConfig(level=logging.WARNING)

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from neoswga.core.coverage import polymerase_extension_reach  # noqa: E402
from neoswga.core.dominating_set_optimizer import (  # noqa: E402
    DominatingSetOptimizer,
    coverage_bin_size,
)
from neoswga.core.position_cache import PositionCache  # noqa: E402
from max_coverage_bound import build_bin_coverage, coverage_bounds  # noqa: E402

WORK = sys.argv[1]
SIZES = [int(s) for s in sys.argv[2].split(",")]
N_RANDOM = int(sys.argv[3]) if len(sys.argv) > 3 else 100

params = json.load(open(f"{WORK}/params.json"))
fg_prefixes = params["fg_prefixes"]
fg_lengths = params["fg_seq_lengths"]

cands = pd.read_csv(f"{WORK}/step3_df.csv")["primer"].astype(str).tolist()
print(f"candidates: {len(cands)}", flush=True)

reach = polymerase_extension_reach(params["polymerase"], "realistic")
bin_size = coverage_bin_size(10000, reach)
print(f"reach={reach} bin_size={bin_size}", flush=True)

t0 = time.time()
cache = PositionCache(fname_prefixes=fg_prefixes, primers=cands)
print(f"cache built in {time.time()-t0:.1f}s", flush=True)

opt = DominatingSetOptimizer(
    cache, fg_prefixes, fg_lengths, bin_size=bin_size, extension_reach=reach
)

primer_to_bins, bin_weight, total_bases = build_bin_coverage(opt, cands)
print(f"primers with bins: {len(primer_to_bins)}  bins: {len(bin_weight)}  "
      f"genome: {total_bases}", flush=True)


def bases_of(selected):
    bins = set()
    for p in selected:
        bins |= primer_to_bins.get(p, set())
    return sum(bin_weight[b] for b in bins)


rng = random.Random(42)
pool = list(primer_to_bins)
rows = []

for S in SIZES:
    t = time.time()
    g = opt.optimize_greedy(cands, max_primers=S, verbose=False)
    greedy_s = time.time() - t
    greedy_cov = bases_of(g.get("primers", [])) / total_bases

    b = coverage_bounds(opt, cands, budget=S, max_seconds=600)
    ilp, lp = b["ilp"], b["lp"]

    rnd = sorted(bases_of(rng.sample(pool, min(S, len(pool)))) / total_bases
                 for _ in range(N_RANDOM))
    rnd_med = rnd[len(rnd) // 2]
    rnd_best = rnd[-1]
    better = sum(1 for r in rnd if r >= greedy_cov)

    row = {
        "budget": S,
        "greedy_coverage": round(greedy_cov, 6),
        "greedy_seconds": round(greedy_s, 2),
        "ilp_optimum": round(ilp.coverage, 6),
        "ilp_proven_optimal": ilp.proven_optimal,
        "ilp_status": ilp.status,
        "ilp_seconds": round(ilp.seconds, 1),
        "lp_bound": round(lp.coverage, 6),
        "gap_pct": round((ilp.coverage - greedy_cov) / ilp.coverage * 100, 3)
        if ilp.coverage else None,
        "random_median": round(rnd_med, 6),
        "random_best_of_n": round(rnd_best, 6),
        "random_beating_greedy": better,
        "n_random": N_RANDOM,
    }
    rows.append(row)
    print(json.dumps(row), flush=True)

json.dump(rows, open(f"{WORK}/gap_results.json", "w"), indent=2)
print("wrote gap_results.json", flush=True)
