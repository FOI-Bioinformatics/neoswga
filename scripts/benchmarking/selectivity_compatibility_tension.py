"""Are the most selective candidates unusually incompatible with each other?

The answer on the Wolbachia pool is: pairwise yes, usefully no. Selective
candidates are GC-richer and their pairwise compatibility is 8 to 10 points
lower than random, but the largest mutually compatible subset is not smaller --
at larger k it is larger. The "6 of 16" that prompted this sits inside the
random range of 6 to 8, so it was small-sample noise.

See docs/validation/no_search_headroom_on_this_pool_2026-09-18.md.

Six of the sixteen most selective candidates on the Wolbachia pool are mutually
dimer-free at max_dimer_bp 3. That could mean selectivity and compatibility are
in tension, or it could be what any sixteen candidates look like. The difference
decides whether more search work on specificity is worth doing, so it needs a
control -- the thing the invalid existence proof lacked.

Usage: tension.py <design_dir> [D]
"""

import importlib.util, json, random, statistics, sys, pathlib

DESIGN = pathlib.Path(sys.argv[1])
D = float(sys.argv[2]) if len(sys.argv) > 2 else 70.0
spec = importlib.util.spec_from_file_location(
    "sb", str(pathlib.Path(__file__).resolve().parent / "selectivity_budget.py")
)
sb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sb)

import pandas as pd
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.dimer_validator import DimerValidator
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_objective import PoolConstraints
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

params = json.loads((DESIGN / "params.json").read_text())
pool = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg, bg = params["fg_prefixes"], params["bg_prefixes"]
_ensure_optimizers_registered()
cache = PositionCache(fg + bg, pool)
opt = OptimizerFactory.create(
    "hybrid",
    cache,
    fg,
    params["fg_seq_lengths"],
    bg,
    params["bg_seq_lengths"],
    config=OptimizerConfig(
        max_dimer_bp=3,
        max_self_dimer_bp=4,
        extension_reach=3000,
        fg_circular=params.get("fg_circular", False),
        verbose=False,
    ),
    conditions=ReactionConditions(temp=30.0),
    polymerase="phi29",
)
budget = sb.SelectivityBudget.build(
    PoolConstraints(coverage_metric="effective", min_selectivity_density=D), opt, pool
)
val = DimerValidator(3, 4)


def stats(subset):
    """Compatible-pair fraction and the largest dimer-free subset a greedy finds."""
    n = len(subset)
    pairs = total = 0
    for i, a in enumerate(subset):
        for b in subset[i + 1 :]:
            total += 1
            if val.is_compatible(a, b):
                pairs += 1
    best = []
    for start in range(min(n, 40)):
        panel = []
        for c in subset[start:] + subset[:start]:
            if not val.has_self_dimer(c) and all(val.is_compatible(c, q) for q in panel):
                panel.append(c)
        if len(panel) > len(best):
            best = panel
    gc = statistics.mean((c.count("G") + c.count("C")) / len(c) for c in subset)
    return pairs / total if total else 1.0, len(best), gc


ranked = sorted(pool, key=budget.slack, reverse=True)
rng = random.Random(20260918)
print(f"pool {len(pool)}, D={D}, max_dimer_bp=3\n")
print(f"{'k':>5} {'set':>12} {'compatible':>11} {'max dimer-free':>15} {'mean GC':>8}")
for k in (16, 32, 64, 128, 256):
    frac, clique, gc = stats(ranked[:k])
    print(f"{k:>5} {'top by slack':>12} {frac:>10.1%} {clique:>15} {gc:>8.3f}")
    controls = [stats(rng.sample(pool, k)) for _ in range(12)]
    cf = [c[0] for c in controls]
    cc = [c[1] for c in controls]
    cg = [c[2] for c in controls]
    print(
        f"{k:>5} {'random x12':>12} {statistics.mean(cf):>10.1%} "
        f"{statistics.mean(cc):>15.1f} {statistics.mean(cg):>8.3f}"
        f"   (range {min(cc)}-{max(cc)})"
    )
