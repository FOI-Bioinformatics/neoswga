"""Where does one objective evaluation spend its time? One variable at a time.

`measure_objective_cost.py` prices one evaluation. This one attributes the
price. It was written because the cost note blamed the 144 Mb background, and
the measurement says the background is 16% of it while the occupancy weighting
is 95%.

Point DESIGN at a directory that has been through count-kmers, filter and
score. Nothing is written; the run takes about a minute once the index is
loaded, most of it in PositionCache construction.

See docs/validation/objective_evaluation_cost_2026-09-17.md and
docs/validation/parallelism_opportunities_2026-09-17.md.
"""

import json, sys, time, pathlib

# The design directory to measure. Override with argv[1].
DESIGN = pathlib.Path(
    sys.argv[1]
    if len(sys.argv) > 1
    else "examples/wolbachia_pool_design/objective_cost_2026-09-17/design"
)

import pandas as pd
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_objective import PoolConstraints, PoolObjective
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

params = json.loads((DESIGN / "params.json").read_text())
candidates = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg, bg = params["fg_prefixes"], params["bg_prefixes"]
cache = PositionCache(fg + bg, candidates)
_ensure_optimizers_registered()


def build(conditions, with_bg=True):
    return OptimizerFactory.create(
        "hybrid",
        cache,
        fg,
        params["fg_seq_lengths"],
        bg if with_bg else [],
        params["bg_seq_lengths"] if with_bg else [],
        config=OptimizerConfig(
            max_dimer_bp=3, max_self_dimer_bp=4, extension_reach=3000, verbose=False
        ),
        conditions=conditions,
        polymerase="phi29",
    )


def timeit(objective, panel, probes, n=200):
    objective.coverage(panel)
    started = time.time()
    for candidate in probes[:n]:
        objective.coverage([*panel, candidate])
    return (time.time() - started) / n * 1e3


panel = candidates[:24]
probes = candidates[24:]

cases = [
    (
        "effective coverage, background present",
        build(ReactionConditions(temp=30.0), True),
        "effective",
    ),
    ("raw coverage, background present", build(None, True), "raw"),
    ("effective coverage, NO background", build(ReactionConditions(temp=30.0), False), "effective"),
]
print(f"{'case':<42} {'ms/call':>8}")
for label, opt, metric in cases:
    obj = PoolObjective(opt.compute_pool_metrics, PoolConstraints(coverage_metric=metric))
    print(f"{label:<42} {timeit(obj, panel, probes):>8.2f}")

# Panel-size scaling: if the cost is the target-length loop it is linear in panel size.
opt = build(ReactionConditions(temp=30.0), True)
obj = PoolObjective(opt.compute_pool_metrics, PoolConstraints(coverage_metric="effective"))
print(f"\n{'panel size':>10} {'ms/call':>8}")
for size in (6, 12, 24, 48):
    print(f"{size:>10} {timeit(obj, candidates[:size], probes, n=100):>8.2f}")
