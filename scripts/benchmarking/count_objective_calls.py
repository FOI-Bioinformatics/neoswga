"""How many objective evaluations does a real plan-pool run actually make?

The cost note prices a full-universe greedy scan. Nothing passes objective= to
optimize_greedy, so the reachable objective-scored searches are the swap loop
and the beam, both budgeted. This counts the calls instead of assuming them.
"""

import json, sys, time, pathlib

# The design directory to measure, and the selectivity density floor. A floor
# nothing violates leaves the repair unattempted and the objective unused, which
# is itself one of the results, so both are arguments.
DESIGN = pathlib.Path(
    sys.argv[2]
    if len(sys.argv) > 2
    else "examples/wolbachia_pool_design/objective_cost_2026-09-17/design"
)

import pandas as pd
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

params = json.loads((DESIGN / "params.json").read_text())
candidates = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg, bg = params["fg_prefixes"], params["bg_prefixes"]
cache = PositionCache(fg + bg, candidates)
_ensure_optimizers_registered()

FLOOR = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0  # 1.0 binds nothing
print(f"selectivity density floor {FLOOR}")
for budget in (10_000, 100_000):
    optimizer = OptimizerFactory.create(
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
            refinement_method="swap",
            swap_max_evaluations=budget,
            allow_dimer_relaxation=False,
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    calls = {"n": 0}
    real = optimizer.compute_pool_metrics

    def counted(primers, _real=real, _calls=calls):
        _calls["n"] += 1
        return _real(primers)

    optimizer.compute_pool_metrics = counted
    started = time.time()
    plan = plan_pool(
        optimizer,
        candidates,
        [12],
        [0.5],
        primer_length=12,
        coverage_metric="effective",
        min_selectivity_density=FLOOR,
    )
    elapsed = time.time() - started
    row = plan["rows"][0]
    stop = (row.get("repair") or {}).get("stop_reason")
    evals = (row.get("repair") or {}).get("evaluations")
    print(
        f"budget {budget:>7,}: {calls['n']:>6} objective calls, "
        f"{elapsed:6.1f}s wall, repair evaluations={evals}, stop={stop}"
    )
