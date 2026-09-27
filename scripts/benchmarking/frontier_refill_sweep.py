"""Does refilling the frontier change the delivered panel on the real pool?

Increment 5's acceptance measurement. The frontier opens at the shortlist and
the inventory holds every candidate that cleared hard QC; each refill doubles,
so four reach the universe.

Usage: frontier_refill_sweep.py [density_floor] [design_dir]

The floor matters more here than anywhere else: a row that qualifies never
refills, so a floor the shortlist already meets measures nothing. On the bundled
Wolbachia design 40 and 60 qualify immediately and 100 does not. The design
directory must have been through count-kmers, filter and score.

See docs/validation/frontier_refill_2026-09-17.md.
"""

import json, sys, time, pathlib

DESIGN = pathlib.Path(
    sys.argv[2]
    if len(sys.argv) > 2
    else "examples/wolbachia_pool_design/objective_cost_2026-09-17/design"
)
import pandas as pd
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.candidate_source import open_candidate_source
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered
from types import SimpleNamespace

FLOOR = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
SIZES = [12]

params = json.loads((DESIGN / "params.json").read_text())
candidates = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg, bg = params["fg_prefixes"], params["bg_prefixes"]
conditions = build_reaction_conditions(SimpleNamespace(**params))
_ensure_optimizers_registered()

print(f"floor {FLOOR}, sizes {SIZES}\n")
for refills in (0, 4):
    source = open_candidate_source(
        str(DESIGN), conditions.fingerprint(), [12], frontier=len(candidates)
    )
    described = source.describe()
    universe = described["universe"]
    # The cache has to cover whatever the frontier can reach.
    everything = source._provider._eligible_in_search_order()
    started = time.time()
    cache = PositionCache(fg + bg, everything, on_missing="warn")
    load_seconds = time.time() - started
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
            fg_circular=params.get("fg_circular", False),
            refinement_method="swap",
            max_frontier_refills=refills,
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    started = time.time()
    plan = plan_pool(
        optimizer,
        source,
        SIZES,
        [0.5],
        primer_length=12,
        coverage_metric="effective",
        min_selectivity_density=FLOOR,
    )
    elapsed = time.time() - started
    row = plan["rows"][0]
    print(
        f"refills allowed {refills}: universe {universe}, cache load {load_seconds:.0f}s, "
        f"design {elapsed:.0f}s"
    )
    print(f"  used {row['frontier_refills']} refills, ended {row['candidates_exhausted']}")
    print(
        f"  coverage {row.get('coverage'):.6f}  density {row.get('selectivity_density'):.3f}  "
        f"feasible {not row.get('failed_constraints')}"
    )
    print(f"  examined {plan['candidate_source']['examined']} of {universe}")
    print(f"  primers {sorted(row.get('primers') or [])[:3]} ... n={len(row.get('primers') or [])}")
    json.dump(
        {
            "refills": refills,
            "row": {k: v for k, v in row.items() if k != "repair"},
            "source": plan["candidate_source"],
        },
        open(f"refill_{refills}.json", "w"),
        indent=1,
    )
