"""Which selectivity floor can each candidate universe actually reach?

Phase 5 of the plan for `docs/validation/pipeline_audit_2026-09-16/`. The
retention benchmark stopped after `filter` and so could not say whether
retaining candidates changes a delivered panel. With Phase 4 landed it can be
asked directly: raise a specificity floor until each universe fails, and compare
where each one stops.

The three universes are the shortlist (frontier pinned at 2,000, no refills),
the post-Gini inventory, and all hard-QC candidates. Budgets are held equal
across universes; only the universe and the floor vary.

Usage: floor_boundary.py <design_dir> <label> [floors]

The design directory must have been through count-kmers, filter and score, and
its `candidate_retention` decides the universe: `post_gini` gives 20,670 on the
Wolbachia pair and `all_qc` gives 491,836. The shortlist rows are produced from
the same directory by pinning the frontier and allowing no refills, so all three
universes are compared under one chemistry and one index.

A full-universe row costs about 13 minutes at 491,836 candidates, plus two
minutes to build the position cache, so pass a short floor list for that mode.

See docs/validation/retention_changes_no_delivered_panel_2026-09-17.md.
"""

import json, pathlib, sys, time
from types import SimpleNamespace

DESIGN = pathlib.Path(sys.argv[1])
LABEL = sys.argv[2]
FLOORS = [
    float(x)
    for x in (sys.argv[3].split(",") if len(sys.argv) > 3 else "40,60,80,100,140".split(","))
]
SIZE = 12
# Held equal across every universe and floor, so only the universe varies.
SWAP_EVALS, SWAP_SECONDS, SCAN_WIDTH = 10_000, 10.0, 64

import pandas as pd
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.candidate_source import open_candidate_source
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered

params = json.loads((DESIGN / "params.json").read_text())
shortlist = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg, bg = params["fg_prefixes"], params["bg_prefixes"]
conditions = build_reaction_conditions(SimpleNamespace(**params))
_ensure_optimizers_registered()

started = time.time()
probe = open_candidate_source(str(DESIGN), conditions.fingerprint(), [12], frontier=len(shortlist))
universe = probe._provider._eligible_in_search_order()
print(f"{LABEL}: universe {len(universe)}, shortlist {len(shortlist)}", flush=True)

cache = PositionCache(fg + bg, universe, on_missing="warn")
print(f"  cache over {len(universe)} candidates in {time.time()-started:.0f}s", flush=True)

rows = []
for refills, tag in ((0, "shortlist"), (8, LABEL)):
    for floor in FLOORS:
        source = open_candidate_source(
            str(DESIGN), conditions.fingerprint(), [12], frontier=len(shortlist)
        )
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
                refinement_method="swap",
                swap_max_evaluations=SWAP_EVALS,
                swap_max_seconds=SWAP_SECONDS,
                objective_scan_width=SCAN_WIDTH,
                max_frontier_refills=refills,
                verbose=False,
            ),
            conditions=ReactionConditions(temp=30.0),
            polymerase="phi29",
        )
        t0 = time.time()
        plan = plan_pool(
            opt,
            source,
            [SIZE],
            [0.5],
            primer_length=12,
            coverage_metric="effective",
            min_selectivity_density=floor,
        )
        r = plan["rows"][0]
        row = dict(
            universe=tag,
            floor=floor,
            feasible=not r.get("failed_constraints"),
            density=r.get("selectivity_density"),
            coverage=r.get("coverage"),
            refills=r.get("frontier_refills"),
            ended=r.get("candidates_exhausted"),
            examined=plan["candidate_source"]["examined"],
            seconds=round(time.time() - t0, 1),
        )
        rows.append(row)
        print(
            f"  {tag:>10} floor {floor:>5}: feasible {str(row['feasible']):>5}  "
            f"density {row['density']:>8.3f}  cov {row['coverage']:.4f}  "
            f"refills {row['refills']}  examined {row['examined']:>7}  {row['seconds']:>5}s",
            flush=True,
        )
        if tag == "shortlist" and refills == 0:
            pass
    if refills == 0 and len(universe) == len(shortlist):
        print("  (universe equals shortlist; the refill rows would repeat)", flush=True)
        break

out = DESIGN.parent / f"floor_boundary_{LABEL}.json"
out.write_text(json.dumps(rows, indent=1))
print(f"\nwrote {out}")
