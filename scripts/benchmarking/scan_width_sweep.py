"""How narrow can the objective-scored scan be before the panel moves?

Phase 4 increment 4's acceptance measurement. The plan asks for widths 16, 64,
256 and unbounded compared on the real pool, with the warning that a Jaccard
against unbounded below about 0.8 at width 64 means the prescreen is deciding
rather than filtering and the width must rise.

Usage: scan_width_sweep.py [density_floor] [design_dir]

A floor that nothing violates leaves the repair unattempted and every width
identical, which measures nothing, so the floor is an argument and defaults to
a value that binds on the bundled Wolbachia design. The design directory must
have been through count-kmers, filter and score.
"""

import json
import pathlib
import sys
import time

DESIGN = pathlib.Path(
    sys.argv[2]
    if len(sys.argv) > 2
    else "examples/wolbachia_pool_design/objective_cost_2026-09-17/design"
)
FLOOR = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
SIZES = [10, 12, 14]
WIDTHS = [16, 64, 256, None]

import pandas as pd

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.pool_planner import plan_pool
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered


def main():
    params = json.loads((DESIGN / "params.json").read_text())
    candidates = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
    fg, bg = params["fg_prefixes"], params["bg_prefixes"]
    cache = PositionCache(fg + bg, candidates)
    _ensure_optimizers_registered()

    print(f"pool {len(candidates)} candidates, density floor {FLOOR}, sizes {SIZES}\n")
    results = {}
    for width in WIDTHS:
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
                objective_scan_width=width,
                verbose=False,
            ),
            conditions=ReactionConditions(temp=30.0),
            polymerase="phi29",
        )
        started = time.time()
        plan = plan_pool(
            optimizer,
            candidates,
            SIZES,
            [0.5],
            primer_length=12,
            coverage_metric="effective",
            min_selectivity_density=FLOOR,
        )
        results[width] = (plan, time.time() - started)
        print(f"width {str(width):>9}: {time.time() - started:6.1f}s")

    reference, _ = results[None]
    ref_rows = {r["requested_size"]: r for r in reference["rows"]}

    print(
        f"\n{'width':>9} {'size':>5} {'jaccard':>8} {'coverage':>11} {'density':>9} "
        f"{'feasible':>9} {'obj evals':>10} {'pairs':>9} {'stop':>18}"
    )
    for width in WIDTHS:
        plan, elapsed = results[width]
        for row in plan["rows"]:
            size = row["requested_size"]
            mine = set(row.get("primers") or [])
            theirs = set(ref_rows[size].get("primers") or [])
            jaccard = len(mine & theirs) / len(mine | theirs) if (mine | theirs) else 1.0
            rep = row.get("repair") or {}
            print(
                f"{str(width):>9} {size:>5} {jaccard:>8.3f} {row.get('coverage'):>11.6f} "
                f"{row.get('selectivity_density'):>9.3f} "
                f"{str(not row.get('failed_constraints')):>9} "
                f"{str(rep.get('objective_evaluations')):>10} "
                f"{str(rep.get('pairs_considered')):>9} {str(rep.get('stop_reason')):>18}"
            )

    print("\nwall clock per width:")
    for width in WIDTHS:
        print(f"  {str(width):>9}: {results[width][1]:6.1f}s")


if __name__ == "__main__":
    main()
