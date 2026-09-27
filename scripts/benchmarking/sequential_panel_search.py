"""SUPERSEDED. Use `equal_allowance_comparison.py` instead.

Kept because `docs/validation/sequential_panel_service_2026-09-19.md` quotes
figures this script produced, and deleting it would leave that record
unreproducible. Do not draw new conclusions from it.

Why it is superseded, in its own original words: "Budgets are per stage, not
equal total compute". An arm running more stages was handed more compute, so a
difference in its panel could not be attributed to the search. It also ran one
unseeded configuration, and it built one proposal OUTSIDE the timed region and
shared it between arms, which understates the cost of the arm that needed it.

Measured consequence, from
docs/validation/equal_allowance_comparison_2026-09-27.md: under a SHARED
allowance the reduction stage can be starved by refinement before it runs at
all, and at 400 evaluations it delivers 12 primers where at 8,000 it delivers 3.
Reading the first as a fact about reduction would be reading the allowance.
"""

import argparse
import json
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from neoswga.core.design_context import design_context_from_params
from neoswga.core.optimization_service import (
    OptimizationRequest,
    panel_violations,
    run_panel_search,
)
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.panel_refinement import objective_for_optimizer
from neoswga.core.pool_objective import PoolConstraints
from neoswga.core.position_cache import PositionCache
from neoswga.core.unified_optimizer import _ensure_optimizers_registered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--rebuild-indexes", action="store_true")
    parser.add_argument("--sizes", type=int, nargs="+", default=[6, 12, 24])
    parser.add_argument("--target", type=float, default=0.7)
    parser.add_argument("--density", type=float, default=20)
    parser.add_argument("--evaluations", type=int, default=1000)
    args = parser.parse_args()
    params = json.loads((args.design / "params.json").read_text())
    context = design_context_from_params(params)
    data_dir = Path(params.get("data_dir", str(args.design)))
    candidates = pd.read_csv(data_dir / "step3_df.csv")["primer"].tolist()
    if args.rebuild_indexes:
        from neoswga.core.string_search import get_positions

        index_dir = args.output.parent / "indexes"
        index_dir.mkdir(parents=True, exist_ok=True)
        for group in ("fg", "bg"):
            prefixes = [
                str((index_dir / f"{group}_{i}").resolve())
                for i in range(len(params[f"{group}_genomes"]))
            ]
            get_positions(
                candidates,
                prefixes,
                params[f"{group}_genomes"],
                circular=params.get(f"{group}_circular", False),
                k_values=sorted({len(p) for p in candidates}),
            )
            params[f"{group}_prefixes"] = prefixes
    cache = PositionCache(
        params["fg_prefixes"] + params["bg_prefixes"], candidates, on_missing="error"
    )
    cache.require_record_metadata(params["fg_prefixes"] + params["bg_prefixes"])
    _ensure_optimizers_registered()
    config = context.optimizer_config(
        verbose=False,
        swap_max_evaluations=args.evaluations,
        swap_max_seconds=30,
        beam_max_evaluations=args.evaluations,
    )
    constraints = PoolConstraints(min_selectivity_density=args.density)
    rows = []
    for size in args.sizes:

        def optimizer():
            return OptimizerFactory.create(
                name="dominating-set",
                position_cache=cache,
                fg_prefixes=params["fg_prefixes"],
                fg_seq_lengths=params["fg_seq_lengths"],
                bg_prefixes=params["bg_prefixes"],
                bg_seq_lengths=params["bg_seq_lengths"],
                config=config,
                conditions=context.conditions,
            )

        proposal_opt = optimizer()
        objective_for_optimizer(proposal_opt, constraints)
        proposal = proposal_opt.optimize(candidates, target_size=size)
        for mode in ("single_pass", "sequential_reduction"):
            opt = optimizer()
            started = time.monotonic()
            result = run_panel_search(
                OptimizationRequest(
                    opt,
                    tuple(candidates),
                    size,
                    constraints=constraints,
                    minimize=mode == "sequential_reduction",
                    target_coverage=args.target,
                ),
                initial_result=proposal,
            )
            row = {
                "requested_size": size,
                "mode": mode,
                "seconds": time.monotonic() - started,
                "size": len(result.primers),
                "effective_coverage": result.metrics.effective_fg_coverage,
                "raw_coverage": result.metrics.fg_coverage,
                "density": result.metrics.selectivity_density,
                "background_sites": result.metrics.total_bg_sites,
                "violations": list(panel_violations(opt, result.primers)),
                "primers": list(result.primers),
                "stages": list(result.stage_history),
            }
            rows.append(row)
            print({k: v for k, v in row.items() if k not in {"primers", "stages"}}, flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "design": str(args.design.resolve()),
                "candidate_count": len(candidates),
                "conditions_fingerprint": context.conditions.fingerprint(),
                "constraints": asdict(constraints),
                "target": args.target,
                "evaluations_per_stage": args.evaluations,
                "interpretation": "Exploratory comparison; per-stage budgets are equal, total compute is not. Coverage is predicted.",
                "rows": rows,
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
