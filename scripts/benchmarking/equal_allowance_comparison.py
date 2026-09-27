"""Compare panel searches under EQUAL declared allowances, across fixed seeds.

This replaces `sequential_panel_search.py`, whose own docstring conceded the
flaw: "Budgets are per stage, not equal total compute". A per-stage budget means
an arm that runs more stages is given more compute, so a difference in the
delivered panel cannot be attributed to the search rather than to the allowance
it was handed. It also ran one unseeded configuration, so a single favourable
panel could not be told from a tendency.

What is equal here, and what is not:

  EQUAL   one shared `SearchBudget` allowance per arm, the same integer for
          every arm and every seed, spent from the moment the arm begins. The
          references, the candidate inventory, the constraints, the reaction
          and the requested size are identical across arms.
  NOT     wall-clock time, which is an outcome rather than an input, and
          reported for that reason.

What the evaluation count covers: uncached evaluations of the shared objective,
plus the `clique` method's own scoring loop since 2026-09-27. It does not cover
one final assessment per stage, which is deliberately uncharged so reporting
cannot consume a search's allowance. `stop_reason` says whether an arm spent
the allowance or finished inside it -- an arm that finished is not being
compared on equal compute, it is being compared on enough compute, and the two
readings are different.

**The position cache is shared across arms, deliberately, and that is a warm
cache.** Building it over the Wolbachia inventory costs about two minutes, and
paying that per arm would dominate every timing while changing no answer: it
holds binding POSITIONS, which are reference data identical for every arm, not
evaluator results. Each arm gets a fresh optimizer from the factory, so any
per-optimizer memoisation starts cold. Proposal generation is INSIDE the timed
and charged region, unlike the script this replaces, which built one proposal
outside both and shared it.

Coverage here is a geometric proxy. Nothing in this repository has measured
sequencing recovery; see docs/validation/design_release_gates.md.

Usage:
    python scripts/benchmarking/equal_allowance_comparison.py \
        --design examples/wolbachia_pool_design \
        --output runs/equal_allowance.json \
        --sizes 6 12 --seeds 1 2 3 --evaluations 400
"""

import argparse
import json
import random
import resource
import statistics
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from neoswga.core.design_context import design_context_from_params  # noqa: E402
from neoswga.core.optimization_service import (  # noqa: E402
    OptimizationRequest,
    panel_violations,
    run_panel_search,
)
from neoswga.core.optimizer_factory import OptimizerFactory  # noqa: E402
from neoswga.core.panel_refinement import objective_for_optimizer  # noqa: E402
from neoswga.core.pool_objective import PoolConstraints  # noqa: E402
from neoswga.core.position_cache import PositionCache  # noqa: E402
from neoswga.core.search_control import SearchBudget  # noqa: E402
from neoswga.core.unified_optimizer import _ensure_optimizers_registered  # noqa: E402

_RSS_SCALE = 1 if sys.platform == "darwin" else 1024

#: The arms. Each is (label, method, minimize), and each receives the identical
#: allowance. `minimize` is the sequential reduction the replaced script
#: compared; the methods are the two that share a Stage-1 greedy, so a
#: difference between them is the refinement rather than the proposal.
ARMS = (
    ("single_pass", "dominating-set", False),
    ("sequential_reduction", "dominating-set", True),
    ("hybrid_single_pass", "hybrid", False),
)


def _peak_rss() -> int:
    """Process high-water mark, which NEVER FALLS.

    Reported once for the whole run and deliberately not per arm. A per-arm
    delta looked informative and was not: the first arm absorbs the peak of
    building the shared position cache, and every later arm then shows a small
    delta because the high-water mark is already set. Measured on a smoke run,
    the first arm read 1,112 MB and the next 54 MB for the same work. A real
    per-arm memory figure needs one process per arm, which is how
    `count_coverage_rss.py` in this directory does it.
    """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * _RSS_SCALE


def _seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def run_arm(
    label,
    method,
    minimize,
    *,
    cache,
    params,
    context,
    config,
    constraints,
    size,
    seed,
    evaluations,
    candidates,
    target_coverage,
):
    """One arm, one seed, one size, under its own fresh allowance."""
    _seed_everything(seed)
    budget = SearchBudget(max_evaluations=evaluations)
    started = time.monotonic()

    optimizer = OptimizerFactory.create(
        name=method,
        position_cache=cache,
        fg_prefixes=params["fg_prefixes"],
        fg_seq_lengths=params["fg_seq_lengths"],
        bg_prefixes=params["bg_prefixes"],
        bg_seq_lengths=params["bg_seq_lengths"],
        config=config,
        conditions=context.conditions,
    )
    objective_for_optimizer(optimizer, constraints)

    failed = None
    try:
        result = run_panel_search(
            OptimizationRequest(
                optimizer,
                tuple(candidates),
                size,
                constraints=constraints,
                minimize=minimize,
                budget=budget,
                # Passed rather than left to the dataclass default of 0.7.
                # `--target` was accepted and recorded in the output while
                # reaching nothing, so every arm of the first sweep ran at the
                # default and the JSON said otherwise. That is the inert-option
                # defect this repository ratchets against, in the harness
                # written to measure it.
                target_coverage=target_coverage,
            )
        )
    except Exception as exc:  # noqa: BLE001 - a failure is a result to report
        return {
            "arm": label,
            "method": method,
            "seed": seed,
            "requested_size": size,
            "failed": f"{type(exc).__name__}: {exc}",
            "seconds": time.monotonic() - started,
        }

    elapsed = time.monotonic() - started
    spent = budget.describe()
    return {
        "arm": label,
        "method": method,
        "seed": seed,
        "requested_size": size,
        "failed": failed,
        "seconds": elapsed,
        "evaluations_spent": spent.get("evaluations"),
        "evaluations_allowed": spent.get("max_evaluations"),
        "stop_reason": spent.get("stop_reason"),
        "uncounted_scopes": list(spent.get("uncounted_scopes") or ()),
        "delivered_size": len(result.primers),
        "effective_coverage": result.metrics.effective_fg_coverage,
        "raw_coverage": result.metrics.fg_coverage,
        "density": result.metrics.selectivity_density,
        "background_sites": result.metrics.total_bg_sites,
        "violations": list(panel_violations(optimizer, result.primers)),
        "candidates_examined": len(candidates),
        # Which stages ran and what each changed. Without this a sweep cannot
        # tell a stage that ran and found nothing from a stage that never ran,
        # and those support opposite conclusions about the arm: the first says
        # there was nothing to find, the second says the arm was not exercised.
        "stages": [
            (
                {k: v for k, v in dict(stage).items() if k != "primers"}
                if isinstance(stage, dict)
                else str(stage)
            )
            for stage in (result.stage_history or [])
        ],
        "primers": sorted(result.primers),
    }


def summarise(rows):
    """Median and range per (size, arm), because one panel is not a tendency."""
    summary = []
    keys = sorted({(r["requested_size"], r["arm"]) for r in rows})
    for size, arm in keys:
        group = [r for r in rows if r["requested_size"] == size and r["arm"] == arm]
        usable = [r for r in group if not r.get("failed")]
        entry = {
            "requested_size": size,
            "arm": arm,
            "seeds": len(group),
            "failures": len(group) - len(usable),
        }
        for field in (
            "delivered_size",
            "effective_coverage",
            "density",
            "background_sites",
            "seconds",
            "evaluations_spent",
        ):
            values = [r[field] for r in usable if r.get(field) is not None]
            if not values:
                continue
            entry[f"{field}_median"] = statistics.median(values)
            entry[f"{field}_min"] = min(values)
            entry[f"{field}_max"] = max(values)
        # Identical panels across seeds mean the arm is deterministic here, and
        # a median over one distinct answer is not an uncertainty estimate.
        entry["distinct_panels"] = len({tuple(r["primers"]) for r in usable})
        entry["stop_reasons"] = sorted({str(r.get("stop_reason")) for r in usable})
        summary.append(entry)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sizes", type=int, nargs="+", default=[6, 12])
    parser.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument(
        "--evaluations",
        type=int,
        required=True,
        help="the shared total allowance every arm receives",
    )
    parser.add_argument("--density", type=float, default=20.0)
    parser.add_argument("--target", type=float, default=0.7)
    parser.add_argument(
        "--candidates", type=int, default=None, help="cap the inventory, for a quicker smoke run"
    )
    args = parser.parse_args()

    params = json.loads((args.design / "params.json").read_text())
    context = design_context_from_params(params)
    work = Path(params.get("data_dir", str(args.design)))
    candidates = pd.read_csv(work / "step3_df.csv")["primer"].tolist()
    if args.candidates:
        candidates = candidates[: args.candidates]

    prefixes = params["fg_prefixes"] + params["bg_prefixes"]
    cache_started = time.monotonic()
    cache = PositionCache(prefixes, candidates, on_missing="error")
    cache.require_record_metadata(prefixes)
    cache_seconds = time.monotonic() - cache_started

    _ensure_optimizers_registered()
    config = context.optimizer_config(verbose=False)
    constraints = PoolConstraints(min_selectivity_density=args.density)

    rows = []
    for size in args.sizes:
        for seed in args.seeds:
            for label, method, minimize in ARMS:
                row = run_arm(
                    label,
                    method,
                    minimize,
                    cache=cache,
                    params=params,
                    context=context,
                    config=config,
                    constraints=constraints,
                    size=size,
                    seed=seed,
                    evaluations=args.evaluations,
                    candidates=candidates,
                    target_coverage=args.target,
                )
                rows.append(row)
                print(
                    {k: v for k, v in row.items() if k not in {"primers", "uncounted_scopes"}},
                    flush=True,
                )

    payload = {
        "design": str(args.design.resolve()),
        "candidate_count": len(candidates),
        "conditions_fingerprint": context.conditions.fingerprint(),
        "constraints": asdict(constraints),
        "coverage_target": args.target,
        "target_coverage_reaches_the_search": True,
        "allowance_per_arm": args.evaluations,
        "seeds": list(args.seeds),
        "position_cache_seconds": cache_seconds,
        "process_peak_rss_mb": _peak_rss() / 1e6,
        "memory_note": (
            "One figure for the whole process. ru_maxrss is a high-water mark "
            "that never falls, so a per-arm delta credits the first arm with "
            "the shared cache build and reports near-zero for the rest. A "
            "per-arm figure needs one process per arm."
        ),
        "what_the_allowance_covers": (
            "uncached evaluations of the shared objective, plus the clique "
            "method's own scoring loop. One final assessment per stage is "
            "deliberately uncharged."
        ),
        "shared_across_arms": (
            "the position cache, which holds binding positions -- reference "
            "data identical for every arm. Each arm builds a fresh optimizer "
            "and generates its own proposal inside the timed region."
        ),
        "interpretation": (
            "Coverage is a geometric proxy at the declared reach and has not "
            "been calibrated against sequencing breadth. An arm whose "
            "stop_reason is null finished inside its allowance, so it was not "
            "constrained by it; comparing such arms on runtime compares "
            "implementations, not search quality under a budget."
        ),
        "rows": rows,
        "summary": summarise(rows),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False, default=str))
    print(f"\nwrote {args.output}")
    for entry in payload["summary"]:
        print(entry, flush=True)


if __name__ == "__main__":
    main()
