"""How long does one objective evaluation take on a real design?

Plan Phase 4: this decides whether bounding the per-step scan is a prerequisite
or an optimisation. Above roughly 200 microseconds per call, scoring every
candidate at every greedy step is not affordable at 491,836 candidates.
"""

import json, pathlib, shutil, subprocess, sys, time

EX = pathlib.Path("examples/wolbachia_pool_design")
OUT = EX / "objective_cost_2026-09-17"
DESIGN = OUT / "design"


def build():
    if (DESIGN / "step3_df.csv").exists():
        return
    if DESIGN.exists():
        shutil.rmtree(DESIGN)
    DESIGN.mkdir(parents=True)
    base = json.loads((EX / "params.json").read_text())
    params = dict(base)
    params["data_dir"] = str(DESIGN.resolve())
    params["candidate_retention"] = "post_gini"
    for key in ("fg_prefixes", "bg_prefixes"):
        params[key] = [str(DESIGN.resolve() / pathlib.Path(p).name) for p in base[key]]
    (DESIGN / "params.json").write_text(json.dumps(params, indent=2))
    for step in ("count-kmers", "filter", "score"):
        started = time.time()
        with open(DESIGN / f"{step}.log", "w") as fh:
            rc = subprocess.run(
                [sys.executable, "-m", "neoswga.cli_unified", step, "-j", "params.json"],
                cwd=DESIGN,
                stdout=fh,
                stderr=subprocess.STDOUT,
            ).returncode
        assert rc == 0, f"{step} failed:\n{(DESIGN / f'{step}.log').read_text()[-1500:]}"
        print(f"{step}: {time.time() - started:.0f}s", flush=True)


def main():
    build()
    sys.path.insert(0, ".")
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

    started = time.time()
    cache = PositionCache(fg + bg, candidates)
    print(f"\nPositionCache over {len(candidates)} candidates: {time.time() - started:.1f}s")

    _ensure_optimizers_registered()
    optimizer = OptimizerFactory.create(
        "hybrid",
        cache,
        fg,
        params["fg_seq_lengths"],
        bg,
        params["bg_seq_lengths"],
        config=OptimizerConfig(
            max_dimer_bp=3, max_self_dimer_bp=4, extension_reach=3000, verbose=False
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    objective = PoolObjective(optimizer.compute_pool_metrics, PoolConstraints())

    panel = candidates[:24]
    objective.coverage(panel)  # warm

    for n in (200, 1000):
        started = time.time()
        for candidate in candidates[24 : 24 + n]:
            objective.coverage([*panel, candidate])
        elapsed = time.time() - started
        per_call = elapsed / n
        print(
            f"{n:>5} calls at panel size 24: {elapsed:6.1f}s  " f"=> {per_call * 1e6:8.0f} us/call"
        )

    print("\n--- what that means for an unbounded scan ---")
    universe = 491_836
    print(f"{'candidates':>12} {'one greedy step':>18} {'24-primer panel':>18}")
    for size in (2_000, 20_670, universe):
        step = size * per_call
        print(f"{size:>12,} {step:>15.0f}s {step * 24 / 3600:>15.1f}h")


main()
