"""Does an interval sweep reproduce _compute_effective_coverage, and how fast?

The occupancy coverage loop is the dominant term in one objective evaluation
(docs/validation/objective_evaluation_cost_2026-09-17.md). It makes two passes
over the whole target per primer. The same quantity can be accumulated over
window endpoints instead, which is independent of the target length.

This is a PROTOTYPE for measurement, not a replacement. It does not pass
record_starts, matching the loop it is compared against, so neither confines a
window to the record holding its site. A real replacement has to decide that
question rather than inherit it.

Run with a design directory that has been through count-kmers, filter and score.
"""

import json, sys, time, pathlib, math
import numpy as np

DESIGN = pathlib.Path(
    sys.argv[1]
    if len(sys.argv) > 1
    else "examples/wolbachia_pool_design/objective_cost_2026-09-17/design"
)

import pandas as pd
from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.unified_optimizer import _ensure_optimizers_registered
from neoswga.core.coverage import _mark_window
from neoswga.core.occupancy import site_occupancy
from neoswga.core.thermodynamics import calculate_enthalpy_entropy

params = json.loads((DESIGN / "params.json").read_text())
candidates = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg = params["fg_prefixes"]
bg = params["bg_prefixes"]
LENGTH = params["fg_seq_lengths"][0]
REACH = 3000
cache = PositionCache(fg + bg, candidates)
_ensure_optimizers_registered()


def sweep(opt, positions_by_primer, total_length, circular):
    """Same quantity, accumulated over window endpoints instead of over bases."""
    events = {}
    for primer, positions in positions_by_primer.items():
        if not positions:
            continue
        tm = opt.conditions.calculate_effective_tm(primer)
        dh, _ = calculate_enthalpy_entropy(primer)
        theta = site_occupancy(dh, tm, total_length and opt.conditions.temp)
        if theta <= 0.0:
            continue
        weight = math.log1p(-theta) if theta < 1.0 else -math.inf
        # Merge this primer's windows first: a primer does not stack with
        # itself, which is the grouping the loop implements.
        spans = []
        for pos in positions:
            a, b = int(pos) - REACH, int(pos) + REACH
            if circular:
                if a < 0:
                    spans.append((0, min(b, total_length)))
                    spans.append((max(0, total_length + a), total_length))
                    continue
                if b > total_length:
                    spans.append((a, total_length))
                    spans.append((0, min(b - total_length, total_length)))
                    continue
            spans.append((max(0, a), min(total_length, b)))
        spans.sort()
        merged = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        for a, b in merged:
            if b <= a:
                continue
            events[a] = events.get(a, 0.0) + weight
            events[b] = events.get(b, 0.0) - weight

    if not events:
        return 0.0
    keys = sorted(events)
    running = 0.0
    covered = 0.0
    prev = keys[0]
    for x in keys:
        if x > prev and running != 0.0:
            covered += (x - prev) * (1.0 - math.exp(running))
        running += events[x]
        prev = x
    return covered / total_length


for circular in (False, True):
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
            extension_reach=REACH,
            fg_circular=circular,
            verbose=False,
        ),
        conditions=ReactionConditions(temp=30.0),
        polymerase="phi29",
    )
    for panel_size in (12, 24):
        panel = candidates[:panel_size]
        pbp = {p: [int(x) for x in cache.get_positions(fg[0], p, "both")] for p in panel}
        t0 = time.time()
        for _ in range(20):
            a = optimizer._compute_effective_coverage(pbp, LENGTH)
        t_loop = (time.time() - t0) / 20
        t0 = time.time()
        for _ in range(20):
            b = sweep(optimizer, pbp, LENGTH, circular)
        t_sweep = (time.time() - t0) / 20
        print(
            f"circular={str(circular):<5} n={panel_size:<3} "
            f"loop {t_loop*1e3:7.2f} ms  sweep {t_sweep*1e3:6.2f} ms  "
            f"speedup {t_loop/t_sweep:5.1f}x  "
            f"loop={a:.10f} sweep={b:.10f} absdiff={abs(a-b):.2e}"
        )
