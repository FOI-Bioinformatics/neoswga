"""Where does a plan-pool run spend its time?

Plan step 204 asks for incremental interval/occupancy updates only where
profiling warrants them. This is the profile.
"""

import cProfile
import io
import pstats
import random
import sys
import time

import h5py
import numpy as np

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.optimizer_factory import OptimizerFactory
from neoswga.core.unified_optimizer import _ensure_optimizers_registered
from neoswga.core.position_cache import PositionCache
from neoswga.core.pool_planner import plan_pool
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 2_000_000
N_CANDIDATES = 300
K = 12
OUT = sys.argv[1] if len(sys.argv) > 1 else "/tmp/profile_pool"


def build(tmp):
    rng = random.Random(4711)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    primers, seen = [], set()
    while len(primers) < N_CANDIDATES:
        p = "".join(rng.choice("ACGT") for _ in range(K))
        if p in seen or reverse_complement(p) in seen:
            continue
        seen.add(p)
        primers.append(p)
        for _ in range(rng.randint(3, 30)):
            pos = rng.randrange(0, GENOME_LENGTH - K)
            seq[pos : pos + K] = list(p)
    genome = "".join(seq)

    fg = f"{tmp}/target"
    bg = f"{tmp}/host"
    for prefix, subject in ((fg, genome), (bg, genome[::-1])):
        with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as f:
            for p in primers:
                for key in {p, reverse_complement(p)}:
                    hits, i = [], subject.find(key)
                    while i != -1:
                        hits.append(i)
                        i = subject.find(key, i + 1)
                    f.create_dataset(key, data=np.array(hits, dtype=np.int64))
    return primers, fg, bg


def main():
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        primers, fg, bg = build(tmp)
        cache = PositionCache([fg, bg], primers)
        _ensure_optimizers_registered()
        config = OptimizerConfig(
            max_dimer_bp=3,
            max_self_dimer_bp=4,
            extension_reach=3_000,
            refinement_method="swap",
            verbose=False,
        )
        opt = OptimizerFactory.create(
            "hybrid",
            cache,
            [fg],
            [GENOME_LENGTH],
            [bg],
            [GENOME_LENGTH],
            config=config,
            conditions=None,
            polymerase="phi29",
        )
        started = time.time()
        profiler = cProfile.Profile()
        profiler.enable()
        plan = plan_pool(
            opt,
            primers,
            range(4, 17, 4),
            [0.8],
            primer_length=K,
            coverage_metric="raw",
            min_selectivity_density=float(sys.argv[2]) if len(sys.argv) > 2 else 1.5,
        )
        profiler.disable()
        elapsed = time.time() - started

        print(f"plan_pool wall clock: {elapsed:.1f}s over {len(plan['rows'])} sizes")
        for r in plan["rows"]:
            rep = r.get("repair", {})
            print(
                f"  size {r['requested_size']:>3}: {r['seconds']:6.1f}s "
                f"eligible={r['eligible']} repair={rep.get('method')} "
                f"swaps={rep.get('swaps')} evals={rep.get('evaluations')} "
                f"beam={rep.get('beam')} density={r.get('selectivity_density')}"
            )
        buf = io.StringIO()
        stats = pstats.Stats(profiler, stream=buf)
        stats.sort_stats("cumulative").print_stats(
            "neoswga/core/(pool_planner|panel_beam|swap_refinement|base_optimizer|coverage)"
        )
        print(buf.getvalue())


main()
