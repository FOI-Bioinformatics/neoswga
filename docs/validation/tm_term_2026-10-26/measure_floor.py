#!/usr/bin/env python3
"""What the Stage-1 floor's metric change does, measured in isolation.

Behind the "Floor" section of docs/validation/2026-10-26-wiring-the-tm-term.md.

Why this is a SEPARATE script from measure.py, and why it works the way it
does. hybrid Stage 2 is not deterministic (unseeded `algebraic_connectivity`
and hash-ordered set iteration), so an A/B of "floor on" against "floor off"
compares two different loop outputs and measures the noise, not the floor. The
honest experiment takes the ONE refinement the loop actually produced and asks
what each floor rule does with that same pair -- the loop never re-runs, so the
nondeterminism cannot leak into the comparison.

For reproducibility of the printed numbers this script pins both noise sources:
run it as `PYTHONHASHSEED=0 python docs/validation/tm_term_2026-10-26/measure_floor.py`
and it seeds `algebraic_connectivity`. NOTE: production was seeded later the
same day (docs/validation/2026-10-26-stage2-determinism.md), so this script's
own seeding is now redundant with the shipped code; it is kept so the script
reproduces identically whether run against that commit or an earlier one.

Run:
  PYTHONHASHSEED=0 PYTHONPATH=. python docs/validation/tm_term_2026-10-26/measure_floor.py
"""

from __future__ import annotations

import functools
import json
import os
import random
import tempfile

import h5py
import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO

# Seed the one unseeded call so the loop output is reproducible for this
# measurement. Production is untouched; see the module docstring.
nx.algebraic_connectivity = functools.partial(nx.algebraic_connectivity, seed=0)

from neoswga.core import hybrid_optimizer as HO  # noqa: E402
from neoswga.core.position_cache import PositionCache  # noqa: E402
from neoswga.core.reaction_conditions import ReactionConditions  # noqa: E402
from neoswga.core.selection_weights import OPTIMAL_TM_OFFSET_C  # noqa: E402
from neoswga.core.stage2_ranking import prefix_floors_refinement  # noqa: E402
from neoswga.core.thermodynamics import reverse_complement  # noqa: E402

EXAMPLE = "examples/plasmid_example"
REACTION_TEMP = 30.0
OPTIMUM = REACTION_TEMP + OPTIMAL_TM_OFFSET_C
PER_K = 150
SEED = 20261026
POOLS = {"single-length k=10": [10], "mixed-length k=9-12": [9, 10, 11, 12]}
WEIGHTS = [0.0, 0.10, 0.25, 0.50, 1.0]


def build(ks):
    seq = str(list(SeqIO.parse(f"{EXAMPLE}/pcDNA.fasta", "fasta"))[0].seq).upper()
    frame = pd.read_csv(f"{EXAMPLE}/step2_df.csv.original")
    frame["k"] = frame["primer"].astype(str).str.len()
    rng = random.Random(SEED)
    prim: list[str] = []
    for k in ks:
        pool = sorted(frame.loc[frame.k == k, "primer"].astype(str).tolist())
        prim.extend(pool if len(pool) <= PER_K else rng.sample(pool, PER_K))
    prim = sorted(set(prim))
    directory = tempfile.mkdtemp(prefix="floor_")
    prefix = os.path.join(directory, "pcDNA")
    by_k: dict[int, list[str]] = {}
    for p in prim:
        by_k.setdefault(len(p), []).append(p)
    for k, group in by_k.items():
        with h5py.File(f"{prefix}_{k}mer_positions.h5", "w") as handle:
            for p in group:
                for key in {p, reverse_complement(p)}:
                    pos, i = [], seq.find(key)
                    while i != -1:
                        pos.append(i)
                        i = seq.find(key, i + 1)
                    handle.create_dataset(key, data=np.array(pos, dtype=np.int64))
    return prim, prefix, len(seq)


def mean_abs_dev(primers, cond):
    return sum(abs(cond.calculate_effective_tm(p) - OPTIMUM) for p in primers) / len(primers)


def run_pool(label, ks):
    prim, prefix, seqlen = build(ks)
    cond = ReactionConditions(temp=REACTION_TEMP, polymerase="phi29")
    rows = []

    captured = {}
    orig = HO.HybridOptimizer._apply_stage1_floor

    def spy(self, current, floor_set, target_count, tm_by_primer, verbose):
        prefix_set = list(floor_set)[:target_count] if floor_set else []
        captured["refined"] = list(current)
        captured["prefix"] = prefix_set
        captured["prefix_axes"] = self._floor_axes(prefix_set, tm_by_primer) if prefix_set else None
        captured["refined_axes"] = self._floor_axes(list(current), tm_by_primer)
        return orig(self, current, floor_set, target_count, tm_by_primer, verbose)

    HO.HybridOptimizer._apply_stage1_floor = spy
    try:
        for weight in WEIGHTS:
            captured.clear()
            optimizer = HO.HybridOptimizer(
                position_cache=PositionCache([prefix], prim),
                fg_prefixes=[prefix],
                fg_seq_lengths=[seqlen],
                bg_prefixes=[],
                bg_seq_lengths=[],
                coverage_reach=500,
                conditions=cond,
                reaction_temp=REACTION_TEMP,
                stage2_tm_weight=weight if weight > 0 else None,
            )
            optimizer.optimize(prim, final_count=6, verbose=False)

            refined, prefix_set = captured["refined"], captured["prefix"]
            cov_r = optimizer._calculate_coverage(refined)
            cov_p = optimizer._calculate_coverage(prefix_set)
            # The two floor rules, on the one refinement the loop produced.
            old_keeps_prefix = cov_p > cov_r + 1e-12
            new_keeps_prefix = prefix_floors_refinement(
                captured["prefix_axes"], captured["refined_axes"]
            )
            old_panel = prefix_set if old_keeps_prefix else refined
            new_panel = prefix_set if new_keeps_prefix else refined
            rows.append(
                {
                    "pool": label,
                    "stage2_tm_weight": weight,
                    "axes": len(captured["refined_axes"]),
                    "refined_cov": round(cov_r, 4),
                    "refined_dev": round(mean_abs_dev(refined, cond), 3),
                    "prefix_cov": round(cov_p, 4),
                    "prefix_dev": round(mean_abs_dev(prefix_set, cond), 3),
                    "old_rule": "prefix" if old_keeps_prefix else "refined",
                    "new_rule": "prefix" if new_keeps_prefix else "refined",
                    "old_delivered_cov": round(optimizer._calculate_coverage(old_panel), 4),
                    "old_delivered_dev": round(mean_abs_dev(old_panel, cond), 3),
                    "new_delivered_cov": round(optimizer._calculate_coverage(new_panel), 4),
                    "new_delivered_dev": round(mean_abs_dev(new_panel, cond), 3),
                }
            )
    finally:
        HO.HybridOptimizer._apply_stage1_floor = orig
    return rows


def main() -> int:
    if os.environ.get("PYTHONHASHSEED") != "0":
        print(
            "WARNING: run with PYTHONHASHSEED=0 for the printed numbers to match "
            "the record; set-iteration order is otherwise unpinned.\n"
        )
    all_rows = []
    for label, ks in POOLS.items():
        rows = run_pool(label, ks)
        all_rows.extend(rows)
        print(f"\n{label}  (phi29 {REACTION_TEMP} C, optimum {OPTIMUM} C, target 6)\n")
        header = (
            f"{'weight':>7}{'axes':>5}{'refined cov/dev':>18}{'prefix cov/dev':>17}"
            f"{'old->':>9}{'new->':>9}{'old deliv dev':>15}{'new deliv dev':>15}"
        )
        print(header)
        print("-" * len(header))
        for r in rows:
            refined = f"{r['refined_cov']:.2f}/{r['refined_dev']:.2f}"
            prefix = f"{r['prefix_cov']:.2f}/{r['prefix_dev']:.2f}"
            print(
                f"{r['stage2_tm_weight']:>7.2f}{r['axes']:>5}{refined:>18}{prefix:>17}"
                f"{r['old_rule']:>9}{r['new_rule']:>9}"
                f"{r['old_delivered_dev']:>15.3f}{r['new_delivered_dev']:>15.3f}"
            )

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "floor_results.json")
    with open(out, "w") as handle:
        json.dump(all_rows, handle, indent=2, default=str)
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
