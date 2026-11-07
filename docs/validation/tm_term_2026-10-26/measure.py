#!/usr/bin/env python3
"""What the Stage-2 Tm term does to a delivered panel, against the noise floor.

Behind docs/validation/2026-10-26-wiring-the-tm-term.md.

Target is the real pcDNA plasmid shipped with the example (6,157 bp, one
record) and the pool is its real step-2 shortlist, so this measures a genuine
candidate pool rather than planted primers. No k-mer counter is needed: the
position indexes are written straight from the sequence, the way
`tests/test_hybrid_optimizer_run.py` does it.

REPEATS ARE THE POINT. `hybrid` Stage 2 is not deterministic. Two independent
sources, both pre-existing and neither introduced by the Tm term:

  * set and dict iteration order, which varies with PYTHONHASHSEED and breaks
    ties in the removal ranking;
  * `nx.algebraic_connectivity(..., method="tracemin_lu")` is called with no
    `seed`, so its random starting basis varies per call.

A single run per weight therefore cannot distinguish the term's effect from
the noise, and an earlier version of this script was fooled exactly that way.
Each weight is run `REPEATS` times and reported as a median with its range.

Two pools, because the term behaves very differently on each:

  single-length  k=10 only. Tm varies with GC alone, over a few degrees.
  mixed-length   k=9..12. One k step moves an effective Tm by 4-7 C, so the
                 term is mostly a length preference. This is the shape the
                 shipped min_k=6/max_k=12 range produces.

Run:
  PYTHONPATH=. python docs/validation/tm_term_2026-10-26/measure.py
"""

from __future__ import annotations

import json
import os
import random
import statistics
import tempfile

import h5py
import numpy as np
import pandas as pd
from Bio import SeqIO

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.selection_weights import OPTIMAL_TM_OFFSET_C
from neoswga.core.thermodynamics import reverse_complement

EXAMPLE = "examples/plasmid_example"
SHORTLIST = f"{EXAMPLE}/step2_df.csv.original"
FASTA = f"{EXAMPLE}/pcDNA.fasta"

REACTION_TEMP = 30.0
POLYMERASE = "phi29"
COVERAGE_REACH = 500
TARGET_SIZE = 6
PER_K = 150
SEED = 20261026
REPEATS = 7

# 0.0 is the shipped default: `stage2_tm_weight` is None, which Stage 2 reads
# as off. The rest are the `tm_weight` values the five `--application`
# profiles carry. They are measured here because a user would reach for them,
# and because whether the profile should drive this term is the question this
# measurement settles.
WEIGHTS = [0.0, 0.10, 0.25, 0.30, 1.0]

POOLS = {"single-length k=10": [10], "mixed-length k=9-12": [9, 10, 11, 12]}


def load_sequence() -> str:
    records = list(SeqIO.parse(FASTA, "fasta"))
    assert len(records) == 1, f"expected one record, got {len(records)}"
    return str(records[0].seq).upper()


def candidates_for(ks: list[int]) -> list[str]:
    frame = pd.read_csv(SHORTLIST)
    frame["k"] = frame["primer"].astype(str).str.len()
    rng = random.Random(SEED)
    out: list[str] = []
    for k in ks:
        pool = sorted(frame.loc[frame.k == k, "primer"].astype(str).tolist())
        out.extend(pool if len(pool) <= PER_K else rng.sample(pool, PER_K))
    return sorted(set(out))


def write_index(prefix: str, sequence: str, primers: list[str]) -> None:
    """One HDF5 per k, in the layout `PositionCache` expects."""
    by_k: dict[int, list[str]] = {}
    for primer in primers:
        by_k.setdefault(len(primer), []).append(primer)
    for k, group in by_k.items():
        with h5py.File(f"{prefix}_{k}mer_positions.h5", "w") as handle:
            for primer in group:
                for key in {primer, reverse_complement(primer)}:
                    positions, i = [], sequence.find(key)
                    while i != -1:
                        positions.append(i)
                        i = sequence.find(key, i + 1)
                    # An absent entry and an empty one are different things, so
                    # the empty one is written rather than skipped.
                    handle.create_dataset(key, data=np.array(positions, dtype=np.int64))


def mean_abs_deviation(primers, conditions) -> float:
    """The quantity the Tm term exists to reduce."""
    optimal = REACTION_TEMP + OPTIMAL_TM_OFFSET_C
    tms = [conditions.calculate_effective_tm(p) for p in primers]
    return sum(abs(t - optimal) for t in tms) / len(tms)


def tm_spread(primers, conditions) -> float:
    tms = [conditions.calculate_effective_tm(p) for p in primers]
    return max(tms) - min(tms)


def one_run(primers, prefix, sequence, conditions, weight):
    optimizer = HybridOptimizer(
        position_cache=PositionCache([prefix], primers),
        fg_prefixes=[prefix],
        fg_seq_lengths=[len(sequence)],
        bg_prefixes=[],
        bg_seq_lengths=[],
        coverage_reach=COVERAGE_REACH,
        conditions=conditions,
        reaction_temp=REACTION_TEMP,
        # None is the shipped default and the only value that leaves Stage 2
        # bit-identical to how it ran before the term existed.
        stage2_tm_weight=weight if weight > 0 else None,
    )
    result = optimizer.optimize(primers, final_count=TARGET_SIZE, verbose=False)
    delivered = tuple(sorted(result.primers))
    return {
        "primers": list(delivered),
        "mean_abs_deviation": mean_abs_deviation(delivered, conditions),
        "tm_spread": tm_spread(delivered, conditions),
        "coverage": float(result.final_coverage or 0.0),
        "lengths": sorted({len(p) for p in delivered}),
    }


def summarise(label, weight, runs, baseline_panels):
    devs = [r["mean_abs_deviation"] for r in runs]
    covs = [r["coverage"] for r in runs]
    panels = {tuple(r["primers"]) for r in runs}
    overlap = []
    for r in runs:
        best = max(
            (len(set(r["primers"]) & set(b)) / len(set(r["primers"]) | set(b)))
            for b in baseline_panels
        )
        overlap.append(best)
    return {
        "pool": label,
        "stage2_tm_weight": weight,
        "repeats": len(runs),
        "distinct_panels": len(panels),
        "mean_abs_deviation_median": round(statistics.median(devs), 3),
        "mean_abs_deviation_min": round(min(devs), 3),
        "mean_abs_deviation_max": round(max(devs), 3),
        "coverage_median": round(statistics.median(covs), 4),
        "coverage_min": round(min(covs), 4),
        "coverage_max": round(max(covs), 4),
        "tm_spread_median": round(statistics.median(r["tm_spread"] for r in runs), 2),
        "best_jaccard_to_a_baseline_panel_median": round(statistics.median(overlap), 3),
        "lengths_seen": sorted({tuple(r["lengths"]) for r in runs}),
        "runs": runs,
    }


def main() -> int:
    sequence = load_sequence()
    rows = []
    for label, ks in POOLS.items():
        primers = candidates_for(ks)
        directory = tempfile.mkdtemp(prefix="tmterm_")
        prefix = os.path.join(directory, "pcDNA")
        write_index(prefix, sequence, primers)
        conditions = ReactionConditions(temp=REACTION_TEMP, polymerase=POLYMERASE)

        by_weight = {
            w: [one_run(primers, prefix, sequence, conditions, w) for _ in range(REPEATS)]
            for w in WEIGHTS
        }
        baseline_panels = [tuple(r["primers"]) for r in by_weight[0.0]]

        print(
            f"\n{label}  ({len(primers)} candidates, target {TARGET_SIZE}, "
            f"{POLYMERASE} {REACTION_TEMP} C, optimum "
            f"{REACTION_TEMP + OPTIMAL_TM_OFFSET_C} C, {REPEATS} repeats)\n"
        )
        header = (
            f"{'weight':>8}{'panels':>8}{'mean |dTm| (med)':>19}{'range':>14}"
            f"{'coverage (med)':>16}{'range':>16}{'jaccard':>9}"
        )
        print(header)
        print("-" * len(header))
        for w in WEIGHTS:
            s = summarise(label, w, by_weight[w], baseline_panels)
            rows.append(s)
            dev_range = f"{s['mean_abs_deviation_min']:.2f}-{s['mean_abs_deviation_max']:.2f}"
            cov_range = f"{s['coverage_min']:.4f}-{s['coverage_max']:.4f}"
            print(
                f"{w:>8.2f}{s['distinct_panels']:>8}"
                f"{s['mean_abs_deviation_median']:>19.3f}{dev_range:>14}"
                f"{s['coverage_median']:>16.4f}{cov_range:>16}"
                f"{s['best_jaccard_to_a_baseline_panel_median']:>9.3f}"
            )

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results.json")
    with open(out, "w") as handle:
        json.dump(rows, handle, indent=2, default=str)
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
