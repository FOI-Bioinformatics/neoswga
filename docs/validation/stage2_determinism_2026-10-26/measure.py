#!/usr/bin/env python3
"""Stage-2 determinism before and after seeding the connectivity solver.

Behind docs/validation/2026-10-26-stage2-determinism.md.

Counts distinct delivered panels over N repeats of ONE fixed input, within a
single process, at a given PYTHONHASHSEED. Before the fix this was 2-3; after
it is 1, and the signature is identical across hash seeds.

The `--unseeded` flag forces `algebraic_connectivity` back to `seed=None` so the
pre-fix behaviour can be reproduced without checking out the old commit.

Run:
  PYTHONPATH=. python docs/validation/stage2_determinism_2026-10-26/measure.py
  PYTHONPATH=. python docs/validation/stage2_determinism_2026-10-26/measure.py --unseeded
"""

from __future__ import annotations

import hashlib
import os
import random
import sys
import tempfile

import h5py
import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO

if "--unseeded" in sys.argv:
    _orig = nx.algebraic_connectivity
    nx.algebraic_connectivity = lambda G, **kw: _orig(G, **{**kw, "seed": None})

from neoswga.core.hybrid_optimizer import HybridOptimizer  # noqa: E402
from neoswga.core.position_cache import PositionCache  # noqa: E402
from neoswga.core.reaction_conditions import ReactionConditions  # noqa: E402
from neoswga.core.thermodynamics import reverse_complement  # noqa: E402

EXAMPLE = "examples/plasmid_example"
REPEATS = 6


def build():
    seq = str(list(SeqIO.parse(f"{EXAMPLE}/pcDNA.fasta", "fasta"))[0].seq).upper()
    frame = pd.read_csv(f"{EXAMPLE}/step2_df.csv.original")
    frame["k"] = frame["primer"].astype(str).str.len()
    rng = random.Random(20261026)
    prim: list[str] = []
    for k in (9, 10, 11, 12):  # mixed-length: the pool that wobbled most
        pool = sorted(frame.loc[frame.k == k, "primer"].astype(str).tolist())
        prim.extend(pool if len(pool) <= 150 else rng.sample(pool, 150))
    prim = sorted(set(prim))
    directory = tempfile.mkdtemp(prefix="det_")
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


def main() -> int:
    prim, prefix, seqlen = build()
    cond = ReactionConditions(temp=30.0, polymerase="phi29")
    seeded = "--unseeded" not in sys.argv
    signatures = set()
    for r in range(REPEATS):
        optimizer = HybridOptimizer(
            position_cache=PositionCache([prefix], prim),
            fg_prefixes=[prefix],
            fg_seq_lengths=[seqlen],
            bg_prefixes=[],
            bg_seq_lengths=[],
            coverage_reach=500,
            conditions=cond,
            reaction_temp=30.0,
        )
        result = optimizer.optimize(prim, final_count=6, verbose=False)
        sig = hashlib.md5("|".join(sorted(result.primers)).encode()).hexdigest()[:8]
        signatures.add(sig)
        print(f"  run {r}: {sig}  coverage={result.final_coverage:.4f}")
    print(
        f"\n{'SEEDED' if seeded else 'UNSEEDED'}: "
        f"{len(signatures)} distinct panel(s) over {REPEATS} runs "
        f"at PYTHONHASHSEED={os.environ.get('PYTHONHASHSEED', 'unset')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
