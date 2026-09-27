"""A density-only ceiling, which bounds nothing deliverable on its own.

READ THIS BEFORE QUOTING THE NUMBER. The ceiling here ignores the coverage
target and the dimer screen, and on the measured pool a panel achieving it has
coverage 0.4042 against a 0.5 target and 30 dimerising pairs out of 66. Treating
it as headroom cost real effort: three Stage 1 search rules were built chasing a
gap that did not exist. See
docs/validation/stage_one_constraint_awareness_2026-09-18.md.

It is kept because the arithmetic is exact and useful for the narrow question it
answers -- does retaining more candidates raise the density a pool can reach at
all, which it does not -- and because the correction is worth having on record.

The highest selectivity density any N-primer panel from a pool can reach.

Exact, not a search result. `occupancy.weighted_site_load` is a sum of
per-primer terms, so a panel's foreground and background loads are additive and

    density(panel) >= D   <=>   SUM_i (f_i * bg_len - D * b_i * fg_len) >= 0

is linear in the per-candidate loads. The best N-subset for a given D is
therefore the N largest of those terms, and the achievable ceiling is found by
binary search on D.

This decides whether Phase 5's negative result is about the search or about the
pool. If the ceiling is near the floor the shortlist already reaches, no search
improvement can do better and Stage 1 is not the problem.

Usage: achievable_density.py <design_dir> [sizes]
"""

import json, pathlib, sys, time
from types import SimpleNamespace

DESIGN = pathlib.Path(sys.argv[1])
SIZES = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["12"])]

import pandas as pd
from neoswga.core.candidate_source import open_candidate_source
from neoswga.core.occupancy import default_mismatch_penalty, mismatch_tm, site_occupancy
from neoswga.core.mismatch_counts import mismatch_class_counts
from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions
from neoswga.core.thermodynamics import calculate_enthalpy_entropy

params = json.loads((DESIGN / "params.json").read_text())
shortlist = list(dict.fromkeys(pd.read_csv(DESIGN / "step3_df.csv")["primer"]))
fg, bg = params["fg_prefixes"], params["bg_prefixes"]
fg_len = float(sum(params["fg_seq_lengths"]))
bg_len = float(sum(params["bg_seq_lengths"]))
conditions = ReactionConditions(temp=30.0, polymerase="phi29")
fingerprint = build_reaction_conditions(SimpleNamespace(**params)).fingerprint()

source = open_candidate_source(str(DESIGN), fingerprint, [12], frontier=len(shortlist))
universe = source._provider._eligible_in_search_order()
penalty = default_mismatch_penalty()
mm = int(params.get("max_mismatches", 1) or 0)


def load(primer, prefixes):
    dh, _ = calculate_enthalpy_entropy(primer)
    tm = conditions.calculate_effective_tm(primer)
    total = 0.0
    for distance, count in mismatch_class_counts(primer, prefixes, mm).items():
        if count:
            total += count * site_occupancy(dh, mismatch_tm(tm, distance, penalty), conditions.temp)
    return total


def profile(pool, name):
    started = time.time()
    rows = []
    for primer in pool:
        rows.append((primer, load(primer, fg), load(primer, bg)))
    print(
        f"  {name}: per-candidate loads for {len(pool)} in {time.time()-started:.0f}s", flush=True
    )
    return rows


def ceiling(rows, n):
    """Highest D for which some n-subset satisfies the density floor."""
    lo, hi = 0.0, 1e7
    for _ in range(200):
        mid = (lo + hi) / 2
        slack = sorted((f * bg_len - mid * b * fg_len for _, f, b in rows), reverse=True)
        if sum(slack[:n]) >= 0:
            lo = mid
        else:
            hi = mid
    return lo


pools = [("shortlist", shortlist), ("universe", universe)]
out = {}
for name, pool in pools:
    rows = profile(pool, name)
    for n in SIZES:
        c = ceiling(rows, n)
        # The panel that achieves it, and its own coverage is a separate matter.
        best = sorted(rows, key=lambda r: r[1] * bg_len - c * r[2] * fg_len, reverse=True)[:n]
        fsum = sum(r[1] for r in best)
        bsum = sum(r[2] for r in best)
        actual = (fsum / fg_len) / (bsum / bg_len) if bsum else float("inf")
        print(
            f"    n={n:<3} achievable density ceiling: {c:10.3f}   "
            f"(check: the argmax subset gives {actual:.3f})",
            flush=True,
        )
        out[f"{name}_n{n}"] = dict(ceiling=c, check=actual, pool=len(pool))

(DESIGN.parent / "achievable_density.json").write_text(json.dumps(out, indent=1))
print(f"\nwrote {DESIGN.parent / 'achievable_density.json'}")
