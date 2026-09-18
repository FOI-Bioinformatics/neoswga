"""What each dimer policy costs, decides, and admits.

Produced the measurements that decided the design of the optional stability
floor (`max_dimer_dg`), recorded in
`docs/validation/dimer_stability_floor_2026-09-18.md`. Needs no design
directory and no external tool.

Usage:
    python scripts/benchmarking/dimer_policy_comparison.py

Three questions:

1. **Panel size.** The run screen is what bounds it: CLAUDE.md records the

   shipped pools supporting 29, 31 and 26 primers at `max_dimer_bp` 3.
   `--allow-dimer-relaxation` is the current escape hatch and it is blunt: it
   admits violating pairs with a warning.
2. **Cost**, per pair, against the run screen.
3. **What each policy ADMITS**, in complementary run length. This is the
   safety-relevant number, because a stability floor bounds stability and not
   length, and this project delivered an 11 bp heterodimer once.

A free-energy floor is a principled alternative to loosening the run length,
and the third measurement is why it is offered as an ADDITIONAL floor rather
than a replacement.
"""

import itertools
import random

from neoswga.core.dimer import is_dimer, is_dimer_thermodynamic
from neoswga.core.reaction_conditions import ReactionConditions

rng = random.Random(20260918)
K, N = 12, 200
COND = ReactionConditions(temp=30.0, polymerase="phi29")


def pool(gc):
    """A GC-matched pool, reseeded per call.

    Reseeded because three sections sample and a shared generator would give
    each one a different pool, so the same policy would report different
    figures in different tables of one run.
    """
    local = random.Random((20260918, round(gc, 4)).__hash__())
    p_gc, p_at = gc / 2, (1 - gc) / 2
    return ["".join(local.choices("ACGT", weights=[p_at, p_gc, p_gc, p_at], k=K)) for _ in range(N)]


def greedy_compatible(primers, conflicts):
    """Largest dimer-free subset a greedy finds, same shape as the optimizers'."""
    chosen = []
    for p in primers:
        if all((p, q) not in conflicts and (q, p) not in conflicts for q in chosen):
            chosen.append(p)
    return chosen


def conflict_set(primers, *, bp=None, dg=None):
    out = set()
    for a, b in itertools.combinations(primers, 2):
        bad = False
        if bp is not None and is_dimer(a, b, max_dimer_bp=bp):
            bad = True
        if (
            not bad
            and dg is not None
            and is_dimer_thermodynamic(a, b, delta_g_threshold=dg, conditions=COND)
        ):
            bad = True
        if bad:
            out.add((a, b))
    return out


POLICIES = [
    ("run <= 3 (shipped default)", dict(bp=3)),
    ("run <= 4", dict(bp=4)),
    ("run <= 5", dict(bp=5)),
    ("dG > -6 only", dict(dg=-6.0)),
    ("run <= 5 and dG > -6", dict(bp=5, dg=-6.0)),
    ("run <= 5 and dG > -4", dict(bp=5, dg=-4.0)),
    ("run <= 6 and dG > -4", dict(bp=6, dg=-4.0)),
    ("run <= 8 and dG > -4", dict(bp=8, dg=-4.0)),
]

for label, gc in (("35% GC", 0.35), ("50% GC", 0.50), ("65% GC", 0.65)):
    primers = pool(gc)
    total = N * (N - 1) // 2
    print(f"\n### {label}, {N} primers, {total:,} pairs")
    print(f"{'policy':<28} {'pairs rejected':>15} {'%':>6} {'greedy panel':>13}")
    for name, kw in POLICIES:
        conflicts = conflict_set(primers, **kw)
        panel = greedy_compatible(primers, conflicts)
        print(
            f"{name:<28} {len(conflicts):>15,} {100*len(conflicts)/total:>5.1f}% {len(panel):>13}"
        )


def admitted_run_lengths():
    """What each policy lets through, in complementary run length.

    The safety-relevant measurement: a stability floor bounds stability, not
    length, so it can admit a long AT-rich run.
    """
    from neoswga.core.dimer import max_complementary_run

    print("\n## Longest complementary run each policy admits")
    print(f"{'pool':<9} {'policy':<24} {'admitted pairs':>15} {'worst run bp':>13} {'>= 8 bp':>9}")
    for label, gc in (("35% GC", 0.35), ("50% GC", 0.50), ("65% GC", 0.65)):
        primers = pool(gc)
        pairs = list(itertools.combinations(primers, 2))
        for name, kw in (
            ("run <= 3 (default)", dict(bp=3)),
            ("run <= 5", dict(bp=5)),
            ("dG > -6 only", dict(dg=-6.0)),
            ("dG > -4 only", dict(dg=-4.0)),
            ("dG > -2.79 only", dict(dg=-2.79)),
        ):
            conflicts = conflict_set(primers, **kw)
            admitted = [
                (a, b) for a, b in pairs if (a, b) not in conflicts and (b, a) not in conflicts
            ]
            runs = [max_complementary_run(a, b) for a, b in admitted]
            worst = max(runs) if runs else 0
            long = sum(1 for r in runs if r >= 8)
            print(f"{label:<9} {name:<24} {len(admitted):>15,} {worst:>13} {long:>9,}")


def per_pair_cost():
    """Whether the floor is affordable. Measured because the audit guessed."""
    import time

    from neoswga.core.dimer import is_dimer

    primers = pool(0.5)[:120]
    pairs = list(itertools.combinations(primers, 2))
    start = time.perf_counter()
    for a, b in pairs:
        is_dimer(a, b, max_dimer_bp=3)
    run = time.perf_counter() - start
    start = time.perf_counter()
    for a, b in pairs:
        is_dimer_thermodynamic(a, b, delta_g_threshold=-6.0, temperature=30.0)
    floor = time.perf_counter() - start
    print(f"\n## Cost over {len(pairs):,} pairs")
    print(f"complementary run : {run / len(pairs) * 1e6:6.1f} us/pair")
    print(f"free energy       : {floor / len(pairs) * 1e6:6.1f} us/pair")
    print(f"ratio             : {floor / run:6.1f}x")


admitted_run_lengths()
per_pair_cost()
