"""What the candidate gate admits, and how informative it is.

Produces every measurement in
`docs/validation/pool_selection_audit_2026-09-18.md`:

  1. Occupancy spread within the admitted pool, which bounds how far Stage 1's
     unweighted bin count can misrank two candidates.
  2. How far each polymerase's Tm floor sits below its reaction temperature,
     and the occupancy of a primer sitting on that floor.
  3. Which bound of the Tm window binds as a function of primer length, which
     decides whether an additive enlarges the candidate pool or shrinks it.
  4. The composition of what an additive gains and loses.

Sampling is over random sequence at a fixed GC fraction rather than over a
genome's distinct k-mers. That answers what fraction of sequence space the gate
admits, which is a property of the gate; it is not a pool composition and should
not be read as one.

Usage:
    python scripts/benchmarking/candidate_gate_audit.py
    python scripts/benchmarking/candidate_gate_audit.py --gc 0.41 --n 40000
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
from typing import Dict, List, Sequence, Tuple

from neoswga.core.occupancy import (
    DEFAULT_MISMATCH_PENALTY_C,
    mismatch_tm,
    site_occupancy,
)
from neoswga.core.parameter import default_tm_range
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import calculate_enthalpy_entropy

SEED = 20260918

# Platform, reaction temperature, additive kwargs, label. The bst row is
# measured at k = 18 in the write-up, the middle of its documented 15-25 bp
# range, because judging it at k = 12 would be judging a configuration nobody
# would run.
PLATFORMS: Sequence[Tuple[str, float, Dict[str, float], str]] = (
    ("phi29", 30.0, {}, "phi29 30 C"),
    ("equiphi29", 42.0, {}, "equiphi29 42 C"),
    ("equiphi29", 42.0, {"dmso_percent": 5.0, "betaine_m": 1.0}, "equiphi29 +DMSO/betaine"),
    ("bst", 63.0, {}, "bst 63 C"),
)

# Representative enthalpies, spanning the range a 12-mer can take.
PROBE_PRIMERS = ("AAAAGGGCGTTA", "CACCGACGACGA", "GCGCGCGCGCGC", "ATATATATATAT")


def sample_primers(k: int, gc: float, n: int, seed: int = SEED) -> List[str]:
    """`n` random k-mers with per-base probabilities matching `gc`."""
    rng = random.Random(seed + k)
    p_gc, p_at = gc / 2.0, (1.0 - gc) / 2.0
    weights = [p_at, p_gc, p_gc, p_at]
    return ["".join(rng.choices("ACGT", weights=weights, k=k)) for _ in range(n)]


def gc_fraction(primer: str) -> float:
    return (primer.count("G") + primer.count("C")) / len(primer)


def occupancy_spread(k: int, gc: float, n: int) -> None:
    """Section: how much Stage 1's occupancy blindness can cost."""
    print(f"\n## Occupancy spread within the admitted pool (k={k}, GC={gc:.2f}, n={n})")
    print(
        f"{'platform':26s} {'admitted':>9s} {'theta p5':>9s} {'median':>8s} "
        f"{'p95':>8s} {'p95/p5':>8s} {'frac<0.5':>9s}"
    )
    primers = sample_primers(k, gc, n)
    for polymerase, temp, additives, label in PLATFORMS:
        low, high = default_tm_range(polymerase)
        conditions = ReactionConditions(temp=temp, polymerase=polymerase, **additives)
        thetas = []
        for primer in primers:
            tm = conditions.calculate_effective_tm(primer)
            if not low <= tm <= high:
                continue
            dh, _ds = calculate_enthalpy_entropy(primer)
            thetas.append(site_occupancy(dh, tm, temp))
        if not thetas:
            print(f"{label:26s} {0:9d}")
            continue
        thetas.sort()
        count = len(thetas)
        p5 = thetas[int(0.05 * count)]
        p95 = thetas[int(0.95 * count)]
        ratio = p95 / p5 if p5 > 0 else float("inf")
        below = sum(1 for t in thetas if t < 0.5) / count
        print(
            f"{label:26s} {count:9d} {p5:9.4f} {statistics.median(thetas):8.4f} "
            f"{p95:8.4f} {ratio:8.1f} {below:9.3f}"
        )


def window_versus_reaction_temperature() -> None:
    """Section: the Tm window is the wrong gate."""
    print("\n## Tm floor against reaction temperature")
    print(f"{'polymerase':12s} {'window':>14s} {'temp':>6s} {'floor - temp':>13s}")
    for polymerase, temp, _additives, _label in PLATFORMS:
        low, high = default_tm_range(polymerase)
        print(f"{polymerase:12s} {f'[{low:.0f},{high:.0f}]':>14s} {temp:6.0f} {low - temp:13.1f}")

    print("\nOccupancy of a primer sitting exactly on the floor:")
    for primer in PROBE_PRIMERS:
        dh, _ds = calculate_enthalpy_entropy(primer)
        line = f"  {primer}  dH={dh:8.1f} kcal/mol"
        for polymerase, temp in (("phi29", 30.0), ("equiphi29", 42.0)):
            low, _high = default_tm_range(polymerase)
            line += f"   theta@{polymerase}-floor={site_occupancy(dh, low, temp):.3f}"
        print(line)


def additive_sign_by_length(gc: float, n: int, k_range: Sequence[int]) -> None:
    """Section: which bound binds decides the sign of the additive's effect."""
    low, high = default_tm_range("equiphi29")
    plain = ReactionConditions(temp=42.0, polymerase="equiphi29")
    treated = ReactionConditions(temp=42.0, polymerase="equiphi29", dmso_percent=5.0, betaine_m=1.0)
    print(
        f"\n## Admitted fraction by primer length (equiphi29 42 C, window [{low:.0f},{high:.0f}])"
    )
    print(f"GC={gc:.2f}, n={n} per length, additive = DMSO 5% + betaine 1 M")
    print(f"{'k':>3s} {'plain':>10s} {'+additive':>10s} {'change':>9s} {'above ceiling':>14s}")
    for k in k_range:
        primers = sample_primers(k, gc, n)
        results = []
        for conditions in (plain, treated):
            kept = above = 0
            for primer in primers:
                tm = conditions.calculate_effective_tm(primer)
                if low <= tm <= high:
                    kept += 1
                elif tm > high:
                    above += 1
            results.append((kept / n, above / n))
        (kept_plain, above_plain), (kept_treated, _) = results
        marker = "  <- crossover" if kept_treated >= kept_plain else ""
        print(
            f"{k:3d} {kept_plain:10.3f} {kept_treated:10.3f} "
            f"{kept_treated - kept_plain:+9.3f} {above_plain:14.3f}{marker}"
        )


def additive_exchange(gc: float, n: int, lengths: Sequence[int]) -> None:
    """Section: the composition of what an additive gains and loses."""
    low, high = default_tm_range("equiphi29")
    plain = ReactionConditions(temp=42.0, polymerase="equiphi29")
    treated = ReactionConditions(temp=42.0, polymerase="equiphi29", dmso_percent=5.0, betaine_m=1.0)
    print("\n## What the additive gains and loses")
    for k in lengths:
        both: List[float] = []
        lost: List[float] = []
        gained: List[float] = []
        for primer in sample_primers(k, gc, n):
            in_plain = low <= plain.calculate_effective_tm(primer) <= high
            in_treated = low <= treated.calculate_effective_tm(primer) <= high
            if in_plain and in_treated:
                both.append(gc_fraction(primer))
            elif in_plain:
                lost.append(gc_fraction(primer))
            elif in_treated:
                gained.append(gc_fraction(primer))

        def describe(values: List[float]) -> str:
            if not values:
                return f"n={0:6d}"
            return f"n={len(values):6d} GC={statistics.mean(values):.3f}"

        print(
            f"  k={k:2d}  both: {describe(both)}   "
            f"lost: {describe(lost)}   gained: {describe(gained)}"
        )


def discrimination_by_gc_class(gc: float, n: int, k: int = 12) -> None:
    """Section: what an additive actually does, per GC class.

    Occupancy and mismatch discrimination move in opposite directions along the
    Tm axis, so the useful candidate is the one near its transition. `disc` is
    the per-candidate ratio theta(perfect) / theta(one mismatch). Their product
    is a reading aid, not a validated figure of merit.
    """
    plain = ReactionConditions(temp=42.0, polymerase="equiphi29")
    treated = ReactionConditions(temp=42.0, polymerase="equiphi29", dmso_percent=5.0, betaine_m=1.0)
    print(f"\n## Per GC class, plain against DMSO 5% + betaine 1 M (k={k}, GC={gc:.2f}, n={n})")
    print(f"one mismatch costed at {DEFAULT_MISMATCH_PENALTY_C:.0f} C")
    print(
        f"{'G+C':>4s} {'n':>6s} | {'Tm':>6s} {'theta':>7s} {'disc':>6s} | "
        f"{'Tm':>6s} {'theta':>7s} {'disc':>6s} | {'product':>16s}"
    )
    primers = sample_primers(k, gc, n)
    for count in range(0, k + 1):
        selected = [p for p in primers if p.count("G") + p.count("C") == count]
        if len(selected) < 30:
            continue
        summary = []
        for conditions in (plain, treated):
            tms, thetas, discs = [], [], []
            for primer in selected:
                dh, _ds = calculate_enthalpy_entropy(primer)
                tm = conditions.calculate_effective_tm(primer)
                theta = site_occupancy(dh, tm, conditions.temp)
                mismatched = site_occupancy(
                    dh, mismatch_tm(tm, 1, DEFAULT_MISMATCH_PENALTY_C), conditions.temp
                )
                tms.append(tm)
                thetas.append(theta)
                discs.append(theta / mismatched if mismatched > 0 else float("inf"))
            summary.append((statistics.mean(tms), statistics.mean(thetas), statistics.mean(discs)))
        (t1, h1, d1), (t2, h2, d2) = summary
        print(
            f"{count:4d} {len(selected):6d} | {t1:6.1f} {h1:7.4f} {d1:6.2f} | "
            f"{t2:6.1f} {h2:7.4f} {d2:6.2f} | {h1 * d1:6.2f} to {h2 * d2:6.2f}"
        )


def routes_to_specificity_conflict(gc: float, n: int, k: int = 12) -> None:
    """Section: the compositional and thermodynamic routes disagree on GC."""
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")
    gcs, discs = [], []
    for primer in sample_primers(k, gc, n):
        dh, _ds = calculate_enthalpy_entropy(primer)
        tm = conditions.calculate_effective_tm(primer)
        theta = site_occupancy(dh, tm, conditions.temp)
        mismatched = site_occupancy(
            dh, mismatch_tm(tm, 1, DEFAULT_MISMATCH_PENALTY_C), conditions.temp
        )
        if mismatched <= 0:
            continue
        gcs.append(gc_fraction(primer))
        discs.append(theta / mismatched)
    mean_gc, mean_disc = statistics.mean(gcs), statistics.mean(discs)
    covariance = sum((g - mean_gc) * (d - mean_disc) for g, d in zip(gcs, discs))
    spread = math.sqrt(
        sum((g - mean_gc) ** 2 for g in gcs) * sum((d - mean_disc) ** 2 for d in discs)
    )
    print(f"\n## The two routes to specificity (k={k}, GC={gc:.2f}, n={len(gcs)})")
    print(
        f"Pearson correlation, GC against thermodynamic discrimination: {covariance / spread:+.3f}"
    )
    print("The compositional route favours GC-rich against an AT-rich host, so they conflict.")


def tm_gate_against_occupancy_band(gc: float, n: int, band=(0.3, 0.8)) -> None:
    """Section: the Tm gate is the wrong parameterisation.

    Compares the shipped window with a band on occupancy at the reaction
    temperature, which is additive- and temperature-aware by construction.
    """
    tm_low, tm_high = default_tm_range("equiphi29")
    plain = ReactionConditions(temp=42.0, polymerase="equiphi29")
    treated = ReactionConditions(temp=42.0, polymerase="equiphi29", dmso_percent=5.0, betaine_m=1.0)
    print(f"\n## Shipped Tm window [{tm_low:.0f},{tm_high:.0f}] against an occupancy band {band}")
    print(f"equiphi29 42 C, GC={gc:.2f}, n={n} per length")
    print(
        f"{'k':>3s} {'reaction':>10s} {'Tm: n':>7s} {'Tm: disc':>9s} {'Tm: GC':>7s} "
        f"{'band: n':>8s} {'band: disc':>11s} {'band: GC':>9s}"
    )
    for k in (12, 14, 16, 18):
        primers = sample_primers(k, gc, n)
        precomputed = [(p, calculate_enthalpy_entropy(p)[0], gc_fraction(p)) for p in primers]
        for label, conditions in (("plain", plain), ("+additive", treated)):
            gates = []
            for gate in ("tm", "band"):
                kept, discs, gcs = 0, [], []
                for primer, dh, primer_gc in precomputed:
                    tm = conditions.calculate_effective_tm(primer)
                    theta = site_occupancy(dh, tm, conditions.temp)
                    admitted = (
                        tm_low <= tm <= tm_high if gate == "tm" else band[0] <= theta <= band[1]
                    )
                    if not admitted:
                        continue
                    mismatched = site_occupancy(
                        dh, mismatch_tm(tm, 1, DEFAULT_MISMATCH_PENALTY_C), conditions.temp
                    )
                    kept += 1
                    discs.append(theta / mismatched if mismatched > 0 else float("inf"))
                    gcs.append(primer_gc)
                gates.append(
                    (
                        kept,
                        statistics.mean(discs) if discs else float("nan"),
                        statistics.mean(gcs) if gcs else float("nan"),
                    )
                )
            (na, da, ga), (nb, db, gb) = gates
            print(
                f"{k:3d} {label:>10s} {na:7d} {da:9.2f} {ga:7.3f} " f"{nb:8d} {db:11.2f} {gb:9.3f}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gc", type=float, default=0.35, help="GC fraction of the sampled space")
    parser.add_argument("--n", type=int, default=20_000, help="primers sampled per length")
    parser.add_argument("--k", type=int, default=12, help="primer length for the spread table")
    args = parser.parse_args()

    window_versus_reaction_temperature()
    occupancy_spread(args.k, args.gc, args.n)
    occupancy_spread(18, args.gc, args.n)
    additive_sign_by_length(args.gc, args.n, range(12, 21))
    additive_exchange(args.gc, args.n, (12, 18))
    discrimination_by_gc_class(args.gc, args.n)
    routes_to_specificity_conflict(args.gc, args.n)
    tm_gate_against_occupancy_band(args.gc, args.n)


if __name__ == "__main__":
    main()
