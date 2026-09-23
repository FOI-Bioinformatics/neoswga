#!/usr/bin/env python3
"""Refit the coverage reach under directional windows.

`docs/validation/reach_calibration.md` fitted 3.0-4.5 kb by asking which reach
reproduces the 60-75% breadth Clarke et al. (2017) measured at 10x depth for
their `TmL/Even` Wolbachia set, using SYMMETRIC windows. The geometry and the
value were fitted together, so correcting one without refitting the other is
unsound in either direction.

This recomputes that fit under directional windows, using the same sets, the
same genome and the same outcome. It needs no new data.

Step 1 rebuilds the SYMMETRIC table first and checks it against the published
one. A table that cannot be rebuilt cannot be refitted, so a mismatch there
stops the run.

    python scripts/benchmarking/refit_reach_directional.py
"""

import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from neoswga.core.coverage import site_spans  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/validation/data/clarke_2017_wolbachia.json"
GENOME = ROOT / "tests/validation/genomes/wolbachia.fna"

REACHES = (1000, 2000, 3000, 4000, 6000, 8000, 10000)

#: The table `reach_calibration.md` publishes, to check the rebuild against.
PUBLISHED = {
    "TmL/Even": (0.28, 0.48, 0.62, 0.72, 0.84, 0.92, 0.96),
    "TmL/Selective": (0.26, 0.45, 0.60, 0.71, 0.85, 0.92, 0.96),
    "Leichty and Brisson 2014": (0.20, 0.36, 0.48, 0.58, 0.71, 0.80, 0.85),
    "TmH/Even": (0.12, 0.24, 0.34, 0.43, 0.57, 0.68, 0.76),
    "TmH/Selective": (0.11, 0.20, 0.28, 0.36, 0.49, 0.59, 0.68),
}

#: Clarke et al. report 10x depth over 60-75% of the genome for `TmL/Even`.
BAND = (0.60, 0.75)
_COMP = str.maketrans("ACGT", "TGCA")


def rc(seq):
    return seq.translate(_COMP)[::-1]


def oriented_sites(genome, primers):
    """Forward and reverse occurrence lists, kept apart.

    `test_reach_calibration_wolbachia._sites` pools them into one set, which is
    all the symmetric convention needs and is exactly what the directional one
    cannot use. Circular, as that helper is: a match spanning the origin counts.
    """
    forward, reverse = [], []
    for primer in primers:
        padded = genome + genome[: len(primer) - 1]
        for pattern, bucket in ((primer, forward), (rc(primer), reverse)):
            bucket.extend(
                m.start() for m in re.finditer(f"(?={pattern})", padded) if m.start() < len(genome)
            )
    return sorted(set(forward)), sorted(set(reverse))


def coverage(forward, reverse, oligo_length, genome_len, reach, geometry):
    occupied = np.zeros(genome_len, dtype=bool)
    for low, high in site_spans(forward, reverse, oligo_length, reach, geometry):
        occupied[max(0, low) : min(genome_len, high)] = True
    return float(occupied.sum()) / genome_len


def band_edges(curve):
    """Where the curve enters and leaves the measured band, linearly.

    Returns `(enter, leave)` in bp, or None for an edge the sampled reaches do
    not bracket.
    """

    def crossing(target):
        for (r0, c0), (r1, c1) in zip(curve, curve[1:], strict=True):
            if c0 <= target <= c1 and c1 > c0:
                return r0 + (r1 - r0) * (target - c0) / (c1 - c0)
        return None

    return crossing(BAND[0]), crossing(BAND[1])


def main():
    if not GENOME.exists():
        sys.exit("missing tests/validation/genomes/wolbachia.fna; see the module docstring")
    genome = "".join(
        line.strip().upper() for line in GENOME.read_text().splitlines() if not line.startswith(">")
    )
    sets = json.loads(DATA.read_text())["sets"]
    n = len(genome)
    print(f"Wolbachia {n:,} bp, {len(sets)} published sets\n")

    tables = {}
    for geometry in ("symmetric", "directional"):
        print(f"== {geometry} ==")
        head = "".join(f"{r // 1000:>7}kb" for r in REACHES)
        print(f"{'set':<28}{head}")
        rows = {}
        for name, entry in sets.items():
            primers = entry["primers"]
            fwd, rev = oriented_sites(genome, primers)
            k = len(primers[0])
            curve = [(r, coverage(fwd, rev, k, n, r, geometry)) for r in REACHES]
            rows[name] = curve
            cells = "".join(f"{c:>9.2f}" for _, c in curve)
            print(f"{name:<28}{cells}")
        tables[geometry] = rows
        print()

    print("== Step 1: does the symmetric rebuild match the published table? ==")
    worst, offenders = 0.0, []
    for name, published in PUBLISHED.items():
        if name not in tables["symmetric"]:
            offenders.append(f"{name}: absent from the data file")
            continue
        got = [c for _, c in tables["symmetric"][name]]
        for reach, want, have in zip(REACHES, published, got, strict=True):
            delta = abs(want - have)
            worst = max(worst, delta)
            if delta > 0.02:
                offenders.append(f"{name} at {reach}: published {want:.2f}, rebuilt {have:.2f}")
    print(f"largest disagreement: {worst:.3f}")
    if offenders:
        print("\nSTOP. The published table does not rebuild:")
        for line in offenders:
            print(f"  {line}")
        print("\nA table that cannot be rebuilt cannot be refitted.")
        return
    print("matches within 0.02 everywhere.\n")

    print("== Step 3: the fit ==")
    print(f"`TmL/Even` reached 10x over {BAND[0]:.0%}-{BAND[1]:.0%} of the genome.\n")
    fits = {}
    for geometry in ("symmetric", "directional"):
        enter, leave = band_edges(tables[geometry]["TmL/Even"])
        fits[geometry] = (enter, leave)
        show = lambda v: "beyond the sampled reaches" if v is None else f"{v:,.0f} bp"  # noqa: E731
        print(f"{geometry:<13} enters at {show(enter):<26} leaves at {show(leave)}")
    s_enter, d_enter = fits["symmetric"][0], fits["directional"][0]
    if s_enter and d_enter:
        print(f"\ndirectional / symmetric at the lower edge: {d_enter / s_enter:.2f}x")
        print("The plan's prediction was about 2x. Anything else changes the decision.")


if __name__ == "__main__":
    main()
