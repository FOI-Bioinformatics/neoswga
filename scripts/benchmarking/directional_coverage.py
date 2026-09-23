#!/usr/bin/env python3
"""Symmetric against directional coverage windows, measured.

`coverage._mark_window` marks `[pos - r, pos + r)` around every site, and the
sites come from `get_positions(..., strand="both")`. A site can only extend one
way: an oligo occurring literally at `i` anneals to the minus strand and extends
toward increasing coordinates, and an occurrence of its reverse complement at
`j` extends toward decreasing coordinates.

Two directional variants are scored, and the distinction matters:

  onesided_r    the window is `r` wide, on the reachable side
  onesided_2r   the window is `2r` wide -- the same total per site as the
                symmetric model credits -- on the reachable side

`onesided_2r` isolates GEOMETRY from MAGNITUDE. It is also the variant that
reproduces the figure recorded in `reach_calibration.md`.

Reads the shipped HDF5 indexes directly, so it needs no `optimize` run and
works on indexes `optimize` refuses.

    python scripts/benchmarking/directional_coverage.py panels
    python scripts/benchmarking/directional_coverage.py trend
"""

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from neoswga.core.thermodynamics import reverse_complement  # noqa: E402

GENOMES = Path("tests/validation/genomes")
PREVOTELLA_BP = 3_168_282


def mark(spans, length):
    occupied = np.zeros(length, dtype=bool)
    for low, high in spans:
        occupied[max(0, low) : min(length, high)] = True
    return occupied


def symmetric(forward, reverse, _oligo_length, reach):
    sites = np.unique(np.concatenate([forward, reverse])) if len(forward) + len(reverse) else []
    return [(int(p) - reach, int(p) + reach) for p in sites]


def directional(forward, reverse, oligo_length, reach, width=None):
    width = reach if width is None else width
    spans = [(int(i), int(i) + width) for i in forward]
    spans += [(int(j) + oligo_length - width, int(j) + oligo_length) for j in reverse]
    return spans


def sites_for(handle, primers):
    forward, reverse, oligo_length = [], [], None
    for primer in primers:
        oligo_length = len(primer)
        if primer in handle:
            forward.append(np.array(handle[primer][:], dtype=np.int64))
        rc = reverse_complement(primer)
        if rc in handle:
            reverse.append(np.array(handle[rc][:], dtype=np.int64))
    empty = np.array([], dtype=np.int64)
    return (
        np.concatenate(forward) if forward else empty,
        np.concatenate(reverse) if reverse else empty,
        oligo_length,
    )


def gaps_of(occupied):
    edges = np.flatnonzero(np.diff(np.concatenate([[0], occupied.view(np.int8), [0]])))
    runs, out, previous = list(zip(edges[::2], edges[1::2], strict=True)), [], 0
    for start, end in runs:
        if start > previous:
            out.append(start - previous)
        previous = end
    if previous < len(occupied):
        out.append(len(occupied) - previous)
    return np.array(out) if out else np.array([0])


def panels(reach=3000):
    """Delivered panels, scored three ways, with map agreement and gaps."""
    print(f"Prevotella {PREVOTELLA_BP:,} bp, reach {reach:,}\n")
    header = f"{'panel':<10}{'k':>3}{'sites':>7}{'symmetric':>12}{'1-sided r':>12}{'1-sided 2r':>13}{'agree':>8}"
    print(header)
    for name in ("out10", "out11", "out12", "out12big"):
        csv = GENOMES / name / "step4_improved_df.csv"
        if not csv.exists():
            continue
        frame = pd.read_csv(csv)
        if "set_index" in frame.columns:
            frame = frame[frame["set_index"] == 0]
        primers = frame["primer"].tolist()
        k = len(primers[0])
        with h5py.File(GENOMES / f"prevotella_{k}mer_positions.h5", "r") as db:
            fwd, rev, length = sites_for(db, primers)
        sym = mark(symmetric(fwd, rev, length, reach), PREVOTELLA_BP)
        one = mark(directional(fwd, rev, length, reach), PREVOTELLA_BP)
        two = mark(directional(fwd, rev, length, reach, 2 * reach), PREVOTELLA_BP)
        agree = 100 * (sym & two).sum() / max(sym.sum(), 1)
        print(
            f"{name:<10}{k:>3}{len(fwd) + len(rev):>7}{sym.sum() / PREVOTELLA_BP:>12.4f}"
            f"{100 * (one.sum() - sym.sum()) / sym.sum():>11.1f}%"
            f"{100 * (two.sum() - sym.sum()) / sym.sum():>12.1f}%{agree:>7.1f}%"
        )
    print("\n`agree` is the share of the symmetric claim that one-sided 2r also claims.")
    print("The totals nearly cancel; the maps do not.")


def trend(sizes=(8, 16, 32, 64), pool=4000, seed=20260923):
    """Does the gap grow with set size, as `reach_calibration.md` records?"""
    with h5py.File(GENOMES / "prevotella_12mer_positions.h5", "r") as db:
        seen, candidates = set(), []
        for key in db:
            if len(key) != 12 or key in seen:
                continue
            rc = reverse_complement(key)
            seen.update({key, rc})
            fwd = np.array(db[key][:], dtype=np.int64)
            rev = np.array(db[rc][:], dtype=np.int64) if rc in db else np.array([], dtype=np.int64)
            if len(fwd) + len(rev) >= 6:
                candidates.append((key, fwd, rev))
    rng = np.random.default_rng(seed)
    if len(candidates) > pool:
        candidates = [candidates[i] for i in rng.choice(len(candidates), pool, replace=False)]

    for reach in (3000, 4500, 6200):
        covered, chosen, rows = np.zeros(PREVOTELLA_BP, dtype=bool), [], {}
        taken = set()
        for step in range(1, max(sizes) + 1):
            best, best_gain, best_mark = None, -1, None
            for key, fwd, rev in candidates:
                if key in taken:
                    continue
                m = mark(symmetric(fwd, rev, 12, reach), PREVOTELLA_BP)
                gain = int((m & ~covered).sum())
                if gain > best_gain:
                    best, best_gain, best_mark = (key, fwd, rev), gain, m
            if best is None or best_gain <= 0:
                break
            chosen.append(best)
            taken.add(best[0])
            covered |= best_mark
            if step in sizes:
                F = np.concatenate([c[1] for c in chosen])
                V = np.concatenate([c[2] for c in chosen])
                s = mark(symmetric(F, V, 12, reach), PREVOTELLA_BP).sum()
                d1 = mark(directional(F, V, 12, reach), PREVOTELLA_BP).sum()
                d2 = mark(directional(F, V, 12, reach, 2 * reach), PREVOTELLA_BP).sum()
                rows[step] = (s / PREVOTELLA_BP, 100 * (s - d1) / s, 100 * (s - d2) / s)
        print(f"\nreach {reach:,} bp   (positive = symmetric claims more)")
        print(f"{'n':>4}{'symmetric':>12}{'vs 1-sided r':>15}{'vs 1-sided 2r':>16}")
        for n in sorted(rows):
            c, o1, o2 = rows[n]
            print(f"{n:>4}{c:>12.4f}{o1:>14.1f}%{o2:>15.1f}%")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "panels"
    {"panels": panels, "trend": trend}[which]()
