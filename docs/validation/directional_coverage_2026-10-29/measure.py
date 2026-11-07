"""Measure how much the symmetric coverage window inflates `fg_coverage`.

Run from the repo root, with a design's work directory:

    PYTHONPATH=. python docs/validation/directional_coverage_2026-10-29/measure.py \
        --prefix <work>/target --length 1270000 --panel <work>/step4_improved_df.csv

What is being compared
----------------------
SHIPPED (symmetric). `coverage._mark_window` marks
``occupied[pos-extension : pos+extension]`` -- a window ``2r`` wide centred on
every site -- and both coverage entry points default to ``strand="both"``, so
the position list holds the oligo's own occurrences and its reverse
complement's together. Every occurrence therefore earns coverage in the
direction it cannot extend in.

DIRECTIONAL. One direction per occurrence, which is what the chemistry allows:

- the oligo occurring literally at ``i`` anneals to the minus strand, the
  nascent strand is plus-sense, and extension runs toward INCREASING
  coordinates, covering ``[i, i + r]``;
- the reverse complement occurring at ``j`` means the oligo anneals to the plus
  strand and extension runs toward DECREASING coordinates, covering
  ``[j + len - r, j + len]``.

`position_cache` separates these: ``strand="forward"`` is the primer's own
entry and ``strand="reverse"`` is its reverse complement's
(`position_cache.py:369-386`), so the distinction is available without
re-scanning.

The two models are read off the SAME position lists, so the difference
reported is the geometry alone and not a different set of sites.

Scope. This measures coverage geometry. It does not re-select a panel under
the directional model, which would be a different and larger question: a
set-cover that knew each site covers ``r`` rather than ``2r`` would choose
differently.
"""

from __future__ import annotations

import argparse
import json

import numpy as np


def _positions(cache, prefix, primers, strand):
    out = {}
    for primer in primers:
        try:
            out[primer] = cache.get_positions(prefix, primer, strand)
        except Exception:
            out[primer] = np.array([], dtype=np.int64)
    return out


def directional_coverage(cache, prefix, primers, length, reach):
    """Covered fraction under one-direction-per-occurrence.

    Marked on a single bool array so overlapping windows are a union, exactly
    as the shipped path does; the only difference is which interval each
    occurrence contributes.
    """
    occupied = np.zeros(length, dtype=bool)
    forward = _positions(cache, prefix, primers, "forward")
    reverse = _positions(cache, prefix, primers, "reverse")
    marked_f = marked_r = 0
    for primer in primers:
        klen = len(primer)
        for pos in np.asarray(forward.get(primer, ()), dtype=np.int64):
            start, end = int(pos), min(length, int(pos) + reach)
            if end > start:
                occupied[start:end] = True
                marked_f += 1
        for pos in np.asarray(reverse.get(primer, ()), dtype=np.int64):
            start = max(0, int(pos) + klen - reach)
            end = min(length, int(pos) + klen)
            if end > start:
                occupied[start:end] = True
                marked_r += 1
    return float(occupied.sum()) / float(length), marked_f, marked_r


def symmetric_coverage(cache, prefix, primers, length, reach):
    """The shipped figure, through the shipped function."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    overall, per_prefix = compute_per_prefix_coverage(
        cache, primers, [prefix], [length], extension=reach, strand="both"
    )
    return float(overall)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--panel", required=True, help="step4_improved_df.csv")
    ap.add_argument("--set-index", type=int, default=0)
    ap.add_argument("--reaches", default="250,500,1000,2000,3000,6000")
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args()

    import pandas as pd

    from neoswga.core.position_cache import PositionCache

    frame = pd.read_csv(args.panel)
    panel = frame[frame["set_index"] == args.set_index]["primer"].astype(str).tolist()
    cache = PositionCache([args.prefix], panel)

    rows = []
    print(f"panel of {len(panel)} primers against {args.length:,} bp\n")
    print(f"{'reach':>7} {'symmetric':>10} {'directional':>12} {'inflation':>10} {'ratio':>7}")
    for reach in [int(r) for r in args.reaches.split(",")]:
        sym = symmetric_coverage(cache, args.prefix, panel, args.length, reach)
        dr, nf, nr = directional_coverage(cache, args.prefix, panel, args.length, reach)
        infl = sym - dr
        ratio = (sym / dr) if dr > 0 else float("nan")
        rows.append(
            {
                "reach": reach,
                "symmetric": sym,
                "directional": dr,
                "inflation_pp": infl * 100,
                "ratio": ratio,
                "forward_sites": nf,
                "reverse_sites": nr,
            }
        )
        print(f"{reach:>7} {sym:>10.4f} {dr:>12.4f} {infl * 100:>9.2f}pp {ratio:>7.3f}")

    if args.json_out:
        with open(args.json_out, "w") as handle:
            json.dump(
                {"panel_size": len(panel), "length": args.length, "rows": rows}, handle, indent=2
            )
        print(f"\nwrote {args.json_out}")


if __name__ == "__main__":
    main()
