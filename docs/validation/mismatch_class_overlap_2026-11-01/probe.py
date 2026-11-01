"""How much background load the mismatch classes counted twice.

A site is one physical duplex and must be weighted at one occupancy. Because a
count file holds only one of each k-mer and its reverse complement, the same
stored k-mer can be a 1-mismatch neighbour of a primer read one way and a
2-mismatch neighbour read the other, so two classes claimed it and
`occupancy.weighted_site_load` weighted it at two different occupancies.

Run from the repository root:

    PYTHONPATH=. python docs/validation/mismatch_class_overlap_2026-11-01/probe.py \
        [--prefix <work>/host] [--k 10] [--primers N]

With no `--prefix` it builds its own 200 kb genome from a fixed seed, so the
whole measurement re-runs from this file alone.

Reports, at `max_mismatches=2`: how often the collision occurs, and what the
double count did to the summed class-1 and class-2 site load -- which is what
feeds `selectivity_ratio` and `selectivity_density`, the figures a panel is
ranked and accepted on.
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import tempfile
from collections import Counter

SEED = 20261101


def _build_table(k: int, size: int = 200_000) -> str:
    """A canonical count file over a synthetic genome, as jellyfish -C writes."""
    from neoswga.core.mismatch_counts import canonical_kmer

    rng = random.Random(SEED)
    seq = "".join(rng.choice("ACGT") for _ in range(size))
    counts = Counter(canonical_kmer(seq[i : i + k]) for i in range(size - k + 1))
    directory = tempfile.mkdtemp()
    prefix = os.path.join(directory, "g")
    with open(f"{prefix}_{k}mer_all.txt", "w") as handle:
        for kmer, count in counts.items():
            handle.write(f"{kmer} {count}\n")
    return prefix


def _classes_without_claiming(primer: str, prefix: str, k: int, depth: int) -> dict[int, int]:
    """The pre-fix arithmetic: each class counted in isolation.

    Reproduced here rather than imported, so this probe keeps measuring the
    defect after the code no longer contains it.
    """
    from neoswga.core.mismatch_counts import (
        _variants_at_distance,
        canonical_kmer,
        load_kmer_counts,
    )

    table = load_kmer_counts(prefix, k)
    primer_canonical = canonical_kmer(primer)
    out = {0: table.get(primer_canonical, 0)}
    for distance in range(1, depth + 1):
        forms = {canonical_kmer(v) for v in _variants_at_distance(primer, distance)}
        forms.discard(primer_canonical)
        out[distance] = sum(table.get(form, 0) for form in forms)
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--primers", type=int, default=300)
    args = parser.parse_args()

    from neoswga.core.mismatch_counts import (
        _variants_at_distance,
        canonical_kmer,
        mismatch_class_counts,
    )

    prefix = args.prefix or _build_table(args.k)
    rng = random.Random(SEED + 1)
    primers = ["".join(rng.choice("ACGT") for _ in range(args.k)) for _ in range(args.primers)]

    collided = 0
    fixed_total = 0
    prefix_total = 0
    worst = (0.0, None)
    for primer in primers:
        at_one = {canonical_kmer(v) for v in _variants_at_distance(primer, 1)}
        at_two = {canonical_kmer(v) for v in _variants_at_distance(primer, 2)}
        shares = bool(at_one & at_two)
        collided += shares

        after = mismatch_class_counts(primer, [prefix], max_mismatches=2)
        before = _classes_without_claiming(primer, prefix, args.k, 2)

        near_after = after[1] + after[2]
        near_before = before[1] + before[2]
        fixed_total += near_after
        prefix_total += near_before
        if shares and near_after:
            inflation = near_before / near_after
            if inflation > worst[0]:
                worst = (inflation, primer)

    print(f"k={args.k}, {len(primers)} random primers, max_mismatches=2")
    print(
        f"  primers with a form shared between class 1 and class 2: "
        f"{collided} ({collided / len(primers):.1%})"
    )
    print(f"  summed class-1 + class-2 sites, corrected : {fixed_total:,}")
    print(f"  summed class-1 + class-2 sites, pre-fix   : {prefix_total:,}")
    if fixed_total:
        over = prefix_total - fixed_total
        print(
            f"  double-counted sites: {over:,} " f"({over / fixed_total:.2%} of the corrected load)"
        )
    if worst[1]:
        print(f"  worst single primer: {worst[1]} at {worst[0]:.3f}x its corrected load")
    return 0


if __name__ == "__main__":
    sys.exit(main())
