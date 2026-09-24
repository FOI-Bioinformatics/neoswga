"""What a clamped depth tail costs `fit_reach`, on the Prevotella geometry.

`count_coverage` clamped its stop to the CONTIG length while `compute_bam_depth`
was handed the prefix's CONCATENATED length, so the tail stayed zero and entered
the fit as observed depth. This measures what that did to the reported
confidence, not to `best_reach`.

Synthetic depth: a true reach of 5,000 through the production kernel, Poisson
noise, four seeds. A direction, not a calibrated magnitude for a real profile.

    python clamped_tail_fit.py
"""
import sys

import numpy as np

from neoswga.core.reach_calibration import fit_reach, predicted_depth

# tests/validation/genomes/prevotella.fna, measured with reference_layout
REC1, REC2 = 1_796_408, 1_371_874
TOTAL = REC1 + REC2
TRUE_REACH = 5_000
N_SITES = 400
DEPTH_SCALE = 30


def main(seeds=(7, 11, 23, 41)) -> None:
    starts = [0, REC1]
    print(f"Prevotella geometry: {REC1:,} + {REC2:,} = {TOTAL:,} bp; "
          f"clamping fabricates {REC2 / TOTAL:.1%} of the array")
    for seed in seeds:
        rng = np.random.default_rng(seed)
        positions = sorted(rng.integers(0, TOTAL, N_SITES).tolist())
        expected = predicted_depth(positions, TOTAL, TRUE_REACH) * DEPTH_SCALE
        honest = rng.poisson(np.maximum(expected, 0.1)).astype(float)
        clamped = honest.copy()
        clamped[REC1:] = 0.0

        a = fit_reach(honest, positions, record_starts=starts)
        b = fit_reach(clamped, positions, record_starts=starts)
        print(
            f"seed {seed:>3}: honest rho={a.correlation:.3f} reach={a.best_reach:>6,} "
            f"plausible={len(a.plausible_reaches)} | "
            f"clamped rho={b.correlation:.3f} reach={b.best_reach:>6,} "
            f"plausible={len(b.plausible_reaches)} informative={b.informative}"
        )


if __name__ == "__main__":
    sys.exit(main())
