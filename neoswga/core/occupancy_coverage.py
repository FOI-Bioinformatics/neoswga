"""Coverage weighted by how much of the time each binding site is occupied.

Extracted from `base_optimizer` on 2026-09-17, when rewriting the computation
pushed that module past its size budget. The answer to that is extraction
rather than a larger budget.

The rewrite is why this module exists. The previous implementation allocated a
boolean window and a float32 accumulator the length of the target and made two
full passes over both per primer, so its cost was linear in the genome and in
the panel and independent of how many sites there actually were. On the real
Wolbachia design it was 95% of one objective evaluation, while the 144 Mb host
genome everyone assumed was responsible was 16%. See
`docs/validation/objective_evaluation_cost_2026-09-17.md`.

`base_optimizer.BaseOptimizer._compute_effective_coverage` delegates here and
supplies the reach and the geometry from its own config, so every existing
caller is unchanged.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Dict, Mapping, Optional, Sequence

from .coverage import merged_window_intervals
from .occupancy import site_occupancy
from .thermodynamics import calculate_enthalpy_entropy


def occupancy_weighted_coverage(
    positions_by_primer: Mapping[str, Sequence[int]],
    total_length: int,
    *,
    extension_reach: int,
    circular: bool,
    conditions,
    reverse_by_primer: Optional[Mapping[str, Sequence[int]]] = None,
    geometry: str = "symmetric",
) -> Optional[float]:
    """Coverage weighted by how much of the time each site is occupied.

    `base_optimizer._compute_coverage` unions binding windows as booleans: a site either
    covers its window or does not. That is the right question for "could
    this primer reach here", and the wrong one for "will it". A primer whose
    effective Tm sits well below the reaction temperature is mostly not
    bound, and its sites contribute far less amplification than a count
    implies.

    A site covers its window with probability theta -- the same two-state
    occupancy the selectivity metric uses -- and the product is taken over
    PRIMERS, not over sites:

        P(covered at x) = 1 - PRODUCT over primers p reaching x of (1 - theta_p)

    The grouping is deliberate and is what the loop below implements: one
    primer's overlapping windows are unioned first and its occupancy applied
    once, because a primer does not stack with itself. Per-site
    independence would multiply (1 - theta) in once per overlapping window
    and report a larger number.

    Corrected 2026-09-14 (audit F5): this formula previously read "PRODUCT
    over sites", describing a model the code does not implement. The two
    differ measurably -- at T = Tm with 100 bp windows on a 1 kb target,
    sites at 500 and 510 give 10.5% under one primer and 15.25% under two
    distinct primers -- and `tests/test_occupancy_grouping_is_specified.py`
    pins the one in use.

    Independence across primers is itself an approximation: sites on one
    template molecule compete for polymerase. Neither model here has been
    compared against a measured reaction, so this is a stated approximation
    rather than a validated one, and the earlier claim that it bounds
    single-molecule recovery from above is not established by this
    arithmetic.

    Returns None when no reaction conditions are attached, since without them
    there is no temperature at which to evaluate occupancy and a fabricated
    number here would be indistinguishable from a measured one.
    """
    if conditions is None:
        return None
    if not positions_by_primer or total_length <= 0:
        return 0.0

    temp = conditions.temp

    # Accumulated in log space over window EDGES, not over bases.
    #
    # This loop used to allocate a boolean window and a float32 accumulator
    # the length of the target and make two full passes over both per
    # primer: `window[:] = False` and `not_covered[window] *= 1 - theta`.
    # Its cost was therefore linear in the genome and in the panel and
    # independent of how many sites there were, which is the wrong way
    # round -- a 1.27 Mb target with 20 sites per primer paid for 1.27
    # million bases either way. Measured, it was 95% of one objective
    # evaluation, and the 144 Mb host everyone blamed was 16%.
    #
    # Taking logs turns the product over primers into a sum, and a sum has
    # a difference array: each primer's merged window contributes its
    # log(1 - theta) at one edge and removes it at another. The answer is
    # then a walk over the edges, which number twice the merged windows
    # rather than once the genome.
    log_weight_at: Dict[int, float] = defaultdict(float)
    saturated_at: Dict[int, int] = defaultdict(int)

    # `reverse_by_primer` is absent under the default geometry, where the
    # caller's `positions_by_primer` is already the pooled union and the
    # orientation would not be used. Under `directional` the two are needed
    # apart, and a primer with sites in only one of them is ordinary.
    reverse_sites = reverse_by_primer or {}

    for primer, positions in positions_by_primer.items():
        backward = list(reverse_sites.get(primer, ()))
        if not positions and not backward:
            continue
        tm = conditions.calculate_effective_tm(primer)
        dh, _ = calculate_enthalpy_entropy(primer)
        theta = site_occupancy(dh, tm, temp)
        if theta <= 0.0:
            continue

        # One primer's windows are merged first, then its occupancy applied
        # once. Applying it per site would multiply (1 - theta) in for every
        # overlapping window of the SAME primer, which understates coverage
        # where a primer binds densely -- the windows overlap, the primer
        # does not stack with itself.
        spans = merged_window_intervals(
            positions,
            extension_reach,
            total_length,
            circular,
            reverse=backward,
            geometry=geometry,
        )
        if not spans:
            continue

        if theta >= 1.0:
            # log(1 - 1) is negative infinity, and adding it at one edge
            # then subtracting it at another gives NaN. A primer that is
            # always bound is counted rather than summed.
            #
            # No valid reaction reaches this: it needs a Tm near 200 C and a
            # 12-mer tops out around 70, so within the polymerase's own
            # temperature band the most extreme real primer gives
            # 0.999999999999. It is handled because it is representable and
            # because the alternative failure is a silent NaN.
            for start, end in spans:
                saturated_at[start] += 1
                saturated_at[end] -= 1
        else:
            weight = math.log1p(-theta)
            for start, end in spans:
                log_weight_at[start] += weight
                log_weight_at[end] -= weight

    if not log_weight_at and not saturated_at:
        return 0.0

    covered = 0.0
    running_log = 0.0
    running_saturated = 0
    edges = sorted(set(log_weight_at) | set(saturated_at))
    previous = edges[0]
    for edge in edges:
        # `running_log < 0.0` rather than `!= 0.0`: adding and removing the
        # same weights in a different order leaves a residue of a few float
        # epsilons, and a positive residue would contribute a negative width.
        if edge > previous and (running_saturated or running_log < 0.0):
            width = edge - previous
            covered += width if running_saturated else width * (1.0 - math.exp(running_log))
        running_log += log_weight_at.get(edge, 0.0)
        running_saturated += saturated_at.get(edge, 0)
        previous = edge

    return covered / total_length
