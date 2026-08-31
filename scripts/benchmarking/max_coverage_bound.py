"""Exact and LP-relaxation bounds for the fixed-budget coverage problem.

NeoSWGA's optimizers all answer the same question: given a budget of S primers,
which S maximise foreground coverage?  The ILP that ships in
`dominating_set_optimizer.optimize_ilp` answers a different one -- the fewest
primers that cover every reachable bin -- so it is infeasible, and returns None,
for every budget below the minimum cover.  It cannot bound the greedy result.

This module builds the max-coverage formulation instead:

    maximise    sum_b w_b * y_b
    subject to  y_b <= sum_{i covers b} x_i     for every bin b
                sum_i x_i <= S
                x_i, y_b binary  (or in [0,1] for the LP relaxation)

`w_b` is the number of genome bases the bin stands for, so the objective is
covered bases rather than covered bins and the result is directly comparable to
the coverage a greedy run reports on the same bins.

Both bounds are reported.  The LP relaxation is always available and is a valid
upper bound even when the integer program is stopped early, which matters at
genome scale where the exact solve will not finish.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence


@dataclass
class BoundResult:
    """One bound on the best achievable coverage at a fixed budget."""

    budget: int
    kind: str  # "ilp" or "lp"
    coverage: float  # fraction of the genome, comparable to fg_coverage
    objective_bases: float
    selected: List[str] = field(default_factory=list)
    proven_optimal: bool = False
    gap: Optional[float] = None  # solver's own relative MIP gap
    seconds: float = 0.0
    status: str = ""


def build_bin_coverage(
    optimizer,
    candidates: Sequence[str],
) -> tuple[Dict[str, set], Dict[int, int], int]:
    """Map each primer to the bins it covers, reusing the shipped BipartiteGraph.

    Going through `BipartiteGraph` rather than re-deriving coverage keeps the
    bound on exactly the same bins the greedy optimizer sees; a bound computed
    against different bins would not be a bound on anything.

    Returns (primer -> bins, bin -> weight in bases, total genome bases).
    """
    from neoswga.core.dominating_set_optimizer import BipartiteGraph

    # The reach must match the optimizer's, or the bound is computed over
    # different bins and bounds nothing.
    reach = getattr(optimizer, "extension_reach", 0)
    graph = BipartiteGraph(bin_size=optimizer.bin_size)
    for primer in candidates:
        for prefix, length in zip(optimizer.fg_prefixes, optimizer.fg_seq_lengths):
            positions = optimizer.cache.get_positions(prefix, primer, "both")
            if len(positions) > 0:
                graph.add_primer_coverage(
                    primer, positions, prefix, length, extension_reach=reach
                )

    primer_to_bins = {p: set(graph.primer_to_regions.get(p, set())) for p in candidates}
    primer_to_bins = {p: bins for p, bins in primer_to_bins.items() if bins}

    # Weight each bin by the bases it actually spans. The last bin of a prefix
    # is short, and `optimize_greedy` divides covered bins by bin *count*, which
    # treats that short bin as full -- on the 6157 bp plasmid at bin_size=500
    # that alone reads as a 5.3% coverage difference with no primer changed.
    # Weighting by span puts both sides in the same units so the residual is a
    # real optimality gap rather than a rounding convention.
    bin_weight: Dict[object, int] = {}
    for bins in primer_to_bins.values():
        for b in bins:
            bin_weight.setdefault(b, b.end - b.start)

    total_bases = int(sum(optimizer.fg_seq_lengths))
    return primer_to_bins, bin_weight, total_bases


def covered_bases(optimizer, candidates, selected) -> tuple[int, int]:
    """Bases covered by `selected`, on the same bins the bound uses."""
    primer_to_bins, bin_weight, total_bases = build_bin_coverage(optimizer, candidates)
    bins = set()
    for p in selected:
        bins |= primer_to_bins.get(p, set())
    return sum(bin_weight[b] for b in bins), total_bases


def _solve(primer_to_bins, bin_weight, total_bases, budget, relax, max_seconds):
    # `sense` is the string 'MAX', not the `maximize` helper. Passing the
    # function silently yields a model the solver reports as empty, which is
    # how `dominating_set_optimizer.optimize_ilp` gets away with
    # `Model(sense=minimize)` -- it happens to fall through to the default.
    from mip import BINARY, CONTINUOUS, MAXIMIZE, Model, xsum

    # CoverageRegion is hashable but not orderable, so keep insertion order.
    bins = list(bin_weight)
    bin_to_primers: Dict[object, List[str]] = {b: [] for b in bins}
    for primer, covered in primer_to_bins.items():
        for b in covered:
            bin_to_primers[b].append(primer)

    vtype = CONTINUOUS if relax else BINARY
    model = Model(sense=MAXIMIZE)
    model.verbose = 0

    x = {p: model.add_var(var_type=vtype, lb=0.0, ub=1.0) for p in primer_to_bins}
    y = {b: model.add_var(var_type=vtype, lb=0.0, ub=1.0) for b in bins}

    model.objective = xsum(bin_weight[b] * y[b] for b in bins)
    for b in bins:
        model += y[b] <= xsum(x[p] for p in bin_to_primers[b])
    model += xsum(x.values()) <= budget

    start = time.time()
    status = model.optimize(max_seconds=max_seconds)
    elapsed = time.time() - start

    objective = float(model.objective_value or 0.0)
    selected = (
        sorted(p for p, v in x.items() if v.x is not None and v.x >= 0.99)
        if not relax
        else []
    )
    return BoundResult(
        budget=budget,
        kind="lp" if relax else "ilp",
        coverage=objective / total_bases if total_bases else 0.0,
        objective_bases=objective,
        selected=selected,
        proven_optimal=(status.name == "OPTIMAL"),
        gap=getattr(model, "gap", None),
        seconds=elapsed,
        status=status.name,
    )


def coverage_bounds(
    optimizer,
    candidates: Sequence[str],
    budget: int,
    max_seconds: int = 300,
) -> Dict[str, BoundResult]:
    """Return the LP and ILP bounds on coverage at `budget` primers."""
    primer_to_bins, bin_weight, total_bases = build_bin_coverage(optimizer, candidates)
    return {
        "lp": _solve(primer_to_bins, bin_weight, total_bases, budget, True, max_seconds),
        "ilp": _solve(primer_to_bins, bin_weight, total_bases, budget, False, max_seconds),
    }
