"""What would a host-coverage field cost the search? (Known Issue 14)

Phase 2.1 of `docs/reports/PLAN_2026-11-02_known_issues.md`. Run from the
repository root against a prepared design directory:

    PYTHONPATH=. python docs/validation/host_coverage_cost_2026-11-02/measure.py \
        --params <work>/params.json --out results.json [--panels 40] [--size 12]

`stage2_background_weight` is REFUSED above zero because the weighted host axis
needs a host-coverage fraction and `PoolMetrics` carries five fields and no
sixth. The stated reason for five is cost: `plan_pool` spent 29.7 s, of which
27.6 s was inside `compute_metrics` across 4,040 calls. So the open question is
not whether a weighted host term is a good idea -- it is what carrying host
coverage through the search would cost, and that figure does not exist.

Three arms, because the alternative to adding the field is not "nothing".
`panel_refinement` already falls back to the FULL `compute_metrics` when a
`max_host_coverage` limit is set, so a host-coverage limit today already pays
the expensive path:

  A  compute_pool_metrics            -- as shipped, no host coverage
  B  A + one `_coverage_over_prefixes` call over the background
  C  compute_metrics                 -- what the search uses today when a
                                        host-coverage limit is set

B's marginal work is exactly one coverage call: `compute_pool_metrics` ALREADY
gathers and sorts every background position (for `total_bg_sites`) and then
throws the positions away. So the extra cost is the union-marking over a
boolean array of the host's length, not the position lookups, which are the
part that would have been expensive.

The marginal call is timed on positions gathered OUTSIDE the timed region,
because that gathering is already inside arm A and counting it twice would
overstate B.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import random
import statistics
import sys
import time

SEED = 20261102


def _load(params_path: str, pool_size: int):
    """A minimal optimizer over the prepared pool, as the search would build."""
    from types import SimpleNamespace

    from neoswga.core.base_optimizer import (
        BaseOptimizer,
        OptimizationResult,
        OptimizerConfig,
    )
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import build_reaction_conditions

    params = json.loads(pathlib.Path(params_path).read_text())
    data_dir = pathlib.Path(params["data_dir"])
    with open(data_dir / "step3_df.csv") as handle:
        rows = list(csv.DictReader(handle))
    key = "primer" if "primer" in rows[0] else list(rows[0])[1]
    primers = [r[key] for r in rows][:pool_size]

    fg = list(params["fg_prefixes"])
    bg = list(params.get("bg_prefixes") or [])
    cache = PositionCache(fg + bg, primers)
    conditions = build_reaction_conditions(SimpleNamespace(**params), from_mapping_only=True)

    class _Minimal(BaseOptimizer):
        def optimize(self, candidates, target_size=10, **kwargs):
            return OptimizationResult.failure("minimal", "not used")

        def _score_set(self, primers):
            return 0.0

        @property
        def name(self) -> str:
            return "minimal"

    optimizer = _Minimal(
        cache,
        fg,
        list(params["fg_seq_lengths"]),
        bg_prefixes=bg,
        bg_seq_lengths=list(params.get("bg_seq_lengths") or []),
        config=OptimizerConfig(target_set_size=params.get("num_primers", 12)),
        conditions=conditions,
    )
    return optimizer, primers, params


def _bg_positions(optimizer, panel):
    """The background positions arm A already gathers, built untimed."""
    by_prefix = {prefix: set() for prefix in optimizer.bg_prefixes}
    for primer in panel:
        for prefix in optimizer.bg_prefixes:
            by_prefix[prefix].update(
                optimizer.get_primer_positions(primer, prefix, "both").tolist()
            )
    return {p: sorted(v) for p, v in by_prefix.items()}


def _best_of(fn, repeats=3):
    """Minimum over repeats: the cleanest estimate of work done."""
    return min(_timed(fn) for _ in range(repeats))


def _timed(fn):
    started = time.perf_counter()
    fn()
    return time.perf_counter() - started


def _scaling(optimizer, sites_per_prefix=400):
    """How the marginal call scales with the host's LENGTH.

    `union_coverage` allocates `np.zeros(total_length, dtype=bool)`, so the
    marginal cost is set by how long the host is, not by how many sites the
    panel has there. Timed over synthetic position lists at a fixed site count
    so length is the only thing moving -- which is what decides whether this
    is affordable on a bacterial host and what it would mean on hg38.
    """
    rng = random.Random(SEED + 7)
    prefix = optimizer.bg_prefixes[0] if optimizer.bg_prefixes else "bg"
    out = {}
    for length in (1_000_000, 4_000_000, 16_000_000, 64_000_000):
        positions = {prefix: sorted(rng.sample(range(length), sites_per_prefix))}
        out[str(length)] = round(
            _best_of(
                lambda p=positions, n=length: optimizer._coverage_over_prefixes(p, [prefix], [n])
            )
            * 1000.0,
            3,
        )
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--panels", type=int, default=40)
    parser.add_argument("--size", type=int, default=12)
    parser.add_argument("--pool", type=int, default=300)
    args = parser.parse_args()

    from neoswga.core.pool_metrics import compute_pool_metrics

    optimizer, primers, params = _load(args.params, args.pool)
    rng = random.Random(SEED)
    panels = [rng.sample(primers, min(args.size, len(primers))) for _ in range(args.panels)]

    rows = []
    for panel in panels:
        bg_by_prefix = _bg_positions(optimizer, panel)

        arm_a = _best_of(lambda p=panel: compute_pool_metrics(optimizer, p))
        marginal = _best_of(
            lambda b=bg_by_prefix: optimizer._coverage_over_prefixes(
                b, optimizer.bg_prefixes, optimizer.bg_seq_lengths
            )
        )
        arm_c = _best_of(lambda p=panel: optimizer.compute_metrics(p))
        rows.append({"a": arm_a, "marginal": marginal, "c": arm_c})

    def ms(key):
        return statistics.median(r[key] for r in rows) * 1000.0

    a, m, c = ms("a"), ms("marginal"), ms("c")
    b = a + m
    out = {
        "panels": len(panels),
        "panel_size": len(panels[0]),
        "pool": len(primers),
        "fg_bp": sum(params["fg_seq_lengths"]),
        "bg_bp": sum(params.get("bg_seq_lengths") or []),
        "median_ms": {
            "A_compute_pool_metrics": round(a, 3),
            "B_with_host_coverage": round(b, 3),
            "C_compute_metrics": round(c, 3),
            "marginal_host_coverage": round(m, 3),
        },
        "B_over_A": round(b / a, 3) if a else None,
        "C_over_A": round(c / a, 3) if a else None,
        "C_over_B": round(c / b, 3) if b else None,
        "marginal_share_of_B": round(m / b, 4) if b else None,
        # Scaling of the marginal call with host length, at a fixed site count.
        "marginal_ms_by_host_bp": _scaling(optimizer),
        # The recorded profile the five-field decision was made against.
        "projected_over_4040_calls_s": {
            "A": round(a * 4040 / 1000.0, 2),
            "B": round(b * 4040 / 1000.0, 2),
            "C": round(c * 4040 / 1000.0, 2),
        },
    }
    pathlib.Path(args.out).write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
