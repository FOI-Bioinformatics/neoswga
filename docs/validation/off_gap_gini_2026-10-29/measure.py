"""Does `off_gap_gini` say anything the figures already reported do not?

Phase 7.2 of `docs/reports/PLAN_2026-10-28_pool_design_trustworthiness.md`.
Run from the repository root:

    PYTHONPATH=. python docs/validation/off_gap_gini_2026-10-29/measure.py \
        --out results.json [--work DIR] [--sizes 6,12,18,24]

`off_gap_gini` is a Gini index of the distances between background binding
sites, computed per strand, averaged across strands and then across genomes.
It carries the second-largest recorded weight in swga 2.0's set-level model,
and nothing in this repository computed it. It is REPORTED here and never
scored; `neoswga/core/background_evenness.py` says why.

The question a report has to answer is whether the figure is informative. A
quantity that moves exactly with `bg_coverage` or with `total_bg_sites` would
add a column and no information. So this records, over several panel sizes on
one pool:

  - `off_gap_gini` against the three background figures already reported,
  - the direction check: a deliberately clustered host against a dispersed one
    of the same site COUNT, through the same function, which is what the
    positive fitted weight claims to see and the counts cannot,
  - the cost, since the figure is computed inside `_compute_metrics`, which
    runs per candidate panel during the search.

The genomes are synthetic and generated here from a fixed seed, so the whole
measurement is re-runnable from this file alone: no reference genome is
reachable from this container (`tests/validation/genomes/` is gitignored and
absent, NCBI is refused by the egress proxy). Site geometry is therefore
Poisson where a real genome clusters, which is a limitation stated in the
record and the reason the direction check is done separately and explicitly.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import random
import statistics
import subprocess
import sys
import time

TARGET_LENGTH = 1_270_000
HOST_LENGTH = 4_000_000
TARGET_GC = 0.35
HOST_GC = 0.42
SEED = 20261029


def _sequence(length: int, gc: float, rng: random.Random) -> str:
    """A sequence of the requested length at the requested GC fraction."""
    weights = [(1.0 - gc) / 2, gc / 2, gc / 2, (1.0 - gc) / 2]
    return "".join(rng.choices("ATGC", weights=weights, k=length))


def _write_fasta(path: pathlib.Path, name: str, sequence: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 70):
            handle.write(sequence[start : start + 70] + "\n")


def build_pool(work: pathlib.Path) -> pathlib.Path:
    """The two genomes and a params.json, written once and reused."""
    work.mkdir(parents=True, exist_ok=True)
    data_dir = work / "work"
    data_dir.mkdir(exist_ok=True)
    target = work / "target.fasta"
    host = work / "host.fasta"
    if not target.exists():
        rng = random.Random(SEED)
        _write_fasta(target, "target", _sequence(TARGET_LENGTH, TARGET_GC, rng))
        _write_fasta(host, "host", _sequence(HOST_LENGTH, HOST_GC, rng))
    params = {
        "schema_version": 1,
        "fg_genomes": [str(target)],
        "fg_prefixes": [str(data_dir / "target")],
        "bg_genomes": [str(host)],
        "bg_prefixes": [str(data_dir / "host")],
        "data_dir": str(data_dir) + "/",
        "fg_seq_lengths": [TARGET_LENGTH],
        "bg_seq_lengths": [HOST_LENGTH],
        "min_k": 10,
        "max_k": 12,
        "polymerase": "equiphi29",
        "reaction_temp": 43.0,
        "min_fg_freq": 1e-5,
        "max_bg_freq": 5e-6,
        "max_gini": 0.7,
        "max_primer": 500,
        "num_primers": 12,
    }
    path = work / "params.json"
    path.write_text(json.dumps(params, indent=2), encoding="utf-8")
    return path


def _neoswga(step: str, params: pathlib.Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", step, "-j", str(params)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        sys.stderr.write(result.stdout[-4000:] + result.stderr[-4000:])
        raise SystemExit(f"{step} failed")


def prepare(params: pathlib.Path) -> None:
    """Steps 1 to 3, which do not depend on the panel size."""
    data_dir = pathlib.Path(json.loads(params.read_text())["data_dir"])
    if (data_dir / "step3_df.csv").exists():
        return
    for step in ("count-kmers", "filter", "prepare-candidates"):
        _neoswga(step, params)


def optimize_at(params: pathlib.Path, size: int) -> dict:
    """One `optimize` run at this panel size; returns the summary metrics."""
    config = json.loads(params.read_text())
    config["num_primers"] = size
    config["target_set_size"] = size
    params.write_text(json.dumps(config, indent=2), encoding="utf-8")
    started = time.time()
    _neoswga("optimize", params)
    elapsed = time.time() - started
    summary = json.loads(
        (pathlib.Path(config["data_dir"]) / "step4_improved_df_summary.json").read_text()
    )
    metrics = summary["metrics"]
    host_prefix = config["bg_prefixes"][0]
    target_prefix = config["fg_prefixes"][0]
    stats = metrics.get("strand_stats") or {}
    return {
        "num_primers": size,
        "delivered": len(summary.get("primers") or []),
        "seconds": round(elapsed, 1),
        "off_gap_gini": metrics.get("off_gap_gini"),
        "host_gap_gini_per_strand": (stats.get(host_prefix) or {}).get("gap_gini_per_strand"),
        "target_gap_gini_per_strand": (stats.get(target_prefix) or {}).get("gap_gini_per_strand"),
        "bg_coverage": metrics.get("bg_coverage"),
        "total_bg_sites": metrics.get("total_bg_sites"),
        "selectivity_density": metrics.get("selectivity_density"),
        "gap_gini_pooled_target": metrics.get("gap_gini"),
        "fg_coverage": metrics.get("fg_coverage"),
    }


def direction_check() -> dict:
    """Same site COUNT, two geometries, through the shipped function.

    This is the claim the positive fitted weight makes and the one the site
    counts cannot see: a host whose sites are clustered leaves most of itself
    out of reach, and the counts are identical either way.
    """
    from neoswga.core.background_evenness import per_strand_gap_gini

    rng = random.Random(SEED)
    count, span = 400, HOST_LENGTH
    dispersed = {
        "forward": sorted(rng.sample(range(span), count)),
        "reverse": sorted(rng.sample(range(span), count)),
    }
    # The same number of sites, packed into 2% of the host in ten clusters.
    clustered: dict[str, list[int]] = {"forward": [], "reverse": []}
    for strand in clustered:
        for cluster in range(10):
            origin = cluster * (span // 10)
            clustered[strand].extend(
                sorted(origin + rng.randrange(span // 500) for _ in range(count // 10))
            )
        clustered[strand] = sorted(set(clustered[strand]))
    return {
        "sites_per_strand": count,
        "dispersed": per_strand_gap_gini(dispersed),
        "clustered": per_strand_gap_gini(clustered),
    }


def cost_check(work: pathlib.Path) -> dict:
    """What the per-strand Gini adds to one `compute_strand_alternation_stats`.

    Timed against the same call with the arithmetic skipped, on the prepared
    index, because the figure is computed inside `_compute_metrics` and that
    runs per candidate panel during the search.
    """
    from neoswga.core import background_evenness
    from neoswga.core.position_cache import PositionCache

    data_dir = pathlib.Path(json.loads((work / "params.json").read_text())["data_dir"])
    rows = (data_dir / "step3_df.csv").read_text().splitlines()[1:21]
    primers = [row.split(",")[0] for row in rows if row]
    prefix = str(data_dir / "target")
    cache = PositionCache([prefix], primers)

    def timed(repeats=5):
        best = []
        for _ in range(repeats):
            started = time.time()
            cache.compute_strand_alternation_stats(prefix, primers, TARGET_LENGTH)
            best.append(time.time() - started)
        return min(best)

    with_gini = timed()
    original = background_evenness.per_strand_gap_gini
    background_evenness.per_strand_gap_gini = lambda _by_strand: None
    try:
        import neoswga.core.position_cache as pc

        patched = pc.per_strand_gap_gini
        pc.per_strand_gap_gini = lambda _by_strand: None
        without = timed()
        pc.per_strand_gap_gini = patched
    finally:
        background_evenness.per_strand_gap_gini = original
    return {
        "primers": len(primers),
        "with_gini_ms": round(with_gini * 1000, 2),
        "without_gini_ms": round(without * 1000, 2),
    }


def _spearman(left: list[float], right: list[float]) -> float | None:
    """Rank correlation, or `None` when there is nothing to correlate."""
    pairs = [(a, b) for a, b in zip(left, right, strict=True) if a is not None and b is not None]
    if len(pairs) < 3:
        return None

    def ranks(values):
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0.0] * len(values)
        for rank, index in enumerate(order):
            out[index] = float(rank)
        return out

    first, second = ranks([p[0] for p in pairs]), ranks([p[1] for p in pairs])
    if len(set(first)) < 2 or len(set(second)) < 2:
        return None
    return round(statistics.correlation(first, second), 3)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--work", default=None)
    parser.add_argument("--sizes", default="6,12,18,24")
    args = parser.parse_args()

    work = pathlib.Path(args.work or "off_gap_gini_work").resolve()
    params = build_pool(work)
    prepare(params)

    rows = [optimize_at(params, int(size)) for size in args.sizes.split(",")]
    gini_values = [row["off_gap_gini"] for row in rows]
    results = {
        "panels": rows,
        "direction_check": direction_check(),
        "cost": cost_check(work),
        "rank_correlation_with_off_gap_gini": {
            name: _spearman(gini_values, [row[name] for row in rows])
            for name in ("bg_coverage", "total_bg_sites", "selectivity_density", "fg_coverage")
        },
    }
    pathlib.Path(args.out).write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
