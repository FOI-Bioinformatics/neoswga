"""Print the tables behind docs/validation/2026-11-08-background-aggregate-measured.md.

Reads WORKDIR/aggregate_control.json, written by
`diversity_baseline.py --aggregate-control`, and prints, for each reading and
requested size: the delivered size, the panel overlap between the readings,
coverage on every target strain, and the selectivity density (exact and
occupancy-weighted) of the reference strain against every host.
"""

import json
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT = os.path.join(
    REPO_ROOT, "tests", "validation", "genomes", "diversity_baseline", "aggregate_control.json"
)


def value(measurement):
    if measurement is None:
        return None
    return measurement.get("value")


def fmt(x, spec):
    return "n/a" if x is None else format(x, spec)


def main(path=DEFAULT):
    with open(path) as handle:
        record = json.load(handle)
    modes = record["modes"]
    print(f"design {record['design']}  git {record['git_sha'][:7]}")
    for mode, block in modes.items():
        print(f"{mode:>10}: step2 {block['step2_candidates']}, step3 {block['step3_candidates']}")
    reference = record["fg"][0]
    sizes = [row["requested_size"] for row in modes["sum"]["sizes"]]
    for size in sizes:
        rows = {m: next(r for r in b["sizes"] if r["requested_size"] == size) for m, b in modes.items()}
        panels = {m: set(r["delivered_primers"] or ()) for m, r in rows.items()}
        shared = len(panels["sum"] & panels["worst-case"])
        print(f"\nn={size}: delivered " + ", ".join(
            f"{m} {r['delivered_size']}" for m, r in rows.items()
        ) + f"; shared {shared}")
        for mode, row in rows.items():
            ev = row["evaluation"]
            cov = {k: value(v["coverage"]) for k, v in ev["per_target"].items()}
            print(f"  {mode:>10} coverage: " + "  ".join(f"{k} {fmt(c, '.4f')}" for k, c in cov.items()))
            print(f"  {mode:>10} worst target coverage {fmt(value(ev['worst_target_coverage']), '.4f')}")
            for pair in ev["target_host_pairs"]:
                if pair["target"] != reference:
                    continue
                host = ev["per_host"][pair["host"]]
                print(
                    f"  {mode:>10} vs {pair['host']:<14} sites {fmt(value(host['sites']), '>8.0f')}"
                    f"  sites/Mb {fmt(value(host['sites_per_mb']), '>9.2f')}"
                    f"  density {fmt(value(pair['selectivity_density']), '>9.2f')}"
                    f"  weighted {fmt(value(pair.get('weighted_selectivity_density')), '>8.2f')}"
                )
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
