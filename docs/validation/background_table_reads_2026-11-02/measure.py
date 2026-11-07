"""How much of a background table does a design actually read? (Known Issue 3)

Phase 1.1 of `docs/reports/PLAN_2026-11-02_known_issues.md`. Run from the
repository root against a prepared design directory:

    PYTHONPATH=. python docs/validation/background_table_reads_2026-11-02/measure.py \
        --params <work>/params.json --out results.json

Known Issue 3 reads "The filter command loads all background k-mers into
memory. Use Bloom filter for large backgrounds." The plan's Phase 1 proposed
moving the callers to `kmer_tables.counts_for`, and `mismatch_counts`' own
docstring already names that as the change that removes its ceiling. Before
writing it: measure which paths still materialise a table, and what fraction
of one a run actually looks up. If the fraction is small the dict is overhead
and the restructuring is justified by the ratio alone; if a run reads most of
the table, it is not.

Two numbers per k:

  table      distinct canonical k-mers the background table holds
  looked up  distinct canonical forms the occupancy path asks for, which is
             the exact match plus the 3k one-mismatch neighbours of every
             shortlisted candidate, deduplicated by canonical form

`max_mismatches` is read from params, since depth 2 raises the second number
by roughly an order of magnitude and that is the regime where the answer could
change.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import time


def _shortlist(data_dir: pathlib.Path) -> list[str]:
    """The candidates the occupancy ranking runs over."""
    for name in ("step3_df.csv", "step2_df.csv"):
        path = data_dir / name
        if not path.exists():
            continue
        with open(path) as handle:
            rows = list(csv.DictReader(handle))
        if not rows:
            continue
        key = "primer" if "primer" in rows[0] else list(rows[0])[1]
        return [r[key] for r in rows]
    raise SystemExit("no step2_df.csv or step3_df.csv in the data directory")


def _looked_up(primers, depth):
    """Distinct canonical forms the mismatch classes ask the table for."""
    from neoswga.core.mismatch_counts import _variants_at_distance, canonical_kmer

    forms = set()
    for primer in primers:
        forms.add(canonical_kmer(primer))
        for distance in range(1, depth + 1):
            forms.update(canonical_kmer(v) for v in _variants_at_distance(primer, distance))
    return forms


def _compare_and_weigh(prefix, k, forms):
    """Whole-table stream against one batched lookup: time, memory, agreement.

    The agreement check is the point. A batched lookup that returned 0 for a
    form the table holds would be a silent wrong answer, and
    `counts_for` guarantees every REQUESTED k-mer appears in its result -- so
    the hazard is asking outside the batch, not the batch itself.
    """
    import resource

    from neoswga.core import kmer_tables

    base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    started = time.perf_counter()
    whole = dict(kmer_tables.iter_table(prefix, k))
    stream_s = time.perf_counter() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    started = time.perf_counter()
    batched = kmer_tables.counts_for(prefix, k, sorted(forms))
    batched_s = time.perf_counter() - started

    agree = all(whole.get(f, 0) == batched.get(f, 0) for f in forms)
    held = len(whole)
    rise_mb = (peak - base) / 1024.0
    del whole
    return {
        "table_kmers": held,
        "stream_s": round(stream_s, 2),
        "batched_s": round(batched_s, 2),
        "speedup": round(stream_s / batched_s, 2) if batched_s else None,
        "whole_table_rss_mb": round(rise_mb, 1),
        "bytes_per_entry": round(rise_mb * 1024 * 1024 / held) if held else None,
        "counts_agree": agree,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--weigh",
        action="store_true",
        help="also time and weigh the whole-table read against one batched lookup",
    )
    args = parser.parse_args()

    from neoswga.core import kmer_tables

    params = json.loads(pathlib.Path(args.params).read_text())
    data_dir = pathlib.Path(params["data_dir"])
    depth = int(params.get("max_mismatches") or 1)
    bg_prefixes = list(params.get("bg_prefixes") or [])
    if not bg_prefixes:
        raise SystemExit("this measurement needs a background genome")

    primers = _shortlist(data_dir)
    by_k: dict[int, list[str]] = {}
    for primer in primers:
        by_k.setdefault(len(primer), []).append(primer)

    rows = []
    for k in sorted(by_k):
        prefix = bg_prefixes[0]
        if not kmer_tables.table_exists(prefix, k):
            continue
        started = time.perf_counter()
        table_size = sum(1 for _ in kmer_tables.iter_table(prefix, k))
        stream_s = time.perf_counter() - started

        forms = _looked_up(by_k[k], depth)
        row = {
            "k": k,
            "candidates": len(by_k[k]),
            "table_kmers": table_size,
            "looked_up": len(forms),
            "fraction_of_table": round(len(forms) / table_size, 6) if table_size else None,
            "stream_whole_table_s": round(stream_s, 2),
        }
        if args.weigh:
            row["batched_vs_stream"] = _compare_and_weigh(prefix, k, forms)
        rows.append(row)

    out = {
        "max_mismatches": depth,
        "bg_bp": sum(params.get("bg_seq_lengths") or []),
        "per_k": rows,
        "totals": {
            "table_kmers": sum(r["table_kmers"] for r in rows),
            "looked_up": sum(r["looked_up"] for r in rows),
        },
    }
    if out["totals"]["table_kmers"]:
        out["totals"]["fraction_of_table"] = round(
            out["totals"]["looked_up"] / out["totals"]["table_kmers"], 6
        )
    pathlib.Path(args.out).write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
