"""Convert a lab-result CSV into the JSON format `neoswga active-learn`
accepts (`--experimental-results`).

CSV schema (columns; those marked * are required):

    primer_set_id*    str - identifier (ignored by active-learn, kept for logs)
    primers*          str - primer sequences, semicolon-separated
    enrichment_fold*  float - target/background ratio measured by sequencing
    uniformity*       float - coverage uniformity 0-1 (higher = flatter)
    coverage          float - overall target coverage 0-1
    total_amplification   float - DNA yield (ng/uL or fold)
    off_target_fraction   float - fraction of reads off-target
    temperature       float - amplification temperature (C); default 30
    time_hours        float - amplification time; default 4
    primer_concentration  float - primer uM; default 1
    replicate_cv      float - CV across replicates
    passed_qc         bool  - 1/0 or true/false; default true
    timestamp         str   - ISO-8601; defaults to file mtime
    notes             str

Usage:

    python scripts/csv_to_experimental_results_json.py \\
        --input results.csv \\
        --output experimental_results.json

Designed to read a spreadsheet export (`File > Save as CSV`) and not
require any extra dependencies beyond the standard library.
"""

from __future__ import annotations

import argparse
import csv
import datetime as _dt
import json
import pathlib
import sys
from typing import Any, Dict, List


REQUIRED = {"primers", "enrichment_fold", "uniformity"}


def _to_float(value: str | None, default: float | None = None) -> float | None:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _to_bool(value: str | None, default: bool = True) -> bool:
    if value is None or value == "":
        return default
    s = str(value).strip().lower()
    if s in {"1", "true", "t", "yes", "y", "pass"}:
        return True
    if s in {"0", "false", "f", "no", "n", "fail"}:
        return False
    return default


def row_to_experimental_result(row: Dict[str, str], default_timestamp: str) -> Dict[str, Any]:
    primers = [p.strip().upper() for p in (row.get("primers") or "").split(";") if p.strip()]
    if not primers:
        raise ValueError(f"row has no primers: {row}")

    enrichment_fold = _to_float(row.get("enrichment_fold"))
    uniformity = _to_float(row.get("uniformity"))
    if enrichment_fold is None or uniformity is None:
        raise ValueError(
            f"enrichment_fold and uniformity are required (got {row})"
        )

    return {
        "primer_set": primers,
        "timestamp": (row.get("timestamp") or default_timestamp).strip(),
        "enrichment_fold": enrichment_fold,
        "uniformity_score": uniformity,
        "total_amplification": _to_float(row.get("total_amplification")),
        "off_target_fraction": _to_float(row.get("off_target_fraction")),
        "temperature": _to_float(row.get("temperature"), 30.0) or 30.0,
        "time_hours": _to_float(row.get("time_hours"), 4.0) or 4.0,
        "primer_concentration": _to_float(row.get("primer_concentration"), 1.0) or 1.0,
        "replicate_cv": _to_float(row.get("replicate_cv")),
        "passed_qc": _to_bool(row.get("passed_qc"), default=True),
        "notes": (row.get("notes") or "").strip(),
    }


def convert(input_path: pathlib.Path) -> List[Dict[str, Any]]:
    default_timestamp = _dt.datetime.fromtimestamp(input_path.stat().st_mtime).isoformat()
    results: List[Dict[str, Any]] = []
    with input_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"CSV is missing required columns: {sorted(missing)}. "
                f"Found: {reader.fieldnames}"
            )
        for row in reader:
            try:
                results.append(row_to_experimental_result(row, default_timestamp))
            except ValueError as e:
                print(f"Skipping row: {e}", file=sys.stderr)
    return results


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "-i", required=True, type=pathlib.Path)
    parser.add_argument("--output", "-o", required=True, type=pathlib.Path)
    args = parser.parse_args(argv)

    if not args.input.is_file():
        print(f"Input CSV not found: {args.input}", file=sys.stderr)
        return 1

    results = convert(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(results)} experimental result(s) to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
