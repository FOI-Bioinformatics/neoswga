#!/usr/bin/env python3
"""Build a training-data CSV for `train_enhanced_rf.py` from real lab measurements.

Inputs:
    --lab-csv       CSV with measured enrichment (per-primer or per-set rows).
    --fg-genome     Target genome FASTA (single record, used for feature engineering).
    --output        Output CSV path (consumed by scripts/train_enhanced_rf.py).
    --params        Optional params.json for reaction conditions.
    --polymerase, --reaction-temp, --na-conc, --mg-conc, --dmso-percent, --betaine-m
                    Optional overrides for reaction conditions.
    --mode          'per-primer' (default) or 'per-set' (primers semicolon-separated).

Lab CSV schema:
    per-primer:  primer, enrichment_fold[, replicate, ...]
    per-set:     primers, enrichment_fold (primers separated by ';')

Output CSV columns:
    primer, [120+ feature columns from AdvancedFeatureEngineer],
    mean_enrichment, log_enrichment, num_sets

The output is consumed by train_enhanced_rf.py with --target-variable log_enrichment.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd
from Bio import SeqIO

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from neoswga.core.advanced_features import AdvancedFeatureEngineer  # noqa: E402
from neoswga.core.reaction_conditions import ReactionConditions  # noqa: E402
from neoswga.core.thermodynamics import reverse_complement  # noqa: E402

logger = logging.getLogger(__name__)

LABEL_COLUMNS = ("mean_enrichment", "log_enrichment", "num_sets", "set_id", "set_size")

PARAMS_TO_CONDITIONS = {
    "reaction_temp": "temp",
    "na_conc": "na_conc",
    "mg_conc": "mg_conc",
    "polymerase": "polymerase",
    "dmso_percent": "dmso_percent",
    "betaine_m": "betaine_m",
    "trehalose_m": "trehalose_m",
    "formamide_percent": "formamide_percent",
    "glycerol_percent": "glycerol_percent",
    "ethanol_percent": "ethanol_percent",
    "urea_m": "urea_m",
    "tmac_m": "tmac_m",
}


def find_primer_positions(genome: str, primer: str) -> Dict[str, List[int]]:
    """Locate forward and reverse-complement occurrences of a primer."""
    fwd, rev = [], []
    rc = reverse_complement(primer)

    start = 0
    while True:
        pos = genome.find(primer, start)
        if pos == -1:
            break
        fwd.append(pos)
        start = pos + 1

    start = 0
    while True:
        pos = genome.find(rc, start)
        if pos == -1:
            break
        rev.append(pos)
        start = pos + 1

    return {"forward": fwd, "reverse": rev}


def aggregate_labels(lab_df: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Return one row per primer with mean_enrichment and num_sets."""
    if "enrichment_fold" not in lab_df.columns:
        raise ValueError(
            "Lab CSV must contain an 'enrichment_fold' column with measured enrichment."
        )

    if mode == "per-primer":
        if "primer" not in lab_df.columns:
            raise ValueError("Per-primer mode requires a 'primer' column.")
        rows = lab_df[["primer", "enrichment_fold"]].copy()
        rows["primer"] = rows["primer"].str.strip().str.upper()

    elif mode == "per-set":
        if "primers" not in lab_df.columns:
            raise ValueError("Per-set mode requires a 'primers' column (semicolon-separated).")
        exploded = []
        for _, row in lab_df.iterrows():
            primers = [p.strip().upper() for p in str(row["primers"]).split(";") if p.strip()]
            for p in primers:
                exploded.append({"primer": p, "enrichment_fold": float(row["enrichment_fold"])})
        rows = pd.DataFrame(exploded)

    else:
        raise ValueError(f"Unknown mode '{mode}'. Use 'per-primer' or 'per-set'.")

    grouped = rows.groupby("primer", as_index=False).agg(
        mean_enrichment=("enrichment_fold", "mean"),
        num_sets=("enrichment_fold", "count"),
    )
    grouped["log_enrichment"] = np.log10(grouped["mean_enrichment"].clip(lower=1e-3))
    return grouped


def build_training_data(
    lab_csv: Path,
    fg_genome_fasta: Path,
    conditions: ReactionConditions,
    mode: str = "per-primer",
) -> pd.DataFrame:
    """Read lab measurements and the target genome, return a training-ready DataFrame."""
    lab_df = pd.read_csv(lab_csv)
    labels = aggregate_labels(lab_df, mode=mode)
    primers = labels["primer"].tolist()

    record = next(SeqIO.parse(str(fg_genome_fasta), "fasta"))
    genome = str(record.seq).upper()

    positions = {p: find_primer_positions(genome, p) for p in primers}
    engineer = AdvancedFeatureEngineer(genome, conditions, positions)
    features = engineer.engineer_features(primers)

    merged = features.merge(labels, on="primer", how="inner")

    cols = ["primer"] + [
        c for c in merged.columns if c != "primer" and c not in LABEL_COLUMNS
    ] + [c for c in LABEL_COLUMNS if c in merged.columns]
    return merged[cols]


def conditions_from_args(args: argparse.Namespace) -> ReactionConditions:
    """Construct ReactionConditions from --params plus CLI overrides."""
    kwargs: Dict[str, object] = {}

    if args.params:
        with open(args.params) as fh:
            params = json.load(fh)
        for src, dst in PARAMS_TO_CONDITIONS.items():
            if src in params:
                kwargs[dst] = params[src]

    overrides = {
        "temp": args.reaction_temp,
        "polymerase": args.polymerase,
        "na_conc": args.na_conc,
        "mg_conc": args.mg_conc,
        "dmso_percent": args.dmso_percent,
        "betaine_m": args.betaine_m,
    }
    for key, val in overrides.items():
        if val is not None:
            kwargs[key] = val

    return ReactionConditions(**kwargs)


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lab-csv", required=True, type=Path)
    parser.add_argument("--fg-genome", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", choices=("per-primer", "per-set"), default="per-primer")
    parser.add_argument("--params", type=Path, default=None,
                        help="Path to params.json for reaction conditions")
    parser.add_argument("--polymerase", default=None,
                        choices=("phi29", "equiphi29", "bst", "klenow"))
    parser.add_argument("--reaction-temp", type=float, default=None)
    parser.add_argument("--na-conc", type=float, default=None)
    parser.add_argument("--mg-conc", type=float, default=None)
    parser.add_argument("--dmso-percent", type=float, default=None)
    parser.add_argument("--betaine-m", type=float, default=None)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args(list(argv) if argv is not None else None)

    logging.basicConfig(level=args.log_level, format="%(asctime)s %(levelname)s %(message)s")

    conditions = conditions_from_args(args)
    logger.info(
        "Conditions: polymerase=%s temp=%.1fC Na=%.1f mM Mg=%.1f mM",
        conditions.polymerase, conditions.temp, conditions.na_conc, conditions.mg_conc,
    )

    df = build_training_data(
        lab_csv=args.lab_csv,
        fg_genome_fasta=args.fg_genome,
        conditions=conditions,
        mode=args.mode,
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    logger.info("Wrote %d primers x %d columns to %s", len(df), len(df.columns), args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
