"""
Primer export functionality for synthesis ordering.

Generates multiple output formats:
- FASTA: Standard sequence format
- CSV: Vendor-ready ordering format (IDT, Twist, Sigma)
- Protocol: Wet-lab protocol markdown

Usage:
    from neoswga.core.export import PrimerExporter

    exporter = PrimerExporter(primers, params)
    exporter.export_fasta("primers.fasta")
    exporter.export_vendor_csv("order.csv", vendor="idt")
    exporter.export_protocol("protocol.md")
"""

import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional, Any
from datetime import datetime

logger = logging.getLogger(__name__)


def calculate_gc(seq: str) -> float:
    """Calculate GC content as fraction."""
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq) if len(seq) > 0 else 0.0


def calculate_simple_tm(seq: str) -> float:
    """Calculate Tm using Wallace rule (for short primers)."""
    seq = seq.upper()
    at_count = seq.count('A') + seq.count('T')
    gc_count = seq.count('G') + seq.count('C')
    return 2 * at_count + 4 * gc_count


def export_to_fasta(
    primers: List[str],
    output_path: str,
    prefix: str = "SWGA",
    include_metadata: bool = False
) -> None:
    """
    Export primers to FASTA format.

    Args:
        primers: List of primer sequences
        output_path: Path for output file
        prefix: Prefix for sequence names (default: SWGA)
        include_metadata: Include Tm and GC in header
    """
    with open(output_path, 'w') as f:
        for i, primer in enumerate(primers, 1):
            header = f">{prefix}_{i:03d}"

            if include_metadata:
                tm = calculate_simple_tm(primer)
                gc = calculate_gc(primer)
                header += f" Tm={tm:.1f}C GC={gc:.1%} len={len(primer)}"

            f.write(f"{header}\n{primer}\n")

    logger.info(f"Exported {len(primers)} primers to {output_path}")


# Vendor format specifications
VENDOR_FORMATS: Dict[str, Dict[str, Any]] = {
    "idt": {
        "columns": ["Name", "Sequence", "Scale", "Purification"],
        "defaults": {"Scale": "25nm", "Purification": "STD"},
    },
    "twist": {
        "columns": ["Name", "Sequence"],
        "defaults": {},
    },
    "sigma": {
        "columns": ["Oligo Name", "Sequence (5' to 3')", "Scale", "Purification"],
        "defaults": {"Scale": "0.025 umol", "Purification": "Desalt"},
    },
    "generic": {
        "columns": ["name", "sequence", "length", "tm", "gc"],
        "defaults": {},
    },
}


def export_to_vendor_csv(
    primers: List[str],
    output_path: str,
    vendor: str = "generic",
    project_name: str = "SWGA",
    scale: Optional[str] = None,
    purification: Optional[str] = None
) -> None:
    """
    Export primers in vendor-specific CSV format.

    Args:
        primers: List of primer sequences.
        output_path: Path for output file.
        vendor: Vendor name ('idt', 'twist', 'sigma', 'generic').
        project_name: Project/order name prefix.
        scale: Override default synthesis scale.
        purification: Override default purification method.
    """
    vendor = vendor.lower()
    if vendor not in VENDOR_FORMATS:
        logger.warning(f"Unknown vendor '{vendor}', using generic format")
        vendor = "generic"

    fmt = VENDOR_FORMATS[vendor]
    columns = fmt["columns"]
    defaults = fmt["defaults"].copy()

    if scale:
        for col in ["Scale", "scale"]:
            if col in columns:
                defaults[col] = scale
    if purification:
        for col in ["Purification", "purification"]:
            if col in columns:
                defaults[col] = purification

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()

        for i, primer in enumerate(primers, 1):
            row = defaults.copy()
            name = f"{project_name}_{i:03d}"

            # Map to vendor-specific column names
            if "Name" in columns:
                row["Name"] = name
            if "Oligo Name" in columns:
                row["Oligo Name"] = name
            if "name" in columns:
                row["name"] = name

            if "Sequence" in columns:
                row["Sequence"] = primer
            if "Sequence (5' to 3')" in columns:
                row["Sequence (5' to 3')"] = primer
            if "sequence" in columns:
                row["sequence"] = primer

            if "length" in columns:
                row["length"] = len(primer)
            if "tm" in columns:
                row["tm"] = f"{calculate_simple_tm(primer):.1f}"
            if "gc" in columns:
                row["gc"] = f"{calculate_gc(primer):.1%}"

            writer.writerow(row)

    logger.info(f"Exported {len(primers)} primers to {output_path} ({vendor} format)")
