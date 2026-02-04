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
