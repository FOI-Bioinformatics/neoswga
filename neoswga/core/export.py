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


# Protocol template for wet-lab use
PROTOCOL_TEMPLATE = """# SWGA Protocol

**Generated:** {date}
**Polymerase:** {polymerase}
**Primers:** {num_primers}

---

## Primer Set

| # | Sequence | Length | Tm | GC |
|---|----------|--------|-----|-----|
{primer_table}

---

## Reaction Conditions

- **Temperature:** {temperature}C
- **Duration:** {duration} hours
- **Mg2+ concentration:** {mg_conc} mM
{additives_section}

---

## Reaction Setup (25 uL)

| Component | Volume | Final Concentration |
|-----------|--------|---------------------|
| 10X Reaction Buffer | 2.5 uL | 1X |
| dNTPs (10 mM each) | 1.0 uL | 0.4 mM each |
| Primer Mix (10 uM each) | 0.5 uL | 200 nM each |
| {polymerase} Polymerase | 0.5 uL | As per manufacturer |
| Template DNA | 1.0 uL | 1-10 ng |
| Nuclease-free water | 19.5 uL | - |

---

## Procedure

1. Thaw all components on ice
2. Prepare master mix (all components except template)
3. Aliquot master mix into tubes
4. Add template DNA
5. Incubate at {temperature}C for {duration} hours
6. Heat-inactivate at 65C for 10 minutes
7. Store product at -20C

---

## Quality Control

- Run 1 uL on 0.8% agarose gel
- Expected: High molecular weight smear (>10 kb)
- Include no-template control (NTC)

---

## Primer Ordering Notes

- Synthesis scale: 25 nmol (standard desalt purification)
- Resuspension: 100 uM in TE buffer (10 mM Tris, 0.1 mM EDTA, pH 8.0)
- Storage: -20C, protect from light
- Working stock: Dilute to 10 uM for primer mix preparation

---

## Primer Mix Preparation

Combine equal volumes of each primer at 10 uM to create working primer mix.
Final concentration in reaction: 200 nM each primer.

---

*Protocol generated by NeoSWGA*
"""


def generate_protocol(
    primers: List[str],
    polymerase: str = "phi29",
    temperature: float = 30.0,
    duration: float = 16.0,
    mg_conc: float = 2.5,
    betaine_m: float = 0.0,
    dmso_percent: float = 0.0,
    trehalose_m: float = 0.0,
) -> str:
    """
    Generate wet-lab protocol markdown.

    Args:
        primers: List of primer sequences.
        polymerase: Polymerase name.
        temperature: Reaction temperature (C).
        duration: Reaction duration (hours).
        mg_conc: Mg2+ concentration (mM).
        betaine_m: Betaine concentration (M).
        dmso_percent: DMSO percentage.
        trehalose_m: Trehalose concentration (M).

    Returns:
        Protocol as markdown string.
    """
    # Build primer table
    table_rows = []
    for i, primer in enumerate(primers, 1):
        tm = calculate_simple_tm(primer)
        gc = calculate_gc(primer)
        table_rows.append(
            f"| {i} | {primer} | {len(primer)} | {tm:.1f}C | {gc:.1%} |"
        )
    primer_table = "\n".join(table_rows)

    # Build additives section
    additives = []
    if betaine_m > 0:
        additives.append(f"- **Betaine:** {betaine_m} M")
    if dmso_percent > 0:
        additives.append(f"- **DMSO:** {dmso_percent}%")
    if trehalose_m > 0:
        additives.append(f"- **Trehalose:** {trehalose_m} M")

    additives_section = "\n".join(additives) if additives else "- No additives"

    return PROTOCOL_TEMPLATE.format(
        date=datetime.now().strftime("%Y-%m-%d"),
        polymerase=polymerase,
        num_primers=len(primers),
        primer_table=primer_table,
        temperature=temperature,
        duration=duration,
        mg_conc=mg_conc,
        additives_section=additives_section,
    )


def export_protocol(
    primers: List[str],
    output_path: str,
    **kwargs
) -> None:
    """
    Export protocol to markdown file.

    Args:
        primers: List of primer sequences.
        output_path: Path for output file.
        **kwargs: Arguments passed to generate_protocol.
    """
    protocol = generate_protocol(primers, **kwargs)

    with open(output_path, 'w') as f:
        f.write(protocol)

    logger.info(f"Exported protocol to {output_path}")
