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
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import List, Dict, Optional, Any
from datetime import datetime

logger = logging.getLogger(__name__)


class ModificationProfile(Enum):
    """Modification profiles for primer synthesis.

    SWGA primers benefit from chemical modifications for Phi29 compatibility:
    - 3' Phosphorothioate (PTO) bonds protect against exonuclease degradation
    - 5' C18 spacer blocks template-independent amplification (tdMDA approach)

    Reference: Wang et al. (2017) BioTechniques 63:21-26
    """
    NONE = "none"
    STANDARD = "standard"
    LOW_INPUT = "low-input"


@dataclass
class PrimerModifications:
    """Primer modification configuration.

    Attributes:
        pto_bonds: Number of 3' phosphorothioate bonds (default: 2)
        five_prime_block: 5' blocking group ('c18', 'c3', or None)
        profile: The modification profile used

    Example:
        Standard: ATCGATCGAT*C*G (2 PTO bonds at 3' end)
        Low-input: /5SpC18/ATCGATCGAT*C*G (C18 spacer + PTO)
    """
    pto_bonds: int = 2
    five_prime_block: Optional[str] = None  # 'c18', 'c3', or None
    profile: ModificationProfile = ModificationProfile.STANDARD

    @classmethod
    def from_profile(cls, profile: str) -> 'PrimerModifications':
        """Create modification config from profile name.

        Args:
            profile: Profile name ('none', 'standard', 'low-input')

        Returns:
            PrimerModifications instance

        Raises:
            ValueError: If profile name is unknown
        """
        profiles = {
            'none': cls(
                pto_bonds=0,
                five_prime_block=None,
                profile=ModificationProfile.NONE
            ),
            'standard': cls(
                pto_bonds=2,
                five_prime_block=None,
                profile=ModificationProfile.STANDARD
            ),
            'low-input': cls(
                pto_bonds=2,
                five_prime_block='c18',
                profile=ModificationProfile.LOW_INPUT
            ),
        }
        if profile not in profiles:
            raise ValueError(
                f"Unknown modification profile: '{profile}'. "
                f"Available profiles: {', '.join(profiles.keys())}"
            )
        return profiles[profile]


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


def export_to_bed(
    primers: List[str],
    positions: Dict[str, List[tuple]],
    genome_name: str,
    output_path: str,
) -> None:
    """Export primer binding sites in BED format for genome browser visualization.

    BED format (6-column):
        chrom  start  end  name  score  strand

    Args:
        primers: List of primer sequences.
        positions: Dict mapping primer to list of (position, strand) tuples
            where strand is 'forward' or 'reverse'.
        genome_name: Chromosome or contig name for BED output.
        output_path: Path for output file.
    """
    with open(output_path, 'w') as f:
        for primer in primers:
            sites = positions.get(primer, [])
            for pos, strand in sites:
                strand_char = '+' if strand == 'forward' else '-'
                end = pos + len(primer)
                f.write(
                    f"{genome_name}\t{pos}\t{end}\t{primer}\t0\t{strand_char}\n"
                )

    total_sites = sum(len(positions.get(p, [])) for p in primers)
    logger.info(f"Exported {total_sites} binding sites to {output_path}")


def export_to_bedgraph(
    primers: List[str],
    positions: Dict[str, List[tuple]],
    genome_name: str,
    genome_length: int,
    output_path: str,
    window_size: int = 1000,
) -> None:
    """Export primer binding density in BedGraph format.

    BedGraph shows coverage depth per genomic window, suitable for
    visualizing primer density across the genome.

    Args:
        primers: List of primer sequences.
        positions: Dict mapping primer to list of (position, strand) tuples.
        genome_name: Chromosome or contig name.
        genome_length: Total genome length in bp.
        output_path: Path for output file.
        window_size: Window size in bp for coverage binning (default: 1000).
    """
    import numpy as np

    num_windows = (genome_length + window_size - 1) // window_size
    counts = np.zeros(num_windows, dtype=int)

    for primer in primers:
        for pos, _strand in positions.get(primer, []):
            window_idx = pos // window_size
            if 0 <= window_idx < num_windows:
                counts[window_idx] += 1

    with open(output_path, 'w') as f:
        for i in range(num_windows):
            if counts[i] > 0:
                start = i * window_size
                end = min((i + 1) * window_size, genome_length)
                f.write(f"{genome_name}\t{start}\t{end}\t{counts[i]}\n")

    logger.info(f"Exported BedGraph with {window_size}bp windows to {output_path}")


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

# Vendor-specific modification syntax
# Reference: IDT ordering guide, Sigma custom oligo guide
VENDOR_MODIFICATION_SYNTAX: Dict[str, Dict[str, Any]] = {
    "idt": {
        "pto_symbol": "*",
        "c18_prefix": "/5SpC18/",
        "c3_prefix": "/5SpC3/",
        "supports_mods": True,
        "min_scale_for_mods": "100nm",
        "notes": "Modifications require 100nm scale minimum",
    },
    "twist": {
        "pto_symbol": "*",
        "c18_prefix": "",
        "c3_prefix": "",
        "supports_mods": False,
        "notes": "Contact Twist directly for modified oligos",
    },
    "sigma": {
        "pto_symbol": "*",
        "c18_prefix": "/5SpC18/",
        "c3_prefix": "/5SpC3/",
        "supports_mods": True,
        "min_scale_for_mods": "0.1 umol",
        "notes": "Uses IDT-compatible syntax",
    },
    "generic": {
        "pto_symbol": "*",
        "c18_prefix": "[C18]-",
        "c3_prefix": "[C3]-",
        "supports_mods": True,
        "notes": "Generic notation for documentation",
    },
}


def apply_modifications(
    sequence: str,
    modifications: PrimerModifications,
    vendor: str = "generic"
) -> str:
    """Apply chemical modifications to primer sequence using vendor syntax.

    Modifications applied:
    - 3' Phosphorothioate (PTO) bonds: Protect against Phi29 3'->5' exonuclease
    - 5' C18 spacer: Block template-independent amplification (tdMDA)

    Args:
        sequence: Raw primer sequence (e.g., 'ATCGATCGATCG')
        modifications: PrimerModifications configuration
        vendor: Target vendor for syntax ('idt', 'twist', 'sigma', 'generic')

    Returns:
        Modified sequence string (e.g., '/5SpC18/ATCGATCGAT*C*G')

    Example:
        >>> mods = PrimerModifications.from_profile('standard')
        >>> apply_modifications('ATCGATCG', mods, vendor='idt')
        'ATCGAT*C*G'

        >>> mods = PrimerModifications.from_profile('low-input')
        >>> apply_modifications('ATCGATCG', mods, vendor='idt')
        '/5SpC18/ATCGAT*C*G'
    """
    if modifications.profile == ModificationProfile.NONE:
        return sequence

    vendor = vendor.lower()
    syntax = VENDOR_MODIFICATION_SYNTAX.get(
        vendor, VENDOR_MODIFICATION_SYNTAX["generic"]
    )

    # Warn if vendor doesn't support modifications
    if not syntax.get("supports_mods", True):
        if modifications.pto_bonds > 0 or modifications.five_prime_block:
            logger.warning(
                f"Vendor '{vendor}' may not support modifications directly. "
                f"{syntax.get('notes', 'Contact vendor for custom synthesis.')}"
            )

    result = sequence

    # Apply 3' PTO bonds
    # PTO bonds are inserted between the last N nucleotides
    # Example: ATCGATCG with 2 PTO -> ATCGAT*C*G
    if modifications.pto_bonds > 0:
        pto = syntax["pto_symbol"]
        n = modifications.pto_bonds
        if len(result) > n:
            # Split into prefix and the last n+1 nucleotides
            prefix = result[:-n-1]
            suffix = result[-n-1:]
            # Insert PTO between each of the last n+1 nucleotides
            modified_suffix = pto.join(suffix)
            result = prefix + modified_suffix
        elif len(result) > 1:
            # Short sequence: insert PTO between all nucleotides
            result = pto.join(result)

    # Apply 5' blocking group
    if modifications.five_prime_block == 'c18':
        prefix = syntax.get("c18_prefix", "[C18]-")
        if prefix:  # Only add if vendor supports it
            result = prefix + result
    elif modifications.five_prime_block == 'c3':
        prefix = syntax.get("c3_prefix", "[C3]-")
        if prefix:
            result = prefix + result

    return result


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


@dataclass
class ExportSummary:
    """Summary of exported primer set."""
    num_primers: int
    mean_length: float
    mean_gc: float
    mean_tm: float
    min_tm: float
    max_tm: float
    estimated_cost: float  # USD


class PrimerExporter:
    """
    Unified primer export handler.

    Provides convenient methods to export primers in multiple formats
    and generate ordering summaries.

    Usage:
        # From primer list
        exporter = PrimerExporter(["ATCGATCG", "GCTAGCTA"])

        # From results directory
        exporter = PrimerExporter.from_results_dir("./results/")

        # Export all formats
        exporter.export_all("./export/", project_name="MyProject")

        # Or individually
        exporter.export_fasta("primers.fasta")
        exporter.export_vendor_csv("order.csv", vendor="idt")
        exporter.export_protocol("protocol.md")
    """

    def __init__(
        self,
        primers: List[str],
        polymerase: str = "phi29",
        temperature: float = 30.0,
        betaine_m: float = 0.0,
        dmso_percent: float = 0.0,
        trehalose_m: float = 0.0,
        mg_conc: float = 2.5,
        modifications: Optional[PrimerModifications] = None,
    ):
        """
        Initialize exporter with primers and reaction conditions.

        Args:
            primers: List of primer sequences
            polymerase: Polymerase type
            temperature: Reaction temperature (C)
            betaine_m: Betaine concentration (M)
            dmso_percent: DMSO percentage
            trehalose_m: Trehalose concentration (M)
            mg_conc: Mg2+ concentration (mM)
            modifications: Primer modification config (default: standard profile)
        """
        self.primers = primers
        self.polymerase = polymerase
        self.temperature = temperature
        self.betaine_m = betaine_m
        self.dmso_percent = dmso_percent
        self.trehalose_m = trehalose_m
        self.mg_conc = mg_conc
        self.modifications = modifications or PrimerModifications.from_profile('standard')

    @classmethod
    def from_results_dir(
        cls,
        results_dir: str,
        params_file: Optional[str] = None
    ) -> 'PrimerExporter':
        """
        Create exporter from pipeline results directory.

        Args:
            results_dir: Path to results directory containing step4_improved_df.csv
            params_file: Optional path to params.json for reaction conditions

        Returns:
            PrimerExporter instance
        """
        import json

        results_path = Path(results_dir)

        # Load primers from step4 output
        step4_file = results_path / "step4_improved_df.csv"
        if not step4_file.exists():
            raise FileNotFoundError(f"Step 4 results not found: {step4_file}")

        primers = []
        with open(step4_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                primer = row.get('primer', row.get('sequence', ''))
                if primer and primer not in primers:
                    primers.append(primer)

        if not primers:
            raise ValueError("No primers found in step4 results")

        # Load reaction conditions from params.json if available
        kwargs = {}
        params_path = params_file or (results_path / "params.json")
        if Path(params_path).exists():
            with open(params_path, 'r') as f:
                params = json.load(f)
                kwargs['polymerase'] = params.get('polymerase', 'phi29')
                kwargs['temperature'] = params.get('reaction_temp', 30.0)
                kwargs['betaine_m'] = params.get('betaine_m', 0.0)
                kwargs['dmso_percent'] = params.get('dmso_percent', 0.0)
                kwargs['trehalose_m'] = params.get('trehalose_m', 0.0)
                kwargs['mg_conc'] = params.get('mg_conc', 2.5)

        return cls(primers, **kwargs)

    def export_fasta(self, output_path: str, **kwargs) -> None:
        """Export primers to FASTA format."""
        export_to_fasta(self.primers, output_path, **kwargs)

    def export_vendor_csv(
        self, output_path: str, vendor: str = "idt", **kwargs
    ) -> None:
        """Export primers in vendor-specific CSV format with modifications applied."""
        # Apply modifications to primers
        modified_primers = [
            apply_modifications(p, self.modifications, vendor)
            for p in self.primers
        ]
        export_to_vendor_csv(modified_primers, output_path, vendor=vendor, **kwargs)

    def export_protocol(self, output_path: str, **kwargs) -> None:
        """Export wet-lab protocol."""
        # Use instance conditions as defaults
        protocol_kwargs = {
            'polymerase': self.polymerase,
            'temperature': self.temperature,
            'betaine_m': self.betaine_m,
            'dmso_percent': self.dmso_percent,
            'trehalose_m': self.trehalose_m,
            'mg_conc': self.mg_conc,
        }
        protocol_kwargs.update(kwargs)
        export_protocol(self.primers, output_path, **protocol_kwargs)

    def export_bed(
        self,
        output_path: str,
        positions: Dict[str, List[tuple]],
        genome_name: str = "genome",
    ) -> None:
        """Export primer binding sites in BED format."""
        export_to_bed(self.primers, positions, genome_name, output_path)

    def export_bedgraph(
        self,
        output_path: str,
        positions: Dict[str, List[tuple]],
        genome_name: str = "genome",
        genome_length: int = 0,
        window_size: int = 1000,
    ) -> None:
        """Export primer binding density in BedGraph format."""
        export_to_bedgraph(
            self.primers, positions, genome_name,
            genome_length, output_path, window_size
        )

    def export_all(
        self,
        output_dir: str,
        project_name: str = "SWGA",
        vendors: Optional[List[str]] = None
    ) -> Dict[str, str]:
        """
        Export all formats to a directory.

        Args:
            output_dir: Output directory path
            project_name: Project name for file naming
            vendors: List of vendors for CSV export (default: ["idt"])

        Returns:
            Dictionary mapping format to output path
        """
        if vendors is None:
            vendors = ["idt"]

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        outputs = {}

        # FASTA
        fasta_path = output_path / f"{project_name}_primers.fasta"
        self.export_fasta(
            str(fasta_path), prefix=project_name, include_metadata=True
        )
        outputs["fasta"] = str(fasta_path)

        # Vendor CSVs
        for vendor in vendors:
            csv_path = output_path / f"{project_name}_order_{vendor}.csv"
            self.export_vendor_csv(
                str(csv_path), vendor=vendor, project_name=project_name
            )
            outputs[f"csv_{vendor}"] = str(csv_path)

        # Protocol
        protocol_path = output_path / f"{project_name}_protocol.md"
        self.export_protocol(str(protocol_path))
        outputs["protocol"] = str(protocol_path)

        # BED and BedGraph (require position data — skip if unavailable)
        try:
            if hasattr(self, '_positions') and self._positions:
                bed_path = output_path / f"{project_name}_primers.bed"
                genome_name = getattr(self, '_genome_name', 'genome')
                self.export_bed(str(bed_path), self._positions, genome_name)
                outputs["bed"] = str(bed_path)

                bedgraph_path = output_path / f"{project_name}_coverage.bedgraph"
                genome_length = getattr(self, '_genome_length', 0)
                self.export_bedgraph(
                    str(bedgraph_path), self._positions, genome_name,
                    genome_length
                )
                outputs["bedgraph"] = str(bedgraph_path)
        except Exception as e:
            logger.debug(f"BED/BedGraph export skipped (no position data): {e}")

        logger.info(f"Exported all formats to {output_dir}")
        return outputs

    def get_summary(self) -> Dict[str, Any]:
        """
        Get summary statistics for the primer set.

        Returns:
            Dictionary with summary statistics
        """
        if not self.primers:
            return {"num_primers": 0}

        lengths = [len(p) for p in self.primers]
        gcs = [calculate_gc(p) for p in self.primers]
        tms = [calculate_simple_tm(p) for p in self.primers]

        # Estimate cost: approximately $5 per primer at 25nmol scale
        estimated_cost = len(self.primers) * 5.0

        return {
            "num_primers": len(self.primers),
            "mean_length": sum(lengths) / len(lengths),
            "mean_gc": sum(gcs) / len(gcs),
            "mean_tm": sum(tms) / len(tms),
            "min_tm": min(tms),
            "max_tm": max(tms),
            "estimated_cost": estimated_cost,
            "polymerase": self.polymerase,
            "temperature": self.temperature,
            "modification_profile": self.modifications.profile.value,
            "pto_bonds": self.modifications.pto_bonds,
            "five_prime_block": self.modifications.five_prime_block,
        }

    def print_summary(self) -> None:
        """Print formatted summary to console."""
        summary = self.get_summary()

        print("\n" + "=" * 50)
        print("PRIMER SET SUMMARY")
        print("=" * 50)
        print(f"  Primers: {summary['num_primers']}")
        print(f"  Mean length: {summary['mean_length']:.1f} bp")
        print(f"  Mean GC: {summary['mean_gc']:.1%}")
        print(f"  Tm range: {summary['min_tm']:.1f} - {summary['max_tm']:.1f} C")
        print(f"  Polymerase: {summary['polymerase']}")
        print(f"  Temperature: {summary['temperature']} C")
        print(f"  Estimated cost: ${summary['estimated_cost']:.2f}")
        print("-" * 50)
        print("  MODIFICATIONS")
        print(f"  Profile: {summary['modification_profile']}")
        print(f"  3' PTO bonds: {summary['pto_bonds']}")
        block = summary['five_prime_block'] or 'none'
        print(f"  5' blocking: {block}")
        print("=" * 50 + "\n")
