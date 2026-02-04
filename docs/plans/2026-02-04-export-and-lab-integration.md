# Export and Lab Integration Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Bridge the gap between pipeline output and laboratory use by adding export functionality and practical lab guidance.

**Architecture:** Add a new `export` module under `neoswga/core/` that handles primer export in multiple formats (FASTA, vendor CSV, protocol). Integrate export into existing report generation and add a standalone `neoswga export` CLI command. Create comprehensive documentation for the results-to-lab workflow.

**Tech Stack:** Python dataclasses, Jinja2-style string templating (already in report module), CSV/FASTA file generation, existing report infrastructure.

---

## Overview

This plan addresses 5 critical gaps:
1. No final primer export file
2. No vendor-ready ordering format
3. Protocol generation not integrated
4. Documentation gap between software and lab
5. No synthesis guidance

**Deliverables:**
- `neoswga/core/export.py` - Export module with FASTA, CSV, protocol generation
- `neoswga export` CLI command
- Enhanced report with primer sequences section
- `docs/FROM_RESULTS_TO_LAB.md` - Complete lab workflow guide
- Updated documentation index

---

## Task 1: Create Export Module Core

**Files:**
- Create: `neoswga/core/export.py`
- Test: `tests/test_export.py`

**Step 1: Write the failing test for FASTA export**

```python
# tests/test_export.py
"""Tests for primer export functionality."""

import pytest
from pathlib import Path
import tempfile


class TestFastaExport:
    """Test FASTA format export."""

    def test_export_primers_to_fasta(self):
        """Export primer list to FASTA format."""
        from neoswga.core.export import export_to_fasta

        primers = ["ATCGATCG", "GCTAGCTA", "TTAATTAA"]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            output_path = f.name

        export_to_fasta(primers, output_path, prefix="SWGA")

        content = Path(output_path).read_text()

        assert ">SWGA_001" in content
        assert "ATCGATCG" in content
        assert ">SWGA_002" in content
        assert ">SWGA_003" in content

        Path(output_path).unlink()

    def test_fasta_with_metadata(self):
        """FASTA export includes Tm and GC in header."""
        from neoswga.core.export import export_to_fasta

        primers = ["ATCGATCGATCG"]  # 12-mer

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            output_path = f.name

        export_to_fasta(primers, output_path, include_metadata=True)

        content = Path(output_path).read_text()

        # Should include Tm and GC in header
        assert "Tm=" in content
        assert "GC=" in content

        Path(output_path).unlink()
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_export.py::TestFastaExport::test_export_primers_to_fasta -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'neoswga.core.export'`

**Step 3: Write minimal implementation**

```python
# neoswga/core/export.py
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
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_export.py::TestFastaExport -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add neoswga/core/export.py tests/test_export.py
git commit -m "feat(export): add FASTA export for primer sequences"
```

---

## Task 2: Add Vendor CSV Export

**Files:**
- Modify: `neoswga/core/export.py`
- Test: `tests/test_export.py`

**Step 1: Write the failing test for vendor CSV**

```python
# Add to tests/test_export.py

class TestVendorCsvExport:
    """Test vendor-ready CSV export."""

    def test_export_idt_format(self):
        """Export in IDT ordering format."""
        from neoswga.core.export import export_to_vendor_csv

        primers = ["ATCGATCG", "GCTAGCTA"]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_path = f.name

        export_to_vendor_csv(primers, output_path, vendor="idt", project_name="MyProject")

        content = Path(output_path).read_text()

        # IDT format has Name, Sequence, Scale, Purification
        assert "Name" in content
        assert "Sequence" in content
        assert "MyProject_001" in content
        assert "ATCGATCG" in content
        assert "25nm" in content.lower() or "25 nm" in content.lower()

        Path(output_path).unlink()

    def test_export_twist_format(self):
        """Export in Twist Bioscience format."""
        from neoswga.core.export import export_to_vendor_csv

        primers = ["ATCGATCGATCG"]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_path = f.name

        export_to_vendor_csv(primers, output_path, vendor="twist")

        content = Path(output_path).read_text()

        # Twist uses different column names
        assert "Sequence" in content or "sequence" in content

        Path(output_path).unlink()

    def test_export_generic_format(self):
        """Export in generic CSV format."""
        from neoswga.core.export import export_to_vendor_csv

        primers = ["ATCGATCG"]

        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            output_path = f.name

        export_to_vendor_csv(primers, output_path, vendor="generic")

        content = Path(output_path).read_text()

        assert "name" in content.lower()
        assert "sequence" in content.lower()
        assert "ATCGATCG" in content

        Path(output_path).unlink()
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_export.py::TestVendorCsvExport -v
```

Expected: FAIL with `ImportError`

**Step 3: Add vendor CSV implementation**

```python
# Add to neoswga/core/export.py

# Vendor format specifications
VENDOR_FORMATS = {
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
        primers: List of primer sequences
        output_path: Path for output file
        vendor: Vendor name ('idt', 'twist', 'sigma', 'generic')
        project_name: Project/order name prefix
        scale: Override default synthesis scale
        purification: Override default purification method
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
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_export.py::TestVendorCsvExport -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add neoswga/core/export.py tests/test_export.py
git commit -m "feat(export): add vendor-specific CSV export (IDT, Twist, Sigma)"
```

---

## Task 3: Add Protocol Generation

**Files:**
- Modify: `neoswga/core/export.py`
- Test: `tests/test_export.py`

**Step 1: Write the failing test**

```python
# Add to tests/test_export.py

class TestProtocolExport:
    """Test wet-lab protocol generation."""

    def test_generate_basic_protocol(self):
        """Generate basic SWGA protocol."""
        from neoswga.core.export import generate_protocol

        primers = ["ATCGATCG", "GCTAGCTA"]

        protocol = generate_protocol(
            primers,
            polymerase="phi29",
            temperature=30.0,
        )

        assert "SWGA Protocol" in protocol
        assert "phi29" in protocol.lower() or "Phi29" in protocol
        assert "30" in protocol
        assert "ATCGATCG" in protocol
        assert "GCTAGCTA" in protocol

    def test_protocol_with_additives(self):
        """Protocol includes additive recommendations."""
        from neoswga.core.export import generate_protocol

        primers = ["ATCGATCGATCG"]

        protocol = generate_protocol(
            primers,
            polymerase="equiphi29",
            temperature=42.0,
            betaine_m=1.0,
            dmso_percent=5.0,
        )

        assert "equiphi29" in protocol.lower() or "EquiPhi29" in protocol
        assert "42" in protocol
        assert "betaine" in protocol.lower()
        assert "DMSO" in protocol or "dmso" in protocol

    def test_protocol_includes_primer_table(self):
        """Protocol includes formatted primer table."""
        from neoswga.core.export import generate_protocol

        primers = ["ATCGATCG", "GCTAGCTA", "TTAATTAA"]

        protocol = generate_protocol(primers, polymerase="phi29")

        # Should have a table with all primers
        for primer in primers:
            assert primer in protocol
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_export.py::TestProtocolExport -v
```

Expected: FAIL with `ImportError`

**Step 3: Add protocol generation implementation**

```python
# Add to neoswga/core/export.py

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
        primers: List of primer sequences
        polymerase: Polymerase name
        temperature: Reaction temperature (C)
        duration: Reaction duration (hours)
        mg_conc: Mg2+ concentration (mM)
        betaine_m: Betaine concentration (M)
        dmso_percent: DMSO percentage
        trehalose_m: Trehalose concentration (M)

    Returns:
        Protocol as markdown string
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
        primers: List of primer sequences
        output_path: Path for output file
        **kwargs: Arguments passed to generate_protocol
    """
    protocol = generate_protocol(primers, **kwargs)

    with open(output_path, 'w') as f:
        f.write(protocol)

    logger.info(f"Exported protocol to {output_path}")
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_export.py::TestProtocolExport -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add neoswga/core/export.py tests/test_export.py
git commit -m "feat(export): add wet-lab protocol generation"
```

---

## Task 4: Add PrimerExporter Class

**Files:**
- Modify: `neoswga/core/export.py`
- Test: `tests/test_export.py`

**Step 1: Write the failing test**

```python
# Add to tests/test_export.py

class TestPrimerExporter:
    """Test unified PrimerExporter class."""

    def test_exporter_from_results_dir(self, tmp_path):
        """Create exporter from results directory."""
        from neoswga.core.export import PrimerExporter

        # Create mock step4 file
        step4_path = tmp_path / "step4_improved_df.csv"
        step4_path.write_text(
            "primer,set_index,score,coverage\n"
            "ATCGATCG,0,0.85,0.75\n"
            "GCTAGCTA,0,0.85,0.75\n"
        )

        exporter = PrimerExporter.from_results_dir(str(tmp_path))

        assert len(exporter.primers) == 2
        assert "ATCGATCG" in exporter.primers
        assert "GCTAGCTA" in exporter.primers

    def test_exporter_export_all(self, tmp_path):
        """Export all formats at once."""
        from neoswga.core.export import PrimerExporter

        primers = ["ATCGATCG", "GCTAGCTA"]
        exporter = PrimerExporter(primers)

        output_dir = tmp_path / "export"
        exporter.export_all(str(output_dir), project_name="TestProject")

        assert (output_dir / "TestProject_primers.fasta").exists()
        assert (output_dir / "TestProject_order_idt.csv").exists()
        assert (output_dir / "TestProject_protocol.md").exists()

    def test_exporter_summary(self):
        """Get export summary."""
        from neoswga.core.export import PrimerExporter

        primers = ["ATCGATCG", "GCTAGCTA", "TTAATTAA"]
        exporter = PrimerExporter(primers)

        summary = exporter.get_summary()

        assert summary["num_primers"] == 3
        assert "mean_length" in summary
        assert "mean_gc" in summary
        assert "estimated_cost" in summary
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_export.py::TestPrimerExporter -v
```

Expected: FAIL with `ImportError`

**Step 3: Add PrimerExporter class**

```python
# Add to neoswga/core/export.py

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
        """
        self.primers = primers
        self.polymerase = polymerase
        self.temperature = temperature
        self.betaine_m = betaine_m
        self.dmso_percent = dmso_percent
        self.trehalose_m = trehalose_m
        self.mg_conc = mg_conc

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

    def export_vendor_csv(self, output_path: str, vendor: str = "idt", **kwargs) -> None:
        """Export primers in vendor-specific CSV format."""
        export_to_vendor_csv(self.primers, output_path, vendor=vendor, **kwargs)

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

    def export_all(
        self,
        output_dir: str,
        project_name: str = "SWGA",
        vendors: List[str] = None
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
        self.export_fasta(str(fasta_path), prefix=project_name, include_metadata=True)
        outputs["fasta"] = str(fasta_path)

        # Vendor CSVs
        for vendor in vendors:
            csv_path = output_path / f"{project_name}_order_{vendor}.csv"
            self.export_vendor_csv(str(csv_path), vendor=vendor, project_name=project_name)
            outputs[f"csv_{vendor}"] = str(csv_path)

        # Protocol
        protocol_path = output_path / f"{project_name}_protocol.md"
        self.export_protocol(str(protocol_path))
        outputs["protocol"] = str(protocol_path)

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

        # Estimate cost: ~$5 per primer at 25nmol scale
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
        print("=" * 50 + "\n")
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_export.py::TestPrimerExporter -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add neoswga/core/export.py tests/test_export.py
git commit -m "feat(export): add PrimerExporter class for unified export"
```

---

## Task 5: Add CLI Export Command

**Files:**
- Modify: `neoswga/cli_unified.py`
- Test: `tests/test_export.py`

**Step 1: Write the failing test**

```python
# Add to tests/test_export.py

class TestExportCLI:
    """Test CLI export command."""

    def test_export_command_exists(self):
        """Export command is registered in CLI."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        # Parse help to check subcommands
        import io
        import sys

        # Check that export subparser exists by trying to parse it
        try:
            args = parser.parse_args(['export', '--help'])
        except SystemExit:
            pass  # --help causes SystemExit, that's OK

    def test_export_from_results_dir(self, tmp_path):
        """CLI exports from results directory."""
        # Create mock results
        step4_path = tmp_path / "step4_improved_df.csv"
        step4_path.write_text(
            "primer,set_index,score,coverage\n"
            "ATCGATCG,0,0.85,0.75\n"
            "GCTAGCTA,0,0.85,0.75\n"
        )

        from neoswga.cli_unified import create_parser

        output_dir = tmp_path / "export"

        parser = create_parser()
        args = parser.parse_args([
            'export',
            '-d', str(tmp_path),
            '-o', str(output_dir),
            '--project', 'TestProject'
        ])

        assert args.dir == str(tmp_path)
        assert args.output == str(output_dir)
        assert args.project == 'TestProject'
```

**Step 2: Run test to verify it fails**

```bash
pytest tests/test_export.py::TestExportCLI -v
```

Expected: FAIL (export command not registered)

**Step 3: Add CLI export command**

Find the `create_parser` function in `cli_unified.py` and add the export subparser. Add after other subparsers (around line 700):

```python
# Add to cli_unified.py in create_parser() function, after other subparsers

    # Export command
    export_parser = subparsers.add_parser(
        'export',
        help='Export primers for synthesis ordering',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
Export optimized primers in formats ready for synthesis ordering.

Generates:
  - FASTA file with primer sequences
  - Vendor-ready CSV (IDT, Twist, Sigma)
  - Wet-lab protocol with reaction conditions

Examples:
  neoswga export -d ./results/
  neoswga export -d ./results/ -o ./order/ --vendor idt --project MyProject
  neoswga export -d ./results/ --format fasta --output primers.fasta
"""
    )
    export_parser.add_argument('-d', '--dir', required=True,
                               help='Results directory containing step4_improved_df.csv')
    export_parser.add_argument('-o', '--output', default='./export',
                               help='Output directory (default: ./export)')
    export_parser.add_argument('--project', default='SWGA',
                               help='Project name for file naming (default: SWGA)')
    export_parser.add_argument('--vendor', default='idt',
                               choices=['idt', 'twist', 'sigma', 'generic'],
                               help='Vendor format for CSV (default: idt)')
    export_parser.add_argument('--format', choices=['all', 'fasta', 'csv', 'protocol'],
                               default='all',
                               help='Export format (default: all)')
    export_parser.add_argument('-j', '--json-file',
                               help='params.json for reaction conditions')
    export_parser.set_defaults(func=run_export)
```

Then add the handler function (add near other run_* functions):

```python
# Add to cli_unified.py

def run_export(args):
    """Export primers for synthesis ordering."""
    from neoswga.core.export import PrimerExporter
    from pathlib import Path

    try:
        # Load exporter from results
        exporter = PrimerExporter.from_results_dir(
            args.dir,
            params_file=getattr(args, 'json_file', None)
        )

        # Print summary
        exporter.print_summary()

        output_dir = Path(args.output)

        if args.format == 'all':
            outputs = exporter.export_all(
                str(output_dir),
                project_name=args.project,
                vendors=[args.vendor]
            )
            print("\nExported files:")
            for fmt, path in outputs.items():
                print(f"  {fmt}: {path}")

        elif args.format == 'fasta':
            output_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = output_dir / f"{args.project}_primers.fasta"
            exporter.export_fasta(str(fasta_path), prefix=args.project, include_metadata=True)
            print(f"Exported: {fasta_path}")

        elif args.format == 'csv':
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_path = output_dir / f"{args.project}_order_{args.vendor}.csv"
            exporter.export_vendor_csv(str(csv_path), vendor=args.vendor, project_name=args.project)
            print(f"Exported: {csv_path}")

        elif args.format == 'protocol':
            output_dir.mkdir(parents=True, exist_ok=True)
            protocol_path = output_dir / f"{args.project}_protocol.md"
            exporter.export_protocol(str(protocol_path))
            print(f"Exported: {protocol_path}")

        print("\nPrimers ready for ordering!")

    except FileNotFoundError as e:
        logger.error(str(e))
        print(f"Error: {e}")
        print("Make sure to run the full pipeline (count-kmers, filter, score, optimize) first.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Export failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
```

**Step 4: Run test to verify it passes**

```bash
pytest tests/test_export.py::TestExportCLI -v
```

Expected: PASS

**Step 5: Commit**

```bash
git add neoswga/cli_unified.py tests/test_export.py
git commit -m "feat(cli): add 'neoswga export' command for synthesis ordering"
```

---

## Task 6: Create Lab Workflow Documentation

**Files:**
- Create: `docs/FROM_RESULTS_TO_LAB.md`

**Step 1: Create the documentation file**

```markdown
# From Results to Lab: Complete Workflow Guide

This guide covers the complete workflow from running NeoSWGA to having oligos ready for your SWGA experiment.

## Overview

After running the NeoSWGA pipeline, you have optimized primer sequences. This guide covers:

1. Understanding your results
2. Exporting primers for ordering
3. Ordering from synthesis vendors
4. Receiving and preparing primers
5. Setting up your SWGA reaction

---

## 1. Understanding Your Results

### Pipeline Output Files

After running the four pipeline steps, your results directory contains:

| File | Description |
|------|-------------|
| `step4_improved_df.csv` | Final optimized primer set |
| `step3_df.csv` | Scored candidates (before optimization) |
| `step2_df.csv` | Filtered candidates |
| `params.json` | Parameters used for the run |

### Checking Quality

Before ordering, verify your primer set quality:

```bash
# Get quality assessment
neoswga interpret -d ./results/

# Generate detailed report
neoswga report -d ./results/
```

**Quality ratings:**
- **EXCELLENT/GOOD**: Ready for ordering
- **ACCEPTABLE**: Consider re-optimization
- **POOR/CRITICAL**: Re-run with different parameters

---

## 2. Exporting Primers for Ordering

### Quick Export (Recommended)

```bash
# Export all formats at once
neoswga export -d ./results/ -o ./order/ --project MyPathogen

# This creates:
#   ./order/MyPathogen_primers.fasta     - Sequences
#   ./order/MyPathogen_order_idt.csv     - IDT order form
#   ./order/MyPathogen_protocol.md       - Lab protocol
```

### Export Options

```bash
# Different vendors
neoswga export -d ./results/ --vendor twist   # Twist Bioscience
neoswga export -d ./results/ --vendor sigma   # Sigma-Aldrich
neoswga export -d ./results/ --vendor generic # Generic CSV

# Single format
neoswga export -d ./results/ --format fasta
neoswga export -d ./results/ --format csv
neoswga export -d ./results/ --format protocol
```

### Manual Export

If you prefer manual extraction:

```python
import pandas as pd

# Load results
df = pd.read_csv('./results/step4_improved_df.csv')

# Get unique primers from set 0 (best set)
primers = df[df['set_index'] == 0]['primer'].tolist()

# Save as FASTA
with open('primers.fasta', 'w') as f:
    for i, seq in enumerate(primers, 1):
        f.write(f">Primer_{i:02d}\n{seq}\n")
```

---

## 3. Ordering from Synthesis Vendors

### Recommended Vendors

| Vendor | Turnaround | Best For |
|--------|------------|----------|
| IDT | 1-2 days | Fast turnaround, reliable |
| Twist | 5-7 days | Bulk orders, pools |
| Sigma | 3-5 days | Budget option |
| Eurofins | 2-3 days | European labs |

### Order Specifications

For SWGA primers, use these specifications:

| Parameter | Recommendation |
|-----------|----------------|
| Scale | 25 nmol (standard) |
| Purification | Desalt (standard) |
| Format | Dry or 100 uM in TE |
| Quantity | 1 tube per primer |

**Cost estimate:** ~$5-8 per primer at standard scale

### IDT-Specific Instructions

1. Go to [IDT Oligo Order](https://www.idtdna.com/site/order)
2. Choose "Upload file" > select `*_order_idt.csv`
3. Review order and adjust scale if needed
4. Add to cart and checkout

### Twist-Specific Instructions

1. Go to [Twist Oligos](https://www.twistbioscience.com/products/genes)
2. Upload `*_order_twist.csv`
3. Select "Oligo Pools" for >10 primers

---

## 4. Receiving and Preparing Primers

### Upon Arrival

1. **Centrifuge** tubes briefly before opening
2. **Resuspend** (if dry):
   - Add TE buffer to 100 uM (see tube label for volume)
   - Vortex gently
   - Let sit 10 min at room temperature
3. **Store** at -20C

### Preparing Working Stocks

```
100 uM stock -> 10 uM working stock
   Add 10 uL stock to 90 uL TE buffer
```

### Preparing Primer Mix

For a 6-primer set at 10 uM each:
```
Combine equal volumes:
  10 uL of Primer 1 (10 uM)
  10 uL of Primer 2 (10 uM)
  10 uL of Primer 3 (10 uM)
  10 uL of Primer 4 (10 uM)
  10 uL of Primer 5 (10 uM)
  10 uL of Primer 6 (10 uM)
  ---
  60 uL total (each primer at 1.67 uM)
```

---

## 5. Setting Up Your SWGA Reaction

### Standard Protocol (Phi29, 30C)

**25 uL reaction:**

| Component | Volume | Final Conc |
|-----------|--------|------------|
| 10X Phi29 Buffer | 2.5 uL | 1X |
| dNTPs (10 mM each) | 1.0 uL | 0.4 mM |
| Primer mix | 3.0 uL | 200 nM each |
| Phi29 Polymerase | 0.5 uL | 10 U |
| Template DNA | 1.0 uL | 1-10 ng |
| Water | 17.0 uL | - |

**Incubation:** 30C for 16-24 hours, then 65C for 10 min

### Enhanced Protocol (EquiPhi29, 42C)

Add to the above:
- Betaine: 1.0 M final
- DMSO: 5% final

**Incubation:** 42C for 8-16 hours, then 65C for 10 min

### Quality Control

1. **Gel check:** Run 1 uL on 0.8% agarose
   - Expected: High MW smear (>10 kb)
   - Bad: No product or ladder pattern

2. **Qubit/NanoDrop:**
   - Expected yield: 1-10 ug from 1 ng input
   - A260/280 > 1.8

3. **qPCR verification:**
   - Test target vs background loci
   - Expected enrichment: 10-1000x

---

## Troubleshooting

| Problem | Cause | Solution |
|---------|-------|----------|
| No amplification | Poor primer quality | Reorder, check Tm |
| Low yield | Insufficient primers | Increase primer conc |
| Ladder pattern | Primer dimers | Reduce primer conc |
| Poor enrichment | Non-specific priming | Use EquiPhi29 + additives |

---

## Quick Reference

```bash
# Complete workflow
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
neoswga interpret -d ./data/
neoswga export -d ./data/ -o ./order/
```

---

*See also: [SWGA Science](SWGA_SCIENCE.md) for thermodynamics and additive details*
```

**Step 2: Commit documentation**

```bash
git add docs/FROM_RESULTS_TO_LAB.md
git commit -m "docs: add complete lab workflow guide (results to oligos)"
```

---

## Task 7: Update Documentation Index

**Files:**
- Modify: `docs/README.md`
- Modify: `README.md`

**Step 1: Update docs/README.md**

Add to the User Guides table:

```markdown
| [From Results to Lab](FROM_RESULTS_TO_LAB.md) | Export, ordering, and lab setup |
```

**Step 2: Update main README.md**

Add to Documentation section:

```markdown
- **[From Results to Lab](docs/FROM_RESULTS_TO_LAB.md)**: Export primers and lab workflow
```

**Step 3: Commit**

```bash
git add docs/README.md README.md
git commit -m "docs: add lab workflow guide to documentation index"
```

---

## Task 8: Integration Test

**Files:**
- Test: `tests/test_export.py`

**Step 1: Write integration test**

```python
# Add to tests/test_export.py

class TestExportIntegration:
    """Integration tests for complete export workflow."""

    def test_full_export_workflow(self, tmp_path):
        """Test complete export from results to all formats."""
        from neoswga.core.export import PrimerExporter
        import json

        # Create realistic mock results
        results_dir = tmp_path / "results"
        results_dir.mkdir()

        # Step 4 results
        step4_content = """primer,set_index,score,coverage,bg_coverage,selectivity,mean_gap,optimizer
ATCGATCGATCG,0,0.847,0.78,0.12,6.5,15000,hybrid
GCTAGCTAGCTA,0,0.847,0.78,0.12,6.5,15000,hybrid
TTAATTAATTAA,0,0.847,0.78,0.12,6.5,15000,hybrid
CCGGCCGGCCGG,0,0.847,0.78,0.12,6.5,15000,hybrid
AACCAACCAACC,0,0.847,0.78,0.12,6.5,15000,hybrid
TTGGTTGGTTGG,0,0.847,0.78,0.12,6.5,15000,hybrid
"""
        (results_dir / "step4_improved_df.csv").write_text(step4_content)

        # Params file
        params = {
            "polymerase": "equiphi29",
            "reaction_temp": 42.0,
            "betaine_m": 1.0,
            "dmso_percent": 5.0,
            "mg_conc": 2.5,
        }
        (results_dir / "params.json").write_text(json.dumps(params))

        # Run export
        exporter = PrimerExporter.from_results_dir(str(results_dir))

        export_dir = tmp_path / "export"
        outputs = exporter.export_all(str(export_dir), project_name="IntegrationTest")

        # Verify all outputs
        assert len(outputs) == 3  # fasta, csv, protocol

        # Check FASTA content
        fasta_content = (export_dir / "IntegrationTest_primers.fasta").read_text()
        assert ">IntegrationTest_001" in fasta_content
        assert "ATCGATCGATCG" in fasta_content
        assert "Tm=" in fasta_content

        # Check CSV content
        csv_content = (export_dir / "IntegrationTest_order_idt.csv").read_text()
        assert "Name" in csv_content
        assert "IntegrationTest_001" in csv_content
        assert "25nm" in csv_content.lower()

        # Check protocol content
        protocol_content = (export_dir / "IntegrationTest_protocol.md").read_text()
        assert "equiphi29" in protocol_content.lower() or "EquiPhi29" in protocol_content
        assert "42" in protocol_content
        assert "Betaine" in protocol_content
        assert "DMSO" in protocol_content

    def test_summary_accuracy(self):
        """Verify summary calculations are accurate."""
        from neoswga.core.export import PrimerExporter

        # Known primers with calculable properties
        primers = [
            "AAAAAAAAAA",  # 10-mer, 0% GC, Tm = 20
            "GGGGGGGGGG",  # 10-mer, 100% GC, Tm = 40
            "ATCGATCGAT",  # 10-mer, 40% GC, Tm = 28
        ]

        exporter = PrimerExporter(primers)
        summary = exporter.get_summary()

        assert summary["num_primers"] == 3
        assert summary["mean_length"] == 10.0
        assert abs(summary["mean_gc"] - 0.467) < 0.01  # (0 + 1 + 0.4) / 3
        assert summary["min_tm"] == 20.0
        assert summary["max_tm"] == 40.0
        assert summary["estimated_cost"] == 15.0  # 3 * $5
```

**Step 2: Run all tests**

```bash
pytest tests/test_export.py -v
```

Expected: All PASS

**Step 3: Commit**

```bash
git add tests/test_export.py
git commit -m "test: add integration tests for export workflow"
```

---

## Summary

This plan delivers:

1. **`neoswga/core/export.py`** - Complete export module with:
   - FASTA export with metadata
   - Vendor-specific CSV (IDT, Twist, Sigma, generic)
   - Wet-lab protocol generation
   - `PrimerExporter` class for unified export

2. **`neoswga export` CLI command** - Easy primer export:
   ```bash
   neoswga export -d ./results/ -o ./order/ --project MyProject
   ```

3. **`docs/FROM_RESULTS_TO_LAB.md`** - Complete guide covering:
   - Understanding results
   - Export commands
   - Vendor ordering instructions
   - Primer preparation
   - SWGA reaction setup

4. **Updated documentation** - Links to new guide in README files

---

**Plan complete and saved to `docs/plans/2026-02-04-export-and-lab-integration.md`. Two execution options:**

**1. Subagent-Driven (this session)** - I dispatch fresh subagent per task, review between tasks, fast iteration

**2. Parallel Session (separate)** - Open new session with executing-plans, batch execution with checkpoints

**Which approach?**
