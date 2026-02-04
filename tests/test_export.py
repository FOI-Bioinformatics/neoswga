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
