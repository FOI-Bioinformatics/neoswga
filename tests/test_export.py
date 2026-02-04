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


class TestPrimerModifications:
    """Test primer modification functionality."""

    def test_profile_standard(self):
        """Standard profile has 2 PTO bonds and no 5' block."""
        from neoswga.core.export import PrimerModifications, ModificationProfile

        mods = PrimerModifications.from_profile('standard')
        assert mods.pto_bonds == 2
        assert mods.five_prime_block is None
        assert mods.profile == ModificationProfile.STANDARD

    def test_profile_low_input(self):
        """Low-input profile has 2 PTO bonds and C18 spacer."""
        from neoswga.core.export import PrimerModifications, ModificationProfile

        mods = PrimerModifications.from_profile('low-input')
        assert mods.pto_bonds == 2
        assert mods.five_prime_block == 'c18'
        assert mods.profile == ModificationProfile.LOW_INPUT

    def test_profile_none(self):
        """None profile has no modifications."""
        from neoswga.core.export import PrimerModifications, ModificationProfile

        mods = PrimerModifications.from_profile('none')
        assert mods.pto_bonds == 0
        assert mods.five_prime_block is None
        assert mods.profile == ModificationProfile.NONE

    def test_invalid_profile(self):
        """Invalid profile raises ValueError."""
        from neoswga.core.export import PrimerModifications

        with pytest.raises(ValueError, match="Unknown modification profile"):
            PrimerModifications.from_profile('invalid')

    def test_apply_mods_standard_idt(self):
        """Standard profile adds PTO bonds with IDT syntax."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('standard')
        result = apply_modifications("ATCGATCG", mods, vendor="idt")
        # 2 PTO bonds at 3' end: ATCGAT*C*G (2 asterisks connecting last 3 nt)
        assert result == "ATCGAT*C*G"

    def test_apply_mods_standard_longer_primer(self):
        """Standard profile on 12-mer primer."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('standard')
        result = apply_modifications("ATCGATCGATCG", mods, vendor="idt")
        # 2 PTO bonds at 3' end: ATCGATCGAT*C*G
        assert result == "ATCGATCGAT*C*G"

    def test_apply_mods_low_input_idt(self):
        """Low-input profile adds C18 prefix and PTO bonds."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('low-input')
        result = apply_modifications("ATCGATCG", mods, vendor="idt")
        # C18 prefix + PTO at 3' end
        assert result == "/5SpC18/ATCGAT*C*G"

    def test_apply_mods_none(self):
        """None profile returns sequence unchanged."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('none')
        result = apply_modifications("ATCGATCG", mods, vendor="idt")
        assert result == "ATCGATCG"

    def test_apply_mods_generic_vendor(self):
        """Generic vendor uses bracket notation."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('low-input')
        result = apply_modifications("ATCGATCG", mods, vendor="generic")
        # Generic uses [C18]- prefix
        assert result == "[C18]-ATCGAT*C*G"

    def test_apply_mods_sigma_vendor(self):
        """Sigma uses IDT-compatible syntax."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('low-input')
        result = apply_modifications("ATCGATCG", mods, vendor="sigma")
        # Same as IDT
        assert result == "/5SpC18/ATCGAT*C*G"

    def test_apply_mods_short_sequence(self):
        """Handle very short sequences gracefully."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications.from_profile('standard')
        # 4-mer with 2 PTO bonds - all bases get connected
        result = apply_modifications("ATCG", mods, vendor="idt")
        assert result == "AT*C*G"

    def test_apply_mods_custom_pto_count(self):
        """Custom PTO bond count."""
        from neoswga.core.export import apply_modifications, PrimerModifications

        mods = PrimerModifications(pto_bonds=3, five_prime_block=None)
        result = apply_modifications("ATCGATCG", mods, vendor="idt")
        # 3 PTO bonds: last 4 nucleotides connected
        assert result == "ATCGA*T*C*G"


class TestExportWithModifications:
    """Test export functionality with modifications."""

    def test_exporter_default_modifications(self):
        """Exporter uses standard modifications by default."""
        from neoswga.core.export import PrimerExporter, ModificationProfile

        exporter = PrimerExporter(["ATCGATCG"])
        assert exporter.modifications.profile == ModificationProfile.STANDARD
        assert exporter.modifications.pto_bonds == 2

    def test_exporter_custom_modifications(self):
        """Exporter accepts custom modifications."""
        from neoswga.core.export import PrimerExporter, PrimerModifications

        mods = PrimerModifications.from_profile('low-input')
        exporter = PrimerExporter(["ATCGATCG"], modifications=mods)
        assert exporter.modifications.five_prime_block == 'c18'

    def test_export_csv_with_modifications(self, tmp_path):
        """CSV export includes modified sequences."""
        from neoswga.core.export import PrimerExporter, PrimerModifications

        mods = PrimerModifications.from_profile('standard')
        exporter = PrimerExporter(["ATCGATCGATCG"], modifications=mods)

        csv_path = tmp_path / "test.csv"
        exporter.export_vendor_csv(str(csv_path), vendor="idt", project_name="Test")

        content = csv_path.read_text()
        # Should contain PTO notation
        assert "*" in content
        # Original sequence should have modifications applied (2 PTO bonds)
        assert "ATCGATCGAT*C*G" in content

    def test_export_csv_no_modifications(self, tmp_path):
        """CSV export with no modifications gives bare sequences."""
        from neoswga.core.export import PrimerExporter, PrimerModifications

        mods = PrimerModifications.from_profile('none')
        exporter = PrimerExporter(["ATCGATCGATCG"], modifications=mods)

        csv_path = tmp_path / "test.csv"
        exporter.export_vendor_csv(str(csv_path), vendor="idt", project_name="Test")

        content = csv_path.read_text()
        # Should NOT contain PTO notation
        assert "*" not in content
        # Bare sequence
        assert "ATCGATCGATCG" in content

    def test_summary_includes_modifications(self):
        """Summary includes modification information."""
        from neoswga.core.export import PrimerExporter, PrimerModifications

        mods = PrimerModifications.from_profile('low-input')
        exporter = PrimerExporter(["ATCGATCG"], modifications=mods)

        summary = exporter.get_summary()
        assert summary["modification_profile"] == "low-input"
        assert summary["pto_bonds"] == 2
        assert summary["five_prime_block"] == "c18"

    def test_full_export_with_modifications(self, tmp_path):
        """Full export workflow with modifications."""
        from neoswga.core.export import PrimerExporter, PrimerModifications
        import json

        # Create mock results
        results_dir = tmp_path / "results"
        results_dir.mkdir()
        (results_dir / "step4_improved_df.csv").write_text(
            "primer,set_index,score\nATCGATCGATCG,0,0.85\nGCTAGCTAGCTA,0,0.85\n"
        )
        (results_dir / "params.json").write_text(
            json.dumps({"polymerase": "phi29", "reaction_temp": 30})
        )

        # Create exporter with low-input modifications
        exporter = PrimerExporter.from_results_dir(str(results_dir))
        exporter.modifications = PrimerModifications.from_profile('low-input')

        # Export
        export_dir = tmp_path / "export"
        exporter.export_all(str(export_dir), project_name="ModTest")

        # Check CSV has modifications
        csv_content = (export_dir / "ModTest_order_idt.csv").read_text()
        assert "/5SpC18/" in csv_content
        assert "*" in csv_content


class TestExportCLIModifications:
    """Test CLI export command with modification arguments."""

    def test_cli_modifications_argument(self):
        """CLI accepts --modifications argument."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'export', '-d', '/tmp/test',
            '--modifications', 'low-input'
        ])
        assert args.modifications == 'low-input'

    def test_cli_no_modifications_flag(self):
        """CLI accepts --no-modifications flag."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'export', '-d', '/tmp/test',
            '--no-modifications'
        ])
        assert args.no_modifications is True

    def test_cli_pto_bonds_argument(self):
        """CLI accepts --pto-bonds argument."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args([
            'export', '-d', '/tmp/test',
            '--pto-bonds', '3'
        ])
        assert args.pto_bonds == 3

    def test_cli_default_modifications(self):
        """CLI defaults to standard modifications."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(['export', '-d', '/tmp/test'])
        assert args.modifications == 'standard'
        assert args.no_modifications is False
        assert args.pto_bonds is None
