"""Reporting / export CLI handlers: report and export (with their private
primer-position and genome-length loaders). Extracted from cli_unified.
"""

import json
import logging
import os
import sys

logger = logging.getLogger(__name__)


def run_report(args):
    """Generate comprehensive report for primer set"""
    from pathlib import Path

    # Import validation module
    from neoswga.core.report.validation import (
        ValidationLevel,
        validate_metrics,
        validate_results_directory,
    )

    results_dir = Path(args.dir)
    quiet = getattr(args, "quiet", False)
    check_only = getattr(args, "check", False)
    interactive = getattr(args, "interactive", False)

    # Check if interactive mode is requested but Plotly is not available
    if interactive:
        from neoswga.core.report.visualizations import is_plotly_available

        if not is_plotly_available():
            print("Warning: Interactive charts requested but Plotly is not installed.")
            print("Install with: pip install plotly  (or: pip install neoswga[interactive])")
            print("Continuing without interactive charts...")
            interactive = False

    def progress(msg):
        """Print progress message unless quiet mode."""
        if not quiet:
            print(msg)

    # Step 1: Validate results directory
    progress("Validating results directory...")
    dir_result = validate_results_directory(args.dir)

    # Display validation issues
    for issue in dir_result.issues:
        if issue.level == ValidationLevel.ERROR:
            print(f"  ERROR: {issue.message}")
        elif issue.level == ValidationLevel.WARNING:
            print(f"  WARNING: {issue.message}")
        elif issue.level == ValidationLevel.INFO and not quiet:
            print(f"  INFO: {issue.message}")

    if not dir_result.is_valid:
        print(f"\nValidation failed: {len(dir_result.errors)} error(s) found.")
        print("Cannot generate report. Please fix the errors above.")
        sys.exit(1)

    # If check-only mode, stop here
    if check_only:
        if dir_result.warnings:
            print(f"\nValidation passed with {len(dir_result.warnings)} warning(s).")
        else:
            print("\nValidation passed. Ready to generate report.")
        return

    # Determine output path
    level = getattr(args, "level", "summary")
    if args.output:
        output_path = args.output
    else:
        if level == "full":
            output_path = str(results_dir / "technical_report.html")
        else:
            output_path = str(results_dir / "report.html")

    try:
        if args.format == "json":
            # JSON export
            import json

            from neoswga.core.report.metrics import collect_pipeline_metrics
            from neoswga.core.report.quality import calculate_quality_grade

            progress("Collecting pipeline metrics...")
            metrics = collect_pipeline_metrics(args.dir)

            # Validate metrics
            metrics_result = validate_metrics(metrics)
            for issue in metrics_result.warnings:
                print(f"  WARNING: {issue.message}")

            progress("Calculating quality grade...")
            quality = calculate_quality_grade(metrics)

            progress("Generating JSON output...")

            # Build JSON output
            report_data = {
                "grade": quality.grade.value,
                "composite_score": quality.composite_score,
                "recommendation": quality.recommendation,
                "recommendation_details": quality.recommendation_details,
                "considerations": quality.considerations,
                "primer_count": metrics.primer_count,
                "components": [
                    {
                        "name": c.name,
                        "weight": c.weight,
                        "raw_value": c.raw_value,
                        "normalized_score": c.normalized_score,
                        "rating": c.rating,
                    }
                    for c in quality.components
                ],
            }

            # Change extension if needed
            if not output_path.endswith(".json"):
                output_path = output_path.rsplit(".", 1)[0] + ".json"

            with open(output_path, "w") as f:
                json.dump(report_data, f, indent=2)

            print(
                f"\nQuality Grade: {quality.grade.value} " f"({quality.composite_score:.2f}/1.00)"
            )
            print(f"Recommendation: {quality.recommendation}")
            print(f"\nJSON report saved to: {output_path}")

        elif level == "full":
            # Full technical report
            try:
                from neoswga.core.report import generate_technical_report
            except ImportError as e:
                logger.error(f"Report module not available: {e}")
                sys.exit(1)

            progress("Collecting pipeline metrics...")
            progress("Calculating quality grade...")
            if interactive:
                progress("Generating technical report with interactive charts...")
            else:
                progress("Generating technical report...")

            data = generate_technical_report(args.dir, output_path, interactive=interactive)
            print(
                f"\nQuality Grade: {data.quality.grade.value} "
                f"({data.quality.composite_score:.2f}/1.00)"
            )
            print(f"Recommendation: {data.quality.recommendation}")
            print(f"Primers analyzed: {data.metrics.primer_count}")
            if interactive:
                print(f"\nTechnical report with interactive charts saved to: {output_path}")
            else:
                print(f"\nTechnical report saved to: {output_path}")

        else:
            # Executive summary (default)
            try:
                from neoswga.core.report import generate_executive_summary
            except ImportError as e:
                logger.error(f"Report module not available: {e}")
                sys.exit(1)

            progress("Collecting pipeline metrics...")
            progress("Calculating quality grade...")
            if interactive:
                progress("Generating executive summary with interactive charts...")
            else:
                progress("Generating executive summary...")

            summary = generate_executive_summary(args.dir, output_path, interactive=interactive)
            print(
                f"\nQuality Grade: {summary.quality.grade.value} "
                f"({summary.quality.composite_score:.2f}/1.00)"
            )
            print(f"Recommendation: {summary.quality.recommendation}")
            if interactive:
                print(f"\nReport with interactive charts saved to: {output_path}")
            else:
                print(f"\nReport saved to: {output_path}")

    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to generate report: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def _load_primer_positions(results_dir, primers):
    """Load primer binding positions from HDF5 files in a results directory.

    Args:
        results_dir: Path to results directory.
        primers: List of primer sequences.

    Returns:
        Dict mapping primer to list of (position, strand) tuples.
    """
    import glob as glob_module

    import h5py

    from neoswga.core.thermodynamics import reverse_complement

    positions = {}
    h5_files = glob_module.glob(os.path.join(results_dir, "*_positions.h5"))

    if not h5_files:
        params_path = os.path.join(results_dir, "params.json")
        if os.path.exists(params_path):
            with open(params_path) as f:
                params = json.load(f)
            fg_prefixes = params.get("fg_prefixes", [])
            for prefix in fg_prefixes:
                h5_files.extend(glob_module.glob(f"{prefix}_*_positions.h5"))

    if not h5_files:
        logger.warning(
            "No position HDF5 files found. "
            "BED/BedGraph export requires position data from the filter step."
        )
        return {p: [] for p in primers}

    rc_map = {p: reverse_complement(p) for p in primers}

    for primer in primers:
        sites = []
        k = len(primer)
        for h5_path in h5_files:
            if f"_{k}mer_" not in h5_path:
                continue
            try:
                with h5py.File(h5_path, "r") as db:
                    if primer in db:
                        for pos in db[primer][:]:
                            sites.append((int(pos), "forward"))
                    rc = rc_map[primer]
                    if rc in db:
                        for pos in db[rc][:]:
                            sites.append((int(pos), "reverse"))
            except Exception as e:
                logger.warning(f"Error reading {h5_path}: {e}")
        positions[primer] = sites

    return positions


def _get_genome_length(results_dir):
    """Get genome length from params.json in results directory.

    Args:
        results_dir: Path to results directory.

    Returns:
        Total genome length in bp, or 0 if not determinable.
    """
    params_path = os.path.join(results_dir, "params.json")
    if os.path.exists(params_path):
        with open(params_path) as f:
            params = json.load(f)
        lengths = params.get("fg_seq_lengths", [])
        if lengths:
            return sum(lengths)
    logger.warning("Could not determine genome length. Using 0.")
    return 0


def run_export(args):
    """Export primers for synthesis ordering."""
    from pathlib import Path

    from neoswga.core.export import PrimerExporter, PrimerModifications

    try:
        # Determine modification profile
        if getattr(args, "no_modifications", False):
            mods = PrimerModifications.from_profile("none")
        else:
            mods = PrimerModifications.from_profile(args.modifications)
            # Override PTO bonds if specified
            if args.pto_bonds is not None:
                mods.pto_bonds = args.pto_bonds

        # Load exporter from results
        exporter = PrimerExporter.from_results_dir(
            args.dir, params_file=getattr(args, "json_file", None)
        )
        exporter.modifications = mods

        # Print summary
        exporter.print_summary()

        output_dir = Path(args.output)

        if args.format == "all":
            outputs = exporter.export_all(
                str(output_dir), project_name=args.project, vendors=[args.vendor]
            )
            print("\nExported files:")
            for fmt, path in outputs.items():
                print(f"  {fmt}: {path}")

        elif args.format == "fasta":
            output_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = output_dir / f"{args.project}_primers.fasta"
            exporter.export_fasta(str(fasta_path), prefix=args.project, include_metadata=True)
            print(f"Exported: {fasta_path}")

        elif args.format == "csv":
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_path = output_dir / f"{args.project}_order_{args.vendor}.csv"
            exporter.export_vendor_csv(str(csv_path), vendor=args.vendor, project_name=args.project)
            print(f"Exported: {csv_path}")

        elif args.format == "protocol":
            output_dir.mkdir(parents=True, exist_ok=True)
            protocol_path = output_dir / f"{args.project}_protocol.md"
            exporter.export_protocol(str(protocol_path))
            print(f"Exported: {protocol_path}")

        elif args.format in ("bed", "bedgraph"):
            output_dir.mkdir(parents=True, exist_ok=True)

            # Load positions from results directory
            positions = _load_primer_positions(args.dir, exporter.primers)

            if args.format == "bed":
                bed_path = output_dir / f"{args.project}_binding_sites.bed"
                exporter.export_bed(str(bed_path), positions, genome_name=args.genome_name)
                print(f"Exported: {bed_path}")
            else:
                genome_length = _get_genome_length(args.dir)
                bg_path = output_dir / f"{args.project}_coverage.bedgraph"
                exporter.export_bedgraph(
                    str(bg_path),
                    positions,
                    genome_name=args.genome_name,
                    genome_length=genome_length,
                    window_size=args.window_size,
                )
                print(f"Exported: {bg_path}")

        print("\nPrimers ready for ordering!")

    except FileNotFoundError as e:
        logger.error(str(e))
        logger.error(
            "Make sure to run the full pipeline (count-kmers, filter, score, optimize) first."
        )
        sys.exit(1)
    except Exception as e:
        logger.error(f"Export failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
