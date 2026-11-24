"""Reporting / export CLI handlers: report and export (with their private
primer-position and genome-length loaders). Extracted from cli_unified.
"""

import json
import logging
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


def _run_targets(results_dir, params_file=None):
    """(params, target prefixes, target FASTAs) for the run being exported.

    The params file is `-j` when given, else `DIR/params.json`; relative paths
    in it resolve against its own directory, as `plan-pool` resolves them.
    """
    from pathlib import Path

    path = Path(params_file) if params_file else Path(results_dir) / "params.json"
    if not path.is_file():
        raise ValueError(
            f"BED/BedGraph export needs the run's params.json to find its target: "
            f"{path} does not exist. Pass it with -j."
        )
    params = json.loads(path.read_text())
    base = path.parent
    prefixes = [str(base / p) for p in params.get("fg_prefixes") or []]
    genomes = [str(base / g) for g in params.get("fg_genomes") or []]
    return params, prefixes, genomes


def _load_primer_positions(results_dir, primers, params_file=None):
    """Each primer's binding sites on the run's one target, by strand.

    The target comes from the params file (see `_run_targets`). The cache scans
    the target FASTA for a primer the position index does not hold, and
    refuses without one. It replaces a glob over every `*_positions.h5` in the
    directory, which wrote a host's sites under the target's name and returned
    no sites, silently, for a primer it could not find (2026-11-24).

    Returns:
        Dict mapping primer to list of (position, strand) tuples.
    """
    from neoswga.cli._common import design_position_cache

    _, prefixes, genomes = _run_targets(results_dir, params_file)
    if len(prefixes) != 1:
        raise ValueError(
            f"BED/BedGraph export writes one target and the run lists "
            f"{len(prefixes)} (fg_prefixes); each would need its own chromosome name."
        )
    prefix = prefixes[0]
    cache = design_position_cache(prefixes, primers, genomes[:1] if genomes else [])
    return {
        primer: [(int(x), "forward") for x in cache.get_positions(prefix, primer, "forward")]
        + [(int(x), "reverse") for x in cache.get_positions(prefix, primer, "reverse")]
        for primer in primers
    }


def _blocking_validator_findings(results_dir):
    """Details of any recorded finding that makes this pool unfit to order.

    Delegates to `core.design_result`, which owns the code set. It lived here
    as a bare literal and again in `results_interpreter`, so a new blocking code
    had to be added twice or the two commands would disagree about one pool.
    Kept as a name because several tests and `run_export` call it.
    """
    from neoswga.core.design_result import blocking_validator_findings

    return blocking_validator_findings(results_dir)


def _get_genome_length(results_dir, params_file=None):
    """The target's length in bp: recorded in params.json, else measured.

    Returned 0 when params.json carried no `fg_seq_lengths`, and BedGraph then
    wrote no windows; an unknown length now refuses (2026-11-24).
    """
    params, _, genomes = _run_targets(results_dir, params_file)
    lengths = params.get("fg_seq_lengths") or []
    if lengths:
        return int(sum(lengths))
    if genomes:
        from neoswga.core.utility import get_all_seq_lengths

        return int(sum(get_all_seq_lengths(fname_genomes=genomes, cpus=1)))
    raise ValueError(
        "BedGraph needs the target length: params.json records no fg_seq_lengths "
        "and names no fg_genomes to measure it from."
    )


def run_export(args):
    """Export primers for synthesis ordering."""
    from pathlib import Path

    from neoswga.core.export import PrimerExporter, PrimerModifications

    if getattr(args, "no_modifications", False) and args.pto_bonds is not None:
        # The bonds were dropped without a word, so an order asked for with
        # them went out bare (2026-11-24).
        logger.error(
            "--no-modifications and --pto-bonds contradict each other: "
            "pass one. Bare sequences: --no-modifications; bonds: --pto-bonds N."
        )
        sys.exit(1)

    try:
        # Determine modification profile
        if getattr(args, "no_modifications", False):
            mods = PrimerModifications.from_profile("none")
        else:
            mods = PrimerModifications.from_profile(args.modifications)
            # Override PTO bonds if specified
            if args.pto_bonds is not None:
                mods.pto_bonds = args.pto_bonds

        # Refuse a directory whose most recent run did not finish with a
        # qualifying panel. The result files here would be from an earlier
        # run: real, stale, and indistinguishable from current without this.
        from neoswga.core.export import export_is_blocked

        # The set about to be WRITTEN, not the directory alone: the findings
        # in it describe set 0, so a gate given only the directory judged a
        # different panel from the one being exported.
        blocked_run = export_is_blocked(args.dir, getattr(args, "set_index", 0) or 0)
        if blocked_run and getattr(args, "allow_unqualified", False):
            # The same shape as --allow-dimer-relaxation: the constraint is
            # traded deliberately, by name, and the run says so loudly.
            logger.warning("Exporting anyway because --allow-unqualified was given.")
            logger.warning("  %s", blocked_run)
            blocked_run = None
        if blocked_run:
            logger.error(blocked_run)
            sys.exit(1)

        # Load exporter from results
        exporter = PrimerExporter.from_results_dir(
            args.dir,
            params_file=getattr(args, "json_file", None),
            set_index=getattr(args, "set_index", 0),
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
            params_file = getattr(args, "json_file", None)
            positions = _load_primer_positions(args.dir, exporter.primers, params_file)

            if args.format == "bed":
                bed_path = output_dir / f"{args.project}_binding_sites.bed"
                exporter.export_bed(str(bed_path), positions, genome_name=args.genome_name)
                print(f"Exported: {bed_path}")
            else:
                genome_length = _get_genome_length(args.dir, params_file)
                bg_path = output_dir / f"{args.project}_coverage.bedgraph"
                exporter.export_bedgraph(
                    str(bg_path),
                    positions,
                    genome_name=args.genome_name,
                    genome_length=genome_length,
                    window_size=args.window_size,
                )
                print(f"Exported: {bg_path}")

        # The pool was checked before any of this was written, by
        # `export_is_blocked`. Reaching here means it passed, or that the user
        # traded the constraint explicitly.
        if getattr(args, "allow_unqualified", False) and _blocking_validator_findings(args.dir):
            print("\nExported, but NOT ready for ordering. The optimizer recorded:")
            for detail in _blocking_validator_findings(args.dir):
                print(f"  - {detail}")
            print(f"\nRun `neoswga interpret -d {args.dir}` for the full assessment.")
        else:
            print("\nPrimers ready for ordering!")

    except FileNotFoundError as e:
        logger.error(str(e))
        logger.error(
            "Make sure to run the full pipeline (count-kmers, filter, prepare-candidates, optimize) first."
        )
        sys.exit(1)
    except Exception as e:
        logger.error(f"Export failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def add_parsers(subparsers):
    """Register this group's subcommands on the shared subparsers object.

    Called by neoswga.cli_unified.create_parser(). Extracted from the former
    monolithic create_parser() so each command group owns its argparse setup
    next to its handlers.
    """
    import argparse  # noqa: F401  (used by some command blocks)

    from neoswga.cli._common import add_common_options  # noqa: F401

    report_parser = subparsers.add_parser(
        "report", help="Generate comprehensive HTML report for primer set"
    )
    report_parser.add_argument(
        "-d", "--dir", required=True, help="Results directory containing pipeline output"
    )
    report_parser.add_argument(
        "-o", "--output", help="Output file path (default: <dir>/report.html)"
    )
    report_parser.add_argument(
        "--format", choices=["html", "json"], default="html", help="Output format (default: html)"
    )
    report_parser.add_argument(
        "--level",
        choices=["summary", "full"],
        default="summary",
        help="Report level: summary (1-page) or full (technical report)",
    )
    report_parser.add_argument(
        "--check", action="store_true", help="Validate only, do not generate report"
    )
    report_parser.add_argument(
        "--interactive",
        action="store_true",
        help="Include interactive Plotly charts (requires plotly)",
    )
    report_parser.add_argument(
        "-q", "--quiet", action="store_true", help="Suppress progress messages"
    )

    # =========================================================================
    # EXPORT: Primer export for synthesis ordering
    # =========================================================================

    export_parser = subparsers.add_parser(
        "export",
        help="Export primers for synthesis ordering",
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
""",
    )
    export_parser.add_argument(
        "-d", "--dir", required=True, help="Results directory containing step4_improved_df.csv"
    )
    export_parser.add_argument(
        "-o", "--output", default="./export", help="Output directory (default: ./export)"
    )
    export_parser.add_argument(
        "--allow-unqualified",
        action="store_true",
        help="Export even when the optimizer recorded a defect in the pool, "
        "such as a delivered pair above max_dimer_bp. The findings are printed "
        "and the files are written anyway. Without this the export refuses and "
        "writes nothing.",
    )
    export_parser.add_argument(
        "--set",
        dest="set_index",
        type=int,
        default=0,
        help="Which of the alternative primer sets to use (default: 0, the set the run summary describes). Alternatives are separate answers to the same design, not additions to it; they are never pooled.",
    )
    export_parser.add_argument(
        "--project", default="SWGA", help="Project name for file naming (default: SWGA)"
    )
    export_parser.add_argument(
        "--vendor",
        default="idt",
        choices=["idt", "twist", "sigma", "generic"],
        help="Vendor format for CSV (default: idt)",
    )
    export_parser.add_argument(
        "--format",
        choices=["all", "fasta", "csv", "protocol", "bed", "bedgraph"],
        default="all",
        help="Export format (default: all)",
    )
    export_parser.add_argument("-j", "--json-file", help="params.json for reaction conditions")
    export_parser.add_argument(
        "--modifications",
        choices=["none", "standard", "low-input"],
        default="standard",
        help="Modification profile: none (bare primers), "
        "standard (3' PTO for exonuclease protection), "
        "low-input (3' PTO + 5' C18 spacer for tdMDA) "
        "(default: standard)",
    )
    export_parser.add_argument(
        "--pto-bonds",
        type=int,
        default=None,
        help="Override number of 3' phosphorothioate bonds (default: 2)",
    )
    export_parser.add_argument(
        "--no-modifications",
        action="store_true",
        help="Export bare sequences without modifications " "(equivalent to --modifications none)",
    )
    export_parser.add_argument(
        "--genome-name",
        default="genome",
        help="Chromosome/contig name for BED/BedGraph output (default: genome)",
    )
    export_parser.add_argument(
        "--window-size",
        type=int,
        default=1000,
        help="Window size in bp for BedGraph binning (default: 1000)",
    )
    export_parser.set_defaults(func=run_export)

    # =========================================================================
    # SETUP: Interactive workflow selector
    # =========================================================================
