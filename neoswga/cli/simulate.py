"""Simulation CLI handlers: ``multi-genome`` pan-design and ``simulate``
(agent-based replication simulation with optional plots/report).

Extracted from cli_unified.py. Shared helpers come from neoswga.cli._common.
"""

import json
import logging
import os
import sys

from neoswga.cli._common import _apply_seed
from neoswga.core.registry import polymerase_names as _polymerase_names

logger = logging.getLogger(__name__)


def run_multi_genome(args):
    """Pan-genome primer design for multiple targets"""
    from pathlib import Path

    from neoswga.core.multi_genome_filter import GenomeSet
    from neoswga.core.multi_genome_pipeline import MultiGenomePipeline

    _apply_seed(args)

    logger.info(f"Designing pan-genome primers for {len(args.genomes)} genomes")

    try:
        # Create genome set from input genomes
        genome_set = GenomeSet()

        # Add all genomes as targets
        for i, genome_path in enumerate(args.genomes):
            genome_name = Path(genome_path).stem
            genome_set.add_genome(name=genome_name, fasta_path=genome_path, role="target")
            logger.info(f"  Added target: {genome_name}")

        # Add background genomes if provided
        if hasattr(args, "background") and args.background:
            for bg_path in args.background:
                bg_name = Path(bg_path).stem
                genome_set.add_genome(name=bg_name, fasta_path=bg_path, role="background")
                logger.info(f"  Added background: {bg_name}")

        # Add blacklist genomes if provided
        if hasattr(args, "blacklist") and args.blacklist:
            for bl_path in args.blacklist:
                bl_name = Path(bl_path).stem
                genome_set.add_genome(
                    name=bl_name, fasta_path=bl_path, role="blacklist", penalty_weight=5.0
                )
                logger.info(f"  Added blacklist: {bl_name}")

        # Get optional parameters
        primer_count = getattr(args, "num_primers", 12)
        min_k = getattr(args, "min_k", None)
        max_k = getattr(args, "max_k", None)
        kmer_range = (min_k, max_k) if min_k and max_k else None
        polymerase = getattr(args, "polymerase", None)

        # Create and run pipeline
        pipeline = MultiGenomePipeline(
            genome_set=genome_set,
            output_dir=args.output,
            kmer_range=kmer_range,
            preferred_polymerase=polymerase,
            primer_count=primer_count,
            validate_with_simulation=getattr(args, "validate_simulation", False),
        )

        result = pipeline.run(verbose=not getattr(args, "quiet", False))
        pipeline.save_results(result)

        logger.info(f"Pan-genome design complete! Results saved to: {args.output}")
        logger.info(f"Designed {result.primer_count} primers")
        logger.info(f"Mean enrichment: {result.mean_enrichment:.1f}x")

    except Exception as e:
        logger.error(f"Multi-genome pipeline failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


# =========================================================================
# CATEGORY 4: Handler function for simulation orphaned features
# =========================================================================


def run_simulate(args):
    """Monte Carlo amplification simulation"""
    import json
    from pathlib import Path

    from Bio import SeqIO

    from neoswga.core import reaction_conditions as rc
    from neoswga.core.genome_io import GenomeLoader
    from neoswga.core.replication_simulator import (
        Phi29Simulator,
        SimulationConfig,
        simulate_primer_set,
    )

    # Resolve primers from --primers or --from-results
    if args.from_results:
        import pandas as pd

        source = args.from_results
        if os.path.isdir(source):
            # Look for step4 CSV in the directory
            candidates = [
                os.path.join(source, "step4_improved_df.csv"),
                os.path.join(source, "step4_df.csv"),
            ]
            csv_path = None
            for c in candidates:
                if os.path.exists(c):
                    csv_path = c
                    break
            if csv_path is None:
                logger.error(f"No step4 CSV found in {source}")
                logger.error("Expected: step4_improved_df.csv or step4_df.csv")
                sys.exit(1)
        else:
            csv_path = source
            if not os.path.exists(csv_path):
                logger.error(f"File not found: {csv_path}")
                sys.exit(1)

        logger.info(f"Loading primers from {csv_path}")
        df = pd.read_csv(csv_path)
        if "primer" not in df.columns:
            logger.error(f"CSV file {csv_path} has no 'primer' column")
            sys.exit(1)
        primers = df["primer"].tolist()
        logger.info(f"Loaded {len(primers)} primers from results")
    else:
        primers = args.primers

    logger.info(f"Running amplification simulation for {len(primers)} primers")
    logger.info(f"Genome: {args.genome}")

    try:
        os.makedirs(args.output, exist_ok=True)

        # Load genome sequence
        logger.info("Loading genome...")
        loader = GenomeLoader()
        genome_sequence = loader.load_genome(args.genome)
        genome_length = len(genome_sequence)
        logger.info(f"  Genome length: {genome_length:,} bp")

        # Find primer positions in genome
        logger.info("Finding primer binding positions...")
        primer_positions = {}
        for primer in primers:
            primer_positions[primer] = {"forward": [], "reverse": []}
            primer_rc = primer.translate(str.maketrans("ACGT", "TGCA"))[::-1]

            # Find forward positions
            pos = 0
            while True:
                pos = genome_sequence.find(primer, pos)
                if pos == -1:
                    break
                primer_positions[primer]["forward"].append(pos)
                pos += 1

            # Find reverse positions
            pos = 0
            while True:
                pos = genome_sequence.find(primer_rc, pos)
                if pos == -1:
                    break
                primer_positions[primer]["reverse"].append(pos)
                pos += 1

            total = len(primer_positions[primer]["forward"]) + len(
                primer_positions[primer]["reverse"]
            )
            logger.info(f"  {primer}: {total} binding sites")

        # Get reaction conditions based on polymerase type
        polymerase = getattr(args, "polymerase", "phi29")

        # Get polymerase characteristics for optimal temperature
        from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS

        if polymerase not in POLYMERASE_CHARACTERISTICS:
            logger.error(f"Unknown polymerase: {polymerase}")
            logger.info(f"Available: {list(POLYMERASE_CHARACTERISTICS.keys())}")
            sys.exit(1)

        poly_info = POLYMERASE_CHARACTERISTICS[polymerase]
        optimal_temp = poly_info["optimal_temp"]
        processivity = poly_info["processivity"]

        logger.info(f"Polymerase: {polymerase}")
        logger.info(f"  Optimal temperature: {optimal_temp}C")
        logger.info(f"  Processivity: {processivity:,} bp")

        # Create conditions with appropriate temperature
        conditions = rc.ReactionConditions(temp=optimal_temp, polymerase=polymerase)

        # Configure simulation
        duration = args.duration * 60.0  # Convert minutes to seconds
        config = SimulationConfig(duration=duration, polymerase_type=polymerase)

        # Run simulation with replicates
        n_replicates = getattr(args, "replicates", 5)
        logger.info(f"\nRunning {n_replicates} simulation replicates...")

        results = simulate_primer_set(
            primers=primers,
            primer_positions=primer_positions,
            genome_length=genome_length,
            genome_sequence=genome_sequence,
            conditions=conditions,
            n_replicates=n_replicates,
            seed=getattr(args, "seed", None),
        )

        logger.info(f"\nSimulation complete!")
        logger.info(f"  Mean coverage: {results['mean_coverage']:.1%}")
        logger.info(f"  Std coverage: {results['std_coverage']:.1%}")
        logger.info(f"  Mean forks created: {results['mean_forks_created']:.0f}")
        logger.info(f"  Mean fork travel: {results['mean_fork_travel']:,.0f} bp")

        # Save results
        results_file = Path(args.output) / "simulation_results.json"
        with open(results_file, "w") as f:
            json.dump(
                {
                    "primers": primers,
                    "genome_length": genome_length,
                    "mean_coverage": results["mean_coverage"],
                    "std_coverage": results["std_coverage"],
                    "mean_forks_created": results["mean_forks_created"],
                    "mean_fork_travel": results["mean_fork_travel"],
                    "n_replicates": n_replicates,
                    "polymerase": polymerase,
                    "duration_seconds": duration,
                },
                f,
                indent=2,
            )
        logger.info(f"Results saved to: {results_file}")

        # Generate visualizations if requested
        if args.visualize:
            logger.info("Generating visualization plots...")
            try:
                from neoswga.core import simulation_plots

                plot_path = str(Path(args.output) / "simulation_plots.png")
                simulation_plots.plot_replication_summary(results, plot_path)
            except ImportError as e:
                logger.warning(f"Visualization unavailable (install neoswga[viz]): {e}")

        # Generate HTML report if requested
        if args.report:
            logger.info("Generating HTML report...")
            from neoswga.core import simulation_report

            report_path = str(Path(args.output) / "simulation_report.html")
            simulation_report.generate_replication_report(results, primers, report_path)

    except Exception as e:
        logger.error(f"Simulation failed: {e}")
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

    multi_genome_parser = subparsers.add_parser(
        "multi-genome", help="Pan-genome primer design for multiple targets"
    )
    multi_genome_parser.add_argument(
        "--genomes", required=True, nargs="+", help="List of target genome FASTA files"
    )
    multi_genome_parser.add_argument(
        "--background", nargs="+", help="Background genome FASTA files (e.g., host DNA)"
    )
    multi_genome_parser.add_argument(
        "--blacklist",
        nargs="+",
        help="Blacklist genome FASTA files (strongly avoid, e.g., other pathogens)",
    )
    multi_genome_parser.add_argument("--output", "-o", required=True, help="Output directory")
    multi_genome_parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducible design."
    )
    multi_genome_parser.add_argument(
        "--num-primers", type=int, default=12, help="Number of primers to design (default: 12)"
    )
    multi_genome_parser.add_argument(
        "--min-k", type=int, default=8, help="Minimum primer length (default: 8)"
    )
    multi_genome_parser.add_argument(
        "--max-k", type=int, default=12, help="Maximum primer length (default: 12)"
    )
    multi_genome_parser.add_argument(
        "--polymerase",
        choices=_polymerase_names(),
        help="Polymerase type: phi29 (30C), equiphi29 (42C), " "bst (63C), klenow (37C)",
    )
    multi_genome_parser.add_argument(
        "--validate-simulation", action="store_true", help="Validate with simulation"
    )
    multi_genome_parser.add_argument("--quiet", action="store_true", help="Minimal output")

    # =========================================================================
    # CATEGORY 4: Simulation command (orphaned features)
    # =========================================================================

    simulate_parser = subparsers.add_parser("simulate", help="Agent-based replication simulation")
    sim_input = simulate_parser.add_mutually_exclusive_group(required=True)
    sim_input.add_argument("--primers", nargs="+", help="Primer sequences to simulate")
    sim_input.add_argument(
        "--from-results",
        type=str,
        metavar="DIR_OR_CSV",
        help="Read primers from pipeline results directory or step4 CSV file "
        "(e.g., --from-results results/ or --from-results step4_improved_df.csv)",
    )
    simulate_parser.add_argument("--genome", required=True, help="Target genome FASTA file")
    simulate_parser.add_argument(
        "--output", "-o", required=True, help="Output directory for simulation results"
    )
    simulate_parser.add_argument(
        "--duration",
        "--cycles",
        type=int,
        default=60,
        help="Simulation duration in minutes (default: 60)",
    )
    simulate_parser.add_argument(
        "--replicates", type=int, default=5, help="Number of simulation replicates (default: 5)"
    )
    simulate_parser.add_argument(
        "--polymerase",
        choices=_polymerase_names(),
        default="phi29",
        help="Polymerase type: phi29 (30C), equiphi29 (42C), "
        "bst (63C), klenow (37C) (default: phi29)",
    )
    simulate_parser.add_argument(
        "--visualize", action="store_true", help="Generate visualization plots"
    )
    simulate_parser.add_argument("--report", action="store_true", help="Generate HTML report")
    simulate_parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible simulation (replicates derive distinct "
        "seeds so a run reproduces exactly).",
    )

    # Expand primers - add new primers to existing validated set
