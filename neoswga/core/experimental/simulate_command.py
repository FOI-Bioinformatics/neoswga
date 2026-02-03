#!/usr/bin/env python3
"""
CLI command for SWGA simulation.

Provides 'neoswga simulate' command to validate primer sets.

Usage:
    neoswga simulate --primers primers.txt --fg target.fasta --bg background.fasta \\
                     --positions positions.h5 --output report.html

    neoswga simulate --config simulation_config.json

    neoswga simulate --primers primers.txt --fg target.fasta --bg background.fasta \\
                     --positions positions.h5 --mode detailed --replicates 10
"""

import argparse
import json
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional
import time

from neoswga.core.swga_simulator import SwgaSimulator, SimulationResult
from neoswga.core.simulation_analysis import SimulationAnalyzer, print_analysis_report

logger = logging.getLogger(__name__)


def add_simulate_parser(subparsers):
    """Add simulate subcommand to parser"""
    sim_parser = subparsers.add_parser(
        'simulate',
        help='Simulate SWGA amplification with primer set',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Quick validation
  %(prog)s simulate --primers primers.txt --fg target.fasta --bg background.fasta \\
           --positions positions.h5

  # Detailed simulation with analysis
  %(prog)s simulate --primers primers.txt --fg target.fasta --bg background.fasta \\
           --positions positions.h5 --mode detailed --analyze --output report.html

  # Use config file
  %(prog)s simulate --config simulation_config.json

  # Compare multiple primer sets
  %(prog)s simulate --config compare_primers.json --compare
        '''
    )

    # Input files
    input_group = sim_parser.add_argument_group('Input Files')
    input_group.add_argument('--primers', '-p',
                            help='Primer sequences file (one per line or tab-delimited)')
    input_group.add_argument('--fg', '--target',
                            help='Target genome FASTA file')
    input_group.add_argument('--bg', '--background',
                            help='Background genome FASTA file')
    input_group.add_argument('--positions', '--pos',
                            help='HDF5 file with primer positions')
    input_group.add_argument('--config', '-c',
                            help='JSON config file (overrides other options)')

    # Simulation options
    sim_group = sim_parser.add_argument_group('Simulation Options')
    sim_group.add_argument('--mode', '-m',
                          choices=['fast', 'detailed', 'validation'],
                          default='fast',
                          help='Simulation mode [default: fast]')
    sim_group.add_argument('--replicates', '-r', type=int, default=5,
                          help='Number of replicates for detailed mode [default: 5]')
    sim_group.add_argument('--bin-size', type=int, default=10000,
                          help='Bin size for coverage analysis (bp) [default: 10000]')
    sim_group.add_argument('--max-extension', type=int, default=70000,
                          help='Max Phi29 extension distance (bp) [default: 70000]')

    # Reaction conditions
    cond_group = sim_parser.add_argument_group('Reaction Conditions')
    cond_group.add_argument('--temperature', '-t', type=float, default=30.0,
                           help='Reaction temperature (°C) [default: 30.0]')
    cond_group.add_argument('--duration', type=int, default=3600,
                           help='Reaction duration (seconds) [default: 3600]')
    cond_group.add_argument('--polymerase', choices=['phi29', 'equiphi29'],
                           default='phi29',
                           help='Polymerase type [default: phi29]')

    # Analysis options
    analysis_group = sim_parser.add_argument_group('Analysis Options')
    analysis_group.add_argument('--analyze', '-a', action='store_true',
                               help='Perform comprehensive analysis')
    analysis_group.add_argument('--plot', action='store_true',
                               help='Generate plots (requires --output)')
    analysis_group.add_argument('--compare', action='store_true',
                               help='Compare multiple primer sets')

    # Output options
    output_group = sim_parser.add_argument_group('Output Options')
    output_group.add_argument('--output', '-o',
                             help='Output file (HTML report or JSON results)')
    output_group.add_argument('--format', choices=['html', 'json', 'text'],
                             default='text',
                             help='Output format [default: text]')
    output_group.add_argument('--quiet', '-q', action='store_true',
                             help='Minimal output')
    output_group.add_argument('--verbose', '-v', action='store_true',
                             help='Verbose output')

    return sim_parser


def load_primers(primer_file: Path) -> List[str]:
    """Load primers from file"""
    primers = []
    with open(primer_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                # Handle tab-delimited format
                parts = line.split('\t')
                primer = parts[0].strip()
                if primer:
                    primers.append(primer)
    return primers


def load_config(config_file: Path) -> Dict:
    """Load configuration from JSON file"""
    with open(config_file) as f:
        return json.load(f)


def save_results(result: SimulationResult, output_path: Path, format: str):
    """Save simulation results"""
    if format == 'json':
        with open(output_path, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)
        logger.info(f"Results saved to: {output_path}")

    elif format == 'text':
        # Text report already printed to stdout
        if output_path:
            # Redirect to file
            pass

    elif format == 'html':
        # HTML report generation
        from neoswga.core.simulation_report import generate_html_report
        generate_html_report(result, output_path)
        logger.info(f"HTML report saved to: {output_path}")


def run_simulation(args) -> int:
    """
    Main simulation execution.

    Returns:
        0 on success, 1 on error
    """
    # Setup logging
    log_level = logging.WARNING if args.quiet else \
                logging.DEBUG if args.verbose else \
                logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

    # Load configuration
    if args.config:
        logger.info(f"Loading configuration from: {args.config}")
        config = load_config(Path(args.config))

        # Override with command-line args if provided
        primers_file = args.primers or config.get('primers')
        fg_genome = args.fg or config.get('target_genome') or config.get('fg_genome')
        bg_genome = args.bg or config.get('background_genome') or config.get('bg_genome')
        positions_file = args.positions or config.get('positions')

        mode = args.mode if args.mode != 'fast' else config.get('mode', 'fast')
        replicates = args.replicates if args.replicates != 5 else config.get('replicates', 5)
    else:
        # Use command-line arguments
        primers_file = args.primers
        fg_genome = args.fg
        bg_genome = args.bg
        positions_file = args.positions
        mode = args.mode
        replicates = args.replicates

    # Validate required arguments
    if not all([primers_file, fg_genome, bg_genome, positions_file]):
        logger.error("Missing required arguments: --primers, --fg, --bg, --positions")
        logger.error("Provide either all arguments or use --config")
        return 1

    # Check files exist
    for file_path, name in [
        (primers_file, "primers"),
        (fg_genome, "target genome"),
        (bg_genome, "background genome"),
        (positions_file, "positions")
    ]:
        if not Path(file_path).exists():
            logger.error(f"{name} file not found: {file_path}")
            return 1

    # Load primers
    logger.info(f"Loading primers from: {primers_file}")
    primers = load_primers(Path(primers_file))
    logger.info(f"Loaded {len(primers)} primers")

    if len(primers) == 0:
        logger.error("No primers found in file")
        return 1

    # Initialize simulator
    logger.info("Initializing simulator...")
    try:
        simulator = SwgaSimulator(
            primers=primers,
            fg_genome=fg_genome,
            bg_genome=bg_genome,
            fg_positions_h5=positions_file,
            bg_positions_h5=positions_file,
            reaction_conditions={
                'temperature': args.temperature,
                'polymerase': args.polymerase,
                'duration': args.duration
            },
            bin_size=args.bin_size,
            max_extension=args.max_extension
        )
    except Exception as e:
        logger.error(f"Failed to initialize simulator: {e}")
        import traceback
        traceback.print_exc()
        return 1

    # Run simulation
    logger.info(f"Running {mode} simulation...")
    start_time = time.time()

    try:
        if mode == 'fast':
            result = simulator.simulate_fast()
        elif mode == 'detailed':
            result = simulator.simulate_detailed(num_replicates=replicates)
        elif mode == 'validation':
            result = simulator.simulate_validation()
        else:
            logger.error(f"Unknown mode: {mode}")
            return 1
    except Exception as e:
        logger.error(f"Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

    elapsed = time.time() - start_time
    logger.info(f"Simulation completed in {elapsed:.1f} seconds")

    # Print basic report
    if not args.quiet:
        simulator.generate_report(result)

    # Comprehensive analysis
    if args.analyze:
        logger.info("Performing comprehensive analysis...")
        analyzer = SimulationAnalyzer(
            result=result,
            fg_positions=simulator.fg_positions,
            bg_positions=simulator.bg_positions,
            fg_length=simulator.fg_length,
            bg_length=simulator.bg_length,
            bin_size=args.bin_size
        )
        analysis = analyzer.analyze()

        if not args.quiet:
            print_analysis_report(analysis)

    # Generate plots
    if args.plot:
        if not args.output:
            logger.error("--plot requires --output to specify output file")
            return 1

        logger.info("Generating plots...")
        from neoswga.core.simulation_plots import generate_plots
        plot_file = Path(args.output).with_suffix('.png')
        generate_plots(result, simulator, analysis if args.analyze else None, plot_file)
        logger.info(f"Plots saved to: {plot_file}")

    # Save results
    if args.output:
        output_path = Path(args.output)
        save_results(result, output_path, args.format)

    # Print summary
    if not args.quiet:
        print("\n" + "=" * 80)
        print("SIMULATION SUMMARY")
        print("=" * 80)
        print(f"Mode: {mode}")
        print(f"Primers: {len(primers)}")
        print(f"Target coverage: {result.target_coverage:.1%}")
        print(f"Enrichment: {result.enrichment:.0f}×")
        print(f"Recommendation: {result.recommendation}")
        print(f"Runtime: {elapsed:.1f}s")
        print("=" * 80)

    # Exit code based on recommendation
    if result.recommendation in ['EXCELLENT', 'GOOD']:
        return 0
    elif result.recommendation == 'FAIR':
        logger.warning("Primer set quality is FAIR - consider improvements")
        return 0
    else:  # POOR
        logger.error("Primer set quality is POOR - redesign recommended")
        return 1


def main(args=None):
    """Main entry point"""
    if args is None:
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(dest='command')
        add_simulate_parser(subparsers)
        args = parser.parse_args()

    if args.command == 'simulate':
        return run_simulation(args)
    else:
        parser.print_help()
        return 1


if __name__ == '__main__':
    sys.exit(main())
