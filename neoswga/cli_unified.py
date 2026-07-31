#!/usr/bin/env python3
"""
Unified NeoSWGA Command Line Interface.

Single entry point for all SWGA primer design functionality:
- Standard pipeline (count-kmers, filter, score, optimize)
- Improved optimization (network-based, 10-100x better)
- Utility commands (build-filter, validate)
- Advanced features (optimize-conditions, analyze-set)

Usage:
    # Standard workflow
    neoswga count-kmers [options]
    neoswga filter [options]
    neoswga score [options]
    neoswga optimize [options] [--optimization-method=hybrid|dominating-set|network|background-aware|ensemble]

    # Utility commands
    neoswga build-filter <genome> <output_dir>
    neoswga validate [--quick]
    neoswga show-presets

    # Advanced features
    neoswga optimize-conditions --fg <genome> --output <dir>
    neoswga analyze-set --primers <seq1> <seq2> --fg <genome> --fg-kmers <prefix> --output <dir>

Examples:
    # Standard workflow
    neoswga count-kmers -j params.json
    neoswga filter -j params.json
    neoswga score -j params.json
    neoswga optimize -j params.json

    # Control optimization method (hybrid|dominating-set|network|background-aware|ensemble)
    neoswga optimize -j params.json --optimization-method=ensemble
    neoswga optimize -j params.json --use-background-filter --background-bloom-path filters/bg_bloom.pkl

    # Pre-build background filter (one-time, 30 min for human genome)
    neoswga build-filter human.fasta ./filters/

    # Show available reaction condition presets
    neoswga show-presets

    # Optimize conditions for a genome
    neoswga optimize-conditions --fg target.fasta --output results/

    # Analyze existing primer set
    neoswga analyze-set --primers ACGTACGTACGT TGCATGCATGCA --fg target.fasta --fg-kmers data/target --output results/
"""

import argparse
import json
import logging
import os
import sys

# Import StepPrerequisiteError for proper exception handling
# This is imported at module level so except clauses can catch it
try:
    from neoswga.core.pipeline import StepPrerequisiteError
except ImportError:
    # Define a fallback if pipeline module isn't available
    class StepPrerequisiteError(Exception):
        """Raised when a pipeline step's prerequisites are not satisfied."""

        pass


# Setup logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Shared CLI helpers now live in neoswga.cli._common. Re-exported here so the
# parser, the remaining handlers, and any external importers keep working.
from neoswga.cli._common import (  # noqa: E402,F401
    FORBIDDEN_PATH_CHARS,
    MAX_PRIMER_LENGTH,
    MIN_PRIMER_LENGTH,
    PRESETS,
    VALID_DNA_BASES,
    _apply_seed,
    _record_run_manifest,
    add_common_options,
    check_jellyfish_available,
    collect_primers_from_args,
    load_preset_conditions,
    load_primers_from_file,
    merge_args_to_parameter,
    setup_gpu_acceleration,
    validate_params_json_file,
    validate_path_security,
    validate_primer_file,
    validate_primer_sequence,
)

COMMAND_GROUPS = [
    (
        "Pipeline",
        [
            "count-kmers",
            "filter",
            "score",
            "optimize",
            "design",
        ],
    ),
    (
        "Setup",
        [
            "init",
            "start",
            "validate",
            "validate-params",
            "schema",
            "doctor",
            "show-presets",
            "suggest",
        ],
    ),
    (
        "Results",
        [
            "interpret",
            "report",
            "export",
        ],
    ),
    (
        "Analysis",
        [
            "analyze-set",
            "analyze-genome",
            "analyze-dimers",
            "analyze-stability",
            "analyze-coverage",
            "evaluate-set",
        ],
    ),
    (
        "Advanced",
        [
            "multi-genome",
            "simulate",
            "optimize-conditions",
            "build-filter",
            "background-list",
            "background-add",
            "genome-add",
            "genome-list",
            "genome-remove",
        ],
    ),
    (
        "Experimental",
        [
            "expand-primers",
            "swap-primer",
            "contract-set",
            "rescore-set",
            "predict-efficiency",
            "validate-model",
        ],
    ),
]


class GroupedHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Formatter that displays subcommands under category headers."""

    def _format_action(self, action):
        # Only customise the subparsers listing
        if not isinstance(action, argparse._SubParsersAction):
            return super()._format_action(action)

        parts = []
        # Build a lookup from command name to its help string
        cmd_help = {}
        for choice_action in action._choices_actions:
            cmd_help[choice_action.dest] = choice_action.help or ""

        for group_name, commands in COMMAND_GROUPS:
            parts.append(f"\n  {group_name}:")
            for cmd in commands:
                help_text = cmd_help.get(cmd, "")
                parts.append(f"    {cmd:<24s}{help_text}")

        parts.append("")
        return "\n".join(parts) + "\n"


def create_parser():
    """Create unified argument parser"""

    from neoswga import __version__

    parser = argparse.ArgumentParser(
        prog="neoswga",
        description="NeoSWGA - Selective Whole Genome Amplification primer design",
        formatter_class=GroupedHelpFormatter,
        epilog="""Quick Start:
  neoswga init --genome target.fasta --background host.fasta
  neoswga count-kmers -j params.json
  neoswga filter -j params.json
  neoswga score -j params.json
  neoswga optimize -j params.json

Run "neoswga <command> --help" for details on a specific command.
        """,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", title="commands", metavar="<command>")

    # Register each command group's subparsers (defined alongside their
    # handlers in neoswga/cli/<group>.py).
    from neoswga.cli import (
        analysis,
        commands,
        evaluate,
        iterate,
        pipeline,
        registry,
        report,
        setup,
        simulate,
    )

    pipeline.add_parsers(subparsers)
    evaluate.add_parsers(subparsers)
    setup.add_parsers(subparsers)
    report.add_parsers(subparsers)
    analysis.add_parsers(subparsers)
    simulate.add_parsers(subparsers)
    iterate.add_parsers(subparsers)
    registry.add_parsers(subparsers)
    commands.add_parsers(subparsers)

    # Validate that all registered commands appear in COMMAND_GROUPS
    grouped_cmds = {cmd for _, cmds in COMMAND_GROUPS for cmd in cmds}
    registered_cmds = set(subparsers.choices.keys()) if hasattr(subparsers, "choices") else set()
    missing = registered_cmds - grouped_cmds
    if missing:
        import warnings

        warnings.warn(
            f"Commands not listed in COMMAND_GROUPS: {', '.join(sorted(missing))}",
            stacklevel=2,
        )

    return parser


# Remaining standalone command handlers live in neoswga.cli.commands.
from neoswga.cli.analysis import (  # noqa: E402,F401
    analyze_primer_set,
    optimize_conditions,
    run_analyze_dimers,
    run_analyze_genome,
    run_analyze_stability,
    show_presets,
)
from neoswga.cli.commands import (  # noqa: E402,F401
    run_analyze_coverage,
    run_design,
    run_predict_efficiency,
    run_start,
    run_suggest,
)
from neoswga.cli.evaluate import run_evaluate_set  # noqa: E402,F401
from neoswga.cli.iterate import (  # noqa: E402,F401
    run_contract_set,
    run_expand_primers,
    run_rescore_set,
    run_swap_primer,
)

# Pipeline-step handlers live in neoswga.cli.pipeline (re-exported for dispatch
# and for external importers, e.g. test_docs_consistency imports
# print_method_guide).
from neoswga.cli.pipeline import (  # noqa: E402,F401
    print_method_guide,
    run_build_filter,
    run_step1,
    run_step2,
    run_step3,
    run_step4,
)
from neoswga.cli.registry import (  # noqa: E402,F401
    run_background_add,
    run_background_list,
    run_genome_add,
    run_genome_list,
    run_genome_remove,
)
from neoswga.cli.report import run_export, run_report  # noqa: E402,F401

# Setup/validation and reporting handlers live in neoswga.cli.* (re-exported
# for the dispatch table).
from neoswga.cli.setup import (  # noqa: E402,F401
    run_doctor,
    run_init,
    run_interpret,
    run_schema,
    run_validate,
    run_validate_model,
    run_validate_params,
)
from neoswga.cli.simulate import run_multi_genome, run_simulate  # noqa: E402,F401


def main():
    """Main entry point"""
    parser = create_parser()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        print(
            "\nTip: New to neoswga? Start with 'neoswga init --genome target.fasta'",
            file=sys.stderr,
        )
        sys.exit(1)

    # Set verbosity
    if getattr(args, "quiet", False):
        logging.getLogger().setLevel(logging.WARNING)
    elif getattr(args, "verbose", False):
        logging.getLogger().setLevel(logging.DEBUG)

    # Dispatch to appropriate command
    commands = {
        # Standard pipeline
        "count-kmers": run_step1,
        "filter": run_step2,
        "score": run_step3,
        "optimize": run_step4,
        # Unified pipeline
        "design": run_design,
        # Utility commands
        "build-filter": run_build_filter,
        "validate": run_validate,
        "show-presets": lambda args: show_presets(),
        "init": run_init,
        "validate-params": run_validate_params,
        "schema": run_schema,
        "doctor": run_doctor,
        "validate-model": run_validate_model,
        "interpret": run_interpret,
        "report": run_report,
        "export": run_export,
        "start": run_start,
        "suggest": run_suggest,
        "optimize-conditions": optimize_conditions,
        "analyze-set": analyze_primer_set,
        # Category 1: Orphaned analysis features (now exposed!)
        "analyze-genome": run_analyze_genome,
        "analyze-dimers": run_analyze_dimers,
        "analyze-stability": run_analyze_stability,
        "analyze-coverage": run_analyze_coverage,
        "evaluate-set": run_evaluate_set,
        # Category 3: Orphaned pipeline features (now exposed!)
        "multi-genome": run_multi_genome,
        # Category 4: Orphaned simulation features (now exposed!)
        "simulate": run_simulate,
        # Category 6: Iterative design
        "expand-primers": run_expand_primers,
        "swap-primer": run_swap_primer,
        "contract-set": run_contract_set,
        "rescore-set": run_rescore_set,
        "predict-efficiency": run_predict_efficiency,
        # Category 7: Background registry
        "background-list": run_background_list,
        "background-add": run_background_add,
        # Category 8: Genome library
        "genome-add": run_genome_add,
        "genome-list": run_genome_list,
        "genome-remove": run_genome_remove,
    }

    command_func = commands.get(args.command)
    if command_func:
        try:
            command_func(args)
        except KeyboardInterrupt:
            logger.info("\nInterrupted by user")
            sys.exit(1)
        except FileNotFoundError as e:
            logger.error(f"File not found: {e}")
            logger.error("Check that all genome files and data paths exist.")
            sys.exit(1)
        except KeyError as e:
            logger.error(f"Missing required parameter: {e}")
            logger.error(
                "Run 'neoswga init --genome target.fasta' to create a complete params.json"
            )
            logger.error(
                "Or run 'neoswga validate-params -j params.json' to check your configuration"
            )
            sys.exit(1)
        except (ValueError, TypeError) as e:
            logger.error(f"Invalid input: {e}")
            if getattr(args, "verbose", False):
                import traceback

                traceback.print_exc()
            sys.exit(1)
        except PermissionError as e:
            logger.error(f"Permission denied: {e}")
            logger.error("Check that the output directory is writable.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Command failed: {e}")
            if getattr(args, "verbose", False):
                import traceback

                traceback.print_exc()
            else:
                logger.error("Run with --verbose for full traceback")
            sys.exit(1)
    else:
        logger.error(f"Unknown command: {args.command}")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
