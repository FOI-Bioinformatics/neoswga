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

    # =========================================================================
    # STEP 1: K-mer preprocessing
    # =========================================================================
    count_kmers_parser = subparsers.add_parser(
        "count-kmers", help="K-mer preprocessing with primer length control"
    )
    add_common_options(count_kmers_parser)
    count_kmers_parser.add_argument("-x", "--fasta-fore", help="Target genome FASTA")
    count_kmers_parser.add_argument("-y", "--fasta-back", help="Background genome FASTA")
    count_kmers_parser.add_argument("-k", "--kmer-fore", help="Target k-mer prefix")
    count_kmers_parser.add_argument("-l", "--kmer-back", help="Background k-mer prefix")

    # Primer length control (defaults from params.json, or 6-12 if not specified)
    count_kmers_parser.add_argument(
        "--min-k",
        type=int,
        default=None,
        help="Minimum primer length (default: from params.json or 6). Supported range: 6-18bp",
    )
    count_kmers_parser.add_argument(
        "--max-k",
        type=int,
        default=None,
        help="Maximum primer length (default: from params.json or 12). Use 15-18bp with DMSO/betaine additives",
    )
    count_kmers_parser.add_argument(
        "--exclusion-genome",
        type=str,
        default=None,
        help="Exclusion genome FASTA file (e.g., mtDNA)",
    )
    count_kmers_parser.add_argument(
        "--blacklist",
        "-bl",
        nargs="+",
        default=None,
        help="Blacklist genome FASTA file(s) (penalty-weighted filtering)",
    )
    count_kmers_parser.add_argument(
        "--bl-penalty", type=float, default=None, help="Blacklist penalty weight (default: 5.0)"
    )
    count_kmers_parser.add_argument(
        "--max-bl-freq",
        type=float,
        default=None,
        help="Maximum blacklist frequency (default: 0.0, any hit rejects)",
    )

    # =========================================================================
    # STEP 2: Candidate filtering
    # =========================================================================
    filter_parser = subparsers.add_parser(
        "filter", help="Candidate primer filtering (adaptive GC + reaction conditions)"
    )
    add_common_options(filter_parser)

    # GC Filtering Options
    step2_gc_group = filter_parser.add_argument_group("GC Content Filtering")
    step2_gc_group.add_argument(
        "--gc-tolerance",
        type=float,
        default=0.15,
        help="GC tolerance for adaptive filter (default: 0.15)",
    )
    step2_gc_group.add_argument(
        "--gc-min", type=float, help="Explicit minimum GC content (overrides adaptive)"
    )
    step2_gc_group.add_argument(
        "--gc-max", type=float, help="Explicit maximum GC content (overrides adaptive)"
    )

    # Reaction Conditions
    step2_rxn_group = filter_parser.add_argument_group("Reaction Conditions")
    step2_rxn_group.add_argument(
        "--reaction-temp",
        type=float,
        help="Reaction temperature in C (default: from params.json or polymerase preset)",
    )
    step2_rxn_group.add_argument(
        "--na-conc",
        type=float,
        help="Sodium concentration in mM (default: from params.json or 50.0)",
    )

    # Additives
    step2_add_group = filter_parser.add_argument_group(
        "Additives (enable longer primers & GC-extreme genomes)"
    )
    step2_add_group.add_argument(
        "--dmso-percent",
        type=float,
        help="DMSO concentration 0-10%% (Tm lowering, secondary structure reduction)",
    )
    step2_add_group.add_argument(
        "--betaine-m",
        type=float,
        help="Betaine concentration 0-2.5 M (equalizes AT/GC, enables longer primers)",
    )
    step2_add_group.add_argument(
        "--trehalose-m", type=float, help="Trehalose concentration 0-1.0 M (Tm lowering)"
    )
    step2_add_group.add_argument(
        "--glycerol-percent", type=float, help="Glycerol concentration 0-15%% (enzyme stabilizer)"
    )
    step2_add_group.add_argument(
        "--bsa", type=float, help="BSA concentration 0-400 ug/mL (inhibitor neutralizer)"
    )
    step2_add_group.add_argument(
        "--peg-percent", type=float, help="PEG concentration 0-15%% (molecular crowding)"
    )
    step2_add_group.add_argument(
        "--mg-conc", type=float, help="Mg2+ concentration in mM (auto-optimized if not specified)"
    )
    step2_add_group.add_argument(
        "--ssb",
        action="store_true",
        help="Use single-strand binding protein (lowers effective annealing temp)",
    )
    step2_add_group.add_argument(
        "--ethanol-percent",
        type=float,
        help="Ethanol concentration 0-5%% (secondary structure reduction)",
    )
    step2_add_group.add_argument(
        "--urea-m", type=float, help="Urea concentration 0-2.0 M (denatures GC-rich regions)"
    )
    step2_add_group.add_argument(
        "--tmac-m", type=float, help="TMAC concentration 0-0.1 M (equalizes AT/GC Tm)"
    )
    step2_add_group.add_argument(
        "--formamide-percent", type=float, help="Formamide concentration 0-10%% (Tm lowering)"
    )

    # Preset Conditions
    filter_parser.add_argument(
        "--preset",
        choices=[
            "standard_phi29",
            "enhanced_equiphi29",
            "high_gc_genome",
            "long_primers_15mer",
            "q_solution",
            "gc_melt",
            "crude_sample",
            "low_temp",
            "bst",
            "klenow",
            "extreme_gc",
        ],
        help="Use predefined reaction conditions preset",
    )

    # Traditional Filtering Parameters
    step2_trad_group = filter_parser.add_argument_group("Traditional Filtering")
    step2_trad_group.add_argument(
        "--min-fg-freq", type=float, help="Minimum target genome frequency (default: 1e-5)"
    )
    step2_trad_group.add_argument(
        "--max-bg-freq", type=float, help="Maximum background genome frequency (default: 5e-6)"
    )
    step2_trad_group.add_argument(
        "--max-gini", type=float, help="Maximum Gini index for evenness (default: 0.6)"
    )
    step2_trad_group.add_argument(
        "--max-primer", type=int, help="Number of top primers to keep (default: 500)"
    )
    step2_trad_group.add_argument(
        "--min-tm", type=float, help="Minimum melting temperature in C (default: 15)"
    )
    step2_trad_group.add_argument(
        "--max-tm", type=float, help="Maximum melting temperature in C (default: 45)"
    )
    step2_trad_group.add_argument(
        "--max-dimer-bp", type=int, help="Maximum heterodimer base pairs (default: 3)"
    )
    step2_trad_group.add_argument(
        "--max-self-dimer-bp", type=int, help="Maximum self-dimer base pairs (default: 4)"
    )

    # Background Bloom Filter (for large background genomes)
    step2_bloom_group = filter_parser.add_argument_group("Background Filtering (for large genomes)")
    step2_bloom_group.add_argument(
        "--use-bloom-filter",
        action="store_true",
        default=False,
        help="Use Bloom filter for background filtering (memory-efficient for human genome)",
    )
    step2_bloom_group.add_argument(
        "--bloom-filter-path",
        type=str,
        help="Path to pre-built Bloom filter (.pkl). Build with: neoswga build-filter",
    )
    step2_bloom_group.add_argument(
        "--sampled-index-path",
        type=str,
        help="Path to pre-built sampled index (.pkl) for count estimation",
    )
    step2_bloom_group.add_argument(
        "--bloom-max-bg-matches",
        type=int,
        default=10,
        help="Maximum background matches via Bloom filter (default: 10)",
    )

    # Exclusion Genome Filtering (zero-tolerance for contaminants)
    step2_excl_group = filter_parser.add_argument_group(
        "Exclusion Genome Filtering (reject primers binding contaminant sequences)"
    )
    step2_excl_group.add_argument(
        "--exclusion-genome",
        type=str,
        default=None,
        help="Exclusion genome FASTA file (e.g., mtDNA, chloroplast). "
        "Primers binding this genome are rejected.",
    )
    step2_excl_group.add_argument(
        "--excl-threshold",
        type=int,
        default=0,
        help="Maximum allowed hits in exclusion genome " "(default: 0 = any hit rejects primer)",
    )

    # Blacklist Genome Filtering (penalty-weighted)
    step2_bl_group = filter_parser.add_argument_group(
        "Blacklist Genome Filtering (penalty-weighted filtering of contaminating sequences)"
    )
    step2_bl_group.add_argument(
        "--blacklist", "-bl", nargs="+", default=None, help="Blacklist genome FASTA file(s)"
    )
    step2_bl_group.add_argument(
        "--bl-penalty", type=float, default=None, help="Blacklist penalty weight (default: 5.0)"
    )
    step2_bl_group.add_argument(
        "--max-bl-freq",
        type=float,
        default=None,
        help="Maximum blacklist frequency (default: 0.0, any hit rejects)",
    )

    # =========================================================================
    # STEP 3: Random forest scoring
    # =========================================================================
    score_parser = subparsers.add_parser(
        "score", help="Amplification efficacy scoring with ML threshold control"
    )
    add_common_options(score_parser)

    score_parser.add_argument(
        "--min-amp-pred", type=float, help="Minimum amplification prediction score (default: 10)"
    )
    score_parser.add_argument(
        "--full-score",
        action="store_true",
        help="Include thermodynamic delta-G histogram features in "
        "random-forest scoring. These features contribute <2%% of "
        "model accuracy but >99%% of scoring compute time, so they "
        "are skipped by default. Pass this flag only if you need "
        "the full histogram output for downstream analysis.",
    )
    score_parser.add_argument(
        "--fast-score",
        action="store_true",
        help="Deprecated alias for the current default behavior "
        "(thermodynamic histogram features are skipped). "
        "Accepted for backwards compatibility.",
    )
    score_parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible scoring. K-mer "
        "sampling (on by default) is otherwise "
        "non-deterministic; set this for repeatable "
        "step3 output.",
    )

    # Enhanced feature engineering options
    score_parser.add_argument(
        "--use-enhanced-features",
        action="store_true",
        help="Use enhanced 120+ feature model (requires enhanced_rf_model.pkl)",
    )
    score_parser.add_argument(
        "--enhanced-model-path", type=str, default=None, help="Path to enhanced random forest model"
    )

    # =========================================================================
    # STEP 4: Primer set optimization
    # =========================================================================
    optimize_parser = subparsers.add_parser(
        "optimize", help="Primer set optimization (network-based + experimental)"
    )
    add_common_options(optimize_parser)

    # Method Selection
    opt_method_group = optimize_parser.add_argument_group("Method Selection")
    opt_method_group.add_argument(
        "-m",
        "--optimization-method",
        choices=["hybrid", "dominating-set", "network", "background-aware", "ensemble"],
        default="hybrid",
        help="Optimization method. "
        "Decision tree: "
        "hybrid (default, general use), "
        "dominating-set (speed-critical, large pools), "
        "background-aware (clinical, 10-20x bg reduction), "
        "network (Tm-weighted, dimer-aware), "
        "ensemble (run several and keep the best by "
        "normalized score). "
        "Use --method-guide for detailed comparison.",
    )
    opt_method_group.add_argument(
        "--ensemble-methods",
        nargs="+",
        default=None,
        metavar="METHOD",
        help="Methods to run when --optimization-method=ensemble "
        "(default: hybrid dominating-set network background-aware). "
        "The best result by application-weighted normalized score "
        "is kept; a per-method comparison table is printed.",
    )
    opt_method_group.add_argument(
        "--ensemble-combine",
        choices=["best", "union"],
        default="best",
        help="How ensemble combines methods: 'best' keeps the single best set; "
        "'union' additionally re-optimizes over the pooled primers from all "
        "methods (can beat any single method; never worsens it).",
    )
    opt_method_group.add_argument(
        "--scoring-weights",
        dest="scoring_weights",
        choices=["clinical", "discovery", "fast", "balanced", "enrichment"],
        default="balanced",
        help="Scoring-weight preset for the network optimizer: "
        "clinical (high specificity), discovery (max coverage), "
        "fast (quick screening), balanced (equal weights), "
        "enrichment (sequencing).",
    )
    opt_method_group.add_argument(
        "--method-guide",
        action="store_true",
        help="Show optimization method selection guide and exit",
    )

    # Primer Strategy
    opt_strategy_group = optimize_parser.add_argument_group("Primer Strategy")
    opt_strategy_group.add_argument(
        "--use-cooperative-binding",
        action="store_true",
        help="[EXPERIMENTAL] Cooperative binding model (not yet fully integrated)",
    )
    opt_strategy_group.add_argument(
        "--primer-strategy",
        choices=["standard", "hybrid"],
        default="standard",
        help="Primer design strategy (standard: uniform length, hybrid: mixed lengths)",
    )

    # Background Filtering
    opt_bg_group = optimize_parser.add_argument_group("Background Filtering")
    opt_bg_group.add_argument(
        "--use-background-filter",
        action="store_true",
        default=False,
        help="Use Bloom filter for background filtering",
    )
    opt_bg_group.add_argument(
        "--no-bg-prefilter",
        action="store_true",
        default=False,
        help="Disable automatic background pre-filtering of candidates "
        "(enabled by default when background genome data is available)",
    )
    opt_bg_group.add_argument(
        "--background-bloom-path", type=str, help="Path to pre-built Bloom filter"
    )
    opt_bg_group.add_argument(
        "--background-sampled-path", type=str, help="Path to pre-built sampled index"
    )
    opt_bg_group.add_argument(
        "--no-background",
        action="store_true",
        default=False,
        help="Host-free mode: optimize without background genome data. "
        "Relies on intrinsic primer quality (Tm, complexity, evenness) "
        "rather than fg/bg selectivity. Use when background genome is "
        "unknown or unavailable.",
    )

    # Performance
    opt_perf_group = optimize_parser.add_argument_group("Performance")
    opt_perf_group.add_argument(
        "--use-position-cache",
        action="store_true",
        default=True,
        help="Use in-memory position cache (default: True, 1000x speedup)",
    )
    opt_perf_group.add_argument(
        "--no-position-cache",
        action="store_false",
        dest="use_position_cache",
        help="Disable position cache (slower)",
    )
    opt_perf_group.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible results. Affects the network "
        "refinement, ensemble re-seeding, RF k-mer sampling, and any "
        "background pre-filtering. Required for clinical/regulated workflows.",
    )
    opt_perf_group.add_argument(
        "-n",
        "--num-primers",
        type=int,
        help="Number of primers to select (default: from params.json or 6)",
    )
    opt_perf_group.add_argument(
        "--max-optimization-time",
        type=int,
        default=300,
        help="Maximum optimization time in seconds (default: 300)",
    )
    opt_perf_group.add_argument(
        "--max-extension",
        type=int,
        default=70000,
        help="Maximum Phi29 extension length in bp (default: 70000)",
    )

    # Set Size & Application
    opt_size_group = optimize_parser.add_argument_group("Set Size & Application")
    opt_size_group.add_argument(
        "--auto-size",
        action="store_true",
        help="Automatically determine optimal primer set size based on "
        "application profile and reaction conditions",
    )
    opt_size_group.add_argument(
        "--application",
        type=str,
        choices=["balanced", "discovery", "clinical", "enrichment", "metagenomics"],
        default="enrichment",
        help="Application profile: "
        "balanced (default scoring), "
        "discovery (maximize sensitivity), "
        "clinical (minimize false positives), "
        "enrichment (balanced, default for auto-size), "
        "metagenomics (capture diversity). "
        "Drives --auto-size AND tunes normalized_score weights.",
    )
    opt_size_group.add_argument(
        "--min-fg-bg-ratio",
        type=float,
        help="Minimum foreground/background binding site ratio. "
        "Overrides application profile default. Higher values = more selective.",
    )
    opt_size_group.add_argument(
        "--show-frontier",
        action="store_true",
        help="Show Pareto frontier of coverage vs fg/bg ratio tradeoffs "
        "(requires matplotlib for plotting)",
    )
    opt_size_group.add_argument(
        "--quick-estimate",
        action="store_true",
        help="Use quick estimation only for auto-size (skip full optimization at multiple sizes)",
    )

    # Mechanistic Model
    opt_mech_group = optimize_parser.add_argument_group("Mechanistic Model")
    opt_mech_group.add_argument(
        "--use-mechanistic-model",
        action="store_true",
        help="Use mechanistic model for primer weighting during optimization",
    )
    opt_mech_group.add_argument(
        "--mechanistic-weight",
        type=float,
        default=None,
        help="Weight for mechanistic model in scoring (0.0-1.0, default 0.3 when "
        "enabled). Setting this > 0 implies --use-mechanistic-model.",
    )
    opt_mech_group.add_argument(
        "--template-gc",
        type=float,
        help="Template genome GC content (0-1). Auto-detected if not specified.",
    )
    opt_mech_group.add_argument(
        "--uniformity-weight",
        type=float,
        default=0.0,
        help="Weight for coverage uniformity in scoring (0.0-1.0, default: 0.0). "
        "Higher values prioritize even genome coverage over raw enrichment.",
    )

    # Post-processing
    opt_post_group = optimize_parser.add_argument_group("Post-processing")
    opt_post_group.add_argument(
        "--minimize-primers",
        action="store_true",
        help="Post-process to minimize primer count while maintaining coverage",
    )
    opt_post_group.add_argument(
        "--target-coverage",
        type=float,
        default=0.70,
        help="Target genome coverage fraction when minimizing primers (default: 0.70)",
    )
    opt_post_group.add_argument(
        "--min-per-target-coverage",
        type=float,
        default=0.0,
        help="Multi-genome runs only: minimum coverage required on every "
        "individual target. Below this, the post-optimization validator "
        "flags a warning. Default 0.0 (disabled).",
    )

    # Validation
    opt_val_group = optimize_parser.add_argument_group("Validation")
    opt_val_group.add_argument(
        "--validate-simulation",
        action="store_true",
        help="Validate results with stochastic simulation (Gillespie algorithm)",
    )
    opt_val_group.add_argument(
        "--simulation-time",
        type=float,
        default=3600.0,
        help="Simulation time in seconds (default: 3600)",
    )

    # =========================================================================
    # UNIFIED: Complete pipeline (all 4 steps)
    # =========================================================================
    design_parser = subparsers.add_parser(
        "design", help="Run complete primer design pipeline (all 4 steps sequentially)"
    )
    add_common_options(design_parser)

    # Multi-genome option (merges multi-genome)
    design_parser.add_argument(
        "--multi-genome",
        nargs="+",
        metavar="GENOME",
        help="Design pan-genome primers for multiple target genomes",
    )
    design_parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        help="Minimum fraction of genomes to cover (default: 0.8)",
    )

    # Step control
    design_parser.add_argument(
        "--start-from",
        type=int,
        choices=[1, 2, 3, 4],
        help="Start from step: 1=count-kmers, 2=filter, 3=score, 4=optimize",
    )
    design_parser.add_argument(
        "--stop-at",
        type=int,
        choices=[1, 2, 3, 4],
        help="Stop at step: 1=count-kmers, 2=filter, 3=score, 4=optimize",
    )

    # =========================================================================
    # UTILITY: Build background filter
    # =========================================================================
    build_filter_parser = subparsers.add_parser(
        "build-filter", help="Build background Bloom filter (one-time setup)"
    )
    build_filter_parser.add_argument(
        "--genome",
        required=True,
        help="Path to background genome FASTA (or k-mer prefix with --from-kmers)",
    )
    build_filter_parser.add_argument(
        "-o",
        "--output",
        "--output-dir",
        required=True,
        dest="output_dir",
        help="Output directory for filter files",
    )
    build_filter_parser.add_argument(
        "--from-kmers",
        action="store_true",
        help="Build from jellyfish k-mer files (MUCH faster for large genomes)",
    )
    build_filter_parser.add_argument(
        "--min-k", type=int, default=6, help="Minimum k-mer length (default: 6)"
    )
    build_filter_parser.add_argument(
        "--max-k", type=int, default=12, help="Maximum k-mer length (default: 12)"
    )
    build_filter_parser.add_argument(
        "--force", action="store_true", help="Rebuild even if filter exists"
    )
    build_filter_parser.add_argument(
        "--capacity", type=int, help="Bloom filter capacity (default: auto from genome size)"
    )
    build_filter_parser.add_argument(
        "--error-rate", type=float, default=0.01, help="Bloom filter error rate (default: 0.01)"
    )
    build_filter_parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    build_filter_parser.add_argument("-q", "--quiet", action="store_true", help="Minimal output")

    # =========================================================================
    # UTILITY: Validate (parent command with install/params/model subcommands)
    # =========================================================================
    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate the installation, a params.json file, or the "
        "mechanistic model. Run without a subcommand for installation "
        "checks (backwards-compatible).",
    )
    # Backwards-compat flags (when invoked as `neoswga validate [--quick|--all]`
    # with no subcommand — preserves the historical installation-check shape).
    validate_parser.add_argument(
        "--quick", action="store_true", help="Run quick installation validation (subset of tests)"
    )
    validate_parser.add_argument(
        "--all", action="store_true", help="Run all installation tests including benchmarks"
    )

    validate_sub = validate_parser.add_subparsers(
        dest="validate_mode",
        metavar="{install,params,model}",
        title="subcommands",
    )

    validate_install_sub = validate_sub.add_parser(
        "install",
        help="Check that Jellyfish, Python dependencies, and ML models load.",
    )
    validate_install_sub.add_argument(
        "--quick", action="store_true", help="Run quick validation (subset of tests)"
    )
    validate_install_sub.add_argument(
        "--all", action="store_true", help="Run all tests including benchmarks"
    )

    validate_params_sub = validate_sub.add_parser(
        "params",
        help="Validate a params.json file against the schema.",
    )
    validate_params_sub.add_argument(
        "-j", "--json-file", required=True, help="Parameters JSON file to validate"
    )

    validate_model_sub = validate_sub.add_parser(
        "model",
        help="Run regression tests on the mechanistic model.",
    )
    validate_model_sub.add_argument(
        "--verbose", "-v", action="store_true", help="Show detailed results for each test"
    )
    validate_model_sub.add_argument(
        "--output-json", action="store_true", help="Output results as JSON"
    )

    # =========================================================================
    # UTILITY: Show presets
    # =========================================================================
    subparsers.add_parser("show-presets", help="Show available reaction condition presets")

    # =========================================================================
    # SETUP: Initialize new project (setup wizard)
    # =========================================================================
    init_parser = subparsers.add_parser(
        "init",
        help="Set up a new primer design project with guided configuration. "
        "Requires a target genome path; use `neoswga start` for a "
        "menu-driven entry point that walks you through choosing a "
        "command first.",
    )
    init_parser.add_argument("--genome", "-g", required=True, help="Target genome FASTA file")
    init_parser.add_argument("--background", "-b", help="Background genome FASTA file (optional)")
    init_parser.add_argument(
        "--output",
        "-o",
        default="params.json",
        help="Output params.json path (default: params.json)",
    )
    init_parser.add_argument(
        "--output-dir",
        default="results",
        help="Directory for pipeline output files (default: results)",
    )
    init_parser.add_argument(
        "--advanced", "-a", action="store_true", help="Show advanced configuration options"
    )
    init_parser.add_argument(
        "--yes", "-y", action="store_true", help="Auto-approve without prompting"
    )
    init_parser.add_argument(
        "--non-interactive", action="store_true", help="Run non-interactively with defaults"
    )
    init_parser.add_argument(
        "--blacklist",
        "-bl",
        nargs="+",
        default=None,
        help="Blacklist genome FASTA file(s) for contamination filtering",
    )

    # =========================================================================
    # SETUP: Validate parameters
    # =========================================================================
    validate_params_parser = subparsers.add_parser(
        "validate-params", help="Validate params.json configuration before running"
    )
    validate_params_parser.add_argument(
        "-j", "--json-file", required=True, help="Parameters JSON file to validate"
    )

    # =========================================================================
    # SETUP: Dump canonical JSON schema
    # =========================================================================
    schema_parser = subparsers.add_parser(
        "schema", help="Inspect or dump the canonical params.json schema"
    )
    schema_parser.add_argument(
        "--dump", action="store_true", help="Print the params.json JSON Schema to stdout"
    )
    schema_parser.add_argument(
        "--output",
        "-o",
        type=str,
        default=None,
        help="Write the schema to the given file instead of stdout",
    )

    # =========================================================================
    # SETUP: Diagnostic health check
    # =========================================================================
    doctor_parser = subparsers.add_parser(
        "doctor",
        help="Run a diagnostic health check on the installation, "
        "params.json, and optimizer capability matrix",
    )
    doctor_parser.add_argument(
        "-j",
        "--json-file",
        type=str,
        default=None,
        help="Optional params.json to include in the diagnostic",
    )
    doctor_parser.add_argument(
        "--json",
        dest="output_json",
        action="store_true",
        help="Emit the diagnostic report as JSON instead of text",
    )

    # =========================================================================
    # SETUP: Validate mechanistic model
    # =========================================================================
    validate_model_parser = subparsers.add_parser(
        "validate-model", help="Validate mechanistic model against expected behavior"
    )
    validate_model_parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show detailed results for each test"
    )
    validate_model_parser.add_argument(
        "--output-json", action="store_true", help="Output results as JSON"
    )

    # =========================================================================
    # SETUP: Interpret results
    # =========================================================================
    interpret_parser = subparsers.add_parser(
        "interpret", help="Interpret pipeline results and provide quality assessment"
    )
    interpret_parser.add_argument(
        "-d", "--dir", required=True, help="Results directory containing step4_improved_df.csv"
    )

    # =========================================================================
    # SETUP: Generate report
    # =========================================================================
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
    subparsers.add_parser(
        "start",
        help="Interactive menu to discover and launch neoswga features. "
        "If you already know you want to create a params.json for a "
        "specific genome, use `neoswga init -g GENOME` directly.",
    )

    # =========================================================================
    # SETUP: Suggest reaction conditions
    # =========================================================================
    suggest_parser = subparsers.add_parser(
        "suggest", help="Suggest optimal reaction conditions based on genome and primer length"
    )
    suggest_parser.add_argument(
        "--genome-gc", type=float, help="Target genome GC content (0-1, e.g., 0.65 for 65%%)"
    )
    suggest_parser.add_argument(
        "--genome", type=str, help="Target genome FASTA file (calculates GC automatically)"
    )
    suggest_parser.add_argument(
        "--primer-length", "-l", type=int, help="Target primer length in bp"
    )
    suggest_parser.add_argument(
        "--context",
        "-c",
        choices=["standard", "clinical", "high_throughput", "low_input"],
        default="standard",
        help="Application context (default: standard)",
    )
    suggest_parser.add_argument(
        "--polymerase",
        choices=["phi29", "equiphi29", "bst", "klenow"],
        default="phi29",
        help="Polymerase type (default: phi29)",
    )
    suggest_parser.add_argument(
        "--optimize-for",
        choices=["amplification", "specificity", "coverage", "processivity"],
        default="amplification",
        help="Optimization goal (default: amplification)",
    )
    suggest_parser.add_argument(
        "--use-optimizer",
        action="store_true",
        help="Use advanced multi-additive optimizer (grid search)",
    )
    suggest_parser.add_argument(
        "--sweep",
        action="store_true",
        help="Sweep a grid of reaction conditions and rank by " "predicted amplification factor",
    )
    suggest_parser.add_argument(
        "--kmer-range",
        action="store_true",
        help="Print a recommended (min_k, max_k) range based on genome size, "
        "GC, and polymerase, then exit",
    )
    suggest_parser.add_argument(
        "--genome-size",
        type=int,
        help="Target genome size in bp (for --kmer-range). Auto-derived from --genome if that flag is given.",
    )
    suggest_parser.add_argument(
        "--output", "-o", type=str, help="Output CSV path for condition sweep results"
    )

    # =========================================================================
    # ADVANCED: Optimize conditions
    # =========================================================================
    optimize_cond_parser = subparsers.add_parser(
        "optimize-conditions", help="Find optimal reaction conditions for a genome"
    )
    optimize_cond_parser.add_argument(
        "--fg", "--foreground", required=True, help="Foreground genome FASTA file"
    )
    optimize_cond_parser.add_argument("--output", "-o", required=True, help="Output directory")
    optimize_cond_parser.add_argument(
        "--target-k", type=int, help="Target k-mer length to optimize for"
    )

    # =========================================================================
    # ADVANCED: Analyze primer set
    # =========================================================================
    analyze_parser = subparsers.add_parser(
        "analyze-set", help="[EXPERIMENTAL] Analyze an existing primer set"
    )
    analyze_parser.add_argument(
        "--primers", required=True, nargs="+", help="Primer sequences to analyze"
    )
    analyze_parser.add_argument("--fg", required=True, help="Foreground genome FASTA")
    analyze_parser.add_argument("--fg-kmers", required=True, help="Foreground k-mer file prefix")
    analyze_parser.add_argument("--output", "-o", required=True, help="Output directory")
    analyze_parser.add_argument(
        "--preset", default="standard_phi29", help="Reaction conditions preset"
    )
    analyze_parser.add_argument(
        "--simulate", action="store_true", help="Run replication simulation"
    )

    # =========================================================================
    # CATEGORY 1: Genome analysis (orphaned feature)
    # =========================================================================
    analyze_genome_parser = subparsers.add_parser(
        "analyze-genome", help="[EXPERIMENTAL] Analyze genome suitability for SWGA"
    )
    analyze_genome_parser.add_argument(
        "--genome", required=True, help="Genome FASTA file to analyze"
    )
    analyze_genome_parser.add_argument(
        "--output", "-o", required=True, help="Output directory for analysis report"
    )
    analyze_genome_parser.add_argument(
        "--window-size",
        type=int,
        default=1000,
        help="Window size for GC profiling (default: 1000bp)",
    )

    # =========================================================================
    # CATEGORY 1: Dimer network analysis (orphaned feature)
    # =========================================================================
    analyze_dimers_parser = subparsers.add_parser(
        "analyze-dimers", help="[EXPERIMENTAL] Analyze primer dimer interaction network"
    )
    analyze_dimers_parser.add_argument(
        "--primers", required=True, nargs="+", help="Primer sequences to analyze"
    )
    analyze_dimers_parser.add_argument(
        "--output", "-o", required=True, help="Output directory for network visualization"
    )
    analyze_dimers_parser.add_argument(
        "--threshold", type=float, default=0.3, help="Dimer severity threshold (default: 0.3)"
    )
    analyze_dimers_parser.add_argument(
        "--visualize", action="store_true", help="Generate network visualization"
    )

    # =========================================================================
    # CATEGORY 1: 3' stability analysis (orphaned feature)
    # =========================================================================
    analyze_stability_parser = subparsers.add_parser(
        "analyze-stability", help="Analyze 3' end stability and specificity"
    )
    analyze_stability_parser.add_argument(
        "--primers", required=True, nargs="+", help="Primer sequences to analyze"
    )
    analyze_stability_parser.add_argument(
        "--output", "-o", required=True, help="Output file for stability report"
    )
    analyze_stability_parser.add_argument(
        "--temp", type=float, default=37.0, help="Reaction temperature in C (default: 37)"
    )

    # =========================================================================
    # CATEGORY 1: Optimal oligo generator (orphaned feature)
    # =========================================================================

    # =========================================================================
    # CATEGORY 3: Multi-genome pipeline (orphaned feature)
    # =========================================================================
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
        choices=["phi29", "equiphi29", "bst", "klenow"],
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
        choices=["phi29", "equiphi29", "bst", "klenow"],
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
    expand_parser = subparsers.add_parser(
        "expand-primers",
        help="Expand existing primer set with additional primers to fill coverage gaps",
    )
    expand_parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    expand_parser.add_argument(
        "--fixed-primers",
        nargs="+",
        required=True,
        help="Primer sequences to keep (already validated)",
    )
    expand_parser.add_argument(
        "--fixed-primers-file", help="File with fixed primers (one per line)"
    )
    expand_parser.add_argument(
        "--failed-primers", nargs="+", help="Primer sequences to exclude (failed in wet lab)"
    )
    expand_parser.add_argument(
        "--failed-primers-file", help="File with failed primers to exclude (one per line)"
    )
    expand_parser.add_argument(
        "--num-new", type=int, default=6, help="Number of new primers to add (default: 6)"
    )
    expand_parser.add_argument(
        "--optimization-method",
        default="hybrid",
        choices=["hybrid", "dominating-set", "network", "background-aware", "ensemble"],
        help="Optimization method (default: hybrid)",
    )
    expand_parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible expansion.",
    )
    expand_parser.add_argument(
        "--bam",
        help="Mapped BAM of real sequencing reads against the "
        "target genome. Low-depth regions are merged with "
        "in-silico coverage gaps and the candidate pool is "
        "focused on primers that bind inside the gaps. "
        "Requires the [bam] extra (pip install 'neoswga[bam]').",
    )
    expand_parser.add_argument(
        "--min-depth",
        type=int,
        default=5,
        help="Sequencing depth below which a base counts as a " "BAM coverage gap (default: 5).",
    )
    expand_parser.add_argument(
        "--min-gap-size",
        type=int,
        default=10000,
        help="Minimum gap length to act on, bp (default: 10000).",
    )
    expand_parser.add_argument(
        "--contig-alias",
        action="append",
        default=None,
        metavar="FG=BAMCONTIG",
        help="Map a foreground prefix/basename to a BAM contig "
        "name when they differ. Repeatable.",
    )
    expand_parser.add_argument(
        "--output", "-o", required=True, help="Output directory for expanded primer set"
    )
    expand_parser.add_argument(
        "--quiet", "-q", action="store_true", help="Suppress progress output"
    )

    # Analyze coverage - read-only inspection of in-silico + BAM gaps
    cov_parser = subparsers.add_parser(
        "analyze-coverage",
        help="Report coverage gaps for a primer set, from in-silico binding "
        "sites and (optionally) real sequencing depth (BAM). Read-only.",
    )
    cov_parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    cov_parser.add_argument(
        "--primers", nargs="+", help="Current primer set (the primers to assess)."
    )
    cov_parser.add_argument("--primers-file", help="File with current primers (one per line).")
    cov_parser.add_argument(
        "--bam", help="Mapped BAM of real reads vs the target genome. " "Requires the [bam] extra."
    )
    cov_parser.add_argument(
        "--min-depth",
        type=int,
        default=5,
        help="Depth below which a base is a BAM gap (default: 5).",
    )
    cov_parser.add_argument(
        "--min-gap-size",
        type=int,
        default=10000,
        help="Minimum gap length to report, bp (default: 10000).",
    )
    cov_parser.add_argument(
        "--contig-alias",
        action="append",
        default=None,
        metavar="FG=BAMCONTIG",
        help="Map a foreground prefix/basename to a BAM contig. " "Repeatable.",
    )
    cov_parser.add_argument(
        "--output", "-o", required=True, help="Output directory for gap BED/JSON."
    )
    cov_parser.add_argument("--quiet", "-q", action="store_true", help="Suppress progress output")

    # Swap primer: replace worst performers with better candidates
    swap_parser = subparsers.add_parser(
        "swap-primer",
        help="Replace under-performing primers in an existing set with better "
        "candidates from a pool (greedy dimer-aware optimization).",
    )
    swap_parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    swap_parser.add_argument("--primers", nargs="+", required=True, help="Current primer set")
    swap_parser.add_argument("--primers-file", help="File with current primers (one per line)")
    swap_parser.add_argument(
        "--candidates-file",
        help="File with candidate pool for replacement, "
        "one primer per line. Defaults to step2_df.csv "
        "primers if present in data_dir.",
    )
    swap_parser.add_argument(
        "--max-swaps", type=int, default=3, help="Maximum replacements to attempt (default: 3)"
    )
    swap_parser.add_argument(
        "--output",
        "-o",
        help="Output JSON with before/after set and per-swap " "deltas (default: stdout)",
    )
    swap_parser.add_argument("--quiet", "-q", action="store_true")
    swap_parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducible swaps."
    )

    # Contract set: remove redundant primers while keeping coverage above threshold
    contract_parser = subparsers.add_parser(
        "contract-set",
        help="Remove redundant primers from a set while keeping coverage above "
        "--min-coverage. Greedy leave-one-out.",
    )
    contract_parser.add_argument("-j", "--json-file", required=True)
    contract_parser.add_argument("--primers", nargs="+", required=True)
    contract_parser.add_argument("--primers-file")
    contract_parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.70,
        help="Minimum foreground coverage to keep (default: 0.70)",
    )
    contract_parser.add_argument("--output", "-o")
    contract_parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducible contraction."
    )
    contract_parser.add_argument("--quiet", "-q", action="store_true")

    # Rescore set under arbitrary conditions
    rescore_parser = subparsers.add_parser(
        "rescore-set",
        help="Rescore an existing primer set under arbitrary reaction "
        "conditions (polymerase + additives) without re-optimizing.",
    )
    rescore_parser.add_argument("-j", "--json-file", required=True)
    rescore_parser.add_argument("--primers", nargs="+", required=True)
    rescore_parser.add_argument("--primers-file")
    rescore_parser.add_argument("--polymerase", choices=["phi29", "equiphi29", "bst", "klenow"])
    rescore_parser.add_argument("--reaction-temp", type=float)
    rescore_parser.add_argument("--mg-conc", type=float)
    rescore_parser.add_argument("--na-conc", type=float)
    rescore_parser.add_argument("--dmso-percent", type=float)
    rescore_parser.add_argument("--betaine-m", type=float)
    rescore_parser.add_argument("--trehalose-m", type=float)
    rescore_parser.add_argument("--formamide-percent", type=float)
    rescore_parser.add_argument("--ethanol-percent", type=float)
    rescore_parser.add_argument("--urea-m", type=float)
    rescore_parser.add_argument("--tmac-m", type=float)
    rescore_parser.add_argument("--output", "-o")
    rescore_parser.add_argument("--quiet", "-q", action="store_true")

    # Predict efficiency - unified confidence score for primer set
    predict_parser = subparsers.add_parser(
        "predict-efficiency", help="Predict efficiency of primer set before synthesis"
    )
    predict_parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    predict_parser.add_argument(
        "--primers", nargs="+", required=True, help="Primer sequences to evaluate"
    )
    predict_parser.add_argument("--primers-file", help="File with primers (one per line)")
    predict_parser.add_argument(
        "--run-simulation",
        action="store_true",
        help="Run simulation for more accurate prediction (slower)",
    )
    predict_parser.add_argument(
        "--track", action="store_true", help="Record prediction for later outcome tracking"
    )
    predict_parser.add_argument("--output", "-o", help="Output JSON file for prediction results")
    predict_parser.add_argument(
        "--quiet", "-q", action="store_true", help="Suppress progress output"
    )

    # Background registry - list available pre-computed backgrounds
    bg_list_parser = subparsers.add_parser(
        "background-list", help="List available pre-computed background genomes"
    )
    bg_list_parser.add_argument(
        "--discover", action="store_true", help="Search for new backgrounds in default directories"
    )
    bg_list_parser.add_argument("--search", help="Search for backgrounds by name or species")
    bg_list_parser.add_argument(
        "--quiet", "-q", action="store_true", help="Suppress verbose output"
    )

    # Background registry - add new background
    bg_add_parser = subparsers.add_parser(
        "background-add", help="Add a pre-computed background to the registry"
    )
    bg_add_parser.add_argument(
        "--name", required=True, help='Human-readable name (e.g., "Human GRCh38")'
    )
    bg_add_parser.add_argument(
        "--species", required=True, help='Species name (e.g., "Homo sapiens")'
    )
    bg_add_parser.add_argument("--bloom-path", help="Path to Bloom filter pickle file")
    bg_add_parser.add_argument(
        "--kmer-prefix", help="Prefix for k-mer files (without _Xmer_all.txt)"
    )
    bg_add_parser.add_argument("--genome-size", type=int, default=0, help="Genome size in bp")
    bg_add_parser.add_argument(
        "--min-k", type=int, default=6, help="Minimum k-mer length (default: 6)"
    )
    bg_add_parser.add_argument(
        "--max-k", type=int, default=12, help="Maximum k-mer length (default: 12)"
    )
    bg_add_parser.add_argument("--description", help="Additional description")
    bg_add_parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing entry with same name"
    )

    # =========================================================================
    # UTILITY: Genome library management
    # =========================================================================
    genome_add_parser = subparsers.add_parser(
        "genome-add", help="Add a genome to the pre-calculated library"
    )
    genome_add_parser.add_argument("name", help="Identifier for the genome (e.g., human-grch38)")
    genome_add_parser.add_argument("fasta", help="Path to genome FASTA file")
    genome_add_parser.add_argument(
        "--role",
        choices=["bg", "bl"],
        default="bg",
        help="Role: bg (background) or bl (blacklist). Default: bg",
    )
    genome_add_parser.add_argument("--species", default="", help="Species name (optional)")
    genome_add_parser.add_argument(
        "--k-ranges", default="6-12,12-18", help="K-mer ranges to compute (default: 6-12,12-18)"
    )
    genome_add_parser.add_argument(
        "--no-bloom", action="store_true", help="Skip Bloom filter construction"
    )

    genome_list_parser = subparsers.add_parser(
        "genome-list", help="List genomes in the pre-calculated library"
    )
    genome_list_parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show detailed information"
    )

    genome_remove_parser = subparsers.add_parser(
        "genome-remove", help="Remove a genome from the pre-calculated library"
    )
    genome_remove_parser.add_argument("name", help="Identifier of genome to remove")

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
