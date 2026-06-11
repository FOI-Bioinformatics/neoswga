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


def run_validate(args):
    """Dispatcher for `neoswga validate [install|params|model]`.

    - No subcommand or `install`: installation / environment check
      (historical behavior of `neoswga validate`).
    - `params`: delegate to :func:`run_validate_params`.
    - `model`: delegate to :func:`run_validate_model`.
    """
    mode = getattr(args, "validate_mode", None)
    if mode == "params":
        return run_validate_params(args)
    if mode == "model":
        return run_validate_model(args)
    # mode is None (bare `validate`) or 'install': installation check.
    return _run_validate_installation(args)


def _run_validate_installation(args):
    """Run installation / environment validation tests."""
    from neoswga.core.gpu_acceleration import get_gpu_info, is_gpu_available
    from neoswga.core.kmer_counter import check_jellyfish_available
    from neoswga.core.validation import ValidationSuite, quick_validation

    # Show system capabilities
    logger.info("System Capabilities:")
    logger.info(f"  Jellyfish: {'Available' if check_jellyfish_available() else 'NOT FOUND'}")
    if is_gpu_available():
        gpu_info = get_gpu_info()
        logger.info(
            f"  GPU: {gpu_info.get('device_name', 'Available')} ({gpu_info.get('memory_total_gb', 0):.1f} GB)"
        )
    else:
        logger.info("  GPU: Not available (install CuPy for GPU acceleration)")

    if args.quick:
        logger.info("Running quick validation...")
        success = quick_validation()
    elif args.all:
        logger.info("Running full validation + benchmarks...")
        suite = ValidationSuite(verbose=True)
        success = suite.run_all_tests()

        # Also run benchmarks if available
        try:
            import importlib.util

            bench_path = os.path.join(
                os.path.dirname(__file__),
                "..",
                "scripts",
                "benchmarking",
                "benchmark_improvements.py",
            )
            spec = importlib.util.spec_from_file_location("benchmark_improvements", bench_path)
            if spec and spec.loader:
                logger.info("\nRunning benchmarks...")
                benchmark_improvements = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(benchmark_improvements)
                benchmark_improvements.main()
            else:
                logger.info("\nBenchmark script not found, skipping.")
        except Exception as e:
            logger.warning(f"Benchmarks unavailable: {e}")
    else:
        logger.info("Running standard validation...")
        suite = ValidationSuite(verbose=True)
        success = suite.run_all_tests()

    sys.exit(0 if success else 1)


def run_doctor(args):
    """Diagnostic health check.

    Reports version, jellyfish location, registered optimizers with
    additive-awareness, params.json summary, and file-existence status.
    Exits non-zero if any blocking issue is found.
    """
    import json as _json
    import os as _os
    import platform as _platform
    import shutil as _shutil
    import subprocess as _subprocess
    import sys as _sys

    blocking: list = []
    warnings: list = []
    summary: dict = {
        "platform": {
            "system": _platform.system(),
            "python": _sys.version.split()[0],
            "neoswga": _import_neoswga_version(),
        },
        "tools": {},
        "optimizers": [],
        "params": None,
    }

    # Jellyfish
    jelly = _shutil.which("jellyfish")
    if jelly:
        try:
            version = _subprocess.run(
                [jelly, "--version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            first = (version.stdout or version.stderr or "").splitlines()
            summary["tools"]["jellyfish"] = {
                "path": jelly,
                "version": first[0] if first else "unknown",
            }
        except Exception as e:
            summary["tools"]["jellyfish"] = {"path": jelly, "error": str(e)}
            warnings.append(f"jellyfish --version failed: {e}")
    else:
        blocking.append(
            "jellyfish is not on PATH. Install with: conda install -c bioconda jellyfish"
        )

    # Optimizer capability matrix. Each entry declares whether the optimizer
    # is registered with the factory and whether its scoring path meaningfully
    # references `self.conditions`. The awareness flag is computed by grepping
    # the source file for ReactionConditions usage plus `self.conditions`.
    try:
        from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry

        # Force registration of built-in optimizers
        try:
            from neoswga.core import unified_optimizer as _uo

            _uo._ensure_optimizers_registered()
        except Exception:
            pass
        registered = OptimizerFactory.list_optimizers()
        for name in sorted(registered.keys()):
            cls = OptimizerRegistry.get(name)
            # Phase 15D: structural class-level flag replaces the old
            # source-grep heuristic. Subclasses that actually use
            # ReactionConditions in selection opt-in by setting
            # ADDITIVE_AWARE = True on the class.
            aware = bool(getattr(cls, "ADDITIVE_AWARE", False))
            summary["optimizers"].append(
                {
                    "name": name,
                    "class": f"{cls.__module__}.{cls.__name__}",
                    "additive_aware": aware,
                    "description": registered[name],
                }
            )
    except Exception as e:
        warnings.append(f"Optimizer discovery failed: {e}")

    # params.json summary (optional)
    if args.json_file:
        if not _os.path.isfile(args.json_file):
            blocking.append(f"params.json not found at {args.json_file}")
        else:
            try:
                with open(args.json_file) as fh:
                    raw = _json.load(fh)
                params_summary = {
                    "path": args.json_file,
                    "schema_version": raw.get("schema_version"),
                    "polymerase": raw.get("polymerase"),
                    "reaction_temp": raw.get("reaction_temp"),
                    "mg_conc": raw.get("mg_conc"),
                    "additives": {
                        k: raw.get(k, 0.0)
                        for k in (
                            "dmso_percent",
                            "betaine_m",
                            "trehalose_m",
                            "formamide_percent",
                            "ethanol_percent",
                            "urea_m",
                            "tmac_m",
                        )
                    },
                    "min_k": raw.get("min_k"),
                    "max_k": raw.get("max_k"),
                    "num_targets": len(raw.get("fg_genomes", []) or []),
                    "num_backgrounds": len(raw.get("bg_genomes", []) or []),
                    "num_blacklists": len(raw.get("bl_genomes", []) or []),
                    "adaptive_gc": raw.get("adaptive_gc", True),
                    "genome_gc": raw.get("genome_gc"),
                }
                # File existence checks for each prefix
                missing_files: list = []
                for prefix in raw.get("fg_prefixes", []) or []:
                    mn = raw.get("min_k", 6)
                    mx = raw.get("max_k", 12)
                    for k in range(mn, mx + 1):
                        candidate = f"{prefix}_{k}mer_all.txt"
                        if not _os.path.isfile(candidate):
                            missing_files.append(candidate)
                            if len(missing_files) > 10:
                                break
                    if len(missing_files) > 10:
                        break
                if missing_files:
                    warnings.append(
                        "Missing k-mer count files (run count-kmers?): "
                        + ", ".join(missing_files[:5])
                        + ("…" if len(missing_files) > 5 else "")
                    )
                    params_summary["missing_kmer_files"] = missing_files[:20]
                summary["params"] = params_summary
            except Exception as e:
                blocking.append(f"Could not parse params.json: {e}")

    summary["ok"] = not blocking
    summary["blocking"] = blocking
    summary["warnings"] = warnings

    if args.output_json:
        print(_json.dumps(summary, indent=2))
    else:
        # Human-readable report
        print("=" * 60)
        print("neoswga doctor")
        print("=" * 60)
        print(f"  platform: {summary['platform']['system']}")
        print(f"  python: {summary['platform']['python']}")
        print(f"  neoswga: {summary['platform']['neoswga']}")
        j = summary["tools"].get("jellyfish")
        if j:
            print(f"  jellyfish: {j.get('version', j.get('error', '?'))} at {j.get('path', '')}")
        else:
            print("  jellyfish: NOT FOUND")
        print()
        print(f"Registered optimizers: {len(summary['optimizers'])}")
        aware_count = sum(1 for o in summary["optimizers"] if o["additive_aware"])
        print(f"  additive-aware: {aware_count}")
        print(f"  additive-blind: {len(summary['optimizers']) - aware_count}")
        for opt in summary["optimizers"]:
            marker = "+" if opt["additive_aware"] else "-"
            print(f"  [{marker}] {opt['name']}")
        print()
        if summary["params"]:
            ps = summary["params"]
            print("params.json summary:")
            print(f"  path: {ps['path']}")
            print(f"  schema_version: {ps['schema_version']}")
            print(f"  polymerase: {ps['polymerase']}  reaction_temp: {ps['reaction_temp']}")
            print(f"  min_k: {ps['min_k']}  max_k: {ps['max_k']}")
            print(
                f"  targets: {ps['num_targets']}  backgrounds: {ps['num_backgrounds']}  blacklists: {ps['num_blacklists']}"
            )
            nz_add = {k: v for k, v in ps["additives"].items() if v}
            if nz_add:
                print(f"  additives: {nz_add}")
            print()
        if warnings:
            print("Warnings:")
            for w in warnings:
                print(f"  ! {w}")
            print()
        if blocking:
            print("BLOCKING:")
            for b in blocking:
                print(f"  X {b}")
            print()
        if summary["ok"]:
            print("Status: READY")
        else:
            print("Status: BLOCKED — fix the items above before running the pipeline")

    sys.exit(0 if summary["ok"] else 1)


def _import_neoswga_version() -> str:
    try:
        import neoswga

        return getattr(neoswga, "__version__", "unknown")
    except Exception:
        return "unknown"


def run_init(args):
    """Run the setup wizard to create params.json"""
    import inspect

    from neoswga.core.wizard import run_wizard

    try:
        # Build kwargs, only pass blacklist_paths if wizard supports it
        wizard_kwargs = dict(
            genome_path=args.genome,
            background_path=args.background,
            output_path=args.output,
            output_dir=args.output_dir,
            interactive=not args.non_interactive,
            advanced=args.advanced,
            auto_approve=args.yes,
        )
        if getattr(args, "blacklist", None):
            sig = inspect.signature(run_wizard)
            if "blacklist_paths" in sig.parameters:
                wizard_kwargs["blacklist_paths"] = args.blacklist
            else:
                logger.warning("Blacklist support for init requires wizard update")
        run_wizard(**wizard_kwargs)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)
    except KeyboardInterrupt:
        logger.info("\nSetup cancelled")
        sys.exit(1)


def run_validate_params(args):
    """Validate params.json configuration.

    Invoked as `neoswga validate params -j params.json` (preferred) or
    `neoswga validate-params -j params.json` (deprecated alias).
    """
    from neoswga.core.param_validator import validate_params_file

    if getattr(args, "command", None) == "validate-params":
        logger.warning(
            "`neoswga validate-params` is deprecated; " "use `neoswga validate params` instead."
        )

    success, _ = validate_params_file(args.json_file, verbose=True)
    sys.exit(0 if success else 1)


def run_schema(args):
    """Dump the canonical params.json JSON Schema.

    Either prints to stdout (default with --dump) or writes to --output path.
    Useful for IDE integration (VS Code can autocomplete against the schema)
    and for human-readable review.
    """
    import json

    try:
        from neoswga.core.schema import load_schema

        schema = load_schema()
    except Exception as e:
        logger.error(f"Could not load params schema: {e}")
        sys.exit(1)

    output = json.dumps(schema, indent=2)
    if getattr(args, "output", None):
        with open(args.output, "w", encoding="utf-8") as fh:
            fh.write(output + "\n")
        logger.info(f"Schema written to {args.output}")
    else:
        print(output)


def run_validate_model(args):
    """Validate mechanistic model against expected behavior.

    Invoked as `neoswga validate model` (preferred) or
    `neoswga validate-model` (deprecated alias).
    """
    import json

    from neoswga.core.model_validation import (
        format_validation_report,
        validate_mechanistic_model,
    )

    if getattr(args, "command", None) == "validate-model":
        logger.warning(
            "`neoswga validate-model` is deprecated; " "use `neoswga validate model` instead."
        )

    logger.info("Validating mechanistic model...")
    results = validate_mechanistic_model()

    if getattr(args, "output_json", False):
        # Output as JSON
        # Convert results to JSON-serializable format
        json_results = []
        for r in results:
            jr = {
                "test": r.get("test", "Unknown"),
                "passed": r.get("passed", False),
                "summary": r.get("summary", ""),
            }
            if "error" in r:
                jr["error"] = r["error"]
            json_results.append(jr)
        print(json.dumps(json_results, indent=2))
    else:
        # Output as formatted report
        report = format_validation_report(results)
        print(report)

    # Exit with appropriate code
    all_passed = all(r.get("passed", False) for r in results)
    sys.exit(0 if all_passed else 1)


def run_interpret(args):
    """Interpret pipeline results"""
    from neoswga.core.results_interpreter import interpret_results

    try:
        interpret_results(args.dir, verbose=True)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)


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


def run_start(args):
    """Run interactive workflow selector"""
    from neoswga.core.workflow_selector import run_workflow_selector

    run_workflow_selector()


def run_suggest(args):
    """Suggest optimal reaction conditions"""
    genome_gc = args.genome_gc

    # Calculate GC from genome file if provided
    if args.genome and not genome_gc:
        from neoswga.core.genome_analysis import calculate_genome_gc

        try:
            genome_gc = calculate_genome_gc(args.genome)
            logger.info(f"Calculated GC content: {genome_gc:.1%}")
        except FileNotFoundError:
            logger.error(f"Genome file not found: {args.genome}")
            sys.exit(1)

    # K-mer range suggestion: emit a recommended (min_k, max_k) based on
    # genome size + GC + polymerase, then exit. Complements the warning in
    # get_params() by giving the user a direct answer.
    if getattr(args, "kmer_range", False):
        from neoswga.core.condition_suggester import suggest_kmer_range
        from neoswga.core.genome_io import GenomeLoader

        genome_size = getattr(args, "genome_size", None)
        if genome_size is None and getattr(args, "genome", None):
            try:
                sequence = GenomeLoader().load_genome(args.genome, return_stats=False)
                genome_size = len(sequence)
            except Exception as e:
                logger.error(f"Could not read genome size from {args.genome}: {e}")
                sys.exit(1)
        if genome_size is None:
            logger.error("--kmer-range requires --genome-size or --genome")
            sys.exit(1)

        polymerase = getattr(args, "polymerase", "phi29")
        min_k, max_k = suggest_kmer_range(
            genome_size,
            gc=genome_gc,
            polymerase=polymerase,
        )
        print(
            f"Recommended k-mer range for {genome_size:,} bp, "
            f"GC={genome_gc if genome_gc is not None else 'unknown'}, "
            f"polymerase={polymerase}:"
        )
        print(f"  min_k = {min_k}")
        print(f"  max_k = {max_k}")
        print(f"\nAdd to params.json:")
        print(f'  "min_k": {min_k},')
        print(f'  "max_k": {max_k}')
        return

    # Condition sweep mode
    if getattr(args, "sweep", False):
        from neoswga.core.condition_suggester import sweep_conditions

        primer_length = args.primer_length
        if primer_length is None:
            # Default primer length by polymerase type
            polymerase = getattr(args, "polymerase", "phi29")
            primer_length_defaults = {
                "phi29": 9,
                "equiphi29": 15,
                "bst": 20,
                "klenow": 12,
            }
            primer_length = primer_length_defaults.get(polymerase, 10)
        if genome_gc is None:
            genome_gc = 0.5
            logger.info("No GC content provided, using default 50%")
        sweep_conditions(
            genome_gc=genome_gc,
            primer_length=primer_length,
            polymerase=getattr(args, "polymerase", "phi29"),
            output_path=getattr(args, "output", None),
            verbose=True,
        )
        return

    if genome_gc is None and args.primer_length is None:
        logger.error("Provide at least --genome-gc or --genome or --primer-length")
        sys.exit(1)

    # Use advanced optimizer if requested or if polymerase/optimize-for specified
    use_optimizer = getattr(args, "use_optimizer", False)
    polymerase = getattr(args, "polymerase", "phi29")
    optimize_for = getattr(args, "optimize_for", "amplification")

    if use_optimizer or polymerase != "phi29" or optimize_for != "amplification":
        # Use new AdditiveOptimizer
        from neoswga.core.additive_optimizer import AdditiveOptimizer
        from neoswga.core.mechanistic_params import get_polymerase_params

        # Determine primer length if not provided
        primer_length = args.primer_length
        if primer_length is None:
            poly_params = get_polymerase_params(polymerase)
            primer_length = int(
                (poly_params["primer_length_range"][0] + poly_params["primer_length_range"][1]) / 2
            )
            logger.info(f"Using default primer length for {polymerase}: {primer_length}bp")

        # Default GC if not provided
        if genome_gc is None:
            genome_gc = 0.5
            logger.info("No GC content provided, using default 50%")

        optimizer = AdditiveOptimizer(polymerase=polymerase)
        recommendation = optimizer.optimize(
            primer_length=primer_length,
            template_gc=genome_gc,
            optimize_for=optimize_for,
        )

        # Print detailed recommendation
        print(recommendation.summary())

    else:
        # Use original suggest_conditions for backward compatibility
        from neoswga.core.condition_suggester import suggest_conditions

        suggest_conditions(
            genome_gc=genome_gc,
            primer_length=args.primer_length,
            context=args.context,
            verbose=True,
        )


# Analysis handlers live in neoswga.cli.analysis (re-exported for dispatch).
from neoswga.cli.analysis import (  # noqa: E402,F401
    analyze_primer_set,
    optimize_conditions,
    run_analyze_dimers,
    run_analyze_genome,
    run_analyze_stability,
    show_presets,
)

# Simulation handlers live in neoswga.cli.simulate (re-exported for dispatch).
from neoswga.cli.simulate import run_multi_genome, run_simulate  # noqa: E402,F401

# =========================================================================
# CATEGORY 3: Handler functions for pipeline orphaned features
# =========================================================================


def run_analyze_coverage(args):
    """Read-only coverage-gap report from in-silico binding sites + optional BAM.

    Computes in-silico gaps for the given primer set, optionally merges low
    sequencing-depth regions from a mapped BAM, and writes the merged gaps as
    BED + JSON (plus an in-silico binding-density BedGraph). No optimizer runs.
    """
    import json as json_module

    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.export import export_gaps_to_bed
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.primer_expansion import PrimerExpander

    quiet = getattr(args, "quiet", False)

    primers = collect_primers_from_args(
        getattr(args, "primers", None),
        getattr(args, "primers_file", None),
        name="primer",
        allow_empty=True,
    )

    validate_params_json_file(args.json_file)
    merge_args_to_parameter(args, parameter, ["json_file"])
    core_pipeline._initialize()

    fg_prefixes = core_pipeline.fg_prefixes
    fg_seq_lengths = core_pipeline.fg_seq_lengths
    bg_prefixes = core_pipeline.bg_prefixes
    bg_seq_lengths = core_pipeline.bg_seq_lengths

    cache = PositionCache(fg_prefixes, list(set(primers)) or ["A" * 8])
    expander = PrimerExpander(
        position_cache=cache,
        fg_prefixes=fg_prefixes,
        fg_seq_lengths=fg_seq_lengths,
        bg_prefixes=bg_prefixes,
        bg_seq_lengths=bg_seq_lengths,
    )

    bam_gaps_list = None
    if getattr(args, "bam", None):
        from neoswga.core.bam_coverage import bam_gaps

        aliases = {}
        for item in getattr(args, "contig_alias", None) or []:
            if "=" in item:
                k, v = item.split("=", 1)
                aliases[k] = v
        fg_circular = bool(getattr(parameter, "fg_circular", False))
        try:
            bam_gaps_list = bam_gaps(
                args.bam,
                fg_prefixes,
                fg_seq_lengths,
                min_depth=args.min_depth,
                min_gap_size=args.min_gap_size,
                circular=fg_circular,
                contig_aliases=aliases or None,
            )
        except RuntimeError as e:
            logger.error(str(e))
            sys.exit(1)
        if not quiet:
            logger.info(f"BAM low-depth gaps: {len(bam_gaps_list)}")

    gaps = expander.identify_gaps(
        primers,
        min_gap_size=args.min_gap_size,
        extra_gaps=bam_gaps_list,
        merge=True,
    )

    os.makedirs(args.output, exist_ok=True)
    export_gaps_to_bed(gaps, os.path.join(args.output, "coverage_gaps.bed"))
    with open(os.path.join(args.output, "coverage_gaps.json"), "w") as f:
        json_module.dump(
            {
                "n_gaps": len(gaps),
                "min_depth": args.min_depth,
                "min_gap_size": args.min_gap_size,
                "used_bam": bool(getattr(args, "bam", None)),
                "gaps": [
                    {"chromosome": g.chromosome, "start": g.start, "end": g.end, "size": g.size}
                    for g in gaps
                ],
            },
            f,
            indent=2,
        )

    # In-silico binding-density BedGraph (per foreground prefix) so the primer
    # binding distribution can be inspected alongside the gaps in a browser.
    if primers:
        from neoswga.core.export import export_to_bedgraph

        for prefix, length in zip(fg_prefixes, fg_seq_lengths):
            positions = {}
            for p in primers:
                sites = []
                for pos in cache.get_positions(prefix, p, "forward"):
                    sites.append((int(pos), "forward"))
                for pos in cache.get_positions(prefix, p, "reverse"):
                    sites.append((int(pos), "reverse"))
                positions[p] = sites
            base = os.path.basename(prefix) or "genome"
            export_to_bedgraph(
                primers,
                positions,
                base,
                length,
                os.path.join(args.output, f"{base}_binding_density.bedgraph"),
            )

    if not quiet:
        logger.info(f"Found {len(gaps)} coverage gap(s); wrote BED + JSON to {args.output}")
        for g in gaps[:10]:
            logger.info(f"  {g.chromosome}: {g.start:,}-{g.end:,} ({g.size:,} bp)")


# Set-editing handlers live in neoswga.cli.iterate (re-exported for dispatch).
from neoswga.cli.iterate import (  # noqa: E402,F401
    run_contract_set,
    run_expand_primers,
    run_rescore_set,
    run_swap_primer,
)


def run_predict_efficiency(args):
    """
    Predict efficiency of primer set before synthesis.

    Provides unified confidence score with SYNTHESIZE/CAUTION/DO_NOT_SYNTHESIZE
    recommendation.
    """
    import json as json_module

    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.efficiency_predictor import EfficiencyPredictor
    from neoswga.core.experimental_tracker import ExperimentalTracker
    from neoswga.core.position_cache import PositionCache

    quiet = getattr(args, "quiet", False)

    if not quiet:
        logger.info("Primer Set Efficiency Prediction")
        logger.info("=" * 60)

    # Load and validate primers using consolidated helper
    primers = collect_primers_from_args(
        cli_primers=args.primers,
        primers_file=args.primers_file,
        name="primer",
        allow_empty=False,  # Must provide at least one primer
    )

    if not quiet:
        logger.info(f"Evaluating {len(primers)} primers")

    # Set json_file for pipeline initialization
    merge_args_to_parameter(args, parameter, ["json_file"])

    try:
        # Initialize pipeline to get genome info
        if not quiet:
            logger.info("Initializing pipeline...")
        core_pipeline._initialize()

        fg_prefixes = core_pipeline.fg_prefixes
        bg_prefixes = core_pipeline.bg_prefixes
        fg_seq_lengths = core_pipeline.fg_seq_lengths
        bg_seq_lengths = core_pipeline.bg_seq_lengths

        # Initialize position cache
        if not quiet:
            logger.info("Loading position data...")
        cache = PositionCache(fg_prefixes, primers)

        # Load genome sequence if simulation requested
        genome_sequence = None
        if args.run_simulation:
            fg_genomes = getattr(parameter, "fg_genomes", [])
            if fg_genomes:
                if not quiet:
                    logger.info("Loading genome for simulation...")
                from neoswga.core.genome_io import read_genome

                genome_sequence = read_genome(fg_genomes[0])

        # Create predictor
        predictor = EfficiencyPredictor(
            position_cache=cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
        )

        # Run prediction
        result = predictor.predict(
            primers=primers,
            run_simulation=args.run_simulation,
            genome_sequence=genome_sequence,
            verbose=not quiet,
        )

        # Track prediction if requested
        if args.track:
            tracker = ExperimentalTracker()
            exp_id = tracker.record_prediction(
                primer_set=primers,
                prediction=result,
            )
            if not quiet:
                logger.info(f"\nPrediction tracked: {exp_id}")
                logger.info("Record experimental outcomes manually for future reference.")

        # Save result if output specified
        if args.output:
            with open(args.output, "w") as f:
                json_module.dump(result.to_dict(), f, indent=2)
            if not quiet:
                logger.info(f"\nResults saved to: {args.output}")

        # Print summary
        if not quiet:
            logger.info("")
            logger.info("=" * 60)
            logger.info("PREDICTION RESULT")
            logger.info("=" * 60)
            logger.info(
                f"Confidence: {result.confidence_score:.2f} ({result.recommendation.value})"
            )
            logger.info(f"")
            logger.info(f"Predicted enrichment: {result.predicted_enrichment:.0f}x")
            logger.info(
                f"  95% CI: {result.enrichment_ci_low:.0f}x - {result.enrichment_ci_high:.0f}x"
            )
            logger.info(f"Predicted coverage: {result.predicted_coverage:.1%}")
            logger.info(f"")
            logger.info(f"Component scores:")
            logger.info(f"  ML score: {result.ml_score:.2f}")
            logger.info(f"  Network score: {result.network_score:.2f}")
            if result.simulation_score is not None:
                logger.info(f"  Simulation score: {result.simulation_score:.2f}")

            if result.limiting_factors:
                logger.info(f"")
                logger.info("Limiting factors:")
                for factor in result.limiting_factors:
                    logger.info(f"  - {factor}")

            if result.improvement_suggestions:
                logger.info(f"")
                logger.info("Suggestions:")
                for suggestion in result.improvement_suggestions:
                    logger.info(f"  - {suggestion}")

    except Exception as e:
        logger.error(f"Prediction failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


# Registry / library handlers live in neoswga.cli.registry (re-exported for
# the dispatch table and external importers).
from neoswga.cli.registry import (  # noqa: E402,F401
    run_background_add,
    run_background_list,
    run_genome_add,
    run_genome_list,
    run_genome_remove,
)

# =========================================================================
# UNIFIED: Complete pipeline handler
# =========================================================================


def run_design(args):
    """Run complete primer design pipeline (all 4 steps)"""
    check_jellyfish_available()
    logger.info("Running complete primer design pipeline")

    # Check for multi-genome mode
    if hasattr(args, "multi_genome") and args.multi_genome:
        logger.info(f"Multi-genome mode enabled for {len(args.multi_genome)} genomes")
        from neoswga.core import multi_genome_pipeline

        try:
            results = multi_genome_pipeline.design_pan_genome_primers(
                genome_paths=args.multi_genome,
                min_coverage=args.min_coverage,
                output_dir=getattr(args, "output", "./multi_genome_output"),
            )
            logger.info(f"Pan-genome design complete! Designed {len(results['primers'])} primers")
            return
        except Exception as e:
            logger.error(f"Multi-genome pipeline failed: {e}")
            import traceback

            traceback.print_exc()
            sys.exit(1)

    # Determine which steps to run
    start_step = getattr(args, "start_from", None) or 1
    stop_step = getattr(args, "stop_at", None) or 4

    if start_step > stop_step:
        logger.error("--start-from must be <= --stop-at")
        sys.exit(1)

    logger.info(f"Running steps {start_step} through {stop_step}")

    # Ensure optimize-specific attributes exist on args namespace
    # (the design subparser doesn't define these but run_step4 needs them)
    optimize_defaults = {
        "optimization_method": "hybrid",
        "num_primers": None,
        "iterations": None,
        "max_sets": None,
        "no_background": False,
        "use_mechanistic_model": False,
        "mechanistic_weight": 0.3,
        "auto_size": False,
        "application": "enrichment",
        "validate_with_simulation": False,
        "method_guide": False,
        "no_bg_prefilter": False,
        "strategy": None,
        "use_position_cache": True,
        "use_background_filter": False,
        "use_cooperative_binding": False,
        "primer_strategy": None,
        "enable_qa": False,
    }
    for attr, default in optimize_defaults.items():
        if not hasattr(args, attr):
            setattr(args, attr, default)

    # Run each step sequentially
    try:
        if start_step <= 1 <= stop_step:
            logger.info("=" * 80)
            logger.info("COUNT-KMERS: K-mer preprocessing")
            logger.info("=" * 80)
            run_step1(args)

        if start_step <= 2 <= stop_step:
            logger.info("=" * 80)
            logger.info("FILTER: Candidate primer filtering")
            logger.info("=" * 80)
            run_step2(args)

        if start_step <= 3 <= stop_step:
            logger.info("=" * 80)
            logger.info("SCORE: Amplification efficacy scoring")
            logger.info("=" * 80)
            run_step3(args)

        if start_step <= 4 <= stop_step:
            logger.info("=" * 80)
            logger.info("OPTIMIZE: Primer set optimization")
            logger.info("=" * 80)
            run_step4(args)

        logger.info("=" * 80)
        logger.info("Pipeline complete!")
        logger.info("=" * 80)

    except Exception as e:
        logger.error(f"Pipeline failed at step: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


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
