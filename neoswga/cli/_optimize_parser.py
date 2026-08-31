"""Argument groups for the `optimize` subcommand.

Split out of `cli/pipeline.py`: three independent fixes landed in that module
at once and pushed it past its size budget. These builders touch nothing but
the parser they are handed, so they carry no coupling back.
"""


def _add_optimize_option_groups(parser):
    """Background-filtering and performance groups for `optimize`."""
    # Background Filtering
    opt_bg_group = parser.add_argument_group("Background Filtering")
    opt_bg_group.add_argument(
        "--use-background-filter",
        action="store_true",
        default=False,
        help="Not implemented: optimize builds no Bloom filter. The "
        "background screening it does perform is controlled by "
        "--no-bg-prefilter",
    )
    opt_bg_group.add_argument(
        "--no-bg-prefilter",
        action="store_true",
        default=False,
        help="Disable automatic background pre-filtering of candidates "
        "(enabled by default when background genome data is available)",
    )
    opt_bg_group.add_argument(
        "--background-bloom-path",
        type=str,
        help="Path to a pre-built Bloom filter. Not implemented: optimize "
        "reads no pre-built background index",
    )
    opt_bg_group.add_argument(
        "--background-sampled-path",
        type=str,
        help="Path to a pre-built sampled index. Not implemented: optimize "
        "reads no pre-built background index",
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
    opt_perf_group = parser.add_argument_group("Performance")
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
        default=None,
        help="Maximum optimization time in seconds. Not implemented: no "
        "optimizer enforces a wall-clock budget, so there is no default to "
        "state",
    )
    opt_perf_group.add_argument(
        "--max-extension",
        type=int,
        default=None,
        help="Polymerase extension length in bp for the Stage-2 "
        "amplification network, which asks whether two binding sites can be "
        "connected by one molecule. Distinct from --coverage-reach, the much "
        "shorter per-primer reach used to select for coverage. Unset takes "
        "the polymerase preset (phi29 70000, bst 2000)",
    )
    opt_perf_group.add_argument(
        "--coverage-reach",
        type=int,
        default=None,
        help="Per-primer extension reach in bp for coverage and set-cover "
        "selection (default: the polymerase's realistic reach, phi29 ~3000). "
        "This is NOT --max-extension, which is the amplification-network "
        "reach. Coverage figures are not comparable across different reaches; "
        "estimate this from data with 'neoswga calibrate-reach --bam'.",
    )

    _add_optimize_selection_groups(parser)


def _add_optimize_selection_groups(parser):
    """Set-size, mechanistic-model, post-processing and validation groups."""
    # Set Size & Application
    opt_size_group = parser.add_argument_group("Set Size & Application")
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
    opt_mech_group = parser.add_argument_group("Mechanistic Model")
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
        default=None,
        help="Weight for coverage uniformity in scoring (0.0-1.0). Higher "
        "values prioritize even genome coverage over raw enrichment. Left "
        "unset, the --application profile supplies it; 0.0 disables the term.",
    )

    # Post-processing
    opt_post_group = parser.add_argument_group("Post-processing")
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
    opt_val_group = parser.add_argument_group("Validation")
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
