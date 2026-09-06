"""Remaining standalone CLI handlers: start (workflow selector), suggest
(reaction conditions), analyze-coverage, predict-efficiency, and design (the
unified count-kmers -> filter -> score -> optimize runner).

Extracted from cli_unified.py. run_design composes the pipeline-step handlers,
which now live in neoswga.cli.pipeline.
"""

import json
import logging
import os
import sys

from neoswga.cli._common import (
    check_jellyfish_available,
    collect_primers_from_args,
    merge_args_to_parameter,
    validate_params_json_file,
)
from neoswga.cli.pipeline import run_step1, run_step2, run_step3, run_step4
from neoswga.core.registry import polymerase_names as _polymerase_names
from neoswga.core.registry import views as _registry_views

logger = logging.getLogger(__name__)


def run_start(args):
    """Run interactive workflow selector"""
    from neoswga.core.workflow_selector import run_workflow_selector

    run_workflow_selector()


def _default_primer_length_for(polymerase: str) -> int:
    """Midpoint of the enzyme's primer-length range, from the registry.

    This read `get_polymerase_params(polymerase)["primer_length_range"]`, which
    returns the *mechanistic enzyme* parameters and has never carried that key,
    so `neoswga suggest --polymerase bst` exited 1 with "Missing required
    parameter". The range lives in the registry, which is where `parameter`,
    `param_validator` and the hybrid optimizer all read it.
    """
    from neoswga.core.registry import views as _views

    low, high = _views.primer_length_ranges().get((polymerase or "").lower(), (6, 12))
    return int((low + high) / 2)


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
            # Scalar suggested length - distinct from the (min_k, max_k) design
            # window in POLYMERASE_PRIMER_LENGTHS.
            primer_length_defaults = _registry_views.default_primer_lengths()
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

        # Determine primer length if not provided
        primer_length = args.primer_length
        if primer_length is None:
            primer_length = _default_primer_length_for(polymerase)
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

# Simulation handlers live in neoswga.cli.simulate (re-exported for dispatch).

# =========================================================================
# CATEGORY 3: Handler functions for pipeline orphaned features
# =========================================================================


def run_calibrate_reach(args):
    """Estimate the per-primer coverage reach from real sequencing depth.

    The reach decides both what coverage is reported and which primers get
    selected, and neither published convention is measured: swga 2.0 uses
    phi29's ~70 kb single-molecule processivity, and NeoSWGA's ~3 kb default is
    derived from inter-primer spacings that designers chose. A BAM from an SWGA
    run using these primers contains the answer, because depth falls away from
    binding sites on the length scale of the reach.

    Read-only. Prints the whole grid rather than one number, because a flat
    optimum means the data cannot separate 3 kb from 20 kb and should not be
    reported as though it could.
    """
    import json as json_module
    import os

    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.bam_coverage import compute_bam_depth, match_contigs
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reach_calibration import fit_reach, format_reach_table

    primers = collect_primers_from_args(
        getattr(args, "primers", None),
        getattr(args, "primers_file", None),
        name="primer",
    )

    validate_params_json_file(args.json_file)
    merge_args_to_parameter(args, parameter, ["json_file"])
    core_pipeline._initialize()

    fg_prefixes = core_pipeline.fg_prefixes
    fg_seq_lengths = core_pipeline.fg_seq_lengths

    aliases = {}
    for item in getattr(args, "contig_alias", None) or []:
        if "=" in item:
            key, value = item.split("=", 1)
            aliases[key] = value

    import pysam  # noqa: F401  (surfaces the [bam] extra's absence early)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        contig_map = match_contigs(
            list(bam.references),
            list(bam.lengths),
            fg_prefixes,
            fg_seq_lengths,
            aliases=aliases or None,
        )
    if not contig_map:
        raise SystemExit(
            f"No BAM contig could be matched to a foreground prefix. "
            f"Map one explicitly with --contig-alias FG=BAMCONTIG."
        )

    # Fit on the longest matched contig: the estimate is a length-scale, and a
    # short contig cannot distinguish reaches comparable to its own size.
    prefix = max(contig_map, key=lambda p: fg_seq_lengths[fg_prefixes.index(p)])
    contig = contig_map[prefix]
    length = fg_seq_lengths[fg_prefixes.index(prefix)]

    cache = PositionCache([prefix], list(dict.fromkeys(primers)))
    positions = sorted(
        {int(x) for p in primers for x in cache.get_positions(prefix, p, strand="both")}
    )
    if not positions:
        raise SystemExit(
            f"None of the {len(primers)} primers bind {prefix}. Reach cannot be "
            f"fitted without binding sites; check the primers and params.json."
        )

    logger.info(f"Fitting reach on {prefix} ({contig}, {length:,} bp, {len(positions)} sites)")
    depth = compute_bam_depth(args.bam, contig, length)
    fit = fit_reach(depth, positions, contig=contig)

    print(format_reach_table(fit))

    output = getattr(args, "output", None)
    if output:
        os.makedirs(os.path.dirname(output) or ".", exist_ok=True)
        with open(output, "w") as handle:
            json_module.dump(fit.to_dict(), handle, indent=2)
        logger.info(f"Wrote {output}")
    return 0 if fit.informative else 1


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
        # None, so `resolve_optimization_method` falls through to
        # params.json. Hardcoding "hybrid" here pinned the one-shot
        # `design` path to it from a second, independent direction.
        "optimization_method": None,
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


def _add_calibrate_reach_parser(subparsers):
    """Parser for `calibrate-reach`, extracted to keep `add_parsers` in budget."""
    parser = subparsers.add_parser(
        "calibrate-reach",
        help="Estimate the per-primer coverage reach from real sequencing "
        "depth (BAM), instead of assuming the polymerase default. Read-only.",
    )
    parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    parser.add_argument(
        "--primers", nargs="+", help="Primers used in the reaction that produced the BAM."
    )
    parser.add_argument("--primers-file", help="File with those primers (one per line).")
    parser.add_argument(
        "--bam",
        required=True,
        help="BAM from an SWGA reaction using these primers, mapped to the "
        "target genome. Requires the [bam] extra.",
    )
    parser.add_argument(
        "--contig-alias",
        action="append",
        default=None,
        metavar="FG=BAMCONTIG",
        help="Map a foreground prefix/basename to a BAM contig. Repeatable.",
    )
    parser.add_argument("--output", "-o", help="Write the fit as JSON to this path.")
    parser.add_argument("--quiet", "-q", action="store_true", help="Suppress progress output")


def add_parsers(subparsers):
    """Register this group's subcommands on the shared subparsers object.

    Called by neoswga.cli_unified.create_parser(). Extracted from the former
    monolithic create_parser() so each command group owns its argparse setup
    next to its handlers.
    """
    import argparse  # noqa: F401  (used by some command blocks)

    from neoswga.cli._common import add_common_options  # noqa: F401

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
        choices=_polymerase_names(),
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

    _add_calibrate_reach_parser(subparsers)

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
