"""Pipeline-step CLI handlers: count-kmers (run_step1), filter (run_step2),
score (run_step3), optimize (run_step4), the optimization-method guide, and
build-filter. Extracted from cli_unified.py.
"""

import json
import logging
import os
import sys

from neoswga.cli._common import (
    _record_run_manifest,
    check_jellyfish_available,
    collect_primers_from_args,
    merge_args_to_parameter,
    setup_gpu_acceleration,
    validate_params_json_file,
)

logger = logging.getLogger(__name__)

# Catchable at module level so the step handlers' except clauses resolve it.
try:
    from neoswga.core.pipeline import StepPrerequisiteError
except ImportError:  # pragma: no cover - defensive

    class StepPrerequisiteError(Exception):
        """Raised when a pipeline step's prerequisites are not satisfied."""

        pass


def run_step1(args):
    """Run step 1: K-mer preprocessing"""
    check_jellyfish_available()
    validate_params_json_file(args.json_file)
    import time as _time

    _t0 = _time.time()
    logger.info("Step 1: K-mer preprocessing starting...")
    try:
        from neoswga.core import parameter, pipeline

        # Set json_file if provided (needed for _initialize to load params)
        merge_args_to_parameter(args, parameter, ["json_file"])

        # Apply CLI k-mer range arguments BEFORE initialization
        # This ensures they're available when _initialize() reads defaults
        merge_args_to_parameter(args, parameter, ["min_k", "max_k"])

        # Initialize pipeline (lazy init, loads params from JSON)
        pipeline._initialize()

        # Re-apply CLI arguments after initialization to ensure they take precedence
        # (initialization may have overwritten them with JSON values)
        merge_args_to_parameter(args, parameter, ["min_k", "max_k"])

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration (boolean flag)
        if hasattr(args, "enable_qa") and args.enable_qa:
            parameter.enable_qa = True

        # Log effective parameters
        if not args.quiet:
            min_k = getattr(parameter, "min_k", 6)
            max_k = getattr(parameter, "max_k", 12)
            logger.info(f"K-mer range: {min_k}-{max_k}bp")

        # Run step1 (QA integration is only available for step2-4)
        if getattr(parameter, "enable_qa", False):
            logger.info(
                "Note: QA integration is available for filter/score/optimize steps, not count-kmers"
            )
        pipeline.step1()

        # Count k-mers for exclusion genome if specified via CLI
        excl_genome = getattr(args, "exclusion_genome", None)
        if excl_genome:
            if not os.path.exists(excl_genome):
                logger.error(f"Exclusion genome not found: {excl_genome}")
                sys.exit(1)
            from neoswga.core.kmer_counter import run_jellyfish as _run_jf

            excl_name = os.path.splitext(os.path.basename(excl_genome))[0]
            excl_prefix = os.path.join(parameter.data_dir, f"excl_{excl_name}")
            min_k = getattr(parameter, "min_k", 6)
            max_k = getattr(parameter, "max_k", 12)
            logger.info(f"Counting k-mers in exclusion genome: {excl_genome}")
            _run_jf(excl_genome, excl_prefix, min_k, max_k)
            # Store exclusion prefix for downstream steps
            parameter.excl_genomes = [excl_genome]
            parameter.excl_prefixes = [excl_prefix]

        # Count k-mers for blacklist genome(s) if specified via CLI
        bl_genomes_cli = getattr(args, "blacklist", None)
        if bl_genomes_cli:
            from neoswga.core.kmer_counter import run_jellyfish as _run_jf_bl

            for bl_genome in bl_genomes_cli:
                if not os.path.exists(bl_genome):
                    logger.error(f"Blacklist genome not found: {bl_genome}")
                    sys.exit(1)
            bl_prefixes_cli = []
            min_k = getattr(parameter, "min_k", 6)
            max_k = getattr(parameter, "max_k", 12)
            for bl_genome in bl_genomes_cli:
                bl_name = os.path.splitext(os.path.basename(bl_genome))[0]
                bl_prefix = os.path.join(parameter.data_dir, f"bl_{bl_name}")
                bl_prefixes_cli.append(bl_prefix)
                logger.info(f"Counting k-mers in blacklist genome: {bl_genome}")
                _run_jf_bl(bl_genome, bl_prefix, min_k, max_k)
            parameter.bl_genomes = bl_genomes_cli
            parameter.bl_prefixes = bl_prefixes_cli
            if args.bl_penalty is not None:
                parameter.bl_penalty = args.bl_penalty
            if args.max_bl_freq is not None:
                parameter.max_bl_freq = args.max_bl_freq

        _step1_inputs = (
            list(getattr(parameter, "fg_genomes", []) or [])
            + list(getattr(parameter, "bg_genomes", []) or [])
            + list(getattr(parameter, "bl_genomes", []) or [])
        )
        _record_run_manifest("count-kmers", args, parameter, input_files=_step1_inputs)

        _elapsed = _time.time() - _t0
        logger.info(f"Step 1 complete in {_elapsed:.1f}s")
        if not args.quiet:
            print("\nNext: neoswga filter -j params.json")

    except ImportError as e:
        logger.error(f"Pipeline not found: {e}")
        sys.exit(1)


def run_step2(args):
    """Run step 2: Candidate filtering with reaction conditions"""
    validate_params_json_file(args.json_file)
    import time as _time

    _t0 = _time.time()
    logger.info("Running step2 (candidate primer filtering)")

    try:
        from neoswga.core import parameter, pipeline

        # Set json_file if provided (needed for _initialize to load params)
        if hasattr(args, "json_file") and args.json_file:
            parameter.json_file = args.json_file

        # Initialize pipeline (lazy init, loads params from JSON)
        pipeline._initialize()

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration
        if hasattr(args, "enable_qa") and args.enable_qa:
            parameter.enable_qa = True

        # Apply preset if specified
        if hasattr(args, "preset") and args.preset:
            logger.info(f"Loading preset: {args.preset}")
            preset_params = load_preset_conditions(args.preset)

            # Apply preset to parameter module (but CLI args can override)
            for key, value in preset_params.items():
                if not hasattr(args, key) or getattr(args, key) is None:
                    setattr(parameter, key, value)
                else:
                    # CLI argument overrides preset
                    setattr(parameter, key, getattr(args, key))

        # Apply CLI arguments to parameter module using merger utility
        # GC filtering
        merge_args_to_parameter(args, parameter, ["gc_min", "gc_max"])

        # Special handling for gc_tolerance (calculates gc_min/gc_max from genome GC)
        if hasattr(args, "gc_tolerance") and args.gc_tolerance is not None:
            if (
                (not hasattr(args, "gc_min") or args.gc_min is None)
                and hasattr(parameter, "genome_gc")
                and parameter.genome_gc
            ):
                parameter.gc_min = max(0.20, parameter.genome_gc - args.gc_tolerance)
                parameter.gc_max = min(0.80, parameter.genome_gc + args.gc_tolerance)

        # Reaction conditions
        merge_args_to_parameter(args, parameter, ["reaction_temp", "na_conc", "mg_conc"])

        # Additives (with bsa -> bsa_ug_ml mapping)
        merge_args_to_parameter(
            args,
            parameter,
            [
                "dmso_percent",
                "betaine_m",
                "trehalose_m",
                "glycerol_percent",
                "peg_percent",
                "ethanol_percent",
                "urea_m",
                "tmac_m",
                "formamide_percent",
            ],
        )
        merge_args_to_parameter(args, parameter, ["bsa"], {"bsa": "bsa_ug_ml"})

        # Boolean flag for SSB (special handling)
        if hasattr(args, "ssb") and args.ssb:
            parameter.ssb = True

        # Traditional filtering parameters
        merge_args_to_parameter(
            args,
            parameter,
            [
                "min_fg_freq",
                "max_bg_freq",
                "max_gini",
                "max_primer",
                "min_tm",
                "max_tm",
                "max_dimer_bp",
                "max_self_dimer_bp",
            ],
        )

        # Polymerase
        merge_args_to_parameter(args, parameter, ["polymerase"])

        # Bloom filter for large background genomes
        use_bloom = getattr(args, "use_bloom_filter", False)
        bloom_path = getattr(args, "bloom_filter_path", None)
        sampled_path = getattr(args, "sampled_index_path", None)
        bloom_max_bg = getattr(args, "bloom_max_bg_matches", 10)

        if use_bloom or bloom_path:
            parameter.use_bloom_filter = True
            parameter.bloom_filter_path = bloom_path
            parameter.sampled_index_path = sampled_path
            parameter.bloom_max_bg_matches = bloom_max_bg
            if not args.quiet:
                logger.info(f"Bloom filter enabled for background filtering")
                if bloom_path:
                    logger.info(f"  Bloom filter: {bloom_path}")
                if sampled_path:
                    logger.info(f"  Sampled index: {sampled_path}")

        # Handle exclusion genome from CLI
        excl_genome = getattr(args, "exclusion_genome", None)
        if excl_genome:
            if not os.path.exists(excl_genome):
                logger.error(f"Exclusion genome not found: {excl_genome}")
                sys.exit(1)
            excl_name = os.path.splitext(os.path.basename(excl_genome))[0]
            excl_prefix = os.path.join(parameter.data_dir, f"excl_{excl_name}")
            parameter.excl_genomes = [excl_genome]
            parameter.excl_prefixes = [excl_prefix]
            excl_threshold = getattr(args, "excl_threshold", 0)
            if excl_threshold is not None:
                parameter.excl_threshold = excl_threshold
            if not args.quiet:
                logger.info(f"Exclusion genome: {excl_genome}")
                logger.info(f"  Threshold: {parameter.excl_threshold} (0 = reject any hit)")

        # Handle blacklist genome(s) from CLI
        bl_genomes_cli = getattr(args, "blacklist", None)
        if bl_genomes_cli:
            for bl_genome in bl_genomes_cli:
                if not os.path.exists(bl_genome):
                    logger.error(f"Blacklist genome not found: {bl_genome}")
                    sys.exit(1)
            bl_prefixes_cli = []
            for bl_genome in bl_genomes_cli:
                bl_name = os.path.splitext(os.path.basename(bl_genome))[0]
                bl_prefix = os.path.join(parameter.data_dir, f"bl_{bl_name}")
                bl_prefixes_cli.append(bl_prefix)
            parameter.bl_genomes = bl_genomes_cli
            parameter.bl_prefixes = bl_prefixes_cli
            # Compute blacklist lengths so the frequency denominator is correct.
            from neoswga.core import utility as _utility

            parameter.bl_seq_lengths = _utility.get_all_seq_lengths(
                fname_genomes=bl_genomes_cli, cpus=getattr(parameter, "cpus", 1)
            )
            if args.bl_penalty is not None:
                parameter.bl_penalty = args.bl_penalty
            if args.max_bl_freq is not None:
                parameter.max_bl_freq = args.max_bl_freq
            if not args.quiet:
                logger.info(f"Blacklist genomes: {bl_genomes_cli}")
                logger.info(f"  Max frequency: {parameter.max_bl_freq}")

        # Log effective parameters
        if not args.quiet:
            logger.info("Effective parameters:")
            if hasattr(parameter, "gc_min"):
                logger.info(f"  GC range: {parameter.gc_min:.3f} - {parameter.gc_max:.3f}")
            if hasattr(parameter, "reaction_temp"):
                logger.info(f"  Reaction temp: {parameter.reaction_temp}C")
            if hasattr(parameter, "polymerase"):
                logger.info(f"  Polymerase: {parameter.polymerase}")
            if hasattr(parameter, "dmso_percent") and parameter.dmso_percent > 0:
                logger.info(f"  DMSO: {parameter.dmso_percent}%")
            if hasattr(parameter, "betaine_m") and parameter.betaine_m > 0:
                logger.info(f"  Betaine: {parameter.betaine_m} M")

        # Run step2 with QA if enabled
        if getattr(parameter, "enable_qa", False):
            from neoswga.core import pipeline_qa_integration

            pipeline_qa_integration.run_step2_with_qa()
        else:
            pipeline.step2()

        if not args.quiet:
            print("\nNext: neoswga score -j params.json")

    except ImportError as e:
        logger.error(f"Failed to import pipeline module: {e}")
        logger.error("This may indicate a corrupted installation.")
        logger.error("Try: pip install -e . --force-reinstall")
        sys.exit(1)
    except StepPrerequisiteError as e:
        # Detailed error with remediation is in the exception message
        logger.error(str(e))
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"Required file not found: {e}")
        logger.error("Ensure Step 1 (count-kmers) has completed successfully.")
        logger.error("Run: neoswga count-kmers -j params.json")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Step 2 failed: {e}")
        if logger.level <= logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)
    _data_dir = getattr(parameter, "data_dir", None)
    _step2_out = os.path.join(_data_dir, "step2_df.csv") if _data_dir else None
    _record_run_manifest(
        "filter", args, parameter, input_files=[_step2_out] if _step2_out else None
    )
    _elapsed = _time.time() - _t0
    logger.info(f"Step 2 complete in {_elapsed:.1f}s")


def run_step3(args):
    """Run step 3: Random forest scoring"""
    validate_params_json_file(args.json_file)
    import time as _time

    _t0 = _time.time()
    try:
        from neoswga.core import parameter, pipeline

        # Set json_file and initialize pipeline
        merge_args_to_parameter(args, parameter, ["json_file"])
        pipeline._initialize()

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # Reproducible scoring: seed the k-mer sampling RNG if requested.
        if getattr(args, "seed", None) is not None:
            from neoswga.core.rf_preprocessing import set_kmer_sampling_seed

            set_kmer_sampling_seed(args.seed)
            if not args.quiet:
                logger.info(f"K-mer sampling seeded with {args.seed} for reproducible scoring")

        # QA integration (boolean flag)
        if hasattr(args, "enable_qa") and args.enable_qa:
            parameter.enable_qa = True

        # Apply CLI arguments to parameter module
        merge_args_to_parameter(args, parameter, ["min_amp_pred"])

        # Log effective parameters
        if not args.quiet:
            min_amp_pred = getattr(parameter, "min_amp_pred", None) or "auto"
            logger.info(f"Minimum amplification prediction score: {min_amp_pred}")

        # Enhanced feature engineering
        use_enhanced = getattr(args, "use_enhanced_features", False)
        enhanced_model_path = getattr(args, "enhanced_model_path", None)

        if use_enhanced:
            from neoswga.core.rf_preprocessing import is_enhanced_model_available

            if is_enhanced_model_available(enhanced_model_path):
                logger.info("Using enhanced 120+ feature model")
                parameter.use_enhanced_features = True
                parameter.enhanced_model_path = enhanced_model_path
            else:
                logger.warning("Enhanced model not found, using standard model")
                logger.warning("To train enhanced model: python scripts/train_enhanced_rf.py")
                use_enhanced = False

        # Scoring mode: fast (skip delta-G histograms) is the default.
        # --full-score opts back in to the full RF feature set.
        if getattr(args, "full_score", False):
            parameter.fast_score = False
            logger.info("Full scoring: computing thermodynamic histogram features (slow)")
        else:
            parameter.fast_score = True
            if getattr(args, "fast_score", False):
                logger.info("--fast-score is now the default; flag is a no-op")

        # Run step3 with QA if enabled
        if getattr(parameter, "enable_qa", False):
            from neoswga.core import pipeline_qa_integration

            pipeline_qa_integration.run_step3_with_qa()
        else:
            pipeline.step3()

        if not args.quiet:
            print("\nNext: neoswga optimize -j params.json")

    except ImportError as e:
        logger.error(f"Failed to import pipeline module: {e}")
        logger.error("This may indicate a corrupted installation.")
        logger.error("Try: pip install -e . --force-reinstall")
        sys.exit(1)
    except StepPrerequisiteError as e:
        # Detailed error with remediation is in the exception message
        logger.error(str(e))
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"Required file not found: {e}")
        logger.error("Ensure Step 2 (filter) has completed successfully.")
        logger.error("Run: neoswga filter -j params.json")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Step 3 failed: {e}")
        if logger.level <= logging.DEBUG:
            import traceback

            traceback.print_exc()
        sys.exit(1)
    _data_dir = getattr(parameter, "data_dir", None)
    _step3_in = os.path.join(_data_dir, "step2_df.csv") if _data_dir else None
    _step3_out = os.path.join(_data_dir, "step3_df.csv") if _data_dir else None
    _record_run_manifest(
        "score", args, parameter, input_files=[p for p in [_step3_in, _step3_out] if p]
    )
    _elapsed = _time.time() - _t0
    logger.info(f"Step 3 complete in {_elapsed:.1f}s")


def print_method_guide():
    """Print optimization method selection guide."""
    guide = """
OPTIMIZATION METHOD SELECTION GUIDE
====================================

Decision Tree:
  Speed critical / large pool?   -> dominating-set
  Clinical / low background?     -> background-aware
  Tm-balanced, dimer-aware?      -> network
  Not sure / want the best?      -> ensemble (runs all, keeps the best)
  Default                        -> hybrid

Methods:
+------------------+--------+-----------+-------------+--------------------------+
| Method           | Speed  | Coverage  | Specificity | Best For                 |
+------------------+--------+-----------+-------------+--------------------------+
| hybrid           | Medium | Excellent | Good        | General use (default)    |
| dominating-set   | Fast   | Excellent | Fair        | Large pools, quick       |
| network          | Medium | Good      | Good        | Tm-balanced, dimer-aware |
| background-aware | Slow   | Good      | Excellent   | Clinical, low background |
| ensemble         | Slow   | Best-of   | Best-of     | Run all, keep best by    |
|                  |        |           |             | normalized score         |
+------------------+--------+-----------+-------------+--------------------------+

Application weighting (--application) tunes how the normalized score and the
selection knobs trade coverage vs specificity:
  balanced | discovery (max coverage) | clinical (max specificity) |
  enrichment (sequencing) | metagenomics (capture diversity)

Usage Examples:
  # Default (hybrid)
  neoswga optimize -j params.json

  # Fast screening on a large candidate pool
  neoswga optimize -j params.json --optimization-method=dominating-set

  # Clinical samples (low background)
  neoswga optimize -j params.json --optimization-method=background-aware --application=clinical

  # Run several methods and keep the best (prints a comparison table)
  neoswga optimize -j params.json --optimization-method=ensemble
  neoswga optimize -j params.json --optimization-method=ensemble \\
      --ensemble-methods hybrid network background-aware

  # Pool the methods' primers and re-optimize over them (never worsens the best)
  neoswga optimize -j params.json --optimization-method=ensemble \\
      --ensemble-combine=union
"""
    print(guide)


def run_step4(args):
    """Run step 4: Primer set optimization (network-based + experimental)"""
    validate_params_json_file(args.json_file)
    import time as _time

    _t0 = _time.time()

    # Handle --method-guide option
    if getattr(args, "method_guide", False):
        print_method_guide()
        return

    logger.info("Primer set optimization")

    try:
        from neoswga.core import parameter

        # Set json_file if provided
        merge_args_to_parameter(args, parameter, ["json_file"])

        # Pass polymerase info to optimizer (use adapted value, not raw JSON,
        # so GC-adaptive strategy recommendations are respected)
        polymerase = getattr(parameter, "polymerase", "phi29")
        if polymerase != "phi29":
            logger.info(f"Polymerase: {polymerase} (config applied to optimizer)")

        logger.info(f"Optimization method: {args.optimization_method}")
        logger.info(f"Position cache: {args.use_position_cache}")
        logger.info(f"Background filter: {args.use_background_filter}")

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration (boolean flag)
        if hasattr(args, "enable_qa") and args.enable_qa:
            parameter.enable_qa = True

        # Set num_primers from CLI or JSON file (special handling needed)
        if hasattr(args, "num_primers") and args.num_primers:
            parameter.num_primers = args.num_primers
            parameter.target_set_size = args.num_primers
            logger.info(f"Target primer count: {args.num_primers}")
        elif not hasattr(parameter, "num_primers") or not parameter.num_primers:
            # Try to read from JSON if loaded
            json_data = getattr(parameter, "_json_data", {})
            num_primers = json_data.get("num_primers", json_data.get("target_set_size", 6))
            parameter.num_primers = num_primers
            parameter.target_set_size = num_primers

        # Automatic set size optimization (overrides num_primers if enabled)
        auto_size = getattr(args, "auto_size", False)
        application = getattr(args, "application", "enrichment")
        min_fg_bg_ratio = getattr(args, "min_fg_bg_ratio", None)
        show_frontier = getattr(args, "show_frontier", False)
        quick_estimate = getattr(args, "quick_estimate", False)

        if auto_size:
            logger.info("=" * 60)
            logger.info("Automatic Set Size Optimization")
            logger.info("=" * 60)

            try:
                from neoswga.core.mechanistic_model import MechanisticModel
                from neoswga.core.reaction_conditions import ReactionConditions
                from neoswga.core.set_size_optimizer import (
                    create_baseline_effects,
                    recommend_set_size,
                )

                # Get genome information
                json_data = getattr(parameter, "_json_data", {})
                fg_seq_lengths = getattr(parameter, "fg_seq_lengths", [])
                genome_length = (
                    sum(fg_seq_lengths)
                    if fg_seq_lengths
                    else json_data.get("fg_genome_length", 1_000_000)
                )
                primer_length = json_data.get("max_k", 12)

                # Get template GC content
                template_gc = getattr(args, "template_gc", None)
                if template_gc is None:
                    # Try to calculate from genome or use default
                    template_gc = json_data.get("fg_gc", 0.5)

                # Get polymerase and create conditions
                polymerase = json_data.get("polymerase", "phi29")
                reaction_temp = json_data.get(
                    "reaction_temp", 30.0 if polymerase == "phi29" else 42.0
                )

                try:
                    conditions = ReactionConditions(
                        temp=reaction_temp,
                        polymerase=polymerase,
                        mg_conc=json_data.get("mg_conc", 2.5),
                        dmso_percent=json_data.get("dmso_percent", 0.0),
                        betaine_m=json_data.get("betaine_m", 0.0),
                    )
                    model = MechanisticModel(conditions)
                    sample_primer = "A" * primer_length  # Neutral sequence
                    mech_effects = model.calculate_effects(sample_primer, template_gc)
                except Exception as e:
                    logger.warning(f"Could not create mechanistic model: {e}")
                    mech_effects = create_baseline_effects()

                # Per-primer reach for the recommender. Uses realistic per-primer
                # reach (Phase 16) so the prediction matches what the optimizer
                # actually measures via base_optimizer.compute_metrics
                # (extension_reach=3000 bp default). Passing theoretical
                # processivity here under-counted required primers ~10-20x.
                from neoswga.core.coverage import polymerase_extension_reach

                processivity = polymerase_extension_reach(polymerase, coverage_metric="realistic")

                # Get recommendation (supports new min_fg_bg_ratio parameter)
                recommendation = recommend_set_size(
                    application=application,
                    genome_length=genome_length,
                    primer_length=primer_length,
                    mech_effects=mech_effects,
                    processivity=processivity,
                    min_fg_bg_ratio=min_fg_bg_ratio,
                )

                auto_num_primers = recommendation["recommended_size"]
                parameter.num_primers = auto_num_primers
                parameter.target_set_size = auto_num_primers

                logger.info(f"Application profile: {application}")
                logger.info(f"  Priority: {recommendation['priority']}")
                logger.info(f"  {recommendation['rationale']}")
                logger.info(f"Genome length: {genome_length:,} bp")
                logger.info(f"Primer length: {primer_length} bp")
                logger.info(f"Target coverage: {recommendation['target_coverage']:.0%}")
                logger.info(f"Min fg/bg ratio: {recommendation['min_fg_bg_ratio']:.1f}")
                logger.info(f"Recommended set size: {auto_num_primers} primers")
                logger.info(
                    f"  (typical range: {recommendation['size_range'][0]}-{recommendation['size_range'][1]})"
                )
                logger.info("=" * 60)

            except Exception as e:
                logger.warning(f"Auto-size failed: {e}. Using default num_primers.")

        # Cooperative binding (experimental - not yet fully integrated)
        if hasattr(args, "use_cooperative_binding") and args.use_cooperative_binding:
            logger.warning("--use-cooperative-binding is experimental and not yet fully integrated")
            parameter.use_cooperative_binding = True

        # Mechanistic model: when enabled, the optimizer's
        # factory kwargs carry mechanistic_weight so NetworkOptimizer (the
        # condition-aware scoring path inside hybrid / network) will
        # construct a MechanisticModel(conditions) and fold a
        # `mechanistic_weight * predicted_amplification_factor` term into
        # its per-primer score. See network_optimizer.py:~630 where
        # self.mech_model is built from self.conditions. Phase 13B.
        _use_mech = getattr(args, "use_mechanistic_model", False)
        # Footgun fix: a non-default --mechanistic-weight implies enabling the
        # model. Previously the weight was silently ignored unless
        # --use-mechanistic-model was also passed ("I set 0.5 and nothing changed").
        _explicit_weight = getattr(args, "mechanistic_weight", None)
        if not _use_mech and _explicit_weight is not None and _explicit_weight > 0:
            _use_mech = True
            logger.info(
                "--mechanistic-weight > 0 implies --use-mechanistic-model; "
                "enabling the mechanistic model."
            )
        _mech_weight = (
            (_explicit_weight if _explicit_weight is not None else 0.3) if _use_mech else 0.0
        )
        if _use_mech:
            logger.info(
                f"Mechanistic model scoring enabled (weight={_mech_weight:.2f}). "
                "Network / hybrid optimizers will apply a weighted mechanistic "
                "amplification-factor term to every primer."
            )

        # Primer strategy
        merge_args_to_parameter(args, parameter, ["primer_strategy"])
        if hasattr(args, "primer_strategy") and args.primer_strategy == "hybrid":
            logger.info("Using hybrid primer strategy (mixed lengths)")

        # Use unified optimizer framework (all methods handled via factory pattern)
        from neoswga.core.unified_optimizer import list_available_optimizers, optimize_step4

        # Get uniformity weight from args (default 0.0 for backward compatibility)
        uniformity_weight = getattr(args, "uniformity_weight", 0.0)
        if uniformity_weight > 0:
            logger.info(f"Uniformity weight: {uniformity_weight:.2f}")

        # Get minimal primer selection options
        minimize_primers = getattr(args, "minimize_primers", False)
        target_coverage = getattr(args, "target_coverage", 0.70)
        if minimize_primers:
            logger.info(
                f"Minimal primer selection enabled (target coverage: {target_coverage:.1%})"
            )

        strategy = getattr(args, "scoring_weights", "balanced")

        # Background pre-filter flag (enabled by default, --no-bg-prefilter disables)
        bg_prefilter = not getattr(args, "no_bg_prefilter", False)

        # Host-free mode: optimize without background genome data
        no_background = getattr(args, "no_background", False)
        if no_background:
            logger.info("Host-free mode: background genome data will be ignored")

        # Run unified optimization
        # Random seed for reproducibility (stochastic optimizers)
        seed = getattr(args, "seed", None)
        if seed is not None:
            logger.info(f"Random seed: {seed}")

        # Pass explicit target_size from parameter module to avoid re-initialization override
        target_size = getattr(parameter, "target_set_size", getattr(parameter, "num_primers", 6))

        results, scores, cache = optimize_step4(
            use_cache=args.use_position_cache,
            use_background_filter=args.use_background_filter,
            optimization_method=args.optimization_method,
            verbose=not args.quiet,
            uniformity_weight=uniformity_weight,
            minimize_primers=minimize_primers,
            target_coverage=target_coverage,
            strategy=strategy,  # Pass strategy for normalized optimizer
            polymerase=polymerase,  # Pass polymerase for hybrid preset config
            bg_prefilter=bg_prefilter,  # Background pre-filtering of candidates
            no_background=no_background,  # Host-free mode
            seed=seed,  # Reproducibility for stochastic optimizers
            target_size=target_size,  # Explicit target from params
            mechanistic_weight=_mech_weight,  # 0 if flag not set; else user value
            # Phase 11D / 14C — per-target coverage floor and application
            # profile. Both are consumed by unified_optimizer.run_optimization.
            min_per_target_coverage=getattr(args, "min_per_target_coverage", 0.0),
            application=getattr(args, "application", "balanced"),
            ensemble_methods=getattr(args, "ensemble_methods", None),
            ensemble_combine=getattr(args, "ensemble_combine", "best"),
        )

        if results:
            target_size = getattr(parameter, "num_primers", 6)
            target_size = getattr(parameter, "target_set_size", target_size)
            num_found = len(results[0])
            logger.info(f"Selected {num_found} primers")
            logger.info(f"Score: {scores[0]:.4f}")
            if num_found < target_size:
                logger.warning(
                    f"WARNING: Found {num_found} primers but target was "
                    f"{target_size} (PARTIAL result: insufficient candidates)"
                )
                logger.warning("To improve, consider:")
                logger.warning("  - Relaxing filter thresholds (max_bg_freq, max_gini)")
                logger.warning("  - Widening k-mer range (min_k / max_k)")
                logger.warning("  - Increasing candidate pool (max_primer)")
                logger.warning(
                    f"  - Trying a different optimizer " f"(current: {args.optimization_method})"
                )
        else:
            logger.error("No primer sets found. Optimization failed.")
            sys.exit(1)

        # MOEA Pareto front text table (Phase 14B). When MOEA has populated
        # OptimizationResult.pareto_front, emit a plain-text tradeoff table
        # in addition to the set-size frontier (which runs below). This
        # surfaces NSGA-III's non-dominated solutions without needing Plotly
        # or pandas.
        if show_frontier and results:
            from neoswga.core import unified_optimizer as _uo

            _last = getattr(_uo, "_LAST_RESULT", None)
            if _last is not None and getattr(_last, "pareto_front", None):
                logger.info("")
                logger.info("=" * 60)
                logger.info("MOEA Pareto Front (non-dominated solutions)")
                logger.info("=" * 60)
                header = (
                    f"{'#':>3}  {'coverage':>9}  {'n_primers':>10}  "
                    f"{'bg_binding':>11}  {'uniformity':>11}  {'score':>7}"
                )
                logger.info(header)
                logger.info("-" * len(header))
                for idx, pm in enumerate(_last.pareto_metrics or [], start=1):
                    row = (
                        f"{idx:>3}  "
                        f"{pm.get('target_coverage', 0.0):>9.3f}  "
                        f"{pm.get('n_primers', 0):>10}  "
                        f"{pm.get('background_binding', 0.0):>11.3f}  "
                        f"{pm.get('uniformity', 0.0):>11.3f}  "
                        f"{pm.get('score', 0.0):>7.3f}"
                    )
                    logger.info(row)

        # Ensemble per-method comparison table. Printed whenever an ensemble
        # run populated OptimizationResult.ensemble_comparison, so the user
        # sees which method won and by how much.
        if results:
            from neoswga.core import unified_optimizer as _uo

            _last = getattr(_uo, "_LAST_RESULT", None)
            _cmp = getattr(_last, "ensemble_comparison", None) if _last is not None else None
            if _cmp:
                logger.info("")
                logger.info("=" * 60)
                logger.info("Ensemble comparison (best by normalized score)")
                logger.info("=" * 60)
                header = (
                    f"{'method':>16}  {'norm_score':>10}  {'n':>3}  "
                    f"{'fg_cov':>7}  {'bg_cov':>7}  {'status':>9}  sel"
                )
                logger.info(header)
                logger.info("-" * len(header))
                for row_d in _cmp:
                    logger.info(
                        f"{row_d.get('method', '?'):>16}  "
                        f"{row_d.get('normalized_score', 0.0):>10.4f}  "
                        f"{row_d.get('n_primers', 0):>3}  "
                        f"{row_d.get('fg_coverage', 0.0):>7.3f}  "
                        f"{row_d.get('bg_coverage', 0.0):>7.3f}  "
                        f"{row_d.get('status', '?'):>9}  "
                        f"{'*' if row_d.get('selected') else ''}"
                    )

        # Show Pareto frontier analysis (optional)
        if show_frontier and results and cache is not None:
            try:
                import pandas as pd

                from neoswga.core.pareto_frontier import (
                    generate_frontier_report,
                    plot_frontier,
                    summarize_frontier_for_cli,
                )
                from neoswga.core.set_size_optimizer import (
                    ParetoFrontierGenerator,
                    select_from_frontier,
                )

                logger.info("")
                logger.info("=" * 60)
                logger.info("Pareto Frontier Analysis")
                logger.info("=" * 60)

                # Load step2 or step3 DataFrame for primer pool
                data_dir = parameter.data_dir
                step3_file = os.path.join(data_dir, "step3_df.csv")
                step2_file = os.path.join(data_dir, "step2_df.csv")

                if os.path.exists(step3_file):
                    primer_pool = pd.read_csv(step3_file)
                elif os.path.exists(step2_file):
                    primer_pool = pd.read_csv(step2_file)
                else:
                    raise FileNotFoundError("No primer pool CSV found")

                # Get genome lengths
                fg_lengths = getattr(parameter, "fg_lengths", [1_000_000])
                bg_lengths = getattr(parameter, "bg_lengths", [])
                fg_prefixes = parameter.fg_prefixes
                bg_prefixes = getattr(parameter, "bg_prefixes", [])

                # Create frontier generator. The generator uses `processivity`
                # as the coverage read-length, so pass the REALISTIC per-primer
                # reach (phi29 ~3 kb), not single-molecule processivity (70 kb),
                # which would inflate the frontier's coverage estimates.
                from neoswga.core.coverage import polymerase_extension_reach

                _frontier_reach = polymerase_extension_reach(
                    getattr(parameter, "polymerase", "phi29") or "phi29",
                    coverage_metric="realistic",
                )
                generator = ParetoFrontierGenerator(
                    primer_pool=primer_pool,
                    position_cache=cache,
                    fg_prefixes=fg_prefixes,
                    bg_prefixes=bg_prefixes,
                    fg_seq_lengths=fg_lengths,
                    bg_seq_lengths=bg_lengths,
                    processivity=_frontier_reach,
                )

                # Generate frontier (quick estimation only if requested)
                frontier_result = generator.generate_frontier(
                    min_size=4,
                    max_size=min(20, len(primer_pool)),
                    quick_only=quick_estimate,
                    verbose=not args.quiet,
                )

                # Select from frontier based on application
                selected, explanation = select_from_frontier(
                    frontier_result.pareto_points,
                    application=application,
                    min_fg_bg_ratio=min_fg_bg_ratio,
                )
                frontier_result.selected_point = selected
                frontier_result.selection_explanation = explanation

                # Display summary
                logger.info(summarize_frontier_for_cli(frontier_result, application))
                logger.info("")
                logger.info(explanation)

                # Try to save plot
                try:
                    fig = plot_frontier(frontier_result, application=application)
                    plot_path = os.path.join(data_dir, "pareto_frontier.png")
                    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                    logger.info(f"Pareto frontier plot saved to: {plot_path}")
                    import matplotlib.pyplot as plt

                    plt.close(fig)
                except Exception as e:
                    logger.debug(f"Could not save frontier plot: {e}")

                logger.info("=" * 60)

            except Exception as e:
                logger.warning(f"Pareto frontier analysis failed: {e}")
                import traceback

                logger.debug(traceback.format_exc())

        # Phase 5: Stochastic validation (optional)
        validate_simulation = getattr(args, "validate_simulation", False)
        if validate_simulation and results:
            simulation_time = getattr(args, "simulation_time", 3600.0)
            logger.info("=" * 60)
            logger.info("Stochastic Validation (Gillespie Algorithm)")
            logger.info("=" * 60)
            logger.info(f"Simulation time: {simulation_time/3600:.1f} hours")

            try:
                from neoswga.core.amplicon_network import AmpliconNetwork
                from neoswga.core.stochastic_simulator import validate_network_predictions

                primers = results[0]

                # Build simplified networks for validation
                # Get genome lengths and positions from cache
                fg_prefixes = parameter.fg_prefixes
                bg_prefixes = getattr(parameter, "bg_prefixes", [])

                if cache is not None:
                    # Build AmpliconNetwork for target
                    fg_genome_length = sum(getattr(parameter, "fg_lengths", [1000000]))
                    bg_genome_length = sum(getattr(parameter, "bg_lengths", [1000000]))

                    # Create networks from cache
                    from neoswga.core.network_optimizer import AmplificationNetwork

                    # Build fg network
                    fg_network = AmplificationNetwork(max_extension=70000)
                    for primer in primers:
                        for prefix in fg_prefixes:
                            try:
                                fw_pos = cache.get_positions(prefix, primer, "forward")
                                rv_pos = cache.get_positions(prefix, primer, "reverse")
                                for pos in fw_pos:
                                    fg_network.add_site(pos, primer, "forward")
                                for pos in rv_pos:
                                    fg_network.add_site(pos, primer, "reverse")
                            except Exception:
                                pass
                    fg_network.build_graph()

                    # Build bg network
                    bg_network = AmplificationNetwork(max_extension=70000)
                    for primer in primers:
                        for prefix in bg_prefixes:
                            try:
                                fw_pos = cache.get_positions(prefix, primer, "forward")
                                rv_pos = cache.get_positions(prefix, primer, "reverse")
                                for pos in fw_pos:
                                    bg_network.add_site(pos, primer, "forward")
                                for pos in rv_pos:
                                    bg_network.add_site(pos, primer, "reverse")
                            except Exception:
                                pass
                    bg_network.build_graph()

                    # Run validation
                    validation_results = validate_network_predictions(
                        primers, fg_network, bg_network, seed=getattr(args, "seed", None)
                    )

                    # Report results
                    predicted = validation_results["predicted"]
                    simulated = validation_results["simulated"]
                    validation = validation_results["validation"]

                    logger.info("\nPrediction vs Simulation:")
                    logger.info(f"  Predicted enrichment: {predicted['enrichment']:.0f}x")
                    logger.info(f"  Simulated enrichment: {simulated['enrichment']:.0f}x")
                    logger.info(f"  Prediction error: {validation['prediction_error']:.1%}")

                    if validation["prediction_accurate"]:
                        logger.info("  Status: VALIDATED (within 50% of simulation)")
                    else:
                        logger.warning("  Status: DIVERGENT (>50% error)")
                        logger.warning("  Consider using stochastic optimizer for better accuracy")
                else:
                    logger.warning("Position cache not available for validation")

            except ImportError as e:
                logger.warning(f"Stochastic validation not available: {e}")
            except Exception as e:
                logger.warning(f"Validation failed: {e}")

        _data_dir = getattr(parameter, "data_dir", None)
        _step4_in = os.path.join(_data_dir, "step3_df.csv") if _data_dir else None
        _step4_out = os.path.join(_data_dir, "step4_improved_df.csv") if _data_dir else None
        _record_run_manifest(
            "optimize", args, parameter, input_files=[p for p in [_step4_in, _step4_out] if p]
        )

        _elapsed = _time.time() - _t0
        logger.info(f"Step 4 complete in {_elapsed:.1f}s")
        if not args.quiet:
            data_dir = getattr(parameter, "data_dir", ".")
            print(f"\nDone! View results:")
            print(f"  neoswga interpret -d {data_dir}")
            print(f"  neoswga report -d {data_dir}")
            print(f"  neoswga export -d {data_dir} --format fasta")

    except StepPrerequisiteError as e:
        # Detailed error with remediation is in the exception message
        logger.error(str(e))
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"Required file not found: {e}")
        logger.error("Ensure Step 3 (score) has completed successfully.")
        logger.error("Run: neoswga score -j params.json")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Step 4 failed: {e}")
        if logger.level <= logging.DEBUG:
            import traceback

            traceback.print_exc()
        else:
            logger.error("Run with --verbose for full traceback")
        sys.exit(1)


def run_build_filter(args):
    """Build background Bloom filter"""
    import pickle

    from neoswga.core.background_filter import BackgroundBloomFilter, SampledGenomeIndex

    logger.info("Building background filter")
    logger.info(f"Input: {args.genome}")
    logger.info(f"Output: {args.output_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    bloom_path = os.path.join(args.output_dir, "bg_bloom.pkl")
    sampled_path = os.path.join(args.output_dir, "bg_sampled.pkl")

    if not args.force and os.path.exists(bloom_path):
        logger.warning("Filter already exists. Use --force to rebuild.")
        return

    try:
        min_k = getattr(args, "min_k", 6)
        max_k = getattr(args, "max_k", 12)

        if getattr(args, "from_kmers", False):
            # Build from pre-computed k-mer files (MUCH faster)
            logger.info(f"Building from k-mer files (prefix: {args.genome})")
            logger.info(f"K-mer range: {min_k}-{max_k}bp")

            # Estimate capacity from k-mer file sizes
            capacity = args.capacity
            if capacity is None:
                total_kmers = 0
                for k in range(min_k, max_k + 1):
                    fpath = f"{args.genome}_{k}mer_all.txt"
                    if os.path.exists(fpath):
                        with open(fpath) as f:
                            total_kmers += sum(1 for _ in f)
                capacity = max(total_kmers * 2, 10000000)  # 2x k-mers or min 10M
                logger.info(f"Auto-detected capacity: {capacity:,} (from {total_kmers:,} k-mers)")

            # Build Bloom filter from k-mer files
            bloom = BackgroundBloomFilter(capacity=capacity, error_rate=args.error_rate)
            bloom.add_from_kmer_files(args.genome, min_k=min_k, max_k=max_k)

            # Build sampled index from k-mer files (simpler - just use the counts)
            logger.info("Building sampled index from k-mer files...")
            sampled = SampledGenomeIndex(
                sample_rate=1
            )  # rate=1 since k-mer files are already unique
            for k in range(min_k, max_k + 1):
                fpath = f"{args.genome}_{k}mer_all.txt"
                if os.path.exists(fpath):
                    with open(fpath) as f:
                        for line in f:
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                kmer, count = parts[0], int(parts[1])
                                sampled.kmers[kmer] = count

            logger.info(f"Sampled index built: {len(sampled.kmers):,} k-mers with counts")

        else:
            # Build from genome FASTA (slower but comprehensive)
            logger.info(f"Building from genome FASTA: {args.genome}")

            # Auto-detect capacity from genome size
            capacity = args.capacity
            if capacity is None:
                from Bio import SeqIO

                total_size = 0
                for record in SeqIO.parse(args.genome, "fasta"):
                    total_size += len(record.seq)
                capacity = total_size * 10  # 10x genome size
                logger.info(f"Auto-detected genome size: {total_size:,} bp")
                logger.info(f"Using capacity: {capacity:,}")

            bloom = BackgroundBloomFilter(capacity=capacity, error_rate=args.error_rate)
            bloom.add_genome(args.genome, include_mismatches=False, min_k=min_k, max_k=max_k)

            sampled = SampledGenomeIndex(sample_rate=100)
            sampled.add_genome(args.genome, min_k=min_k, max_k=max_k)

        # Save filters
        bloom.save(bloom_path)
        sampled.save(sampled_path)

        logger.info("Filter built successfully!")
        logger.info(f"Bloom filter: {bloom_path} ({bloom.memory_usage_mb():.1f} MB)")
        logger.info(f"Sampled index: {sampled_path}")
        logger.info(f"Total k-mers indexed: {bloom.kmer_count:,}")

    except Exception as e:
        logger.error(f"Failed to build filter: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
