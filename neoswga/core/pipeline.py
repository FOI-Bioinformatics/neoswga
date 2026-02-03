import warnings
import argparse
import neoswga.core.parameter as parameter
import neoswga.core.utility as utility
from neoswga.core.kmer_counter import run_jellyfish, get_primer_list_from_kmers
import neoswga.core.filter as filter_module
import neoswga.core.string_search as string_search
import neoswga.core.rf_preprocessing as rf_preprocessing
from neoswga.core.progress import progress_context
import multiprocessing
import json
import sys
import pandas as pd
import pickle
import numpy as np
import os
import logging
from typing import List, Tuple, Optional
from dataclasses import dataclass

logger = logging.getLogger(__name__)


# =============================================================================
# Step Prerequisite Validation
# =============================================================================

@dataclass
class StepValidationResult:
    """Result of step prerequisite validation."""
    valid: bool
    missing_files: List[str]
    error_message: str
    remediation: str


class StepPrerequisiteError(Exception):
    """Raised when step prerequisites are not met."""

    def __init__(self, step: int, validation: StepValidationResult):
        self.step = step
        self.validation = validation
        message = f"\n{'='*60}\n"
        message += f"STEP {step} PREREQUISITE ERROR\n"
        message += f"{'='*60}\n"
        message += f"\n{validation.error_message}\n"
        if validation.missing_files:
            message += f"\nMissing files:\n"
            for f in validation.missing_files[:5]:  # Show first 5
                message += f"  - {f}\n"
            if len(validation.missing_files) > 5:
                message += f"  ... and {len(validation.missing_files) - 5} more\n"
        message += f"\nTo fix this:\n  {validation.remediation}\n"
        message += f"{'='*60}\n"
        super().__init__(message)


def validate_step1_prerequisites(data_dir: str, fg_genomes: List[str], bg_genomes: List[str]) -> StepValidationResult:
    """
    Validate prerequisites for Step 1 (k-mer counting).

    Checks:
    - Foreground genome files exist
    - Background genome files exist (if specified)
    - Data directory is writable
    """
    missing = []

    # Check foreground genomes
    for genome in fg_genomes:
        if not os.path.exists(genome):
            missing.append(genome)

    # Check background genomes
    for genome in bg_genomes:
        if not os.path.exists(genome):
            missing.append(genome)

    if missing:
        return StepValidationResult(
            valid=False,
            missing_files=missing,
            error_message="Genome files not found. Step 1 requires FASTA files for foreground and background genomes.",
            remediation="Check paths in params.json (fasta_fore, fasta_back) or use --fasta-fore and --fasta-back CLI options."
        )

    # Check data directory
    if data_dir and not os.path.exists(data_dir):
        try:
            os.makedirs(data_dir)
        except OSError as e:
            return StepValidationResult(
                valid=False,
                missing_files=[data_dir],
                error_message=f"Cannot create data directory: {e}",
                remediation="Check write permissions or specify a different --data-dir."
            )

    return StepValidationResult(valid=True, missing_files=[], error_message="", remediation="")


def validate_step2_prerequisites(data_dir: str, fg_prefixes: List[str], bg_prefixes: List[str],
                                  min_k: int = 6, max_k: int = 12) -> StepValidationResult:
    """
    Validate prerequisites for Step 2 (filtering).

    Checks:
    - K-mer count files exist from Step 1
    """
    missing = []

    # Check k-mer files for each prefix and k value
    for prefix in fg_prefixes + bg_prefixes:
        for k in range(min_k, max_k + 1):
            kmer_file = f"{prefix}_{k}mer_all.txt"
            if not os.path.exists(kmer_file):
                missing.append(kmer_file)

    if missing:
        return StepValidationResult(
            valid=False,
            missing_files=missing,
            error_message="K-mer count files not found. Step 2 requires output from Step 1.",
            remediation="Run 'neoswga count-kmers -j params.json' (Step 1) first."
        )

    return StepValidationResult(valid=True, missing_files=[], error_message="", remediation="")


def validate_step3_prerequisites(data_dir: str) -> StepValidationResult:
    """
    Validate prerequisites for Step 3 (scoring).

    Checks:
    - step2_df.csv exists from Step 2
    """
    step2_file = os.path.join(data_dir, "step2_df.csv")

    if not os.path.exists(step2_file):
        return StepValidationResult(
            valid=False,
            missing_files=[step2_file],
            error_message="Step 2 output not found. Step 3 requires filtered primers from Step 2.",
            remediation="Run 'neoswga filter -j params.json' (Step 2) first."
        )

    # Check file is not empty
    try:
        df = pd.read_csv(step2_file)
        if len(df) == 0:
            return StepValidationResult(
                valid=False,
                missing_files=[],
                error_message="Step 2 output is empty (no primers passed filtering).",
                remediation="Relax filtering parameters (--min-fg-freq, --max-bg-freq, --max-gini) and re-run Step 2."
            )
    except Exception as e:
        return StepValidationResult(
            valid=False,
            missing_files=[step2_file],
            error_message=f"Cannot read Step 2 output: {e}",
            remediation="Re-run 'neoswga filter -j params.json' (Step 2)."
        )

    return StepValidationResult(valid=True, missing_files=[], error_message="", remediation="")


def validate_step4_prerequisites(data_dir: str, fg_prefixes: List[str]) -> StepValidationResult:
    """
    Validate prerequisites for Step 4 (optimization).

    Checks:
    - step3_df.csv exists from Step 3
    - Position HDF5 files exist for optimization
    """
    step3_file = os.path.join(data_dir, "step3_df.csv")

    if not os.path.exists(step3_file):
        return StepValidationResult(
            valid=False,
            missing_files=[step3_file],
            error_message="Step 3 output not found. Step 4 requires scored primers from Step 3.",
            remediation="Run 'neoswga score -j params.json' (Step 3) first."
        )

    # Check file has primers
    try:
        df = pd.read_csv(step3_file)
        if len(df) == 0:
            return StepValidationResult(
                valid=False,
                missing_files=[],
                error_message="Step 3 output is empty (no primers passed scoring threshold).",
                remediation="Lower --min-amp-pred threshold and re-run Step 3."
            )
    except Exception as e:
        return StepValidationResult(
            valid=False,
            missing_files=[step3_file],
            error_message=f"Cannot read Step 3 output: {e}",
            remediation="Re-run 'neoswga score -j params.json' (Step 3)."
        )

    # Check position files (at least one should exist)
    missing_positions = []
    found_any = False
    for prefix in fg_prefixes:
        for k in range(6, 13):
            pos_file = f"{prefix}_{k}mer_positions.h5"
            if os.path.exists(pos_file):
                found_any = True
                break
        if found_any:
            break

    if not found_any:
        return StepValidationResult(
            valid=False,
            missing_files=missing_positions,
            error_message="Position files not found. Step 4 requires binding positions from Step 2.",
            remediation="Position files are generated during Step 2. Re-run 'neoswga filter -j params.json'."
        )

    return StepValidationResult(valid=True, missing_files=[], error_message="", remediation="")


defaults = {
    "min_fg_freq": float(1 / 100000),
    "max_bg_freq": float(1 / 200000),
    "max_gini": 0.6,
    "max_primer": 500,
    "min_amp_pred": 10,
    "min_tm": 15,
    "max_tm": 45,
    "max_dimer_bp": 3,
    "max_self_dimer_bp": 4,
    "selection_metric": "deterministic",
    "iterations": 8,
    "top_set_count": 10,
    "retries": 5,
    "max_sets": 5,
    "fg_circular": True,
    "bg_circular": False,
    "drop_iterations": 5,
    "verbose": False,
    "cpus": int(multiprocessing.cpu_count()),
}

# Module-level variables (will be initialized lazily)
_initialized = False
params = None
fg_prefixes = None
bg_prefixes = None
fg_genomes = None
bg_genomes = None
fg_seq_lengths = None
bg_seq_lengths = None
fg_circular = None
bg_circular = None

def _initialize():
    """Lazy initialization - only parse args when actually needed"""
    global _initialized, params, fg_prefixes, bg_prefixes, fg_genomes, bg_genomes
    global fg_seq_lengths, bg_seq_lengths, fg_circular, bg_circular

    if _initialized:
        return

    # Create empty options object that returns None for any attribute
    # (CLI arguments are handled by cli_unified.py)
    class EmptyOptions:
        def __init__(self):
            # Get json_file from parameter module if it was set by CLI
            self.json_file = getattr(parameter, 'json_file', None)
        def __getattr__(self, name):
            return None
    options = EmptyOptions()
    params = parameter.get_params(options)

    fg_prefixes = params["fg_prefixes"]
    bg_prefixes = params["bg_prefixes"]

    fg_genomes = params["fg_genomes"]
    bg_genomes = params["bg_genomes"]

    fg_seq_lengths = params["fg_seq_lengths"]
    bg_seq_lengths = params["bg_seq_lengths"]

    fg_circular = params["fg_circular"]
    bg_circular = params["bg_circular"]

    # Apply GC-adaptive strategy if genome_gc is available
    _apply_gc_adaptive_defaults()

    _initialized = True


def _apply_gc_adaptive_defaults():
    """
    Apply GC-adaptive parameter defaults if genome_gc is set.

    Only applies defaults for parameters that weren't explicitly set by the user.
    This enables automatic optimization for different genome types without
    requiring manual configuration.
    """
    genome_gc = getattr(parameter, 'genome_gc', None)
    if genome_gc is None:
        return  # No genome GC specified, use explicit parameters

    try:
        from neoswga.core.gc_adaptive_strategy import GCAdaptiveStrategy

        strategy = GCAdaptiveStrategy(genome_gc_content=genome_gc)
        adaptive_params = strategy.get_parameters()

        # Apply polymerase if not explicitly set
        current_polymerase = getattr(parameter, 'polymerase', None)
        if current_polymerase is None or current_polymerase == 'phi29':
            # Only override if using default phi29 and adaptive recommends different
            if adaptive_params.recommended_polymerase != 'phi29':
                logger.info(f"GC-adaptive: Using {adaptive_params.recommended_polymerase} "
                           f"for {genome_gc:.1%} GC genome")
                parameter.polymerase = adaptive_params.recommended_polymerase

        # Apply reaction temperature if not explicitly set
        current_temp = getattr(parameter, 'reaction_temp', None)
        if current_temp is None:
            parameter.reaction_temp = adaptive_params.reaction_temp
            logger.info(f"GC-adaptive: Setting reaction temp to {adaptive_params.reaction_temp}C")

        # Apply k-mer range if using defaults
        current_min_k = getattr(parameter, 'min_k', None)
        current_max_k = getattr(parameter, 'max_k', None)
        if current_min_k is None or current_max_k is None:
            if current_min_k is None:
                parameter.min_k = adaptive_params.kmer_range[0]
            if current_max_k is None:
                parameter.max_k = adaptive_params.kmer_range[1]
            logger.info(f"GC-adaptive: Setting k-mer range to "
                       f"{parameter.min_k}-{parameter.max_k}bp")

        # Apply betaine if not explicitly set and recommended
        current_betaine = getattr(parameter, 'betaine_m', 0.0)
        if current_betaine == 0.0 and adaptive_params.betaine_m > 0:
            parameter.betaine_m = adaptive_params.betaine_m
            logger.info(f"GC-adaptive: Setting betaine to {adaptive_params.betaine_m}M")

        # Apply DMSO if not explicitly set and recommended
        current_dmso = getattr(parameter, 'dmso_percent', 0.0)
        if current_dmso == 0.0 and adaptive_params.dmso_percent > 0:
            parameter.dmso_percent = adaptive_params.dmso_percent
            logger.info(f"GC-adaptive: Setting DMSO to {adaptive_params.dmso_percent}%")

        # Log overall strategy
        logger.info(f"GC-adaptive strategy: {adaptive_params.genome_class.value} genome, "
                   f"confidence {adaptive_params.confidence:.0%}")

    except ImportError as e:
        logger.warning(f"Could not import GCAdaptiveStrategy: {e}")
    except Exception as e:
        logger.warning(f"Error applying GC-adaptive defaults: {e}")


def step1():
    """
    Creates files of all k-mers of length 6 to 12 which is located in the path specificed by --kmer_fore and --kmer_b .
    """
    _initialize()  # Lazy initialization
    for prefix in fg_prefixes + bg_prefixes:
        if not os.path.exists(os.path.dirname(prefix)):
            os.makedirs(os.path.dirname(prefix))

    min_k = getattr(parameter, 'min_k', 6)
    max_k = getattr(parameter, 'max_k', 12)

    # Progress tracking for k-mer counting
    total_genomes = len(fg_prefixes) + len(bg_prefixes)
    if total_genomes > 0:
        logger.info(f"Counting k-mers ({min_k}-{max_k}bp) for {total_genomes} genome(s)...")

    logger.info("Running jellyfish for foreground...")
    for i, fg_prefix in enumerate(fg_prefixes):
        genome_name = os.path.basename(fg_genomes[i])
        with progress_context(f"  Foreground {i+1}/{len(fg_prefixes)}: {genome_name}"):
            run_jellyfish(fg_genomes[i], fg_prefix, min_k, max_k)

    if bg_prefixes:
        logger.info("Running jellyfish for background...")
        for i, bg_prefix in enumerate(bg_prefixes):
            genome_name = os.path.basename(bg_genomes[i])
            with progress_context(f"  Background {i+1}/{len(bg_prefixes)}: {genome_name}"):
                run_jellyfish(bg_genomes[i], bg_prefix, min_k, max_k)

    logger.info("Done running jellyfish")

    return fg_seq_lengths, bg_seq_lengths


def step2(all_primers=None, validate_prerequisites=True):
    """
    Filters all candidate primers according to primer design principles (http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html)
    and minimum foreground frequency, maximum background frequency, and maximum Gini index of distances between binding sites.

    Args:
        all_primers: The list of candidate primers to consider. Defaults to all k-mers from 6 to 12 of the target genome if not input.
        validate_prerequisites: If True, validate that Step 1 outputs exist before running.

    Returns:
        filtered_gini_df: Pandas dataframe containing sequences which pass the foreground frequency, background frequency, and Gini index filters.
    """
    _initialize()  # Lazy initialization

    # Validate prerequisites
    if validate_prerequisites:
        min_k = getattr(parameter, 'min_k', 6)
        max_k = getattr(parameter, 'max_k', 12)
        validation = validate_step2_prerequisites(
            parameter.data_dir, fg_prefixes, bg_prefixes, min_k, max_k
        )
        if not validation.valid:
            raise StepPrerequisiteError(2, validation)
    # Report adaptive GC filtering parameters
    genome_gc = getattr(parameter, 'genome_gc', None)
    gc_min = getattr(parameter, 'gc_min', 0.375)
    gc_max = getattr(parameter, 'gc_max', 0.625)
    if genome_gc is not None:
        logger.info(f"Adaptive GC filtering: genome GC={genome_gc:.1%}, primer range={gc_min:.1%}-{gc_max:.1%}")
    else:
        logger.info(f"GC filtering: primer range={gc_min:.1%}-{gc_max:.1%} (default)")

    # Report Tm filtering parameters
    tm_min = getattr(parameter, 'min_tm', 15)
    tm_max = getattr(parameter, 'max_tm', 55)
    reaction_temp = getattr(parameter, 'reaction_temp', 30)
    logger.info(f"Tm filtering: {tm_min}C-{tm_max}C (reaction temp: {reaction_temp}C)")

    kwargs = {
        "fg_prefixes": fg_prefixes,
        "bg_prefixes": bg_prefixes,
        "fg_total_length": sum(fg_seq_lengths),
        "bg_total_length": sum(bg_seq_lengths),
    }
    if all_primers is None:
        # Use stored min_k and max_k values from parameter module (with safe defaults)
        min_k = getattr(parameter, 'min_k', 6)
        max_k = getattr(parameter, 'max_k', 12)
        kmer_lengths = range(min_k, max_k + 1)
        with progress_context("Loading candidate k-mers"):
            all_primers = get_primer_list_from_kmers(fg_prefixes, kmer_lengths=kmer_lengths)
        logger.info(f"Loaded {len(all_primers)} candidate primers")

    with progress_context("Computing foreground/background rates"):
        rate_df = filter_module.get_all_rates(all_primers, **kwargs)
    filtered_rate_df = rate_df[(rate_df["fg_bool"]) & (rate_df["bg_bool"])]
    filtered_rate_df = filtered_rate_df.drop(["fg_bool", "bg_bool"], axis=1)
    logger.info(
        f"Filtered {len(rate_df) - len(filtered_rate_df)} primers based on foreground/background rate"
    )

    if parameter.verbose:
        logger.debug(f"Rate dataframe:\n{rate_df}")

    # Create position files BEFORE Gini calculation (Gini needs these files to exist)
    # Check if position files exist for speedup
    import os
    k = parameter.min_k  # Use first k-mer size for check
    fg_position_file = f"{fg_prefixes[0]}_{k}mer_positions.h5"
    position_files_exist = os.path.exists(fg_position_file)

    if position_files_exist:
        logger.info(f"Reusing existing position files (incremental update only)")

    with progress_context("Creating position files"):
        string_search.get_positions(
            filtered_rate_df["primer"], fg_prefixes, fg_genomes, circular=parameter.fg_circular
        )
        if len(bg_prefixes) > 0 and len(bg_genomes) > 0:
            string_search.get_positions(
                filtered_rate_df["primer"], bg_prefixes, bg_genomes, circular=parameter.bg_circular
            )

    with progress_context("Computing Gini index"):
        gini_df = filter_module.get_gini(fg_prefixes, fg_genomes, fg_seq_lengths, filtered_rate_df, fg_circular)
    logger.info(f"Filtered {len(filtered_rate_df) - len(gini_df)} primers based on Gini index")
    # Calculate ratio with division-by-zero protection
    # When fg_count is 0, set ratio to infinity (primer never binds target = worst case)
    gini_df["ratio"] = gini_df["bg_count"] / gini_df["fg_count"].replace(0, np.nan)
    gini_df["ratio"] = gini_df["ratio"].fillna(float('inf'))
    filtered_gini_df = gini_df.sort_values(by=["ratio"], ascending=False)[: parameter.max_primer]

    filtered_gini_df.to_csv(os.path.join(parameter.data_dir, "step2_df.csv"))
    logger.info(f"Number of remaining primers: {len(filtered_gini_df['primer'])}")
    return filtered_gini_df


# RANK BY RANDOM FOREST
def step3(validate_prerequisites=True):
    """
    Filters primers according to primer efficacy. To adjust the threshold, use option -a or --min_amp_pred.

    Args:
        validate_prerequisites: If True, validate that Step 2 outputs exist before running.

    Returns:
        joined_step3_df: Pandas dataframe of sequences passing step 3.
    """
    _initialize()  # Lazy initialization

    # Validate prerequisites
    if validate_prerequisites:
        validation = validate_step3_prerequisites(parameter.data_dir)
        if not validation.valid:
            raise StepPrerequisiteError(3, validation)

    # Check for sampling parameters in params
    # Default to 10% sampling for 10x speedup with <5% accuracy loss
    sample_rate = getattr(parameter, 'sample_rate', 0.1)  # Default 10% sampling
    disable_sampling = getattr(parameter, 'disable_kmer_sampling', False)

    if not disable_sampling and sample_rate < 1.0:
        min_count = getattr(parameter, 'min_sample_count', 5)
        rf_preprocessing.enable_kmer_sampling(sample_rate=sample_rate, min_count=min_count)
        logger.info(f"K-mer sampling enabled: {sample_rate*100:.0f}% sample rate (10x speedup)")
    else:
        rf_preprocessing.disable_kmer_sampling()

    step2_df = pd.read_csv(os.path.join(parameter.data_dir, "step2_df.csv"))

    primer_list = step2_df["primer"]
    logger.info(f"Scoring {len(primer_list)} primers...")
    fg_scale = sum(fg_seq_lengths) / 6200

    with progress_context("Computing random forest features"):
        df_pred = rf_preprocessing.create_augmented_df(fg_prefixes, primer_list)
        df_pred = rf_preprocessing.scale_delta_Gs(df_pred, on_scale=fg_scale)
        df_pred["molarity"] = 2.5

    with progress_context("Predicting amplification efficacy"):
        results = rf_preprocessing.predict_new_primers(df_pred)
    results.sort_values(by=["on.target.pred"], ascending=[False], inplace=True)

    step3_df = results[results["on.target.pred"] >= parameter.min_amp_pred]

    step2_df = step2_df.set_index("primer")
    step3_df = step3_df.rename({"sequence": "primer"}, axis="columns")
    step3_df = step3_df.set_index("primer")

    joined_step3_df = step3_df.join(step2_df[["ratio", "gini", "fg_count", "bg_count"]], how="left").sort_values(
        by="gini"
    )

    joined_step3_df.to_csv(os.path.join(parameter.data_dir, "step3_df.csv"))

    logger.info(f"Filtered {step2_df.shape[0] - joined_step3_df.shape[0]} primers based on efficacy")

    if parameter.verbose:
        logger.debug(f"Step 3 results:\n{joined_step3_df}")

    return joined_step3_df


# optimize AND SEARCH
def step4(primer_list=None, scores=None, initial_primer_sets=None, validate_prerequisites=True):
    """
    DEPRECATED: This function is no longer supported.

    Use ImprovedPipeline from neoswga.core.improved_pipeline instead:

        from neoswga.core.improved_pipeline import ImprovedPipeline, PipelineConfig

        config = PipelineConfig(optimization_method='hybrid')
        pipeline = ImprovedPipeline(config=config)
        results = pipeline.run_step4(primer_df, fg_prefixes, bg_prefixes, ...)

    The legacy optimize module has been removed. The improved_pipeline module
    provides superior optimization algorithms including:
    - hybrid: Combined network and greedy optimization (default)
    - dominating-set: Fast graph-based coverage optimization
    - background-aware: Clinical-grade background minimization
    - genetic: Multi-objective genetic algorithm
    - moea: Pareto optimization

    Raises:
        NotImplementedError: Always raised, use ImprovedPipeline instead.
    """
    raise NotImplementedError(
        "step4() is deprecated. Use ImprovedPipeline from "
        "neoswga.core.improved_pipeline instead. See docstring for details."
    )


