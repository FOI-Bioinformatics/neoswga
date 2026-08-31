"""Shared CLI helpers extracted from cli_unified.py.

Pure leaf utilities (input validation, params loading, GPU setup, reaction
presets, RNG seeding, run-manifest writing). Handler submodules and
cli_unified both import from here, so this module must not import either of
them (no circular dependency).
"""

import argparse
import functools
import json
import logging
import os
import sys

from neoswga.core.registry import polymerase_names as _polymerase_names

logger = logging.getLogger(__name__)


def params_command(_func=None, *, merge=("json_file",), seed=False, validate=True):
    """Decorator for params.json-based CLI handlers' common entry preamble.

    Runs, in order, before the wrapped handler body:
      1. ``validate_params_json_file(args.json_file)`` -- when ``validate`` and
         ``args`` carries a non-None ``json_file``;
      2. ``merge_args_to_parameter(args, parameter, merge)`` -- when ``merge`` is
         truthy (default: just the standard ``["json_file"]`` merge);
      3. ``_apply_seed(args)`` -- when ``seed`` is True.

    Handler-specific extras (additional ``merge_args_to_parameter`` calls,
    ``_record_run_manifest`` writes, timing, try/except) stay in the body; this
    only factors out the mechanical, zero-variation opener. Validating before
    seeding means a bad params path fails fast (the previous hand-written order
    seeded first, which was wasted work on the error path).
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(args):
            if validate and getattr(args, "json_file", None) is not None:
                validate_params_json_file(args.json_file)
            if merge:
                from neoswga.core import parameter

                merge_args_to_parameter(args, parameter, list(merge))
            if seed:
                _apply_seed(args)
            return func(args)

        return wrapper

    return decorator(_func) if _func is not None else decorator


def check_jellyfish_available():
    """Check that Jellyfish is installed and raise a clear error if not.

    Called before commands that require Jellyfish (count-kmers, design).
    """
    from neoswga.core.kmer_counter import check_jellyfish_available as _check
    from neoswga.core.kmer_counter import get_jellyfish_version

    if not _check():
        logger.error(
            "Jellyfish is required but not found in PATH.\n"
            "Jellyfish is an external k-mer counting tool that must be installed separately.\n\n"
            "Install via conda:  conda install -c bioconda jellyfish\n"
            "Install via brew:   brew install jellyfish\n"
            "Build from source:  https://github.com/gmarcais/Jellyfish\n\n"
            "After installation, verify with:  jellyfish --version"
        )
        sys.exit(1)
    version = get_jellyfish_version()
    if version and version.split(".")[0] == "1":
        logger.error(
            f"Jellyfish version {version} detected. NeoSWGA requires Jellyfish 2.x.\n"
            "Please upgrade: https://github.com/gmarcais/Jellyfish"
        )
        sys.exit(1)


# =============================================================================
# Input Validation Utilities
# =============================================================================

# =============================================================================
# Validation Constants
# =============================================================================

# Valid DNA bases (IUPAC nucleotide codes)
VALID_DNA_BASES = frozenset("ACGTRYSWKMBDHVN")

# Primer length bounds
MAX_PRIMER_LENGTH = 50
MIN_PRIMER_LENGTH = 4

# Forbidden characters in file paths (shell metacharacters, null bytes)
FORBIDDEN_PATH_CHARS = frozenset(";&|`$\x00")


def validate_primer_sequence(primer: str, name: str = "primer") -> str:
    """
    Validate a primer sequence for safety and correctness.

    Args:
        primer: The primer sequence to validate
        name: Name for error messages (e.g., "fixed primer")

    Returns:
        Validated and normalized (uppercase) primer sequence

    Raises:
        ValueError: If primer is invalid
    """
    if not primer:
        raise ValueError(f"Empty {name} sequence")

    # Normalize to uppercase
    primer = primer.strip().upper()

    # Check length
    if len(primer) < MIN_PRIMER_LENGTH:
        raise ValueError(
            f"Invalid {name} '{primer}': too short " f"(minimum {MIN_PRIMER_LENGTH} bp)"
        )
    if len(primer) > MAX_PRIMER_LENGTH:
        raise ValueError(
            f"Invalid {name} '{primer}': too long " f"(maximum {MAX_PRIMER_LENGTH} bp)"
        )

    # Check for valid DNA bases only
    invalid_chars = set(primer) - VALID_DNA_BASES
    if invalid_chars:
        raise ValueError(
            f"Invalid {name} '{primer}': contains invalid characters "
            f"{invalid_chars}. Only IUPAC nucleotide codes allowed."
        )

    return primer


def validate_path_security(path: str, context: str = "path") -> None:
    """
    Validate a path for security issues (traversal, injection).

    This is the shared validation logic used by both primer file
    and log file path validation.

    Args:
        path: The path string to validate (before normalization)
        context: Description for error messages (e.g., "primer file", "log")

    Raises:
        ValueError: If path contains security issues
    """
    # Check for path traversal BEFORE normalization (normpath resolves '..')
    if ".." in path:
        raise ValueError(f"Invalid {context}: directory traversal not allowed")

    # Check for forbidden characters (shell metacharacters, null bytes)
    if FORBIDDEN_PATH_CHARS.intersection(path):
        raise ValueError(f"Invalid {context}: contains suspicious characters")


def validate_primer_file(filepath: str) -> str:
    """
    Validate a primer file path for safety.

    Args:
        filepath: Path to validate

    Returns:
        Validated absolute path

    Raises:
        ValueError: If path is invalid or suspicious
    """
    validate_path_security(filepath, "primer file path")

    # Now safe to normalize
    filepath = os.path.expanduser(filepath)
    filepath = os.path.normpath(filepath)
    filepath = os.path.abspath(filepath)

    return filepath


def load_primers_from_file(filepath: str, name: str = "primer") -> list:
    """
    Safely load and validate primers from a file.

    Args:
        filepath: Path to primer file (one primer per line)
        name: Name for error messages

    Returns:
        List of validated primer sequences
    """
    filepath = validate_primer_file(filepath)

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Primer file not found: {filepath}")

    primers = []
    with open(filepath) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue
            try:
                primer = validate_primer_sequence(line, f"{name} (line {line_num})")
                if primer not in primers:
                    primers.append(primer)
            except ValueError as e:
                logger.warning(f"Skipping invalid line {line_num}: {e}")

    return primers


def collect_primers_from_args(
    cli_primers: list,
    primers_file: str,
    name: str = "primer",
    allow_empty: bool = False,
) -> list:
    """
    Collect and validate primers from CLI arguments and/or file.

    Consolidates the repeated pattern of loading primers from both
    command-line arguments and files, with deduplication.

    Args:
        cli_primers: List of primer sequences from CLI arguments (or None)
        primers_file: Path to primer file (or None)
        name: Name for error messages (e.g., "fixed primer")
        allow_empty: If False, raises error when no primers found

    Returns:
        List of validated, deduplicated primer sequences

    Raises:
        SystemExit: On validation errors (logs error and exits)
    """
    primers = []

    # Load from CLI arguments
    if cli_primers:
        for primer in cli_primers:
            try:
                validated = validate_primer_sequence(primer, name)
                if validated not in primers:
                    primers.append(validated)
            except ValueError as e:
                logger.error(str(e))
                sys.exit(1)

    # Load from file
    if primers_file:
        try:
            file_primers = load_primers_from_file(primers_file, name)
            for primer in file_primers:
                if primer not in primers:
                    primers.append(primer)
        except (ValueError, FileNotFoundError) as e:
            logger.error(str(e))
            sys.exit(1)

    # Check if we got any primers
    if not allow_empty and not primers:
        logger.error(f"No {name}s provided")
        sys.exit(1)

    return primers


# =============================================================================
# Parameter File Validation
# =============================================================================


def validate_params_json_file(path):
    """Validate that a parameter file exists and contains valid JSON.

    Prints a user-friendly error and exits if the file is missing or malformed.
    """
    if path is None:
        print("Error: -j/--json-file is required.", file=sys.stderr)
        print("Usage: neoswga <command> -j params.json", file=sys.stderr)
        print("To create a params.json: neoswga init --genome target.fasta", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(path):
        print(f"Error: parameter file '{path}' not found.", file=sys.stderr)
        sys.exit(1)
    try:
        with open(path) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: '{path}' is not valid JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Validate that essential fields are present
    if not data:
        print(
            f"Error: '{path}' is empty. At minimum, provide 'fg_genome' and 'data_dir'.",
            file=sys.stderr,
        )
        print(
            "Run 'neoswga init --genome target.fasta' to create a valid configuration.",
            file=sys.stderr,
        )
        sys.exit(1)
    # Check for data_dir
    if "data_dir" not in data:
        print(f"Error: '{path}' is missing required field 'data_dir'.", file=sys.stderr)
        print(
            "Run 'neoswga init --genome target.fasta' to create a valid configuration.",
            file=sys.stderr,
        )
        sys.exit(1)
    # Validate data_dir exists (create if needed)
    data_dir = data["data_dir"]
    if data_dir and data_dir != "./" and data_dir != ".":
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir, exist_ok=True)
            logger.info(f"Created data directory: {data_dir}")
    # Check for fg_prefixes (common missing field — #1 new-user crash)
    if "fg_genomes" in data and "fg_prefixes" not in data:
        logger.warning(
            f"'{path}' is missing 'fg_prefixes'. "
            "Consider running 'neoswga init --genome target.fasta' for a complete config."
        )
    # Check for genome specification (multiple accepted conventions)
    has_genome = any(k in data for k in ("fg_genome", "fg_genomes", "fg_prefixes"))
    if not has_genome:
        print(f"Error: '{path}' is missing genome specification.", file=sys.stderr)
        print(
            "Provide 'fg_genome', 'fg_genomes', or 'fg_prefixes' in your params.json.",
            file=sys.stderr,
        )
        print(
            "Run 'neoswga init --genome target.fasta' to create a valid configuration.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Run automatic parameter validation (ERROR-level issues only)
    try:
        from neoswga.core.param_validator import ParamValidator, ValidationLevel

        validator = ParamValidator()
        results = validator.validate_params(data)
        # Skip errors already handled by the checks above (required keys, file existence)
        # to avoid false positives — only catch value/type errors (ranges, polymerase, method)
        already_checked = {"fg_genomes", "bg_genomes", "fg_prefixes", "data_dir", "params_file"}
        errors = [
            r
            for r in results
            if r.level == ValidationLevel.ERROR
            and r.parameter not in already_checked
            and "File not found" not in r.message
        ]
        if errors:
            print(
                f"\nParameter validation found {len(errors)} error(s) in '{path}':", file=sys.stderr
            )
            for err in errors:
                print(f"  - {err.message}", file=sys.stderr)
            print(f"\nRun 'neoswga validate-params -j {path}' for full details.", file=sys.stderr)
            sys.exit(1)
    except ImportError:
        pass  # ParamValidator not available — skip validation
    except Exception as e:
        # Don't block pipeline on validator bugs
        logger.debug(f"Parameter validation skipped: {e}")


# =============================================================================
# Parameter Merger Utility
# =============================================================================


def merge_args_to_parameter(args, parameter, param_names, mapping=None):
    """
    Merge CLI arguments to the parameter module.

    Reduces repetitive hasattr/setattr patterns by handling parameter
    transfer in a generic way.

    Args:
        args: argparse Namespace with CLI arguments
        parameter: The parameter module to update
        param_names: List of parameter names to transfer
        mapping: Optional dict mapping arg names to parameter names
                 (e.g., {'bsa': 'bsa_ug_ml'})

    Example:
        merge_args_to_parameter(args, parameter, [
            'gc_min', 'gc_max', 'min_tm', 'max_tm',
            'min_fg_freq', 'max_bg_freq', 'max_gini'
        ])

        # With mapping
        merge_args_to_parameter(args, parameter, ['bsa'], {'bsa': 'bsa_ug_ml'})
    """
    if mapping is None:
        mapping = {}

    for param in param_names:
        if hasattr(args, param):
            value = getattr(args, param)
            if value is not None:
                target_name = mapping.get(param, param)
                setattr(parameter, target_name, value)


UNIMPLEMENTED_MESSAGE = (
    "{flag} is accepted but not implemented: {detail}. "
    "The run continues with the option ignored."
)

# Options the parsers accept that no module consumes. They stay in the
# interface because removing a published argument turns a working script into a
# parse error, which is the repository owner's call; what they no longer do is
# log a message that reads as confirmation. An entry here is a promise the tool
# does not keep, so wiring an option means deleting its entry, not editing it.
UNIMPLEMENTED_OPTIONS = {
    "use_enhanced_features": (
        "--use-enhanced-features",
        "the scoring step runs the standard random-forest model in "
        "neoswga/core/models/, and rf_preprocessing.predict_with_enhanced_"
        "features has no caller in the pipeline",
    ),
    "enhanced_model_path": (
        "--enhanced-model-path",
        "no scoring path loads an enhanced model file",
    ),
    "use_background_filter": (
        "--use-background-filter",
        "optimize builds no Bloom filter; candidates are screened against the "
        "background by the separate pre-filter that --no-bg-prefilter disables",
    ),
    "background_bloom_path": (
        "--background-bloom-path",
        "a pre-built Bloom filter is read only by improved_pipeline, which the "
        "optimize command does not use; see --no-bg-prefilter for the screening "
        "optimize does perform",
    ),
    "background_sampled_path": (
        "--background-sampled-path",
        "a pre-built sampled index is read only by improved_pipeline, which the "
        "optimize command does not use; see --no-bg-prefilter for the screening "
        "optimize does perform",
    ),
    "max_optimization_time": (
        "--max-optimization-time",
        "no optimizer enforces a wall-clock budget, so a long run is not " "interrupted",
    ),
}


def report_unimplemented_options(args):
    """Warn about each unconsumed option the caller actually set.

    Returns the flags reported, so a caller can act on the list rather than
    parse the log. An option left at its default is silent: only a request the
    tool cannot honour is worth interrupting the user for.
    """
    reported = []

    for attr, (flag, detail) in UNIMPLEMENTED_OPTIONS.items():
        value = getattr(args, attr, None)
        if value is None or value is False:
            continue
        logger.warning(UNIMPLEMENTED_MESSAGE.format(flag=flag, detail=detail))
        reported.append(flag)

    return reported


def setup_gpu_acceleration(args, parameter, quiet=False):
    """Record the GPU request and report that it does not change the run.

    `neoswga.core.gpu_acceleration` exists, but no pipeline stage dispatches to
    it: neither `--use-gpu` nor CuPy being importable changes a result or a
    runtime. The attributes are still written because they are part of the
    parameter object's surface and cost nothing to keep.

    `quiet` is accepted for call-site compatibility and does not suppress the
    warning -- an option that cannot do what it says is not routine output.
    """
    no_gpu = getattr(args, "no_gpu", False)
    use_gpu = getattr(args, "use_gpu", False)

    if no_gpu:
        parameter.use_gpu = False
        return

    if use_gpu:
        parameter.use_gpu = True
        parameter.gpu_device = getattr(args, "gpu_device", 0)
        logger.warning(
            UNIMPLEMENTED_MESSAGE.format(
                flag="--use-gpu",
                detail=(
                    "no pipeline stage dispatches to neoswga.core.gpu_acceleration, "
                    "so results and runtime are unchanged whether or not CuPy is "
                    "installed"
                ),
            )
        )
        return

    parameter.use_gpu = False
    logger.debug("GPU acceleration is not connected to any pipeline stage")


# Preset configurations for reaction conditions
# The named presets, each an optimizer-facing config derived from the SAME
# reaction-conditions factory that `--preset` applies (see PRESET_FACTORIES and
# load_preset_conditions below). It used to be a hand-maintained literal that had
# drifted: `show-presets` printed 37 C / 0 mM Mg for standard_phi29 while
# `--preset standard_phi29` applied 30 C / 2.5 mM, and long_primers_15mer
# advertised 7% DMSO / 1.5 M betaine while applying 5% / 1.0 M. Three of four
# presets described something the tool would not do. Deriving them removes the
# divergence by construction.
_PRESET_OPTIMIZER_HINTS = {
    "standard_phi29": {"optimization_method": "hybrid", "target_set_size": 6},
    "enhanced_equiphi29": {"optimization_method": "hybrid", "target_set_size": 6},
    "long_primers_15mer": {"optimization_method": "hybrid", "target_set_size": 6},
    "high_gc_genome": {"optimization_method": "hybrid", "target_set_size": 8},
}


def _build_presets():
    out = {}
    for name, hints in _PRESET_OPTIMIZER_HINTS.items():
        applied = load_preset_conditions(name)
        out[name] = {
            "temperature": applied["reaction_temp"],
            "polymerase": applied["polymerase"],
            "dmso_percent": applied["dmso_percent"],
            "betaine_m": applied["betaine_m"],
            "na_conc": applied["na_conc"],
            "mg_conc": applied["mg_conc"],
            "ssb": applied["ssb"],
            **hints,
        }
    return out


# Command groups for organized help display


# Concentration additives carried by a reaction-conditions preset. Kept as an
# explicit tuple rather than introspected so the set is greppable and reviewable;
# tests/test_preset_round_trip.py asserts it stays complete against
# ReactionConditions.__init__.
PRESET_ADDITIVE_FIELDS = (
    "dmso_percent",
    "betaine_m",
    "trehalose_m",
    "formamide_percent",
    "glycerol_percent",
    "bsa_ug_ml",
    "peg_percent",
    "ethanol_percent",
    "urea_m",
    "tmac_m",
    "propanediol_m",
    # Buffer species added with the schema v2 chemistry work. Defaults of 0
    # preserve behaviour; presets may set them to describe a real buffer.
    "k_conc",
    "nh4_conc",
    "dntp_conc",
    "dtt_mm",
)


def load_preset_conditions(preset_name):
    """
    Load predefined reaction conditions preset.

    Args:
        preset_name: Name of preset to load

    Returns:
        Dictionary of reaction condition parameters
    """
    from neoswga.core import reaction_conditions as rc

    preset_map = {
        "standard_phi29": rc.get_standard_conditions(),
        "enhanced_equiphi29": rc.get_enhanced_conditions(),
        "high_gc_genome": rc.get_high_gc_conditions(),
        "long_primers_15mer": rc.get_enhanced_conditions(),  # Same as enhanced but user intent is longer primers
        "q_solution": rc.get_q_solution_equivalent(),
        "gc_melt": rc.get_gc_melt_conditions(),
        "crude_sample": rc.get_crude_sample_conditions(),
        "low_temp": rc.get_low_temp_conditions(),
        "bst": rc.get_bst_conditions(),
        "klenow": rc.get_klenow_conditions(),
        "extreme_gc": rc.get_extreme_gc_conditions(),
    }

    if preset_name not in preset_map:
        logger.error(f"Unknown preset: {preset_name}")
        return {}

    conditions = preset_map[preset_name]

    # Convert ReactionConditions to dict.
    #
    # Every concentration knob on ReactionConditions must appear here. This list
    # previously omitted urea_m, tmac_m, ethanol_percent and formamide_percent, so
    # `--preset extreme_gc` silently dropped the two additives that define it
    # (urea_m=0.5, tmac_m=0.05) and applied a preset that was no longer that preset.
    # tests/test_preset_round_trip.py enforces coverage against the dataclass.
    out = {
        "reaction_temp": conditions.temp,
        "polymerase": conditions.polymerase,
        "na_conc": conditions.na_conc,
        "mg_conc": conditions.mg_conc,
        "ssb": conditions.ssb,
    }
    for name in PRESET_ADDITIVE_FIELDS:
        out[name] = getattr(conditions, name)
    return out


# Built after load_preset_conditions so the two cannot disagree.
PRESETS = _build_presets()


def add_common_options(parser):
    """Add options common to all steps"""
    parser.add_argument("-j", "--json-file", type=str, help="Parameters JSON file")
    parser.add_argument("-z", "--data-dir", type=str, help="Data directory")
    parser.add_argument(
        "--polymerase",
        choices=_polymerase_names(),
        help="Polymerase type: phi29 (30-40C), equiphi29 (42-45C), bst (60-65C), klenow (25-40C)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("-q", "--quiet", action="store_true", help="Minimal output")

    # GPU acceleration. Accepted, and not implemented: see
    # UNIMPLEMENTED_OPTIONS and setup_gpu_acceleration.
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="Not implemented: no pipeline stage dispatches to the CuPy "
        "backend, so this changes neither results nor runtime",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Accepted for compatibility. No stage uses a GPU, so this is " "already the behaviour",
    )
    parser.add_argument(
        "--gpu-device",
        type=int,
        default=0,
        help="GPU device ID. Recorded and not used; see --use-gpu",
    )

    # Quality assurance (Category 3)
    parser.add_argument(
        "--enable-qa", action="store_true", help="Enable quality assurance checks at each step"
    )


def _apply_seed(args):
    """Seed all RNGs from --seed for the post-optimize set-producing commands.

    The `optimize` step seeds inside unified_optimizer.run_optimization, but
    expand-primers / swap-primer / contract-set / multi-genome run their
    greedy/stochastic logic directly, so seed here for reproducibility.
    """
    seed = getattr(args, "seed", None)
    if seed is None:
        return
    import random

    import numpy as np

    random.seed(seed)
    np.random.seed(seed)
    try:
        from neoswga.core.rf_preprocessing import set_kmer_sampling_seed

        set_kmer_sampling_seed(seed)
    except Exception as e:  # pragma: no cover - defensive
        logger.debug(f"k-mer sampling seed skipped: {e}")
    logger.info(f"Random seed set to {seed} for reproducibility")


def _effective_conditions(parameter):
    """The reaction the run is actually using, read off the parameter module.

    params.json is not the answer: `retune_for_polymerase` and the GC-adaptive
    strategy both write back here, so a GC-rich target designed with 2 M
    betaine leaves no trace in the file the user wrote. Snapshotting the
    resolved state is what lets `export` and `report` correct a Tm for the
    buffer the design was optimized under instead of a different one.
    """
    from neoswga.core.parameter import REACTION_PARAM_DEFAULTS

    fields = ("polymerase", "reaction_temp", "na_conc", "mg_conc", *REACTION_PARAM_DEFAULTS)
    conditions = {}
    for name in fields:
        value = getattr(parameter, name, None)
        if value is not None:
            conditions[name] = value
    return conditions or None


def _record_run_manifest(step: str, args, parameter, input_files=None):
    """Best-effort wrapper around run_manifest.write_manifest.

    Failures are swallowed so manifest issues never break a pipeline that
    otherwise succeeded.
    """
    try:
        from neoswga.core.run_manifest import write_manifest

        write_manifest(
            effective_conditions=_effective_conditions(parameter),
            step=step,
            data_dir=getattr(parameter, "data_dir", None),
            params_path=getattr(args, "json_file", None),
            input_files=input_files,
            # The CLI seed lives on args (--seed), not parameter; the previous
            # getattr(parameter, "seed") recorded None even when --seed was set.
            seed=getattr(args, "seed", None) or getattr(parameter, "seed", None),
        )
    except Exception as e:
        logger.debug(f"run_manifest write skipped: {e}")


def bootstrap_params_from_genome(genome_paths, output_dir, circular=True):
    """Populate the parameter module from bare genome FASTAs.

    The iterative commands all require a params.json carrying `fg_prefixes`,
    `fg_genomes` and `fg_seq_lengths`, which in practice means a prior
    `neoswga init` plus `count-kmers`. Someone holding an oligo list and a
    reference genome has none of that, and the missing piece is bookkeeping
    rather than anything scientific -- so synthesise it.

    K-mer index files are NOT created here. Commands that need binding
    positions should pass `genome_paths` with `on_missing='scan'` so
    PositionCache finds the primers by scanning the FASTA directly.

    Args:
        genome_paths: One or more foreground FASTA paths.
        output_dir: Directory used for `data_dir` and as the k-mer prefix root.
        circular: Whether the targets are circular.

    Returns:
        The dict of parameters that were applied.
    """
    import os

    from neoswga.core import parameter
    from neoswga.core import utility as _utility

    genome_paths = [genome_paths] if isinstance(genome_paths, str) else list(genome_paths)
    for path in genome_paths:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Genome file not found: {path}")

    os.makedirs(output_dir, exist_ok=True)
    prefixes = [
        os.path.join(output_dir, os.path.splitext(os.path.basename(g))[0]) for g in genome_paths
    ]
    lengths = _utility.get_all_seq_lengths(
        fname_genomes=genome_paths, cpus=getattr(parameter, "cpus", 1) or 1
    )

    applied = {
        "fg_genomes": genome_paths,
        "fg_prefixes": prefixes,
        "fg_seq_lengths": lengths,
        "data_dir": output_dir,
        "src_dir": output_dir,
        "fg_circular": circular,
    }
    for key, value in applied.items():
        setattr(parameter, key, value)

    logger.info(
        "Bootstrapped parameters from %d genome(s): %s bp total",
        len(genome_paths),
        f"{sum(lengths):,}",
    )
    return applied


def add_position_source_options(parser):
    """Add the flags that let a command find binding positions without an index.

    Shared so `--genome` means the same thing everywhere: scan this FASTA for
    any primer that is not already in the HDF5 position files, instead of
    silently reporting zero coverage for it.
    """
    parser.add_argument(
        "--genome",
        nargs="+",
        metavar="FASTA",
        help="Target genome FASTA(s). Primers missing from the k-mer index are "
        "found by scanning these directly, so an externally-designed oligo set "
        "can be evaluated without a prior count-kmers run. Also lets -j be "
        "omitted entirely.",
    )
    parser.add_argument(
        "--scan-background",
        action="store_true",
        help="Also scan background genomes for missing primers. Off by default: "
        "background specificity is a frequency question the k-mer counts already "
        "answer, and scanning a host-sized genome holds it in memory.",
    )
    return parser


def set_size_shortfall_advice(num_found, target_size, method=None):
    """What to say when fewer primers came back than were asked for.

    This check compares two integers. It cannot see why they differ, and it
    used to assert one: "PARTIAL result: insufficient candidates", followed by
    four remedies all about enlarging the pool.

    At least one common cause is not that, and the advice points away from it.
    Stage 1 stops once coverage is satisfied, so a set-cover that reaches the
    coverage target with fewer primers than requested returns early -- the
    optimizer succeeding, reported as a partial failure. Loosening filters in
    response degrades a design that was not deficient. The same misdirection is
    on record in docs/validation/additive_specificity.md, where the advice was
    followed from a different root cause.

    So the shortfall is stated, the benign explanation is offered first, and
    the pool remedies are kept but demoted to what they are: possibilities.
    """
    lines = [
        f"Found {num_found} primers but the target was {target_size}.",
        "This is not necessarily a shortfall. Stage 1 stops once the coverage "
        "target is met, so a set that covers the genome with fewer primers "
        "returns early. Check the reported coverage before treating it as one.",
        "If coverage is also below target, the pool may genuinely be too small:",
        "  - Relaxing filter thresholds (max_bg_freq, max_gini)",
        "  - Widening k-mer range (min_k / max_k)",
        "  - Increasing candidate pool (max_primer)",
    ]
    if method:
        lines.append(f"  - Trying a different optimizer (current: {method})")
    return lines
