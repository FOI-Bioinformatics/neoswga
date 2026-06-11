"""Shared CLI helpers extracted from cli_unified.py.

Pure leaf utilities (input validation, params loading, GPU setup, reaction
presets, RNG seeding, run-manifest writing). Handler submodules and
cli_unified both import from here, so this module must not import either of
them (no circular dependency).
"""

import argparse
import json
import logging
import os
import sys

logger = logging.getLogger(__name__)


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


def setup_gpu_acceleration(args, parameter, quiet=False):
    """
    Configure GPU acceleration with auto-detection.

    Auto-enables GPU if CuPy is available unless --no-gpu is specified.
    Can be explicitly enabled with --use-gpu.
    """
    # Check if user explicitly disabled GPU
    no_gpu = getattr(args, "no_gpu", False)
    use_gpu = getattr(args, "use_gpu", False)

    if no_gpu:
        parameter.use_gpu = False
        return

    if use_gpu:
        # User explicitly requested GPU
        parameter.use_gpu = True
        parameter.gpu_device = getattr(args, "gpu_device", 0)
        if not quiet:
            logger.info(f"GPU acceleration enabled (device {parameter.gpu_device})")
        return

    # Auto-detect GPU availability
    try:
        from neoswga.core.gpu_acceleration import is_gpu_available

        if is_gpu_available():
            parameter.use_gpu = True
            parameter.gpu_device = getattr(args, "gpu_device", 0)
            if not quiet:
                logger.info(f"GPU auto-detected and enabled (device {parameter.gpu_device})")
    except ImportError:
        pass  # GPU module not available


# Preset configurations for reaction conditions
PRESETS = {
    "standard_phi29": {
        "temperature": 37.0,
        "dmso_percent": 0.0,
        "betaine_m": 0.0,
        "polymerase": "phi29",
        "na_conc": 50.0,
        "mg_conc": 0.0,
        "ssb": False,
        "optimization_method": "hybrid",
        "target_set_size": 6,
    },
    "enhanced_equiphi29": {
        "temperature": 42.0,
        "dmso_percent": 5.0,
        "betaine_m": 1.0,
        "polymerase": "equiphi29",
        "na_conc": 50.0,
        "mg_conc": 0.0,
        "ssb": True,
        "optimization_method": "hybrid",
        "target_set_size": 6,
    },
    "long_primers_15mer": {
        "temperature": 45.0,
        "dmso_percent": 7.0,
        "betaine_m": 1.5,
        "polymerase": "equiphi29",
        "na_conc": 50.0,
        "mg_conc": 0.0,
        "ssb": True,
        "optimization_method": "hybrid",
        "target_set_size": 6,
    },
    "high_gc_genome": {
        "temperature": 45.0,
        "dmso_percent": 10.0,
        "betaine_m": 2.0,
        "polymerase": "equiphi29",
        "na_conc": 50.0,
        "mg_conc": 0.0,
        "ssb": True,
        "optimization_method": "hybrid",
        "target_set_size": 8,
    },
}


# Command groups for organized help display


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

    # Convert ReactionConditions to dict
    return {
        "reaction_temp": conditions.temp,
        "polymerase": conditions.polymerase,
        "dmso_percent": conditions.dmso_percent,
        "betaine_m": conditions.betaine_m,
        "trehalose_m": conditions.trehalose_m,
        "glycerol_percent": conditions.glycerol_percent,
        "bsa_ug_ml": conditions.bsa_ug_ml,
        "peg_percent": conditions.peg_percent,
        "na_conc": conditions.na_conc,
        "mg_conc": conditions.mg_conc,
        "ssb": conditions.ssb,
    }


def add_common_options(parser):
    """Add options common to all steps"""
    parser.add_argument("-j", "--json-file", type=str, help="Parameters JSON file")
    parser.add_argument("-z", "--data-dir", type=str, help="Data directory")
    parser.add_argument(
        "--polymerase",
        choices=["phi29", "equiphi29", "bst", "klenow"],
        help="Polymerase type: phi29 (30-40C), equiphi29 (42-45C), bst (60-65C), klenow (25-40C)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("-q", "--quiet", action="store_true", help="Minimal output")

    # GPU acceleration (Category 1 orphaned feature)
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="Use GPU acceleration (requires CuPy). Provides 10-100x speedup for large genomes",
    )
    parser.add_argument(
        "--no-gpu",
        action="store_true",
        help="Disable auto-GPU detection (GPU is auto-enabled when CuPy is available)",
    )
    parser.add_argument(
        "--gpu-device", type=int, default=0, help="GPU device ID to use (default: 0)"
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


def _record_run_manifest(step: str, args, parameter, input_files=None):
    """Best-effort wrapper around run_manifest.write_manifest.

    Failures are swallowed so manifest issues never break a pipeline that
    otherwise succeeded.
    """
    try:
        from neoswga.core.run_manifest import write_manifest

        write_manifest(
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
