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
    neoswga optimize [options] [--optimization-method=hybrid|greedy|milp]

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

    # Control optimization method
    neoswga optimize -j params.json --optimization-method=milp
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

import sys
import os
import argparse
import logging
import json

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
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)


def check_jellyfish_available():
    """Check that Jellyfish is installed and raise a clear error if not.

    Called before commands that require Jellyfish (count-kmers, design).
    """
    from neoswga.core.kmer_counter import check_jellyfish_available as _check, get_jellyfish_version
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
    if version and version.split('.')[0] == '1':
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
VALID_DNA_BASES = frozenset('ACGTRYSWKMBDHVN')

# Primer length bounds
MAX_PRIMER_LENGTH = 50
MIN_PRIMER_LENGTH = 4

# Forbidden characters in file paths (shell metacharacters, null bytes)
FORBIDDEN_PATH_CHARS = frozenset(';&|`$\x00')


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
            f"Invalid {name} '{primer}': too short "
            f"(minimum {MIN_PRIMER_LENGTH} bp)"
        )
    if len(primer) > MAX_PRIMER_LENGTH:
        raise ValueError(
            f"Invalid {name} '{primer}': too long "
            f"(maximum {MAX_PRIMER_LENGTH} bp)"
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
    if '..' in path:
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
            if not line or line.startswith('#'):
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
        print(f"Error: '{path}' is empty. At minimum, provide 'fg_genome' and 'data_dir'.", file=sys.stderr)
        print("Run 'neoswga init --genome target.fasta' to create a valid configuration.", file=sys.stderr)
        sys.exit(1)
    # Check for data_dir
    if 'data_dir' not in data:
        print(f"Error: '{path}' is missing required field 'data_dir'.", file=sys.stderr)
        print("Run 'neoswga init --genome target.fasta' to create a valid configuration.", file=sys.stderr)
        sys.exit(1)
    # Validate data_dir exists (create if needed)
    data_dir = data['data_dir']
    if data_dir and data_dir != './' and data_dir != '.':
        if not os.path.isdir(data_dir):
            os.makedirs(data_dir, exist_ok=True)
            logger.info(f"Created data directory: {data_dir}")
    # Check for fg_prefixes (common missing field — #1 new-user crash)
    if 'fg_genomes' in data and 'fg_prefixes' not in data:
        logger.warning(
            f"'{path}' is missing 'fg_prefixes'. "
            "Consider running 'neoswga init --genome target.fasta' for a complete config."
        )
    # Check for genome specification (multiple accepted conventions)
    has_genome = any(k in data for k in ('fg_genome', 'fg_genomes', 'fg_prefixes'))
    if not has_genome:
        print(f"Error: '{path}' is missing genome specification.", file=sys.stderr)
        print("Provide 'fg_genome', 'fg_genomes', or 'fg_prefixes' in your params.json.", file=sys.stderr)
        print("Run 'neoswga init --genome target.fasta' to create a valid configuration.", file=sys.stderr)
        sys.exit(1)

    # Run automatic parameter validation (ERROR-level issues only)
    try:
        from neoswga.core.param_validator import ParamValidator, ValidationLevel
        validator = ParamValidator()
        results = validator.validate_params(data)
        # Skip errors already handled by the checks above (required keys, file existence)
        # to avoid false positives — only catch value/type errors (ranges, polymerase, method)
        already_checked = {'fg_genomes', 'bg_genomes', 'fg_prefixes', 'data_dir', 'params_file'}
        errors = [
            r for r in results
            if r.level == ValidationLevel.ERROR
            and r.parameter not in already_checked
            and 'File not found' not in r.message
        ]
        if errors:
            print(f"\nParameter validation found {len(errors)} error(s) in '{path}':",
                  file=sys.stderr)
            for err in errors:
                print(f"  - {err.message}", file=sys.stderr)
            print(f"\nRun 'neoswga validate-params -j {path}' for full details.",
                  file=sys.stderr)
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
    no_gpu = getattr(args, 'no_gpu', False)
    use_gpu = getattr(args, 'use_gpu', False)

    if no_gpu:
        parameter.use_gpu = False
        return

    if use_gpu:
        # User explicitly requested GPU
        parameter.use_gpu = True
        parameter.gpu_device = getattr(args, 'gpu_device', 0)
        if not quiet:
            logger.info(f"GPU acceleration enabled (device {parameter.gpu_device})")
        return

    # Auto-detect GPU availability
    try:
        from neoswga.core.gpu_acceleration import is_gpu_available
        if is_gpu_available():
            parameter.use_gpu = True
            parameter.gpu_device = getattr(args, 'gpu_device', 0)
            if not quiet:
                logger.info(f"GPU auto-detected and enabled (device {parameter.gpu_device})")
    except ImportError:
        pass  # GPU module not available


# Preset configurations for reaction conditions
PRESETS = {
    'standard_phi29': {
        'temperature': 37.0,
        'dmso_percent': 0.0,
        'betaine_m': 0.0,
        'polymerase': 'phi29',
        'na_conc': 50.0,
        'mg_conc': 0.0,
        'ssb': False,
        'optimization_method': 'genetic_algorithm',
        'target_set_size': 6
    },
    'enhanced_equiphi29': {
        'temperature': 42.0,
        'dmso_percent': 5.0,
        'betaine_m': 1.0,
        'polymerase': 'equiphi29',
        'na_conc': 50.0,
        'mg_conc': 0.0,
        'ssb': True,
        'optimization_method': 'genetic_algorithm',
        'target_set_size': 6
    },
    'long_primers_15mer': {
        'temperature': 45.0,
        'dmso_percent': 7.0,
        'betaine_m': 1.5,
        'polymerase': 'equiphi29',
        'na_conc': 50.0,
        'mg_conc': 0.0,
        'ssb': True,
        'optimization_method': 'genetic_algorithm',
        'target_set_size': 6
    },
    'high_gc_genome': {
        'temperature': 45.0,
        'dmso_percent': 10.0,
        'betaine_m': 2.0,
        'polymerase': 'equiphi29',
        'na_conc': 50.0,
        'mg_conc': 0.0,
        'ssb': True,
        'optimization_method': 'genetic_algorithm',
        'target_set_size': 8
    }
}


# Command groups for organized help display
COMMAND_GROUPS = [
    ('Pipeline', [
        'count-kmers', 'filter', 'score', 'optimize', 'design',
    ]),
    ('Setup', [
        'init', 'start', 'validate', 'validate-params', 'schema', 'show-presets', 'suggest',
    ]),
    ('Results', [
        'interpret', 'report', 'export',
    ]),
    ('Analysis', [
        'analyze-set', 'analyze-genome', 'analyze-dimers', 'analyze-stability',
    ]),
    ('Advanced', [
        'auto-pipeline', 'multi-genome', 'simulate', 'optimize-conditions',
        'build-filter', 'background-list', 'background-add',
        'genome-add', 'genome-list', 'genome-remove',
    ]),
    ('Experimental', [
        'active-learn', 'expand-primers', 'predict-efficiency', 'ml-predict',
        'design-oligos', 'validate-model',
    ]),
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
            cmd_help[choice_action.dest] = choice_action.help or ''

        for group_name, commands in COMMAND_GROUPS:
            parts.append(f'\n  {group_name}:')
            for cmd in commands:
                help_text = cmd_help.get(cmd, '')
                parts.append(f'    {cmd:<24s}{help_text}')

        parts.append('')
        return '\n'.join(parts) + '\n'


def create_parser():
    """Create unified argument parser"""

    from neoswga import __version__

    parser = argparse.ArgumentParser(
        prog='neoswga',
        description='NeoSWGA - Selective Whole Genome Amplification primer design',
        formatter_class=GroupedHelpFormatter,
        epilog='''Quick Start:
  neoswga init --genome target.fasta --background host.fasta
  neoswga count-kmers -j params.json
  neoswga filter -j params.json
  neoswga score -j params.json
  neoswga optimize -j params.json

Run "neoswga <command> --help" for details on a specific command.
        '''
    )
    parser.add_argument('--version', action='version',
                        version=f'%(prog)s {__version__}')

    subparsers = parser.add_subparsers(
        dest='command', title='commands', metavar='<command>'
    )

    # =========================================================================
    # STEP 1: K-mer preprocessing
    # =========================================================================
    count_kmers_parser = subparsers.add_parser('count-kmers',
                                               help='K-mer preprocessing with primer length control')
    add_common_options(count_kmers_parser)
    count_kmers_parser.add_argument('-x', '--fasta-fore', help='Target genome FASTA')
    count_kmers_parser.add_argument('-y', '--fasta-back', help='Background genome FASTA')
    count_kmers_parser.add_argument('-k', '--kmer-fore', help='Target k-mer prefix')
    count_kmers_parser.add_argument('-l', '--kmer-back', help='Background k-mer prefix')

    # Primer length control (defaults from params.json, or 6-12 if not specified)
    count_kmers_parser.add_argument('--min-k', type=int, default=None,
                             help='Minimum primer length (default: from params.json or 6). Supported range: 6-18bp')
    count_kmers_parser.add_argument('--max-k', type=int, default=None,
                             help='Maximum primer length (default: from params.json or 12). Use 15-18bp with DMSO/betaine additives')
    count_kmers_parser.add_argument('--exclusion-genome', type=str, default=None,
                             help='Exclusion genome FASTA file (e.g., mtDNA)')
    count_kmers_parser.add_argument('--blacklist', '-bl', nargs='+', default=None,
                             help='Blacklist genome FASTA file(s) (penalty-weighted filtering)')
    count_kmers_parser.add_argument('--bl-penalty', type=float, default=None,
                             help='Blacklist penalty weight (default: 5.0)')
    count_kmers_parser.add_argument('--max-bl-freq', type=float, default=None,
                             help='Maximum blacklist frequency (default: 0.0, any hit rejects)')

    # =========================================================================
    # STEP 2: Candidate filtering
    # =========================================================================
    filter_parser = subparsers.add_parser('filter',
                                          help='Candidate primer filtering (adaptive GC + reaction conditions)')
    add_common_options(filter_parser)

    # GC Filtering Options
    step2_gc_group = filter_parser.add_argument_group('GC Content Filtering')
    step2_gc_group.add_argument('--gc-tolerance', type=float, default=0.15,
                                help='GC tolerance for adaptive filter (default: 0.15)')
    step2_gc_group.add_argument('--gc-min', type=float,
                                help='Explicit minimum GC content (overrides adaptive)')
    step2_gc_group.add_argument('--gc-max', type=float,
                                help='Explicit maximum GC content (overrides adaptive)')

    # Reaction Conditions
    step2_rxn_group = filter_parser.add_argument_group('Reaction Conditions')
    step2_rxn_group.add_argument('--reaction-temp', type=float,
                                 help='Reaction temperature in C (default: from params.json or polymerase preset)')
    step2_rxn_group.add_argument('--na-conc', type=float,
                                 help='Sodium concentration in mM (default: from params.json or 50.0)')

    # Additives
    step2_add_group = filter_parser.add_argument_group('Additives (enable longer primers & GC-extreme genomes)')
    step2_add_group.add_argument('--dmso-percent', type=float,
                                 help='DMSO concentration 0-10%% (Tm lowering, secondary structure reduction)')
    step2_add_group.add_argument('--betaine-m', type=float,
                                 help='Betaine concentration 0-2.5 M (equalizes AT/GC, enables longer primers)')
    step2_add_group.add_argument('--trehalose-m', type=float,
                                 help='Trehalose concentration 0-1.0 M (Tm lowering)')
    step2_add_group.add_argument('--glycerol-percent', type=float,
                                 help='Glycerol concentration 0-15%% (enzyme stabilizer)')
    step2_add_group.add_argument('--bsa', type=float,
                                 help='BSA concentration 0-400 ug/mL (inhibitor neutralizer)')
    step2_add_group.add_argument('--peg-percent', type=float,
                                 help='PEG concentration 0-15%% (molecular crowding)')
    step2_add_group.add_argument('--mg-conc', type=float,
                                 help='Mg2+ concentration in mM (auto-optimized if not specified)')
    step2_add_group.add_argument('--ssb', action='store_true',
                                 help='Use single-strand binding protein (lowers effective annealing temp)')
    step2_add_group.add_argument('--ethanol-percent', type=float,
                                 help='Ethanol concentration 0-5%% (secondary structure reduction)')
    step2_add_group.add_argument('--urea-m', type=float,
                                 help='Urea concentration 0-2.0 M (denatures GC-rich regions)')
    step2_add_group.add_argument('--tmac-m', type=float,
                                 help='TMAC concentration 0-0.1 M (equalizes AT/GC Tm)')
    step2_add_group.add_argument('--formamide-percent', type=float,
                                 help='Formamide concentration 0-10%% (Tm lowering)')

    # Preset Conditions
    filter_parser.add_argument('--preset', choices=[
        'standard_phi29', 'enhanced_equiphi29', 'high_gc_genome',
        'long_primers_15mer', 'q_solution', 'gc_melt', 'crude_sample', 'low_temp',
        'bst', 'klenow', 'extreme_gc'
    ], help='Use predefined reaction conditions preset')

    # Traditional Filtering Parameters
    step2_trad_group = filter_parser.add_argument_group('Traditional Filtering')
    step2_trad_group.add_argument('--min-fg-freq', type=float,
                                  help='Minimum target genome frequency (default: 1e-5)')
    step2_trad_group.add_argument('--max-bg-freq', type=float,
                                  help='Maximum background genome frequency (default: 5e-6)')
    step2_trad_group.add_argument('--max-gini', type=float,
                                  help='Maximum Gini index for evenness (default: 0.6)')
    step2_trad_group.add_argument('--max-primer', type=int,
                                  help='Number of top primers to keep (default: 500)')
    step2_trad_group.add_argument('--min-tm', type=float,
                                  help='Minimum melting temperature in C (default: 15)')
    step2_trad_group.add_argument('--max-tm', type=float,
                                  help='Maximum melting temperature in C (default: 45)')
    step2_trad_group.add_argument('--max-dimer-bp', type=int,
                                  help='Maximum heterodimer base pairs (default: 3)')
    step2_trad_group.add_argument('--max-self-dimer-bp', type=int,
                                  help='Maximum self-dimer base pairs (default: 4)')

    # Background Bloom Filter (for large background genomes)
    step2_bloom_group = filter_parser.add_argument_group('Background Filtering (for large genomes)')
    step2_bloom_group.add_argument('--use-bloom-filter', action='store_true', default=False,
                                   help='Use Bloom filter for background filtering (memory-efficient for human genome)')
    step2_bloom_group.add_argument('--bloom-filter-path', type=str,
                                   help='Path to pre-built Bloom filter (.pkl). Build with: neoswga build-filter')
    step2_bloom_group.add_argument('--sampled-index-path', type=str,
                                   help='Path to pre-built sampled index (.pkl) for count estimation')
    step2_bloom_group.add_argument('--bloom-max-bg-matches', type=int, default=10,
                                   help='Maximum background matches via Bloom filter (default: 10)')

    # Exclusion Genome Filtering (zero-tolerance for contaminants)
    step2_excl_group = filter_parser.add_argument_group(
        'Exclusion Genome Filtering (reject primers binding contaminant sequences)')
    step2_excl_group.add_argument('--exclusion-genome', type=str, default=None,
                                  help='Exclusion genome FASTA file (e.g., mtDNA, chloroplast). '
                                       'Primers binding this genome are rejected.')
    step2_excl_group.add_argument('--excl-threshold', type=int, default=0,
                                  help='Maximum allowed hits in exclusion genome '
                                       '(default: 0 = any hit rejects primer)')

    # Blacklist Genome Filtering (penalty-weighted)
    step2_bl_group = filter_parser.add_argument_group(
        'Blacklist Genome Filtering (penalty-weighted filtering of contaminating sequences)')
    step2_bl_group.add_argument('--blacklist', '-bl', nargs='+', default=None,
                                help='Blacklist genome FASTA file(s)')
    step2_bl_group.add_argument('--bl-penalty', type=float, default=None,
                                help='Blacklist penalty weight (default: 5.0)')
    step2_bl_group.add_argument('--max-bl-freq', type=float, default=None,
                                help='Maximum blacklist frequency (default: 0.0, any hit rejects)')

    # =========================================================================
    # STEP 3: Random forest scoring
    # =========================================================================
    score_parser = subparsers.add_parser('score',
                                         help='Amplification efficacy scoring with ML threshold control')
    add_common_options(score_parser)

    # Scoring method (merges ml-predict)
    score_parser.add_argument('--method', choices=['rf', 'ml', 'both'],
                             default='rf',
                             help='Scoring method: rf (random forest, default), ml (deep learning), both (compare methods)')
    score_parser.add_argument('--min-amp-pred', type=float,
                             help='Minimum amplification prediction score (default: 10)')
    score_parser.add_argument('--fast-score', action='store_true',
                             help='Skip thermodynamic histogram features for ~100x speedup. '
                                  'These features contribute <2%% of model accuracy but >99%% '
                                  'of computation time. Recommended for large primer pools.')

    # ML-specific options (when --method ml or both)
    score_parser.add_argument('--model-path', type=str,
                             help='Path to pre-trained deep learning model')
    score_parser.add_argument('--embedding-dim', type=int, default=128,
                             help='Embedding dimension for deep learning (default: 128)')

    # Enhanced feature engineering options
    score_parser.add_argument('--use-enhanced-features', action='store_true',
                             help='Use enhanced 120+ feature model (requires enhanced_rf_model.pkl)')
    score_parser.add_argument('--enhanced-model-path', type=str, default=None,
                             help='Path to enhanced random forest model')

    # Blacklist metadata (filtering happens at step2; this is a passthrough so
    # params.json contains the full record for audit trails).
    score_parser.add_argument('--blacklist', '-bl', nargs='+', default=None,
                             help='Blacklist genome FASTA file(s); informational at this step, filtering happens in step2.')

    # =========================================================================
    # STEP 4: Primer set optimization
    # =========================================================================
    optimize_parser = subparsers.add_parser('optimize',
                                            help='Primer set optimization (network-based + experimental)')
    add_common_options(optimize_parser)

    # Method Selection
    opt_method_group = optimize_parser.add_argument_group('Method Selection')
    opt_method_group.add_argument('-m', '--optimization-method',
                             choices=['hybrid', 'greedy', 'dominating-set', 'weighted-set-cover',
                                     'network', 'genetic', 'background-aware',
                                     'milp', 'equiphi29', 'moea', 'normalized', 'tiling',
                                     'clique', 'coverage-then-dimerfree', 'dimerfree-scored',
                                     'bg-prefilter', 'bg-prefilter-hybrid'],
                             default='hybrid',
                             help='Optimization method. '
                                  'Decision tree: '
                                  'hybrid (default, general use), '
                                  'dominating-set (speed-critical, large pools), '
                                  'background-aware (clinical, 10-20x bg reduction), '
                                  'clique (dimer-free requirement), '
                                  'equiphi29 (GC-rich genomes >65%%), '
                                  'milp (exact solution, small pools). '
                                  'Pipelines: '
                                  'coverage-then-dimerfree (DS->clique cascade), '
                                  'dimerfree-scored (clique->network scoring). '
                                  'Use --method-guide for detailed comparison.')
    opt_method_group.add_argument('--strategy',
                             choices=['clinical', 'discovery', 'fast', 'balanced', 'enrichment'],
                             default='balanced',
                             help='Strategy preset for normalized optimizer: '
                                  'clinical (high specificity), discovery (max coverage), '
                                  'fast (quick screening), balanced (equal weights), '
                                  'enrichment (sequencing). Only used with --optimization-method=normalized')
    opt_method_group.add_argument('--method-guide', action='store_true',
                             help='Show optimization method selection guide and exit')

    # Primer Strategy
    opt_strategy_group = optimize_parser.add_argument_group('Primer Strategy')
    opt_strategy_group.add_argument('--use-cooperative-binding', action='store_true',
                             help='[EXPERIMENTAL] Cooperative binding model (not yet fully integrated)')
    opt_strategy_group.add_argument('--primer-strategy', choices=['standard', 'hybrid'],
                             default='standard',
                             help='Primer design strategy (standard: uniform length, hybrid: mixed lengths)')

    # Background Filtering
    opt_bg_group = optimize_parser.add_argument_group('Background Filtering')
    opt_bg_group.add_argument('--use-background-filter', action='store_true', default=False,
                             help='Use Bloom filter for background filtering')
    opt_bg_group.add_argument('--no-bg-prefilter', action='store_true', default=False,
                             help='Disable automatic background pre-filtering of candidates '
                                  '(enabled by default when background genome data is available)')
    opt_bg_group.add_argument('--background-bloom-path', type=str,
                             help='Path to pre-built Bloom filter')
    opt_bg_group.add_argument('--background-sampled-path', type=str,
                             help='Path to pre-built sampled index')
    opt_bg_group.add_argument('--no-background', action='store_true', default=False,
                             help='Host-free mode: optimize without background genome data. '
                                  'Relies on intrinsic primer quality (Tm, complexity, evenness) '
                                  'rather than fg/bg selectivity. Use when background genome is '
                                  'unknown or unavailable.')

    # Performance
    opt_perf_group = optimize_parser.add_argument_group('Performance')
    opt_perf_group.add_argument('--use-position-cache', action='store_true', default=True,
                             help='Use in-memory position cache (default: True, 1000x speedup)')
    opt_perf_group.add_argument('--no-position-cache', action='store_false', dest='use_position_cache',
                             help='Disable position cache (slower)')
    opt_perf_group.add_argument('--seed', type=int, default=None,
                             help='Random seed for reproducible results. '
                                  'Affects stochastic optimizers (genetic, moea). '
                                  'Required for clinical/regulated workflows.')
    opt_perf_group.add_argument('-n', '--num-primers', type=int,
                             help='Number of primers to select (default: from params.json or 6)')
    opt_perf_group.add_argument('--max-optimization-time', type=int, default=300,
                             help='Maximum optimization time in seconds (default: 300)')
    opt_perf_group.add_argument('--max-extension', type=int, default=70000,
                             help='Maximum Phi29 extension length in bp (default: 70000)')

    # Set Size & Application
    opt_size_group = optimize_parser.add_argument_group('Set Size & Application')
    opt_size_group.add_argument('--auto-size', action='store_true',
                             help='Automatically determine optimal primer set size based on '
                                  'application profile and reaction conditions')
    opt_size_group.add_argument('--application', type=str,
                             choices=['discovery', 'clinical', 'enrichment', 'metagenomics'],
                             default='enrichment',
                             help='Application profile for auto-sizing: '
                                  'discovery (maximize sensitivity), '
                                  'clinical (minimize false positives), '
                                  'enrichment (balanced, default), '
                                  'metagenomics (capture diversity)')
    opt_size_group.add_argument('--min-fg-bg-ratio', type=float,
                             help='Minimum foreground/background binding site ratio. '
                                  'Overrides application profile default. Higher values = more selective.')
    opt_size_group.add_argument('--show-frontier', action='store_true',
                             help='Show Pareto frontier of coverage vs fg/bg ratio tradeoffs '
                                  '(requires matplotlib for plotting)')
    opt_size_group.add_argument('--quick-estimate', action='store_true',
                             help='Use quick estimation only for auto-size (skip full optimization at multiple sizes)')

    # Mechanistic Model
    opt_mech_group = optimize_parser.add_argument_group('Mechanistic Model')
    opt_mech_group.add_argument('--use-mechanistic-model', action='store_true',
                             help='Use mechanistic model for primer weighting during optimization')
    opt_mech_group.add_argument('--mechanistic-weight', type=float, default=0.3,
                             help='Weight for mechanistic model in scoring (0.0-1.0, default: 0.3)')
    opt_mech_group.add_argument('--template-gc', type=float,
                             help='Template genome GC content (0-1). Auto-detected if not specified.')
    opt_mech_group.add_argument('--uniformity-weight', type=float, default=0.0,
                             help='Weight for coverage uniformity in scoring (0.0-1.0, default: 0.0). '
                                  'Higher values prioritize even genome coverage over raw enrichment.')

    # Post-processing
    opt_post_group = optimize_parser.add_argument_group('Post-processing')
    opt_post_group.add_argument('--minimize-primers', action='store_true',
                             help='Post-process to minimize primer count while maintaining coverage')
    opt_post_group.add_argument('--target-coverage', type=float, default=0.70,
                             help='Target genome coverage fraction when minimizing primers (default: 0.70)')

    # Validation
    opt_val_group = optimize_parser.add_argument_group('Validation')
    opt_val_group.add_argument('--validate-simulation', action='store_true',
                             help='Validate results with stochastic simulation (Gillespie algorithm)')
    opt_val_group.add_argument('--simulation-time', type=float, default=3600.0,
                             help='Simulation time in seconds (default: 3600)')

    # Blacklist passthrough: filtering already happened at step2; included for
    # CLI surface consistency and audit-trail parity with params.json.
    optimize_parser.add_argument('--blacklist', '-bl', nargs='+', default=None,
                             help='Blacklist genome FASTA file(s); informational at this step, filtering happens in step2.')

    # =========================================================================
    # UNIFIED: Complete pipeline (all 4 steps)
    # =========================================================================
    design_parser = subparsers.add_parser('design',
                                          help='Run complete primer design pipeline (all 4 steps sequentially)')
    add_common_options(design_parser)

    # Auto-tuning option (merges auto-pipeline)
    design_parser.add_argument('--auto', action='store_true',
                              help='Enable automatic parameter optimization')
    design_parser.add_argument('--auto-iterations', type=int, default=5,
                              help='Number of auto-tuning iterations (default: 5)')

    # Multi-genome option (merges multi-genome)
    design_parser.add_argument('--multi-genome', nargs='+', metavar='GENOME',
                              help='Design pan-genome primers for multiple target genomes')
    design_parser.add_argument('--min-coverage', type=float, default=0.8,
                              help='Minimum fraction of genomes to cover (default: 0.8)')

    # Step control
    design_parser.add_argument('--start-from', type=int, choices=[1, 2, 3, 4],
                              help='Start from step: 1=count-kmers, 2=filter, 3=score, 4=optimize')
    design_parser.add_argument('--stop-at', type=int, choices=[1, 2, 3, 4],
                              help='Stop at step: 1=count-kmers, 2=filter, 3=score, 4=optimize')

    # =========================================================================
    # UTILITY: Build background filter
    # =========================================================================
    build_filter_parser = subparsers.add_parser('build-filter',
                                               help='Build background Bloom filter (one-time setup)')
    build_filter_parser.add_argument('--genome', required=True, help='Path to background genome FASTA (or k-mer prefix with --from-kmers)')
    build_filter_parser.add_argument('-o', '--output', '--output-dir', required=True,
                                    dest='output_dir',
                                    help='Output directory for filter files')
    build_filter_parser.add_argument('--from-kmers', action='store_true',
                                    help='Build from jellyfish k-mer files (MUCH faster for large genomes)')
    build_filter_parser.add_argument('--min-k', type=int, default=6,
                                    help='Minimum k-mer length (default: 6)')
    build_filter_parser.add_argument('--max-k', type=int, default=12,
                                    help='Maximum k-mer length (default: 12)')
    build_filter_parser.add_argument('--force', action='store_true',
                                    help='Rebuild even if filter exists')
    build_filter_parser.add_argument('--capacity', type=int,
                                    help='Bloom filter capacity (default: auto from genome size)')
    build_filter_parser.add_argument('--error-rate', type=float, default=0.01,
                                    help='Bloom filter error rate (default: 0.01)')
    build_filter_parser.add_argument('-v', '--verbose', action='store_true',
                                    help='Verbose output')
    build_filter_parser.add_argument('-q', '--quiet', action='store_true',
                                    help='Minimal output')

    # =========================================================================
    # UTILITY: Validate installation
    # =========================================================================
    validate_parser = subparsers.add_parser('validate',
                                           help='Validate installation and run tests')
    validate_parser.add_argument('--quick', action='store_true',
                                help='Run quick validation (subset of tests)')
    validate_parser.add_argument('--all', action='store_true',
                                help='Run all tests including benchmarks')

    # =========================================================================
    # UTILITY: Show presets
    # =========================================================================
    subparsers.add_parser('show-presets',
                         help='Show available reaction condition presets')

    # =========================================================================
    # SETUP: Initialize new project (setup wizard)
    # =========================================================================
    init_parser = subparsers.add_parser('init',
                                        help='Set up a new primer design project with guided configuration')
    init_parser.add_argument('--genome', '-g', required=True,
                            help='Target genome FASTA file')
    init_parser.add_argument('--background', '-b',
                            help='Background genome FASTA file (optional)')
    init_parser.add_argument('--output', '-o', default='params.json',
                            help='Output params.json path (default: params.json)')
    init_parser.add_argument('--output-dir', default='results',
                            help='Directory for pipeline output files (default: results)')
    init_parser.add_argument('--advanced', '-a', action='store_true',
                            help='Show advanced configuration options')
    init_parser.add_argument('--yes', '-y', action='store_true',
                            help='Auto-approve without prompting')
    init_parser.add_argument('--non-interactive', action='store_true',
                            help='Run non-interactively with defaults')
    init_parser.add_argument('--blacklist', '-bl', nargs='+', default=None,
                            help='Blacklist genome FASTA file(s) for contamination filtering')

    # =========================================================================
    # SETUP: Validate parameters
    # =========================================================================
    validate_params_parser = subparsers.add_parser('validate-params',
                                                   help='Validate params.json configuration before running')
    validate_params_parser.add_argument('-j', '--json-file', required=True,
                                       help='Parameters JSON file to validate')

    # =========================================================================
    # SETUP: Dump canonical JSON schema
    # =========================================================================
    schema_parser = subparsers.add_parser('schema',
                                          help='Inspect or dump the canonical params.json schema')
    schema_parser.add_argument('--dump', action='store_true',
                               help='Print the params.json JSON Schema to stdout')
    schema_parser.add_argument('--output', '-o', type=str, default=None,
                               help='Write the schema to the given file instead of stdout')

    # =========================================================================
    # SETUP: Validate mechanistic model
    # =========================================================================
    validate_model_parser = subparsers.add_parser('validate-model',
                                                   help='Validate mechanistic model against expected behavior')
    validate_model_parser.add_argument('--verbose', '-v', action='store_true',
                                       help='Show detailed results for each test')
    validate_model_parser.add_argument('--output-json', action='store_true',
                                       help='Output results as JSON')

    # =========================================================================
    # SETUP: Interpret results
    # =========================================================================
    interpret_parser = subparsers.add_parser('interpret',
                                             help='Interpret pipeline results and provide quality assessment')
    interpret_parser.add_argument('-d', '--dir', required=True,
                                 help='Results directory containing step4_improved_df.csv')

    # =========================================================================
    # SETUP: Generate report
    # =========================================================================
    report_parser = subparsers.add_parser('report',
                                          help='Generate comprehensive HTML report for primer set')
    report_parser.add_argument('-d', '--dir', required=True,
                               help='Results directory containing pipeline output')
    report_parser.add_argument('-o', '--output',
                               help='Output file path (default: <dir>/report.html)')
    report_parser.add_argument('--format', choices=['html', 'json'], default='html',
                               help='Output format (default: html)')
    report_parser.add_argument('--level', choices=['summary', 'full'], default='summary',
                               help='Report level: summary (1-page) or full (technical report)')
    report_parser.add_argument('--check', action='store_true',
                               help='Validate only, do not generate report')
    report_parser.add_argument('--interactive', action='store_true',
                               help='Include interactive Plotly charts (requires plotly)')
    report_parser.add_argument('-q', '--quiet', action='store_true',
                               help='Suppress progress messages')

    # =========================================================================
    # EXPORT: Primer export for synthesis ordering
    # =========================================================================
    export_parser = subparsers.add_parser(
        'export',
        help='Export primers for synthesis ordering',
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
"""
    )
    export_parser.add_argument('-d', '--dir', required=True,
                               help='Results directory containing step4_improved_df.csv')
    export_parser.add_argument('-o', '--output', default='./export',
                               help='Output directory (default: ./export)')
    export_parser.add_argument('--project', default='SWGA',
                               help='Project name for file naming (default: SWGA)')
    export_parser.add_argument('--vendor', default='idt',
                               choices=['idt', 'twist', 'sigma', 'generic'],
                               help='Vendor format for CSV (default: idt)')
    export_parser.add_argument('--format', choices=['all', 'fasta', 'csv', 'protocol', 'bed', 'bedgraph'],
                               default='all',
                               help='Export format (default: all)')
    export_parser.add_argument('-j', '--json-file',
                               help='params.json for reaction conditions')
    export_parser.add_argument('--modifications',
                               choices=['none', 'standard', 'low-input'],
                               default='standard',
                               help="Modification profile: none (bare primers), "
                                    "standard (3' PTO for exonuclease protection), "
                                    "low-input (3' PTO + 5' C18 spacer for tdMDA) "
                                    "(default: standard)")
    export_parser.add_argument('--pto-bonds', type=int, default=None,
                               help="Override number of 3' phosphorothioate bonds (default: 2)")
    export_parser.add_argument('--no-modifications', action='store_true',
                               help='Export bare sequences without modifications '
                                    '(equivalent to --modifications none)')
    export_parser.add_argument('--genome-name', default='genome',
                               help='Chromosome/contig name for BED/BedGraph output (default: genome)')
    export_parser.add_argument('--window-size', type=int, default=1000,
                               help='Window size in bp for BedGraph binning (default: 1000)')
    export_parser.set_defaults(func=run_export)

    # =========================================================================
    # SETUP: Interactive workflow selector
    # =========================================================================
    subparsers.add_parser('start',
                         help='Interactive menu to discover and launch neoswga features')

    # =========================================================================
    # SETUP: Suggest reaction conditions
    # =========================================================================
    suggest_parser = subparsers.add_parser('suggest',
                                           help='Suggest optimal reaction conditions based on genome and primer length')
    suggest_parser.add_argument('--genome-gc', type=float,
                               help='Target genome GC content (0-1, e.g., 0.65 for 65%%)')
    suggest_parser.add_argument('--genome', type=str,
                               help='Target genome FASTA file (calculates GC automatically)')
    suggest_parser.add_argument('--primer-length', '-l', type=int,
                               help='Target primer length in bp')
    suggest_parser.add_argument('--context', '-c', choices=['standard', 'clinical', 'high_throughput', 'low_input'],
                               default='standard',
                               help='Application context (default: standard)')
    suggest_parser.add_argument('--polymerase', choices=['phi29', 'equiphi29', 'bst', 'klenow'],
                               default='phi29',
                               help='Polymerase type (default: phi29)')
    suggest_parser.add_argument('--optimize-for', choices=['amplification', 'specificity', 'coverage', 'processivity'],
                               default='amplification',
                               help='Optimization goal (default: amplification)')
    suggest_parser.add_argument('--use-optimizer', action='store_true',
                               help='Use advanced multi-additive optimizer (grid search)')
    suggest_parser.add_argument('--sweep', action='store_true',
                               help='Sweep a grid of reaction conditions and rank by '
                                    'predicted amplification factor')
    suggest_parser.add_argument('--kmer-range', action='store_true',
                               help='Print a recommended (min_k, max_k) range based on genome size, '
                                    'GC, and polymerase, then exit')
    suggest_parser.add_argument('--genome-size', type=int,
                               help='Target genome size in bp (for --kmer-range). Auto-derived from --genome if that flag is given.')
    suggest_parser.add_argument('--output', '-o', type=str,
                               help='Output CSV path for condition sweep results')

    # =========================================================================
    # ADVANCED: Optimize conditions
    # =========================================================================
    optimize_cond_parser = subparsers.add_parser('optimize-conditions',
                                                 help='Find optimal reaction conditions for a genome')
    optimize_cond_parser.add_argument('--fg', '--foreground', required=True,
                                     help='Foreground genome FASTA file')
    optimize_cond_parser.add_argument('--output', '-o', required=True,
                                     help='Output directory')
    optimize_cond_parser.add_argument('--target-k', type=int,
                                     help='Target k-mer length to optimize for')

    # =========================================================================
    # ADVANCED: Analyze primer set
    # =========================================================================
    analyze_parser = subparsers.add_parser('analyze-set',
                                          help='[EXPERIMENTAL] Analyze an existing primer set')
    analyze_parser.add_argument('--primers', required=True, nargs='+',
                               help='Primer sequences to analyze')
    analyze_parser.add_argument('--fg', required=True,
                               help='Foreground genome FASTA')
    analyze_parser.add_argument('--fg-kmers', required=True,
                               help='Foreground k-mer file prefix')
    analyze_parser.add_argument('--output', '-o', required=True,
                               help='Output directory')
    analyze_parser.add_argument('--preset', default='standard_phi29',
                               help='Reaction conditions preset')
    analyze_parser.add_argument('--simulate', action='store_true',
                               help='Run replication simulation')

    # =========================================================================
    # CATEGORY 1: Genome analysis (orphaned feature)
    # =========================================================================
    analyze_genome_parser = subparsers.add_parser('analyze-genome',
                                                   help='[EXPERIMENTAL] Analyze genome suitability for SWGA')
    analyze_genome_parser.add_argument('--genome', required=True,
                                      help='Genome FASTA file to analyze')
    analyze_genome_parser.add_argument('--output', '-o', required=True,
                                      help='Output directory for analysis report')
    analyze_genome_parser.add_argument('--window-size', type=int, default=1000,
                                      help='Window size for GC profiling (default: 1000bp)')

    # =========================================================================
    # CATEGORY 1: Dimer network analysis (orphaned feature)
    # =========================================================================
    analyze_dimers_parser = subparsers.add_parser('analyze-dimers',
                                                   help='[EXPERIMENTAL] Analyze primer dimer interaction network')
    analyze_dimers_parser.add_argument('--primers', required=True, nargs='+',
                                      help='Primer sequences to analyze')
    analyze_dimers_parser.add_argument('--output', '-o', required=True,
                                      help='Output directory for network visualization')
    analyze_dimers_parser.add_argument('--threshold', type=float, default=0.3,
                                      help='Dimer severity threshold (default: 0.3)')
    analyze_dimers_parser.add_argument('--visualize', action='store_true',
                                      help='Generate network visualization')

    # =========================================================================
    # CATEGORY 1: 3' stability analysis (orphaned feature)
    # =========================================================================
    analyze_stability_parser = subparsers.add_parser('analyze-stability',
                                                     help="Analyze 3' end stability and specificity")
    analyze_stability_parser.add_argument('--primers', required=True, nargs='+',
                                         help='Primer sequences to analyze')
    analyze_stability_parser.add_argument('--output', '-o', required=True,
                                         help='Output file for stability report')
    analyze_stability_parser.add_argument('--temp', type=float, default=37.0,
                                         help='Reaction temperature in C (default: 37)')

    # =========================================================================
    # CATEGORY 1: ML-based prediction (orphaned feature)
    # =========================================================================
    ml_predict_parser = subparsers.add_parser('ml-predict',
                                              help='Deep learning-based amplification prediction')
    ml_predict_parser.add_argument('--primers-file', '--primers', required=True,
                                  dest='primers_file',
                                  help='File containing primer sequences (one per line)')
    ml_predict_parser.add_argument('--output', '-o', required=True,
                                  help='Output directory for predictions')
    ml_predict_parser.add_argument('--model', type=str,
                                  help='Path to pre-trained model (optional)')
    ml_predict_parser.add_argument('--use-gpu', action='store_true',
                                  help='Use GPU acceleration for inference')

    # =========================================================================
    # CATEGORY 1: Optimal oligo generator (orphaned feature)
    # =========================================================================
    design_oligos_parser = subparsers.add_parser('design-oligos',
                                                  help='Alternative comprehensive primer design system')
    design_oligos_parser.add_argument('--genome', required=True,
                                     help='Target genome FASTA file')
    design_oligos_parser.add_argument('--output', '-o', required=True,
                                     help='Output directory for designed primers')
    design_oligos_parser.add_argument('--num-primers', type=int, default=10,
                                     help='Number of primers to design (default: 10)')
    design_oligos_parser.add_argument('--min-k', type=int, default=6,
                                     help='Minimum primer length (default: 6)')
    design_oligos_parser.add_argument('--max-k', type=int, default=12,
                                     help='Maximum primer length (default: 12)')

    # =========================================================================
    # CATEGORY 3: Auto-tuning pipeline (orphaned feature)
    # =========================================================================
    auto_pipeline_parser = subparsers.add_parser('auto-pipeline',
                                                  help='Automatic parameter optimization pipeline (experimental, use standard pipeline instead)')
    auto_pipeline_parser.add_argument('-j', '--json-file', required=True,
                                     help='Base parameters JSON file')
    auto_pipeline_parser.add_argument('--iterations', type=int, default=5,
                                     help='Number of optimization iterations (default: 5)')
    auto_pipeline_parser.add_argument('--output', '-o', required=True,
                                     help='Output directory for optimized parameters')

    # =========================================================================
    # CATEGORY 3: Multi-genome pipeline (orphaned feature)
    # =========================================================================
    multi_genome_parser = subparsers.add_parser('multi-genome',
                                                help='Pan-genome primer design for multiple targets')
    multi_genome_parser.add_argument('--genomes', required=True, nargs='+',
                                    help='List of target genome FASTA files')
    multi_genome_parser.add_argument('--background', nargs='+',
                                    help='Background genome FASTA files (e.g., host DNA)')
    multi_genome_parser.add_argument('--blacklist', nargs='+',
                                    help='Blacklist genome FASTA files (strongly avoid, e.g., other pathogens)')
    multi_genome_parser.add_argument('--output', '-o', required=True,
                                    help='Output directory')
    multi_genome_parser.add_argument('--num-primers', type=int, default=12,
                                    help='Number of primers to design (default: 12)')
    multi_genome_parser.add_argument('--min-k', type=int, default=8,
                                    help='Minimum primer length (default: 8)')
    multi_genome_parser.add_argument('--max-k', type=int, default=12,
                                    help='Maximum primer length (default: 12)')
    multi_genome_parser.add_argument('--polymerase',
                                    choices=['phi29', 'equiphi29', 'bst', 'klenow'],
                                    help='Polymerase type: phi29 (30C), equiphi29 (42C), '
                                         'bst (63C), klenow (37C)')
    multi_genome_parser.add_argument('--validate-simulation', action='store_true',
                                    help='Validate with simulation')
    multi_genome_parser.add_argument('--quiet', action='store_true',
                                    help='Minimal output')

    # =========================================================================
    # CATEGORY 4: Simulation command (orphaned features)
    # =========================================================================
    simulate_parser = subparsers.add_parser('simulate',
                                           help='Agent-based replication simulation')
    sim_input = simulate_parser.add_mutually_exclusive_group(required=True)
    sim_input.add_argument('--primers', nargs='+',
                                help='Primer sequences to simulate')
    sim_input.add_argument('--from-results', type=str, metavar='DIR_OR_CSV',
                                help='Read primers from pipeline results directory or step4 CSV file '
                                     '(e.g., --from-results results/ or --from-results step4_improved_df.csv)')
    simulate_parser.add_argument('--genome', required=True,
                                help='Target genome FASTA file')
    simulate_parser.add_argument('--output', '-o', required=True,
                                help='Output directory for simulation results')
    simulate_parser.add_argument('--duration', '--cycles', type=int, default=60,
                                help='Simulation duration in minutes (default: 60)')
    simulate_parser.add_argument('--replicates', type=int, default=5,
                                help='Number of simulation replicates (default: 5)')
    simulate_parser.add_argument('--polymerase',
                                choices=['phi29', 'equiphi29', 'bst', 'klenow'],
                                default='phi29',
                                help='Polymerase type: phi29 (30C), equiphi29 (42C), '
                                     'bst (63C), klenow (37C) (default: phi29)')
    simulate_parser.add_argument('--visualize', action='store_true',
                                help='Generate visualization plots')
    simulate_parser.add_argument('--report', action='store_true',
                                help='Generate HTML report')

    # =========================================================================
    # CATEGORY 5: Active learning command (experimental)
    # =========================================================================
    active_learn_parser = subparsers.add_parser('active-learn',
                                                help='[EXPERIMENTAL] Active learning for iterative primer optimization')
    active_learn_parser.add_argument('-j', '--json-file', required=True,
                                     help='Parameters JSON file (from optimize step)')
    active_learn_parser.add_argument('--experimental-results',
                                     help='CSV file with experimental outcomes (optional)')
    active_learn_parser.add_argument('--num-candidates', type=int, default=10,
                                     help='Number of candidate primer sets to generate (default: 10)')
    active_learn_parser.add_argument('--max-primers', type=int, default=15,
                                     help='Maximum primers per set (default: 15)')
    active_learn_parser.add_argument('--exploration-weight', type=float, default=2.0,
                                     help='Balance exploration vs exploitation (default: 2.0)')
    active_learn_parser.add_argument('--output', '-o', required=True,
                                     help='Output directory for active learning results')
    active_learn_parser.add_argument('--quiet', '-q', action='store_true',
                                     help='Suppress progress output')

    # Expand primers - add new primers to existing validated set
    expand_parser = subparsers.add_parser('expand-primers',
        help='Expand existing primer set with additional primers to fill coverage gaps')
    expand_parser.add_argument('-j', '--json-file', required=True,
                               help='Parameters JSON file')
    expand_parser.add_argument('--fixed-primers', nargs='+', required=True,
                               help='Primer sequences to keep (already validated)')
    expand_parser.add_argument('--fixed-primers-file',
                               help='File with fixed primers (one per line)')
    expand_parser.add_argument('--failed-primers', nargs='+',
                               help='Primer sequences to exclude (failed in wet lab)')
    expand_parser.add_argument('--failed-primers-file',
                               help='File with failed primers to exclude (one per line)')
    expand_parser.add_argument('--num-new', type=int, default=6,
                               help='Number of new primers to add (default: 6)')
    expand_parser.add_argument('--optimization-method', default='hybrid',
                               choices=['hybrid', 'dominating-set'],
                               help='Optimization method (default: hybrid)')
    expand_parser.add_argument('--output', '-o', required=True,
                               help='Output directory for expanded primer set')
    expand_parser.add_argument('--quiet', '-q', action='store_true',
                               help='Suppress progress output')

    # Predict efficiency - unified confidence score for primer set
    predict_parser = subparsers.add_parser('predict-efficiency',
        help='Predict efficiency of primer set before synthesis')
    predict_parser.add_argument('-j', '--json-file', required=True,
                                help='Parameters JSON file')
    predict_parser.add_argument('--primers', nargs='+', required=True,
                                help='Primer sequences to evaluate')
    predict_parser.add_argument('--primers-file',
                                help='File with primers (one per line)')
    predict_parser.add_argument('--run-simulation', action='store_true',
                                help='Run simulation for more accurate prediction (slower)')
    predict_parser.add_argument('--track', action='store_true',
                                help='Record prediction for later outcome tracking')
    predict_parser.add_argument('--output', '-o',
                                help='Output JSON file for prediction results')
    predict_parser.add_argument('--quiet', '-q', action='store_true',
                                help='Suppress progress output')

    # Background registry - list available pre-computed backgrounds
    bg_list_parser = subparsers.add_parser('background-list',
        help='List available pre-computed background genomes')
    bg_list_parser.add_argument('--discover', action='store_true',
                                help='Search for new backgrounds in default directories')
    bg_list_parser.add_argument('--search',
                                help='Search for backgrounds by name or species')
    bg_list_parser.add_argument('--quiet', '-q', action='store_true',
                                help='Suppress verbose output')

    # Background registry - add new background
    bg_add_parser = subparsers.add_parser('background-add',
        help='Add a pre-computed background to the registry')
    bg_add_parser.add_argument('--name', required=True,
                               help='Human-readable name (e.g., "Human GRCh38")')
    bg_add_parser.add_argument('--species', required=True,
                               help='Species name (e.g., "Homo sapiens")')
    bg_add_parser.add_argument('--bloom-path',
                               help='Path to Bloom filter pickle file')
    bg_add_parser.add_argument('--kmer-prefix',
                               help='Prefix for k-mer files (without _Xmer_all.txt)')
    bg_add_parser.add_argument('--genome-size', type=int, default=0,
                               help='Genome size in bp')
    bg_add_parser.add_argument('--min-k', type=int, default=6,
                               help='Minimum k-mer length (default: 6)')
    bg_add_parser.add_argument('--max-k', type=int, default=12,
                               help='Maximum k-mer length (default: 12)')
    bg_add_parser.add_argument('--description',
                               help='Additional description')
    bg_add_parser.add_argument('--overwrite', action='store_true',
                               help='Overwrite existing entry with same name')

    # =========================================================================
    # UTILITY: Genome library management
    # =========================================================================
    genome_add_parser = subparsers.add_parser('genome-add',
                                               help='Add a genome to the pre-calculated library')
    genome_add_parser.add_argument('name', help='Identifier for the genome (e.g., human-grch38)')
    genome_add_parser.add_argument('fasta', help='Path to genome FASTA file')
    genome_add_parser.add_argument('--role', choices=['bg', 'bl'], default='bg',
                                   help='Role: bg (background) or bl (blacklist). Default: bg')
    genome_add_parser.add_argument('--species', default='',
                                   help='Species name (optional)')
    genome_add_parser.add_argument('--k-ranges', default='6-12,12-18',
                                   help='K-mer ranges to compute (default: 6-12,12-18)')
    genome_add_parser.add_argument('--no-bloom', action='store_true',
                                   help='Skip Bloom filter construction')

    genome_list_parser = subparsers.add_parser('genome-list',
                                                help='List genomes in the pre-calculated library')
    genome_list_parser.add_argument('--verbose', '-v', action='store_true',
                                    help='Show detailed information')

    genome_remove_parser = subparsers.add_parser('genome-remove',
                                                  help='Remove a genome from the pre-calculated library')
    genome_remove_parser.add_argument('name', help='Identifier of genome to remove')

    # Validate that all registered commands appear in COMMAND_GROUPS
    grouped_cmds = {cmd for _, cmds in COMMAND_GROUPS for cmd in cmds}
    registered_cmds = set(subparsers.choices.keys()) if hasattr(subparsers, 'choices') else set()
    missing = registered_cmds - grouped_cmds
    if missing:
        import warnings
        warnings.warn(
            f"Commands not listed in COMMAND_GROUPS: {', '.join(sorted(missing))}",
            stacklevel=2,
        )

    return parser


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
        'standard_phi29': rc.get_standard_conditions(),
        'enhanced_equiphi29': rc.get_enhanced_conditions(),
        'high_gc_genome': rc.get_high_gc_conditions(),
        'long_primers_15mer': rc.get_enhanced_conditions(),  # Same as enhanced but user intent is longer primers
        'q_solution': rc.get_q_solution_equivalent(),
        'gc_melt': rc.get_gc_melt_conditions(),
        'crude_sample': rc.get_crude_sample_conditions(),
        'low_temp': rc.get_low_temp_conditions(),
        'bst': rc.get_bst_conditions(),
        'klenow': rc.get_klenow_conditions(),
        'extreme_gc': rc.get_extreme_gc_conditions()
    }

    if preset_name not in preset_map:
        logger.error(f"Unknown preset: {preset_name}")
        return {}

    conditions = preset_map[preset_name]

    # Convert ReactionConditions to dict
    return {
        'reaction_temp': conditions.temp,
        'polymerase': conditions.polymerase,
        'dmso_percent': conditions.dmso_percent,
        'betaine_m': conditions.betaine_m,
        'trehalose_m': conditions.trehalose_m,
        'glycerol_percent': conditions.glycerol_percent,
        'bsa_ug_ml': conditions.bsa_ug_ml,
        'peg_percent': conditions.peg_percent,
        'na_conc': conditions.na_conc,
        'mg_conc': conditions.mg_conc,
        'ssb': conditions.ssb
    }


def add_common_options(parser):
    """Add options common to all steps"""
    parser.add_argument('-j', '--json-file', type=str,
                       help='Parameters JSON file')
    parser.add_argument('-z', '--data-dir', type=str, help='Data directory')
    parser.add_argument('--polymerase', choices=['phi29', 'equiphi29', 'bst', 'klenow'],
                       help='Polymerase type: phi29 (30-40C), equiphi29 (42-45C), bst (60-65C), klenow (25-40C)')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Verbose output')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='Minimal output')

    # GPU acceleration (Category 1 orphaned feature)
    parser.add_argument('--use-gpu', action='store_true',
                       help='Use GPU acceleration (requires CuPy). Provides 10-100x speedup for large genomes')
    parser.add_argument('--no-gpu', action='store_true',
                       help='Disable auto-GPU detection (GPU is auto-enabled when CuPy is available)')
    parser.add_argument('--gpu-device', type=int, default=0,
                       help='GPU device ID to use (default: 0)')

    # Quality assurance (Category 3)
    parser.add_argument('--enable-qa', action='store_true',
                       help='Enable quality assurance checks at each step')


def run_step1(args):
    """Run step 1: K-mer preprocessing"""
    check_jellyfish_available()
    validate_params_json_file(args.json_file)
    import time as _time
    _t0 = _time.time()
    logger.info("Step 1: K-mer preprocessing starting...")
    try:
        from neoswga.core import pipeline, parameter

        # Set json_file if provided (needed for _initialize to load params)
        merge_args_to_parameter(args, parameter, ['json_file'])

        # Apply CLI k-mer range arguments BEFORE initialization
        # This ensures they're available when _initialize() reads defaults
        merge_args_to_parameter(args, parameter, ['min_k', 'max_k'])

        # Initialize pipeline (lazy init, loads params from JSON)
        pipeline._initialize()

        # Re-apply CLI arguments after initialization to ensure they take precedence
        # (initialization may have overwritten them with JSON values)
        merge_args_to_parameter(args, parameter, ['min_k', 'max_k'])

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration (boolean flag)
        if hasattr(args, 'enable_qa') and args.enable_qa:
            parameter.enable_qa = True

        # Log effective parameters
        if not args.quiet:
            min_k = getattr(parameter, 'min_k', 6)
            max_k = getattr(parameter, 'max_k', 12)
            logger.info(f"K-mer range: {min_k}-{max_k}bp")

        # Run step1 (QA integration is only available for step2-4)
        if getattr(parameter, 'enable_qa', False):
            logger.info("Note: QA integration is available for filter/score/optimize steps, not count-kmers")
        pipeline.step1()

        # Count k-mers for exclusion genome if specified via CLI
        excl_genome = getattr(args, 'exclusion_genome', None)
        if excl_genome:
            if not os.path.exists(excl_genome):
                logger.error(f"Exclusion genome not found: {excl_genome}")
                sys.exit(1)
            from neoswga.core.kmer_counter import run_jellyfish as _run_jf
            excl_name = os.path.splitext(os.path.basename(excl_genome))[0]
            excl_prefix = os.path.join(parameter.data_dir, f"excl_{excl_name}")
            min_k = getattr(parameter, 'min_k', 6)
            max_k = getattr(parameter, 'max_k', 12)
            logger.info(f"Counting k-mers in exclusion genome: {excl_genome}")
            _run_jf(excl_genome, excl_prefix, min_k, max_k)
            # Store exclusion prefix for downstream steps
            parameter.excl_genomes = [excl_genome]
            parameter.excl_prefixes = [excl_prefix]

        # Count k-mers for blacklist genome(s) if specified via CLI
        bl_genomes_cli = getattr(args, 'blacklist', None)
        if bl_genomes_cli:
            from neoswga.core.kmer_counter import run_jellyfish as _run_jf_bl
            for bl_genome in bl_genomes_cli:
                if not os.path.exists(bl_genome):
                    logger.error(f"Blacklist genome not found: {bl_genome}")
                    sys.exit(1)
            bl_prefixes_cli = []
            min_k = getattr(parameter, 'min_k', 6)
            max_k = getattr(parameter, 'max_k', 12)
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
        from neoswga.core import pipeline, parameter

        # Set json_file if provided (needed for _initialize to load params)
        if hasattr(args, 'json_file') and args.json_file:
            parameter.json_file = args.json_file

        # Initialize pipeline (lazy init, loads params from JSON)
        pipeline._initialize()

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration
        if hasattr(args, 'enable_qa') and args.enable_qa:
            parameter.enable_qa = True

        # Apply preset if specified
        if hasattr(args, 'preset') and args.preset:
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
        merge_args_to_parameter(args, parameter, ['gc_min', 'gc_max'])

        # Special handling for gc_tolerance (calculates gc_min/gc_max from genome GC)
        if hasattr(args, 'gc_tolerance') and args.gc_tolerance is not None:
            if (not hasattr(args, 'gc_min') or args.gc_min is None) and hasattr(parameter, 'genome_gc') and parameter.genome_gc:
                parameter.gc_min = max(0.20, parameter.genome_gc - args.gc_tolerance)
                parameter.gc_max = min(0.80, parameter.genome_gc + args.gc_tolerance)

        # Reaction conditions
        merge_args_to_parameter(args, parameter, ['reaction_temp', 'na_conc', 'mg_conc'])

        # Additives (with bsa -> bsa_ug_ml mapping)
        merge_args_to_parameter(args, parameter, [
            'dmso_percent', 'betaine_m', 'trehalose_m', 'glycerol_percent',
            'peg_percent', 'ethanol_percent', 'urea_m', 'tmac_m', 'formamide_percent'
        ])
        merge_args_to_parameter(args, parameter, ['bsa'], {'bsa': 'bsa_ug_ml'})

        # Boolean flag for SSB (special handling)
        if hasattr(args, 'ssb') and args.ssb:
            parameter.ssb = True

        # Traditional filtering parameters
        merge_args_to_parameter(args, parameter, [
            'min_fg_freq', 'max_bg_freq', 'max_gini', 'max_primer',
            'min_tm', 'max_tm', 'max_dimer_bp', 'max_self_dimer_bp'
        ])

        # Polymerase
        merge_args_to_parameter(args, parameter, ['polymerase'])

        # Bloom filter for large background genomes
        use_bloom = getattr(args, 'use_bloom_filter', False)
        bloom_path = getattr(args, 'bloom_filter_path', None)
        sampled_path = getattr(args, 'sampled_index_path', None)
        bloom_max_bg = getattr(args, 'bloom_max_bg_matches', 10)

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
        excl_genome = getattr(args, 'exclusion_genome', None)
        if excl_genome:
            if not os.path.exists(excl_genome):
                logger.error(f"Exclusion genome not found: {excl_genome}")
                sys.exit(1)
            excl_name = os.path.splitext(os.path.basename(excl_genome))[0]
            excl_prefix = os.path.join(parameter.data_dir, f"excl_{excl_name}")
            parameter.excl_genomes = [excl_genome]
            parameter.excl_prefixes = [excl_prefix]
            excl_threshold = getattr(args, 'excl_threshold', 0)
            if excl_threshold is not None:
                parameter.excl_threshold = excl_threshold
            if not args.quiet:
                logger.info(f"Exclusion genome: {excl_genome}")
                logger.info(f"  Threshold: {parameter.excl_threshold} (0 = reject any hit)")

        # Handle blacklist genome(s) from CLI
        bl_genomes_cli = getattr(args, 'blacklist', None)
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
                fname_genomes=bl_genomes_cli, cpus=getattr(parameter, 'cpus', 1))
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
            if hasattr(parameter, 'gc_min'):
                logger.info(f"  GC range: {parameter.gc_min:.3f} - {parameter.gc_max:.3f}")
            if hasattr(parameter, 'reaction_temp'):
                logger.info(f"  Reaction temp: {parameter.reaction_temp}C")
            if hasattr(parameter, 'polymerase'):
                logger.info(f"  Polymerase: {parameter.polymerase}")
            if hasattr(parameter, 'dmso_percent') and parameter.dmso_percent > 0:
                logger.info(f"  DMSO: {parameter.dmso_percent}%")
            if hasattr(parameter, 'betaine_m') and parameter.betaine_m > 0:
                logger.info(f"  Betaine: {parameter.betaine_m} M")

        # Run step2 with QA if enabled
        if getattr(parameter, 'enable_qa', False):
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
    _elapsed = _time.time() - _t0
    logger.info(f"Step 2 complete in {_elapsed:.1f}s")


def run_step3(args):
    """Run step 3: Random forest scoring"""
    validate_params_json_file(args.json_file)
    import time as _time
    _t0 = _time.time()
    try:
        from neoswga.core import pipeline, parameter

        # Set json_file and initialize pipeline
        merge_args_to_parameter(args, parameter, ['json_file'])
        pipeline._initialize()

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration (boolean flag)
        if hasattr(args, 'enable_qa') and args.enable_qa:
            parameter.enable_qa = True

        # Apply CLI arguments to parameter module
        merge_args_to_parameter(args, parameter, ['min_amp_pred'])

        # Log effective parameters
        if not args.quiet:
            min_amp_pred = getattr(parameter, 'min_amp_pred', None) or 'auto'
            logger.info(f"Minimum amplification prediction score: {min_amp_pred}")

        # Enhanced feature engineering
        use_enhanced = getattr(args, 'use_enhanced_features', False)
        enhanced_model_path = getattr(args, 'enhanced_model_path', None)

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

        # Fast scoring mode: skip expensive delta-G histogram features
        if getattr(args, 'fast_score', False):
            parameter.fast_score = True
            logger.info("Fast scoring: skipping thermodynamic histogram features (~100x speedup)")

        # Run step3 with QA if enabled
        if getattr(parameter, 'enable_qa', False):
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
    _elapsed = _time.time() - _t0
    logger.info(f"Step 3 complete in {_elapsed:.1f}s")


def print_method_guide():
    """Print optimization method selection guide."""
    guide = """
OPTIMIZATION METHOD SELECTION GUIDE
====================================

Decision Tree:
  Speed critical?           -> dominating-set
  Need dimer-free?          -> clique (or coverage-then-dimerfree)
  Clinical/diagnostic?      -> background-aware
  GC-rich genome (>65%)?    -> equiphi29
  Exact solution needed?    -> milp (pools <100)
  Want coverage + no dimers -> coverage-then-dimerfree
  Default                   -> hybrid

Single-Stage Methods:
+------------------+--------+-----------+-------------+------------------------+
| Method           | Speed  | Coverage  | Specificity | Best For               |
+------------------+--------+-----------+-------------+------------------------+
| hybrid           | Medium | Excellent | Good        | General use (default)  |
| dominating-set   | Fast   | Excellent | Fair        | Large pools, quick     |
| background-aware | Slow   | Good      | Excellent   | Clinical, low bg       |
| clique           | Var.   | Good      | Good        | Dimer-free guarantee   |
| equiphi29        | Medium | Good      | Good        | GC-rich, 42-45C        |
| normalized       | Medium | Varies    | Varies      | Configurable weights   |
| network          | Medium | Good      | Good        | Tm-balanced sets       |
| genetic          | Slow   | Good      | Good        | Multi-objective        |
| moea             | Slow   | Good      | Good        | Pareto optimization    |
| milp             | Var.   | Optimal   | Good        | Exact (small sets)     |
| greedy           | Fast   | Fair      | Fair        | Simple baseline        |
| tiling           | Fast   | Excellent | Fair        | Uniform spacing        |
+------------------+--------+-----------+-------------+------------------------+

Pipeline Methods (serial cascades):
+-------------------------+----------+--------------------------------------+
| Method                  | Speed    | Description                          |
+-------------------------+----------+--------------------------------------+
| coverage-then-dimerfree | Moderate | DS -> Clique (coverage + dimer-free) |
| dimerfree-scored        | Moderate | Clique -> Network (dimer-free +      |
|                         |          | connectivity scoring)                |
| bg-prefilter            | Moderate | Background pruning pre-filter alone  |
| bg-prefilter-hybrid     | Moderate | BG pruning -> Hybrid (bg reduction   |
|                         |          | + general optimization)              |
+-------------------------+----------+--------------------------------------+

Strategy Presets (for --optimization-method=normalized):
+-------------+--------------------------------------------------+
| Strategy    | Description                                      |
+-------------+--------------------------------------------------+
| clinical    | High specificity (40% background weight)         |
| discovery   | Max coverage (40% coverage weight)               |
| fast        | Quick screening (no background penalty)          |
| balanced    | Equal weights across all objectives              |
| enrichment  | Sequencing enrichment (balanced coverage/amp)    |
+-------------+--------------------------------------------------+

Usage Examples:
  # Default (hybrid)
  neoswga optimize -j params.json

  # Fast screening
  neoswga optimize -j params.json --optimization-method=dominating-set

  # Clinical samples (low background)
  neoswga optimize -j params.json --optimization-method=background-aware

  # Guaranteed dimer-free sets (pools <200 primers)
  neoswga optimize -j params.json --optimization-method=clique

  # Coverage-first, then dimer-free filtering
  neoswga optimize -j params.json --optimization-method=coverage-then-dimerfree

  # Dimer-free with network connectivity scoring
  neoswga optimize -j params.json --optimization-method=dimerfree-scored

  # Strategy-based with clinical preset
  neoswga optimize -j params.json --optimization-method=normalized --strategy=clinical

For detailed documentation: docs/optimization_guide.md
"""
    print(guide)


def run_step4(args):
    """Run step 4: Primer set optimization (network-based + experimental)"""
    validate_params_json_file(args.json_file)
    import time as _time
    _t0 = _time.time()

    # Handle --method-guide option
    if getattr(args, 'method_guide', False):
        print_method_guide()
        return

    logger.info("Primer set optimization")

    try:
        from neoswga.core import parameter

        # Set json_file if provided
        merge_args_to_parameter(args, parameter, ['json_file'])

        # Pass polymerase info to optimizer (use adapted value, not raw JSON,
        # so GC-adaptive strategy recommendations are respected)
        polymerase = getattr(parameter, 'polymerase', 'phi29')
        if polymerase != 'phi29':
            logger.info(f"Polymerase: {polymerase} (config applied to optimizer)")

        logger.info(f"Optimization method: {args.optimization_method}")
        logger.info(f"Position cache: {args.use_position_cache}")
        logger.info(f"Background filter: {args.use_background_filter}")

        # GPU acceleration (with auto-detection)
        setup_gpu_acceleration(args, parameter, quiet=args.quiet)

        # QA integration (boolean flag)
        if hasattr(args, 'enable_qa') and args.enable_qa:
            parameter.enable_qa = True

        # Set num_primers from CLI or JSON file (special handling needed)
        if hasattr(args, 'num_primers') and args.num_primers:
            parameter.num_primers = args.num_primers
            parameter.target_set_size = args.num_primers
            logger.info(f"Target primer count: {args.num_primers}")
        elif not hasattr(parameter, 'num_primers') or not parameter.num_primers:
            # Try to read from JSON if loaded
            json_data = getattr(parameter, '_json_data', {})
            num_primers = json_data.get('num_primers', json_data.get('target_set_size', 6))
            parameter.num_primers = num_primers
            parameter.target_set_size = num_primers

        # Automatic set size optimization (overrides num_primers if enabled)
        auto_size = getattr(args, 'auto_size', False)
        application = getattr(args, 'application', 'enrichment')
        min_fg_bg_ratio = getattr(args, 'min_fg_bg_ratio', None)
        show_frontier = getattr(args, 'show_frontier', False)
        quick_estimate = getattr(args, 'quick_estimate', False)

        if auto_size:
            logger.info("=" * 60)
            logger.info("Automatic Set Size Optimization")
            logger.info("=" * 60)

            try:
                from neoswga.core.set_size_optimizer import recommend_set_size, create_baseline_effects
                from neoswga.core.mechanistic_model import MechanisticModel
                from neoswga.core.reaction_conditions import ReactionConditions

                # Get genome information
                json_data = getattr(parameter, '_json_data', {})
                fg_seq_lengths = getattr(parameter, 'fg_seq_lengths', [])
                genome_length = sum(fg_seq_lengths) if fg_seq_lengths else json_data.get('fg_genome_length', 1_000_000)
                primer_length = json_data.get('max_k', 12)

                # Get template GC content
                template_gc = getattr(args, 'template_gc', None)
                if template_gc is None:
                    # Try to calculate from genome or use default
                    template_gc = json_data.get('fg_gc', 0.5)

                # Get polymerase and create conditions
                polymerase = json_data.get('polymerase', 'phi29')
                reaction_temp = json_data.get('reaction_temp', 30.0 if polymerase == 'phi29' else 42.0)

                try:
                    conditions = ReactionConditions(
                        temp=reaction_temp,
                        polymerase=polymerase,
                        mg_conc=json_data.get('mg_conc', 2.5),
                        dmso_percent=json_data.get('dmso_percent', 0.0),
                        betaine_m=json_data.get('betaine_m', 0.0),
                    )
                    model = MechanisticModel(conditions)
                    sample_primer = 'A' * primer_length  # Neutral sequence
                    mech_effects = model.calculate_effects(sample_primer, template_gc)
                except Exception as e:
                    logger.warning(f"Could not create mechanistic model: {e}")
                    mech_effects = create_baseline_effects()

                # Get processivity based on polymerase
                processivity_map = {'phi29': 70000, 'equiphi29': 80000, 'bst': 2000, 'klenow': 10000}
                processivity = processivity_map.get(polymerase, 70000)

                # Get recommendation (supports new min_fg_bg_ratio parameter)
                recommendation = recommend_set_size(
                    application=application,
                    genome_length=genome_length,
                    primer_length=primer_length,
                    mech_effects=mech_effects,
                    processivity=processivity,
                    min_fg_bg_ratio=min_fg_bg_ratio,
                )

                auto_num_primers = recommendation['recommended_size']
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
                logger.info(f"  (typical range: {recommendation['size_range'][0]}-{recommendation['size_range'][1]})")
                logger.info("=" * 60)

            except Exception as e:
                logger.warning(f"Auto-size failed: {e}. Using default num_primers.")

        # Cooperative binding (experimental - not yet fully integrated)
        if hasattr(args, 'use_cooperative_binding') and args.use_cooperative_binding:
            logger.warning("--use-cooperative-binding is experimental and not yet fully integrated")
            parameter.use_cooperative_binding = True

        # Mechanistic model (not yet integrated into optimizer scoring)
        if getattr(args, 'use_mechanistic_model', False):
            logger.warning(
                "--use-mechanistic-model is not yet integrated into the optimization step. "
                "The mechanistic model is used in simulation and auto-size features."
            )

        # Primer strategy
        merge_args_to_parameter(args, parameter, ['primer_strategy'])
        if hasattr(args, 'primer_strategy') and args.primer_strategy == 'hybrid':
            logger.info("Using hybrid primer strategy (mixed lengths)")

        # Use unified optimizer framework (all methods handled via factory pattern)
        from neoswga.core.unified_optimizer import optimize_step4, list_available_optimizers

        # Get uniformity weight from args (default 0.0 for backward compatibility)
        uniformity_weight = getattr(args, 'uniformity_weight', 0.0)
        if uniformity_weight > 0:
            logger.info(f"Uniformity weight: {uniformity_weight:.2f}")

        # Get minimal primer selection options
        minimize_primers = getattr(args, 'minimize_primers', False)
        target_coverage = getattr(args, 'target_coverage', 0.70)
        if minimize_primers:
            logger.info(f"Minimal primer selection enabled (target coverage: {target_coverage:.1%})")

        # Get strategy for normalized optimizer
        strategy = getattr(args, 'strategy', 'balanced')
        if args.optimization_method == 'normalized':
            logger.info(f"Using normalized optimizer with strategy: {strategy}")

        # Background pre-filter flag (enabled by default, --no-bg-prefilter disables)
        bg_prefilter = not getattr(args, 'no_bg_prefilter', False)

        # Host-free mode: optimize without background genome data
        no_background = getattr(args, 'no_background', False)
        if no_background:
            logger.info("Host-free mode: background genome data will be ignored")

        # Run unified optimization
        # Random seed for reproducibility (stochastic optimizers)
        seed = getattr(args, 'seed', None)
        if seed is not None:
            logger.info(f"Random seed: {seed}")

        # Pass explicit target_size from parameter module to avoid re-initialization override
        target_size = getattr(parameter, 'target_set_size',
                              getattr(parameter, 'num_primers', 6))

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
        )

        if results:
            target_size = getattr(parameter, 'num_primers', 6)
            target_size = getattr(parameter, 'target_set_size', target_size)
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
                    f"  - Trying a different optimizer "
                    f"(current: {args.optimization_method})"
                )
        else:
            logger.error("No primer sets found. Optimization failed.")
            sys.exit(1)

        # Show Pareto frontier analysis (optional)
        if show_frontier and results and cache is not None:
            try:
                import pandas as pd
                from neoswga.core.set_size_optimizer import (
                    ParetoFrontierGenerator, select_from_frontier
                )
                from neoswga.core.pareto_frontier import (
                    generate_frontier_report, summarize_frontier_for_cli, plot_frontier
                )

                logger.info("")
                logger.info("=" * 60)
                logger.info("Pareto Frontier Analysis")
                logger.info("=" * 60)

                # Load step2 or step3 DataFrame for primer pool
                data_dir = parameter.data_dir
                step3_file = os.path.join(data_dir, 'step3_df.csv')
                step2_file = os.path.join(data_dir, 'step2_df.csv')

                if os.path.exists(step3_file):
                    primer_pool = pd.read_csv(step3_file)
                elif os.path.exists(step2_file):
                    primer_pool = pd.read_csv(step2_file)
                else:
                    raise FileNotFoundError("No primer pool CSV found")

                # Get genome lengths
                fg_lengths = getattr(parameter, 'fg_lengths', [1_000_000])
                bg_lengths = getattr(parameter, 'bg_lengths', [])
                fg_prefixes = parameter.fg_prefixes
                bg_prefixes = getattr(parameter, 'bg_prefixes', [])

                # Create frontier generator
                generator = ParetoFrontierGenerator(
                    primer_pool=primer_pool,
                    position_cache=cache,
                    fg_prefixes=fg_prefixes,
                    bg_prefixes=bg_prefixes,
                    fg_seq_lengths=fg_lengths,
                    bg_seq_lengths=bg_lengths,
                    processivity=70000,
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
                    plot_path = os.path.join(data_dir, 'pareto_frontier.png')
                    fig.savefig(plot_path, dpi=150, bbox_inches='tight')
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
        validate_simulation = getattr(args, 'validate_simulation', False)
        if validate_simulation and results:
            simulation_time = getattr(args, 'simulation_time', 3600.0)
            logger.info("=" * 60)
            logger.info("Stochastic Validation (Gillespie Algorithm)")
            logger.info("=" * 60)
            logger.info(f"Simulation time: {simulation_time/3600:.1f} hours")

            try:
                from neoswga.core.stochastic_simulator import validate_network_predictions
                from neoswga.core.amplicon_network import AmpliconNetwork

                primers = results[0]

                # Build simplified networks for validation
                # Get genome lengths and positions from cache
                fg_prefixes = parameter.fg_prefixes
                bg_prefixes = getattr(parameter, 'bg_prefixes', [])

                if cache is not None:
                    # Build AmpliconNetwork for target
                    fg_genome_length = sum(getattr(parameter, 'fg_lengths', [1000000]))
                    bg_genome_length = sum(getattr(parameter, 'bg_lengths', [1000000]))

                    # Create networks from cache
                    from neoswga.core.network_optimizer import AmplificationNetwork

                    # Build fg network
                    fg_network = AmplificationNetwork(max_extension=70000)
                    for primer in primers:
                        for prefix in fg_prefixes:
                            try:
                                fw_pos = cache.get_positions(prefix, primer, 'forward')
                                rv_pos = cache.get_positions(prefix, primer, 'reverse')
                                for pos in fw_pos:
                                    fg_network.add_site(pos, primer, 'forward')
                                for pos in rv_pos:
                                    fg_network.add_site(pos, primer, 'reverse')
                            except Exception:
                                pass
                    fg_network.build_graph()

                    # Build bg network
                    bg_network = AmplificationNetwork(max_extension=70000)
                    for primer in primers:
                        for prefix in bg_prefixes:
                            try:
                                fw_pos = cache.get_positions(prefix, primer, 'forward')
                                rv_pos = cache.get_positions(prefix, primer, 'reverse')
                                for pos in fw_pos:
                                    bg_network.add_site(pos, primer, 'forward')
                                for pos in rv_pos:
                                    bg_network.add_site(pos, primer, 'reverse')
                            except Exception:
                                pass
                    bg_network.build_graph()

                    # Run validation
                    validation_results = validate_network_predictions(
                        primers, fg_network, bg_network
                    )

                    # Report results
                    predicted = validation_results['predicted']
                    simulated = validation_results['simulated']
                    validation = validation_results['validation']

                    logger.info("\nPrediction vs Simulation:")
                    logger.info(f"  Predicted enrichment: {predicted['enrichment']:.0f}x")
                    logger.info(f"  Simulated enrichment: {simulated['enrichment']:.0f}x")
                    logger.info(f"  Prediction error: {validation['prediction_error']:.1%}")

                    if validation['prediction_accurate']:
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

        _elapsed = _time.time() - _t0
        logger.info(f"Step 4 complete in {_elapsed:.1f}s")
        if not args.quiet:
            data_dir = getattr(parameter, 'data_dir', '.')
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
    from neoswga.core.background_filter import BackgroundBloomFilter, SampledGenomeIndex
    import pickle

    logger.info("Building background filter")
    logger.info(f"Input: {args.genome}")
    logger.info(f"Output: {args.output_dir}")

    os.makedirs(args.output_dir, exist_ok=True)

    bloom_path = os.path.join(args.output_dir, 'bg_bloom.pkl')
    sampled_path = os.path.join(args.output_dir, 'bg_sampled.pkl')

    if not args.force and os.path.exists(bloom_path):
        logger.warning("Filter already exists. Use --force to rebuild.")
        return

    try:
        min_k = getattr(args, 'min_k', 6)
        max_k = getattr(args, 'max_k', 12)

        if getattr(args, 'from_kmers', False):
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
            sampled = SampledGenomeIndex(sample_rate=1)  # rate=1 since k-mer files are already unique
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


def run_validate(args):
    """Run validation tests"""
    from neoswga.core.validation import ValidationSuite, quick_validation
    from neoswga.core.gpu_acceleration import is_gpu_available, get_gpu_info
    from neoswga.core.kmer_counter import check_jellyfish_available

    # Show system capabilities
    logger.info("System Capabilities:")
    logger.info(f"  Jellyfish: {'Available' if check_jellyfish_available() else 'NOT FOUND'}")
    if is_gpu_available():
        gpu_info = get_gpu_info()
        logger.info(f"  GPU: {gpu_info.get('device_name', 'Available')} ({gpu_info.get('memory_total_gb', 0):.1f} GB)")
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
            bench_path = os.path.join(os.path.dirname(__file__), '..', 'scripts', 'benchmarking', 'benchmark_improvements.py')
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


def run_init(args):
    """Run the setup wizard to create params.json"""
    from neoswga.core.wizard import run_wizard
    import inspect

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
        if getattr(args, 'blacklist', None):
            sig = inspect.signature(run_wizard)
            if 'blacklist_paths' in sig.parameters:
                wizard_kwargs['blacklist_paths'] = args.blacklist
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
    """Validate params.json configuration"""
    from neoswga.core.param_validator import validate_params_file

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
    if getattr(args, 'output', None):
        with open(args.output, 'w', encoding='utf-8') as fh:
            fh.write(output + "\n")
        logger.info(f"Schema written to {args.output}")
    else:
        print(output)


def run_validate_model(args):
    """Validate mechanistic model against expected behavior."""
    from neoswga.core.model_validation import (
        validate_mechanistic_model,
        format_validation_report,
    )
    import json

    logger.info("Validating mechanistic model...")
    results = validate_mechanistic_model()

    if getattr(args, 'output_json', False):
        # Output as JSON
        # Convert results to JSON-serializable format
        json_results = []
        for r in results:
            jr = {
                'test': r.get('test', 'Unknown'),
                'passed': r.get('passed', False),
                'summary': r.get('summary', ''),
            }
            if 'error' in r:
                jr['error'] = r['error']
            json_results.append(jr)
        print(json.dumps(json_results, indent=2))
    else:
        # Output as formatted report
        report = format_validation_report(results)
        print(report)

    # Exit with appropriate code
    all_passed = all(r.get('passed', False) for r in results)
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
        validate_results_directory,
        validate_metrics,
        ValidationLevel,
    )

    results_dir = Path(args.dir)
    quiet = getattr(args, 'quiet', False)
    check_only = getattr(args, 'check', False)
    interactive = getattr(args, 'interactive', False)

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
    level = getattr(args, 'level', 'summary')
    if args.output:
        output_path = args.output
    else:
        if level == 'full':
            output_path = str(results_dir / 'technical_report.html')
        else:
            output_path = str(results_dir / 'report.html')

    try:
        if args.format == 'json':
            # JSON export
            from neoswga.core.report.metrics import collect_pipeline_metrics
            from neoswga.core.report.quality import calculate_quality_grade
            import json

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
                'grade': quality.grade.value,
                'composite_score': quality.composite_score,
                'recommendation': quality.recommendation,
                'recommendation_details': quality.recommendation_details,
                'considerations': quality.considerations,
                'primer_count': metrics.primer_count,
                'components': [
                    {
                        'name': c.name,
                        'weight': c.weight,
                        'raw_value': c.raw_value,
                        'normalized_score': c.normalized_score,
                        'rating': c.rating,
                    }
                    for c in quality.components
                ],
            }

            # Change extension if needed
            if not output_path.endswith('.json'):
                output_path = output_path.rsplit('.', 1)[0] + '.json'

            with open(output_path, 'w') as f:
                json.dump(report_data, f, indent=2)

            print(f"\nQuality Grade: {quality.grade.value} "
                  f"({quality.composite_score:.2f}/1.00)")
            print(f"Recommendation: {quality.recommendation}")
            print(f"\nJSON report saved to: {output_path}")

        elif level == 'full':
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
            print(f"\nQuality Grade: {data.quality.grade.value} "
                  f"({data.quality.composite_score:.2f}/1.00)")
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
            print(f"\nQuality Grade: {summary.quality.grade.value} "
                  f"({summary.quality.composite_score:.2f}/1.00)")
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
    import h5py
    import glob as glob_module
    from neoswga.core.thermodynamics import reverse_complement

    positions = {}
    h5_files = glob_module.glob(os.path.join(results_dir, '*_positions.h5'))

    if not h5_files:
        params_path = os.path.join(results_dir, 'params.json')
        if os.path.exists(params_path):
            with open(params_path) as f:
                params = json.load(f)
            fg_prefixes = params.get('fg_prefixes', [])
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
                with h5py.File(h5_path, 'r') as db:
                    if primer in db:
                        for pos in db[primer][:]:
                            sites.append((int(pos), 'forward'))
                    rc = rc_map[primer]
                    if rc in db:
                        for pos in db[rc][:]:
                            sites.append((int(pos), 'reverse'))
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
    params_path = os.path.join(results_dir, 'params.json')
    if os.path.exists(params_path):
        with open(params_path) as f:
            params = json.load(f)
        lengths = params.get('fg_seq_lengths', [])
        if lengths:
            return sum(lengths)
    logger.warning("Could not determine genome length. Using 0.")
    return 0


def run_export(args):
    """Export primers for synthesis ordering."""
    from neoswga.core.export import PrimerExporter, PrimerModifications
    from pathlib import Path

    try:
        # Determine modification profile
        if getattr(args, 'no_modifications', False):
            mods = PrimerModifications.from_profile('none')
        else:
            mods = PrimerModifications.from_profile(args.modifications)
            # Override PTO bonds if specified
            if args.pto_bonds is not None:
                mods.pto_bonds = args.pto_bonds

        # Load exporter from results
        exporter = PrimerExporter.from_results_dir(
            args.dir,
            params_file=getattr(args, 'json_file', None)
        )
        exporter.modifications = mods

        # Print summary
        exporter.print_summary()

        output_dir = Path(args.output)

        if args.format == 'all':
            outputs = exporter.export_all(
                str(output_dir),
                project_name=args.project,
                vendors=[args.vendor]
            )
            print("\nExported files:")
            for fmt, path in outputs.items():
                print(f"  {fmt}: {path}")

        elif args.format == 'fasta':
            output_dir.mkdir(parents=True, exist_ok=True)
            fasta_path = output_dir / f"{args.project}_primers.fasta"
            exporter.export_fasta(str(fasta_path), prefix=args.project, include_metadata=True)
            print(f"Exported: {fasta_path}")

        elif args.format == 'csv':
            output_dir.mkdir(parents=True, exist_ok=True)
            csv_path = output_dir / f"{args.project}_order_{args.vendor}.csv"
            exporter.export_vendor_csv(str(csv_path), vendor=args.vendor, project_name=args.project)
            print(f"Exported: {csv_path}")

        elif args.format == 'protocol':
            output_dir.mkdir(parents=True, exist_ok=True)
            protocol_path = output_dir / f"{args.project}_protocol.md"
            exporter.export_protocol(str(protocol_path))
            print(f"Exported: {protocol_path}")

        elif args.format in ('bed', 'bedgraph'):
            output_dir.mkdir(parents=True, exist_ok=True)

            # Load positions from results directory
            positions = _load_primer_positions(args.dir, exporter.primers)

            if args.format == 'bed':
                bed_path = output_dir / f"{args.project}_binding_sites.bed"
                exporter.export_bed(
                    str(bed_path), positions,
                    genome_name=args.genome_name
                )
                print(f"Exported: {bed_path}")
            else:
                genome_length = _get_genome_length(args.dir)
                bg_path = output_dir / f"{args.project}_coverage.bedgraph"
                exporter.export_bedgraph(
                    str(bg_path), positions,
                    genome_name=args.genome_name,
                    genome_length=genome_length,
                    window_size=args.window_size
                )
                print(f"Exported: {bg_path}")

        print("\nPrimers ready for ordering!")

    except FileNotFoundError as e:
        logger.error(str(e))
        logger.error("Make sure to run the full pipeline (count-kmers, filter, score, optimize) first.")
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
    if getattr(args, 'kmer_range', False):
        from neoswga.core.condition_suggester import suggest_kmer_range
        from neoswga.core.genome_io import GenomeLoader

        genome_size = getattr(args, 'genome_size', None)
        if genome_size is None and getattr(args, 'genome', None):
            try:
                sequence = GenomeLoader().load_genome(args.genome, return_stats=False)
                genome_size = len(sequence)
            except Exception as e:
                logger.error(f"Could not read genome size from {args.genome}: {e}")
                sys.exit(1)
        if genome_size is None:
            logger.error("--kmer-range requires --genome-size or --genome")
            sys.exit(1)

        polymerase = getattr(args, 'polymerase', 'phi29')
        min_k, max_k = suggest_kmer_range(
            genome_size, gc=genome_gc, polymerase=polymerase,
        )
        print(f"Recommended k-mer range for {genome_size:,} bp, "
              f"GC={genome_gc if genome_gc is not None else 'unknown'}, "
              f"polymerase={polymerase}:")
        print(f"  min_k = {min_k}")
        print(f"  max_k = {max_k}")
        print(f"\nAdd to params.json:")
        print(f'  "min_k": {min_k},')
        print(f'  "max_k": {max_k}')
        return

    # Condition sweep mode
    if getattr(args, 'sweep', False):
        from neoswga.core.condition_suggester import sweep_conditions
        primer_length = args.primer_length
        if primer_length is None:
            # Default primer length by polymerase type
            polymerase = getattr(args, 'polymerase', 'phi29')
            primer_length_defaults = {
                'phi29': 9, 'equiphi29': 15, 'bst': 20, 'klenow': 12,
            }
            primer_length = primer_length_defaults.get(polymerase, 10)
        if genome_gc is None:
            genome_gc = 0.5
            logger.info("No GC content provided, using default 50%")
        sweep_conditions(
            genome_gc=genome_gc,
            primer_length=primer_length,
            polymerase=getattr(args, 'polymerase', 'phi29'),
            output_path=getattr(args, 'output', None),
            verbose=True,
        )
        return

    if genome_gc is None and args.primer_length is None:
        logger.error("Provide at least --genome-gc or --genome or --primer-length")
        sys.exit(1)

    # Use advanced optimizer if requested or if polymerase/optimize-for specified
    use_optimizer = getattr(args, 'use_optimizer', False)
    polymerase = getattr(args, 'polymerase', 'phi29')
    optimize_for = getattr(args, 'optimize_for', 'amplification')

    if use_optimizer or polymerase != 'phi29' or optimize_for != 'amplification':
        # Use new AdditiveOptimizer
        from neoswga.core.additive_optimizer import AdditiveOptimizer
        from neoswga.core.mechanistic_params import get_polymerase_params

        # Determine primer length if not provided
        primer_length = args.primer_length
        if primer_length is None:
            poly_params = get_polymerase_params(polymerase)
            primer_length = int((poly_params['primer_length_range'][0] +
                                poly_params['primer_length_range'][1]) / 2)
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
            verbose=True
        )


def show_presets():
    """Show available reaction condition presets"""
    print("\n" + "="*70)
    print("Available Reaction Condition Presets")
    print("="*70 + "\n")

    # Show presets from the PRESETS dict
    for name, config in PRESETS.items():
        print(f"{name}:")
        print(f"  Temperature: {config['temperature']}C")
        print(f"  Polymerase: {config['polymerase']}")
        print(f"  DMSO: {config['dmso_percent']}%")
        print(f"  Betaine: {config['betaine_m']} M")
        print(f"  SSB: {config['ssb']}")
        print(f"  Optimization: {config['optimization_method']}")
        print()

    # Show additional presets available via --preset but not in PRESETS dict
    additional = ['q_solution', 'gc_melt', 'crude_sample', 'low_temp', 'bst', 'klenow', 'extreme_gc']
    try:
        from neoswga.core import reaction_conditions as rc
        preset_funcs = {
            'q_solution': ('Q-solution equivalent', rc.get_q_solution_equivalent),
            'gc_melt': ('GC-melt conditions', rc.get_gc_melt_conditions),
            'crude_sample': ('Crude sample conditions', rc.get_crude_sample_conditions),
            'low_temp': ('Low temperature conditions', rc.get_low_temp_conditions),
            'bst': ('Bst polymerase conditions', rc.get_bst_conditions),
            'klenow': ('Klenow polymerase conditions', rc.get_klenow_conditions),
            'extreme_gc': ('Extreme GC content conditions', rc.get_extreme_gc_conditions),
        }
        for name in additional:
            desc, func = preset_funcs[name]
            try:
                cond = func()
                print(f"{name}:")
                print(f"  {desc}")
                print(f"  Temperature: {cond.temp}C, Na+: {cond.na_conc} mM, Mg2+: {cond.mg_conc} mM")
                if hasattr(cond, 'dmso_percent') and cond.dmso_percent:
                    print(f"  DMSO: {cond.dmso_percent}%")
                if hasattr(cond, 'betaine_m') and cond.betaine_m:
                    print(f"  Betaine: {cond.betaine_m} M")
                print()
            except Exception:
                print(f"{name}: {desc}")
                print()
    except ImportError:
        for name in additional:
            print(f"{name}: (use --preset {name} with filter command)")
        print()


def optimize_conditions(args):
    """Optimize reaction conditions for a genome"""
    from neoswga.core import reaction_conditions as rc
    from Bio import SeqIO

    logger.info("Analyzing genome and recommending optimal conditions...")

    # Load genome
    records = list(SeqIO.parse(args.fg, 'fasta'))
    genome_seq = str(records[0].seq).upper()
    genome_length = len(genome_seq)

    # Calculate GC content
    gc_content = (genome_seq.count('G') + genome_seq.count('C')) / genome_length

    print(f"\nGenome properties:")
    print(f"  Length: {genome_length:,} bp")
    print(f"  GC content: {gc_content:.1%}")

    # Recommend conditions
    recommendations = rc.recommend_conditions(genome_seq, target_k=args.target_k)

    print(f"\nRecommended conditions:")
    print(f"  Optimal k-mer length: {recommendations['optimal_k']}")
    print(f"  Recommended temperature: {recommendations['temperature']:.1f}C")
    print(f"  Polymerase: {recommendations['polymerase']}")
    print(f"  DMSO: {recommendations['dmso_percent']:.1f}%")
    print(f"  Betaine: {recommendations['betaine_m']:.1f} M")
    if recommendations['ssb']:
        print(f"  SSB: Recommended")

    # Save recommendations
    os.makedirs(args.output, exist_ok=True)
    rec_file = os.path.join(args.output, 'recommended_conditions.json')
    with open(rec_file, 'w') as f:
        json.dump(recommendations, f, indent=2)

    print(f"\nRecommendations saved to: {rec_file}")


def analyze_primer_set(args):
    """Analyze an existing primer set"""
    from neoswga.core import reaction_conditions as rc
    from neoswga.core import thermodynamics as thermo
    from neoswga.core import secondary_structure as ss

    logger.info("Analyzing primer set...")

    # Load preset config
    config_dict = PRESETS[args.preset]
    conditions = rc.ReactionConditions(
        temp=config_dict['temperature'],
        dmso_percent=config_dict['dmso_percent'],
        betaine_m=config_dict['betaine_m'],
        polymerase=config_dict['polymerase']
    )

    primers = args.primers

    print(f"\nAnalyzing {len(primers)} primers:")
    for primer in primers:
        print(f"  {primer}")
    print()

    # Thermodynamic analysis
    print("Thermodynamic Properties:")
    for primer in primers:
        tm = conditions.calculate_effective_tm(primer)
        dg = thermo.calculate_free_energy(primer, conditions.temp)
        gc = thermo.gc_content(primer)
        print(f"  {primer}:")
        print(f"    Tm: {tm:.1f}C")
        print(f"    dG: {dg:.2f} kcal/mol")
        print(f"    GC: {gc:.1%}")

    # Secondary structure analysis
    print("\nSecondary Structure Analysis:")
    for primer in primers:
        hairpins = ss.check_hairpins(primer, conditions)
        homodimer = ss.check_homodimer(primer, conditions)
        print(f"  {primer}:")
        if hairpins:
            worst = min(hairpins, key=lambda h: h.energy)
            print(f"    Hairpin: dG={worst.energy:.2f} kcal/mol, stem={worst.stem_length}bp")
        else:
            print(f"    Hairpin: none detected")
        print(f"    Homodimer: dG={homodimer.energy:.2f} kcal/mol")

    # Check heterodimers
    print("\nHeterodimer Analysis:")
    for i, p1 in enumerate(primers):
        for j, p2 in enumerate(primers):
            if i < j:
                result = ss.check_heterodimer(p1, p2, conditions)
                if result.energy < -6.0:
                    print(f"  {p1} x {p2}: dG={result.energy:.2f} kcal/mol (warning)")

    print("\nAnalysis complete!")


# =========================================================================
# CATEGORY 1: Handler functions for orphaned features
# =========================================================================

def run_analyze_genome(args):
    """Analyze genome suitability for SWGA"""
    from neoswga.core import genome_analysis

    logger.info(f"Analyzing genome: {args.genome}")

    try:
        results = genome_analysis.analyze_genome(
            genome_path=args.genome,
            window_size=args.window_size,
            output_dir=args.output
        )

        logger.info(f"Analysis complete! Report saved to: {args.output}")

    except Exception as e:
        logger.error(f"Genome analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_analyze_dimers(args):
    """Analyze primer dimer interaction network"""
    from neoswga.core.dimer_network_analyzer import DimerNetworkAnalyzer

    # Validate primer sequences
    for primer in args.primers:
        validate_primer_sequence(primer, name="--primers")

    logger.info(f"Analyzing dimer network for {len(args.primers)} primers")

    try:
        analyzer = DimerNetworkAnalyzer(severity_threshold=args.threshold)
        metrics, profiles, matrix = analyzer.analyze_primer_set(args.primers)

        os.makedirs(args.output, exist_ok=True)

        # Write summary
        summary_path = os.path.join(args.output, 'dimer_summary.txt')
        with open(summary_path, 'w') as f:
            f.write(str(metrics) + '\n')
            for primer, profile in profiles.items():
                f.write(f"\n{primer}: {profile}\n")

        logger.info(f"Dimer analysis complete! Results saved to: {args.output}")

        if getattr(metrics, 'num_hub_primers', 0) > 0:
            logger.warning(f"Found {metrics.num_hub_primers} hub primers with many interactions")

    except Exception as e:
        logger.error(f"Dimer analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_analyze_stability(args):
    """Analyze 3' end stability and specificity"""
    from neoswga.core.three_prime_stability import ThreePrimeStabilityAnalyzer

    # Validate primer sequences
    for primer in args.primers:
        validate_primer_sequence(primer, name="--primers")

    logger.info(f"Analyzing 3' stability for {len(args.primers)} primers")

    try:
        analyzer = ThreePrimeStabilityAnalyzer()
        results = []
        for primer in args.primers:
            stability = analyzer.analyze_primer(primer)
            results.append(stability)

        # Ensure output directory exists
        output_dir = os.path.dirname(args.output)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # Write report
        with open(args.output, 'w') as f:
            for primer, stability in zip(args.primers, results):
                f.write(f"{primer}\t{stability}\n")

        logger.info(f"Stability analysis complete! Report saved to: {args.output}")

    except Exception as e:
        logger.error(f"Stability analysis failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_ml_predict(args):
    """Deep learning-based amplification prediction"""
    from neoswga.core import deep_learning

    logger.info(f"Running ML prediction on primers from: {args.primers_file}")

    try:
        # Load primers from file
        with open(args.primers_file, 'r') as f:
            primer_list = [line.strip() for line in f if line.strip()]

        logger.info(f"Loaded {len(primer_list)} primers")

        # Run prediction
        results = deep_learning.predict_amplification(
            primers=primer_list,
            model_path=args.model,
            use_gpu=args.use_gpu,
            output_dir=args.output
        )

        logger.info(f"ML prediction complete! Results saved to: {args.output}")

    except Exception as e:
        logger.error(f"ML prediction failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_design_oligos(args):
    """Alternative comprehensive primer design system"""
    from neoswga.core import optimal_oligo_generator

    logger.info(f"Designing {args.num_primers} primers for: {args.genome}")
    logger.info(f"K-mer range: {args.min_k}-{args.max_k}")

    try:
        generator = optimal_oligo_generator.OligoGenerator(
            genome_path=args.genome,
            k_min=args.min_k,
            k_max=args.max_k
        )

        primers = generator.design_primers(
            num_primers=args.num_primers,
            output_dir=args.output
        )

        logger.info(f"Primer design complete! {len(primers)} primers saved to: {args.output}")

    except Exception as e:
        logger.error(f"Primer design failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


# =========================================================================
# CATEGORY 3: Handler functions for pipeline orphaned features
# =========================================================================

def run_auto_pipeline(args):
    """Automatic parameter optimization pipeline"""
    from neoswga.core import auto_swga_pipeline

    logger.info(f"Running auto-tuning pipeline with {args.iterations} iterations")
    logger.info(f"Base parameters: {args.json_file}")

    try:
        optimized_params = auto_swga_pipeline.auto_optimize(
            base_params_file=args.json_file,
            iterations=args.iterations,
            output_dir=args.output
        )

        logger.info(f"Auto-tuning complete! Optimized parameters saved to: {args.output}")

    except Exception as e:
        logger.error(f"Auto-tuning failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_multi_genome(args):
    """Pan-genome primer design for multiple targets"""
    from neoswga.core.multi_genome_pipeline import MultiGenomePipeline
    from neoswga.core.multi_genome_filter import GenomeSet
    from pathlib import Path

    logger.info(f"Designing pan-genome primers for {len(args.genomes)} genomes")

    try:
        # Create genome set from input genomes
        genome_set = GenomeSet()

        # Add all genomes as targets
        for i, genome_path in enumerate(args.genomes):
            genome_name = Path(genome_path).stem
            genome_set.add_genome(
                name=genome_name,
                fasta_path=genome_path,
                role="target"
            )
            logger.info(f"  Added target: {genome_name}")

        # Add background genomes if provided
        if hasattr(args, 'background') and args.background:
            for bg_path in args.background:
                bg_name = Path(bg_path).stem
                genome_set.add_genome(
                    name=bg_name,
                    fasta_path=bg_path,
                    role="background"
                )
                logger.info(f"  Added background: {bg_name}")

        # Add blacklist genomes if provided
        if hasattr(args, 'blacklist') and args.blacklist:
            for bl_path in args.blacklist:
                bl_name = Path(bl_path).stem
                genome_set.add_genome(
                    name=bl_name,
                    fasta_path=bl_path,
                    role="blacklist",
                    penalty_weight=5.0
                )
                logger.info(f"  Added blacklist: {bl_name}")

        # Get optional parameters
        primer_count = getattr(args, 'num_primers', 12)
        min_k = getattr(args, 'min_k', None)
        max_k = getattr(args, 'max_k', None)
        kmer_range = (min_k, max_k) if min_k and max_k else None
        polymerase = getattr(args, 'polymerase', None)

        # Create and run pipeline
        pipeline = MultiGenomePipeline(
            genome_set=genome_set,
            output_dir=args.output,
            kmer_range=kmer_range,
            preferred_polymerase=polymerase,
            primer_count=primer_count,
            validate_with_simulation=getattr(args, 'validate_simulation', False)
        )

        result = pipeline.run(verbose=not getattr(args, 'quiet', False))
        pipeline.save_results(result)

        logger.info(f"Pan-genome design complete! Results saved to: {args.output}")
        logger.info(f"Designed {result.primer_count} primers")
        logger.info(f"Mean enrichment: {result.mean_enrichment:.1f}x")

    except Exception as e:
        logger.error(f"Multi-genome pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


# =========================================================================
# CATEGORY 4: Handler function for simulation orphaned features
# =========================================================================

def run_simulate(args):
    """Monte Carlo amplification simulation"""
    from neoswga.core.replication_simulator import Phi29Simulator, SimulationConfig, simulate_primer_set
    from neoswga.core import reaction_conditions as rc
    from neoswga.core.genome_io import GenomeLoader
    from Bio import SeqIO
    from pathlib import Path
    import json

    # Resolve primers from --primers or --from-results
    if args.from_results:
        import pandas as pd
        source = args.from_results
        if os.path.isdir(source):
            # Look for step4 CSV in the directory
            candidates = [
                os.path.join(source, 'step4_improved_df.csv'),
                os.path.join(source, 'step4_df.csv'),
            ]
            csv_path = None
            for c in candidates:
                if os.path.exists(c):
                    csv_path = c
                    break
            if csv_path is None:
                logger.error(f"No step4 CSV found in {source}")
                logger.error("Expected: step4_improved_df.csv or step4_df.csv")
                sys.exit(1)
        else:
            csv_path = source
            if not os.path.exists(csv_path):
                logger.error(f"File not found: {csv_path}")
                sys.exit(1)

        logger.info(f"Loading primers from {csv_path}")
        df = pd.read_csv(csv_path)
        if 'primer' not in df.columns:
            logger.error(f"CSV file {csv_path} has no 'primer' column")
            sys.exit(1)
        primers = df['primer'].tolist()
        logger.info(f"Loaded {len(primers)} primers from results")
    else:
        primers = args.primers

    logger.info(f"Running amplification simulation for {len(primers)} primers")
    logger.info(f"Genome: {args.genome}")

    try:
        os.makedirs(args.output, exist_ok=True)

        # Load genome sequence
        logger.info("Loading genome...")
        loader = GenomeLoader()
        genome_sequence = loader.load_genome(args.genome)
        genome_length = len(genome_sequence)
        logger.info(f"  Genome length: {genome_length:,} bp")

        # Find primer positions in genome
        logger.info("Finding primer binding positions...")
        primer_positions = {}
        for primer in primers:
            primer_positions[primer] = {'forward': [], 'reverse': []}
            primer_rc = primer.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

            # Find forward positions
            pos = 0
            while True:
                pos = genome_sequence.find(primer, pos)
                if pos == -1:
                    break
                primer_positions[primer]['forward'].append(pos)
                pos += 1

            # Find reverse positions
            pos = 0
            while True:
                pos = genome_sequence.find(primer_rc, pos)
                if pos == -1:
                    break
                primer_positions[primer]['reverse'].append(pos)
                pos += 1

            total = len(primer_positions[primer]['forward']) + len(primer_positions[primer]['reverse'])
            logger.info(f"  {primer}: {total} binding sites")

        # Get reaction conditions based on polymerase type
        polymerase = getattr(args, 'polymerase', 'phi29')

        # Get polymerase characteristics for optimal temperature
        from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
        if polymerase not in POLYMERASE_CHARACTERISTICS:
            logger.error(f"Unknown polymerase: {polymerase}")
            logger.info(f"Available: {list(POLYMERASE_CHARACTERISTICS.keys())}")
            sys.exit(1)

        poly_info = POLYMERASE_CHARACTERISTICS[polymerase]
        optimal_temp = poly_info['optimal_temp']
        processivity = poly_info['processivity']

        logger.info(f"Polymerase: {polymerase}")
        logger.info(f"  Optimal temperature: {optimal_temp}C")
        logger.info(f"  Processivity: {processivity:,} bp")

        # Create conditions with appropriate temperature
        conditions = rc.ReactionConditions(
            temp=optimal_temp,
            polymerase=polymerase
        )

        # Configure simulation
        duration = args.duration * 60.0  # Convert minutes to seconds
        config = SimulationConfig(
            duration=duration,
            polymerase_type=polymerase
        )

        # Run simulation with replicates
        n_replicates = getattr(args, 'replicates', 5)
        logger.info(f"\nRunning {n_replicates} simulation replicates...")

        results = simulate_primer_set(
            primers=primers,
            primer_positions=primer_positions,
            genome_length=genome_length,
            genome_sequence=genome_sequence,
            conditions=conditions,
            n_replicates=n_replicates
        )

        logger.info(f"\nSimulation complete!")
        logger.info(f"  Mean coverage: {results['mean_coverage']:.1%}")
        logger.info(f"  Std coverage: {results['std_coverage']:.1%}")
        logger.info(f"  Mean forks created: {results['mean_forks_created']:.0f}")
        logger.info(f"  Mean fork travel: {results['mean_fork_travel']:,.0f} bp")

        # Save results
        results_file = Path(args.output) / "simulation_results.json"
        with open(results_file, 'w') as f:
            json.dump({
                'primers': primers,
                'genome_length': genome_length,
                'mean_coverage': results['mean_coverage'],
                'std_coverage': results['std_coverage'],
                'mean_forks_created': results['mean_forks_created'],
                'mean_fork_travel': results['mean_fork_travel'],
                'n_replicates': n_replicates,
                'polymerase': polymerase,
                'duration_seconds': duration
            }, f, indent=2)
        logger.info(f"Results saved to: {results_file}")

        # Generate visualizations if requested
        if args.visualize:
            logger.info("Generating visualization plots...")
            try:
                from neoswga.core import simulation_plots
                simulation_plots.generate_plots(results, args.output)
            except ImportError:
                logger.warning("Visualization module not available")

        # Generate HTML report if requested
        if args.report:
            logger.info("Generating HTML report...")
            try:
                from neoswga.core import simulation_report
                simulation_report.generate_report(results, args.output)
            except ImportError:
                logger.warning("Report module not available")

    except Exception as e:
        logger.error(f"Simulation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_active_learn(args):
    """
    Active learning for iterative primer optimization.

    Uses Bayesian optimization to select the most informative primer sets
    to test experimentally, learning from feedback to improve future designs.
    """
    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    import pandas as pd
    import json as json_module

    quiet = getattr(args, 'quiet', False)

    if not quiet:
        logger.info("Active Learning for Primer Optimization")
        logger.info("=" * 60)

    # Set json_file for pipeline initialization
    merge_args_to_parameter(args, parameter, ['json_file'])

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    try:
        # Initialize pipeline to get genome info
        if not quiet:
            logger.info("Initializing pipeline...")
        core_pipeline._initialize()

        fg_prefixes = core_pipeline.fg_prefixes
        bg_prefixes = core_pipeline.bg_prefixes
        fg_seq_lengths = core_pipeline.fg_seq_lengths
        bg_seq_lengths = core_pipeline.bg_seq_lengths

        if not quiet:
            logger.info(f"  Target genomes: {len(fg_prefixes)}")
            logger.info(f"  Background genomes: {len(bg_prefixes)}")

        # Load candidate primers from scored step3 output
        data_dir = getattr(parameter, 'data_dir', '.')
        step3_file = os.path.join(data_dir, 'step3_df.csv')

        if not os.path.exists(step3_file):
            logger.error(f"Step 3 output not found: {step3_file}")
            logger.error("Run 'neoswga score' first to generate candidate primers.")
            sys.exit(1)

        if not quiet:
            logger.info(f"Loading candidates from {step3_file}...")
        df = pd.read_csv(step3_file)

        primer_col = 'primer' if 'primer' in df.columns else 'seq'
        if primer_col not in df.columns:
            logger.error(f"Invalid step3 file: missing 'primer' or 'seq' column. Found: {list(df.columns)}")
            sys.exit(1)

        candidates = df[primer_col].tolist()
        if not quiet:
            logger.info(f"  Loaded {len(candidates)} candidate primers")

        # Initialize position cache
        from neoswga.core.position_cache import PositionCache
        if not quiet:
            logger.info("Loading primer binding positions...")
        cache = PositionCache(fg_prefixes, candidates[:500])  # Limit for memory

        # Import active learning module
        try:
            from neoswga.core.active_learning import ActiveLearningLoop, HAS_GP
            if not HAS_GP:
                logger.error("Active learning requires scikit-learn with Gaussian Processes.")
                logger.error("Install with: pip install scikit-learn")
                sys.exit(1)
        except ImportError as e:
            logger.error(f"Failed to import active learning module: {e}")
            sys.exit(1)

        # Initialize active learning loop
        if not quiet:
            logger.info("Initializing active learning loop...")
        al_loop = ActiveLearningLoop(
            cache=cache,
            fg_prefixes=fg_prefixes,
            bg_prefixes=bg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_seq_lengths=bg_seq_lengths,
            results_dir=args.output
        )

        # Load experimental results if provided
        if args.experimental_results:
            if not quiet:
                logger.info(f"Loading experimental results from {args.experimental_results}...")
            try:
                from neoswga.core.active_learning import ExperimentalResult
                exp_df = pd.read_csv(args.experimental_results)
                for _, row in exp_df.iterrows():
                    result = ExperimentalResult(
                        primer_set=row['primer_set'].split(';') if isinstance(row['primer_set'], str) else row['primer_set'],
                        timestamp=row.get('timestamp', ''),
                        enrichment_fold=row['enrichment_fold'],
                        uniformity_score=row.get('uniformity_score', 0.5),
                        temperature=row.get('temperature', 30.0),
                        time_hours=row.get('time_hours', 4.0)
                    )
                    al_loop.add_experimental_result(result)
                if not quiet:
                    logger.info(f"  Loaded {len(exp_df)} experimental results")
            except Exception as e:
                logger.warning(f"Could not load experimental results: {e}")

        # Recommend next experiment
        if not quiet:
            logger.info("")
            logger.info("Generating candidate primer sets...")
        recommendation = al_loop.recommend_next_experiment(
            candidates=candidates[:500],  # Limit for performance
            n_candidates=args.num_candidates,
            max_primers=args.max_primers,
            exploration_weight=args.exploration_weight
        )

        # Output results
        output_file = os.path.join(args.output, 'recommendation.json')
        with open(output_file, 'w') as f:
            output_data = {
                'primer_set': recommendation['primer_set'],
                'predicted_enrichment': recommendation['predicted_enrichment'],
                'prediction_uncertainty': recommendation['prediction_uncertainty'],
                'n_experiments_so_far': recommendation['n_experiments_so_far'],
                'features': {
                    'n_primers': recommendation['features'].n_primers,
                    'coverage_fraction': recommendation['features'].coverage_fraction,
                    'mean_gap': recommendation['features'].mean_gap,
                    'mean_tm': recommendation['features'].mean_tm
                }
            }
            json_module.dump(output_data, f, indent=2)

        # Print summary
        if not quiet:
            logger.info("")
            logger.info("=" * 60)
            logger.info("RECOMMENDATION")
            logger.info("=" * 60)
            logger.info(f"Primer set ({len(recommendation['primer_set'])} primers):")
            for primer in recommendation['primer_set']:
                logger.info(f"  {primer}")
            logger.info("")
            logger.info(f"Predicted enrichment: {recommendation['predicted_enrichment']:.2f}x")
            logger.info(f"Prediction uncertainty: {recommendation['prediction_uncertainty']:.2f}")
            logger.info(f"Previous experiments: {recommendation['n_experiments_so_far']}")
            logger.info("")
            logger.info(f"Results saved to: {output_file}")
            logger.info("")
            logger.info("Next steps:")
            logger.info("  1. Test this primer set experimentally")
            logger.info("  2. Record enrichment and uniformity results")
            logger.info("  3. Add results to experimental_results.csv")
            logger.info("  4. Re-run active-learn to get next recommendation")

    except Exception as e:
        logger.error(f"Active learning failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_expand_primers(args):
    """
    Expand existing primer set with additional primers.

    Useful for iterative design where some primers have been validated
    experimentally and user wants to add more to fill coverage gaps.
    """
    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.primer_expansion import PrimerExpander
    from neoswga.core.position_cache import PositionCache
    import pandas as pd
    import json as json_module

    quiet = getattr(args, 'quiet', False)

    if not quiet:
        logger.info("Primer Set Expansion")
        logger.info("=" * 60)

    # Load and validate primers using consolidated helper
    fixed_primers = collect_primers_from_args(
        cli_primers=args.fixed_primers,
        primers_file=args.fixed_primers_file,
        name="fixed primer",
        allow_empty=True,  # Fixed primers can be empty for initial design
    )

    failed_primers = collect_primers_from_args(
        cli_primers=args.failed_primers,
        primers_file=args.failed_primers_file,
        name="failed primer",
        allow_empty=True,  # Failed primers are optional
    )

    if not quiet:
        logger.info(f"Fixed primers: {len(fixed_primers)}")
        logger.info(f"Failed primers (excluded): {len(failed_primers)}")
        logger.info(f"Target new primers: {args.num_new}")

    # Set json_file for pipeline initialization
    merge_args_to_parameter(args, parameter, ['json_file'])

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    try:
        # Initialize pipeline to get genome info
        if not quiet:
            logger.info("Initializing pipeline...")
        core_pipeline._initialize()

        fg_prefixes = core_pipeline.fg_prefixes
        bg_prefixes = core_pipeline.bg_prefixes
        fg_seq_lengths = core_pipeline.fg_seq_lengths
        bg_seq_lengths = core_pipeline.bg_seq_lengths

        # Load candidate primers from scored step3 output
        data_dir = getattr(parameter, 'data_dir', '.')
        step3_file = os.path.join(data_dir, 'step3_df.csv')

        if not os.path.exists(step3_file):
            logger.error(f"Step 3 output not found: {step3_file}")
            logger.error("Run 'neoswga score' first to generate candidate primers.")
            sys.exit(1)

        if not quiet:
            logger.info(f"Loading candidates from {step3_file}...")

        df = pd.read_csv(step3_file)
        primer_col = 'primer' if 'primer' in df.columns else 'seq'
        if primer_col not in df.columns:
            logger.error(f"Invalid step3 file: missing 'primer' or 'seq' column. Found: {list(df.columns)}")
            sys.exit(1)
        candidates = df[primer_col].tolist()

        if not quiet:
            logger.info(f"Candidate pool: {len(candidates)} primers")

        # Initialize position cache with candidates and fixed primers
        all_primers = list(set(candidates + fixed_primers))
        if not quiet:
            logger.info("Loading position data...")

        cache = PositionCache(fg_prefixes, all_primers)

        # Create expander
        expander = PrimerExpander(
            position_cache=cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
        )

        # Run expansion
        result = expander.expand(
            candidates=candidates,
            fixed_primers=fixed_primers,
            failed_primers=failed_primers,
            target_new=args.num_new,
            optimization_method=args.optimization_method,
            verbose=not quiet,
        )

        # Save results
        result_file = os.path.join(args.output, 'expansion_result.json')
        with open(result_file, 'w') as f:
            json_module.dump(result.to_dict(), f, indent=2)

        # Save primers as CSV
        primers_file = os.path.join(args.output, 'expanded_primers.csv')
        pd.DataFrame({
            'primer': result.combined_set,
            'type': ['fixed'] * result.n_fixed + ['new'] * result.n_new,
        }).to_csv(primers_file, index=False)

        # Print summary
        if not quiet:
            logger.info("")
            logger.info("=" * 60)
            logger.info("EXPANSION COMPLETE")
            logger.info("=" * 60)
            logger.info(f"Fixed primers: {result.n_fixed}")
            logger.info(f"New primers: {result.n_new}")
            logger.info(f"Total set: {result.n_total}")
            logger.info("")
            logger.info(f"Coverage: {result.coverage_before:.1%} -> {result.coverage_after:.1%}")
            logger.info(f"Improvement: {result.predicted_improvement:.2f}x")
            logger.info(f"Gaps remaining: {result.gaps_remaining}")
            logger.info("")
            logger.info("New primers:")
            for primer in result.new_primers:
                logger.info(f"  {primer}")
            logger.info("")
            logger.info(f"Results saved to: {args.output}")

    except Exception as e:
        logger.error(f"Primer expansion failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def run_predict_efficiency(args):
    """
    Predict efficiency of primer set before synthesis.

    Provides unified confidence score with SYNTHESIZE/CAUTION/DO_NOT_SYNTHESIZE
    recommendation.
    """
    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.efficiency_predictor import EfficiencyPredictor
    from neoswga.core.experimental_tracker import ExperimentalTracker
    from neoswga.core.position_cache import PositionCache
    import json as json_module

    quiet = getattr(args, 'quiet', False)

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
    merge_args_to_parameter(args, parameter, ['json_file'])

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
            fg_genomes = getattr(parameter, 'fg_genomes', [])
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
            with open(args.output, 'w') as f:
                json_module.dump(result.to_dict(), f, indent=2)
            if not quiet:
                logger.info(f"\nResults saved to: {args.output}")

        # Print summary
        if not quiet:
            logger.info("")
            logger.info("=" * 60)
            logger.info("PREDICTION RESULT")
            logger.info("=" * 60)
            logger.info(f"Confidence: {result.confidence_score:.2f} ({result.recommendation.value})")
            logger.info(f"")
            logger.info(f"Predicted enrichment: {result.predicted_enrichment:.0f}x")
            logger.info(f"  95% CI: {result.enrichment_ci_low:.0f}x - {result.enrichment_ci_high:.0f}x")
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


def run_background_list(args):
    """List available pre-computed background genomes."""
    from neoswga.core.background_registry import BackgroundRegistry

    quiet = getattr(args, 'quiet', False)

    registry = BackgroundRegistry(auto_discover=False)

    # Run discovery if requested
    if args.discover:
        if not quiet:
            logger.info("Discovering backgrounds...")
        n_discovered = registry.discover(verbose=not quiet)
        if not quiet:
            logger.info(f"Discovered {n_discovered} new backgrounds")
            logger.info("")

    # Search or list all
    if args.search:
        entries = registry.search(args.search)
        if not quiet:
            logger.info(f"Search results for '{args.search}':")
    else:
        entries = registry.list_all()
        if not quiet:
            logger.info("Available backgrounds:")

    if not entries:
        logger.info("  No backgrounds found")
        logger.info("")
        logger.info("Add backgrounds with:")
        logger.info("  neoswga background-add --name 'Human GRCh38' --species 'Homo sapiens' --bloom-path /path/to/bloom.pkl")
        logger.info("")
        logger.info("Or run discovery to find existing files:")
        logger.info("  neoswga background-list --discover")
        return

    logger.info("")
    for entry in entries:
        logger.info(f"  {entry.name}")
        if not quiet:
            logger.info(f"    Species: {entry.species}")
            if entry.genome_size > 0:
                size_gb = entry.genome_size / 1e9
                logger.info(f"    Size: {size_gb:.2f} Gbp")
            if entry.has_bloom:
                logger.info(f"    Bloom filter: {entry.bloom_path}")
            if entry.has_kmers:
                logger.info(f"    K-mer files: {entry.kmer_prefix}_*mer_all.txt")
                logger.info(f"    K range: {entry.k_range[0]}-{entry.k_range[1]}")
            logger.info("")

    if not quiet:
        logger.info(f"Total: {len(entries)} backgrounds")


def run_background_add(args):
    """Add a pre-computed background to the registry."""
    from neoswga.core.background_registry import BackgroundRegistry

    registry = BackgroundRegistry(auto_discover=False)

    # Validate paths
    if args.bloom_path and not os.path.exists(args.bloom_path):
        logger.error(f"Bloom filter not found: {args.bloom_path}")
        sys.exit(1)

    if args.kmer_prefix:
        # Check if at least one k-mer file exists
        found = False
        for k in range(args.min_k, args.max_k + 1):
            kmer_file = f"{args.kmer_prefix}_{k}mer_all.txt"
            if os.path.exists(kmer_file):
                found = True
                break
        if not found:
            logger.warning(f"No k-mer files found at {args.kmer_prefix}_*mer_all.txt")

    # Add to registry
    entry = registry.add(
        name=args.name,
        species=args.species,
        genome_size=args.genome_size,
        bloom_path=args.bloom_path,
        kmer_prefix=args.kmer_prefix,
        k_range=(args.min_k, args.max_k),
        description=args.description or "",
        overwrite=args.overwrite,
    )

    logger.info(f"Added background: {entry.name}")
    logger.info("")
    logger.info("Background genome registered. It will be used in subsequent filter runs.")


# =========================================================================
# UTILITY: Genome library management
# =========================================================================

def run_genome_add(args):
    """Add a genome to the pre-calculated library."""
    from neoswga.core.genome_library import GenomeLibrary

    # Parse k-ranges
    k_ranges = []
    for part in args.k_ranges.split(','):
        parts = part.strip().split('-')
        if len(parts) == 2:
            k_ranges.append((int(parts[0]), int(parts[1])))
        else:
            logger.error(f"Invalid k-range format: {part}. Use min-max (e.g., 6-12)")
            sys.exit(1)

    fasta_path = args.fasta
    if not os.path.exists(fasta_path):
        logger.error(f"FASTA file not found: {fasta_path}")
        sys.exit(1)

    library = GenomeLibrary()
    bloom = 'no' if args.no_bloom else 'auto'
    entry = library.add(
        name=args.name,
        fasta_path=fasta_path,
        role=args.role,
        species=args.species,
        k_ranges=k_ranges,
        build_bloom=bloom,
    )
    print(f"Added '{entry.name}' to genome library")
    print(f"  Size: {entry.genome_size:,} bp")
    print(f"  GC: {entry.gc_content:.1%}")
    print(f"  K-mer ranges: {entry.computed_k_ranges}")
    if entry.bloom_path:
        print(f"  Bloom filter: {entry.bloom_path}")


def run_genome_list(args):
    """List genomes in the pre-calculated library."""
    from neoswga.core.genome_library import GenomeLibrary

    library = GenomeLibrary()
    entries = library.list()

    if not entries:
        print("No genomes in library.")
        print("Add one with: neoswga genome-add <name> <fasta>")
        return

    print(f"Genome library ({len(entries)} entries):")
    print("-" * 60)
    for entry in entries:
        size_str = f"{entry.genome_size / 1e6:.1f} Mbp" if entry.genome_size >= 1e6 else f"{entry.genome_size:,} bp"
        print(f"  {entry.name:<25} {size_str:>12}  GC={entry.gc_content:.1%}  role={entry.role}")
        if getattr(args, 'verbose', False):
            print(f"    FASTA: {entry.fasta_path}")
            print(f"    K-mer ranges: {entry.computed_k_ranges}")
            if entry.bloom_path:
                print(f"    Bloom: {entry.bloom_path}")
            print(f"    Created: {entry.created_date}")


def run_genome_remove(args):
    """Remove a genome from the pre-calculated library."""
    from neoswga.core.genome_library import GenomeLibrary

    library = GenomeLibrary()
    if library.remove(args.name):
        print(f"Removed '{args.name}' from genome library")
    else:
        logger.error(f"Genome '{args.name}' not found in library")
        sys.exit(1)


# =========================================================================
# UNIFIED: Complete pipeline handler
# =========================================================================

def run_design(args):
    """Run complete primer design pipeline (all 4 steps)"""
    check_jellyfish_available()
    logger.info("Running complete primer design pipeline")

    # Check for auto-tuning mode
    if hasattr(args, 'auto') and args.auto:
        logger.info(f"Auto-tuning mode enabled ({args.auto_iterations} iterations)")
        from neoswga.core import auto_swga_pipeline
        try:
            optimized_params = auto_swga_pipeline.auto_optimize(
                base_params_file=args.json_file,
                iterations=args.auto_iterations,
                output_dir=getattr(args, 'output', './auto_optimized')
            )
            logger.info("Auto-tuning complete! Using optimized parameters for pipeline run.")
        except Exception as e:
            logger.error(f"Auto-tuning failed: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # Check for multi-genome mode
    if hasattr(args, 'multi_genome') and args.multi_genome:
        logger.info(f"Multi-genome mode enabled for {len(args.multi_genome)} genomes")
        from neoswga.core import multi_genome_pipeline
        try:
            results = multi_genome_pipeline.design_pan_genome_primers(
                genome_paths=args.multi_genome,
                min_coverage=args.min_coverage,
                output_dir=getattr(args, 'output', './multi_genome_output')
            )
            logger.info(f"Pan-genome design complete! Designed {len(results['primers'])} primers")
            return
        except Exception as e:
            logger.error(f"Multi-genome pipeline failed: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

    # Determine which steps to run
    start_step = getattr(args, 'start_from', None) or 1
    stop_step = getattr(args, 'stop_at', None) or 4

    if start_step > stop_step:
        logger.error("--start-from must be <= --stop-at")
        sys.exit(1)

    logger.info(f"Running steps {start_step} through {stop_step}")

    # Ensure optimize-specific attributes exist on args namespace
    # (the design subparser doesn't define these but run_step4 needs them)
    optimize_defaults = {
        'optimization_method': 'hybrid', 'num_primers': None,
        'iterations': None, 'max_sets': None, 'no_background': False,
        'use_mechanistic_model': False, 'mechanistic_weight': 0.3,
        'auto_size': False, 'application': 'enrichment',
        'validate_with_simulation': False, 'method_guide': False,
        'no_bg_prefilter': False, 'strategy': None,
        'use_position_cache': True, 'use_background_filter': False,
        'use_cooperative_binding': False, 'primer_strategy': None,
        'enable_qa': False,
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
        print("\nTip: New to neoswga? Start with 'neoswga init --genome target.fasta'", file=sys.stderr)
        sys.exit(1)

    # Set verbosity
    if getattr(args, 'quiet', False):
        logging.getLogger().setLevel(logging.WARNING)
    elif getattr(args, 'verbose', False):
        logging.getLogger().setLevel(logging.DEBUG)

    # Dispatch to appropriate command
    commands = {
        # Standard pipeline
        'count-kmers': run_step1,
        'filter': run_step2,
        'score': run_step3,
        'optimize': run_step4,

        # Unified pipeline
        'design': run_design,

        # Utility commands
        'build-filter': run_build_filter,
        'validate': run_validate,
        'show-presets': lambda args: show_presets(),
        'init': run_init,
        'validate-params': run_validate_params,
        'schema': run_schema,
        'validate-model': run_validate_model,
        'interpret': run_interpret,
        'report': run_report,
        'export': run_export,
        'start': run_start,
        'suggest': run_suggest,
        'optimize-conditions': optimize_conditions,
        'analyze-set': analyze_primer_set,

        # Category 1: Orphaned analysis features (now exposed!)
        'analyze-genome': run_analyze_genome,
        'analyze-dimers': run_analyze_dimers,
        'analyze-stability': run_analyze_stability,
        'ml-predict': run_ml_predict,
        'design-oligos': run_design_oligos,

        # Category 3: Orphaned pipeline features (now exposed!)
        'auto-pipeline': run_auto_pipeline,
        'multi-genome': run_multi_genome,

        # Category 4: Orphaned simulation features (now exposed!)
        'simulate': run_simulate,

        # Category 5: Active learning (experimental)
        'active-learn': run_active_learn,

        # Category 6: Iterative design
        'expand-primers': run_expand_primers,
        'predict-efficiency': run_predict_efficiency,

        # Category 7: Background registry
        'background-list': run_background_list,
        'background-add': run_background_add,

        # Category 8: Genome library
        'genome-add': run_genome_add,
        'genome-list': run_genome_list,
        'genome-remove': run_genome_remove,
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
            logger.error("Run 'neoswga init --genome target.fasta' to create a complete params.json")
            logger.error("Or run 'neoswga validate-params -j params.json' to check your configuration")
            sys.exit(1)
        except (ValueError, TypeError) as e:
            logger.error(f"Invalid input: {e}")
            if getattr(args, 'verbose', False):
                import traceback
                traceback.print_exc()
            sys.exit(1)
        except PermissionError as e:
            logger.error(f"Permission denied: {e}")
            logger.error("Check that the output directory is writable.")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Command failed: {e}")
            if getattr(args, 'verbose', False):
                import traceback
                traceback.print_exc()
            else:
                logger.error("Run with --verbose for full traceback")
            sys.exit(1)
    else:
        logger.error(f"Unknown command: {args.command}")
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
