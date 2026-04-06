import os
import json
import logging
from dataclasses import dataclass, field
from typing import List, Optional
from neoswga.core import utility as _utility

logger = logging.getLogger(__name__)

src_dir=os.path.dirname(os.path.abspath(__file__))

# Optimization parameters (loaded from JSON)
num_primers = 6
target_set_size = 6
_json_data = {}  # Store raw JSON data for access by CLI


# =============================================================================
# Typed Configuration Class
# =============================================================================
# This dataclass provides type-safe access to pipeline parameters.
# New code should use get_current_config() instead of accessing globals directly.

@dataclass
class PipelineParameters:
    """
    Typed configuration for the SWGA pipeline.

    This class provides type-safe access to all pipeline parameters.
    Use get_current_config() to get the current parameter state.

    Example:
        config = get_current_config()
        if config.min_fg_freq > 1e-4:
            print("Frequency threshold is high")
    """
    # K-mer range (supports 6-18bp, 15-18bp requires long_primer_mode optimizations)
    min_k: int = 6
    max_k: int = 12

    # Long primer mode: enables optimizations for 15-18bp primers
    # Auto-enabled when max_k >= 15, includes:
    # - K-mer sampling (default 2% for 18bp)
    # - GPU acceleration (if available)
    # - 2-bit k-mer encoding (memory reduction)
    # - Hybrid prefix filtering (optional)
    long_primer_mode: bool = False

    # Filtering thresholds
    min_fg_freq: float = 1e-5
    max_bg_freq: float = 5e-6
    min_tm: float = 15.0
    max_tm: float = 45.0
    max_gini: float = 0.6
    max_primer: int = 500
    # Unitless amplification prediction score (scale ~0-20). Combines Tm
    # optimality, GC content, 3' stability, and binding energy features.
    # Higher values indicate stronger predicted amplification; the default
    # of 10.0 retains above-average primers.
    min_amp_pred: float = 10.0

    # Dimer constraints
    max_dimer_bp: int = 3
    max_self_dimer_bp: int = 4

    # GC content filtering
    gc_min: float = 0.375
    gc_max: float = 0.625
    gc_tolerance: float = 0.15
    genome_gc: Optional[float] = None

    # Polymerase and reaction conditions
    polymerase: str = 'phi29'
    reaction_temp: Optional[float] = None
    na_conc: float = 50.0
    mg_conc: float = 2.0
    primer_conc: float = 0.5e-6

    # Thermodynamic additives
    dmso_percent: float = 0.0
    betaine_m: float = 0.0
    trehalose_m: float = 0.0
    formamide_percent: float = 0.0
    ethanol_percent: float = 0.0
    urea_m: float = 0.0
    tmac_m: float = 0.0

    # Optimization settings
    iterations: int = 8
    max_sets: int = 5
    retries: int = 5
    top_set_count: int = 10
    drop_iterations: int = 5
    selection_metric: str = 'deterministic'

    # Genome topology
    fg_circular: bool = True
    bg_circular: bool = False

    # K-mer sampling (for performance)
    sample_rate: Optional[float] = None
    min_sample_count: int = 5

    # Bloom filter (for large backgrounds)
    use_bloom_filter: bool = False
    bloom_filter_path: Optional[str] = None

    # Paths
    data_dir: str = ''
    src_dir: str = ''

    # Genome data
    fg_genomes: List[str] = field(default_factory=list)
    bg_genomes: List[str] = field(default_factory=list)
    fg_prefixes: List[str] = field(default_factory=list)
    bg_prefixes: List[str] = field(default_factory=list)
    fg_seq_lengths: List[int] = field(default_factory=list)
    bg_seq_lengths: List[int] = field(default_factory=list)

    # Exclusion genomes (zero-tolerance filtering)
    excl_genomes: List[str] = field(default_factory=list)
    excl_prefixes: List[str] = field(default_factory=list)
    excl_threshold: int = 0  # Max allowed hits in exclusion genome (0 = any hit rejects)

    # Blacklist genomes (penalty-weighted filtering)
    bl_genomes: List[str] = field(default_factory=list)
    bl_prefixes: List[str] = field(default_factory=list)
    bl_seq_lengths: List[int] = field(default_factory=list)
    bl_penalty: float = 5.0
    max_bl_freq: float = 0.0

    # Runtime
    cpus: int = 1
    verbose: bool = False


def get_current_config() -> PipelineParameters:
    """
    Get current pipeline parameters as a typed object.

    Returns a PipelineParameters instance with all current global values.
    This is the recommended way to access parameters in new code.

    Returns:
        PipelineParameters with current global state
    """
    return PipelineParameters(
        min_k=globals().get('min_k', 6),
        max_k=globals().get('max_k', 12),
        long_primer_mode=globals().get('long_primer_mode', False),
        min_fg_freq=globals().get('min_fg_freq', 1e-5),
        max_bg_freq=globals().get('max_bg_freq', 5e-6),
        min_tm=globals().get('min_tm', 15.0),
        max_tm=globals().get('max_tm', 45.0),
        max_gini=globals().get('max_gini', 0.6),
        max_primer=globals().get('max_primer', 500),
        min_amp_pred=globals().get('min_amp_pred', 10.0),
        max_dimer_bp=globals().get('max_dimer_bp', 3),
        max_self_dimer_bp=globals().get('max_self_dimer_bp', 4),
        gc_min=globals().get('gc_min', 0.375),
        gc_max=globals().get('gc_max', 0.625),
        gc_tolerance=globals().get('gc_tolerance', 0.15),
        genome_gc=globals().get('genome_gc', None),
        polymerase=globals().get('polymerase', 'phi29'),
        reaction_temp=globals().get('reaction_temp', None),
        na_conc=globals().get('na_conc', 50.0),
        mg_conc=globals().get('mg_conc', 2.0),
        primer_conc=globals().get('primer_conc', 0.5e-6),
        dmso_percent=globals().get('dmso_percent', 0.0),
        betaine_m=globals().get('betaine_m', 0.0),
        trehalose_m=globals().get('trehalose_m', 0.0),
        formamide_percent=globals().get('formamide_percent', 0.0),
        ethanol_percent=globals().get('ethanol_percent', 0.0),
        urea_m=globals().get('urea_m', 0.0),
        tmac_m=globals().get('tmac_m', 0.0),
        iterations=globals().get('iterations', 8),
        max_sets=globals().get('max_sets', 5),
        retries=globals().get('retries', 5),
        top_set_count=globals().get('top_set_count', 10),
        drop_iterations=globals().get('drop_iterations', 5),
        selection_metric=globals().get('selection_metric', 'deterministic'),
        fg_circular=globals().get('fg_circular', True),
        bg_circular=globals().get('bg_circular', False),
        sample_rate=globals().get('sample_rate', None),
        min_sample_count=globals().get('min_sample_count', 5),
        use_bloom_filter=globals().get('use_bloom_filter', False),
        bloom_filter_path=globals().get('bloom_filter_path', None),
        data_dir=globals().get('data_dir', ''),
        src_dir=globals().get('src_dir', ''),
        fg_genomes=globals().get('fg_genomes', []) or [],
        bg_genomes=globals().get('bg_genomes', []) or [],
        fg_prefixes=globals().get('fg_prefixes', []) or [],
        bg_prefixes=globals().get('bg_prefixes', []) or [],
        fg_seq_lengths=globals().get('fg_seq_lengths', []) or [],
        bg_seq_lengths=globals().get('bg_seq_lengths', []) or [],
        excl_genomes=globals().get('excl_genomes', []) or [],
        excl_prefixes=globals().get('excl_prefixes', []) or [],
        excl_threshold=globals().get('excl_threshold', 0),
        bl_genomes=globals().get('bl_genomes', []) or [],
        bl_prefixes=globals().get('bl_prefixes', []) or [],
        bl_seq_lengths=globals().get('bl_seq_lengths', []) or [],
        bl_penalty=globals().get('bl_penalty', 5.0),
        max_bl_freq=globals().get('max_bl_freq', 0.0),
        cpus=globals().get('cpus', 1),
        verbose=globals().get('verbose', False),
    )


def set_from_config(config: PipelineParameters) -> None:
    """
    Set module globals from a PipelineParameters instance.

    This is the inverse of get_current_config(). Use it to apply
    a typed configuration to the global state for backward compatibility.

    Args:
        config: PipelineParameters instance with desired settings

    Example:
        config = PipelineParameters(min_k=8, max_k=15, polymerase='equiphi29')
        set_from_config(config)
        # Now all modules using parameter globals will see new values
    """
    g = globals()

    # K-mer range
    g['min_k'] = config.min_k
    g['max_k'] = config.max_k
    g['long_primer_mode'] = config.long_primer_mode

    # Filtering thresholds
    g['min_fg_freq'] = config.min_fg_freq
    g['max_bg_freq'] = config.max_bg_freq
    g['min_tm'] = config.min_tm
    g['max_tm'] = config.max_tm
    g['max_gini'] = config.max_gini
    g['max_primer'] = config.max_primer
    g['min_amp_pred'] = config.min_amp_pred

    # Dimer constraints
    g['max_dimer_bp'] = config.max_dimer_bp
    g['max_self_dimer_bp'] = config.max_self_dimer_bp

    # GC content filtering
    g['gc_min'] = config.gc_min
    g['gc_max'] = config.gc_max
    g['gc_tolerance'] = config.gc_tolerance
    g['genome_gc'] = config.genome_gc

    # Polymerase and reaction conditions
    g['polymerase'] = config.polymerase
    g['reaction_temp'] = config.reaction_temp
    g['na_conc'] = config.na_conc
    g['mg_conc'] = config.mg_conc
    g['primer_conc'] = config.primer_conc

    # Thermodynamic additives
    g['dmso_percent'] = config.dmso_percent
    g['betaine_m'] = config.betaine_m
    g['trehalose_m'] = config.trehalose_m
    g['formamide_percent'] = config.formamide_percent
    g['ethanol_percent'] = config.ethanol_percent
    g['urea_m'] = config.urea_m
    g['tmac_m'] = config.tmac_m

    # Optimization settings
    g['iterations'] = config.iterations
    g['max_sets'] = config.max_sets
    g['retries'] = config.retries
    g['top_set_count'] = config.top_set_count
    g['drop_iterations'] = config.drop_iterations
    g['selection_metric'] = config.selection_metric

    # Genome topology
    g['fg_circular'] = config.fg_circular
    g['bg_circular'] = config.bg_circular

    # K-mer sampling
    g['sample_rate'] = config.sample_rate
    g['min_sample_count'] = config.min_sample_count

    # Bloom filter
    g['use_bloom_filter'] = config.use_bloom_filter
    g['bloom_filter_path'] = config.bloom_filter_path

    # Paths
    g['data_dir'] = config.data_dir
    g['src_dir'] = config.src_dir

    # Genome data
    g['fg_genomes'] = config.fg_genomes
    g['bg_genomes'] = config.bg_genomes
    g['fg_prefixes'] = config.fg_prefixes
    g['bg_prefixes'] = config.bg_prefixes
    g['fg_seq_lengths'] = config.fg_seq_lengths
    g['bg_seq_lengths'] = config.bg_seq_lengths

    # Exclusion genomes
    g['excl_genomes'] = config.excl_genomes
    g['excl_prefixes'] = config.excl_prefixes
    g['excl_threshold'] = config.excl_threshold

    # Blacklist genomes
    g['bl_genomes'] = config.bl_genomes
    g['bl_prefixes'] = config.bl_prefixes
    g['bl_seq_lengths'] = config.bl_seq_lengths
    g['bl_penalty'] = config.bl_penalty
    g['max_bl_freq'] = config.max_bl_freq

    # Runtime
    g['cpus'] = config.cpus
    g['verbose'] = config.verbose


def reset_to_defaults() -> None:
    """
    Reset all parameters to their default values.

    Useful for testing or when starting a fresh pipeline run.
    """
    set_from_config(PipelineParameters())


# =============================================================================
# Legacy Global Variables
# =============================================================================
# These are preserved for backward compatibility with existing code.
# New code should use get_current_config() instead.

# Default k-mer range parameters (can be overridden by CLI or JSON)
# These are module-level so they can be set before get_params() is called
# Supported range: 6-18bp (15-18bp requires long_primer_mode optimizations)
min_k = 6
max_k = 12

# Long primer mode: enables optimizations for 15-18bp primers
# Auto-enabled when max_k >= 15. Includes k-mer sampling, GPU acceleration,
# and memory-efficient encoding for the larger k-mer space.
long_primer_mode = False

# K-mer sampling parameters for faster scoring
# sample_rate: fraction of k-mers to use (1.0 = all k-mers, 0.1 = 10%)
# min_sample_count: always include k-mers with count >= this value
sample_rate = None  # None means disabled (use all k-mers)
min_sample_count = 5

# Bloom filter parameters for large background genome filtering
use_bloom_filter = False
bloom_filter_path = None

# Polymerase and reaction condition parameters
# polymerase: phi29 (30C), equiphi29 (42-45C), bst (60-65C), klenow (37C)
polymerase = 'phi29'
reaction_temp = None  # Auto-set based on polymerase if not specified
na_conc = 50.0   # mM
mg_conc = 2.0    # mM
primer_conc = 0.5e-6  # M (0.5 uM)

# Thermodynamic additives
dmso_percent = 0.0
betaine_m = 0.0
trehalose_m = 0.0
formamide_percent = 0.0
ethanol_percent = 0.0
urea_m = 0.0
tmac_m = 0.0

# Exclusion genome parameters (e.g., mitochondrial DNA, chloroplast sequences)
excl_genomes = []
excl_prefixes = []
excl_threshold = 0

# Blacklist genome parameters (penalty-weighted filtering)
bl_genomes = []
bl_prefixes = []
bl_seq_lengths = []
bl_penalty = 5.0
max_bl_freq = 0.0

# GC content filtering parameters
# gc_tolerance: adaptive filtering uses genome_gc +/- gc_tolerance
gc_tolerance = 0.15
gc_min = 0.375
gc_max = 0.625

def get_all_files(prefix_or_file):
    """
    If passed a directory path, this returns all the fasta files in a directory. Otherwise, if given a file it
    passes back the file if it is a fasta file.

    Args:
        prefix_or_file: Directory path or path to a fasta file.

    Returns:
        total_files: List of fasta file(s).
    """
    if type(prefix_or_file) is str:
        if 'fasta' == prefix_or_file.split('.')[-1] or 'fa' == prefix_or_file.split('.')[-1]:
            return [prefix_or_file]
        else:
            prefix_or_file = [prefix_or_file]
    total_files = []
    for individual_prefix in prefix_or_file:
        list_of_files = os.listdir(os.path.dirname(individual_prefix))  # list of files in the current directory
        for each_file in list_of_files:
            if each_file.startswith(individual_prefix.split('/')[-1]):
                if each_file.endswith('fasta') or each_file.endswith('fa'):
                    total_files.append(os.path.join(os.path.dirname(individual_prefix), each_file))
    return total_files

def get_value_or_default(arg_value, data, key):
    """
    Helper function for getting the input values. If a command line argument value was given, just return arg_value.
    Else we grab the value from the dictionary parsed from the json. This is to ensure command line inputs overwrite
    json inputs.

    Args:
        arg_value: Argument value to be returned
        data: dictionary of parameter to parameter values read from the input json file.
        key: The name of the parameter.

    Returns:
        val: Returns the command line input value or the json input value if command line one wasn't entered.
    """
    OPTIONAL_PARAMS = {
        'src_dir', 'min_amp_pred', 'max_dimer_bp', 'max_self_dimer_bp',
        'verbose', 'drop_iterations', 'top_set_count', 'retries',
        'selection_metric', 'fg_circular', 'bg_circular', 'mismatch_penalty'
    }
    if arg_value is not None:
        return arg_value
    if key in data:
        return data[key]
    else:
        if key in OPTIONAL_PARAMS:
            logger.debug("Missing optional parameter '%s' in params.json, using default.", key)
        else:
            logger.warning("Missing parameter '%s' in params.json. Please provide a value.", key)
        return None

def get_params(args):
    """
    Writes the arguments of a pipeline instance to a json file for future use.

    Args:
        args: All the command line parameter inputs.

    Returns:
        data: Returns all the parameter inputs in dictionary form.
    """
    global min_fg_freq
    global max_bg_freq
    global min_tm
    global max_tm
    global max_gini
    global max_primer
    global min_amp_pred
    global cpus
    global max_dimer_bp
    global max_self_dimer_bp
    global mismatch_penalty
    global data_dir
    global selection_metric
    global iterations
    global max_sets
    global fg_circular
    global bg_circular
    global drop_iterations
    global verbose
    global top_set_count
    global retries
    global src_dir
    global genome_gc
    global fg_genomes
    global bg_genomes
    global fg_prefixes
    global bg_prefixes
    global gc_min
    global gc_max
    global gc_tolerance
    global min_k
    global max_k
    global long_primer_mode
    global sample_rate
    global min_sample_count
    global use_bloom_filter
    global bloom_filter_path
    global polymerase
    global reaction_temp
    global na_conc
    global mg_conc
    global primer_conc
    global dmso_percent
    global betaine_m
    global trehalose_m
    global formamide_percent
    global ethanol_percent
    global urea_m
    global tmac_m
    global num_primers
    global target_set_size
    global _json_data
    global excl_genomes
    global excl_prefixes
    global excl_threshold
    global bl_genomes
    global bl_prefixes
    global bl_seq_lengths
    global bl_penalty
    global max_bl_freq


    data = {}

    # data['json_file'] = get_value_or_default(args.json_file, data, 'json_file', os.path.join(data['data_dir'], 'params')) if out_fname is None else out_fname

    if args.json_file is not None:
        data['json_file'] = args.json_file
        json_dir = os.path.dirname(data['json_file'])
        # Only create directory if it's not empty (i.e., json_file has a path component)
        if json_dir and not os.path.exists(json_dir):
            os.makedirs(json_dir)
        # Check if json file exists and read it
        if os.path.isfile(data['json_file']):
            data_extra = read_args_from_json(args.json_file)
            _json_data.update(data_extra)  # Store raw JSON for CLI access

            # Schema versioning
            CURRENT_SCHEMA_VERSION = 1
            schema_version = data_extra.get('schema_version', None) if isinstance(data_extra, dict) else None
            if schema_version is None:
                logger.warning(
                    "params.json has no 'schema_version' field. "
                    "Defaults may differ between NeoSWGA versions. "
                    "Add '\"schema_version\": 1' to your params.json for reproducibility."
                )
            elif schema_version > CURRENT_SCHEMA_VERSION:
                logger.warning(
                    f"params.json schema_version {schema_version} is newer than "
                    f"this NeoSWGA version supports (max: {CURRENT_SCHEMA_VERSION}). "
                    f"Some parameters may not be recognized."
                )

            for k, v in data_extra.items():
                if k not in data:
                    data[k] = v

    # Load optimization parameters from JSON
    num_primers = data.get('num_primers', data.get('target_set_size', 6))
    target_set_size = data.get('target_set_size', num_primers)

    data_dir = data['data_dir'] = get_value_or_default(args.data_dir, data, 'data_dir')
    src_dir = data['src_dir'] = get_value_or_default(args.src_dir, data, 'src_dir')
    min_fg_freq = data['min_fg_freq'] = get_value_or_default(args.min_fg_freq, data, 'min_fg_freq')
    max_bg_freq = data['max_bg_freq'] = get_value_or_default(args.max_bg_freq, data, 'max_bg_freq')

    min_tm = data['min_tm'] = get_value_or_default(args.min_tm, data, 'min_tm')
    max_tm = data['max_tm'] = get_value_or_default(args.max_tm, data, 'max_tm')



    max_gini = data['max_gini'] = get_value_or_default(args.max_gini, data, 'max_gini')
    max_primer = data['max_primer'] = get_value_or_default(args.max_primer, data, 'max_primer')
    min_amp_pred = data['min_amp_pred'] = get_value_or_default(args.min_amp_pred, data, 'min_amp_pred')
    cpus = data['cpus'] = get_value_or_default(args.cpus, data, 'cpus')
    max_dimer_bp = data['max_dimer_bp'] = get_value_or_default(args.max_dimer_bp, data, 'max_dimer_bp')
    max_self_dimer_bp = data['max_self_dimer_bp'] = get_value_or_default(args.max_self_dimer_bp, data, 'max_self_dimer_bp')
    verbose = data['verbose'] = get_value_or_default(args.verbose, data, 'verbose')

    drop_iterations = data['drop_iterations'] = get_value_or_default(args.drop_iterations, data, 'drop_iterations')
    iterations = data['iterations'] = get_value_or_default(args.iterations, data, 'iterations')
    top_set_count = data['top_set_count'] = get_value_or_default(args.top_set_count, data, 'top_set_count')
    retries = data['retries'] = get_value_or_default(args.retries, data, 'retries')
    max_sets = data['max_sets'] = get_value_or_default(args.max_sets, data, 'max_sets')
    selection_metric = data['selection_metric'] = get_value_or_default(args.selection_metric, data, 'selection_metric')
    fg_circular = data['fg_circular'] = get_value_or_default(args.fg_circular, data, 'fg_circular')
    bg_circular = data['bg_circular'] = get_value_or_default(args.bg_circular, data, 'bg_circular')

    # K-mer range parameters (with defaults)
    # Supported range: 6-18bp
    if 'min_k' not in data:
        data['min_k'] = 6
    if 'max_k' not in data:
        data['max_k'] = 12
    min_k = data['min_k'] = get_value_or_default(args.min_k if hasattr(args, 'min_k') else None, data, 'min_k')
    max_k = data['max_k'] = get_value_or_default(args.max_k if hasattr(args, 'max_k') else None, data, 'max_k')

    # Validate k-mer range
    if min_k < 4:
        logger.warning("min_k=%d is very small. Minimum recommended is 6.", min_k)
    elif min_k < 6:
        logger.warning("min_k=%d is below recommended minimum of 6bp.", min_k)
    if max_k > 31:
        logger.error("max_k=%d exceeds jellyfish maximum of 31. Setting to 31.", max_k)
        max_k = 31
        data['max_k'] = 31
    elif max_k > 18:
        logger.warning("max_k=%d exceeds supported maximum of 18bp.", max_k)
    if min_k > max_k:
        raise ValueError(f"min_k ({min_k}) must be <= max_k ({max_k})")

    # Auto-scale dimer thresholds based on max_k: a fixed 3bp match is
    # meaningful for 6bp primers (50%) but negligible for 18bp primers (17%).
    # Only auto-scale if the user did not explicitly set these values.
    if 'max_dimer_bp' not in _json_data and args.max_dimer_bp is None:
        max_dimer_bp = data['max_dimer_bp'] = max(3, int(max_k * 0.25))
    if 'max_self_dimer_bp' not in _json_data and args.max_self_dimer_bp is None:
        max_self_dimer_bp = data['max_self_dimer_bp'] = max(4, int(max_k * 0.28))

    # Long primer mode: auto-enable for k >= 15
    # Enables k-mer sampling, GPU acceleration, and memory optimizations
    long_primer_mode = data.get('long_primer_mode', False)
    if max_k >= 15 and not long_primer_mode:
        long_primer_mode = True
        if verbose:
            print(f"Auto-enabling long_primer_mode for {max_k}bp primers")
    data['long_primer_mode'] = long_primer_mode

    # K-mer sampling parameters for faster scoring
    # Auto-configure for long primers if not explicitly set
    sample_rate = data.get('sample_rate', None)
    min_sample_count = data.get('min_sample_count', 5)

    # Auto-enable sampling for long primers (15-18bp)
    if long_primer_mode and sample_rate is None:
        # Scale sampling rate with k-mer length:
        # 15bp: 5%, 16bp: 3%, 17bp: 2.5%, 18bp: 2%
        if max_k >= 18:
            sample_rate = 0.02
        elif max_k >= 17:
            sample_rate = 0.025
        elif max_k >= 16:
            sample_rate = 0.03
        else:
            sample_rate = 0.05
        if verbose:
            print(f"Auto-setting sample_rate={sample_rate:.1%} for {max_k}bp primers")
        data['sample_rate'] = sample_rate

    # Bloom filter parameters for large background genome filtering (optional)
    use_bloom_filter = data.get('use_bloom_filter', False)
    bloom_filter_path = data.get('bloom_filter_path', None)

    # Polymerase and reaction condition parameters (optional)
    # Load from JSON file if present, otherwise use module defaults
    polymerase = data.get('polymerase', 'phi29')
    na_conc = data.get('na_conc', 50.0)
    mg_conc = data.get('mg_conc', 2.0)
    primer_conc = data.get('primer_conc', 0.5e-6)

    # Auto-set reaction temperature based on polymerase if not specified
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
    if 'reaction_temp' in data and data['reaction_temp'] is not None:
        reaction_temp = data['reaction_temp']
    elif polymerase.lower() in POLYMERASE_CHARACTERISTICS:
        reaction_temp = POLYMERASE_CHARACTERISTICS[polymerase.lower()]['optimal_temp']
    else:
        reaction_temp = 30.0  # Default to phi29 temperature

    # Thermodynamic additives (optional)
    dmso_percent = data.get('dmso_percent', 0.0)
    betaine_m = data.get('betaine_m', 0.0)
    trehalose_m = data.get('trehalose_m', 0.0)
    formamide_percent = data.get('formamide_percent', 0.0)
    ethanol_percent = data.get('ethanol_percent', 0.0)
    urea_m = data.get('urea_m', 0.0)
    tmac_m = data.get('tmac_m', 0.0)

    # Store all reaction parameters in data dict
    data['polymerase'] = polymerase
    data['reaction_temp'] = reaction_temp
    data['na_conc'] = na_conc
    data['mg_conc'] = mg_conc
    data['primer_conc'] = primer_conc
    data['dmso_percent'] = dmso_percent
    data['betaine_m'] = betaine_m
    data['trehalose_m'] = trehalose_m
    data['formamide_percent'] = formamide_percent
    data['ethanol_percent'] = ethanol_percent
    data['urea_m'] = urea_m
    data['tmac_m'] = tmac_m

    if args.fasta_fore is not None:
        data['fg_genomes'] = get_all_files(args.fasta_fore)
    if args.fasta_back is not None:
        data['bg_genomes'] = get_all_files(args.fasta_back)

    if args.kmer_fore is not None:
        if len(data['fg_genomes']) > 1:
            data['fg_prefixes'] = [args.kmer_fore + os.path.splitext(os.path.basename(genome_fname))[0] for genome_fname
                                   in data['fg_genomes']]
        else:
            data['fg_prefixes'] = [args.kmer_fore]
    if args.kmer_back is not None:
        if len(data['bg_genomes']) > 1:
            data['bg_prefixes'] = [args.kmer_back + '_' + os.path.splitext(os.path.basename(genome_fname))[0] for
                                   genome_fname in data['bg_genomes']]
        else:
            data['bg_prefixes'] = [args.kmer_back]

    # Set global variables for genomes and prefixes
    fg_genomes = data.get('fg_genomes', [])
    bg_genomes = data.get('bg_genomes', [])
    fg_prefixes = data.get('fg_prefixes', [])
    bg_prefixes = data.get('bg_prefixes', [])

    # Exclusion genome support
    excl_genomes = data.get('excl_genomes', [])
    excl_prefixes = data.get('excl_prefixes', [])
    excl_threshold = data.get('excl_threshold', 0)

    # Blacklist genome support
    bl_genomes = data.get('bl_genomes', [])
    bl_prefixes = data.get('bl_prefixes', [])
    bl_seq_lengths = data.get('bl_seq_lengths', [])
    bl_penalty = data.get('bl_penalty', 5.0)
    max_bl_freq = data.get('max_bl_freq', 0.0)

    # Resolve genome library names for bg_genomes and bl_genomes
    try:
        from neoswga.core.genome_library import GenomeLibrary
        _library = GenomeLibrary()
    except Exception as e:
        logger.debug(f"Ignored error loading genome library: {e}")
        _library = None

    if _library is not None:
        # Resolve background genome names
        resolved_bg = []
        resolved_bg_prefixes = list(bg_prefixes) if bg_prefixes else []
        for bg in bg_genomes:
            if not os.path.isfile(bg):
                entry = _library.get(bg)
                if entry:
                    logger.info(f"Resolved '{bg}' from genome library: {entry.fasta_path}")
                    resolved_bg.append(entry.fasta_path)
                    if entry.kmer_prefix and entry.kmer_prefix not in resolved_bg_prefixes:
                        resolved_bg_prefixes.append(entry.kmer_prefix)
                else:
                    resolved_bg.append(bg)  # Keep as-is, will fail later with clear error
            else:
                resolved_bg.append(bg)
        if resolved_bg != list(bg_genomes):
            bg_genomes = resolved_bg
            data['bg_genomes'] = bg_genomes
            bg_prefixes = resolved_bg_prefixes
            data['bg_prefixes'] = bg_prefixes

        # Resolve blacklist genome names
        resolved_bl = []
        resolved_bl_prefixes = list(bl_prefixes) if bl_prefixes else []
        for bl in bl_genomes:
            if not os.path.isfile(bl):
                entry = _library.get(bl)
                if entry:
                    logger.info(f"Resolved '{bl}' from genome library: {entry.fasta_path}")
                    resolved_bl.append(entry.fasta_path)
                    if entry.kmer_prefix and entry.kmer_prefix not in resolved_bl_prefixes:
                        resolved_bl_prefixes.append(entry.kmer_prefix)
                else:
                    resolved_bl.append(bl)
            else:
                resolved_bl.append(bl)
        if resolved_bl != list(bl_genomes):
            bl_genomes = resolved_bl
            data['bl_genomes'] = bl_genomes
            bl_prefixes = resolved_bl_prefixes
            data['bl_prefixes'] = bl_prefixes

    if 'fg_genomes' in data and ('fg_seq_lengths' not in data or len(data['fg_seq_lengths']) != len(data['fg_genomes'])):
        data['fg_seq_lengths'] = _utility.get_all_seq_lengths(fname_genomes=data['fg_genomes'], cpus=data['cpus'])

    # In Bloom filter mode, bg_prefixes is empty - ensure bg_seq_lengths matches
    # This allows optimization to proceed without background position data
    if use_bloom_filter and len(data.get('bg_prefixes', [])) == 0:
        data['bg_seq_lengths'] = []
    elif 'bg_genomes' in data and ('bg_seq_lengths' not in data or len(data['bg_seq_lengths']) != len(data['bg_genomes'])):
        data['bg_seq_lengths'] = _utility.get_all_seq_lengths(fname_genomes=data['bg_genomes'], cpus=data['cpus'])

    # Calculate genome GC content from foreground genomes if not already provided
    # This enables adaptive GC filtering for extreme GC genomes
    if 'genome_gc' not in data or data['genome_gc'] is None:
        if 'fg_genomes' in data and len(data['fg_genomes']) > 0:
            try:
                gc_count = 0
                total_length = 0

                # Use genome_io for automatic gzip/zip detection
                import neoswga.core.genome_io as genome_io
                loader = genome_io.GenomeLoader()

                for fg_genome in data['fg_genomes']:
                    # Load genome with automatic compression detection
                    sequence = loader.load_genome(fg_genome, return_stats=False)
                    seq_upper = sequence.upper()
                    gc_count += seq_upper.count('G') + seq_upper.count('C')
                    total_length += len(seq_upper)

                if total_length > 0:
                    data['genome_gc'] = gc_count / total_length
                else:
                    data['genome_gc'] = None
            except Exception as e:
                logger.debug(f"Ignored error calculating genome GC: {e}")
                data['genome_gc'] = None
        else:
            data['genome_gc'] = None

    genome_gc = data.get('genome_gc', None)

    # Get gc_tolerance for adaptive filtering (default 0.15 = +/-15%)
    gc_tolerance = data.get('gc_tolerance', 0.15)

    # Adaptive GC filtering: compute gc_min and gc_max based on genome_gc
    # If genome_gc is known, use genome_gc +/- gc_tolerance
    # This enables extreme GC genomes like Francisella (33%) and Burkholderia (67%)
    if genome_gc is not None and genome_gc > 0:
        # Adaptive: primer GC should match target genome GC +/- tolerance
        gc_min = max(0.15, genome_gc - gc_tolerance)
        gc_max = min(0.85, genome_gc + gc_tolerance)
    else:
        # Fallback to explicit gc_min/gc_max if specified, otherwise use defaults
        gc_min = data.get('gc_min', 0.375)
        gc_max = data.get('gc_max', 0.625)

    # Store computed values in data dict for downstream use
    data['gc_min'] = gc_min
    data['gc_max'] = gc_max
    data['gc_tolerance'] = gc_tolerance

    return data

def read_args_from_json(in_fname):
    """
    Returns the parameters inputs in json form as a dictionary.

    Args:
        in_fname: Path to the json file.

    Returns:
        data: Parameter inputs in dictionary form.

    """
    with open(in_fname, 'r') as json_file:
        data = json.load(json_file)
    return data