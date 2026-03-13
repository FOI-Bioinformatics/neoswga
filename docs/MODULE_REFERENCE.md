# NeoSWGA Module Reference

Complete reference for all modules in the NeoSWGA core package.

## Table of Contents

1. [Pipeline Modules](#pipeline-modules)
2. [Filtering Modules](#filtering-modules)
3. [Scoring Modules](#scoring-modules)
4. [Optimization Modules](#optimization-modules)
5. [Thermodynamic Modules](#thermodynamic-modules)
6. [Performance Modules](#performance-modules)
7. [Simulation Modules](#simulation-modules)
8. [Analysis Modules](#analysis-modules)
9. [Utility Modules](#utility-modules)

---

## Pipeline Modules

### pipeline.py

Step prerequisite validation and workflow control.

**Classes:**
- `StepPrerequisiteError`: Exception for missing prerequisites
- `StepValidationResult`: Validation result dataclass

**Functions:**
```python
def validate_step_prerequisites(step: str, params: dict) -> StepValidationResult:
    """Validate that prerequisites for a step are met."""
```

### improved_pipeline.py

Performance-optimized pipeline with modern features.

**Key Features:**
- Position cache integration
- Adaptive GC filtering
- Bloom filter support
- Network-based optimization

### auto_swga_pipeline.py

Automatic parameter optimization pipeline.

**Functions:**
```python
def run_auto_pipeline(
    fg_genome: str,
    bg_genome: str,
    output_dir: str,
    **kwargs
) -> dict:
    """Run pipeline with automatic parameter tuning."""
```

### multi_genome_pipeline.py

Pan-genome primer design for multiple targets.

**Functions:**
```python
def design_pan_primers(
    target_genomes: List[str],
    background_genomes: List[str],
    output_dir: str,
    **kwargs
) -> dict:
    """Design primers that work across multiple target genomes."""
```

### pipeline_qa_integration.py

Quality assurance hooks for the filtering step.

**Classes:**
- `QAValidator`: Validates primer candidates
- `QAReport`: Quality report dataclass

---

## Filtering Modules

### filter.py

Core filtering logic with five sequence rules.

**Functions:**
```python
def filter_extra(primer: str) -> bool:
    """
    Apply all five filtering rules:
    1. No 3+ consecutive G/C at 3' end
    2. GC content 40-60%
    3. GC clamp: max 3 G/C in last 5 bases
    4. Max 4 consecutive di-nucleotide repeats
    5. Max 4 bp homopolymer runs
    """

def get_bg_rates_via_bloom(
    primer_list: List[str],
    bloom_path: str
) -> Dict[str, int]:
    """Get background rates using Bloom filter."""
```

### adaptive_filters.py

Adaptive GC filtering for extreme GC genomes.

**Functions:**
```python
def get_adaptive_gc_bounds(
    genome_gc: float,
    stringency: str = 'normal'
) -> Tuple[float, float]:
    """
    Calculate GC bounds based on target genome composition.

    Args:
        genome_gc: Target genome GC content (0-1)
        stringency: 'relaxed', 'normal', or 'strict'

    Returns:
        (min_gc, max_gc) bounds for primer filtering
    """
```

### kmer.py

K-mer counting via Jellyfish subprocess.

**Functions:**
```python
def count_kmers(
    genome_path: str,
    k: int,
    output_prefix: str,
    threads: int = 8
) -> str:
    """
    Count k-mers using Jellyfish.

    Returns:
        Path to output k-mer file
    """

def generate_position_file(
    genome_path: str,
    kmer_file: str,
    output_h5: str
) -> None:
    """Generate HDF5 file with k-mer binding positions."""
```

### kmer_counter.py

Multi-genome k-mer counter with batch processing.

**Classes:**
```python
class MultiGenomeKmerCounter:
    """Count k-mers across multiple genomes efficiently."""

    def count_all(
        self,
        genomes: List[str],
        k_range: Tuple[int, int],
        output_dir: str
    ) -> Dict[str, str]:
        """Count k-mers for all genomes and k values."""
```

### thermodynamic_filter.py

Filter by thermodynamic properties.

**Functions:**
```python
def filter_by_tm(
    primers: List[str],
    min_tm: float,
    max_tm: float,
    conditions: ReactionConditions
) -> List[str]:
    """Filter primers by melting temperature range."""

def filter_by_secondary_structure(
    primers: List[str],
    max_hairpin_tm: float,
    max_dimer_severity: float,
    conditions: ReactionConditions
) -> List[str]:
    """Filter primers by secondary structure stability."""
```

### thermo_background_filter.py

Thermodynamic-aware background filtering.

**Functions:**
```python
def filter_background_thermodynamic(
    primers: List[str],
    bg_genome: str,
    conditions: ReactionConditions,
    max_bg_binding_prob: float = 0.01
) -> List[str]:
    """Filter by thermodynamic binding probability to background."""
```

### multi_genome_filter.py

Filter primers across multiple target genomes.

**Functions:**
```python
def filter_multi_genome(
    primers: List[str],
    target_genomes: List[str],
    min_coverage_fraction: float = 0.8
) -> List[str]:
    """Keep primers that bind to specified fraction of targets."""
```

---

## Scoring Modules

### rf_preprocessing.py

Random Forest feature engineering.

**Functions:**
```python
def load_model_safely(
    model_path: str,
    verify_hash: bool = True
) -> Any:
    """
    Load pickle model with SHA-256 hash verification.

    Raises:
        ModelIntegrityError: If hash doesn't match
    """

def compute_features(
    primer: str,
    fg_positions: np.ndarray,
    bg_positions: np.ndarray,
    fg_length: int,
    bg_length: int,
    conditions: ReactionConditions
) -> np.ndarray:
    """
    Compute feature vector for primer.

    Returns:
        Array of 52-120 features depending on mode
    """
```

**Feature Categories:**
1. Thermodynamic histograms (bins of Tm distribution)
2. Positional gap statistics
3. GC content distribution
4. Gini index components
5. Secondary structure scores

### primer_attributes.py

Primer property calculations.

**Functions:**
```python
def calculate_tm(
    primer: str,
    conditions: ReactionConditions
) -> float:
    """Calculate effective Tm under reaction conditions."""

def calculate_self_complementarity(
    primer: str
) -> Tuple[int, float]:
    """
    Calculate self-complementarity.

    Returns:
        (max_consecutive_bp, severity_score)
    """

def calculate_complexity(primer: str) -> float:
    """
    Calculate sequence complexity (0-1).

    Low complexity indicates repetitive sequence.
    """
```

### integrated_quality_scorer.py

Multi-criteria quality scoring.

**Classes:**
```python
class IntegratedQualityScorer:
    """Combine multiple scoring criteria."""

    def score_primer(
        self,
        primer: str,
        positions: np.ndarray,
        genome_length: int
    ) -> float:
        """
        Score primer quality (0-1).

        Combines: coverage, specificity, Tm, secondary structure
        """

    def score_primer_set(
        self,
        primers: List[str],
        position_cache: PositionCache,
        fg_prefix: str,
        genome_length: int
    ) -> Dict[str, float]:
        """
        Score primer set with multiple metrics.

        Returns:
            Dict with coverage, enrichment, uniformity, etc.
        """
```

---

## Optimization Modules

### base_optimizer.py

Abstract base class for all optimizers.

**Classes:**
```python
@dataclass
class OptimizerConfig:
    """Common optimizer configuration."""
    target_size: int = 6
    max_iterations: int = 100
    timeout_seconds: float = 300.0
    min_coverage: float = 0.5

@dataclass(frozen=True)
class OptimizationResult:
    """Optimization result container."""
    primer_set: List[str]
    metrics: PrimerSetMetrics
    status: OptimizationStatus
    message: str

class BaseOptimizer(ABC):
    """Abstract base for optimizers."""

    @abstractmethod
    def optimize(
        self,
        candidates: List[str],
        target_size: int,
        **kwargs
    ) -> OptimizationResult:
        """Find optimal primer set."""
```

### optimizer_factory.py

Factory for optimizer instantiation.

**Classes:**
```python
class OptimizerFactory:
    """Create optimizer instances."""

    @staticmethod
    def create(
        name: str,
        cache: PositionCache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        **kwargs
    ) -> BaseOptimizer:
        """Create optimizer by name."""

    @staticmethod
    def list_optimizers() -> Dict[str, str]:
        """List available optimizers."""

class OptimizerRegistry:
    """Thread-safe optimizer class registry."""

    @classmethod
    def register(
        cls,
        name: str,
        aliases: Optional[List[str]] = None,
        description: Optional[str] = None
    ) -> Callable:
        """Decorator to register optimizer class."""
```

### greedy_optimizer.py

Breadth-first search optimization.

**Algorithm:**
1. Start with highest-scoring primer
2. Iteratively add primer that maximizes coverage gain
3. Check dimer compatibility at each step
4. Stop at target size

### network_optimizer.py

Network-based optimization with Tm weighting.

**Key Insight:** SWGA efficiency depends on network connectivity between binding sites, not just frequency.

**Classes:**
```python
@dataclass
class BindingSite:
    position: int
    strand: str
    primer: str
    affinity: float

class AmplificationNetwork:
    """Model primer binding network."""

    def build(self, primers: List[str], positions: Dict) -> None:
        """Build bipartite graph of primer-target positions."""

    def score_connectivity(self) -> float:
        """Score network connectivity for amplification."""
```

### genetic_algorithm.py

Evolutionary optimization.

**Classes:**
```python
@dataclass
class GAConfig:
    population_size: int = 200
    generations: int = 100
    mutation_rate: float = 0.15
    crossover_rate: float = 0.8
    elitism_fraction: float = 0.10
    tournament_size: int = 5
    min_set_size: int = 4
    max_set_size: int = 8

class PrimerSetGA:
    """Genetic algorithm for primer set selection."""

    def evolve(self, verbose: bool = True) -> Individual:
        """Run evolution and return best individual."""
```

**Operators:**
- Selection: Tournament (k=5)
- Crossover: Uniform with dimer checking
- Mutation: Add/remove/replace
- Elitism: Preserve top 10%

### hybrid_optimizer.py

Combines network and greedy approaches.

**Phases:**
1. Network coverage optimization
2. Greedy refinement
3. Dimer optimization

### dominating_set_optimizer.py

Graph-based set cover optimization (8x faster).

**Algorithm:**
1. Divide genome into bins
2. Build coverage graph (primer -> bins covered)
3. Greedy set cover with ln(n) approximation
4. Post-process for dimer avoidance

### background_aware_optimizer.py

Three-stage clinical optimization.

**Stages:**
1. Maximize target coverage
2. Minimize background binding
3. Avoid primer dimers

**Result:** 10-20x background reduction compared to standard methods.

### milp_optimizer.py

Mixed-integer linear programming (exact solutions).

**Requirements:** `mip` package

**Formulation:**
- Binary variables: x_i (primer i selected)
- Objective: minimize sum(x_i) subject to coverage constraint
- Constraints: coverage >= threshold, dimers avoided

### moea_optimizer.py

Multi-objective evolutionary algorithm.

**Objectives:**
1. Maximize coverage
2. Maximize specificity
3. Minimize primer count

**Output:** Pareto frontier of non-dominated solutions

### equiphi29_optimizer.py

EquiPhi29-specific optimization at 42-45C.

**Features:**
- Longer primer support (12-18 bp)
- Enhanced thermodynamic modeling
- Higher temperature constraints

### normalized_optimizer.py

Normalized scoring with strategy presets.

**Purpose:** Provides a normalized [0,1] composite score for comparing results across optimizers. Includes strategy presets (discovery, clinical, enrichment, metagenomics) that adjust scoring weights.

### tiling_optimizer.py

Interval-based tiling coverage optimization.

**Algorithm:**
1. Model primer binding as genomic intervals
2. Select primers to tile the genome with minimal gaps
3. Merge overlapping intervals to compute uncovered regions

### clique_optimizer.py

Clique-based dimer-free primer set selection.

**Purpose:** Builds a primer compatibility graph (edges between primers that do not form dimers) and finds maximum cliques to guarantee dimer-free sets.

**Classes:**
```python
class CliqueOptimizer(BaseOptimizer):
    """Find dimer-free primer sets using maximum clique enumeration."""
```

### dimer_validator.py

Post-optimization dimer validation.

**Classes:**
```python
class DimerValidator:
    """Validate primer sets for dimer interactions and suggest replacements."""
```

### background_prefilter.py

Background-aware candidate pruning optimizer.

**Purpose:** Wraps another optimizer, first pruning candidates with poor foreground/background ratios before delegating to the inner optimizer.

### serial_cascade_optimizer.py

Serial pipeline combinations of optimizers.

**Classes:**
```python
class CoverageThenDimerFreeOptimizer(SerialCascadeOptimizer):
    """Dominating-set coverage followed by clique dimer removal."""

class DimerFreeScoredOptimizer(SerialCascadeOptimizer):
    """Clique dimer-free selection followed by network scoring."""

class BgPrefilterHybridOptimizer(SerialCascadeOptimizer):
    """Background pre-filter followed by hybrid optimization."""
```

### multi_agent_optimizer.py

Multi-agent parallel optimizer execution.

**Purpose:** Runs multiple optimizer strategies concurrently and aggregates results using configurable strategies (best, union, voting, pareto).

### dominating_set_adapter.py

BaseOptimizer adapter for dominating set and weighted set cover.

**Classes:**
```python
class DominatingSetAdapter(BaseOptimizer):
    """Adapter wrapping DominatingSetOptimizer to BaseOptimizer interface."""

class WeightedSetCoverOptimizer(DominatingSetAdapter):
    """Weighted set cover variant prioritizing high-scoring primers."""
```

### unified_optimizer.py

Unified entry point for running any registered optimizer from CLI or programmatic use.

**Functions:**
```python
def run_optimization(method, candidates, fg_prefixes, ...) -> OptimizationResult
def list_available_optimizers() -> Dict[str, str]
def optimize_step4(**kwargs) -> Dict
```

---

## Thermodynamic Modules

### thermodynamics.py

SantaLucia nearest-neighbor calculations.

**Constants:**
```python
R = 1.987  # cal/(mol*K)
ENTHALPY_NN = {...}  # 16 dinucleotide stacks
ENTROPY_NN = {...}
INIT_ENTHALPY = 0.2  # kcal/mol
INIT_ENTROPY = -5.7  # cal/(mol*K)
TERMINAL_AT_ENTHALPY = 2.2
TERMINAL_AT_ENTROPY = 6.9
SYMMETRY_ENTROPY = -1.4
```

**Functions:**
```python
def calculate_tm_with_salt(seq, na_conc, mg_conc, ...) -> float
def calculate_free_energy(seq, temp, ...) -> float
def gc_content(seq) -> float
def reverse_complement(seq) -> str
def is_palindrome(seq) -> bool
```

### reaction_conditions.py

Reaction condition modeling.

**Polymerase Database:**
| Name | Temp | Processivity | Use Case |
|------|------|--------------|----------|
| phi29 | 30C | 70 kb | Standard |
| equiphi29 | 42C | 80 kb | High specificity |
| bst | 63C | 2 kb | LAMP |
| klenow | 37C | 10 kb | Budget |

**Classes:**
```python
class ReactionConditions:
    temp: float
    dmso_percent: float
    betaine_m: float
    na_conc: float
    polymerase: str
    # ... other attributes

    def calculate_effective_tm(self, seq: str) -> float
    def max_primer_length(self) -> int
```

### secondary_structure.py

Hairpin and dimer prediction.

**Functions:**
```python
def check_heterodimer(
    seq1: str,
    seq2: str,
    conditions: ReactionConditions
) -> Dict:
    """
    Check for primer-primer interaction.

    Returns:
        {'energy': float, 'severity': float, 'forms_dimer': bool}
    """

def check_hairpins(seq: str) -> List[Dict]:
    """
    Find hairpin structures.

    Returns:
        List of {'stem_length': int, 'loop_size': int, 'stable': bool}
    """

def calculate_dimer_matrix(
    primers: List[str],
    conditions: ReactionConditions
) -> np.ndarray:
    """Calculate pairwise dimer severity matrix."""
```

### dimer.py

Primer-dimer calculations.

**Functions:**
```python
def is_dimer(
    primer1: str,
    primer2: str,
    max_bp: int = 4
) -> bool:
    """Check if primers form dimer with >= max_bp consecutive base pairs."""

def find_dimer_sites(
    primer1: str,
    primer2: str
) -> List[Tuple[int, int, int]]:
    """
    Find all dimer binding sites.

    Returns:
        List of (pos1, pos2, length) tuples
    """
```

### three_prime_stability.py

3' end stability analysis.

**Functions:**
```python
def calculate_3prime_stability(
    seq: str,
    n_bases: int = 5
) -> float:
    """
    Calculate 3' end stability score.

    Higher = more stable = better extension
    """

def check_3prime_mismatch_tolerance(
    primer: str,
    target: str,
    position: int
) -> float:
    """Calculate mismatch tolerance at 3' end."""
```

---

## Performance Modules

### position_cache.py

In-memory position cache (1000x speedup).

**Classes:**
```python
class PositionCache:
    def __init__(self, fname_prefixes: List[str], primers: List[str])
    def get_positions(self, fname_prefix: str, primer: str, strand: str = 'both') -> np.ndarray
    def get_all_positions(self, fname_prefix: str, primers: List[str]) -> Dict
    def compute_coverage_vectorized(self, fname_prefix: str, primers: List[str], genome_length: int) -> np.ndarray
```

**Memory:** ~4 MB for 500 primers x 1000 sites

### background_filter.py

Bloom filter for large backgrounds.

**Classes:**
```python
class BackgroundBloomFilter:
    """Probabilistic set membership."""

    @classmethod
    def build(cls, kmer_file: str, error_rate: float = 0.01) -> 'BackgroundBloomFilter'

    def contains(self, primer: str) -> bool

    def save(self, path: str) -> None

    @classmethod
    def load(cls, path: str) -> 'BackgroundBloomFilter'

class SampledGenomeIndex:
    """Sample-based count estimation."""

    def estimate_count(self, primer: str) -> int
```

### gpu_acceleration.py

Optional CuPy-based GPU acceleration.

**Functions:**
```python
def is_gpu_available() -> bool:
    """Check if GPU acceleration is available."""

def get_gpu_info() -> Dict[str, Any]:
    """Get GPU device information."""

class GPUThermodynamics:
    """GPU-accelerated thermodynamic calculations."""

    def batch_calculate_tm(self, primers: List[str]) -> np.ndarray:
        """Calculate Tm for many primers in parallel (10-100x faster)."""
```

---

## Simulation Modules

### replication_simulator.py

Agent-based DNA replication simulation.

**Classes:**
```python
class ReplicationSimulator:
    """Simulate phi29 DNA replication."""

    def simulate(
        self,
        primers: List[str],
        genome: str,
        time_points: List[float],
        conditions: ReactionConditions
    ) -> SimulationResult:
        """
        Run replication simulation.

        Returns:
            SimulationResult with yield curves and statistics
        """
```

### swga_simulator.py

SWGA reaction simulation.

**Classes:**
```python
class SWGASimulator:
    """Simulate complete SWGA reaction."""

    def simulate(
        self,
        primers: List[str],
        target_genome: str,
        background_genome: str,
        conditions: ReactionConditions
    ) -> Dict:
        """
        Simulate SWGA with target and background.

        Returns:
            {'target_yield': float, 'background_yield': float, 'enrichment': float}
        """
```

### stochastic_simulator.py

Gillespie algorithm simulation.

**Classes:**
```python
class GillespieSimulator:
    """Exact stochastic chemical kinetics."""

    def simulate(
        self,
        initial_state: Dict,
        reactions: List[Reaction],
        t_end: float
    ) -> Trajectory:
        """Run Gillespie simulation."""
```

### simulation_analysis.py

Analysis of simulation results.

**Functions:**
```python
def analyze_yield_curve(trajectory: Trajectory) -> Dict:
    """Analyze yield over time."""

def calculate_enrichment_factor(
    target_yield: float,
    background_yield: float,
    target_fraction: float
) -> float:
    """Calculate enrichment factor."""
```

### simulation_plots.py

Visualization of simulation results.

**Functions:**
```python
def plot_amplification_curve(result: SimulationResult, ax=None) -> plt.Axes
def plot_coverage_heatmap(coverage: np.ndarray, ax=None) -> plt.Axes
def plot_enrichment_comparison(results: List[SimulationResult]) -> plt.Figure
```

### simulation_report.py

HTML report generation.

**Functions:**
```python
def generate_simulation_report(
    results: SimulationResult,
    output_path: str,
    include_plots: bool = True
) -> str:
    """Generate HTML report with embedded plots."""
```

---

## Analysis Modules

### genome_analysis.py

Genome suitability analysis.

**Functions:**
```python
def analyze_genome(genome_path: str) -> Dict:
    """
    Analyze genome for SWGA suitability.

    Returns:
        Dict with gc_content, complexity, repeat_fraction, recommendation
    """

def identify_difficult_regions(genome: str) -> List[Tuple[int, int]]:
    """Find regions difficult to amplify."""
```

### genome_io.py

Genome file I/O utilities.

**Functions:**
```python
def read_fasta(path: str) -> Dict[str, str]:
    """Read FASTA file to dict of name -> sequence."""

def read_genbank(path: str) -> Dict:
    """Read GenBank file with annotations."""

def write_fasta(sequences: Dict[str, str], path: str) -> None:
    """Write sequences to FASTA file."""
```

### strand_bias_analyzer.py

Strand bias detection.

**Functions:**
```python
def analyze_strand_bias(
    primers: List[str],
    positions: Dict
) -> Dict:
    """
    Analyze forward/reverse strand balance.

    Returns:
        {'forward_fraction': float, 'reverse_fraction': float, 'bias_score': float}
    """
```

### dimer_network_analyzer.py

Primer dimer interaction network analysis.

**Classes:**
```python
class DimerNetworkAnalyzer:
    """Analyze primer-dimer interaction network."""

    def build_network(self, primers: List[str]) -> nx.Graph:
        """Build dimer interaction graph."""

    def find_problematic_clusters(self) -> List[Set[str]]:
        """Find groups of primers with many interactions."""

    def visualize(self, output_path: str) -> None:
        """Generate network visualization."""
```

### amplicon_network.py

Amplicon network analysis.

**Functions:**
```python
def analyze_amplicon_distribution(
    primers: List[str],
    positions: Dict,
    genome_length: int
) -> Dict:
    """
    Analyze expected amplicon size distribution.

    Returns:
        {'mean_size': float, 'std_size': float, 'histogram': np.ndarray}
    """
```

---

## Utility Modules

### utility.py

Multiprocessing and file utilities.

**Functions:**
```python
def create_pool(
    func: Callable,
    input_list: List,
    cpus: int,
    initializer: Callable = None,
    initargs: tuple = ()
) -> List:
    """Create multiprocessing pool and map function."""

def create_pool_with_progress(
    func: Callable,
    input_list: List,
    cpus: int,
    desc: str = "Processing"
) -> List:
    """Pool with tqdm progress bar."""

def flatten(l: List[List]) -> List:
    """Flatten nested list."""

def gini_exact(array: np.ndarray) -> float:
    """Calculate exact Gini coefficient."""

def get_positional_gap_lengths(
    positions: np.ndarray,
    circular: bool = True,
    seq_length: int = None
) -> np.ndarray:
    """Calculate gap lengths between positions."""
```

### parameter.py

Global parameter handling.

**Functions:**
```python
def get_params(json_path: str) -> dict:
    """Load parameters from JSON file."""

def get_current_config() -> PipelineParameters:
    """Get current configuration as dataclass."""

@dataclass
class PipelineParameters:
    min_k: int = 6
    max_k: int = 12
    min_fg_freq: float = 1e-5
    max_bg_freq: float = 5e-6
    # ... 50+ parameters
```

### validation.py

Installation validation.

**Functions:**
```python
def validate_installation(quick: bool = False) -> bool:
    """
    Validate NeoSWGA installation.

    Checks:
    - Python version
    - Required packages
    - Jellyfish availability
    - Model file integrity

    Args:
        quick: Skip slow validation tests

    Returns:
        True if all checks pass
    """
```

### exceptions.py

Custom exception hierarchy.

**Classes:**
```python
class NeoSWGAError(Exception):
    """Base exception for all NeoSWGA errors."""

class OptimizerNotFoundError(NeoSWGAError):
    """Raised when optimizer name is not registered."""

class ModelIntegrityError(NeoSWGAError):
    """Raised when model file hash doesn't match."""

class StepPrerequisiteError(NeoSWGAError):
    """Raised when pipeline step prerequisites are not met."""

class ConfigurationError(NeoSWGAError):
    """Raised for invalid configuration."""

class PositionFileNotFoundError(FileError):
    """Raised when HDF5 position file is missing."""

class GenomeFileError(FileError):
    """Raised for genome file I/O problems."""
```

### export.py

Primer export for synthesis ordering.

**Formats:** FASTA, CSV (vendor-ready for IDT/Twist/Sigma), BED, BedGraph, lab protocol.

**Classes:**
```python
class PrimerExporter:
    """High-level exporter supporting all output formats."""

class ModificationProfile(Enum):
    """Modification profiles: NONE, STANDARD, LOW_INPUT."""
```

**Functions:**
```python
def export_to_fasta(primers, output_path, ...) -> None
def export_to_bed(primers, output_path, genome_name, ...) -> None
def export_to_bedgraph(primers, output_path, ...) -> None
def export_to_vendor_csv(primers, output_path, vendor, ...) -> None
def generate_protocol(primers, conditions, ...) -> str
```

### progress.py

Progress bar utilities for long-running pipeline steps.

**Functions:**
```python
def progress_context(desc: str, disable: bool = False)
def progress_bar(iterable, **kwargs)
```

### search_context.py

Data structures for BFS-based optimization search.

**Classes:**
```python
class GenomeInfo: ...
class PositionData: ...
class DimerConstraints: ...
class SearchState: ...
class BFSConfig: ...
```

### pareto_frontier.py

Pareto frontier visualization and reporting.

**Functions:**
```python
def plot_frontier(results, output_path, ...) -> None
def generate_frontier_report(results, ...) -> str
def summarize_frontier_for_cli(results, ...) -> str
```

### primer_expansion.py

Expand primer sets to improve coverage in under-represented regions.

**Classes:**
```python
class PrimerExpander:
    """Identify coverage gaps and suggest additional primers."""
```

### efficiency_predictor.py

Predict amplification efficiency and enrichment fold-change.

**Classes:**
```python
class EfficiencyPredictor:
    """Predict enrichment based on primer set metrics and mechanistic model."""
```

### experimental_tracker.py

Track experimental outcomes for calibrating computational predictions.

**Classes:**
```python
class ExperimentalTracker:
    """Log experimental results and compute calibration reports."""
```

---

## Mechanistic Modeling Modules

### mechanistic_model.py

Four-pathway mechanistic model for SWGA amplification.

**Classes:**
```python
class MechanisticEffects:
    """Container for per-pathway factor values."""

class MechanisticModel:
    """Calculate amplification predictions from four pathways:
    Tm modification, secondary structure accessibility,
    enzyme activity, and binding kinetics."""
```

### mechanistic_params.py

Literature-derived parameters for the mechanistic model.

### additives.py

Additive effects on Tm with sigmoid GC normalization.

### additive_interactions.py

Registry of pairwise additive interactions (synergistic, antagonistic).

**Classes:**
```python
class AdditiveInteractionRegistry:
    """Registry of known additive-additive interactions."""
```

### additive_optimizer.py

Optimize additive cocktail for given genome and primer properties.

**Classes:**
```python
class AdditiveOptimizer:
    """Recommend optimal additive concentrations."""
```

### set_size_optimizer.py

Automatic primer set size recommendation by application profile.

**Functions:**
```python
def recommend_set_size(application, genome_length, ...) -> Dict
def quick_size_estimate(application, genome_length) -> int
```

### model_validation.py

Validation tests for the mechanistic model against expected literature behavior.

### background_registry.py

Registry of common background genomes with metadata.

**Classes:**
```python
class BackgroundRegistry:
    """Registry of pre-characterized background genomes (human, mouse, etc.)."""
```

---

## User Experience Modules

### wizard.py

Interactive setup wizard for guided params.json creation.

### param_validator.py

Parameter validation with error/warning/info levels.

### condition_suggester.py

Reaction condition recommendations based on genome properties.

### results_interpreter.py

Quality assessment and go/no-go recommendations for pipeline output.

### workflow_selector.py

Interactive menu for feature discovery.

---

## See Also

- [API Reference](API_REFERENCE.md) - Public API documentation
- [Architecture Diagrams](ARCHITECTURE_DIAGRAMS.md) - Visual architecture
- [User Guide](user-guide.md) - Usage tutorials
