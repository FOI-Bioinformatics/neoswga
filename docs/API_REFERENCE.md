# NeoSWGA API Reference

Complete API documentation for the NeoSWGA Python package. This reference covers all public interfaces, classes, and functions.

## Table of Contents

1. [CLI Commands](#cli-commands)
2. [Core Modules](#core-modules)
   - [thermodynamics](#thermodynamics)
   - [reaction_conditions](#reaction_conditions)
   - [position_cache](#position_cache)
   - [filter](#filter)
   - [genetic_algorithm](#genetic_algorithm)
   - [optimizer_factory](#optimizer_factory)
3. [Data Structures](#data-structures)
4. [Configuration](#configuration)

---

## CLI Commands

### Entry Point

```bash
neoswga <command> [options]
```

### Pipeline Commands

#### count-kmers

Generate k-mer counts for target and background genomes.

```bash
neoswga count-kmers -j params.json [options]
```

**Options:**
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-j, --json` | path | required | Parameter JSON file |
| `-x, --fasta-fore` | path | - | Target genome FASTA |
| `-y, --fasta-back` | path | - | Background genome FASTA |
| `--min-k` | int | 6 | Minimum k-mer length |
| `--max-k` | int | 12 | Maximum k-mer length |
| `--cpus` | int | 8 | Number of CPUs |

**Output:**
- `{prefix}_{k}mer_all.txt`: Tab-separated k-mer counts
- `{prefix}_{k}mer_positions.h5`: HDF5 file with binding positions

---

#### filter

Filter candidate primers based on sequence properties.

```bash
neoswga filter -j params.json [options]
```

**Options:**
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-j, --json` | path | required | Parameter JSON file |
| `--min-fg-freq` | float | 1e-5 | Minimum foreground frequency |
| `--max-bg-freq` | float | 5e-6 | Maximum background frequency |
| `--max-gini` | float | 0.6 | Maximum Gini index |
| `--max-primer` | int | 500 | Maximum primers to keep |
| `--use-background-filter` | flag | - | Use Bloom filter for background |

**Filtering Rules:**
1. No 3+ consecutive G/C at 3' end
2. GC content 40-60%
3. GC clamp: max 3 G/C in last 5 bases
4. Max 4 consecutive di-nucleotide repeats
5. Max 4 bp homopolymer runs

**Output:**
- `step2_df.csv`: Filtered primers with properties

---

#### score

Score primers for amplification efficacy using machine learning.

```bash
neoswga score -j params.json [options]
```

**Options:**
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-j, --json` | path | required | Parameter JSON file |
| `--cpus` | int | 8 | Number of CPUs |

**Output:**
- `step3_df.csv`: Primers with `amp_pred` score column

---

#### optimize

Find optimal primer sets using specified optimization method.

```bash
neoswga optimize -j params.json [options]
```

**Options:**
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-j, --json` | path | required | Parameter JSON file |
| `--optimization-method` | str | hybrid | Optimization algorithm |
| `--num-primers` | int | 6 | Target set size |
| `--iterations` | int | 8 | Search iterations |
| `--max-sets` | int | 5 | Parallel sets to build |
| `--no-background` | flag | false | Host-free optimization mode (no background genome) |
| `--use-mechanistic-model` | flag | false | Use mechanistic model for primer weighting |
| `--mechanistic-weight` | float | 0.3 | Weight for mechanistic model scoring |
| `--auto-size` | flag | false | Auto-size primer set based on application |
| `--application` | str | enrichment | Application profile (discovery, clinical, enrichment, metagenomics) |
| `--validate-with-simulation` | flag | false | Post-hoc simulation validation of results |

**Optimization Methods:**
| Method | Speed | Use Case |
|--------|-------|----------|
| `hybrid` | Medium | General use (default) |
| `greedy` | Fast | Simple optimization |
| `network` | Medium | Tm-weighted selection |
| `dominating-set` | Fast | Large primer pools, set-cover |
| `weighted-set-cover` | Fast | Score-weighted set cover |
| `background-aware` | Slow | Clinical applications |
| `genetic` | Moderate | Evolutionary multi-criteria |
| `moea` | Slow | Pareto optimization (requires pymoo) |
| `milp` | Variable | Exact solutions (requires mip) |
| `equiphi29` | Medium | EquiPhi29-specific at 42-45C |
| `tiling` | Fast | Interval-based coverage tiling |
| `normalized` | Medium | Strategy-preset scoring |
| `clique` | Moderate | Dimer-free set via max-clique |
| `multi-agent` | Slow | Parallel multi-optimizer ensemble |
| `bg-prefilter` | Medium | Background pruning + inner optimizer |
| `coverage-then-dimerfree` | Medium | Dominating-set then clique cascade |
| `dimerfree-scored` | Medium | Clique then network scoring cascade |
| `bg-prefilter-hybrid` | Medium | Background pre-filter then hybrid |

**Output:**
- `step4_improved_df.csv`: Optimized primer sets

---

### Utility Commands

#### build-filter

Pre-build Bloom filter for large background genomes.

```bash
neoswga build-filter <genome.fasta> <output_dir>
```

**Arguments:**
- `genome`: Path to background genome FASTA
- `output_dir`: Directory for filter files

**Output:**
- `background_bloom.pkl`: Serialized Bloom filter

---

#### validate

Validate installation and dependencies.

```bash
neoswga validate [--quick]
```

**Options:**
| Option | Description |
|--------|-------------|
| `--quick` | Run minimal validation |

---

#### show-presets

Display available reaction condition presets.

```bash
neoswga show-presets
```

**Available Presets:**
- `standard_phi29`: 30C, standard conditions
- `enhanced_equiphi29`: 42C, 5% DMSO, 1M betaine
- `long_primers_15mer`: 45C, 7% DMSO, 1.5M betaine
- `high_gc_genome`: 45C, 10% DMSO, 2M betaine

---

### Setup and Analysis Commands

#### init

Interactive setup wizard for creating `params.json`.

```bash
neoswga init --genome target.fna [--background host.fna] [-o params.json]
```

---

#### validate-params

Validate a `params.json` file before running the pipeline.

```bash
neoswga validate-params -j params.json
```

Reports errors, warnings, and informational messages about parameter values.

---

#### suggest

Suggest reaction conditions based on target genome properties.

```bash
neoswga suggest --genome-gc 0.5 [--primer-length 12] [--polymerase phi29]
neoswga suggest --genome target.fna   # Auto-calculates GC content
```

**Options:**
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--genome-gc` | float | required | Target genome GC content (0-1) |
| `--genome` | path | - | Calculate GC from FASTA (alternative to --genome-gc) |
| `--primer-length` | int | 10 | Expected primer length |
| `--polymerase` | str | phi29 | Polymerase (phi29, equiphi29, bst, klenow) |
| `--sweep` | flag | false | Run condition grid search (108 combinations) |
| `--output` | path | - | Write sweep results to CSV |

---

#### interpret

Interpret pipeline results with quality assessment and enrichment prediction.

```bash
neoswga interpret -d results/
```

Output includes primer count, quality ratings, recommendations, and (when `params.json` is present) predicted enrichment fold-change based on the mechanistic model.

---

#### report

Generate quality reports.

```bash
neoswga report -d results/                        # Executive summary
neoswga report -d results/ --level full            # Full technical report
neoswga report -d results/ --interactive           # With Plotly charts
neoswga report -d results/ --check                 # Validate only
```

---

#### export

Export primer results in various formats.

```bash
neoswga export -d results/ [--format FORMAT] [-o output_path]
```

**Options:**
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-d, --data-dir` | path | required | Results directory |
| `--format` | str | fasta | Output format |
| `-o, --output` | path | - | Output file path |
| `--genome-name` | str | genome | Chromosome/contig name for BED/BedGraph |
| `--window-size` | int | 1000 | Window size for BedGraph coverage |

**Export Formats:**
| Format | Description |
|--------|-------------|
| `fasta` | FASTA sequences (default) |
| `csv` | CSV with primer attributes |
| `bed` | 6-column BED for genome browsers |
| `bedgraph` | Windowed coverage density |
| `protocol` | Lab protocol with ordering info |
| `all` | All formats at once |

---

### Advanced Commands

#### analyze-set

Analyze an existing primer set.

```bash
neoswga analyze-set --primers SEQ1 SEQ2 ... --fg <genome> --fg-kmers <prefix> --output <dir>
```

#### analyze-genome

Analyze genome suitability for SWGA.

```bash
neoswga analyze-genome --genome <genome.fasta> --output <dir>
```

#### analyze-dimers

Analyze primer-dimer interactions.

```bash
neoswga analyze-dimers --primers SEQ1 SEQ2 ... --output <dir> [--visualize]
```

#### simulate

Run replication simulation.

```bash
neoswga simulate --primers SEQ1 SEQ2 --genome <genome.fasta> --output <dir>
```

#### multi-genome

Pan-genome primer design across multiple targets.

```bash
neoswga multi-genome --genomes target1.fna target2.fna --output <dir>
```

---

## Core Modules

### thermodynamics

DNA hybridization thermodynamics using the SantaLucia nearest-neighbor model.

```python
from neoswga.core import thermodynamics as thermo
```

#### Functions

##### calculate_tm_with_salt

Calculate melting temperature with salt correction.

```python
def calculate_tm_with_salt(
    seq: str,
    na_conc: float = 50.0,
    mg_conc: float = 0.0,
    dmso_percent: float = 0.0,
    betaine_m: float = 0.0,
    trehalose_m: float = 0.0
) -> float:
    """
    Calculate Tm using SantaLucia nearest-neighbor model with Owczarzy salt correction.

    Args:
        seq: DNA sequence (5' to 3')
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM
        dmso_percent: DMSO percentage (0-10)
        betaine_m: Betaine concentration in M
        trehalose_m: Trehalose concentration in M

    Returns:
        Melting temperature in Celsius

    Example:
        >>> thermo.calculate_tm_with_salt('ATCGATCGATCG', na_conc=50)
        32.5
    """
```

##### calculate_free_energy

Calculate Gibbs free energy at specified temperature.

```python
def calculate_free_energy(
    seq: str,
    temp: float = 37.0,
    na_conc: float = 50.0,
    mg_conc: float = 0.0
) -> float:
    """
    Calculate DeltaG using NN model.

    Args:
        seq: DNA sequence
        temp: Temperature in Celsius
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM

    Returns:
        Free energy in kcal/mol (negative = favorable)

    Example:
        >>> thermo.calculate_free_energy('ATCGATCG', temp=30)
        -8.2
    """
```

##### gc_content

Calculate GC content fraction.

```python
def gc_content(seq: str) -> float:
    """
    Calculate GC content as fraction.

    Args:
        seq: DNA sequence

    Returns:
        GC fraction (0.0 to 1.0)

    Example:
        >>> thermo.gc_content('ATCGATCG')
        0.5
    """
```

##### reverse_complement

Return reverse complement of DNA sequence.

```python
def reverse_complement(seq: str) -> str:
    """
    Return reverse complement.

    Args:
        seq: DNA sequence

    Returns:
        Reverse complement sequence

    Example:
        >>> thermo.reverse_complement('ATCG')
        'CGAT'
    """
```

##### is_palindrome

Check if sequence is palindromic.

```python
def is_palindrome(seq: str) -> bool:
    """
    Check if sequence equals its reverse complement.

    Args:
        seq: DNA sequence

    Returns:
        True if palindromic

    Example:
        >>> thermo.is_palindrome('GAATTC')
        True
    """
```

#### Constants

```python
# Universal gas constant
R = 1.987  # cal/(mol*K)

# SantaLucia nearest-neighbor parameters
ENTHALPY_NN = {
    'AA/TT': -7.9,
    'AT/TA': -7.2,
    'TA/AT': -7.2,
    'CA/GT': -8.5,
    'GT/CA': -8.4,
    'CT/GA': -7.8,
    'GA/CT': -8.2,
    'CG/GC': -10.6,
    'GC/CG': -9.8,
    'GG/CC': -8.0,
    # ... reverse complements
}

ENTROPY_NN = {
    'AA/TT': -22.2,
    'AT/TA': -20.4,
    # ... all 16 dinucleotide stacks
}

# Helix initiation
INIT_ENTHALPY = 0.2   # kcal/mol
INIT_ENTROPY = -5.7   # cal/(mol*K)

# Terminal AT penalty
TERMINAL_AT_ENTHALPY = 2.2  # kcal/mol
TERMINAL_AT_ENTROPY = 6.9   # cal/(mol*K)

# Symmetry correction
SYMMETRY_ENTROPY = -1.4  # cal/(mol*K)
```

---

### reaction_conditions

Model reaction conditions for SWGA experiments.

```python
from neoswga.core import reaction_conditions as rc
```

#### Classes

##### ReactionConditions

Container for all reaction parameters.

```python
class ReactionConditions:
    """
    Comprehensive reaction condition model.

    Attributes:
        temp: Reaction temperature (Celsius)
        dmso_percent: DMSO concentration (0-10%)
        betaine_m: Betaine concentration (0-2.5 M)
        trehalose_m: Trehalose concentration (0-1.0 M)
        formamide_percent: Formamide concentration (0-10%)
        glycerol_percent: Glycerol concentration (0-15%)
        bsa_ug_ml: BSA concentration (0-400 ug/mL)
        peg_percent: PEG concentration (0-15%)
        ethanol_percent: Ethanol concentration (0-5%)
        urea_m: Urea concentration (0-2.0 M)
        tmac_m: TMAC concentration (0-0.1 M)
        na_conc: Na+ concentration (mM)
        mg_conc: Mg2+ concentration (mM)
        polymerase: Polymerase name ('phi29', 'equiphi29', 'bst', 'klenow')
        ssb: Enable SSB proteins (bool)
    """

    def calculate_effective_tm(self, seq: str) -> float:
        """Calculate Tm with all additive corrections."""

    def max_primer_length(self) -> int:
        """Get maximum recommended primer length for conditions."""

    def is_primer_compatible(self, seq: str) -> bool:
        """Check if primer is compatible with conditions."""
```

#### Functions

##### get_standard_conditions

Get standard Phi29 conditions.

```python
def get_standard_conditions() -> ReactionConditions:
    """
    Standard Phi29 conditions at 30C.

    Returns:
        ReactionConditions with default Phi29 parameters

    Example:
        >>> conditions = rc.get_standard_conditions()
        >>> conditions.temp
        30.0
        >>> conditions.max_primer_length()
        12
    """
```

##### get_enhanced_conditions

Get enhanced EquiPhi29 conditions.

```python
def get_enhanced_conditions() -> ReactionConditions:
    """
    Enhanced EquiPhi29 conditions with additives.

    Enables longer primers (up to 15bp) through:
    - Temperature: 42C
    - DMSO: 5%
    - Betaine: 1M

    Returns:
        ReactionConditions for EquiPhi29

    Example:
        >>> conditions = rc.get_enhanced_conditions()
        >>> conditions.max_primer_length()
        15
    """
```

##### get_polymerase_processivity

Get processivity for a polymerase.

```python
def get_polymerase_processivity(polymerase: str) -> int:
    """
    Get maximum extension length in base pairs.

    Args:
        polymerase: 'phi29', 'equiphi29', 'bst', or 'klenow'

    Returns:
        Processivity in bp

    Example:
        >>> rc.get_polymerase_processivity('phi29')
        70000
    """
```

##### list_polymerases

List all supported polymerases.

```python
def list_polymerases() -> Dict[str, str]:
    """
    List polymerases with descriptions.

    Returns:
        Dict mapping name to description
    """
```

#### Polymerase Database

| Polymerase | Temp Range | Optimal | Processivity | Use Case |
|------------|------------|---------|--------------|----------|
| `phi29` | 30-40C | 30C | 70 kb | Standard SWGA |
| `equiphi29` | 42-45C | 42C | 80 kb | Higher specificity |
| `bst` | 60-65C | 63C | 1-2 kb | LAMP-compatible |
| `klenow` | 25-40C | 37C | 10 kb | Budget alternative |

#### Additive Effects

| Additive | Effect on Tm | Typical Range |
|----------|-------------|---------------|
| DMSO | -0.6C per % | 0-10% |
| Betaine | -2.3C per M | 0-2.5 M |
| Trehalose | -5C per M | 0-1.0 M |
| Formamide | -0.72C per % | 0-10% |
| Ethanol | -0.5C per % | 0-5% |
| Urea | -2.0C per M | 0-2.0 M |
| TMAC | Equalizes AT/GC | 0-0.1 M |

---

### position_cache

In-memory cache for primer binding positions with O(1) lookup.

```python
from neoswga.core.position_cache import PositionCache
```

#### Classes

##### PositionCache

High-performance position lookup cache.

```python
class PositionCache:
    """
    In-memory cache of all primer binding positions.

    Provides 1000x speedup over repeated HDF5 reads.
    Memory usage: ~4 MB for 500 primers x 1000 sites.

    Example:
        cache = PositionCache(fg_prefixes, candidate_primers)
        positions = cache.get_positions('fg_prefix', 'ATCGATCG')
    """

    def __init__(self, fname_prefixes: List[str], primers: List[str]):
        """
        Load all positions for given primers from HDF5 files.

        Args:
            fname_prefixes: List of path prefixes (e.g., ['data/ecoli'])
            primers: List of primer sequences to cache

        Raises:
            FileNotFoundError: If HDF5 files don't exist
        """

    def get_positions(
        self,
        fname_prefix: str,
        primer: str,
        strand: str = 'both'
    ) -> np.ndarray:
        """
        Get positions for a primer. O(1) lookup.

        Args:
            fname_prefix: Genome identifier
            primer: Primer sequence
            strand: 'forward', 'reverse', or 'both'

        Returns:
            Array of positions (empty if not found)

        Note:
            Strand must be 'forward', 'reverse', or 'both'.
            Do NOT use '+' or '-'.
        """

    def get_all_positions(
        self,
        fname_prefix: str,
        primers: List[str]
    ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Get positions for multiple primers (batched).

        Returns:
            Dict mapping primer -> (forward_positions, reverse_positions)
        """

    def compute_coverage_vectorized(
        self,
        fname_prefix: str,
        primers: List[str],
        genome_length: int
    ) -> np.ndarray:
        """
        Compute coverage array for primer set (vectorized).

        Args:
            fname_prefix: Genome identifier
            primers: List of primers
            genome_length: Length of genome

        Returns:
            Boolean array of length genome_length (True = covered)
        """
```

##### BindingSite

Single primer binding site dataclass.

```python
@dataclass
class BindingSite:
    """Single primer binding site."""
    position: int
    strand: str  # 'forward' or 'reverse'
    primer: str

    def __hash__(self):
        return hash((self.position, self.strand, self.primer))
```

---

### filter

Primer filtering based on sequence properties.

```python
from neoswga.core import filter
```

#### Functions

##### filter_extra

Apply all five filtering rules.

```python
def filter_extra(primer: str) -> bool:
    """
    Filter primer based on sequence rules.

    Rules:
    1. No 3+ consecutive G/C at 3' end
    2. GC content 40-60%
    3. GC clamp: max 3 G/C in last 5 bases
    4. Max 4 consecutive di-nucleotide repeats
    5. Max 4 bp homopolymer runs

    Args:
        primer: Primer sequence (5' to 3')

    Returns:
        True if passes all filters

    Example:
        >>> filter.filter_extra('ATCGATCGATCG')
        True
        >>> filter.filter_extra('GGGGGGATCG')
        False  # homopolymer run
    """
```

##### get_bg_rates_via_bloom

Check background presence using Bloom filter.

```python
def get_bg_rates_via_bloom(
    primer_list: List[str],
    bloom_path: str
) -> Dict[str, int]:
    """
    Get background rates using pre-built Bloom filter.

    Memory-efficient for large genomes (human 3 Gbp).
    O(1) lookup per primer.

    Args:
        primer_list: Primers to check
        bloom_path: Path to Bloom filter file

    Returns:
        Dict mapping primer -> estimated count (0 = not in background)
    """
```

---

### genetic_algorithm

Evolutionary optimization for primer set selection.

```python
from neoswga.core.genetic_algorithm import PrimerSetGA, GAConfig, Individual
```

#### Classes

##### GAConfig

Configuration dataclass.

```python
@dataclass
class GAConfig:
    """Configuration for genetic algorithm."""
    population_size: int = 200
    generations: int = 100
    mutation_rate: float = 0.15
    crossover_rate: float = 0.8
    elitism_fraction: float = 0.10
    tournament_size: int = 5
    min_set_size: int = 4
    max_set_size: int = 8
    max_dimer_severity: float = 0.5
    n_processes: int = None  # None = use all CPUs
```

##### Individual

Represents a primer set in population.

```python
@dataclass
class Individual:
    """Represents a primer set (individual in population)."""
    primers: List[str]
    fitness: float = None
    metrics: Dict = None
```

##### PrimerSetGA

Genetic algorithm optimizer.

```python
class PrimerSetGA:
    """
    Genetic Algorithm for primer set optimization.

    Evolutionary operators:
    - Selection: Tournament selection
    - Crossover: Uniform crossover with dimer checking
    - Mutation: Add/remove/replace primers
    - Elitism: Preserve top performers
    """

    def __init__(
        self,
        primer_pool: List[str],
        fg_prefixes: List[str],
        bg_prefixes: List[str],
        fg_lengths: List[int],
        bg_lengths: List[int],
        conditions: rc.ReactionConditions,
        config: Optional[GAConfig] = None,
        position_cache=None
    ):
        """
        Initialize genetic algorithm.

        Args:
            primer_pool: Available primers to choose from
            fg_prefixes: Foreground HDF5 file prefixes
            bg_prefixes: Background HDF5 file prefixes
            fg_lengths: Foreground genome lengths
            bg_lengths: Background genome lengths
            conditions: Reaction conditions
            config: GA configuration
            position_cache: Optional PositionCache
        """

    def evolve(self, verbose: bool = True) -> Individual:
        """
        Run genetic algorithm evolution.

        Args:
            verbose: Print progress

        Returns:
            Best individual found

        Example:
            ga = PrimerSetGA(primers, fg, bg, fg_lens, bg_lens, conditions)
            best = ga.evolve()
            print(f"Best primers: {best.primers}")
            print(f"Fitness: {best.fitness}")
        """
```

---

### optimizer_factory

Factory for creating optimizer instances.

```python
from neoswga.core.optimizer_factory import OptimizerFactory, OptimizerRegistry
```

#### Classes

##### OptimizerFactory

Factory for optimizer instantiation.

```python
class OptimizerFactory:
    """
    Factory for creating optimizer instances.

    Example:
        optimizer = OptimizerFactory.create('greedy', cache=cache, ...)
    """

    @staticmethod
    def create(
        name: str,
        cache: PositionCache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        **kwargs
    ) -> BaseOptimizer:
        """
        Create optimizer by name.

        Args:
            name: Optimizer name (e.g., 'greedy', 'network', 'genetic')
            cache: Position cache
            fg_prefixes: Foreground file prefixes
            fg_seq_lengths: Foreground genome lengths
            **kwargs: Additional optimizer-specific arguments

        Returns:
            Configured optimizer instance

        Raises:
            OptimizerNotFoundError: If name not registered
        """

    @staticmethod
    def list_optimizers() -> Dict[str, str]:
        """
        List available optimizers with descriptions.

        Returns:
            Dict mapping name to description
        """
```

##### OptimizerRegistry

Registry for optimizer classes.

```python
class OptimizerRegistry:
    """
    Thread-safe singleton registry for optimizer classes.

    Example:
        @OptimizerRegistry.register('my-optimizer', aliases=['my', 'mo'])
        class MyOptimizer(BaseOptimizer):
            ...
    """

    @classmethod
    def register(
        cls,
        name: str,
        aliases: Optional[List[str]] = None,
        description: Optional[str] = None
    ) -> Callable:
        """
        Decorator to register an optimizer class.

        Args:
            name: Primary name (kebab-case preferred)
            aliases: Alternative names
            description: Optional description

        Returns:
            Decorator function
        """
```

---

## Data Structures

### File Formats

#### K-mer Counts (`*_Xmer_all.txt`)

Tab-separated text file:

```
ATCGATCG    1523    45
GCTAGCTA    892     12
...
```

Columns: primer | foreground_count | background_count

#### Filtered Primers (`step2_df.csv`)

CSV with columns:

| Column | Type | Description |
|--------|------|-------------|
| sequence | str | Primer sequence |
| fg_freq | float | Foreground frequency |
| bg_freq | float | Background frequency |
| gini | float | Gini index (binding evenness) |
| tm | float | Melting temperature |
| gc | float | GC content |
| self_dimer_bp | int | Self-dimer base pairs |
| het_dimer_bp | int | Hetero-dimer base pairs |
| complexity | float | Sequence complexity |

#### Scored Primers (`step3_df.csv`)

Adds column:

| Column | Type | Description |
|--------|------|-------------|
| amp_pred | float | Amplification prediction score (0-1) |

#### Optimized Sets (`step4_improved_df.csv`)

Final output with set metrics:

| Column | Type | Description |
|--------|------|-------------|
| set_id | int | Primer set identifier |
| primers | str | Comma-separated primer list |
| coverage | float | Target genome coverage |
| enrichment | float | Enrichment factor |
| mean_tm | float | Mean Tm of set |
| max_dimer_score | float | Maximum dimer severity |

#### Position Data (`*_positions.h5`)

HDF5 file structure:

```
/ATCGATCG    -> int32[] (binding positions)
/GCTAGCTA    -> int32[] (binding positions)
...
```

---

## Configuration

### Parameter JSON Schema

```json
{
  "fg_genomes": ["path/to/target.fasta"],
  "bg_genomes": ["path/to/background.fasta"],
  "fg_prefixes": ["path/to/target_kmers"],
  "bg_prefixes": ["path/to/background_kmers"],
  "data_dir": "results/",

  "min_k": 6,
  "max_k": 12,
  "min_fg_freq": 1e-5,
  "max_bg_freq": 5e-6,
  "max_gini": 0.6,
  "max_primer": 500,

  "polymerase": "phi29",
  "reaction_temp": 30.0,
  "na_conc": 50.0,
  "mg_conc": 0.0,
  "dmso_percent": 0.0,
  "betaine_m": 0.0,
  "min_tm": 15.0,
  "max_tm": 45.0,

  "optimization_method": "hybrid",
  "num_primers": 6,
  "iterations": 8,
  "max_sets": 5,

  "cpus": 8,
  "verbose": false
}
```

### Environment Variables

| Variable | Description |
|----------|-------------|
| `NEOSWGA_DATA_DIR` | Default data directory |
| `NEOSWGA_CPUS` | Default CPU count |
| `CUDA_VISIBLE_DEVICES` | GPU device selection |

---

## See Also

- [User Guide](user-guide.md) - Usage tutorials
- [Architecture](development/architecture.md) - Design documentation
- [Algorithms](development/algorithms.md) - Algorithm details
