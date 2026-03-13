# NeoSWGA Architecture Diagrams

Visual documentation of the NeoSWGA system architecture, data flows, and component relationships.

## Table of Contents

1. [System Overview](#system-overview)
2. [Pipeline Data Flow](#pipeline-data-flow)
3. [Module Dependencies](#module-dependencies)
4. [Optimizer Architecture](#optimizer-architecture)
5. [Thermodynamic Calculations](#thermodynamic-calculations)
6. [Performance Optimization](#performance-optimization)

---

## System Overview

High-level architecture showing major components and their interactions.

```mermaid
graph TB
    subgraph "User Interface"
        CLI[CLI - cli_unified.py]
        JSON[params.json]
    end

    subgraph "Pipeline Steps"
        S1[Step 1: count-kmers]
        S2[Step 2: filter]
        S3[Step 3: score]
        S4[Step 4: optimize]
    end

    subgraph "Core Modules"
        THERMO[Thermodynamics]
        RC[Reaction Conditions]
        CACHE[Position Cache]
        FILTER[Filter Engine]
        RF[Random Forest Model]
        OPT[Optimizer Factory]
    end

    subgraph "Data Layer"
        HDF5[(HDF5 Positions)]
        CSV[(CSV Results)]
        BLOOM[(Bloom Filter)]
    end

    subgraph "External"
        JF[Jellyfish]
        GPU[GPU - CuPy]
    end

    CLI --> S1
    CLI --> S2
    CLI --> S3
    CLI --> S4
    JSON --> CLI

    S1 --> JF
    S1 --> HDF5

    S2 --> FILTER
    S2 --> THERMO
    S2 --> BLOOM
    FILTER --> CSV

    S3 --> RF
    S3 --> CACHE
    RF --> CSV

    S4 --> OPT
    S4 --> CACHE
    OPT --> RC
    OPT --> THERMO
    OPT --> CSV

    CACHE --> HDF5
    THERMO --> GPU
```

---

## Pipeline Data Flow

Detailed view of data transformations through the pipeline.

```mermaid
flowchart LR
    subgraph Input
        FG[Target Genome<br/>FASTA]
        BG[Background Genome<br/>FASTA]
    end

    subgraph "Step 1: count-kmers"
        JF1[Jellyfish<br/>k-mer counting]
        KMER[K-mer Files<br/>*_Xmer_all.txt]
        POS[Position Files<br/>*_positions.h5]
    end

    subgraph "Step 2: filter"
        FREQ[Frequency Filter<br/>fg_freq, bg_freq]
        GINI[Gini Filter<br/>binding evenness]
        TM[Tm Filter<br/>melting temperature]
        SEQ[Sequence Filter<br/>5 rules]
        S2OUT[step2_df.csv<br/>500 primers]
    end

    subgraph "Step 3: score"
        FEAT[Feature Engineering<br/>52-120 features]
        RFMOD[Random Forest<br/>Classifier]
        S3OUT[step3_df.csv<br/>+ amp_pred]
    end

    subgraph "Step 4: optimize"
        OPTIM[Optimizer<br/>hybrid/GA/MILP]
        DIMER[Dimer Check<br/>secondary structure]
        S4OUT[step4_improved_df.csv<br/>6 primers]
    end

    FG --> JF1
    BG --> JF1
    JF1 --> KMER
    JF1 --> POS

    KMER --> FREQ
    POS --> GINI
    FREQ --> TM
    GINI --> TM
    TM --> SEQ
    SEQ --> S2OUT

    S2OUT --> FEAT
    POS --> FEAT
    FEAT --> RFMOD
    RFMOD --> S3OUT

    S3OUT --> OPTIM
    POS --> OPTIM
    OPTIM --> DIMER
    DIMER --> S4OUT
```

---

## Module Dependencies

Import relationships between core modules.

```mermaid
graph TD
    subgraph "Entry Point"
        CLI[cli_unified.py]
    end

    subgraph "Pipeline"
        PIPE[pipeline.py]
        IMPIPE[improved_pipeline.py]
        AUTOPIPE[auto_swga_pipeline.py]
        MULTIPIPE[multi_genome_pipeline.py]
    end

    subgraph "Filtering"
        FILT[filter.py]
        ADAPT[adaptive_filters.py]
        KMER[kmer.py]
        TFILT[thermodynamic_filter.py]
    end

    subgraph "Scoring"
        RF[rf_preprocessing.py]
        PATTR[primer_attributes.py]
        IQS[integrated_quality_scorer.py]
    end

    subgraph "Optimization"
        OPTF[optimizer_factory.py]
        BASE[base_optimizer.py]
        GREEDY[greedy_optimizer.py]
        NETWORK[network_optimizer.py]
        GA[genetic_algorithm.py]
        MILP[milp_optimizer.py]
        DS[dominating_set_optimizer.py]
        HYBRID[hybrid_optimizer.py]
        TILE[tiling_optimizer.py]
        CLIQUE[clique_optimizer.py]
        NORM[normalized_optimizer.py]
        EQOPT[equiphi29_optimizer.py]
        CASCADE[serial_cascade_optimizer.py]
        MAGENT[multi_agent_optimizer.py]
    end

    subgraph "Thermodynamics"
        THERMO[thermodynamics.py]
        RC[reaction_conditions.py]
        SS[secondary_structure.py]
        DIMER[dimer.py]
    end

    subgraph "Performance"
        PCACHE[position_cache.py]
        BGFILT[background_filter.py]
        GPU[gpu_acceleration.py]
    end

    subgraph "Utilities"
        UTIL[utility.py]
        PARAM[parameter.py]
        EXCEPT[exceptions.py]
    end

    CLI --> PIPE
    CLI --> IMPIPE

    PIPE --> FILT
    PIPE --> RF
    PIPE --> OPTF

    FILT --> THERMO
    FILT --> DIMER
    FILT --> PARAM

    RF --> THERMO
    RF --> PATTR

    OPTF --> BASE
    OPTF --> GREEDY
    OPTF --> NETWORK
    OPTF --> GA
    OPTF --> HYBRID

    GREEDY --> BASE
    NETWORK --> BASE
    GA --> BASE
    MILP --> BASE
    DS --> BASE
    HYBRID --> BASE
    TILE --> BASE
    CLIQUE --> BASE
    NORM --> BASE
    EQOPT --> BASE
    CASCADE --> BASE
    MAGENT --> BASE

    BASE --> PCACHE
    BASE --> THERMO
    BASE --> RC

    THERMO --> UTIL
    RC --> THERMO

    SS --> THERMO
    DIMER --> THERMO

    PCACHE --> UTIL
    GPU --> THERMO
```

---

## Optimizer Architecture

Factory pattern and optimizer inheritance hierarchy.

```mermaid
classDiagram
    class BaseOptimizer {
        <<abstract>>
        +cache: PositionCache
        +fg_prefixes: List[str]
        +fg_seq_lengths: List[int]
        +optimize(candidates, target_size) OptimizationResult
        #_calculate_metrics(primers) PrimerSetMetrics
    }

    class OptimizerFactory {
        +create(name, cache, ...) BaseOptimizer
        +list_optimizers() Dict[str, str]
    }

    class OptimizerRegistry {
        -_registry: Dict[str, Type]
        -_aliases: Dict[str, str]
        +register(name, aliases) Decorator
        +get(name) Type[BaseOptimizer]
    }

    class GreedyOptimizer {
        +optimize(candidates, target_size) OptimizationResult
        -_greedy_select(candidates) List[str]
    }

    class NetworkOptimizer {
        +network: AmplificationNetwork
        +optimize(candidates, target_size) OptimizationResult
        -_build_network(primers) None
        -_score_connectivity() float
    }

    class GeneticOptimizer {
        +config: GAConfig
        +population: List[Individual]
        +optimize(candidates, target_size) OptimizationResult
        -_evolve() Individual
        -_select_parents() List[Individual]
        -_crossover(p1, p2) Individual
        -_mutate(ind) Individual
    }

    class HybridOptimizer {
        +network_opt: NetworkOptimizer
        +greedy_opt: GreedyOptimizer
        +optimize(candidates, target_size) OptimizationResult
    }

    class DominatingSetOptimizer {
        +bin_size: int
        +optimize(candidates, target_size) OptimizationResult
        -_build_coverage_graph() Graph
        -_greedy_set_cover() Set[str]
    }

    class MILPOptimizer {
        +solver: mip.Model
        +optimize(candidates, target_size) OptimizationResult
        -_formulate_ilp() None
        -_solve() List[str]
    }

    class BackgroundAwareOptimizer {
        +stages: int
        +optimize(candidates, target_size) OptimizationResult
        -_stage1_coverage() List[str]
        -_stage2_background() List[str]
        -_stage3_dimers() List[str]
    }

    class OptimizationResult {
        +primer_set: List[str]
        +metrics: PrimerSetMetrics
        +status: OptimizationStatus
        +message: str
    }

    class PrimerSetMetrics {
        +fg_coverage: float
        +bg_coverage: float
        +selectivity_ratio: float
        +mean_tm: float
        +dimer_risk_score: float
        +gap_gini: float
        +normalized_score() float
    }

    class CompositeOptimizer {
        +optimizers: List[BaseOptimizer]
        +optimize(candidates, target_size) OptimizationResult
    }

    class TilingOptimizer {
        +optimize(candidates, target_size) OptimizationResult
    }

    class NormalizedOptimizer {
        +strategy: str
        +optimize(candidates, target_size) OptimizationResult
    }

    class CliqueOptimizer {
        +optimize(candidates, target_size) OptimizationResult
    }

    class EquiPhi29Optimizer {
        +optimize(candidates, target_size) OptimizationResult
    }

    class MultiAgentOrchestrator {
        +strategies: List
        +optimize(candidates, target_size) OptimizationResult
    }

    class SerialCascadeOptimizer {
        +stages: List[BaseOptimizer]
        +optimize(candidates, target_size) OptimizationResult
    }

    OptimizerFactory --> OptimizerRegistry
    OptimizerRegistry --> BaseOptimizer
    BaseOptimizer <|-- GreedyOptimizer
    BaseOptimizer <|-- NetworkOptimizer
    BaseOptimizer <|-- GeneticOptimizer
    BaseOptimizer <|-- HybridOptimizer
    BaseOptimizer <|-- DominatingSetOptimizer
    BaseOptimizer <|-- MILPOptimizer
    BaseOptimizer <|-- BackgroundAwareOptimizer
    BaseOptimizer <|-- TilingOptimizer
    BaseOptimizer <|-- NormalizedOptimizer
    BaseOptimizer <|-- CliqueOptimizer
    BaseOptimizer <|-- EquiPhi29Optimizer
    BaseOptimizer <|-- MultiAgentOrchestrator
    BaseOptimizer <|-- SerialCascadeOptimizer
    BaseOptimizer <|-- CompositeOptimizer

    HybridOptimizer --> NetworkOptimizer
    HybridOptimizer --> GreedyOptimizer

    BaseOptimizer --> OptimizationResult
    OptimizationResult --> PrimerSetMetrics
```

---

## Thermodynamic Calculations

Nearest-neighbor model calculation flow.

```mermaid
flowchart TB
    subgraph Input
        SEQ[DNA Sequence<br/>5'-ATCGATCG-3']
        SALT[Salt Concentrations<br/>Na+, Mg2+]
        ADD[Additives<br/>DMSO, Betaine]
        TEMP[Temperature]
    end

    subgraph "Nearest-Neighbor Model"
        NN[Dinucleotide Stacks<br/>16 parameters]
        INIT[Helix Initiation<br/>DH=0.2, DS=-5.7]
        TERM[Terminal AT Penalty<br/>DH=2.2, DS=6.9]
        SYM[Symmetry Correction<br/>DS=-1.4]
    end

    subgraph "Calculations"
        DH[Total DeltaH<br/>sum of stack enthalpies]
        DS[Total DeltaS<br/>sum of stack entropies]
        SALTCORR[Salt Correction<br/>Owczarzy formula]
        ADDCORR[Additive Correction<br/>empirical factors]
    end

    subgraph Output
        TM[Melting Temperature<br/>Tm = DH / DS]
        DG[Free Energy<br/>DG = DH - T*DS]
        PROB[Binding Probability<br/>Boltzmann weighted]
    end

    SEQ --> NN
    NN --> DH
    NN --> DS
    INIT --> DH
    INIT --> DS
    TERM --> DH
    TERM --> DS
    SYM --> DS

    SALT --> SALTCORR
    DH --> SALTCORR
    DS --> SALTCORR

    ADD --> ADDCORR
    SALTCORR --> ADDCORR

    ADDCORR --> TM
    ADDCORR --> DG
    TEMP --> DG
    DG --> PROB
```

### SantaLucia Parameters Table

| Stack | DeltaH (kcal/mol) | DeltaS (cal/mol*K) |
|-------|-------------------|---------------------|
| AA/TT | -7.9 | -22.2 |
| AT/TA | -7.2 | -20.4 |
| TA/AT | -7.2 | -21.3 |
| CA/GT | -8.5 | -22.7 |
| GT/CA | -8.4 | -22.4 |
| CT/GA | -7.8 | -21.0 |
| GA/CT | -8.2 | -22.2 |
| CG/GC | -10.6 | -27.2 |
| GC/CG | -9.8 | -24.4 |
| GG/CC | -8.0 | -19.9 |

---

## Performance Optimization

Caching and acceleration strategies.

```mermaid
flowchart TB
    subgraph "Without Cache"
        HDF5A[(HDF5 File)]
        READ1[Read Position 1]
        READ2[Read Position 2]
        READN[Read Position N]
        SLOW[~1000 reads/sec]
    end

    subgraph "With PositionCache"
        HDF5B[(HDF5 File)]
        LOAD[Single Load<br/>all positions]
        MEM[(In-Memory<br/>Dict Cache)]
        FAST[~1M lookups/sec]
    end

    subgraph "Memory Layout"
        KEY["(prefix, primer, strand)"]
        VAL[np.ndarray int32]
    end

    HDF5A --> READ1
    HDF5A --> READ2
    HDF5A --> READN
    READ1 --> SLOW
    READ2 --> SLOW
    READN --> SLOW

    HDF5B --> LOAD
    LOAD --> MEM
    MEM --> FAST

    KEY --> MEM
    VAL --> MEM
```

### GPU Acceleration

```mermaid
flowchart LR
    subgraph CPU
        PRIMERS[1000 Primers]
        CPU_CALC[Sequential Tm<br/>1000 iterations]
        CPU_TIME[~100ms]
    end

    subgraph GPU
        BATCH[Batch Transfer<br/>to GPU memory]
        GPU_CALC[Parallel Tm<br/>1000 threads]
        RESULT[Batch Transfer<br/>back to CPU]
        GPU_TIME[~1ms]
    end

    PRIMERS --> CPU_CALC
    CPU_CALC --> CPU_TIME

    PRIMERS --> BATCH
    BATCH --> GPU_CALC
    GPU_CALC --> RESULT
    RESULT --> GPU_TIME
```

### Bloom Filter for Background

```mermaid
flowchart TB
    subgraph "Traditional Approach"
        BGFILE[(Background k-mers<br/>3 GB file)]
        BGMEM[Load to Memory<br/>~8 GB RAM]
        LOOKUP1[Dict Lookup<br/>O(1)]
    end

    subgraph "Bloom Filter Approach"
        BUILD[One-time Build<br/>~30 min]
        BLOOM[(Bloom Filter<br/>~100 MB)]
        LOOKUP2[Membership Test<br/>O(1)]
        FP[False Positive<br/>~1%]
    end

    BGFILE --> BGMEM
    BGMEM --> LOOKUP1

    BGFILE --> BUILD
    BUILD --> BLOOM
    BLOOM --> LOOKUP2
    LOOKUP2 --> FP
```

---

## Sequence Diagrams

### Pipeline Execution

```mermaid
sequenceDiagram
    participant User
    participant CLI
    participant Jellyfish
    participant Filter
    participant RF as RandomForest
    participant Optimizer
    participant Cache

    User->>CLI: neoswga count-kmers -j params.json
    CLI->>Jellyfish: Count k-mers (fg genome)
    Jellyfish-->>CLI: fg_Xmer_all.txt
    CLI->>Jellyfish: Count k-mers (bg genome)
    Jellyfish-->>CLI: bg_Xmer_all.txt
    CLI->>CLI: Build position HDF5

    User->>CLI: neoswga filter -j params.json
    CLI->>Filter: Load k-mer files
    Filter->>Filter: Apply 5 filtering rules
    Filter->>Filter: Tm filtering
    Filter->>Filter: Gini calculation
    Filter-->>CLI: step2_df.csv

    User->>CLI: neoswga score -j params.json
    CLI->>Cache: Load positions
    CLI->>RF: Feature engineering
    RF->>RF: Predict amplification
    RF-->>CLI: step3_df.csv

    User->>CLI: neoswga optimize -j params.json
    CLI->>Cache: Load positions
    CLI->>Optimizer: Select method (hybrid)
    Optimizer->>Optimizer: Build network
    Optimizer->>Optimizer: Greedy refinement
    Optimizer->>Optimizer: Dimer check
    Optimizer-->>CLI: step4_improved_df.csv
    CLI-->>User: Final primer set
```

### Optimizer Selection

```mermaid
sequenceDiagram
    participant CLI
    participant Factory as OptimizerFactory
    participant Registry as OptimizerRegistry
    participant Optimizer

    CLI->>Factory: create('hybrid', cache, ...)
    Factory->>Registry: get('hybrid')
    Registry->>Registry: Check aliases
    Registry-->>Factory: HybridOptimizer class
    Factory->>Optimizer: __init__(cache, fg_prefixes, ...)
    Optimizer-->>Factory: optimizer instance
    Factory-->>CLI: optimizer

    CLI->>Optimizer: optimize(candidates, target_size=6)
    Optimizer->>Optimizer: Phase 1: Network coverage
    Optimizer->>Optimizer: Phase 2: Greedy refinement
    Optimizer->>Optimizer: Phase 3: Dimer optimization
    Optimizer-->>CLI: OptimizationResult
```

---

## Component Summary

| Component | Purpose | Key Files |
|-----------|---------|-----------|
| CLI | User interface | cli_unified.py |
| Pipeline | Workflow orchestration | pipeline.py, improved_pipeline.py |
| Filtering | Primer candidate selection | filter.py, adaptive_filters.py |
| Scoring | ML-based ranking | rf_preprocessing.py, primer_attributes.py |
| Optimization | Set selection (18 strategies) | optimizer_factory.py, base_optimizer.py, unified_optimizer.py |
| Thermodynamics | NN calculations | thermodynamics.py, reaction_conditions.py |
| Mechanistic | Four-pathway amplification model | mechanistic_model.py, additives.py |
| Performance | Speed optimization | position_cache.py, gpu_acceleration.py |
| Analysis | Post-processing | genome_analysis.py, dimer_network_analyzer.py |
| Simulation | Replication modeling | replication_simulator.py, swga_simulator.py |
| Reporting | Quality reports (A-F grades) | report/quality.py, report/executive_summary.py |
| Export | Synthesis ordering formats | export.py |

---

## See Also

- [API Reference](API_REFERENCE.md) - Complete API documentation
- [User Guide](user-guide.md) - Usage tutorials
- [Algorithms](development/algorithms.md) - Algorithm details
