"""The typed pipeline configuration object.

Extracted from `core/parameter.py` on 2026-10-28, which sat at 1596 lines
against a 1601 budget -- five lines of headroom, while the last new params.json
key cost six. `parameter` re-exports this name, so
`from neoswga.core.parameter import PipelineParameters` still resolves and
`isinstance` still holds: it is the same class object.

It moved rather than its neighbours because it is the one large top-level item
in that module that touches no module global. `get_current_config` reads
`globals()` and `set_from_config` writes into it, so both are bound to the
module they live in and cannot be moved this way.
"""

from __future__ import annotations

from dataclasses import dataclass, field


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
    max_gini: float = 0.7
    min_gini_sites: int = 3
    max_primer: int = 500
    # Unitless amplification prediction score (scale ~0-20). Combines Tm
    # optimality, GC content, 3' stability, and binding energy features.
    # Higher values indicate stronger predicted amplification; the default
    # of 10.0 retains above-average primers.
    min_amp_pred: float = 10.0

    # Dimer constraints
    max_dimer_bp: int = 3
    allow_dimer_relaxation: bool = False
    refinement_method: str = "network"

    swap_max_evaluations: int = 10000
    swap_max_seconds: float = 10.0
    max_self_dimer_bp: int = 4

    # GC content filtering
    gc_min: float = 0.375
    gc_max: float = 0.625
    gc_tolerance: float = 0.15
    genome_gc: float | None = None

    # Per-primer extension reach for set-cover selection and reported coverage.
    # None means "use the polymerase's realistic reach"; an explicit value here
    # or on --coverage-reach overrides it.
    coverage_reach: int | None = None

    # Candidate ranking in step 2. None means "use the read site's default",
    # which keeps one source of truth for it (pipeline.py).
    occupancy_ranking: bool = True
    occupancy_shortlist: int | None = None
    max_mismatches: int | None = None

    # Polymerase and reaction conditions
    polymerase: str = "phi29"
    reaction_temp: float | None = None
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
    propanediol_m: float = 0.0

    # Enzyme-acting additives (no Tm model; recorded and passed through)
    glycerol_percent: float = 0.0
    peg_percent: float = 0.0
    bsa_ug_ml: float = 0.0
    ssb: bool = False

    # Buffer species beyond Na+/Mg2+
    k_conc: float = 0.0
    nh4_conc: float = 0.0
    dntp_conc: float = 0.0
    dtt_mm: float = 0.0

    # Optimization settings
    iterations: int = 8
    max_sets: int = 5
    retries: int = 5
    top_set_count: int = 10
    drop_iterations: int = 5
    selection_metric: str = "deterministic"

    # Genome topology
    fg_circular: bool = True
    bg_circular: bool = False

    # K-mer sampling (for performance)
    sample_rate: float | None = None
    min_sample_count: int = 5

    # Bloom filter (for large backgrounds)
    use_bloom_filter: bool = False
    bloom_filter_path: str | None = None
    sampled_index_path: str | None = None

    # Paths
    data_dir: str = ""
    src_dir: str = ""

    # Genome data
    fg_genomes: list[str] = field(default_factory=list)
    bg_genomes: list[str] = field(default_factory=list)
    fg_prefixes: list[str] = field(default_factory=list)
    bg_prefixes: list[str] = field(default_factory=list)
    fg_seq_lengths: list[int] = field(default_factory=list)
    bg_seq_lengths: list[int] = field(default_factory=list)

    # Exclusion genomes (zero-tolerance filtering)
    excl_genomes: list[str] = field(default_factory=list)
    excl_prefixes: list[str] = field(default_factory=list)
    excl_threshold: int = 0  # Max allowed hits in exclusion genome (0 = any hit rejects)

    # Blacklist genomes (penalty-weighted filtering)
    bl_genomes: list[str] = field(default_factory=list)
    bl_prefixes: list[str] = field(default_factory=list)
    bl_seq_lengths: list[int] = field(default_factory=list)
    bl_penalty: float = 5.0
    max_bl_freq: float = 0.0

    # Runtime
    cpus: int = 1
    verbose: bool = False
