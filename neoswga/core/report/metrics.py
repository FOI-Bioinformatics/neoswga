"""
Metrics collection for SWGA pipeline reports.

Standardizes metrics from pipeline output files into a consistent format
for report generation.
"""

import csv
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)


def _get_version() -> str:
    """Get NeoSWGA version from package metadata."""
    try:
        from neoswga import __version__

        return __version__
    except ImportError:
        return "unknown"


def _normalize_amp_pred(raw_score: float) -> float:
    """Normalize RF amplification prediction from 0-20 scale to 0-1.

    The RF model (RandomForestRegressor) outputs scores on a ~0-20 scale.
    Report modules expect 0-1 for quality display (stars, heatmaps, grades).
    """
    if raw_score <= 0:
        return 0.0
    # Model outputs ~0-20; normalize to 0-1
    return min(1.0, raw_score / 20.0)


def _safe_float(value: Any, default: float = 0.0) -> float:
    """
    Safely convert a value to float.

    Args:
        value: Value to convert (str, int, float, or None)
        default: Default value if conversion fails

    Returns:
        Float value or default (NaN and Inf are treated as invalid)
    """
    import math

    if value is None or value == "":
        return default

    # Check for string representations of invalid values
    if isinstance(value, str):
        value_lower = value.lower().strip()
        if value_lower in ("nan", "inf", "-inf", "infinity", "-infinity", "none", "null"):
            return default

    try:
        result = float(value)
        # Check for NaN or Inf
        if math.isnan(result) or math.isinf(result):
            return default
        return result
    except (ValueError, TypeError) as e:
        logger.debug(f"Failed to convert '{value}' to float: {e}")
        return default


def _safe_int(value: Any, default: int = 0) -> int:
    """
    Safely convert a value to int.

    Args:
        value: Value to convert
        default: Default value if conversion fails

    Returns:
        Int value or default
    """
    if value is None or value == "":
        return default
    try:
        # Handle float strings like "123.0"
        return int(float(value))
    except (ValueError, TypeError) as e:
        logger.debug(f"Failed to convert '{value}' to int: {e}")
        return default


@dataclass
class GenomeInfo:
    """Information about a genome."""

    name: str
    size: int
    gc_content: float
    n_chromosomes: int = 1
    classification: str = "unknown"

    @property
    def size_mbp(self) -> float:
        """Size in megabase pairs."""
        return self.size / 1_000_000


@dataclass
class PrimerMetrics:
    """Metrics for a single primer."""

    sequence: str
    length: int
    gc_content: float
    tm: float
    fg_freq: float
    bg_freq: float
    fg_sites: int
    bg_sites: int
    gini: float
    specificity: float
    amp_pred: float = 0.0
    dimer_score: float = 0.0
    hairpin_dg: float = 0.0
    self_dimer_dg: float = 0.0
    three_prime_stability: float = 0.0
    strand_ratio: float = 1.0
    quality_rank: int = 0

    @classmethod
    def from_row(cls, row: Dict[str, str]) -> "PrimerMetrics":
        """
        Create PrimerMetrics from a CSV row.

        Uses safe conversions to handle malformed or missing data gracefully.
        """
        seq = str(row.get("primer", row.get("sequence", "")))

        # Calculate GC content if not provided
        gc = _safe_float(row.get("gc", row.get("gc_content", 0)))
        if gc == 0 and seq:
            # Use uppercase for case-insensitive GC counting
            seq_upper = seq.upper()
            gc = (seq_upper.count("G") + seq_upper.count("C")) / len(seq) if len(seq) > 0 else 0

        # Calculate specificity with safe division
        fg_freq = _safe_float(row.get("fg_freq", 0))
        bg_freq = _safe_float(row.get("bg_freq", 0))
        # Use max to avoid division by zero; 1e-10 represents "essentially no background"
        specificity = fg_freq / max(bg_freq, 1e-10)

        return cls(
            sequence=seq,
            length=len(seq),
            gc_content=gc,
            tm=_safe_float(row.get("tm", row.get("Tm", 0))),
            fg_freq=fg_freq,
            bg_freq=bg_freq,
            fg_sites=_safe_int(row.get("fg_count", row.get("fg_sites", 0))),
            bg_sites=_safe_int(row.get("bg_count", row.get("bg_sites", 0))),
            gini=_safe_float(row.get("gini", row.get("gini_index", 0))),
            specificity=specificity,
            amp_pred=_normalize_amp_pred(
                _safe_float(row.get("amp_pred", row.get("on.target.pred", 0)))
            ),
            dimer_score=_safe_float(row.get("dimer_score", row.get("dimer_risk_score", 0))),
            hairpin_dg=_safe_float(row.get("hairpin_dg", 0)),
            self_dimer_dg=_safe_float(row.get("self_dimer_dg", 0)),
            three_prime_stability=_safe_float(row.get("three_prime_stability", 0)),
            strand_ratio=_safe_float(row.get("strand_ratio"), default=1.0),
        )


@dataclass
class FilteringStats:
    """Statistics from the filtering step."""

    total_kmers: int = 0
    after_frequency: int = 0
    after_background: int = 0
    after_gini: int = 0
    after_thermodynamic: int = 0
    after_complexity: int = 0
    final_candidates: int = 0

    def as_funnel(self) -> List[tuple]:
        """Return filtering stages as a funnel list, in actual pipeline order.

        Order matches how the filter step applies stages (frequency ->
        thermodynamic/sequence-quality -> exclusion/blacklist -> Gini -> cap),
        so the rendered funnel is monotonically non-increasing. Intermediate
        stages with a 0 count (not recorded) are skipped; the total and final
        are always shown.
        """
        ordered = [
            ("Total k-mers", self.total_kmers),
            ("After frequency filter", self.after_frequency),
            ("After thermodynamic filter", self.after_thermodynamic),
            ("After background/blacklist", self.after_background),
            ("After Gini filter", self.after_gini),
            ("After complexity filter", self.after_complexity),
            ("Final candidates", self.final_candidates),
        ]
        # Keep the endpoints; drop zero-count intermediate stages.
        return [
            (label, count)
            for i, (label, count) in enumerate(ordered)
            if i in (0, len(ordered) - 1) or count > 0
        ]


@dataclass
class CoverageMetrics:
    """Coverage analysis metrics."""

    overall_coverage: float = 0.0
    covered_bases: int = 0
    total_bases: int = 0
    n_gaps: int = 0
    largest_gap: int = 0
    critical_gaps: int = 0  # >100kb
    high_gaps: int = 0  # 50-100kb
    medium_gaps: int = 0  # 20-50kb
    low_gaps: int = 0  # <20kb
    gap_locations: List[Dict] = field(default_factory=list)
    mean_gap: float = 0.0
    max_gap: float = 0.0
    gap_gini: float = 0.0
    gap_entropy: float = 0.0
    from_optimizer: bool = False
    # Per-primer extension reach (bp) used by the optimizer to compute
    # overall_coverage. None when the value comes from the fallback
    # estimate (n_primers * 30 kb) instead of a measured optimizer run.
    extension_reach: Optional[int] = None
    # Per-target coverage (multi-genome runs): maps fg prefix -> coverage.
    # Empty for single-genome runs. Read from the optimizer summary.
    per_target_coverage: Dict[str, float] = field(default_factory=dict)
    # Measured coverage uniformity (Gini of gap sizes) from the optimizer.
    coverage_uniformity: Optional[float] = None


@dataclass
class SpecificityMetrics:
    """Specificity analysis metrics."""

    enrichment_ratio: float = 0.0
    target_sites: int = 0
    background_sites: int = 0
    target_density: float = 0.0  # sites per Mbp
    background_density: float = 0.0
    # Measured values from the optimizer summary (preferred over the
    # density-ratio estimate above when available).
    selectivity_ratio: Optional[float] = None
    bg_coverage: Optional[float] = None
    from_optimizer: bool = False


@dataclass
class ThermodynamicMetrics:
    """Thermodynamic analysis metrics."""

    mean_tm: float = 0.0
    min_tm: float = 0.0
    max_tm: float = 0.0
    tm_range: float = 0.0
    reaction_temp: float = 30.0
    polymerase: str = "phi29"
    max_heterodimer_dg: float = 0.0
    dimer_risk_level: str = "low"


@dataclass
class UniformityMetrics:
    """Binding uniformity metrics."""

    mean_gini: float = 0.0
    max_gini: float = 0.0
    mean_sites_per_bin: float = 0.0
    sites_std: float = 0.0
    coefficient_of_variation: float = 0.0
    forward_sites: int = 0
    reverse_sites: int = 0
    strand_ratio: float = 1.0
    # Measured strand-interleaving metrics from the optimizer (None = not
    # available). strand_alternation_score: fraction of adjacent binding sites
    # that alternate strands; strand_coverage_ratio: min/max fwd-vs-rev balance.
    strand_alternation_score: Optional[float] = None
    strand_coverage_ratio: Optional[float] = None
    # Set-level binding-evenness Gini measured by the optimizer (None = not
    # available). This is the quantity `quality_thresholds.GINI` is calibrated
    # on -- published sets, not individual primers -- so grading prefers it
    # over `max_gini`, which is the worst single primer and always worse.
    measured_gini: Optional[float] = None
    from_optimizer: bool = False


@dataclass
class PipelineMetrics:
    """Complete metrics from a pipeline run."""

    # Metadata
    results_dir: str = ""
    generated_at: str = ""
    pipeline_version: str = field(default_factory=_get_version)

    # Genome info
    target_genome: Optional[GenomeInfo] = None
    background_genome: Optional[GenomeInfo] = None

    # Parameters
    parameters: Dict[str, Any] = field(default_factory=dict)

    # Primers
    primers: List[PrimerMetrics] = field(default_factory=list)
    primer_count: int = 0

    # Step-by-step metrics
    filtering: Optional[FilteringStats] = None
    coverage: Optional[CoverageMetrics] = None
    specificity: Optional[SpecificityMetrics] = None
    thermodynamics: Optional[ThermodynamicMetrics] = None
    uniformity: Optional[UniformityMetrics] = None

    # Optimizer extras read from step4_improved_df_summary.json (authoritative).
    # ensemble_comparison: per-method rows from an --optimization-method=ensemble
    # run; pareto_metrics: set-level Pareto solutions from MOEA / frontier.
    ensemble_comparison: List[Dict] = field(default_factory=list)
    pareto_metrics: List[Dict] = field(default_factory=list)
    # Coverage gaps from the analyze-coverage command (coverage_gaps.json /
    # merged_gaps.bed), when present in the results dir. None = not run.
    coverage_gaps: Optional[Dict] = None
    # Reaction conditions / additives used (from params.json), for the report.
    reaction_conditions: Dict[str, Any] = field(default_factory=dict)

    # Runtime
    total_runtime_seconds: float = 0.0
    step_runtimes: Dict[str, float] = field(default_factory=dict)

    # Post-optimization validator report loaded from
    # step4_improved_df_validation.json when available. Each issue is a
    # dict {"level": "error"|"warning"|"info", "code": str, "detail": str}.
    # The report module surfaces these in the HTML so users do not have to
    # open the JSON to see per_target_coverage_below_threshold,
    # blacklist_primer_in_set, duplicate-primer, etc. warnings.
    validation_issues: List[Dict[str, str]] = field(default_factory=list)
    validation_ok: bool = True


def _load_csv(filepath: Path) -> List[Dict[str, str]]:
    """
    Load CSV file as list of dictionaries.

    Args:
        filepath: Path to CSV file

    Returns:
        List of dictionaries (one per row), empty list on error
    """
    if not filepath.exists():
        return []

    try:
        with open(filepath, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            return list(reader)
    except (csv.Error, UnicodeDecodeError, PermissionError, OSError) as e:
        logger.warning(f"Failed to load CSV file {filepath}: {e}")
        return []


def _load_params(results_dir: Path) -> Dict[str, Any]:
    """
    Load params.json if present.

    Searches in results_dir first, then parent directory.

    Args:
        results_dir: Path to results directory

    Returns:
        Dictionary of parameters, empty dict on error
    """
    params_file = results_dir / "params.json"
    if not params_file.exists():
        # Try parent directory
        params_file = results_dir.parent / "params.json"

    if not params_file.exists():
        return {}

    try:
        with open(params_file, encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, UnicodeDecodeError, PermissionError, OSError) as e:
        logger.warning(f"Failed to load params file {params_file}: {e}")
        return {}


def _load_optimizer_summary(results_path: Path) -> Optional[Dict]:
    """Load step4 summary JSON if available."""
    summary_file = results_path / "step4_improved_df_summary.json"
    if not summary_file.exists():
        return None
    try:
        with open(summary_file, encoding="utf-8") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to load optimizer summary: {e}")
        return None


def _load_validation_report(results_path: Path) -> tuple:
    """Load step4 validator report written by OptimizationResult.validate().

    Returns:
        (issues, ok) where issues is a list of {"level","code","detail"}
        dicts and ok is True when the validator reported no errors. Missing
        or unreadable files yield ([], True) so downstream report code can
        treat the absence as "no issues to surface".
    """
    report_file = results_path / "step4_improved_df_validation.json"
    if not report_file.exists():
        return [], True
    try:
        with open(report_file, encoding="utf-8") as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to load validation report: {e}")
        return [], True
    issues = data.get("issues", []) if isinstance(data, dict) else []
    ok = bool(data.get("ok", True)) if isinstance(data, dict) else True
    # Normalise: ensure each issue has level/code/detail strings.
    normalised = []
    for it in issues:
        if not isinstance(it, dict):
            continue
        normalised.append(
            {
                "level": str(it.get("level", "warning")),
                "code": str(it.get("code", "unknown")),
                "detail": str(it.get("detail", "")),
            }
        )
    return normalised, ok


def _genome_size(params: Dict, prefix: str) -> int:
    """Total length of the `fg`/`bg` genomes, in bp, or 0 if params omit it.

    `fg_size` is not a key the pipeline writes; `fg_seq_lengths` is. Reading
    only the former left coverage with a genome of length 0 and the specificity
    density dividing by a placeholder 1.
    """
    legacy = {"fg": "foreground_size", "bg": "background_size"}.get(prefix, "")
    size = _safe_int(params.get(f"{prefix}_size") or params.get(legacy) or 0)
    if size:
        return size
    lengths = params.get(f"{prefix}_seq_lengths") or []
    if isinstance(lengths, (list, tuple)):
        return sum(_safe_int(length) for length in lengths)
    return 0


def _extract_genome_info(params: Dict, prefix: str) -> Optional[GenomeInfo]:
    """Extract genome info from params or file metadata.

    The plural keys are the ones the pipeline actually writes: params.json
    carries `fg_genomes` and `fg_seq_lengths`, never `fg_genome` or `fg_size`.
    Reading only the singular forms left the report with no genome -- no name,
    no size, and `covered_bases` pinned at 0 because it is guarded on
    `total_bases`.
    """
    genome_file = params.get(f"{prefix}_genome") or params.get(prefix) or ""
    if not genome_file:
        genomes = params.get(f"{prefix}_genomes") or params.get(f"{prefix}_prefixes") or []
        if isinstance(genomes, (list, tuple)) and genomes:
            genome_file = genomes[0]
        elif isinstance(genomes, str):
            genome_file = genomes
    if not genome_file:
        return None

    size = _genome_size(params, prefix)

    gc = _safe_float(params.get(f"{prefix}_gc", 0))

    n_chromosomes = 1
    genomes = params.get(f"{prefix}_genomes")
    if isinstance(genomes, (list, tuple)) and genomes:
        n_chromosomes = len(genomes)

    return GenomeInfo(
        name=Path(str(genome_file)).stem,
        size=size,
        gc_content=gc,
        n_chromosomes=n_chromosomes,
    )


def _calculate_coverage_metrics(
    primers: List[PrimerMetrics],
    params: Dict,
) -> CoverageMetrics:
    """Calculate coverage metrics from primer data."""
    coverage = CoverageMetrics()

    # Note: PrimerMetrics does not currently have a coverage attribute.
    # Coverage would need to be calculated from position data (HDF5 files)
    # or provided in the step4 output. For now, we estimate from primer count
    # and genome size using a simplified model.

    # Get genome size
    fg_size = _genome_size(params, "fg")
    if fg_size > 0:
        coverage.total_bases = fg_size

        # Estimate coverage based on number of primers and expected amplification
        # Typical phi29 amplification range is ~10-50 kb per primer
        # This is a rough estimate; actual coverage requires position data
        if primers:
            n_primers = len(primers)
            avg_amp_range = 30000  # 30 kb average amplification range
            estimated_covered = min(n_primers * avg_amp_range, fg_size)
            coverage.overall_coverage = estimated_covered / fg_size
            coverage.covered_bases = int(estimated_covered)
            logger.info("Coverage estimated from primer count (no optimizer summary available)")

    return coverage


def _calculate_specificity_metrics(
    primers: List[PrimerMetrics],
    params: Dict,
) -> SpecificityMetrics:
    """Calculate specificity metrics from primer data."""
    if not primers:
        return SpecificityMetrics()

    # Sum up sites
    target_sites = sum(p.fg_sites for p in primers)
    bg_sites = sum(p.bg_sites for p in primers)

    # Calculate densities
    fg_size = _genome_size(params, "fg")
    bg_size = _genome_size(params, "bg")

    target_density = (target_sites / fg_size) * 1_000_000 if fg_size > 0 else 0
    bg_density = (bg_sites / bg_size) * 1_000_000 if bg_size > 0 else 0

    # Calculate enrichment
    enrichment = target_density / max(bg_density, 1e-10)

    return SpecificityMetrics(
        enrichment_ratio=enrichment,
        target_sites=target_sites,
        background_sites=bg_sites,
        target_density=target_density,
        background_density=bg_density,
    )


def _calculate_thermodynamic_metrics(
    primers: List[PrimerMetrics],
    params: Dict,
) -> ThermodynamicMetrics:
    """Calculate thermodynamic metrics from primer data."""
    if not primers:
        return ThermodynamicMetrics()

    tms = [p.tm for p in primers if p.tm > 0]
    if not tms:
        return ThermodynamicMetrics()

    # Find worst (most negative) heterodimer dG
    # Note: dimer scores are negative, more negative = worse
    worst_dimer_dg = min((p.dimer_score for p in primers), default=0)

    # Assess dimer risk level based on worst dG
    if worst_dimer_dg < -8.0:
        risk_level = "high"
    elif worst_dimer_dg < -5.0:
        risk_level = "moderate"
    else:
        risk_level = "low"

    return ThermodynamicMetrics(
        mean_tm=sum(tms) / len(tms),
        min_tm=min(tms),
        max_tm=max(tms),
        tm_range=max(tms) - min(tms),
        reaction_temp=_safe_float(params.get("reaction_temp", 30.0), default=30.0),
        polymerase=str(params.get("polymerase", "phi29")),
        max_heterodimer_dg=worst_dimer_dg,
        dimer_risk_level=risk_level,
    )


def _calculate_uniformity_metrics(
    primers: List[PrimerMetrics],
) -> UniformityMetrics:
    """Calculate uniformity metrics from primer data."""
    if not primers:
        return UniformityMetrics()

    ginis = [p.gini for p in primers if p.gini > 0]

    mean_gini = sum(ginis) / len(ginis) if ginis else 0
    max_gini = max(ginis) if ginis else 0

    # Strand balance
    # Count primers favoring forward vs reverse strand
    forward = sum(1 for p in primers if p.strand_ratio >= 1.0)
    reverse = len(primers) - forward

    # Calculate strand ratio with meaningful interpretation:
    # - 1.0 = balanced (equal forward and reverse)
    # - >1.0 = forward bias
    # - <1.0 = reverse bias
    # Clamp to prevent extreme values when one count is 0
    if forward == 0 and reverse == 0:
        strand_ratio = 1.0  # No data
    elif reverse == 0:
        strand_ratio = min(forward, 10.0)  # Cap at 10 for all-forward
    elif forward == 0:
        strand_ratio = max(1.0 / reverse, 0.1)  # Floor at 0.1 for all-reverse
    else:
        strand_ratio = forward / reverse

    return UniformityMetrics(
        mean_gini=mean_gini,
        max_gini=max_gini,
        strand_ratio=strand_ratio,
    )


def _load_coverage_gaps(results_path: Path) -> Optional[Dict]:
    """Load coverage gaps written by `neoswga analyze-coverage`, if present.

    Reads ``coverage_gaps.json`` (preferred) which records the merged in-silico
    and BAM-derived gaps. Returns a dict with ``n_gaps``, ``used_bam``,
    ``min_depth``, ``min_gap_size`` and a ``gaps`` list, or None when the
    analyze-coverage step was not run.
    """
    gaps_json = results_path / "coverage_gaps.json"
    if not gaps_json.exists():
        return None
    try:
        with open(gaps_json, encoding="utf-8") as f:
            data = json.load(f)
        if isinstance(data, dict) and "gaps" in data:
            return data
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to load coverage_gaps.json: {e}")
    return None


_CONDITION_ADDITIVES = (
    "dmso_percent",
    "betaine_m",
    "trehalose_m",
    "formamide_percent",
    "ethanol_percent",
    "urea_m",
    "tmac_m",
)


def _recorded_conditions(results_path: Optional[Path]) -> Optional[Dict[str, Any]]:
    """Reaction conditions the run manifest recorded, if there is one."""
    if results_path is None:
        return None
    try:
        from neoswga.core.run_manifest import read_effective_conditions
    except ImportError:  # pragma: no cover - report can render without core
        return None
    return read_effective_conditions(str(results_path))


def _conditions_from_params(params: Dict[str, Any], results_path: Optional[Path] = None):
    """The reaction the design was run under, or None if params cannot name one.

    Used to recompute quantities the results CSV does not carry, so the report
    states them in the same terms the optimizer did.

    The run manifest wins over params.json where it has an opinion: the
    GC-adaptive strategy resolves betaine and DMSO from the target's GC at run
    time and never writes them back to the file, so params.json alone names a
    reaction the design was not optimized under.
    """
    recorded = _recorded_conditions(results_path)
    if recorded:
        params = {**(params or {}), **recorded}

    if not params:
        return None
    try:
        from neoswga.core.parameter import default_reaction_temp
        from neoswga.core.reaction_conditions import ReactionConditions
    except ImportError:  # pragma: no cover - report can render without core
        return None

    polymerase = str(params.get("polymerase", "phi29"))
    additives = {name: params[name] for name in _CONDITION_ADDITIVES if params.get(name)}
    try:
        return ReactionConditions(
            temp=_safe_float(
                params.get("reaction_temp"), default=default_reaction_temp(polymerase)
            ),
            polymerase=polymerase,
            mg_conc=params.get("mg_conc"),
            **additives,
        )
    except (ValueError, TypeError) as e:
        logger.warning(f"Could not reconstruct reaction conditions from params: {e}")
        return None


def _backfill_primer_tm(
    primers: List["PrimerMetrics"],
    params: Dict[str, Any],
    results_path: Optional[Path] = None,
) -> None:
    """Give each primer the Tm its row does not carry.

    `step4_improved_df.csv` holds a set-level `mean_tm` and no per-primer
    column, so `PrimerMetrics.from_row` fell through to its `0` default and the
    report printed a Tm of 0.0 C for every primer in the table. The sequence and
    the reaction conditions are both in hand here, so the value is computable --
    and computing it from `calculate_effective_tm` keeps it the same quantity
    the optimizer and the exporter report.

    Only a physically meaningful result is written back. Nearest-neighbour on a
    sequence far below primer length returns a Tm below freezing, and swapping
    one meaningless number for another is not an improvement.
    """
    missing = [p for p in primers if not p.tm > 0 and p.sequence]
    if not missing:
        return

    conditions = _conditions_from_params(params, results_path)
    if conditions is not None:
        estimate = conditions.calculate_effective_tm
    else:
        try:
            from neoswga.core.thermodynamics import calculate_tm_basic as estimate
        except ImportError:  # pragma: no cover - leave the zeros rather than guess
            return

    for primer in missing:
        try:
            tm = float(estimate(primer.sequence))
        except (ValueError, TypeError, KeyError):
            continue
        if tm > 0:
            primer.tm = tm


def _backfill_primer_sites(
    primers: List["PrimerMetrics"], results_path: Path, params: Dict[str, Any]
) -> None:
    """Give each primer the site counts its row does not carry.

    `step4_improved_df.csv` records the set-level result, not per-primer
    measurements, so `fg_sites`, `bg_sites` and `gini` all fell through to 0
    and per-primer specificity came out `0 / 1e-10 = 0`. The report then showed
    "0x" against every primer, including ones with no background sites at all
    -- the best case rendered as the worst.

    The filter step already measured these, per primer, in `step2_df.csv`. This
    reads them from there rather than inventing them.
    """
    missing = [p for p in primers if p.sequence and not p.fg_sites and not p.bg_sites]
    if not missing:
        return

    step2 = results_path / "step2_df.csv"
    if not step2.exists():
        return
    by_primer = {str(row.get("primer", row.get("sequence", ""))): row for row in _load_csv(step2)}
    if not by_primer:
        return

    fg_len = sum(params.get("fg_seq_lengths") or []) or 0
    bg_len = sum(params.get("bg_seq_lengths") or []) or 0

    for primer in missing:
        row = by_primer.get(primer.sequence)
        if row is None:
            continue
        primer.fg_sites = _safe_int(row.get("fg_count", row.get("fg_sites", 0)))
        primer.bg_sites = _safe_int(row.get("bg_count", row.get("bg_sites", 0)))
        if not primer.gini:
            primer.gini = _safe_float(row.get("gini", row.get("gini_index", 0)))

        if fg_len > 0 and bg_len > 0 and primer.fg_sites:
            # The same density ratio enrichment is graded on. A primer with no
            # observed background is floored at one site: that is the detection
            # limit of this background, and claiming more specificity than the
            # data can show is how a zero became a certainty elsewhere.
            fg_density = primer.fg_sites / fg_len
            bg_density = max(primer.bg_sites, 1) / bg_len
            primer.specificity = fg_density / bg_density


def _extract_reaction_conditions(params: Dict[str, Any]) -> Dict[str, Any]:
    """Pull the reaction conditions / additives actually used from params.

    Surfaced in the report so a reader can see the polymerase, temperature,
    salts, and additive concentrations behind the design. Only non-default
    additive values are kept so the block stays focused.
    """
    if not params:
        return {}
    conditions: Dict[str, Any] = {}
    for key in ("polymerase", "reaction_temp", "na_conc", "mg_conc", "primer_conc"):
        if params.get(key) is not None:
            conditions[key] = params[key]
    additives = {}
    for key in (
        "dmso_percent",
        "betaine_m",
        "trehalose_m",
        "formamide_percent",
        "ethanol_percent",
        "urea_m",
        "tmac_m",
    ):
        val = params.get(key)
        if val:  # only non-zero additives
            additives[key] = val
    if additives:
        conditions["additives"] = additives
    return conditions


def collect_pipeline_metrics(results_dir: str) -> PipelineMetrics:
    """
    Collect all metrics from a pipeline results directory.

    Args:
        results_dir: Path to results directory containing step2_df.csv,
                    step3_df.csv, step4_improved_df.csv, etc.

    Returns:
        PipelineMetrics with all collected data
    """
    results_path = Path(results_dir)

    if not results_path.exists():
        raise FileNotFoundError(f"Results directory not found: {results_dir}")

    logger.info(f"Collecting metrics from {results_dir}")

    # Initialize metrics
    metrics = PipelineMetrics(
        results_dir=str(results_path.absolute()),
        generated_at=datetime.now().isoformat(),
    )

    # Load parameters
    metrics.parameters = _load_params(results_path)

    # Extract genome info
    metrics.target_genome = _extract_genome_info(metrics.parameters, "fg")
    metrics.background_genome = _extract_genome_info(metrics.parameters, "bg")

    # Load primer data (prefer step4, fall back to step3, then step2)
    primer_rows = []
    for step_file in ["step4_improved_df.csv", "step3_df.csv", "step2_df.csv"]:
        filepath = results_path / step_file
        if filepath.exists():
            primer_rows = _load_csv(filepath)
            logger.info(f"Loaded {len(primer_rows)} primers from {step_file}")
            break

    # Convert to PrimerMetrics
    metrics.primers = [PrimerMetrics.from_row(row) for row in primer_rows]
    metrics.primer_count = len(metrics.primers)
    _backfill_primer_tm(metrics.primers, metrics.parameters, results_path)
    _backfill_primer_sites(metrics.primers, results_path, metrics.parameters)

    # Calculate derived metrics
    metrics.coverage = _calculate_coverage_metrics(metrics.primers, metrics.parameters)
    metrics.specificity = _calculate_specificity_metrics(metrics.primers, metrics.parameters)
    metrics.thermodynamics = _calculate_thermodynamic_metrics(metrics.primers, metrics.parameters)
    metrics.uniformity = _calculate_uniformity_metrics(metrics.primers)

    # Load optimizer summary for real metrics
    optimizer_summary = _load_optimizer_summary(results_path)
    if optimizer_summary and "metrics" in optimizer_summary:
        opt_metrics = optimizer_summary["metrics"]
        # Override estimated coverage with real optimizer data
        if "fg_coverage" in opt_metrics:
            metrics.coverage.overall_coverage = opt_metrics["fg_coverage"]
            metrics.coverage.from_optimizer = True
            metrics.coverage.extension_reach = opt_metrics.get("extension_reach")
            if metrics.coverage.total_bases > 0:
                metrics.coverage.covered_bases = int(
                    opt_metrics["fg_coverage"] * metrics.coverage.total_bases
                )
        # Add gap metrics
        metrics.coverage.mean_gap = opt_metrics.get("mean_gap", 0.0)
        metrics.coverage.max_gap = opt_metrics.get("max_gap", 0.0)
        metrics.coverage.gap_gini = opt_metrics.get("gap_gini", 0.0)
        metrics.coverage.gap_entropy = opt_metrics.get("gap_entropy", 0.0)
        # Per-target coverage (multi-genome) and measured coverage uniformity.
        if isinstance(opt_metrics.get("per_target_coverage"), dict):
            metrics.coverage.per_target_coverage = opt_metrics["per_target_coverage"]
        if "coverage_uniformity" in opt_metrics:
            metrics.coverage.coverage_uniformity = _safe_float(opt_metrics["coverage_uniformity"])

        # Prefer the optimizer's MEASURED selectivity / background coverage / site
        # totals over the density-ratio estimate derived from the CSV.
        if "selectivity_ratio" in opt_metrics:
            metrics.specificity.selectivity_ratio = _safe_float(opt_metrics["selectivity_ratio"])
            metrics.specificity.from_optimizer = True
        # Enrichment as graded here is a site-DENSITY ratio (target sites per Mb
        # over background sites per Mb), which is what `selectivity_density`
        # measures. Deriving it from the CSV needed per-primer fg_count/bg_count
        # columns that step4 does not have, so it came out 0.0 and a set with no
        # background sites -- the best possible result -- was graded "Critical"
        # on 30% of the composite score.
        if "selectivity_density" in opt_metrics:
            metrics.specificity.enrichment_ratio = _safe_float(opt_metrics["selectivity_density"])
        if "bg_coverage" in opt_metrics:
            metrics.specificity.bg_coverage = _safe_float(opt_metrics["bg_coverage"])
        if "total_fg_sites" in opt_metrics:
            metrics.specificity.target_sites = _safe_int(opt_metrics["total_fg_sites"])
        if "total_bg_sites" in opt_metrics:
            metrics.specificity.background_sites = _safe_int(opt_metrics["total_bg_sites"])

        # Measured strand-interleaving metrics.
        if metrics.uniformity is not None:
            if "strand_alternation_score" in opt_metrics:
                metrics.uniformity.strand_alternation_score = _safe_float(
                    opt_metrics["strand_alternation_score"]
                )
                metrics.uniformity.from_optimizer = True
            if "strand_coverage_ratio" in opt_metrics:
                metrics.uniformity.strand_coverage_ratio = _safe_float(
                    opt_metrics["strand_coverage_ratio"]
                )
            # The set-level evenness the grade is calibrated against.
            for key in ("gap_gini", "coverage_uniformity"):
                if key in opt_metrics:
                    metrics.uniformity.measured_gini = _safe_float(opt_metrics[key])
                    break
        logger.info("Using real optimizer metrics for coverage/specificity/strand data")

    # Top-level optimizer extras: ensemble per-method comparison and Pareto front.
    if optimizer_summary:
        if isinstance(optimizer_summary.get("ensemble_comparison"), list):
            metrics.ensemble_comparison = optimizer_summary["ensemble_comparison"]
        if isinstance(optimizer_summary.get("pareto_metrics"), list):
            metrics.pareto_metrics = optimizer_summary["pareto_metrics"]

    # Coverage gaps from the analyze-coverage command, if present.
    metrics.coverage_gaps = _load_coverage_gaps(results_path)

    # Reaction conditions / additives used (from params), for the report.
    metrics.reaction_conditions = _extract_reaction_conditions(
        {**(metrics.parameters or {}), **(_recorded_conditions(results_path) or {})}
    )

    # Load validator report so the HTML can surface warnings
    # (per_target_coverage_below_threshold, blacklist_primer_in_set,
    # set_size_mismatch, duplicate_primers, ...) that the optimizer
    # writes to JSON.
    metrics.validation_issues, metrics.validation_ok = _load_validation_report(results_path)

    # Try to load filtering stats if available
    filter_stats_file = results_path / "filter_stats.json"
    if filter_stats_file.exists():
        try:
            with open(filter_stats_file, encoding="utf-8") as f:
                stats = json.load(f)
                # Only use fields that FilteringStats accepts
                valid_fields = {
                    "total_kmers",
                    "after_frequency",
                    "after_background",
                    "after_gini",
                    "after_thermodynamic",
                    "after_complexity",
                    "final_candidates",
                }
                filtered_stats = {k: v for k, v in stats.items() if k in valid_fields}
                metrics.filtering = FilteringStats(**filtered_stats)
        except (json.JSONDecodeError, TypeError, KeyError, OSError) as e:
            logger.warning(f"Failed to load filter stats: {e}")
            # Continue without filtering stats

    return metrics
