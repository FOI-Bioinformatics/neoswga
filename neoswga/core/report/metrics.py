"""
Metrics collection for SWGA pipeline reports.

Standardizes metrics from pipeline output files into a consistent format
for report generation.
"""

import csv
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from datetime import datetime

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

    if value is None or value == '':
        return default

    # Check for string representations of invalid values
    if isinstance(value, str):
        value_lower = value.lower().strip()
        if value_lower in ('nan', 'inf', '-inf', 'infinity', '-infinity', 'none', 'null'):
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
    if value is None or value == '':
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
    def from_row(cls, row: Dict[str, str]) -> 'PrimerMetrics':
        """
        Create PrimerMetrics from a CSV row.

        Uses safe conversions to handle malformed or missing data gracefully.
        """
        seq = str(row.get('primer', row.get('sequence', '')))

        # Calculate GC content if not provided
        gc = _safe_float(row.get('gc', row.get('gc_content', 0)))
        if gc == 0 and seq:
            # Use uppercase for case-insensitive GC counting
            seq_upper = seq.upper()
            gc = (seq_upper.count('G') + seq_upper.count('C')) / len(seq) if len(seq) > 0 else 0

        # Calculate specificity with safe division
        fg_freq = _safe_float(row.get('fg_freq', 0))
        bg_freq = _safe_float(row.get('bg_freq', 0))
        # Use max to avoid division by zero; 1e-10 represents "essentially no background"
        specificity = fg_freq / max(bg_freq, 1e-10)

        return cls(
            sequence=seq,
            length=len(seq),
            gc_content=gc,
            tm=_safe_float(row.get('tm', row.get('Tm', 0))),
            fg_freq=fg_freq,
            bg_freq=bg_freq,
            fg_sites=_safe_int(row.get('fg_count', row.get('fg_sites', 0))),
            bg_sites=_safe_int(row.get('bg_count', row.get('bg_sites', 0))),
            gini=_safe_float(row.get('gini', row.get('gini_index', 0))),
            specificity=specificity,
            amp_pred=_normalize_amp_pred(
                _safe_float(row.get('amp_pred', row.get('on.target.pred', 0)))
            ),
            dimer_score=_safe_float(row.get('dimer_score', row.get('dimer_risk_score', 0))),
            hairpin_dg=_safe_float(row.get('hairpin_dg', 0)),
            self_dimer_dg=_safe_float(row.get('self_dimer_dg', 0)),
            three_prime_stability=_safe_float(row.get('three_prime_stability', 0)),
            strand_ratio=_safe_float(row.get('strand_ratio'), default=1.0),
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
        """Return filtering stages as a funnel list."""
        return [
            ("Total k-mers", self.total_kmers),
            ("After frequency filter", self.after_frequency),
            ("After background filter", self.after_background),
            ("After Gini filter", self.after_gini),
            ("After thermodynamic filter", self.after_thermodynamic),
            ("After complexity filter", self.after_complexity),
            ("Final candidates", self.final_candidates),
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
    high_gaps: int = 0      # 50-100kb
    medium_gaps: int = 0    # 20-50kb
    low_gaps: int = 0       # <20kb
    gap_locations: List[Dict] = field(default_factory=list)
    mean_gap: float = 0.0
    max_gap: float = 0.0
    gap_gini: float = 0.0
    gap_entropy: float = 0.0
    from_optimizer: bool = False


@dataclass
class SpecificityMetrics:
    """Specificity analysis metrics."""
    enrichment_ratio: float = 0.0
    target_sites: int = 0
    background_sites: int = 0
    target_density: float = 0.0  # sites per Mbp
    background_density: float = 0.0


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
        with open(filepath, newline='', encoding='utf-8') as f:
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
    params_file = results_dir / 'params.json'
    if not params_file.exists():
        # Try parent directory
        params_file = results_dir.parent / 'params.json'

    if not params_file.exists():
        return {}

    try:
        with open(params_file, encoding='utf-8') as f:
            return json.load(f)
    except (json.JSONDecodeError, UnicodeDecodeError, PermissionError, OSError) as e:
        logger.warning(f"Failed to load params file {params_file}: {e}")
        return {}


def _load_optimizer_summary(results_path: Path) -> Optional[Dict]:
    """Load step4 summary JSON if available."""
    summary_file = results_path / 'step4_improved_df_summary.json'
    if not summary_file.exists():
        return None
    try:
        with open(summary_file, encoding='utf-8') as f:
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
    report_file = results_path / 'step4_improved_df_validation.json'
    if not report_file.exists():
        return [], True
    try:
        with open(report_file, encoding='utf-8') as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        logger.warning(f"Failed to load validation report: {e}")
        return [], True
    issues = data.get('issues', []) if isinstance(data, dict) else []
    ok = bool(data.get('ok', True)) if isinstance(data, dict) else True
    # Normalise: ensure each issue has level/code/detail strings.
    normalised = []
    for it in issues:
        if not isinstance(it, dict):
            continue
        normalised.append({
            'level': str(it.get('level', 'warning')),
            'code': str(it.get('code', 'unknown')),
            'detail': str(it.get('detail', '')),
        })
    return normalised, ok


def _extract_genome_info(params: Dict, prefix: str) -> Optional[GenomeInfo]:
    """Extract genome info from params or file metadata."""
    genome_file = params.get(f'{prefix}_genome', params.get(f'{prefix}', ''))
    if not genome_file:
        return None

    # Try to get size and GC from params
    size = params.get(f'{prefix}_size', 0)
    gc = params.get(f'{prefix}_gc', 0)

    name = Path(genome_file).stem if genome_file else prefix

    return GenomeInfo(
        name=name,
        size=size,
        gc_content=gc,
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
    fg_size = params.get('fg_size', params.get('foreground_size', 0))
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
    fg_size = params.get('fg_size', params.get('foreground_size', 1))
    bg_size = params.get('bg_size', params.get('background_size', 1))

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
        reaction_temp=_safe_float(params.get('reaction_temp', 30.0), default=30.0),
        polymerase=str(params.get('polymerase', 'phi29')),
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
    metrics.target_genome = _extract_genome_info(metrics.parameters, 'fg')
    metrics.background_genome = _extract_genome_info(metrics.parameters, 'bg')

    # Load primer data (prefer step4, fall back to step3, then step2)
    primer_rows = []
    for step_file in ['step4_improved_df.csv', 'step3_df.csv', 'step2_df.csv']:
        filepath = results_path / step_file
        if filepath.exists():
            primer_rows = _load_csv(filepath)
            logger.info(f"Loaded {len(primer_rows)} primers from {step_file}")
            break

    # Convert to PrimerMetrics
    metrics.primers = [PrimerMetrics.from_row(row) for row in primer_rows]
    metrics.primer_count = len(metrics.primers)

    # Calculate derived metrics
    metrics.coverage = _calculate_coverage_metrics(
        metrics.primers, metrics.parameters
    )
    metrics.specificity = _calculate_specificity_metrics(
        metrics.primers, metrics.parameters
    )
    metrics.thermodynamics = _calculate_thermodynamic_metrics(
        metrics.primers, metrics.parameters
    )
    metrics.uniformity = _calculate_uniformity_metrics(metrics.primers)

    # Load optimizer summary for real metrics
    optimizer_summary = _load_optimizer_summary(results_path)
    if optimizer_summary and 'metrics' in optimizer_summary:
        opt_metrics = optimizer_summary['metrics']
        # Override estimated coverage with real optimizer data
        if 'fg_coverage' in opt_metrics:
            metrics.coverage.overall_coverage = opt_metrics['fg_coverage']
            metrics.coverage.from_optimizer = True
            if metrics.coverage.total_bases > 0:
                metrics.coverage.covered_bases = int(
                    opt_metrics['fg_coverage'] * metrics.coverage.total_bases
                )
        # Add gap metrics
        metrics.coverage.mean_gap = opt_metrics.get('mean_gap', 0.0)
        metrics.coverage.max_gap = opt_metrics.get('max_gap', 0.0)
        metrics.coverage.gap_gini = opt_metrics.get('gap_gini', 0.0)
        metrics.coverage.gap_entropy = opt_metrics.get('gap_entropy', 0.0)
        logger.info("Using real optimizer metrics for coverage data")

    # Load validator report so the HTML can surface warnings
    # (per_target_coverage_below_threshold, blacklist_primer_in_set,
    # set_size_mismatch, duplicate_primers, ...) that the optimizer
    # writes to JSON.
    metrics.validation_issues, metrics.validation_ok = _load_validation_report(
        results_path
    )

    # Try to load filtering stats if available
    filter_stats_file = results_path / 'filter_stats.json'
    if filter_stats_file.exists():
        try:
            with open(filter_stats_file, encoding='utf-8') as f:
                stats = json.load(f)
                # Only use fields that FilteringStats accepts
                valid_fields = {
                    'total_kmers', 'after_frequency', 'after_background',
                    'after_gini', 'after_thermodynamic', 'after_complexity',
                    'final_candidates'
                }
                filtered_stats = {k: v for k, v in stats.items() if k in valid_fields}
                metrics.filtering = FilteringStats(**filtered_stats)
        except (json.JSONDecodeError, TypeError, KeyError, OSError) as e:
            logger.warning(f"Failed to load filter stats: {e}")
            # Continue without filtering stats

    return metrics
