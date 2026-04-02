"""
Base optimizer interface for SWGA primer set optimization.

All optimizer implementations should inherit from BaseOptimizer to ensure
consistent interfaces, enable polymorphism, and simplify testing.

Design principles:
- Strategy pattern: swap optimizers without changing client code
- Dependency injection: pass caches and conditions, don't create internally
- Immutable results: OptimizationResult is a frozen dataclass
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set, Tuple, Any
from enum import Enum
import numpy as np
import logging

logger = logging.getLogger(__name__)


class OptimizationStatus(Enum):
    """Status of optimization run."""
    SUCCESS = "success"
    PARTIAL = "partial"  # Converged but below target
    NO_CONVERGENCE = "no_convergence"
    ERROR = "error"


@dataclass(frozen=True)
class PrimerSetMetrics:
    """
    Metrics for evaluating a primer set.

    Immutable to prevent accidental modification after optimization.
    """
    # Coverage metrics
    fg_coverage: float  # Fraction of foreground genome covered
    bg_coverage: float  # Fraction of background genome covered (lower is better)
    coverage_uniformity: float  # Gini coefficient of gap sizes (lower is better)

    # Binding metrics
    total_fg_sites: int  # Total binding sites in foreground
    total_bg_sites: int  # Total binding sites in background
    selectivity_ratio: float  # fg_sites / bg_sites (higher is better)

    # Thermodynamic metrics
    mean_tm: float  # Mean melting temperature
    tm_range: Tuple[float, float]  # (min_tm, max_tm)
    dimer_risk_score: float  # 0-1, fraction of primer pairs with dimer risk

    # Gap statistics
    mean_gap: float  # Mean gap between binding sites
    max_gap: float  # Maximum gap (coverage hole)
    gap_gini: float  # Gini coefficient of gaps
    gap_entropy: float  # Shannon entropy of gap distribution (bits)

    # Strand alternation metrics
    strand_alternation_score: float  # Fraction of adjacent pairs alternating strands (0-1)
    strand_coverage_ratio: float  # Balance between forward and reverse strand sites (0-1)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'fg_coverage': self.fg_coverage,
            'bg_coverage': self.bg_coverage,
            'coverage_uniformity': self.coverage_uniformity,
            'total_fg_sites': self.total_fg_sites,
            'total_bg_sites': self.total_bg_sites,
            'selectivity_ratio': self.selectivity_ratio,
            'mean_tm': self.mean_tm,
            'tm_range': list(self.tm_range),
            'dimer_risk_score': self.dimer_risk_score,
            'mean_gap': self.mean_gap,
            'max_gap': self.max_gap,
            'gap_gini': self.gap_gini,
            'gap_entropy': self.gap_entropy,
            'strand_alternation_score': self.strand_alternation_score,
            'strand_coverage_ratio': self.strand_coverage_ratio,
        }

    def normalized_score(
        self,
        coverage_w: float = 0.35,
        selectivity_w: float = 0.30,
        dimer_w: float = 0.15,
        evenness_w: float = 0.10,
        tm_w: float = 0.10,
    ) -> float:
        """Compute a [0,1] composite score from metrics.

        Provides a single comparable value regardless of which optimizer
        produced the result.  All component metrics are clamped to [0,1]
        before weighting.

        Args:
            coverage_w: Weight for foreground coverage (default 0.35).
            selectivity_w: Weight for selectivity (default 0.30).
            dimer_w: Weight for dimer safety (default 0.15).
            evenness_w: Weight for gap evenness (default 0.10).
            tm_w: Weight for Tm tightness (default 0.10).
        """
        # Coverage: already [0,1]
        cov = min(self.fg_coverage, 1.0)

        # Selectivity: cap ratio at 100, normalize
        sel = min(self.selectivity_ratio / 100.0, 1.0) if self.selectivity_ratio > 0 else 0.0

        # Dimer safety: 1 - risk (lower risk is better)
        dimer_safe = 1.0 - min(self.dimer_risk_score, 1.0)

        # Evenness: 1 - Gini (lower Gini is more uniform)
        even = 1.0 - min(self.gap_gini, 1.0)

        # Tm tightness: narrow Tm spread is better
        tm_spread = abs(self.tm_range[1] - self.tm_range[0]) if self.tm_range else 0.0
        tm_tight = max(0.0, 1.0 - tm_spread / 20.0)  # 0C spread -> 1.0, 20C+ -> 0.0

        return (
            coverage_w * cov
            + selectivity_w * sel
            + dimer_w * dimer_safe
            + evenness_w * even
            + tm_w * tm_tight
        )

    @classmethod
    def empty(cls) -> 'PrimerSetMetrics':
        """Create empty metrics for failed optimization."""
        return cls(
            fg_coverage=0.0,
            bg_coverage=0.0,
            coverage_uniformity=1.0,
            total_fg_sites=0,
            total_bg_sites=0,
            selectivity_ratio=0.0,
            mean_tm=0.0,
            tm_range=(0.0, 0.0),
            dimer_risk_score=1.0,
            mean_gap=float('inf'),
            max_gap=float('inf'),
            gap_gini=1.0,
            gap_entropy=0.0,
            strand_alternation_score=0.0,
            strand_coverage_ratio=0.0,
        )


@dataclass(frozen=True)
class OptimizationResult:
    """
    Result of primer set optimization.

    Immutable to ensure results are not accidentally modified.
    Contains the selected primers, score, and detailed metrics.
    """
    primers: Tuple[str, ...]  # Selected primer sequences (tuple for immutability)
    score: float  # Overall optimization score
    status: OptimizationStatus
    metrics: PrimerSetMetrics
    iterations: int  # Number of iterations/generations used
    optimizer_name: str  # Name of optimizer that produced this result

    # Optional detailed results
    all_scores: Optional[Tuple[float, ...]] = None  # Score history
    message: str = ""  # Status message or error description

    @property
    def num_primers(self) -> int:
        """Number of primers in the set."""
        return len(self.primers)

    @property
    def is_success(self) -> bool:
        """Whether optimization succeeded."""
        return self.status == OptimizationStatus.SUCCESS

    @property
    def normalized_score(self) -> float:
        """Comparable [0,1] score derived from metrics.

        Unlike ``score`` (which varies in scale per optimizer), this value
        is always in [0,1] and comparable across different optimizer types.
        Returns 0.0 only for ERROR or NO_CONVERGENCE status; PARTIAL
        results with valid metrics are scored normally.
        """
        if self.status in (OptimizationStatus.ERROR, OptimizationStatus.NO_CONVERGENCE):
            return 0.0
        return self.metrics.normalized_score()

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            'primers': list(self.primers),
            'score': self.score,
            'status': self.status.value,
            'metrics': self.metrics.to_dict(),
            'iterations': self.iterations,
            'optimizer_name': self.optimizer_name,
            'num_primers': self.num_primers,
            'message': self.message,
        }

    @classmethod
    def failure(cls, optimizer_name: str, message: str) -> 'OptimizationResult':
        """Create a failure result."""
        return cls(
            primers=tuple(),
            score=float('-inf'),
            status=OptimizationStatus.ERROR,
            metrics=PrimerSetMetrics.empty(),
            iterations=0,
            optimizer_name=optimizer_name,
            message=message,
        )


@dataclass
class OptimizerConfig:
    """
    Base configuration for optimizers.

    Subclasses can extend with optimizer-specific parameters.
    """
    target_set_size: int = 6
    max_iterations: int = 100
    max_dimer_bp: int = 4
    min_tm: float = 20.0
    max_tm: float = 50.0
    verbose: bool = True

    # Polymerase-specific extension reach (bp) for coverage computation.
    # Default 70000 (Phi29). Set from POLYMERASE_CHARACTERISTICS.
    extension_reach: int = 70000

    # Convergence criteria
    convergence_threshold: float = 0.001
    patience: int = 10  # Iterations without improvement before stopping

    def validate(self) -> None:
        """Validate configuration parameters."""
        if self.target_set_size < 1:
            raise ValueError(f"target_set_size must be >= 1, got {self.target_set_size}")
        if self.max_iterations < 1:
            raise ValueError(f"max_iterations must be >= 1, got {self.max_iterations}")
        if self.min_tm >= self.max_tm:
            raise ValueError(f"min_tm ({self.min_tm}) must be < max_tm ({self.max_tm})")


class BaseOptimizer(ABC):
    """
    Abstract base class for primer set optimizers.

    All optimizer implementations must inherit from this class and implement
    the optimize() method. This ensures consistent interfaces across different
    optimization strategies (greedy, genetic, network-based, etc.).

    Usage:
        optimizer = SomeOptimizer(cache, fg_prefixes, fg_lengths, ...)
        result = optimizer.optimize(candidates, target_size=10)

        if result.is_success:
            print(f"Selected primers: {result.primers}")
            print(f"Coverage: {result.metrics.fg_coverage:.1%}")

    Subclassing:
        class MyOptimizer(BaseOptimizer):
            @property
            def name(self) -> str:
                return "my-optimizer"

            def optimize(self, candidates, target_size=10, **kwargs):
                # Implementation
                return OptimizationResult(...)
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
    ):
        """
        Initialize optimizer with position data.

        Args:
            position_cache: PositionCache instance for primer position lookups
            fg_prefixes: Foreground (target) genome HDF5 prefixes
            fg_seq_lengths: Foreground genome lengths
            bg_prefixes: Background (off-target) genome HDF5 prefixes (optional)
            bg_seq_lengths: Background genome lengths (optional)
            config: Optimizer configuration
        """
        self.cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.config = config or OptimizerConfig()

        # Computed properties
        self.fg_total_length = sum(fg_seq_lengths)
        self.bg_total_length = sum(bg_seq_lengths) if bg_seq_lengths else 0

        # Validate
        self.config.validate()
        self._validate_inputs()

    def _validate_inputs(self) -> None:
        """Validate constructor inputs."""
        if len(self.fg_prefixes) != len(self.fg_seq_lengths):
            raise ValueError(
                f"fg_prefixes ({len(self.fg_prefixes)}) and fg_seq_lengths "
                f"({len(self.fg_seq_lengths)}) must have same length"
            )
        if self.bg_prefixes and len(self.bg_prefixes) != len(self.bg_seq_lengths):
            raise ValueError(
                f"bg_prefixes ({len(self.bg_prefixes)}) and bg_seq_lengths "
                f"({len(self.bg_seq_lengths)}) must have same length"
            )

    @property
    @abstractmethod
    def name(self) -> str:
        """
        Human-readable optimizer name.

        Used for logging, result attribution, and factory registration.
        Should be kebab-case (e.g., 'greedy-bfs', 'dominating-set').
        """
        pass

    @property
    def description(self) -> str:
        """
        One-line description of the optimizer.

        Override in subclasses to provide helpful description.
        """
        return f"{self.name} optimizer"

    @property
    def supports_background(self) -> bool:
        """
        Whether this optimizer uses background genome data.

        Some optimizers (e.g., dominating-set) only consider foreground
        coverage and don't use background data for optimization.
        """
        return True

    @abstractmethod
    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs
    ) -> OptimizationResult:
        """
        Find optimal primer set from candidates.

        This is the main entry point for optimization. Implementations should:
        1. Validate inputs
        2. Run optimization algorithm
        3. Compute metrics for result
        4. Return OptimizationResult

        Args:
            candidates: Pool of candidate primer sequences
            target_size: Desired number of primers (overrides config if provided)
            fixed_primers: Optional list of primers that must be included in the
                result. These primers are pre-selected and the optimizer will
                select additional primers from candidates to complement them.
                Useful for iterative design where some primers have already been
                validated experimentally.
            **kwargs: Optimizer-specific parameters

        Returns:
            OptimizationResult with selected primers and metrics

        Raises:
            ValueError: If candidates is empty or invalid
        """
        pass

    def _validate_candidates(self, candidates: List[str]) -> List[str]:
        """
        Validate and deduplicate candidate primers.

        Args:
            candidates: Raw candidate list

        Returns:
            Validated, deduplicated candidate list

        Raises:
            ValueError: If no valid candidates
        """
        if not candidates:
            raise ValueError("candidates list cannot be empty")

        # Deduplicate while preserving order
        seen = set()
        valid = []
        for primer in candidates:
            if primer not in seen:
                seen.add(primer)
                # Basic validation
                if primer and all(c in 'ATCG' for c in primer.upper()):
                    valid.append(primer)
                else:
                    logger.warning(f"Skipping invalid primer: {primer}")

        if not valid:
            raise ValueError("No valid primers in candidates list")

        return valid

    def get_primer_positions(
        self,
        primer: str,
        prefix: str,
        strand: str = 'both'
    ) -> np.ndarray:
        """
        Get binding positions for a primer in a genome.

        Wrapper around position cache with consistent interface.

        Args:
            primer: Primer sequence
            prefix: Genome HDF5 prefix
            strand: 'forward', 'reverse', or 'both'

        Returns:
            Array of binding positions (empty array if cache unavailable)
        """
        if self.cache is None:
            logger.warning("Position cache not available")
            return np.array([], dtype=np.int64)

        try:
            return self.cache.get_positions(prefix, primer, strand)
        except Exception as e:
            logger.warning(f"Failed to get positions for {primer}: {e}")
            return np.array([], dtype=np.int64)

    def compute_metrics(
        self,
        primers: List[str],
    ) -> PrimerSetMetrics:
        """
        Compute detailed metrics for a primer set.

        Args:
            primers: List of primer sequences

        Returns:
            PrimerSetMetrics with coverage, binding, and thermodynamic stats
        """
        if not primers:
            return PrimerSetMetrics.empty()

        # Collect all positions
        fg_positions = []
        bg_positions = []

        for primer in primers:
            for prefix in self.fg_prefixes:
                positions = self.get_primer_positions(primer, prefix, 'both')
                fg_positions.extend(positions.tolist())

            for prefix in self.bg_prefixes:
                positions = self.get_primer_positions(primer, prefix, 'both')
                bg_positions.extend(positions.tolist())

        fg_positions = sorted(set(fg_positions))
        bg_positions = sorted(set(bg_positions))

        # Coverage
        fg_coverage = self._compute_coverage(fg_positions, self.fg_total_length)
        bg_coverage = self._compute_coverage(bg_positions, self.bg_total_length) if self.bg_total_length > 0 else 0.0

        # Gap statistics
        gaps = self._compute_gaps(fg_positions, self.fg_total_length)
        mean_gap = np.mean(gaps) if gaps else float('inf')
        max_gap = max(gaps) if gaps else float('inf')
        gap_gini = self._gini(gaps) if gaps else 1.0

        # Shannon entropy of gap distribution
        if len(gaps) >= 2:
            from scipy.stats import entropy as scipy_entropy
            counts, _ = np.histogram(gaps, bins=min(50, len(gaps)))
            counts = counts[counts > 0]
            probs = counts / counts.sum()
            gap_entropy = float(scipy_entropy(probs, base=2))
        else:
            gap_entropy = 0.0

        # Selectivity
        total_fg = len(fg_positions)
        total_bg = len(bg_positions)
        selectivity = total_fg / max(total_bg, 1)

        # Tm calculation (simplified)
        tms = [self._estimate_tm(p) for p in primers]
        mean_tm = np.mean(tms)
        tm_range = (min(tms), max(tms))

        # Dimer risk (simplified - would use dimer module in full implementation)
        dimer_risk = self._estimate_dimer_risk(primers)

        # Strand alternation metrics
        strand_alt_score = 0.0
        strand_cov_ratio = 0.0
        if self.cache is not None and hasattr(self.cache, 'compute_strand_alternation_stats'):
            for prefix, length in zip(self.fg_prefixes, self.fg_seq_lengths):
                try:
                    strand_stats = self.cache.compute_strand_alternation_stats(
                        prefix, primers, length
                    )
                    strand_alt_score = strand_stats['strand_alternation_score']
                    strand_cov_ratio = strand_stats['strand_coverage_ratio']
                    break  # Use first fg genome stats
                except Exception:
                    pass

        return PrimerSetMetrics(
            fg_coverage=fg_coverage,
            bg_coverage=bg_coverage,
            coverage_uniformity=gap_gini,
            total_fg_sites=total_fg,
            total_bg_sites=total_bg,
            selectivity_ratio=selectivity,
            mean_tm=mean_tm,
            tm_range=tm_range,
            dimer_risk_score=dimer_risk,
            mean_gap=mean_gap,
            max_gap=max_gap,
            gap_gini=gap_gini,
            gap_entropy=gap_entropy,
            strand_alternation_score=strand_alt_score,
            strand_coverage_ratio=strand_cov_ratio,
        )

    def _compute_coverage(self, positions: List[int], total_length: int) -> float:
        """Compute coverage fraction from positions.

        Uses ``config.extension_reach`` (polymerase processivity in bp) to
        determine how far each binding site extends along the genome.
        Phi29 extends in both directions from the primer binding site.
        """
        if not positions or total_length == 0:
            return 0.0
        extension_reach = self.config.extension_reach

        # Build intervals: each binding site covers [pos - reach, pos + reach]
        intervals = []
        for pos in sorted(positions):
            start = max(0, pos - extension_reach)
            end = min(total_length, pos + extension_reach)
            intervals.append((start, end))

        # Merge overlapping intervals
        merged = [intervals[0]]
        for start, end in intervals[1:]:
            if start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))

        covered = sum(end - start for start, end in merged)
        return min(1.0, covered / total_length)

    def _compute_gaps(self, positions: List[int], total_length: int) -> List[float]:
        """Compute gaps between adjacent binding sites."""
        if len(positions) < 2:
            return [float(total_length)]
        positions = sorted(positions)
        gaps = []
        for i in range(1, len(positions)):
            gaps.append(positions[i] - positions[i-1])
        # Add wrap-around gap for circular genomes
        gaps.append(total_length - positions[-1] + positions[0])
        return gaps

    def _gini(self, values: List[float]) -> float:
        """Compute Gini coefficient."""
        if not values:
            return 1.0
        values = sorted(values)
        n = len(values)
        total = sum(values)
        if total == 0:
            return 0.0
        cumsum = 0
        gini_sum = 0
        for i, v in enumerate(values):
            cumsum += v
            gini_sum += (2 * (i + 1) - n - 1) * v
        return gini_sum / (n * total)

    def _estimate_tm(self, primer: str) -> float:
        """Estimate melting temperature using Wallace rule."""
        gc = primer.upper().count('G') + primer.upper().count('C')
        at = primer.upper().count('A') + primer.upper().count('T')
        return 4 * gc + 2 * at

    def _estimate_dimer_risk(self, primers: List[str]) -> float:
        """
        Estimate dimer risk for primer set.

        Returns fraction of primer pairs with potential dimer formation.
        Uses the centralized DimerValidator for consistent results across
        all optimizers.
        """
        if len(primers) < 2:
            return 0.0

        try:
            from .dimer_validator import DimerValidator
            validator = DimerValidator(max_dimer_bp=self.config.max_dimer_bp)
            return validator.dimer_risk(primers)
        except ImportError:
            logger.warning("Dimer module not available, using estimate")
            return 0.1  # Fallback estimate


class CompositeOptimizer(BaseOptimizer):
    """
    Optimizer that combines results from multiple sub-optimizers.

    Useful for ensemble approaches or multi-stage optimization.
    """

    def __init__(
        self,
        optimizers: List[BaseOptimizer],
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config
        )
        self.optimizers = optimizers

    @property
    def name(self) -> str:
        return "composite"

    @property
    def description(self) -> str:
        names = [opt.name for opt in self.optimizers]
        return f"Composite optimizer combining: {', '.join(names)}"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        fixed_primers: Optional[List[str]] = None,
        **kwargs
    ) -> OptimizationResult:
        """Run all sub-optimizers and return best result."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        best_result = None

        for optimizer in self.optimizers:
            try:
                result = optimizer.optimize(
                    candidates, target, fixed_primers=fixed_primers, **kwargs
                )
                if best_result is None or result.score > best_result.score:
                    best_result = result
            except Exception as e:
                logger.warning(f"Optimizer {optimizer.name} failed: {e}")

        if best_result is None:
            return OptimizationResult.failure(self.name, "All sub-optimizers failed")

        return best_result
