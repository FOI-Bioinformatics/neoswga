"""
Experimental outcome tracking for SWGA primer sets.

Tracks predictions vs actual experimental results to:
1. Validate prediction accuracy
2. Calibrate prediction models
3. Learn from experimental feedback

Usage:
    from neoswga.core.experimental_tracker import ExperimentalTracker

    tracker = ExperimentalTracker('experiment_log.json')

    # Record a prediction
    tracker.record_prediction(
        primer_set=['ATCGATCG', 'GCTAGCTA'],
        prediction=efficiency_prediction,
    )

    # Later, record the outcome
    tracker.record_outcome(
        primer_set=['ATCGATCG', 'GCTAGCTA'],
        actual_enrichment=120.5,
        actual_coverage=0.85,
        notes='Good amplification, some background'
    )

    # Get calibration report
    report = tracker.get_calibration_report()
    print(f"Mean absolute error: {report['mae']:.1f}x")
"""

import json
import logging
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Any
import os

# fcntl is Unix-only; provide fallback for Windows
try:
    import fcntl
    HAS_FCNTL = True
except ImportError:
    HAS_FCNTL = False

logger = logging.getLogger(__name__)


# Forbidden characters in file paths (shell metacharacters, null bytes)
FORBIDDEN_PATH_CHARS = frozenset(';&|`$\x00')


def _validate_log_path(path: str) -> str:
    """
    Validate and sanitize log file path.

    Args:
        path: Path to validate

    Returns:
        Validated absolute path

    Raises:
        ValueError: If path is invalid or suspicious
    """
    # Security checks BEFORE normalization (normpath resolves '..')
    if '..' in path:
        raise ValueError("Invalid log path: directory traversal not allowed")

    if FORBIDDEN_PATH_CHARS.intersection(path):
        raise ValueError("Invalid log path: contains suspicious characters")

    # Ensure it's a JSON file
    if not path.rstrip().lower().endswith('.json'):
        raise ValueError("Invalid log path: must be a .json file")

    # Now safe to normalize
    path = os.path.expanduser(path)
    path = os.path.normpath(path)
    path = os.path.abspath(path)

    # Ensure parent directory exists or can be created
    parent = os.path.dirname(path)
    if parent and not os.path.exists(parent):
        try:
            os.makedirs(parent, exist_ok=True)
        except OSError as e:
            raise ValueError(f"Cannot create log directory: {e}")

    return path


@dataclass
class ExperimentalOutcome:
    """
    Record of an experimental outcome for a primer set.

    Attributes:
        primer_set: List of primer sequences
        predicted_enrichment: Predicted enrichment from model
        predicted_coverage: Predicted coverage from model
        confidence_score: Prediction confidence (0-1)
        actual_enrichment: Measured enrichment (optional)
        actual_coverage: Measured coverage (optional)
        experiment_date: Date of experiment
        notes: Additional notes
        experiment_id: Unique identifier
    """
    primer_set: List[str]
    predicted_enrichment: float
    predicted_coverage: float
    confidence_score: float
    actual_enrichment: Optional[float] = None
    actual_coverage: Optional[float] = None
    experiment_date: Optional[str] = None
    notes: str = ""
    experiment_id: str = ""

    def __post_init__(self):
        if not self.experiment_date:
            self.experiment_date = datetime.now().isoformat()
        if not self.experiment_id:
            self.experiment_id = f"exp_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    @property
    def has_outcome(self) -> bool:
        """Whether experimental outcome has been recorded."""
        return self.actual_enrichment is not None

    @property
    def enrichment_error(self) -> Optional[float]:
        """Absolute error in enrichment prediction."""
        if self.actual_enrichment is None:
            return None
        return abs(self.predicted_enrichment - self.actual_enrichment)

    @property
    def coverage_error(self) -> Optional[float]:
        """Absolute error in coverage prediction."""
        if self.actual_coverage is None:
            return None
        return abs(self.predicted_coverage - self.actual_coverage)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ExperimentalOutcome':
        """Create from dictionary."""
        return cls(**data)


@dataclass
class CalibrationReport:
    """
    Report on prediction calibration based on experimental outcomes.

    Attributes:
        n_experiments: Number of experiments with outcomes
        mae_enrichment: Mean absolute error for enrichment
        mae_coverage: Mean absolute error for coverage
        correlation: Correlation between predicted and actual
        overestimate_rate: Fraction of overestimates
        underestimate_rate: Fraction of underestimates
        recommendation_accuracy: Fraction of correct recommendations
    """
    n_experiments: int
    mae_enrichment: float
    mae_coverage: float
    correlation: float
    overestimate_rate: float
    underestimate_rate: float
    recommendation_accuracy: float
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'n_experiments': self.n_experiments,
            'mae_enrichment': self.mae_enrichment,
            'mae_coverage': self.mae_coverage,
            'correlation': self.correlation,
            'overestimate_rate': self.overestimate_rate,
            'underestimate_rate': self.underestimate_rate,
            'recommendation_accuracy': self.recommendation_accuracy,
            'details': self.details,
        }

    def __str__(self) -> str:
        return f"""Calibration Report ({self.n_experiments} experiments)
  Enrichment MAE: {self.mae_enrichment:.1f}x
  Coverage MAE: {self.mae_coverage:.2%}
  Correlation: {self.correlation:.2f}
  Overestimate rate: {self.overestimate_rate:.0%}
  Underestimate rate: {self.underestimate_rate:.0%}
  Recommendation accuracy: {self.recommendation_accuracy:.0%}"""


class ExperimentalTracker:
    """
    Tracks experimental outcomes for primer sets.

    Maintains a log of predictions and outcomes to:
    - Validate prediction accuracy over time
    - Provide calibration feedback for models
    - Guide iterative primer design
    """

    def __init__(self, log_path: Optional[str] = None):
        """
        Initialize tracker.

        Args:
            log_path: Path to experiment log JSON file.
                     If None, uses ~/.neoswga/experiment_log.json

        Raises:
            ValueError: If log_path is invalid
        """
        if log_path is None:
            log_dir = Path.home() / '.neoswga'
            log_dir.mkdir(exist_ok=True)
            log_path = str(log_dir / 'experiment_log.json')

        # Validate and sanitize path
        self.log_path = _validate_log_path(log_path)
        self.experiments: List[ExperimentalOutcome] = []
        self._load()

    def _load(self) -> None:
        """Load experiments from log file with file locking."""
        if os.path.exists(self.log_path):
            try:
                with open(self.log_path) as f:
                    # Acquire shared lock for reading (Unix only)
                    if HAS_FCNTL:
                        fcntl.flock(f.fileno(), fcntl.LOCK_SH)
                    try:
                        data = json.load(f)
                        self.experiments = [
                            ExperimentalOutcome.from_dict(exp)
                            for exp in data.get('experiments', [])
                        ]
                    finally:
                        if HAS_FCNTL:
                            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
                logger.debug(f"Loaded {len(self.experiments)} experiments")
            except Exception as e:
                logger.warning(f"Failed to load experiment log: {e}")
                self.experiments = []
        else:
            self.experiments = []

    def _save(self) -> None:
        """Save experiments to log file with file locking.

        Uses 'a+' mode to avoid truncation before acquiring lock,
        preventing race conditions when multiple processes write.
        On Windows (no fcntl), falls back to simple write without locking.
        """
        try:
            data = {
                'version': '1.0',
                'updated': datetime.now().isoformat(),
                'experiments': [exp.to_dict() for exp in self.experiments],
            }

            if HAS_FCNTL:
                # Unix: Use 'a+' mode to avoid truncation before lock
                # This prevents race condition where 'w' truncates before lock
                with open(self.log_path, 'a+') as f:
                    # Acquire exclusive lock for writing
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                    try:
                        # Now safe to truncate and write with lock held
                        f.seek(0)
                        f.truncate()
                        json.dump(data, f, indent=2)
                        f.flush()  # Ensure data written before releasing lock
                    finally:
                        fcntl.flock(f.fileno(), fcntl.LOCK_UN)
            else:
                # Windows: No file locking available, use simple write
                with open(self.log_path, 'w') as f:
                    json.dump(data, f, indent=2)

        except Exception as e:
            logger.error(f"Failed to save experiment log: {e}")

    def record_prediction(
        self,
        primer_set: List[str],
        prediction: 'EfficiencyPrediction',
        notes: str = "",
    ) -> str:
        """
        Record a new prediction.

        Args:
            primer_set: Primer sequences
            prediction: EfficiencyPrediction object
            notes: Optional notes

        Returns:
            Experiment ID for later reference
        """
        outcome = ExperimentalOutcome(
            primer_set=primer_set,
            predicted_enrichment=prediction.predicted_enrichment,
            predicted_coverage=prediction.predicted_coverage,
            confidence_score=prediction.confidence_score,
            notes=notes,
        )

        self.experiments.append(outcome)
        self._save()

        logger.info(f"Recorded prediction: {outcome.experiment_id}")
        return outcome.experiment_id

    def record_outcome(
        self,
        primer_set: Optional[List[str]] = None,
        experiment_id: Optional[str] = None,
        actual_enrichment: Optional[float] = None,
        actual_coverage: Optional[float] = None,
        notes: str = "",
    ) -> bool:
        """
        Record experimental outcome for a prediction.

        Args:
            primer_set: Primer sequences (to find matching prediction)
            experiment_id: Experiment ID (alternative to primer_set)
            actual_enrichment: Measured enrichment
            actual_coverage: Measured coverage
            notes: Additional notes

        Returns:
            True if matching experiment found and updated
        """
        # Find matching experiment
        exp = None

        if experiment_id:
            for e in self.experiments:
                if e.experiment_id == experiment_id:
                    exp = e
                    break
        elif primer_set:
            primer_set_sorted = sorted(primer_set)
            for e in self.experiments:
                if sorted(e.primer_set) == primer_set_sorted:
                    exp = e
                    break

        if exp is None:
            logger.warning("No matching experiment found")
            return False

        # Update with outcome
        if actual_enrichment is not None:
            exp.actual_enrichment = actual_enrichment
        if actual_coverage is not None:
            exp.actual_coverage = actual_coverage
        if notes:
            exp.notes = f"{exp.notes}\nOutcome: {notes}".strip()

        self._save()
        logger.info(f"Recorded outcome for {exp.experiment_id}")
        return True

    def get_experiment(self, experiment_id: str) -> Optional[ExperimentalOutcome]:
        """Get experiment by ID."""
        for exp in self.experiments:
            if exp.experiment_id == experiment_id:
                return exp
        return None

    def get_experiments_with_outcomes(self) -> List[ExperimentalOutcome]:
        """Get experiments that have recorded outcomes."""
        return [exp for exp in self.experiments if exp.has_outcome]

    def get_calibration_report(self) -> CalibrationReport:
        """
        Generate calibration report comparing predictions to outcomes.

        Returns:
            CalibrationReport with statistics
        """
        with_outcomes = self.get_experiments_with_outcomes()
        n = len(with_outcomes)

        if n == 0:
            return CalibrationReport(
                n_experiments=0,
                mae_enrichment=0.0,
                mae_coverage=0.0,
                correlation=0.0,
                overestimate_rate=0.0,
                underestimate_rate=0.0,
                recommendation_accuracy=0.0,
            )

        # Calculate errors
        enrichment_errors = []
        coverage_errors = []
        predicted_enrichments = []
        actual_enrichments = []
        overestimates = 0
        underestimates = 0
        correct_recommendations = 0

        for exp in with_outcomes:
            if exp.enrichment_error is not None:
                enrichment_errors.append(exp.enrichment_error)
                predicted_enrichments.append(exp.predicted_enrichment)
                actual_enrichments.append(exp.actual_enrichment)

                if exp.predicted_enrichment > exp.actual_enrichment * 1.2:
                    overestimates += 1
                elif exp.predicted_enrichment < exp.actual_enrichment * 0.8:
                    underestimates += 1

                # Check if recommendation was correct
                # SYNTHESIZE (>0.7 confidence) -> actual enrichment > 50x
                if exp.confidence_score >= 0.7:
                    if exp.actual_enrichment > 50:
                        correct_recommendations += 1
                elif exp.confidence_score >= 0.4:
                    if 10 < exp.actual_enrichment <= 50:
                        correct_recommendations += 1
                else:
                    if exp.actual_enrichment <= 10:
                        correct_recommendations += 1

            if exp.coverage_error is not None:
                coverage_errors.append(exp.coverage_error)

        # Calculate statistics
        import numpy as np

        mae_enrichment = np.mean(enrichment_errors) if enrichment_errors else 0.0
        mae_coverage = np.mean(coverage_errors) if coverage_errors else 0.0

        # Correlation
        if len(predicted_enrichments) >= 2:
            correlation = np.corrcoef(predicted_enrichments, actual_enrichments)[0, 1]
            if np.isnan(correlation):
                correlation = 0.0
        else:
            correlation = 0.0

        return CalibrationReport(
            n_experiments=n,
            mae_enrichment=float(mae_enrichment),
            mae_coverage=float(mae_coverage),
            correlation=float(correlation),
            overestimate_rate=overestimates / n if n > 0 else 0.0,
            underestimate_rate=underestimates / n if n > 0 else 0.0,
            recommendation_accuracy=correct_recommendations / n if n > 0 else 0.0,
            details={
                'enrichment_errors': enrichment_errors,
                'coverage_errors': coverage_errors,
            }
        )

    def list_experiments(self, limit: int = 10) -> List[Dict]:
        """
        List recent experiments.

        Args:
            limit: Maximum number to return

        Returns:
            List of experiment summaries
        """
        results = []
        for exp in sorted(
            self.experiments,
            key=lambda e: e.experiment_date,
            reverse=True
        )[:limit]:
            results.append({
                'id': exp.experiment_id,
                'date': exp.experiment_date,
                'n_primers': len(exp.primer_set),
                'predicted_enrichment': exp.predicted_enrichment,
                'actual_enrichment': exp.actual_enrichment,
                'has_outcome': exp.has_outcome,
            })
        return results

    def clear(self) -> None:
        """Clear all experiments (use with caution)."""
        self.experiments = []
        self._save()
        logger.info("Cleared all experiments")
