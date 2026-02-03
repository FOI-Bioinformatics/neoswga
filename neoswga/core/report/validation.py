"""
Validation utilities for report generation.

Provides validation of results directories and metrics before report generation.
Uses a three-level severity system: ERROR, WARNING, INFO.
"""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

from neoswga.core.report.metrics import PipelineMetrics


class ValidationLevel(Enum):
    """Severity level for validation issues."""
    ERROR = "error"      # Blocks report generation
    WARNING = "warning"  # Should be noted but doesn't block
    INFO = "info"        # Informational message


@dataclass
class ValidationIssue:
    """A single validation issue."""
    level: ValidationLevel
    message: str
    field: str = ""

    def __str__(self) -> str:
        prefix = f"[{self.level.value.upper()}]"
        if self.field:
            return f"{prefix} {self.field}: {self.message}"
        return f"{prefix} {self.message}"


class ValidationResult:
    """Result of validation with issues by severity."""

    def __init__(self, issues: Optional[List[ValidationIssue]] = None):
        self.issues = issues or []

    @property
    def errors(self) -> List[ValidationIssue]:
        """Get all ERROR level issues."""
        return [i for i in self.issues if i.level == ValidationLevel.ERROR]

    @property
    def warnings(self) -> List[ValidationIssue]:
        """Get all WARNING level issues."""
        return [i for i in self.issues if i.level == ValidationLevel.WARNING]

    @property
    def infos(self) -> List[ValidationIssue]:
        """Get all INFO level issues."""
        return [i for i in self.issues if i.level == ValidationLevel.INFO]

    @property
    def is_valid(self) -> bool:
        """True if no errors (warnings/info allowed)."""
        return len(self.errors) == 0

    def add(self, issue: ValidationIssue) -> None:
        """Add an issue to the result."""
        self.issues.append(issue)

    def add_error(self, message: str, field: str = "") -> None:
        """Add an error issue."""
        self.add(ValidationIssue(ValidationLevel.ERROR, message, field))

    def add_warning(self, message: str, field: str = "") -> None:
        """Add a warning issue."""
        self.add(ValidationIssue(ValidationLevel.WARNING, message, field))

    def add_info(self, message: str, field: str = "") -> None:
        """Add an info issue."""
        self.add(ValidationIssue(ValidationLevel.INFO, message, field))


def validate_results_directory(path: str) -> ValidationResult:
    """
    Validate a results directory before report generation.

    Checks:
    - Directory exists
    - Required files present (step4_improved_df.csv)
    - Optional files present (step3_df.csv, step2_df.csv)

    Args:
        path: Path to results directory

    Returns:
        ValidationResult with any issues found
    """
    result = ValidationResult()
    results_path = Path(path)

    # Check directory exists
    if not results_path.exists():
        result.add_error(
            f"Results directory not found: {path}",
            "results_dir"
        )
        return result

    if not results_path.is_dir():
        result.add_error(
            f"Path is not a directory: {path}",
            "results_dir"
        )
        return result

    # Check required files
    required_files = ["step4_improved_df.csv"]
    for fname in required_files:
        if not (results_path / fname).exists():
            result.add_error(
                f"Required file missing: {fname}",
                fname
            )

    # Check optional files
    optional_files = [
        ("step3_df.csv", "Some scoring metrics may be unavailable"),
        ("step2_df.csv", "Filtering statistics may be incomplete"),
        ("filter_stats.json", "Filter funnel data will not be shown"),
        ("params.json", "Parameter information will not be shown"),
    ]
    for fname, note in optional_files:
        # Check in results dir and parent
        if not (results_path / fname).exists():
            parent_file = results_path.parent / fname
            if fname == "params.json" and parent_file.exists():
                continue  # params.json in parent is OK
            result.add_warning(
                f"Optional file missing: {fname} - {note}",
                fname
            )

    return result


def validate_metrics(metrics: PipelineMetrics) -> ValidationResult:
    """
    Validate collected metrics for completeness.

    Checks:
    - Primers list is not empty
    - Target genome size is known
    - Background genome size is known (if applicable)
    - Thermodynamic data is present

    Args:
        metrics: Collected pipeline metrics

    Returns:
        ValidationResult with any issues found
    """
    result = ValidationResult()

    # Check primers
    if not metrics.primers:
        result.add_error(
            "No primers found in results",
            "primers"
        )
    elif len(metrics.primers) < 3:
        result.add_warning(
            f"Only {len(metrics.primers)} primer(s) found - SWGA typically needs 6+",
            "primers"
        )

    # Check genome sizes
    if metrics.target_genome is None or metrics.target_genome.size == 0:
        result.add_warning(
            "Target genome size unknown - coverage estimate unavailable",
            "target_genome"
        )

    if metrics.background_genome is None or metrics.background_genome.size == 0:
        result.add_warning(
            "Background genome size unknown - specificity calculation affected",
            "background_genome"
        )

    # Check thermodynamics
    if metrics.thermodynamics is None:
        result.add_warning(
            "Thermodynamic data missing - Tm metrics unavailable",
            "thermodynamics"
        )
    elif metrics.primers:
        tms = [p.tm for p in metrics.primers if p.tm > 0]
        if not tms:
            result.add_warning(
                "No Tm values found for primers",
                "thermodynamics"
            )

    # Check for potentially problematic values
    if metrics.primers:
        # Check for missing GC content
        missing_gc = sum(1 for p in metrics.primers if p.gc_content == 0)
        if missing_gc > 0:
            result.add_info(
                f"{missing_gc} primer(s) have GC content = 0",
                "primers"
            )

        # Check for very high dimer scores (very negative)
        worst_dimer = min((p.dimer_score for p in metrics.primers), default=0)
        if worst_dimer < -8.0:
            result.add_warning(
                f"High dimer risk detected (dG = {worst_dimer:.1f} kcal/mol)",
                "dimer_risk"
            )

    return result


def format_validation_result(result: ValidationResult) -> str:
    """
    Format validation result for display.

    Args:
        result: Validation result

    Returns:
        Formatted string with all issues
    """
    if not result.issues:
        return "Validation passed: No issues found."

    lines = []

    if result.errors:
        lines.append("ERRORS:")
        for issue in result.errors:
            lines.append(f"  - {issue.message}")

    if result.warnings:
        lines.append("WARNINGS:")
        for issue in result.warnings:
            lines.append(f"  - {issue.message}")

    if result.infos:
        lines.append("INFO:")
        for issue in result.infos:
            lines.append(f"  - {issue.message}")

    if result.is_valid:
        lines.append("")
        lines.append("Validation passed with warnings.")
    else:
        lines.append("")
        lines.append(f"Validation failed: {len(result.errors)} error(s).")

    return "\n".join(lines)
