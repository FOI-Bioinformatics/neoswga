"""
Custom exceptions for NeoSWGA.

Provides a hierarchy of exceptions for different error conditions,
enabling proper error handling and informative error messages.
"""

from dataclasses import dataclass
from typing import Any


class NeoSWGAError(Exception):
    """
    Base exception for all NeoSWGA errors.

    All custom exceptions inherit from this class, allowing callers to
    catch all NeoSWGA-specific errors with a single except clause.

    Usage:
        try:
            result = optimizer.optimize(candidates)
        except NeoSWGAError as e:
            logger.error(f"NeoSWGA error: {e}")
    """

    pass


# =============================================================================
# File and I/O Errors
# =============================================================================


class FileError(NeoSWGAError):
    """Base class for file-related errors."""

    pass


class PositionFileNotFoundError(FileError):
    """
    HDF5 position file not found.

    Raised when attempting to load primer positions from a non-existent file.
    """

    def __init__(self, filepath: str, primer_length: int | None = None):
        self.filepath = filepath
        self.primer_length = primer_length
        if primer_length:
            msg = f"Position file not found: {filepath} (k={primer_length})"
        else:
            msg = f"Position file not found: {filepath}"
        super().__init__(msg)


class PositionFileCorruptError(FileError):
    """
    HDF5 position file is corrupted or unreadable.

    Raised when an HDF5 file exists but cannot be read properly.
    """

    def __init__(self, filepath: str, reason: str = ""):
        self.filepath = filepath
        self.reason = reason
        msg = f"Corrupted position file: {filepath}"
        if reason:
            msg += f" ({reason})"
        super().__init__(msg)


class GenomeFileError(FileError):
    """
    Error reading genome FASTA file.
    """

    def __init__(self, filepath: str, reason: str = ""):
        self.filepath = filepath
        self.reason = reason
        msg = f"Cannot read genome file: {filepath}"
        if reason:
            msg += f" ({reason})"
        super().__init__(msg)


class KmerFileError(FileError):
    """
    Error with k-mer count file.
    """

    def __init__(self, filepath: str, reason: str = ""):
        self.filepath = filepath
        self.reason = reason
        msg = f"K-mer file error: {filepath}"
        if reason:
            msg += f" ({reason})"
        super().__init__(msg)


# =============================================================================
# Primer Validation Errors
# =============================================================================


class PrimerError(NeoSWGAError):
    """Base class for primer-related errors."""

    pass


class InvalidPrimerError(PrimerError):
    """
    Invalid primer sequence.

    Raised when a primer contains invalid characters or has invalid properties.
    """

    def __init__(self, primer: str, reason: str):
        self.primer = primer
        self.reason = reason
        super().__init__(f"Invalid primer '{primer}': {reason}")


class PrimerLengthError(PrimerError):
    """
    Primer length outside valid range.
    """

    def __init__(self, primer: str, length: int, min_k: int, max_k: int):
        self.primer = primer
        self.length = length
        self.min_k = min_k
        self.max_k = max_k
        super().__init__(
            f"Primer '{primer}' length {length} outside valid range [{min_k}, {max_k}]"
        )


class PrimerDimerError(PrimerError):
    """
    Primer dimer formation detected.
    """

    def __init__(self, primer1: str, primer2: str, dimer_length: int):
        self.primer1 = primer1
        self.primer2 = primer2
        self.dimer_length = dimer_length
        super().__init__(
            f"Dimer formation between '{primer1}' and '{primer2}' "
            f"({dimer_length} bp complementary)"
        )


class NoCandidatesError(PrimerError):
    """
    No candidate primers available after filtering.

    Raised when all primers are filtered out before optimization.
    """

    def __init__(self, original_count: int, filter_stage: str = ""):
        self.original_count = original_count
        self.filter_stage = filter_stage
        msg = f"No candidates remaining from {original_count} primers"
        if filter_stage:
            msg += f" after {filter_stage}"
        super().__init__(msg)


# =============================================================================
# Optimization Errors
# =============================================================================


class OptimizationError(NeoSWGAError):
    """Base class for optimization-related errors."""

    pass


class OptimizerConvergenceError(OptimizationError):
    """
    Optimizer failed to converge.

    Raised when optimization does not find a satisfactory solution
    within the allowed iterations.
    """

    def __init__(
        self,
        optimizer_name: str,
        iterations: int,
        best_score: float,
        target_score: float | None = None,
    ):
        self.optimizer_name = optimizer_name
        self.iterations = iterations
        self.best_score = best_score
        self.target_score = target_score
        msg = f"Optimizer '{optimizer_name}' failed to converge after {iterations} iterations"
        msg += f" (best score: {best_score:.4f}"
        if target_score is not None:
            msg += f", target: {target_score:.4f}"
        msg += ")"
        super().__init__(msg)


class InsufficientCoverageError(OptimizationError):
    """
    Insufficient genome coverage achieved.

    Raised when the optimized primer set does not meet coverage requirements.
    """

    def __init__(self, achieved_coverage: float, required_coverage: float, num_primers: int):
        self.achieved_coverage = achieved_coverage
        self.required_coverage = required_coverage
        self.num_primers = num_primers
        super().__init__(
            f"Insufficient coverage: {achieved_coverage:.1%} achieved with "
            f"{num_primers} primers (required: {required_coverage:.1%})"
        )


class OptimizerNotFoundError(OptimizationError):
    """
    Unknown optimizer type requested.
    """

    def __init__(self, optimizer_name: str, available: list[str]):
        self.optimizer_name = optimizer_name
        self.available = available
        super().__init__(
            f"Unknown optimizer '{optimizer_name}'. " f"Available: {', '.join(sorted(available))}"
        )


# =============================================================================
# Configuration Errors
# =============================================================================


class ConfigurationError(NeoSWGAError):
    """Base class for configuration-related errors."""

    pass


class InvalidParameterError(ConfigurationError):
    """
    Invalid parameter value in configuration.
    """

    def __init__(self, param_name: str, value: Any, reason: str):
        self.param_name = param_name
        self.value = value
        self.reason = reason
        super().__init__(f"Invalid value for '{param_name}': {value} ({reason})")


class MissingParameterError(ConfigurationError):
    """
    Required parameter missing from configuration.
    """

    def __init__(self, param_name: str, context: str = ""):
        self.param_name = param_name
        self.context = context
        msg = f"Missing required parameter: '{param_name}'"
        if context:
            msg += f" (in {context})"
        super().__init__(msg)


class IncompatibleParametersError(ConfigurationError):
    """
    Incompatible parameter combination.
    """

    def __init__(self, params: list[str], reason: str):
        self.params = params
        self.reason = reason
        super().__init__(f"Incompatible parameters {params}: {reason}")


# =============================================================================
# Thermodynamic Errors
# =============================================================================


class ThermodynamicError(NeoSWGAError):
    """Base class for thermodynamic calculation errors."""

    pass


class TmOutOfRangeError(ThermodynamicError):
    """
    Melting temperature outside valid range.
    """

    def __init__(self, primer: str, tm: float, min_tm: float, max_tm: float):
        self.primer = primer
        self.tm = tm
        self.min_tm = min_tm
        self.max_tm = max_tm
        super().__init__(f"Tm for '{primer}' is {tm:.1f}C, outside range [{min_tm}, {max_tm}]C")


class SecondaryStructureError(ThermodynamicError):
    """
    Problematic secondary structure detected.
    """

    def __init__(self, primer: str, structure_type: str, delta_g: float):
        self.primer = primer
        self.structure_type = structure_type
        self.delta_g = delta_g
        super().__init__(f"Problematic {structure_type} in '{primer}' (dG={delta_g:.1f} kcal/mol)")


# =============================================================================
# Runtime Errors
# =============================================================================


class ResourceError(NeoSWGAError):
    """Base class for resource-related errors."""

    pass


class MemoryLimitError(ResourceError):
    """
    Memory limit exceeded.
    """

    def __init__(self, operation: str, required_mb: float, available_mb: float):
        self.operation = operation
        self.required_mb = required_mb
        self.available_mb = available_mb
        super().__init__(
            f"Memory limit exceeded for {operation}: "
            f"requires {required_mb:.0f}MB, only {available_mb:.0f}MB available"
        )


class OperationTimeoutError(ResourceError):
    """
    Operation timed out.

    Note: Named OperationTimeoutError to avoid shadowing the built-in TimeoutError.
    """

    def __init__(self, operation: str, timeout_seconds: float):
        self.operation = operation
        self.timeout_seconds = timeout_seconds
        super().__init__(f"Operation '{operation}' timed out after {timeout_seconds:.0f} seconds")


class JellyfishError(ResourceError):
    """
    Error running Jellyfish k-mer counter.
    """

    def __init__(self, command: str, return_code: int, stderr: str = ""):
        self.command = command
        self.return_code = return_code
        self.stderr = stderr
        msg = f"Jellyfish command failed (exit code {return_code}): {command}"
        if stderr:
            msg += f"\nStderr: {stderr}"
        super().__init__(msg)


# =============================================================================
# Validation Errors
# =============================================================================


class ValidationError(NeoSWGAError):
    """Base class for validation errors."""

    pass


class DataIntegrityError(ValidationError):
    """
    Data integrity check failed.
    """

    def __init__(self, data_type: str, expected: Any, actual: Any):
        self.data_type = data_type
        self.expected = expected
        self.actual = actual
        super().__init__(f"Data integrity error in {data_type}: expected {expected}, got {actual}")


class PipelineStateError(ValidationError):
    """
    Pipeline state is invalid for requested operation.

    For example, trying to run optimize before score.
    """

    def __init__(self, operation: str, required_state: str, current_state: str):
        self.operation = operation
        self.required_state = required_state
        self.current_state = current_state
        super().__init__(
            f"Cannot run '{operation}': requires '{required_state}' state, "
            f"but current state is '{current_state}'"
        )


# =============================================================================
# Pipeline step prerequisites
# =============================================================================
#
# These lived in `core/pipeline.py` until 2026-09-06. They moved here because
# `cli_unified.py` and `cli/pipeline.py` imported that module at module scope
# for the single purpose of making the exception catchable, and importing it
# pulls in `rf_preprocessing` and therefore scikit-learn: 0.60 s of a 1.14 s
# warm import, paid by `--help` and by every step of a pipeline that no longer
# runs the amplification model. `core/pipeline.py` re-exports both names, so
# existing importers are unaffected.


@dataclass
class StepValidationResult:
    """Result of step prerequisite validation."""

    valid: bool
    missing_files: list[str]
    error_message: str
    remediation: str


class StepPrerequisiteError(NeoSWGAError):
    """Raised when step prerequisites are not met."""

    def __init__(self, step: int, validation: StepValidationResult):
        self.step = step
        self.validation = validation
        message = f"\n{'='*60}\n"
        message += f"STEP {step} PREREQUISITE ERROR\n"
        message += f"{'='*60}\n"
        message += f"\n{validation.error_message}\n"
        if validation.missing_files:
            message += "\nMissing files:\n"
            for f in validation.missing_files[:5]:  # Show first 5
                message += f"  - {f}\n"
            if len(validation.missing_files) > 5:
                message += f"  ... and {len(validation.missing_files) - 5} more\n"
        message += f"\nTo fix this:\n  {validation.remediation}\n"
        message += f"{'='*60}\n"
        super().__init__(message)


# =============================================================================
# Design-path errors
# =============================================================================
#
# Added 2026-09-21 for the valid-design contract. These exist because the
# design path used to answer a failed calculation with a plausible substitute:
# zero free energy for a dimer check that raised, NaN for a Tm that could not
# be computed, a default reach for an unknown polymerase, an empty position
# array for a prefix nobody indexed. Each substitute is the most permissive
# value available, so the run continues and delivers a panel that reads exactly
# like a measured one.
#
# The distinction this family draws is between a MEASUREMENT and its ABSENCE.
# A candidate that misses a Tm window has been measured and rejected; a
# candidate whose Tm raised has not been measured at all. `DesignError` is
# never a QC reason code, and `describe_failure` in `design_result.py` records
# it as a failed run rather than as a screened candidate.
#
# `SearchBudgetExhausted` is deliberately NOT in this family: spending an
# allowance is control flow and a recorded termination reason, not a failure.
#
# Keep this module's imports to `dataclasses`, `typing` and `enum`. It is
# imported at CLI startup and Known Issue 12 is what a heavier import costs.


class DesignError(NeoSWGAError):
    """A design run cannot continue and must not deliver a panel.

    Catch this at the command boundary, write a structured failure record and
    exit nonzero. Do not catch it to try a different model or algorithm: an
    unrequested substitution is what this family exists to prevent.
    """

    #: Never set on this family. Present so that code asking "was this a QC
    #: rejection?" gets a definite no rather than an AttributeError.
    qc_reason = None


class InvalidDesignRequest(DesignError):
    """The requested design cannot be resolved into a valid configuration.

    Unknown or retired settings, contradictory ones, non-finite values, and
    missing required fields. Raised before any search begins.
    """

    def __init__(self, field: str, reason: str, value: Any = None):
        self.field = field
        self.reason = reason
        self.value = value
        detail = f"Invalid design request field '{field}': {reason}"
        if value is not None:
            detail += f" (got {value!r})"
        super().__init__(detail)


class ReferenceDataError(DesignError):
    """A reference answer is missing, stale, corrupt or of an unknown schema.

    Distinct from a verified zero. An index that records "this primer occurs
    nowhere on this reference" is a measurement; an index that was never asked
    the question is not, and returning an empty array for the second is the
    silent-zero shape recorded in Known Issues 5, 6, 13 and 15.
    """

    def __init__(self, artifact: str, reason: str, remediation: str = ""):
        self.artifact = artifact
        self.reason = reason
        self.remediation = remediation
        detail = f"Reference data unusable ({artifact}): {reason}"
        if remediation:
            detail += f"\nTo fix this: {remediation}"
        super().__init__(detail)


class UnsupportedModelError(DesignError):
    """A requested computation falls outside the supported domain of its model.

    An unknown polymerase, an oligo length no parameter set covers, an additive
    combination with no coefficient, a modified base with no terms. The absence
    of an effect model is not a known zero effect.
    """

    def __init__(self, model: str, requested: Any, supported: str = ""):
        self.model = model
        self.requested = requested
        self.supported = supported
        detail = f"Model '{model}' does not support {requested!r}"
        if supported:
            detail += f"; supported domain: {supported}"
        super().__init__(detail)


class ModelEvaluationError(DesignError):
    """A supported computation failed or returned a non-finite value.

    The quantity is named so the failure record can say what could not be
    computed, for which input, rather than reporting that "the command failed".
    """

    def __init__(self, quantity: str, subject: Any, reason: str):
        self.quantity = quantity
        self.subject = subject
        self.reason = reason
        super().__init__(f"Could not compute {quantity} for {subject!r}: {reason}")
