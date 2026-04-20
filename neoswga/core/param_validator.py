"""
Parameter validation for neoswga configuration.

Validates params.json files for:
1. Required parameters
2. Type correctness
3. Value ranges
4. Parameter interdependencies
5. File existence

Provides three levels of feedback:
- ERROR: Configuration will fail
- WARNING: Suboptimal settings that may affect results
- INFO: Suggestions for improvement
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class ValidationLevel(Enum):
    """Severity level for validation messages."""
    ERROR = "ERROR"
    WARNING = "WARNING"
    INFO = "INFO"


@dataclass
class ValidationMessage:
    """A single validation message."""
    level: ValidationLevel
    parameter: str
    message: str
    suggestion: Optional[str] = None

    def __str__(self):
        result = f"{self.level.value}: [{self.parameter}] {self.message}"
        if self.suggestion:
            result += f"\n    Suggestion: {self.suggestion}"
        return result


# Parameter specifications
REQUIRED_PARAMS = ['fg_genomes', 'fg_prefixes', 'data_dir']

PARAM_RANGES = {
    'min_k': (4, 30),
    'max_k': (4, 30),
    'min_fg_freq': (0.0, 1.0),
    'max_bg_freq': (0.0, 1.0),
    'max_gini': (0.0, 1.0),
    'max_primer': (1, 10000),
    'reaction_temp': (20.0, 70.0),
    'dmso_percent': (0.0, 10.0),
    'betaine_m': (0.0, 2.5),
    'trehalose_m': (0.0, 1.0),
    # Additives below now validated. Upper bounds mirror the authoritative
    # ReactionConditions._validate() checks in reaction_conditions.py:235-264
    # so users get a clear error from ParamValidator before ReactionConditions
    # raises a ValueError deep in the pipeline.
    'formamide_percent': (0.0, 10.0),
    'ethanol_percent': (0.0, 5.0),
    'urea_m': (0.0, 2.0),
    'tmac_m': (0.0, 0.1),
    'glycerol_percent': (0.0, 15.0),
    'peg_percent': (0.0, 15.0),
    'bsa_ug_ml': (0.0, 400.0),
    'na_conc': (0.0, 1000.0),
    'mg_conc': (0.0, 20.0),
    'primer_conc': (1e-9, 1e-4),
    'gc_min': (0.0, 1.0),
    'gc_max': (0.0, 1.0),
    'gc_tolerance': (0.0, 0.5),
    'genome_gc': (0.0, 1.0),
    'min_tm': (0.0, 100.0),
    'max_tm': (0.0, 100.0),
    'num_primers': (1, 50),
    'target_set_size': (1, 50),
    'iterations': (1, 100),
    'max_sets': (1, 100),
    'bl_penalty': (0.0, 100.0),
    'max_bl_freq': (0.0, 1.0),
    'cpus': (1, 128),
}

VALID_POLYMERASES = ['phi29', 'equiphi29', 'bst', 'klenow']
VALID_OPTIMIZATION_METHODS = [
    'hybrid', 'greedy', 'network', 'genetic', 'milp',
    'dominating-set', 'background-aware', 'moea'
]

# Polymerase temperature ranges
POLYMERASE_TEMP_RANGES = {
    'phi29': (20.0, 40.0),
    'equiphi29': (40.0, 47.0),
    'bst': (55.0, 70.0),
    'klenow': (20.0, 42.0),
}

# Recommended primer lengths by polymerase
POLYMERASE_PRIMER_LENGTHS = {
    'phi29': (6, 12),
    'equiphi29': (10, 18),
    'bst': (15, 25),
    'klenow': (8, 15),
}


class ParamValidator:
    """
    Validates neoswga configuration parameters.

    Usage:
        validator = ParamValidator()
        messages = validator.validate_file('params.json')
        for msg in messages:
            print(msg)
    """

    def __init__(self):
        self.messages: List[ValidationMessage] = []

    def validate_file(self, params_path: str) -> List[ValidationMessage]:
        """
        Validate a params.json file.

        Args:
            params_path: Path to params.json

        Returns:
            List of validation messages
        """
        self.messages = []
        path = Path(params_path)

        if not path.exists():
            self.messages.append(ValidationMessage(
                level=ValidationLevel.ERROR,
                parameter='params_file',
                message=f"File not found: {params_path}"
            ))
            return self.messages

        try:
            with open(path) as f:
                params = json.load(f)
        except json.JSONDecodeError as e:
            self.messages.append(ValidationMessage(
                level=ValidationLevel.ERROR,
                parameter='params_file',
                message=f"Invalid JSON: {e}"
            ))
            return self.messages

        return self.validate_params(params)

    def validate_params(self, params: Dict[str, Any]) -> List[ValidationMessage]:
        """
        Validate a parameters dictionary.

        Args:
            params: Parameters dictionary

        Returns:
            List of validation messages
        """
        self.messages = []

        # Check required parameters
        self._check_required(params)

        # Structural validation against the shipped JSON schema. Deliberately
        # kept soft: if jsonschema is not installed we skip (no hard optional
        # dependency) and fall through to the range / interdependency checks.
        self._check_schema(params)

        # Check types and ranges
        self._check_ranges(params)

        # Check file existence
        self._check_files(params)

        # Check parameter interdependencies
        self._check_interdependencies(params)

        # Check for suboptimal settings
        self._check_optimization(params)

        return self.messages

    def _check_required(self, params: Dict) -> None:
        """Check required parameters exist."""
        for param in REQUIRED_PARAMS:
            if param not in params:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.ERROR,
                    parameter=param,
                    message="Required parameter missing"
                ))

    def _check_schema(self, params: Dict) -> None:
        """Validate params against the shipped JSON Schema.

        If jsonschema is not installed this silently no-ops so users don't need
        the extra dependency just to run the pipeline. The schema drift remains
        caught by the range / interdependency checks below.
        """
        try:
            import jsonschema  # type: ignore
        except Exception:
            return
        try:
            from neoswga.core.schema import load_schema
            schema = load_schema()
        except Exception as e:
            logger.debug(f"Skipped schema check (could not load schema): {e}")
            return

        validator = jsonschema.Draft202012Validator(schema)
        for err in validator.iter_errors(params):
            loc = ".".join(str(x) for x in err.absolute_path) or err.validator
            self.messages.append(ValidationMessage(
                level=ValidationLevel.ERROR,
                parameter=str(loc),
                message=f"Schema violation: {err.message}",
            ))

    def _check_ranges(self, params: Dict) -> None:
        """Check parameter values are within valid ranges."""
        for param, (min_val, max_val) in PARAM_RANGES.items():
            if param in params:
                value = params[param]
                if not isinstance(value, (int, float)):
                    self.messages.append(ValidationMessage(
                        level=ValidationLevel.ERROR,
                        parameter=param,
                        message=f"Expected numeric value, got {type(value).__name__}"
                    ))
                elif value < min_val or value > max_val:
                    self.messages.append(ValidationMessage(
                        level=ValidationLevel.ERROR,
                        parameter=param,
                        message=f"Value {value} outside valid range [{min_val}, {max_val}]"
                    ))

        # Check polymerase
        if 'polymerase' in params:
            poly = params['polymerase']
            if poly not in VALID_POLYMERASES:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.ERROR,
                    parameter='polymerase',
                    message=f"Unknown polymerase: {poly}",
                    suggestion=f"Use one of: {', '.join(VALID_POLYMERASES)}"
                ))

        # Check optimization method
        if 'optimization_method' in params:
            method = params['optimization_method']
            if method not in VALID_OPTIMIZATION_METHODS:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.ERROR,
                    parameter='optimization_method',
                    message=f"Unknown optimization method: {method}",
                    suggestion=f"Use one of: {', '.join(VALID_OPTIMIZATION_METHODS)}"
                ))

    def _check_files(self, params: Dict) -> None:
        """Check that genome files exist."""
        for param in ['fg_genomes', 'bg_genomes']:
            if param in params:
                genomes = params[param]
                if not isinstance(genomes, list):
                    genomes = [genomes]
                for genome in genomes:
                    if not Path(genome).exists():
                        self.messages.append(ValidationMessage(
                            level=ValidationLevel.ERROR,
                            parameter=param,
                            message=f"File not found: {genome}"
                        ))

        # Check data directory can be created
        if 'data_dir' in params:
            data_dir = Path(params['data_dir'])
            if data_dir.exists() and not data_dir.is_dir():
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.ERROR,
                    parameter='data_dir',
                    message=f"Path exists but is not a directory: {data_dir}"
                ))

        # Check k-mer files (INFO level - they're created by step1)
        self._check_kmer_files(params)

    def _check_kmer_files(self, params: Dict) -> None:
        """
        Check if k-mer files from step1 exist.

        This is an INFO-level check because k-mer files are created by step1,
        and the user may be validating before running the pipeline.
        """
        fg_prefixes = params.get('fg_prefixes', [])
        bg_prefixes = params.get('bg_prefixes', [])
        min_k = params.get('min_k', 6)
        max_k = params.get('max_k', 12)

        if not fg_prefixes:
            return

        # Check if at least one k-mer file exists per prefix
        missing_fg = []
        missing_bg = []

        for prefix in fg_prefixes:
            has_any = False
            for k in range(min_k, max_k + 1):
                kmer_file = f"{prefix}_{k}mer_all.txt"
                if Path(kmer_file).exists():
                    has_any = True
                    break
            if not has_any:
                missing_fg.append(prefix)

        for prefix in bg_prefixes:
            has_any = False
            for k in range(min_k, max_k + 1):
                kmer_file = f"{prefix}_{k}mer_all.txt"
                if Path(kmer_file).exists():
                    has_any = True
                    break
            if not has_any:
                missing_bg.append(prefix)

        # Report missing k-mer files (one message per category)
        if missing_fg:
            self.messages.append(ValidationMessage(
                level=ValidationLevel.INFO,
                parameter='fg_prefixes',
                message=f"K-mer files not found for {len(missing_fg)} foreground prefix(es)",
                suggestion="Run 'neoswga count-kmers' (Step 1) to generate k-mer files"
            ))

        if missing_bg:
            self.messages.append(ValidationMessage(
                level=ValidationLevel.INFO,
                parameter='bg_prefixes',
                message=f"K-mer files not found for {len(missing_bg)} background prefix(es)",
                suggestion="Run 'neoswga count-kmers' (Step 1) to generate k-mer files"
            ))

    def _check_interdependencies(self, params: Dict) -> None:
        """Check parameter interdependencies and compatibility."""
        polymerase = params.get('polymerase', 'phi29')
        temp = params.get('reaction_temp')
        min_k = params.get('min_k')
        max_k = params.get('max_k')
        betaine = params.get('betaine_m', 0.0)
        dmso = params.get('dmso_percent', 0.0)
        gc_min = params.get('gc_min')
        gc_max = params.get('gc_max')
        genome_gc = params.get('genome_gc')

        # Temperature vs polymerase
        if temp is not None and polymerase in POLYMERASE_TEMP_RANGES:
            min_temp, max_temp = POLYMERASE_TEMP_RANGES[polymerase]
            if temp < min_temp or temp > max_temp:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.WARNING,
                    parameter='reaction_temp',
                    message=f"Temperature {temp}C outside optimal range for {polymerase} ({min_temp}-{max_temp}C)",
                    suggestion=f"Consider using a different polymerase or adjusting temperature"
                ))

        # Primer length vs polymerase
        if min_k is not None and polymerase in POLYMERASE_PRIMER_LENGTHS:
            rec_min, rec_max = POLYMERASE_PRIMER_LENGTHS[polymerase]
            if min_k < rec_min:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.INFO,
                    parameter='min_k',
                    message=f"Primer length {min_k}bp is shorter than typical for {polymerase} ({rec_min}-{rec_max}bp)"
                ))
            if max_k is not None and max_k > rec_max:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.INFO,
                    parameter='max_k',
                    message=f"Primer length {max_k}bp is longer than typical for {polymerase} ({rec_min}-{rec_max}bp)"
                ))

        # Long primers need additives
        if max_k is not None and max_k > 14:
            if betaine < 0.5 and dmso < 3.0:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.WARNING,
                    parameter='max_k',
                    message=f"Long primers ({max_k}bp) work better with additives",
                    suggestion="Consider adding betaine (1.0M) and/or DMSO (5%) for primers >14bp"
                ))

        # Very long primers need more additives
        if max_k is not None and max_k >= 17:
            if betaine < 1.5 or dmso < 3.0:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.WARNING,
                    parameter='max_k',
                    message=f"Very long primers ({max_k}bp) require higher additive concentrations",
                    suggestion="Use betaine >= 1.5M and DMSO >= 3% for 17+ bp primers"
                ))

        # GC bounds vs genome GC
        if gc_min is not None and gc_max is not None and genome_gc is not None:
            if gc_min > genome_gc + 0.1 or gc_max < genome_gc - 0.1:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.WARNING,
                    parameter='gc_min/gc_max',
                    message=f"GC bounds [{gc_min:.2f}-{gc_max:.2f}] may not match genome GC ({genome_gc:.2f})",
                    suggestion="Center GC bounds around genome GC content"
                ))

        # min_k should be <= max_k
        if min_k is not None and max_k is not None and min_k > max_k:
            self.messages.append(ValidationMessage(
                level=ValidationLevel.ERROR,
                parameter='min_k/max_k',
                message=f"min_k ({min_k}) cannot be greater than max_k ({max_k})"
            ))

        # min_tm should be <= max_tm
        min_tm = params.get('min_tm')
        max_tm = params.get('max_tm')
        if min_tm is not None and max_tm is not None and min_tm > max_tm:
            self.messages.append(ValidationMessage(
                level=ValidationLevel.ERROR,
                parameter='min_tm/max_tm',
                message=f"min_tm ({min_tm}) cannot be greater than max_tm ({max_tm})"
            ))

    def _check_optimization(self, params: Dict) -> None:
        """Check for optimization opportunities."""
        polymerase = params.get('polymerase', 'phi29')
        genome_gc = params.get('genome_gc')
        betaine = params.get('betaine_m', 0.0)
        dmso = params.get('dmso_percent', 0.0)

        # High GC genomes benefit from additives
        if genome_gc is not None and genome_gc > 0.60:
            if betaine < 1.5:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.INFO,
                    parameter='betaine_m',
                    message=f"GC-rich genome ({genome_gc:.1%}) may benefit from higher betaine",
                    suggestion="Consider betaine >= 1.5M for GC-rich genomes"
                ))
            if dmso < 3.0:
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.INFO,
                    parameter='dmso_percent',
                    message=f"GC-rich genome ({genome_gc:.1%}) may benefit from DMSO",
                    suggestion="Consider DMSO >= 3% for GC-rich genomes"
                ))

        # AT-rich genomes might use phi29 instead of equiphi29
        if genome_gc is not None and genome_gc < 0.35:
            if polymerase == 'equiphi29':
                self.messages.append(ValidationMessage(
                    level=ValidationLevel.INFO,
                    parameter='polymerase',
                    message=f"AT-rich genome ({genome_gc:.1%}) may work well with phi29",
                    suggestion="phi29 at 30C works well for AT-rich genomes with shorter primers"
                ))


def validate_params_file(params_path: str, verbose: bool = True) -> Tuple[bool, List[ValidationMessage]]:
    """
    Validate a params.json file.

    Args:
        params_path: Path to params.json
        verbose: Print messages to console

    Returns:
        Tuple of (success, messages)
        success is True if no errors, False otherwise
    """
    validator = ParamValidator()
    messages = validator.validate_file(params_path)

    errors = [m for m in messages if m.level == ValidationLevel.ERROR]
    warnings = [m for m in messages if m.level == ValidationLevel.WARNING]
    infos = [m for m in messages if m.level == ValidationLevel.INFO]

    if verbose:
        print("\n" + "=" * 60)
        print("PARAMETER VALIDATION REPORT")
        print("=" * 60)

        if errors:
            print(f"\nErrors ({len(errors)}):")
            for msg in errors:
                print(f"  {msg}")

        if warnings:
            print(f"\nWarnings ({len(warnings)}):")
            for msg in warnings:
                print(f"  {msg}")

        if infos:
            print(f"\nSuggestions ({len(infos)}):")
            for msg in infos:
                print(f"  {msg}")

        print()
        if not errors:
            if not warnings:
                print("Configuration is valid and optimized.")
            else:
                print("Configuration is valid but may be suboptimal.")
        else:
            print(f"Configuration has {len(errors)} error(s) that must be fixed.")

    return len(errors) == 0, messages


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print("Usage: python param_validator.py <params.json>")
        sys.exit(1)

    success, _ = validate_params_file(sys.argv[1])
    sys.exit(0 if success else 1)
