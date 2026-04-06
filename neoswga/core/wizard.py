"""
Setup wizard for neoswga primer design.

Provides guided configuration for new users by:
1. Analyzing genome composition (GC content, size)
2. Recommending polymerase and reaction conditions
3. Suggesting primer length and additive concentrations
4. Generating validated params.json

Supports two modes:
- Simple: Answer 2-3 questions, get working configuration
- Advanced: Full parameter control with explanations
"""

import json
import logging
import os
import sys
from pathlib import Path
from typing import Optional, Dict, Any, List

from neoswga.core.genome_analysis import (
    calculate_genome_stats,
    get_gc_class,
    recommend_adaptive_qa
)
from neoswga.core.gc_adaptive_strategy import GCAdaptiveStrategy

logger = logging.getLogger(__name__)


# =============================================================================
# Input Validation Helpers
# =============================================================================

def _get_valid_int(prompt: str, default: int, min_val: int, max_val: int,
                   max_attempts: int = 3) -> int:
    """
    Get validated integer input with retry.

    Args:
        prompt: Input prompt
        default: Default value if empty input
        min_val: Minimum allowed value
        max_val: Maximum allowed value
        max_attempts: Maximum retry attempts

    Returns:
        Validated integer value
    """
    for attempt in range(max_attempts):
        response = input(prompt).strip()
        if not response:
            return default
        try:
            value = int(response)
            if min_val <= value <= max_val:
                return value
            print(f"  Please enter a value between {min_val} and {max_val}")
        except ValueError:
            print("  Please enter a valid integer")

    print(f"  Using default value: {default}")
    return default


def _get_valid_float(prompt: str, default: float, min_val: float, max_val: float,
                     max_attempts: int = 3) -> float:
    """
    Get validated float input with retry.

    Args:
        prompt: Input prompt
        default: Default value if empty input
        min_val: Minimum allowed value
        max_val: Maximum allowed value
        max_attempts: Maximum retry attempts

    Returns:
        Validated float value
    """
    for attempt in range(max_attempts):
        response = input(prompt).strip()
        if not response:
            return default
        try:
            value = float(response)
            if min_val <= value <= max_val:
                return value
            print(f"  Please enter a value between {min_val} and {max_val}")
        except ValueError:
            print("  Please enter a valid number")

    print(f"  Using default value: {default}")
    return default


def _get_valid_range(prompt: str, default_min: int, default_max: int,
                     abs_min: int = 4, abs_max: int = 30,
                     max_attempts: int = 3) -> tuple:
    """
    Get validated range input (two integers).

    Args:
        prompt: Input prompt
        default_min: Default minimum
        default_max: Default maximum
        abs_min: Absolute minimum allowed
        abs_max: Absolute maximum allowed
        max_attempts: Maximum retry attempts

    Returns:
        Tuple of (min_val, max_val)
    """
    for attempt in range(max_attempts):
        response = input(prompt).strip()
        if not response:
            return (default_min, default_max)
        parts = response.split()
        if len(parts) != 2:
            print("  Please enter two numbers separated by space (e.g., '10 14')")
            continue
        try:
            val_min, val_max = int(parts[0]), int(parts[1])
            if val_min > val_max:
                print("  First value must be less than or equal to second")
                continue
            if val_min < abs_min or val_max > abs_max:
                print(f"  Values must be between {abs_min} and {abs_max}")
                continue
            return (val_min, val_max)
        except ValueError:
            print("  Please enter valid integers")

    print(f"  Using default range: {default_min}-{default_max}")
    return (default_min, default_max)


def check_jellyfish_available() -> bool:
    """Check if Jellyfish is available in PATH."""
    import shutil
    return shutil.which('jellyfish') is not None


# GC class to human-readable description
GC_CLASS_DESCRIPTIONS = {
    'extreme_at': 'Extreme AT-rich (like Plasmodium)',
    'at_rich': 'AT-rich (like Francisella, Wolbachia)',
    'balanced': 'Balanced (like E. coli)',
    'gc_rich': 'GC-rich (like Mycobacterium, Burkholderia)',
    'extreme_gc': 'Extreme GC-rich (like Streptomyces)'
}

# Polymerase descriptions for users
POLYMERASE_DESCRIPTIONS = {
    'phi29': 'Standard SWGA at 30C with short primers (6-12 bp)',
    'equiphi29': 'Higher specificity at 42C with longer primers (12-18 bp)'
}


def format_bp(length: int) -> str:
    """Format base pair count."""
    if length >= 1_000_000:
        return f"{length/1_000_000:.2f} Mb"
    elif length >= 1000:
        return f"{length/1000:.1f} kb"
    return f"{length} bp"


class SetupWizard:
    """
    Interactive setup wizard for neoswga configuration.

    Usage:
        wizard = SetupWizard()
        wizard.analyze_genome('/path/to/target.fasta')
        wizard.analyze_background('/path/to/background.fasta')  # optional
        config = wizard.generate_config()
        wizard.write_params('params.json')
    """

    def __init__(self, interactive: bool = True, advanced: bool = False):
        """
        Initialize wizard.

        Args:
            interactive: Whether to prompt for user input
            advanced: Whether to show advanced parameters
        """
        self.interactive = interactive
        self.advanced = advanced

        # Genome analysis results
        self.target_stats: Optional[Dict] = None
        self.background_stats: Optional[Dict] = None
        self.target_path: Optional[Path] = None
        self.background_path: Optional[Path] = None
        self.blacklist_paths: List[Path] = []
        self.blacklist_stats: List[Dict] = []

        # Recommendations
        self.gc_class: Optional[str] = None
        self.recommended_polymerase: Optional[str] = None
        self.recommended_kmer_range: Optional[tuple] = None
        self.recommended_additives: Optional[Dict] = None
        self.recommended_temp: Optional[float] = None

        # User overrides
        self.user_overrides: Dict[str, Any] = {}

    def analyze_genome(self, genome_path: str) -> Dict:
        """
        Analyze target genome and generate recommendations.

        Args:
            genome_path: Path to target genome FASTA

        Returns:
            Dictionary with genome analysis and recommendations
        """
        self.target_path = Path(genome_path)

        if not self.target_path.exists():
            raise FileNotFoundError(f"Genome file not found: {genome_path}")

        print(f"\nAnalyzing genome: {self.target_path.name}")
        print("-" * 50)

        # Calculate stats
        self.target_stats = calculate_genome_stats(self.target_path)
        gc = self.target_stats['gc_content']
        length = self.target_stats['length']
        n_contigs = self.target_stats['n_contigs']
        gc_std = self.target_stats.get('gc_std', 0)

        # Classify
        self.gc_class = get_gc_class(gc)

        # Get adaptive QA recommendation
        qa_rec = recommend_adaptive_qa(gc)

        # Use GCAdaptiveStrategy for parameter recommendations
        strategy = GCAdaptiveStrategy(genome_gc_content=gc)
        params = strategy.get_parameters()

        # Store recommendations
        self.recommended_polymerase = params.recommended_polymerase
        self.recommended_kmer_range = params.kmer_range
        self.recommended_temp = params.reaction_temp
        self.recommended_additives = {
            'betaine_m': params.betaine_concentration,
            'dmso_percent': params.dmso_concentration
        }

        # Collect warnings
        self.warnings: List[str] = []

        # Check for potential issues
        if gc < 0.20:
            self.warnings.append(
                f"Extreme AT-rich genome ({gc:.1%} GC) - may have limited primer candidates. "
                "Consider using multi-genome pipeline with closely related species."
            )
        elif gc > 0.75:
            self.warnings.append(
                f"Extreme GC-rich genome ({gc:.1%} GC) - requires high additive concentrations. "
                "Consider betaine >= 2.0M and DMSO >= 5%."
            )

        if length < 500_000:
            self.warnings.append(
                f"Small genome ({format_bp(length)}) - may benefit from longer primers "
                "to reduce non-specific binding."
            )
        elif length > 10_000_000:
            self.warnings.append(
                f"Large genome ({format_bp(length)}) - pipeline may take longer. "
                "Consider using --cpus to parallelize."
            )

        if n_contigs > 100:
            self.warnings.append(
                f"Highly fragmented genome ({n_contigs} contigs) - "
                "primer binding positions may be incomplete at contig boundaries."
            )

        if gc_std > 0.10:
            self.warnings.append(
                f"High GC variability (std={gc_std:.3f}) - genome has regions with "
                "varying GC content. Consider using adaptive GC filtering."
            )

        # Display analysis
        print(f"  Length:     {format_bp(length)}")
        print(f"  Contigs:    {n_contigs}")
        print(f"  GC content: {gc:.1%}")
        print(f"  GC class:   {GC_CLASS_DESCRIPTIONS.get(self.gc_class, self.gc_class)}")
        if gc_std > 0:
            print(f"  GC std:     {gc_std:.4f}")
        print()
        print("Recommendations:")
        print(f"  Polymerase:    {self.recommended_polymerase.upper()}")
        print(f"  Temperature:   {self.recommended_temp}C")
        print(f"  Primer length: {self.recommended_kmer_range[0]}-{self.recommended_kmer_range[1]} bp")
        if self.recommended_additives['betaine_m'] > 0:
            print(f"  Betaine:       {self.recommended_additives['betaine_m']}M")
        if self.recommended_additives['dmso_percent'] > 0:
            print(f"  DMSO:          {self.recommended_additives['dmso_percent']}%")
        print()
        print(f"Adaptive QA: {qa_rec['reason']}")

        # Display warnings
        if self.warnings:
            print()
            print("Warnings:")
            for w in self.warnings:
                print(f"  ! {w}")

        return {
            'stats': self.target_stats,
            'gc_class': self.gc_class,
            'polymerase': self.recommended_polymerase,
            'kmer_range': self.recommended_kmer_range,
            'temperature': self.recommended_temp,
            'additives': self.recommended_additives,
            'qa_recommendation': qa_rec,
            'warnings': self.warnings
        }

    def analyze_background(self, genome_path: str) -> Dict:
        """
        Analyze background genome (optional).

        Args:
            genome_path: Path to background genome FASTA

        Returns:
            Dictionary with background genome stats
        """
        self.background_path = Path(genome_path)

        if not self.background_path.exists():
            raise FileNotFoundError(f"Background genome not found: {genome_path}")

        print(f"\nAnalyzing background: {self.background_path.name}")
        print("-" * 50)

        self.background_stats = calculate_genome_stats(self.background_path)

        print(f"  Length:     {format_bp(self.background_stats['length'])}")
        print(f"  Contigs:    {self.background_stats['n_contigs']}")
        print(f"  GC content: {self.background_stats['gc_content']:.1%}")

        return {'stats': self.background_stats}

    def analyze_blacklist(self, genome_paths: List[str]) -> List[Dict]:
        """Analyze blacklist genome(s).

        Args:
            genome_paths: Paths to blacklist genome FASTA files.

        Returns:
            List of dictionaries with blacklist genome stats.
        """
        results = []
        for genome_path in genome_paths:
            path = Path(genome_path)
            if not path.exists():
                raise FileNotFoundError(f"Blacklist genome not found: {genome_path}")

            print(f"\nAnalyzing blacklist genome: {path.name}")
            print("-" * 50)

            stats = calculate_genome_stats(path)
            self.blacklist_paths.append(path)
            self.blacklist_stats.append(stats)

            print(f"  Length:     {format_bp(stats['length'])}")
            print(f"  Contigs:    {stats['n_contigs']}")
            print(f"  GC content: {stats['gc_content']:.1%}")

            results.append({'stats': stats})

        return results

    def suggest_from_library(self) -> None:
        """Check if genomes match library entries and suggest pre-calculated data."""
        try:
            from neoswga.core.genome_library import GenomeLibrary
            library = GenomeLibrary()
            entries = library.list()
        except Exception:
            return

        if not entries:
            return

        print("\nGenome library check:")
        found = False

        # Check background
        if self.background_path:
            bg_stem = self.background_path.stem
            for entry in entries:
                if entry.name == bg_stem or bg_stem in entry.name:
                    print(f"  [OK] Background '{bg_stem}' found in library as '{entry.name}'")
                    print(f"       Pre-calculated k-mers available (saves computation time)")
                    found = True
                    break

        # Check blacklist
        for bl_path in self.blacklist_paths:
            bl_stem = bl_path.stem
            for entry in entries:
                if entry.name == bl_stem or bl_stem in entry.name:
                    print(f"  [OK] Blacklist '{bl_stem}' found in library as '{entry.name}'")
                    found = True
                    break

        if not found:
            print("  No matching pre-calculated genomes found")
            print("  Tip: Use 'neoswga genome-add <name> <fasta>' to pre-calculate for reuse")

    def prompt_overrides(self) -> None:
        """
        Prompt user for parameter overrides (interactive mode).
        """
        if not self.interactive:
            return

        print("\n" + "=" * 50)
        print("CONFIGURATION OPTIONS")
        print("=" * 50)

        # Polymerase selection
        print(f"\nPolymerase [recommended: {self.recommended_polymerase}]:")
        print("  1. phi29    - " + POLYMERASE_DESCRIPTIONS['phi29'])
        print("  2. equiphi29 - " + POLYMERASE_DESCRIPTIONS['equiphi29'])

        choice = input("Select (1/2) or press Enter for recommended: ").strip()
        if choice == '1':
            self.user_overrides['polymerase'] = 'phi29'
        elif choice == '2':
            self.user_overrides['polymerase'] = 'equiphi29'

        # Primer length with validation
        rec_min, rec_max = self.recommended_kmer_range
        print(f"\nPrimer length range [recommended: {rec_min}-{rec_max} bp]:")
        min_k, max_k = _get_valid_range(
            f"Enter range (e.g., {rec_min} {rec_max}) or press Enter: ",
            default_min=rec_min, default_max=rec_max,
            abs_min=6, abs_max=25
        )
        if min_k != rec_min or max_k != rec_max:
            self.user_overrides['min_k'] = min_k
            self.user_overrides['max_k'] = max_k

        # Advanced options
        if self.advanced:
            self._prompt_advanced_options()

    def _prompt_advanced_options(self) -> None:
        """Prompt for advanced configuration options."""
        print("\n--- Advanced Options ---")

        # Additives with validation
        rec_betaine = self.recommended_additives['betaine_m']
        print(f"\nBetaine concentration [recommended: {rec_betaine}M]:")
        betaine = _get_valid_float(
            "Enter concentration (M) or press Enter: ",
            default=rec_betaine, min_val=0.0, max_val=2.5
        )
        if betaine != rec_betaine:
            self.user_overrides['betaine_m'] = betaine

        rec_dmso = self.recommended_additives['dmso_percent']
        print(f"\nDMSO concentration [recommended: {rec_dmso}%]:")
        dmso = _get_valid_float(
            "Enter concentration (%) or press Enter: ",
            default=rec_dmso, min_val=0.0, max_val=10.0
        )
        if dmso != rec_dmso:
            self.user_overrides['dmso_percent'] = dmso

        # Optimization method
        print("\nOptimization method:")
        print("  1. hybrid (default) - Network + greedy combination")
        print("  2. dominating-set   - Fast graph-based (8x faster)")
        print("  3. background-aware - Clinical applications (10-20x bg reduction)")
        method = input("Select (1/2/3) or press Enter for default: ").strip()
        methods = {'1': 'hybrid', '2': 'dominating-set', '3': 'background-aware'}
        if method in methods:
            self.user_overrides['optimization_method'] = methods[method]

        # Number of primers with validation
        print("\nTarget primer set size [default: 6]:")
        num = _get_valid_int(
            "Enter number or press Enter: ",
            default=6, min_val=1, max_val=50
        )
        if num != 6:
            self.user_overrides['num_primers'] = num

    def generate_config(self, output_dir: str = 'results') -> Dict[str, Any]:
        """
        Generate complete configuration dictionary.

        Args:
            output_dir: Directory for output files

        Returns:
            Complete params.json configuration
        """
        if self.target_stats is None:
            raise ValueError("Must analyze genome first")

        # Start with recommended values
        polymerase = self.user_overrides.get('polymerase', self.recommended_polymerase)
        min_k = self.user_overrides.get('min_k', self.recommended_kmer_range[0])
        max_k = self.user_overrides.get('max_k', self.recommended_kmer_range[1])
        betaine = self.user_overrides.get('betaine_m', self.recommended_additives['betaine_m'])
        dmso = self.user_overrides.get('dmso_percent', self.recommended_additives['dmso_percent'])

        # Calculate GC bounds
        gc = self.target_stats['gc_content']
        gc_tolerance = 0.15
        gc_min = max(0.0, gc - gc_tolerance)
        gc_max = min(1.0, gc + gc_tolerance)

        # Build config
        config = {
            # Schema version for reproducibility
            'schema_version': 1,

            # Genome paths
            'fg_genomes': [str(self.target_path.resolve())],
            'fg_prefixes': [str(Path(output_dir) / self.target_path.stem)],
            'data_dir': output_dir,

            # Genome characteristics
            'genome_gc': round(gc, 4),

            # Polymerase and temperature (match temp to selected polymerase)
            'polymerase': polymerase,
            'reaction_temp': 42.0 if polymerase == 'equiphi29' else (63.0 if polymerase == 'bst' else 30.0),

            # Primer length
            'min_k': min_k,
            'max_k': max_k,

            # GC content filtering
            'gc_tolerance': gc_tolerance,
            'gc_min': round(gc_min, 2),
            'gc_max': round(gc_max, 2),

            # Additives
            'dmso_percent': dmso,
            'betaine_m': betaine,
            'trehalose_m': 0.0,

            # Salt concentrations
            'na_conc': 50.0,
            'mg_conc': 2.0 if self.gc_class in ('extreme_at', 'at_rich') else 0.0,

            # Filtering thresholds
            'min_fg_freq': 1e-5,
            'max_bg_freq': 5e-5,
            'max_gini': 0.6,
            'max_primer': 500,

            # Thermodynamic filters
            'min_tm': 10.0 if polymerase == 'phi29' else 20.0,
            'max_tm': 45.0 if polymerase == 'phi29' else 55.0,
            'max_dimer_bp': 3,
            'max_self_dimer_bp': 4,

            # Optimization
            'num_primers': self.user_overrides.get('num_primers', 6),
            'target_set_size': self.user_overrides.get('num_primers', 6),
            'optimization_method': self.user_overrides.get('optimization_method', 'hybrid'),
            'iterations': 8,
            'max_sets': 5,

            # Execution
            'cpus': 4,
            'verbose': True,

            # Genome topology
            'fg_circular': True,
            'bg_circular': False
        }

        # Add background if provided
        if self.background_path is not None:
            config['bg_genomes'] = [str(self.background_path.resolve())]
            config['bg_prefixes'] = [str(Path(output_dir) / self.background_path.stem)]

        # Add blacklist if provided
        if self.blacklist_paths:
            config['bl_genomes'] = [str(p.resolve()) for p in self.blacklist_paths]
            config['bl_prefixes'] = [
                str(Path(output_dir) / f"bl_{p.stem}") for p in self.blacklist_paths
            ]
            config['bl_penalty'] = 5.0
            config['max_bl_freq'] = 0.0

        return config

    def write_params(self, output_path: str, output_dir: str = 'results') -> Path:
        """
        Write configuration to params.json.

        Args:
            output_path: Path for params.json
            output_dir: Directory for pipeline output files

        Returns:
            Path to written file
        """
        config = self.generate_config(output_dir)
        output_path = Path(output_path)

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w') as f:
            json.dump(config, f, indent=2)

        print(f"\nConfiguration written to: {output_path}")
        return output_path

    def display_summary(self) -> None:
        """Display configuration summary."""
        if self.target_stats is None:
            return

        polymerase = self.user_overrides.get('polymerase', self.recommended_polymerase)
        min_k = self.user_overrides.get('min_k', self.recommended_kmer_range[0])
        max_k = self.user_overrides.get('max_k', self.recommended_kmer_range[1])

        print("\n" + "=" * 50)
        print("CONFIGURATION SUMMARY")
        print("=" * 50)
        print(f"Target genome:    {self.target_path.name}")
        print(f"GC content:       {self.target_stats['gc_content']:.1%} ({self.gc_class})")
        if self.background_path:
            print(f"Background:       {self.background_path.name}")
        if self.blacklist_paths:
            bl_names = ', '.join(p.name for p in self.blacklist_paths)
            print(f"Blacklist:        {bl_names}")
        print()
        print(f"Polymerase:       {polymerase.upper()}")
        print(f"Temperature:      {self.recommended_temp if polymerase == 'equiphi29' else 30.0}C")
        print(f"Primer length:    {min_k}-{max_k} bp")
        print()
        print("Next steps:")
        print("  1. neoswga count-kmers -j params.json")
        print("  2. neoswga filter -j params.json")
        print("  3. neoswga score -j params.json")
        print("  4. neoswga optimize -j params.json")


def run_wizard(
    genome_path: str,
    background_path: Optional[str] = None,
    blacklist_paths: Optional[List[str]] = None,
    output_path: str = 'params.json',
    output_dir: str = 'results',
    interactive: bool = True,
    advanced: bool = False,
    auto_approve: bool = False
) -> Dict[str, Any]:
    """
    Run the setup wizard.

    Args:
        genome_path: Path to target genome FASTA
        background_path: Optional path to background genome
        blacklist_paths: Optional list of paths to blacklist genomes
        output_path: Path for params.json output
        output_dir: Directory for pipeline outputs
        interactive: Whether to prompt for input
        advanced: Whether to show advanced options
        auto_approve: Skip confirmation prompts

    Returns:
        Generated configuration dictionary
    """
    print("\n" + "=" * 50)
    print("NEOSWGA SETUP WIZARD")
    print("=" * 50)

    # Pre-flight checks
    print("\nPre-flight checks:")

    # Check Jellyfish
    if check_jellyfish_available():
        print("  [OK] Jellyfish found in PATH")
    else:
        print("  [!!] Jellyfish NOT found in PATH")
        print("       Install with: conda install -c bioconda jellyfish")
        print("       Or: brew install jellyfish (macOS)")
        if interactive and not auto_approve:
            cont = input("\nContinue anyway? [y/N]: ").strip().lower()
            if cont != 'y':
                print("Aborted. Please install Jellyfish first.")
                raise RuntimeError("Jellyfish not found in PATH")

    # Check genome file
    if Path(genome_path).exists():
        print(f"  [OK] Target genome found: {genome_path}")
    else:
        raise FileNotFoundError(f"Genome file not found: {genome_path}")

    # Check background file
    if background_path:
        if Path(background_path).exists():
            print(f"  [OK] Background genome found: {background_path}")
        else:
            raise FileNotFoundError(f"Background genome not found: {background_path}")

    # Check blacklist files
    if blacklist_paths:
        for bl_path in blacklist_paths:
            if Path(bl_path).exists():
                print(f"  [OK] Blacklist genome found: {bl_path}")
            else:
                raise FileNotFoundError(f"Blacklist genome not found: {bl_path}")

    # Check output path
    output_path_obj = Path(output_path)
    if output_path_obj.exists():
        print(f"  [!!] Output file exists: {output_path}")
        if interactive and not auto_approve:
            overwrite = input("       Overwrite? [y/N]: ").strip().lower()
            if overwrite != 'y':
                print("Aborted. Choose a different output path.")
                return {}

    wizard = SetupWizard(interactive=interactive, advanced=advanced)

    # Analyze genomes
    wizard.analyze_genome(genome_path)

    if background_path:
        wizard.analyze_background(background_path)

    if blacklist_paths:
        wizard.analyze_blacklist(blacklist_paths)

    # Check genome library for suggestions
    wizard.suggest_from_library()

    # Get user input
    if interactive:
        wizard.prompt_overrides()

    # Display summary
    wizard.display_summary()

    # Confirm
    if interactive and not auto_approve:
        print()
        confirm = input("Write configuration? [Y/n]: ").strip().lower()
        if confirm == 'n':
            print("Aborted.")
            return wizard.generate_config(output_dir)

    # Write config
    wizard.write_params(output_path, output_dir)

    # Create output directory
    output_dir_path = Path(output_dir)
    if not output_dir_path.exists():
        output_dir_path.mkdir(parents=True)
        print(f"Created output directory: {output_dir}")

    # Suggest next steps
    print("\n" + "-" * 50)
    print("Setup complete! Next steps:")
    print()
    print(f"  1. Review configuration: cat {output_path}")
    print(f"  2. Validate settings:    neoswga validate-params -j {output_path}")
    print(f"  3. Run pipeline:")
    print(f"       neoswga count-kmers -j {output_path}")
    print(f"       neoswga filter -j {output_path}")
    print(f"       neoswga score -j {output_path}")
    print(f"       neoswga optimize -j {output_path}")
    print(f"  4. Interpret results:    neoswga interpret -d {output_dir}")
    print("-" * 50)

    return wizard.generate_config(output_dir)


if __name__ == '__main__':
    # Example usage
    if len(sys.argv) < 2:
        print("Usage: python wizard.py <genome.fasta> [background.fasta]")
        sys.exit(1)

    genome = sys.argv[1]
    background = sys.argv[2] if len(sys.argv) > 2 else None

    run_wizard(genome, background)
