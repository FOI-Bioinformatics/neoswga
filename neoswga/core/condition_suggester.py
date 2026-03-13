"""
Reaction condition suggester for neoswga.

Recommends optimal reaction conditions (polymerase, additives, temperature)
based on:
1. Target genome GC content
2. Desired primer length
3. Application context (standard vs clinical)

Provides detailed rationale for each recommendation.
"""

import logging
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from enum import Enum

logger = logging.getLogger(__name__)


class ApplicationContext(Enum):
    """Application context for primer design."""
    STANDARD = "standard"           # General SWGA
    CLINICAL = "clinical"           # Clinical diagnostics (minimize background)
    HIGH_THROUGHPUT = "high_throughput"  # Many samples, speed matters
    LOW_INPUT = "low_input"         # Minimal template DNA


@dataclass
class ConditionRecommendation:
    """Complete condition recommendation with rationale."""
    # Primary settings
    polymerase: str
    reaction_temp: float
    primer_length_range: Tuple[int, int]

    # Additives
    betaine_m: float
    dmso_percent: float
    trehalose_m: float
    mg_conc: float

    # Rationale
    rationale: List[str]
    warnings: List[str]
    alternatives: List[Dict[str, Any]]

    # Confidence
    confidence: float  # 0-1

    # Preset name if applicable
    preset_name: Optional[str] = None


# Decision matrix for polymerase selection
# Format: (gc_min, gc_max, primer_min, primer_max) -> polymerase
POLYMERASE_MATRIX = [
    # AT-rich genomes
    ((0.0, 0.30), (6, 12), 'phi29'),
    ((0.0, 0.30), (12, 18), 'equiphi29'),

    # Moderately AT-rich
    ((0.30, 0.40), (6, 12), 'phi29'),
    ((0.30, 0.40), (12, 18), 'equiphi29'),

    # Balanced genomes
    ((0.40, 0.60), (6, 12), 'phi29'),
    ((0.40, 0.60), (10, 15), 'equiphi29'),
    ((0.40, 0.60), (15, 20), 'equiphi29'),

    # GC-rich genomes
    ((0.60, 0.70), (10, 15), 'equiphi29'),
    ((0.60, 0.70), (15, 20), 'equiphi29'),

    # Extreme GC
    ((0.70, 1.0), (12, 18), 'equiphi29'),
]

# Additive recommendations by GC content and primer length
ADDITIVE_MATRIX = {
    # (gc_class, primer_length_class) -> (betaine_m, dmso_percent, trehalose_m)
    ('extreme_at', 'short'): (0.5, 0.0, 0.0),
    ('extreme_at', 'medium'): (0.8, 3.0, 0.0),
    ('extreme_at', 'long'): (1.0, 5.0, 0.0),

    ('at_rich', 'short'): (0.5, 0.0, 0.0),
    ('at_rich', 'medium'): (1.0, 3.0, 0.0),
    ('at_rich', 'long'): (1.5, 5.0, 0.0),

    ('balanced', 'short'): (0.0, 0.0, 0.0),
    ('balanced', 'medium'): (1.0, 3.0, 0.0),
    ('balanced', 'long'): (1.5, 5.0, 0.0),

    ('gc_rich', 'short'): (1.0, 3.0, 0.0),
    ('gc_rich', 'medium'): (1.5, 5.0, 0.0),
    ('gc_rich', 'long'): (2.0, 5.0, 0.2),

    ('extreme_gc', 'short'): (1.5, 5.0, 0.2),
    ('extreme_gc', 'medium'): (2.0, 5.0, 0.3),
    ('extreme_gc', 'long'): (2.5, 7.0, 0.5),
}

# Mg2+ recommendations
MG_RECOMMENDATIONS = {
    'extreme_at': 2.5,   # Extra Mg2+ stabilizes AT-rich
    'at_rich': 2.0,
    'balanced': 1.5,
    'gc_rich': 1.0,      # Less Mg2+ for GC-rich (already stable)
    'extreme_gc': 0.5,
}


def classify_gc(gc_content: float) -> str:
    """Classify genome by GC content."""
    if gc_content < 0.25:
        return 'extreme_at'
    elif gc_content < 0.40:
        return 'at_rich'
    elif gc_content < 0.60:
        return 'balanced'
    elif gc_content < 0.70:
        return 'gc_rich'
    else:
        return 'extreme_gc'


def classify_primer_length(length: int) -> str:
    """Classify primer length."""
    if length <= 10:
        return 'short'
    elif length <= 14:
        return 'medium'
    else:
        return 'long'


class ConditionSuggester:
    """
    Suggests optimal reaction conditions for SWGA.

    Usage:
        suggester = ConditionSuggester(genome_gc=0.65, primer_length=15)
        recommendation = suggester.suggest()
        suggester.print_recommendation(recommendation)
    """

    def __init__(
        self,
        genome_gc: Optional[float] = None,
        primer_length: Optional[int] = None,
        context: ApplicationContext = ApplicationContext.STANDARD
    ):
        """
        Initialize suggester.

        Args:
            genome_gc: Target genome GC content (0-1)
            primer_length: Target primer length (bp)
            context: Application context
        """
        self.genome_gc = genome_gc
        self.primer_length = primer_length
        self.context = context

        # Derived classifications
        self.gc_class = classify_gc(genome_gc) if genome_gc else 'balanced'
        self.length_class = classify_primer_length(primer_length) if primer_length else 'medium'

    def suggest(self) -> ConditionRecommendation:
        """
        Generate condition recommendation.

        Returns:
            ConditionRecommendation with all settings and rationale
        """
        rationale = []
        warnings = []
        alternatives = []

        # Determine polymerase
        polymerase, poly_rationale = self._select_polymerase()
        rationale.append(poly_rationale)

        # Determine temperature
        if polymerase == 'phi29':
            temp = 30.0
            rationale.append(f"Temperature: 30C (phi29 optimal range)")
        else:
            temp = 42.0
            rationale.append(f"Temperature: 42C (equiphi29 optimal range)")

        # Determine primer length range
        primer_range = self._suggest_primer_range(polymerase)
        rationale.append(f"Primer length: {primer_range[0]}-{primer_range[1]}bp "
                        f"(optimal for {polymerase} with {self.gc_class.replace('_', ' ')} genome)")

        # Determine additives
        additives, add_rationale = self._select_additives()
        rationale.extend(add_rationale)

        # Determine Mg2+
        mg_conc = MG_RECOMMENDATIONS.get(self.gc_class, 1.5)
        if mg_conc != 1.5:
            rationale.append(f"Mg2+: {mg_conc}mM (adjusted for {self.gc_class.replace('_', ' ')} genome)")

        # Add warnings for edge cases
        if self.genome_gc and self.genome_gc < 0.25:
            warnings.append("Extreme AT-rich genome may have limited primer candidates")
        if self.genome_gc and self.genome_gc > 0.70:
            warnings.append("Extreme GC-rich genome requires careful additive optimization")
        if self.primer_length and self.primer_length > 16 and additives[0] < 1.5:
            warnings.append("Long primers (>16bp) typically require betaine >= 1.5M")

        # Context-specific adjustments
        if self.context == ApplicationContext.CLINICAL:
            rationale.append("Clinical context: Prioritizing specificity over coverage")
            alternatives.append({
                'method': 'background-aware',
                'description': 'Use --optimization-method=background-aware for 10-20x background reduction'
            })

        # Suggest alternatives
        if polymerase == 'equiphi29':
            alternatives.append({
                'polymerase': 'phi29',
                'temp': 30.0,
                'description': 'Simpler setup if shorter primers acceptable'
            })
        else:
            alternatives.append({
                'polymerase': 'equiphi29',
                'temp': 42.0,
                'description': 'Higher specificity with longer primers'
            })

        # Determine preset name
        preset_name = self._match_preset(polymerase, additives)

        # Calculate confidence
        confidence = self._calculate_confidence()

        return ConditionRecommendation(
            polymerase=polymerase,
            reaction_temp=temp,
            primer_length_range=primer_range,
            betaine_m=additives[0],
            dmso_percent=additives[1],
            trehalose_m=additives[2],
            mg_conc=mg_conc,
            rationale=rationale,
            warnings=warnings,
            alternatives=alternatives,
            confidence=confidence,
            preset_name=preset_name
        )

    def _select_polymerase(self) -> Tuple[str, str]:
        """Select polymerase with rationale."""
        # Default selections based on GC class
        if self.gc_class in ('extreme_at', 'at_rich'):
            if self.primer_length and self.primer_length > 12:
                return 'equiphi29', f"Polymerase: EquiPhi29 (longer primers for {self.gc_class.replace('_', ' ')} genome)"
            return 'phi29', f"Polymerase: Phi29 (simpler for {self.gc_class.replace('_', ' ')} genome)"

        elif self.gc_class in ('gc_rich', 'extreme_gc'):
            return 'equiphi29', f"Polymerase: EquiPhi29 (better GC tolerance for {self.gc_class.replace('_', ' ')} genome)"

        else:  # balanced
            if self.primer_length and self.primer_length > 14:
                return 'equiphi29', "Polymerase: EquiPhi29 (longer primers benefit from higher temperature)"
            return 'phi29', "Polymerase: Phi29 (standard choice for balanced GC genomes)"

    def _suggest_primer_range(self, polymerase: str) -> Tuple[int, int]:
        """Suggest primer length range."""
        if self.primer_length:
            # Center range around requested length
            margin = 2
            return (max(6, self.primer_length - margin),
                    min(20, self.primer_length + margin))

        # Default ranges by polymerase and GC
        if polymerase == 'phi29':
            if self.gc_class in ('extreme_at', 'at_rich'):
                return (8, 11)
            return (8, 12)
        else:  # equiphi29
            if self.gc_class in ('gc_rich', 'extreme_gc'):
                return (12, 16)
            return (11, 15)

    def _select_additives(self) -> Tuple[Tuple[float, float, float], List[str]]:
        """Select additives with rationale."""
        key = (self.gc_class, self.length_class)
        additives = ADDITIVE_MATRIX.get(key, (1.0, 3.0, 0.0))

        rationale = []

        if additives[0] > 0:
            rationale.append(f"Betaine: {additives[0]}M "
                           f"({'GC equalization' if self.gc_class in ('gc_rich', 'extreme_gc') else 'primer stability'})")

        if additives[1] > 0:
            rationale.append(f"DMSO: {additives[1]}% (secondary structure reduction)")

        if additives[2] > 0:
            rationale.append(f"Trehalose: {additives[2]}M (Tm lowering for GC-rich)")

        if not rationale:
            rationale.append("Additives: None required (balanced conditions)")

        return additives, rationale

    def _match_preset(self, polymerase: str, additives: Tuple[float, float, float]) -> Optional[str]:
        """Match to a named preset if applicable."""
        betaine, dmso, trehalose = additives

        if polymerase == 'phi29' and betaine == 0 and dmso == 0:
            return 'standard_phi29'
        elif polymerase == 'equiphi29':
            if betaine >= 2.0 and dmso >= 5:
                return 'high_gc_genome'
            elif betaine >= 1.0 and dmso >= 5:
                return 'enhanced_equiphi29'
            elif betaine == 0 and dmso == 0:
                return 'equiphi29_baseline'
        return None

    def _calculate_confidence(self) -> float:
        """Calculate recommendation confidence."""
        confidence = 0.85  # Base confidence

        # Higher confidence for balanced genomes
        if self.gc_class == 'balanced':
            confidence += 0.10

        # Lower confidence for extreme genomes
        if self.gc_class in ('extreme_at', 'extreme_gc'):
            confidence -= 0.15

        # Higher confidence when both GC and length provided
        if self.genome_gc and self.primer_length:
            confidence += 0.05

        return min(0.95, max(0.5, confidence))

    def print_recommendation(self, rec: ConditionRecommendation) -> None:
        """Print formatted recommendation."""
        print("\n" + "=" * 60)
        print("REACTION CONDITION RECOMMENDATIONS")
        print("=" * 60)

        if self.genome_gc:
            print(f"\nTarget genome: {self.genome_gc:.1%} GC ({self.gc_class.replace('_', ' ')})")
        if self.primer_length:
            print(f"Target primer length: {self.primer_length}bp")

        print(f"\n--- Recommended Settings ---")
        print(f"  Polymerase:     {rec.polymerase.upper()}")
        print(f"  Temperature:    {rec.reaction_temp}C")
        print(f"  Primer range:   {rec.primer_length_range[0]}-{rec.primer_length_range[1]}bp")
        print(f"  Betaine:        {rec.betaine_m}M")
        print(f"  DMSO:           {rec.dmso_percent}%")
        if rec.trehalose_m > 0:
            print(f"  Trehalose:      {rec.trehalose_m}M")
        print(f"  Mg2+:           {rec.mg_conc}mM")

        if rec.preset_name:
            print(f"\n  Matching preset: {rec.preset_name}")
            print(f"  Use: --preset {rec.preset_name}")

        print(f"\n--- Rationale ---")
        for r in rec.rationale:
            print(f"  - {r}")

        if rec.warnings:
            print(f"\n--- Warnings ---")
            for w in rec.warnings:
                print(f"  ! {w}")

        if rec.alternatives:
            print(f"\n--- Alternatives ---")
            for alt in rec.alternatives:
                if 'polymerase' in alt:
                    print(f"  - {alt['polymerase'].upper()} at {alt['temp']}C: {alt['description']}")
                elif 'method' in alt:
                    print(f"  - {alt['description']}")

        print(f"\nConfidence: {rec.confidence:.0%}")
        print()


def suggest_conditions(
    genome_gc: Optional[float] = None,
    primer_length: Optional[int] = None,
    context: str = 'standard',
    verbose: bool = True
) -> ConditionRecommendation:
    """
    Suggest optimal reaction conditions.

    Args:
        genome_gc: Genome GC content (0-1)
        primer_length: Target primer length
        context: 'standard', 'clinical', 'high_throughput', or 'low_input'
        verbose: Print recommendation

    Returns:
        ConditionRecommendation
    """
    ctx = ApplicationContext(context) if context in [e.value for e in ApplicationContext] else ApplicationContext.STANDARD

    suggester = ConditionSuggester(
        genome_gc=genome_gc,
        primer_length=primer_length,
        context=ctx
    )

    rec = suggester.suggest()

    if verbose:
        suggester.print_recommendation(rec)

    return rec


def sweep_conditions(
    genome_gc: float,
    primer_length: int = 10,
    polymerase: str = 'phi29',
    output_path: Optional[str] = None,
    verbose: bool = True,
) -> List[Dict[str, Any]]:
    """Sweep a grid of reaction conditions and rank by predicted amplification.

    Evaluates combinations of temperature, DMSO, betaine, and Mg2+
    using the mechanistic model to identify conditions that maximize
    predicted amplification for the given genome GC content.

    Args:
        genome_gc: Target genome GC content (0-1).
        primer_length: Representative primer length in bp.
        polymerase: Polymerase type.
        output_path: Optional CSV file path for results.
        verbose: Print results to console.

    Returns:
        List of condition dictionaries sorted by predicted amplification
        (best first).
    """
    from neoswga.core.mechanistic_model import MechanisticModel
    from neoswga.core.reaction_conditions import ReactionConditions
    from neoswga.core.mechanistic_params import get_polymerase_params

    poly_params = get_polymerase_params(polymerase)
    base_temp = poly_params.get('optimal_temp', 30.0)

    # Define sweep grid
    temp_offsets = [-2, 0, 2, 4]
    dmso_values = [0.0, 2.5, 5.0]
    betaine_values = [0.0, 0.5, 1.0]
    mg_values = [1.5, 2.5, 4.0]

    # Representative primer with mixed GC content matching target genome
    # Using alternating bases provides a more realistic Tm estimate
    # than homopolymer sequences
    gc_bases = 'GC'
    at_bases = 'AT'
    gc_count = int(primer_length * genome_gc)
    at_count = primer_length - gc_count
    sample_primer = (gc_bases * gc_count + at_bases * at_count)[:primer_length]

    results = []
    for t_off in temp_offsets:
        temp = base_temp + t_off
        for dmso in dmso_values:
            for betaine in betaine_values:
                for mg in mg_values:
                    try:
                        conditions = ReactionConditions(
                            temp=temp,
                            polymerase=polymerase,
                            mg_conc=mg,
                            dmso_percent=dmso,
                            betaine_m=betaine,
                        )
                        model = MechanisticModel(conditions)
                        effects = model.calculate_effects(sample_primer, genome_gc)

                        results.append({
                            'temperature': temp,
                            'dmso_percent': dmso,
                            'betaine_m': betaine,
                            'mg_conc': mg,
                            'amplification_factor': effects.predicted_amplification_factor,
                            'processivity_factor': effects.processivity_factor,
                            'accessibility_factor': effects.accessibility_factor,
                            'effective_tm': effects.effective_tm,
                        })
                    except Exception:
                        continue

    # Sort by amplification factor (best first)
    results.sort(key=lambda x: x['amplification_factor'], reverse=True)

    if verbose and results:
        print(f"\nCondition sweep for {polymerase} (GC={genome_gc:.1%}, primer={primer_length}bp)")
        print("=" * 80)
        print(f"{'Rank':<5} {'Temp':>5} {'DMSO%':>6} {'Bet.M':>6} {'Mg2+':>5} "
              f"{'Ampl.':>7} {'Proc.':>6} {'Access':>7} {'Eff.Tm':>7}")
        print("-" * 80)
        for i, r in enumerate(results[:15], 1):
            print(f"{i:<5} {r['temperature']:>5.1f} {r['dmso_percent']:>6.1f} "
                  f"{r['betaine_m']:>6.1f} {r['mg_conc']:>5.1f} "
                  f"{r['amplification_factor']:>7.3f} {r['processivity_factor']:>6.2f} "
                  f"{r['accessibility_factor']:>7.2f} {r['effective_tm']:>7.1f}")
        if len(results) > 15:
            print(f"... ({len(results)} total combinations evaluated)")

    if output_path and results:
        import csv
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        if verbose:
            print(f"\nFull results saved to {output_path}")

    return results


if __name__ == '__main__':
    import sys

    # Example usage
    if len(sys.argv) >= 2:
        gc = float(sys.argv[1])
        length = int(sys.argv[2]) if len(sys.argv) >= 3 else None
        suggest_conditions(genome_gc=gc, primer_length=length)
    else:
        print("Usage: python condition_suggester.py <gc_content> [primer_length]")
        print("\nExamples:")
        print("  python condition_suggester.py 0.32        # AT-rich genome")
        print("  python condition_suggester.py 0.65 16     # GC-rich, 16bp primers")
        print("  python condition_suggester.py 0.50 12     # Balanced, 12bp primers")
