"""
Multi-additive optimization for SWGA reaction conditions.

Finds optimal additive concentrations (DMSO, betaine, trehalose, etc.)
based on primer characteristics, template properties, and application goals.

This module combines:
1. Mechanistic model predictions
2. Literature-based constraints
3. Grid search optimization
4. Rule-based recommendations

Usage:
    from neoswga.core.additive_optimizer import AdditiveOptimizer

    optimizer = AdditiveOptimizer(polymerase='equiphi29')
    recommendation = optimizer.optimize(
        primer_length=15,
        template_gc=0.65,
        optimize_for='amplification'
    )
    print(recommendation.summary())

Literature basis:
- Musso et al. (2006): SWGA additive optimization
- Henke et al. (1997): Betaine effects
- Varadharajan et al. (2017): EquiPhi29 conditions
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import itertools
import math

from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
from neoswga.core.mechanistic_params import get_polymerase_params


@dataclass
class AdditiveRecommendation:
    """
    Recommended additive concentrations and predicted effects.

    Attributes:
        dmso_percent: Recommended DMSO concentration (%)
        betaine_m: Recommended betaine concentration (M)
        trehalose_m: Recommended trehalose concentration (M)
        formamide_percent: Recommended formamide concentration (%)
        glycerol_percent: Recommended glycerol concentration (%)
        mg_conc: Recommended Mg2+ concentration (mM)
        ssb: Whether to use single-strand binding protein

        predicted_amplification: Expected amplification factor (0-1)
        predicted_processivity: Expected processivity factor
        predicted_accessibility: Expected template accessibility
        predicted_binding: Expected primer binding efficiency

        confidence: Confidence level ('high', 'medium', 'low')
        optimization_score: Combined optimization score
        warnings: List of warnings about the recommendation
        rationale: Human-readable explanation
    """
    # Recommended concentrations
    dmso_percent: float = 0.0
    betaine_m: float = 0.0
    trehalose_m: float = 0.0
    formamide_percent: float = 0.0
    glycerol_percent: float = 0.0
    mg_conc: float = 2.5
    ssb: bool = False

    # Predicted effects
    predicted_amplification: float = 0.0
    predicted_processivity: float = 1.0
    predicted_accessibility: float = 1.0
    predicted_binding: float = 0.0

    # Metadata
    confidence: str = 'medium'
    optimization_score: float = 0.0
    warnings: List[str] = field(default_factory=list)
    rationale: str = ""

    # Input parameters (for reference)
    primer_length: int = 0
    template_gc: float = 0.5
    polymerase: str = 'phi29'
    optimize_for: str = 'amplification'

    def to_conditions(self, temp: Optional[float] = None) -> ReactionConditions:
        """
        Convert recommendation to ReactionConditions.

        Args:
            temp: Optional temperature override (uses polymerase default if None)

        Returns:
            ReactionConditions with recommended additives
        """
        poly_params = get_polymerase_params(self.polymerase)
        reaction_temp = temp if temp is not None else poly_params['optimal_temp']

        return ReactionConditions(
            temp=reaction_temp,
            polymerase=self.polymerase,
            dmso_percent=self.dmso_percent,
            betaine_m=self.betaine_m,
            trehalose_m=self.trehalose_m,
            formamide_percent=self.formamide_percent,
            glycerol_percent=self.glycerol_percent,
            mg_conc=self.mg_conc,
            ssb=self.ssb,
        )

    def summary(self) -> str:
        """Generate human-readable summary of recommendation."""
        lines = [
            "=" * 60,
            "ADDITIVE OPTIMIZATION RECOMMENDATION",
            "=" * 60,
            "",
            f"Target: {self.primer_length}bp primers, {self.template_gc:.0%} GC template",
            f"Polymerase: {self.polymerase}",
            f"Optimization goal: {self.optimize_for}",
            "",
            "RECOMMENDED CONCENTRATIONS:",
            "-" * 30,
        ]

        # Only show non-zero additives
        if self.dmso_percent > 0:
            lines.append(f"  DMSO: {self.dmso_percent:.1f}%")
        if self.betaine_m > 0:
            lines.append(f"  Betaine: {self.betaine_m:.1f} M")
        if self.trehalose_m > 0:
            lines.append(f"  Trehalose: {self.trehalose_m:.2f} M")
        if self.formamide_percent > 0:
            lines.append(f"  Formamide: {self.formamide_percent:.1f}%")
        if self.glycerol_percent > 0:
            lines.append(f"  Glycerol: {self.glycerol_percent:.1f}%")
        lines.append(f"  Mg2+: {self.mg_conc:.1f} mM")
        if self.ssb:
            lines.append("  SSB: Yes (recommended)")

        lines.extend([
            "",
            "PREDICTED EFFECTS:",
            "-" * 30,
            f"  Amplification factor: {self.predicted_amplification:.2f}",
            f"  Processivity: {self.predicted_processivity:.2f}",
            f"  Template accessibility: {self.predicted_accessibility:.2f}",
            f"  Primer binding: {self.predicted_binding:.2f}",
            "",
            f"Optimization score: {self.optimization_score:.3f}",
            f"Confidence: {self.confidence}",
        ])

        if self.warnings:
            lines.extend([
                "",
                "WARNINGS:",
                "-" * 30,
            ])
            for warning in self.warnings:
                lines.append(f"  - {warning}")

        if self.rationale:
            lines.extend([
                "",
                "RATIONALE:",
                "-" * 30,
                self.rationale,
            ])

        lines.append("=" * 60)
        return "\n".join(lines)


class AdditiveOptimizer:
    """
    Optimize additive concentrations for SWGA.

    Uses a combination of grid search and rule-based constraints
    to find optimal additive combinations for specific conditions.

    Optimization goals:
    - 'amplification': Maximize overall amplification (default)
    - 'specificity': Maximize fg/bg discrimination
    - 'coverage': Maximize genome coverage
    - 'processivity': Maximize polymerase processivity

    Example:
        optimizer = AdditiveOptimizer(polymerase='equiphi29')
        rec = optimizer.optimize(primer_length=15, template_gc=0.65)
        conditions = rec.to_conditions()
    """

    # Additive search ranges (based on literature)
    SEARCH_RANGES = {
        'dmso_percent': [0.0, 2.0, 4.0, 5.0, 6.0, 8.0],
        'betaine_m': [0.0, 0.5, 1.0, 1.5, 2.0],
        'trehalose_m': [0.0, 0.1, 0.2, 0.3, 0.5],
        'mg_conc': [1.5, 2.0, 2.5, 3.0, 4.0],
    }

    # Quick search for faster optimization
    QUICK_RANGES = {
        'dmso_percent': [0.0, 3.0, 5.0],
        'betaine_m': [0.0, 1.0, 1.5],
        'trehalose_m': [0.0, 0.2],
        'mg_conc': [2.5],
    }

    # Constraints based on polymerase
    POLYMERASE_CONSTRAINTS = {
        'phi29': {
            'max_dmso': 8.0,
            'max_betaine': 2.5,
            'optimal_temp': 30.0,
            'max_primer_length_base': 12,
        },
        'equiphi29': {
            'max_dmso': 6.0,  # More sensitive at higher temp
            'max_betaine': 2.0,
            'optimal_temp': 42.0,
            'max_primer_length_base': 15,
        },
        'bst': {
            'max_dmso': 5.0,
            'max_betaine': 1.5,
            'optimal_temp': 60.0,
            'max_primer_length_base': 18,
        },
        'klenow': {
            'max_dmso': 10.0,
            'max_betaine': 2.0,
            'optimal_temp': 37.0,
            'max_primer_length_base': 12,
        },
    }

    def __init__(self, polymerase: str = 'phi29'):
        """
        Initialize optimizer for a specific polymerase.

        Args:
            polymerase: Polymerase type ('phi29', 'equiphi29', 'bst', 'klenow')
        """
        self.polymerase = polymerase
        self.poly_params = get_polymerase_params(polymerase)
        self.constraints = self.POLYMERASE_CONSTRAINTS.get(
            polymerase, self.POLYMERASE_CONSTRAINTS['phi29']
        )

    def optimize(
        self,
        primer_length: int,
        template_gc: float,
        optimize_for: str = 'amplification',
        primer_gc: Optional[float] = None,
        quick: bool = False,
        include_ssb: bool = True,
        include_formamide: bool = False,
        include_glycerol: bool = False,
    ) -> AdditiveRecommendation:
        """
        Find optimal additive concentrations.

        Args:
            primer_length: Length of primers in bp
            template_gc: Template genome GC content (0-1)
            optimize_for: Optimization goal:
                - 'amplification': Overall amplification efficiency
                - 'specificity': Maximize discrimination
                - 'coverage': Maximize genome coverage
                - 'processivity': Maximize extension length
            primer_gc: Primer GC content (estimated from template if not provided)
            quick: Use quick search (fewer combinations, faster)
            include_ssb: Consider SSB in optimization
            include_formamide: Include formamide in search
            include_glycerol: Include glycerol in search

        Returns:
            AdditiveRecommendation with optimal concentrations
        """
        # Estimate primer GC if not provided
        if primer_gc is None:
            # Primers typically have GC close to template
            primer_gc = template_gc

        # Get search ranges
        ranges = self.QUICK_RANGES if quick else self.SEARCH_RANGES

        # Filter ranges by polymerase constraints
        filtered_ranges = self._filter_ranges_by_constraints(ranges)

        # Generate all combinations
        combinations = self._generate_combinations(
            filtered_ranges,
            include_formamide=include_formamide,
            include_glycerol=include_glycerol,
        )

        # Evaluate each combination
        best_score = -float('inf')
        best_combo = None
        best_effects = None

        for combo in combinations:
            score, effects, warnings = self._evaluate_combination(
                combo,
                primer_length=primer_length,
                template_gc=template_gc,
                primer_gc=primer_gc,
                optimize_for=optimize_for,
            )

            if score > best_score:
                best_score = score
                best_combo = combo
                best_effects = effects
                best_warnings = warnings

        # Determine if SSB would help
        ssb_recommended = False
        if include_ssb and template_gc > 0.55:
            ssb_recommended = True

        # Build recommendation
        recommendation = self._build_recommendation(
            combo=best_combo,
            effects=best_effects,
            score=best_score,
            warnings=best_warnings,
            primer_length=primer_length,
            template_gc=template_gc,
            optimize_for=optimize_for,
            ssb=ssb_recommended,
        )

        return recommendation

    def _filter_ranges_by_constraints(
        self, ranges: Dict[str, List[float]]
    ) -> Dict[str, List[float]]:
        """Filter additive ranges by polymerase constraints."""
        filtered = {}

        for additive, values in ranges.items():
            if additive == 'dmso_percent':
                max_val = self.constraints['max_dmso']
                filtered[additive] = [v for v in values if v <= max_val]
            elif additive == 'betaine_m':
                max_val = self.constraints['max_betaine']
                filtered[additive] = [v for v in values if v <= max_val]
            else:
                filtered[additive] = values

        return filtered

    def _generate_combinations(
        self,
        ranges: Dict[str, List[float]],
        include_formamide: bool = False,
        include_glycerol: bool = False,
    ) -> List[Dict[str, float]]:
        """Generate all additive combinations."""
        # Base additives
        keys = ['dmso_percent', 'betaine_m', 'trehalose_m', 'mg_conc']
        value_lists = [ranges.get(k, [0.0]) for k in keys]

        # Optional additives
        if include_formamide:
            keys.append('formamide_percent')
            value_lists.append([0.0, 1.0, 2.0])
        if include_glycerol:
            keys.append('glycerol_percent')
            value_lists.append([0.0, 5.0, 10.0])

        # Generate combinations
        combinations = []
        for values in itertools.product(*value_lists):
            combo = dict(zip(keys, values))
            combinations.append(combo)

        return combinations

    def _evaluate_combination(
        self,
        combo: Dict[str, float],
        primer_length: int,
        template_gc: float,
        primer_gc: float,
        optimize_for: str,
    ) -> Tuple[float, MechanisticEffects, List[str]]:
        """
        Evaluate an additive combination using mechanistic model.

        Returns:
            Tuple of (score, effects, warnings)
        """
        warnings = []

        # Create conditions
        conditions = ReactionConditions(
            temp=self.constraints['optimal_temp'],
            polymerase=self.polymerase,
            dmso_percent=combo.get('dmso_percent', 0.0),
            betaine_m=combo.get('betaine_m', 0.0),
            trehalose_m=combo.get('trehalose_m', 0.0),
            formamide_percent=combo.get('formamide_percent', 0.0),
            glycerol_percent=combo.get('glycerol_percent', 0.0),
            mg_conc=combo.get('mg_conc', 2.5),
        )

        # Check primer length safety
        max_safe = conditions.max_primer_length(primer_gc)
        if primer_length > max_safe:
            warnings.append(
                f"Primer length {primer_length}bp exceeds safe maximum "
                f"({max_safe}bp) for these conditions"
            )

        # Create mechanistic model and evaluate
        model = MechanisticModel(conditions)

        # Create representative primer sequence
        # (approximation based on length and GC)
        gc_bases = int(primer_length * primer_gc)
        at_bases = primer_length - gc_bases
        primer = 'G' * (gc_bases // 2) + 'C' * (gc_bases - gc_bases // 2)
        primer += 'A' * (at_bases // 2) + 'T' * (at_bases - at_bases // 2)

        effects = model.calculate_effects(primer, template_gc)

        # Calculate score based on optimization goal
        if optimize_for == 'amplification':
            # Balance all factors
            score = effects.predicted_amplification_factor
        elif optimize_for == 'specificity':
            # Prioritize binding stability (lower koff)
            score = effects.effective_binding_rate * (1.0 / effects.koff_factor)
        elif optimize_for == 'coverage':
            # Prioritize processivity and accessibility
            score = effects.processivity_factor * effects.accessibility_factor
        elif optimize_for == 'processivity':
            # Just processivity
            score = effects.processivity_factor
        else:
            score = effects.predicted_amplification_factor

        # Penalize extreme conditions
        if combo.get('dmso_percent', 0) > 6.0:
            score *= 0.9  # High DMSO penalty
            warnings.append("High DMSO (>6%) may inhibit enzyme")
        if combo.get('betaine_m', 0) > 1.5:
            score *= 0.95
            warnings.append("High betaine (>1.5M) may affect viscosity")

        return score, effects, warnings

    def _build_recommendation(
        self,
        combo: Dict[str, float],
        effects: MechanisticEffects,
        score: float,
        warnings: List[str],
        primer_length: int,
        template_gc: float,
        optimize_for: str,
        ssb: bool,
    ) -> AdditiveRecommendation:
        """Build final recommendation from optimization results."""
        # Determine confidence
        if score > 0.5:
            confidence = 'high'
        elif score > 0.2:
            confidence = 'medium'
        else:
            confidence = 'low'

        # Build rationale
        rationale_parts = []

        dmso = combo.get('dmso_percent', 0)
        betaine = combo.get('betaine_m', 0)
        trehalose = combo.get('trehalose_m', 0)

        if dmso > 0:
            rationale_parts.append(
                f"DMSO at {dmso:.0f}% helps melt secondary structure "
                f"(accessibility +{(effects.accessibility_factor-0.5)*100:.0f}%)"
            )
        if betaine > 0:
            rationale_parts.append(
                f"Betaine at {betaine:.1f}M normalizes GC effects and "
                f"stabilizes enzyme"
            )
        if trehalose > 0:
            rationale_parts.append(
                f"Trehalose at {trehalose:.2f}M provides additional "
                f"enzyme stabilization"
            )
        if template_gc > 0.55:
            rationale_parts.append(
                f"High-GC template ({template_gc:.0%}) benefits from "
                f"GC-normalizing additives"
            )
        if ssb:
            rationale_parts.append(
                "SSB recommended for high-GC templates to maintain "
                "single-stranded regions"
            )

        rationale = " ".join(rationale_parts) if rationale_parts else (
            "Standard conditions without additives are suitable for this target."
        )

        return AdditiveRecommendation(
            dmso_percent=combo.get('dmso_percent', 0.0),
            betaine_m=combo.get('betaine_m', 0.0),
            trehalose_m=combo.get('trehalose_m', 0.0),
            formamide_percent=combo.get('formamide_percent', 0.0),
            glycerol_percent=combo.get('glycerol_percent', 0.0),
            mg_conc=combo.get('mg_conc', 2.5),
            ssb=ssb,
            predicted_amplification=effects.predicted_amplification_factor,
            predicted_processivity=effects.processivity_factor,
            predicted_accessibility=effects.accessibility_factor,
            predicted_binding=effects.effective_binding_rate,
            confidence=confidence,
            optimization_score=score,
            warnings=warnings,
            rationale=rationale,
            primer_length=primer_length,
            template_gc=template_gc,
            polymerase=self.polymerase,
            optimize_for=optimize_for,
        )

    def recommend_for_application(
        self,
        application: str,
        template_gc: float,
        primer_length: Optional[int] = None,
    ) -> AdditiveRecommendation:
        """
        Get additive recommendations for standard application profiles.

        Args:
            application: Application type:
                - 'discovery': Pathogen discovery (maximize coverage)
                - 'clinical': Diagnostics (maximize specificity)
                - 'enrichment': Sequencing enrichment (balanced)
                - 'metagenomics': Capture diversity (high coverage)
            template_gc: Template GC content (0-1)
            primer_length: Primer length (uses polymerase default if None)

        Returns:
            AdditiveRecommendation optimized for the application
        """
        # Default primer lengths by application and polymerase
        if primer_length is None:
            base_length = self.constraints['max_primer_length_base']
            if application == 'discovery':
                primer_length = base_length - 2  # Shorter for coverage
            elif application == 'clinical':
                primer_length = base_length + 2  # Longer for specificity
            elif application == 'metagenomics':
                primer_length = base_length - 3  # Shortest for diversity
            else:
                primer_length = base_length

        # Map application to optimization goal
        goal_map = {
            'discovery': 'coverage',
            'clinical': 'specificity',
            'enrichment': 'amplification',
            'metagenomics': 'coverage',
        }
        optimize_for = goal_map.get(application, 'amplification')

        return self.optimize(
            primer_length=primer_length,
            template_gc=template_gc,
            optimize_for=optimize_for,
        )


def optimize_conditions(
    primer_length: int,
    template_gc: float,
    polymerase: str = 'phi29',
    optimize_for: str = 'amplification',
) -> AdditiveRecommendation:
    """
    Convenience function for quick additive optimization.

    Args:
        primer_length: Length of primers in bp
        template_gc: Template genome GC content (0-1)
        polymerase: Polymerase type
        optimize_for: Optimization goal

    Returns:
        AdditiveRecommendation

    Example:
        >>> rec = optimize_conditions(15, 0.65, polymerase='equiphi29')
        >>> print(rec.summary())
    """
    optimizer = AdditiveOptimizer(polymerase=polymerase)
    return optimizer.optimize(
        primer_length=primer_length,
        template_gc=template_gc,
        optimize_for=optimize_for,
    )
