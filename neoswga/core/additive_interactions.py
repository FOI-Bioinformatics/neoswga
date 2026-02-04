"""
Systematic framework for modeling additive interactions in SWGA.

Replaces ad-hoc hardcoded interactions with a configurable rule-based system.
Each interaction is defined by:
- Additives involved (2 or more)
- Pathway affected (processivity, stability, speed, tm, accessibility, kon, koff)
- Effect type (synergy or antagonism)
- Concentration thresholds
- Temperature dependence (optional)
- Literature reference

This framework supports:
1. Pairwise interactions (2 additives)
2. Multi-way interactions (3+ additives)
3. Temperature-dependent effect strengths
4. GC-dependent activation
5. Easy addition of new interactions from literature

Literature basis:
- Henke et al. (1997): Betaine effects on PCR
- Sarkar et al. (1990): DMSO in PCR
- Musso et al. (2006): SWGA additive combinations
- Varadharajan et al. (2017): EquiPhi29 optimization
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, TYPE_CHECKING, Callable
from enum import Enum
import math

if TYPE_CHECKING:
    from neoswga.core.reaction_conditions import ReactionConditions


class Pathway(Enum):
    """Pathways that can be affected by additive interactions."""
    PROCESSIVITY = "processivity"
    STABILITY = "stability"
    SPEED = "speed"
    TM = "tm"
    ACCESSIBILITY = "accessibility"
    KON = "kon"
    KOFF = "koff"


class EffectType(Enum):
    """Type of interaction effect."""
    SYNERGY = "synergy"           # Combined effect > sum of parts
    ANTAGONISM = "antagonism"     # Additives counteract each other
    CONDITIONAL = "conditional"   # Effect depends on specific conditions


@dataclass
class AdditiveInteraction:
    """
    Defines an interaction between two or more additives.

    Each interaction specifies:
    - Which additives are involved
    - What pathway is affected
    - Whether it's synergistic or antagonistic
    - Concentration thresholds for activation
    - The magnitude of the effect

    Attributes:
        name: Unique identifier for this interaction
        additives: List of additive names involved (e.g., ['dmso', 'betaine'])
        pathway: Which pathway is affected
        effect_type: Synergy, antagonism, or conditional
        base_coefficient: Base effect magnitude (multiplied by concentrations)
        min_concentrations: Minimum concentrations for activation
        max_concentrations: Maximum concentrations to consider (caps effect)
        temperature_coefficient: How effect changes with temperature deviation
        reference_temp: Reference temperature for temperature coefficient
        requires_high_gc: Only active for high-GC templates (>55%)
        requires_low_gc: Only active for low-GC templates (<45%)
        description: Human-readable description of the mechanism
        reference: Literature reference
    """
    name: str
    additives: List[str]
    pathway: Pathway
    effect_type: EffectType
    base_coefficient: float

    # Concentration thresholds
    min_concentrations: Dict[str, float] = field(default_factory=dict)
    max_concentrations: Dict[str, float] = field(default_factory=dict)

    # Temperature dependence
    temperature_coefficient: float = 0.0  # Effect per degree from reference
    reference_temp: float = 30.0  # Reference temperature (Celsius)

    # Template GC dependence
    requires_high_gc: bool = False
    requires_low_gc: bool = False
    gc_threshold: float = 0.55  # Threshold for high/low GC

    # Documentation
    description: str = ""
    reference: str = ""

    def get_additive_concentration(self, conditions: 'ReactionConditions',
                                    additive: str) -> float:
        """Get concentration of an additive from conditions."""
        # Map additive names to condition attributes
        mapping = {
            'dmso': conditions.dmso_percent,
            'betaine': conditions.betaine_m,
            'trehalose': conditions.trehalose_m,
            'formamide': conditions.formamide_percent,
            'mg': conditions.mg_conc,
            'na': conditions.na_conc,
            'glycerol': conditions.glycerol_percent,
            'peg': conditions.peg_percent,
            'bsa': conditions.bsa_ug_ml / 100.0,  # Normalize to useful scale
            'ssb': 1.0 if conditions.ssb else 0.0,
            'ethanol': conditions.ethanol_percent,
            'urea': conditions.urea_m,
            'tmac': conditions.tmac_m,
        }
        return mapping.get(additive.lower(), 0.0)

    def is_active(self, conditions: 'ReactionConditions',
                  template_gc: float = 0.5) -> bool:
        """
        Check if this interaction is active under given conditions.

        Args:
            conditions: Reaction conditions
            template_gc: Template GC content (0-1)

        Returns:
            True if interaction is active
        """
        # Check GC requirements
        if self.requires_high_gc and template_gc < self.gc_threshold:
            return False
        if self.requires_low_gc and template_gc > (1 - self.gc_threshold + 0.1):
            return False

        # Check minimum concentrations for all additives
        for additive in self.additives:
            conc = self.get_additive_concentration(conditions, additive)
            min_conc = self.min_concentrations.get(additive, 0.0)
            if conc < min_conc:
                return False

        return True

    def calculate_effect(self, conditions: 'ReactionConditions',
                         template_gc: float = 0.5) -> float:
        """
        Calculate the effect multiplier for this interaction.

        For synergies, returns a value > 1.0 (boost)
        For antagonisms, returns a value < 1.0 (reduction)

        Args:
            conditions: Reaction conditions
            template_gc: Template GC content (0-1)

        Returns:
            Effect multiplier (1.0 = no effect)
        """
        if not self.is_active(conditions, template_gc):
            return 1.0

        # Calculate concentration-based effect magnitude
        effect_magnitude = self.base_coefficient

        for additive in self.additives:
            conc = self.get_additive_concentration(conditions, additive)
            min_conc = self.min_concentrations.get(additive, 0.0)
            max_conc = self.max_concentrations.get(additive, float('inf'))

            # Use concentration above threshold, capped at max
            effective_conc = min(conc, max_conc) - min_conc
            effective_conc = max(0.0, effective_conc)

            # Scale effect by concentration
            effect_magnitude *= effective_conc

        # Apply temperature dependence
        if self.temperature_coefficient != 0.0:
            temp_delta = conditions.temp - self.reference_temp
            temp_modifier = 1.0 + self.temperature_coefficient * temp_delta
            effect_magnitude *= max(0.0, temp_modifier)

        # Convert to multiplier based on effect type
        if self.effect_type == EffectType.SYNERGY:
            return 1.0 + effect_magnitude
        elif self.effect_type == EffectType.ANTAGONISM:
            return max(0.1, 1.0 - effect_magnitude)
        else:  # CONDITIONAL
            return 1.0 + effect_magnitude  # Treat as synergy by default


class AdditiveInteractionRegistry:
    """
    Registry of all known additive interactions.

    Provides methods to:
    - Register new interactions
    - Query active interactions
    - Calculate combined effects on pathways
    - Validate interaction consistency

    Usage:
        registry = AdditiveInteractionRegistry()
        registry.load_defaults()

        # Get combined effect on processivity
        proc_modifier = registry.calculate_pathway_modifier(
            Pathway.PROCESSIVITY, conditions, template_gc=0.5
        )
    """

    def __init__(self):
        self.interactions: Dict[str, AdditiveInteraction] = {}
        self._pathway_cache: Dict[str, List[AdditiveInteraction]] = {}

    def register(self, interaction: AdditiveInteraction) -> None:
        """Register a new interaction."""
        self.interactions[interaction.name] = interaction
        self._invalidate_cache()

    def unregister(self, name: str) -> None:
        """Remove an interaction by name."""
        if name in self.interactions:
            del self.interactions[name]
            self._invalidate_cache()

    def get(self, name: str) -> Optional[AdditiveInteraction]:
        """Get an interaction by name."""
        return self.interactions.get(name)

    def _invalidate_cache(self) -> None:
        """Clear the pathway cache."""
        self._pathway_cache.clear()

    def _get_pathway_interactions(self, pathway: Pathway) -> List[AdditiveInteraction]:
        """Get all interactions affecting a pathway (cached)."""
        key = pathway.value
        if key not in self._pathway_cache:
            self._pathway_cache[key] = [
                i for i in self.interactions.values()
                if i.pathway == pathway
            ]
        return self._pathway_cache[key]

    def get_active_interactions(self, conditions: 'ReactionConditions',
                                 template_gc: float = 0.5) -> List[AdditiveInteraction]:
        """Get all interactions that are active under given conditions."""
        return [
            i for i in self.interactions.values()
            if i.is_active(conditions, template_gc)
        ]

    def calculate_pathway_modifier(self, pathway: Pathway,
                                    conditions: 'ReactionConditions',
                                    template_gc: float = 0.5) -> float:
        """
        Calculate combined effect of all active interactions on a pathway.

        Effects are multiplicative: if two synergies each give 1.1x,
        the combined effect is 1.1 * 1.1 = 1.21x

        Args:
            pathway: The pathway to calculate
            conditions: Reaction conditions
            template_gc: Template GC content (0-1)

        Returns:
            Combined multiplier (1.0 = no effect)
        """
        modifier = 1.0
        for interaction in self._get_pathway_interactions(pathway):
            modifier *= interaction.calculate_effect(conditions, template_gc)
        return modifier

    def calculate_all_modifiers(self, conditions: 'ReactionConditions',
                                 template_gc: float = 0.5) -> Dict[str, float]:
        """
        Calculate modifiers for all pathways.

        Args:
            conditions: Reaction conditions
            template_gc: Template GC content (0-1)

        Returns:
            Dictionary mapping pathway names to modifiers
        """
        return {
            pathway.value: self.calculate_pathway_modifier(pathway, conditions, template_gc)
            for pathway in Pathway
        }

    def get_interaction_report(self, conditions: 'ReactionConditions',
                                template_gc: float = 0.5) -> str:
        """
        Generate a human-readable report of active interactions.

        Args:
            conditions: Reaction conditions
            template_gc: Template GC content

        Returns:
            Formatted report string
        """
        active = self.get_active_interactions(conditions, template_gc)

        if not active:
            return "No additive interactions active under current conditions."

        lines = ["Active Additive Interactions:", "=" * 40]

        for interaction in active:
            effect = interaction.calculate_effect(conditions, template_gc)
            effect_str = f"+{(effect-1)*100:.1f}%" if effect > 1 else f"{(effect-1)*100:.1f}%"

            lines.append(f"\n{interaction.name}:")
            lines.append(f"  Additives: {', '.join(interaction.additives)}")
            lines.append(f"  Pathway: {interaction.pathway.value}")
            lines.append(f"  Effect: {effect_str}")
            if interaction.description:
                lines.append(f"  Mechanism: {interaction.description}")

        # Summary
        lines.append("\n" + "=" * 40)
        lines.append("Combined effects by pathway:")
        modifiers = self.calculate_all_modifiers(conditions, template_gc)
        for pathway, mod in modifiers.items():
            if abs(mod - 1.0) > 0.001:
                effect_str = f"+{(mod-1)*100:.1f}%" if mod > 1 else f"{(mod-1)*100:.1f}%"
                lines.append(f"  {pathway}: {effect_str}")

        return "\n".join(lines)

    def load_defaults(self) -> None:
        """Load default interactions based on literature."""

        # =====================================================================
        # Existing interactions (migrated from mechanistic_model.py)
        # =====================================================================

        # Note: DMSO-Mg chelation is handled directly in MechanisticModel._calculate_base_enzyme_activity()
        # because it has an inverse relationship with Mg (higher Mg = less chelation effect)
        # that doesn't fit the standard interaction model.
        # Original: processivity *= (1 - 0.08 * (dmso - 3.0) * (3.0 - mg) / 3.0)

        # Betaine-trehalose synergy
        # Both stabilize proteins through different mechanisms
        # Combined effect is greater than sum of individual effects
        self.register(AdditiveInteraction(
            name='betaine_trehalose_synergy',
            additives=['betaine', 'trehalose'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.15,
            min_concentrations={'betaine': 0.5, 'trehalose': 0.1},
            max_concentrations={'betaine': 1.5, 'trehalose': 0.5},
            description="Combined osmolyte protection enhances enzyme stability",
            reference="Vasilescu et al. (2008)"
        ))

        # Betaine-trehalose also boosts processivity slightly
        self.register(AdditiveInteraction(
            name='betaine_trehalose_processivity',
            additives=['betaine', 'trehalose'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.05,  # 0.3 * 0.15 from original
            min_concentrations={'betaine': 0.5, 'trehalose': 0.1},
            max_concentrations={'betaine': 1.5, 'trehalose': 0.5},
            description="Enzyme stabilization indirectly improves processivity",
            reference="Vasilescu et al. (2008)"
        ))

        # DMSO-formamide antagonism
        # Both destabilize DNA but compete for similar binding sites on enzyme
        # Combined destabilization is less than expected
        self.register(AdditiveInteraction(
            name='dmso_formamide_antagonism',
            additives=['dmso', 'formamide'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,  # Actually reduces penalty, so synergy on stability
            base_coefficient=0.015,  # 0.03 * 0.5 from original
            min_concentrations={'dmso': 2.0, 'formamide': 1.0},
            max_concentrations={'dmso': 8.0, 'formamide': 5.0},
            description="Competing destabilizers partially neutralize each other",
            reference="Sarkar et al. (1990)"
        ))

        # =====================================================================
        # New interactions (previously missing)
        # =====================================================================

        # Betaine-DMSO interaction
        # Both affect GC normalization but through different mechanisms
        # Betaine reduces GC Tm penalty, DMSO destabilizes all duplexes
        # Combined: enhanced ability to handle high-GC templates
        self.register(AdditiveInteraction(
            name='betaine_dmso_gc_synergy',
            additives=['betaine', 'dmso'],
            pathway=Pathway.ACCESSIBILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.08,
            min_concentrations={'betaine': 0.5, 'dmso': 2.0},
            max_concentrations={'betaine': 2.0, 'dmso': 8.0},
            requires_high_gc=True,
            gc_threshold=0.55,
            description="Combined GC normalization improves high-GC template accessibility",
            reference="Henke et al. (1997)"
        ))

        # Betaine-DMSO also affects Tm
        self.register(AdditiveInteraction(
            name='betaine_dmso_tm_enhancement',
            additives=['betaine', 'dmso'],
            pathway=Pathway.TM,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.02,
            min_concentrations={'betaine': 0.5, 'dmso': 2.0},
            max_concentrations={'betaine': 2.0, 'dmso': 8.0},
            requires_high_gc=True,
            gc_threshold=0.55,
            description="Combined effect reduces Tm depression for high-GC primers",
            reference="Henke et al. (1997)"
        ))

        # Trehalose-DMSO interaction
        # Trehalose stabilizes enzyme while DMSO improves DNA accessibility
        # Combined: can use higher DMSO with less enzyme inhibition
        self.register(AdditiveInteraction(
            name='trehalose_dmso_protection',
            additives=['trehalose', 'dmso'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.12,
            min_concentrations={'trehalose': 0.1, 'dmso': 3.0},
            max_concentrations={'trehalose': 0.5, 'dmso': 8.0},
            description="Trehalose protects enzyme from DMSO destabilization",
            reference="Musso et al. (2006)"
        ))

        # SSB-betaine synergy for binding
        # SSB keeps DNA single-stranded, betaine helps primers access
        self.register(AdditiveInteraction(
            name='ssb_betaine_binding',
            additives=['ssb', 'betaine'],
            pathway=Pathway.KON,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.3,
            min_concentrations={'ssb': 1.0, 'betaine': 0.5},
            max_concentrations={'betaine': 2.0},
            description="SSB maintains ssDNA while betaine enhances primer binding",
            reference="Varadharajan et al. (2017)"
        ))

        # PEG-betaine interaction
        # Both are crowding agents - excessive combined crowding can inhibit
        self.register(AdditiveInteraction(
            name='peg_betaine_crowding',
            additives=['peg', 'betaine'],
            pathway=Pathway.SPEED,
            effect_type=EffectType.ANTAGONISM,
            base_coefficient=0.02,
            min_concentrations={'peg': 2.0, 'betaine': 1.0},
            max_concentrations={'peg': 10.0, 'betaine': 2.0},
            description="Excessive molecular crowding can slow polymerase",
            reference="Ralser et al. (2006)"
        ))

        # =====================================================================
        # Temperature-dependent interactions
        # =====================================================================

        # High-temperature betaine effectiveness
        # Betaine is more effective at higher temperatures (EquiPhi29 conditions)
        self.register(AdditiveInteraction(
            name='betaine_high_temp_boost',
            additives=['betaine'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.05,
            min_concentrations={'betaine': 0.5},
            max_concentrations={'betaine': 2.0},
            temperature_coefficient=0.02,  # Effect increases with temp
            reference_temp=30.0,
            description="Betaine stabilization is more important at elevated temperatures",
            reference="Varadharajan et al. (2017)"
        ))

        # =====================================================================
        # Three-way interactions
        # =====================================================================

        # Betaine-DMSO-trehalose triple
        # Optimal SWGA additive combination
        self.register(AdditiveInteraction(
            name='optimal_swga_triple',
            additives=['betaine', 'dmso', 'trehalose'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.03,
            min_concentrations={'betaine': 1.0, 'dmso': 3.0, 'trehalose': 0.2},
            max_concentrations={'betaine': 2.0, 'dmso': 6.0, 'trehalose': 0.5},
            description="Optimal three-additive combination for enhanced SWGA",
            reference="Musso et al. (2006)"
        ))

        # Betaine-DMSO-trehalose also enhances accessibility
        self.register(AdditiveInteraction(
            name='optimal_swga_triple_accessibility',
            additives=['betaine', 'dmso', 'trehalose'],
            pathway=Pathway.ACCESSIBILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.04,
            min_concentrations={'betaine': 1.0, 'dmso': 3.0, 'trehalose': 0.2},
            max_concentrations={'betaine': 2.0, 'dmso': 6.0, 'trehalose': 0.5},
            requires_high_gc=True,
            gc_threshold=0.50,
            description="Enhanced template accessibility with optimal additive mix",
            reference="Musso et al. (2006)"
        ))


# Global default registry
_default_registry: Optional[AdditiveInteractionRegistry] = None


def get_default_registry() -> AdditiveInteractionRegistry:
    """Get the default interaction registry (lazy-loaded singleton)."""
    global _default_registry
    if _default_registry is None:
        _default_registry = AdditiveInteractionRegistry()
        _default_registry.load_defaults()
    return _default_registry


def reset_default_registry() -> None:
    """Reset the default registry (for testing)."""
    global _default_registry
    _default_registry = None


def calculate_interaction_modifiers(conditions: 'ReactionConditions',
                                     template_gc: float = 0.5,
                                     registry: Optional[AdditiveInteractionRegistry] = None
                                     ) -> Dict[str, float]:
    """
    Convenience function to calculate all interaction modifiers.

    Args:
        conditions: Reaction conditions
        template_gc: Template GC content (0-1)
        registry: Optional custom registry (uses default if None)

    Returns:
        Dictionary mapping pathway names to modifiers
    """
    if registry is None:
        registry = get_default_registry()
    return registry.calculate_all_modifiers(conditions, template_gc)
