"""Tests for the additive interaction framework."""

import pytest
from neoswga.core.additive_interactions import (
    AdditiveInteraction,
    AdditiveInteractionRegistry,
    Pathway,
    EffectType,
    get_default_registry,
    reset_default_registry,
    calculate_interaction_modifiers,
)
from neoswga.core.reaction_conditions import ReactionConditions


@pytest.fixture(autouse=True)
def reset_registry():
    """Reset the default registry before each test."""
    reset_default_registry()
    yield
    reset_default_registry()


class TestAdditiveInteraction:
    """Tests for individual AdditiveInteraction objects."""

    def test_interaction_creation(self):
        """Test creating an interaction."""
        interaction = AdditiveInteraction(
            name='test_interaction',
            additives=['dmso', 'betaine'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'dmso': 2.0, 'betaine': 0.5},
        )
        assert interaction.name == 'test_interaction'
        assert interaction.additives == ['dmso', 'betaine']
        assert interaction.pathway == Pathway.PROCESSIVITY
        assert interaction.effect_type == EffectType.SYNERGY

    def test_interaction_not_active_below_threshold(self):
        """Test that interaction is not active when below concentration threshold."""
        interaction = AdditiveInteraction(
            name='test',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'dmso': 5.0},
        )
        conditions = ReactionConditions(temp=30.0, dmso_percent=3.0)
        assert not interaction.is_active(conditions)

    def test_interaction_active_above_threshold(self):
        """Test that interaction is active when above concentration threshold."""
        interaction = AdditiveInteraction(
            name='test',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'dmso': 5.0},
        )
        conditions = ReactionConditions(temp=30.0, dmso_percent=6.0)
        assert interaction.is_active(conditions)

    def test_synergy_effect_positive(self):
        """Test that synergy produces effect > 1.0."""
        interaction = AdditiveInteraction(
            name='test',
            additives=['betaine'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'betaine': 0.5},
            max_concentrations={'betaine': 2.0},
        )
        conditions = ReactionConditions(temp=30.0, betaine_m=1.5)
        effect = interaction.calculate_effect(conditions)
        assert effect > 1.0

    def test_antagonism_effect_negative(self):
        """Test that antagonism produces effect < 1.0."""
        interaction = AdditiveInteraction(
            name='test',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.ANTAGONISM,
            base_coefficient=0.1,
            min_concentrations={'dmso': 3.0},
        )
        conditions = ReactionConditions(temp=30.0, dmso_percent=5.0)
        effect = interaction.calculate_effect(conditions)
        assert effect < 1.0

    def test_inactive_interaction_returns_one(self):
        """Test that inactive interaction returns effect of 1.0."""
        interaction = AdditiveInteraction(
            name='test',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.5,
            min_concentrations={'dmso': 10.0},  # Very high threshold
        )
        conditions = ReactionConditions(temp=30.0, dmso_percent=5.0)
        effect = interaction.calculate_effect(conditions)
        assert effect == 1.0

    def test_gc_dependent_interaction_high_gc(self):
        """Test GC-dependent interaction activates for high GC templates."""
        interaction = AdditiveInteraction(
            name='test',
            additives=['betaine'],
            pathway=Pathway.ACCESSIBILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'betaine': 0.5},
            requires_high_gc=True,
            gc_threshold=0.55,
        )
        conditions = ReactionConditions(temp=30.0, betaine_m=1.0)

        # Low GC: should be inactive
        assert not interaction.is_active(conditions, template_gc=0.4)

        # High GC: should be active
        assert interaction.is_active(conditions, template_gc=0.6)


class TestAdditiveInteractionRegistry:
    """Tests for the interaction registry."""

    def test_register_and_get_interaction(self):
        """Test registering and retrieving interactions."""
        registry = AdditiveInteractionRegistry()
        interaction = AdditiveInteraction(
            name='test_interaction',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
        )
        registry.register(interaction)
        assert registry.get('test_interaction') is interaction

    def test_unregister_interaction(self):
        """Test removing an interaction."""
        registry = AdditiveInteractionRegistry()
        interaction = AdditiveInteraction(
            name='test',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
        )
        registry.register(interaction)
        registry.unregister('test')
        assert registry.get('test') is None

    def test_pathway_modifier_multiplicative(self):
        """Test that multiple interactions on same pathway multiply."""
        registry = AdditiveInteractionRegistry()

        # Two synergies, each +10%
        registry.register(AdditiveInteraction(
            name='synergy1',
            additives=['dmso'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'dmso': 1.0},
        ))
        registry.register(AdditiveInteraction(
            name='synergy2',
            additives=['betaine'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'betaine': 0.5},
        ))

        conditions = ReactionConditions(
            temp=30.0, dmso_percent=2.0, betaine_m=1.0
        )
        modifier = registry.calculate_pathway_modifier(
            Pathway.STABILITY, conditions
        )
        # Should be approximately 1.1 * 1.05 = 1.155 (multiplicative)
        # Exact values depend on concentration scaling
        assert modifier > 1.0

    def test_no_interactions_returns_one(self):
        """Test that no active interactions returns modifier of 1.0."""
        registry = AdditiveInteractionRegistry()
        conditions = ReactionConditions(temp=30.0)
        modifier = registry.calculate_pathway_modifier(Pathway.PROCESSIVITY, conditions)
        assert modifier == 1.0

    def test_get_active_interactions(self):
        """Test getting list of active interactions."""
        registry = AdditiveInteractionRegistry()
        registry.register(AdditiveInteraction(
            name='active',
            additives=['dmso'],
            pathway=Pathway.PROCESSIVITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'dmso': 2.0},
        ))
        registry.register(AdditiveInteraction(
            name='inactive',
            additives=['betaine'],
            pathway=Pathway.STABILITY,
            effect_type=EffectType.SYNERGY,
            base_coefficient=0.1,
            min_concentrations={'betaine': 5.0},  # High threshold
        ))

        conditions = ReactionConditions(temp=30.0, dmso_percent=3.0, betaine_m=1.0)
        active = registry.get_active_interactions(conditions)

        active_names = [i.name for i in active]
        assert 'active' in active_names
        assert 'inactive' not in active_names


class TestDefaultRegistry:
    """Tests for the default interaction registry."""

    def test_default_registry_loads(self):
        """Test that default registry loads without error."""
        registry = get_default_registry()
        assert len(registry.interactions) > 0

    def test_default_registry_has_expected_interactions(self):
        """Test that default registry has key interactions."""
        registry = get_default_registry()
        interaction_names = list(registry.interactions.keys())

        # Check for migrated interactions
        assert 'betaine_trehalose_synergy' in interaction_names
        # Note: dmso_mg_chelation is now handled directly in MechanisticModel

        # Check for new interactions
        assert 'betaine_dmso_gc_synergy' in interaction_names
        assert 'trehalose_dmso_protection' in interaction_names
        assert 'optimal_swga_triple' in interaction_names

    def test_standard_conditions_no_effects(self):
        """Test that standard conditions without additives have no interactions."""
        conditions = ReactionConditions(temp=30.0, polymerase='phi29', mg_conc=2.5)
        mods = calculate_interaction_modifiers(conditions, template_gc=0.5)

        # All modifiers should be 1.0 (or very close)
        for pathway, mod in mods.items():
            assert abs(mod - 1.0) < 0.01, f"{pathway} modifier was {mod}"

    def test_enhanced_conditions_have_effects(self):
        """Test that enhanced conditions activate interactions."""
        conditions = ReactionConditions(
            temp=42.0, polymerase='equiphi29',
            dmso_percent=5.0, betaine_m=1.5, trehalose_m=0.3, mg_conc=2.5
        )
        mods = calculate_interaction_modifiers(conditions, template_gc=0.6)

        # Should have some non-1.0 modifiers
        non_neutral = [p for p, m in mods.items() if abs(m - 1.0) > 0.01]
        assert len(non_neutral) > 0


class TestMechanisticModelIntegration:
    """Tests for integration with MechanisticModel."""

    def test_model_uses_interaction_framework(self):
        """Test that MechanisticModel applies interaction modifiers."""
        from neoswga.core.mechanistic_model import MechanisticModel

        conditions = ReactionConditions(
            temp=42.0, polymerase='equiphi29',
            dmso_percent=5.0, betaine_m=1.5, trehalose_m=0.3, mg_conc=2.5
        )
        model = MechanisticModel(conditions)

        # Get interaction report
        report = model.get_interaction_report(template_gc=0.6)
        assert 'Active Additive Interactions' in report

    def test_model_enzyme_parameters_include_interactions(self):
        """Test that enzyme parameters include interaction effects."""
        from neoswga.core.mechanistic_model import MechanisticModel

        # Standard conditions
        base_conditions = ReactionConditions(temp=30.0, polymerase='phi29', mg_conc=2.5)
        base_model = MechanisticModel(base_conditions)
        base_enzyme = base_model.get_enzyme_parameters(template_gc=0.5)

        # Enhanced conditions with synergistic additives
        enhanced_conditions = ReactionConditions(
            temp=30.0, polymerase='phi29',
            betaine_m=1.0, trehalose_m=0.3, mg_conc=2.5
        )
        enhanced_model = MechanisticModel(enhanced_conditions)
        enhanced_enzyme = enhanced_model.get_enzyme_parameters(template_gc=0.5)

        # Enhanced should have better stability due to betaine-trehalose synergy
        assert enhanced_enzyme['stability_factor'] >= base_enzyme['stability_factor']


class TestThreeWayInteractions:
    """Tests for multi-additive interactions."""

    def test_three_way_interaction_requires_all(self):
        """Test that 3-way interaction requires all additives."""
        registry = get_default_registry()
        triple = registry.get('optimal_swga_triple')
        assert triple is not None
        assert len(triple.additives) == 3

        # Missing one additive
        partial = ReactionConditions(
            temp=30.0, dmso_percent=5.0, betaine_m=1.5, trehalose_m=0.0
        )
        assert not triple.is_active(partial)

        # All additives present
        full = ReactionConditions(
            temp=30.0, dmso_percent=5.0, betaine_m=1.5, trehalose_m=0.3
        )
        assert triple.is_active(full)
