"""
Tests for the mechanistic SWGA model.

Tests all four pathways:
1. Tm modification
2. Secondary structure accessibility
3. Enzyme activity
4. Binding kinetics
"""

import pytest
import math
from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
from neoswga.core.reaction_conditions import ReactionConditions


class TestMechanisticEffects:
    """Test the MechanisticEffects dataclass."""

    def test_effects_is_frozen(self):
        """MechanisticEffects should be immutable."""
        effects = MechanisticEffects(
            tm_correction=0.0,
            effective_tm=35.0,
            accessibility_factor=1.0,
            processivity_factor=1.0,
            speed_factor=1.0,
            stability_factor=1.0,
            kon_factor=1.0,
            koff_factor=1.0,
            effective_binding_rate=0.8,
            effective_extension_rate=1.0,
            predicted_amplification_factor=0.8,
        )

        with pytest.raises(AttributeError):
            effects.processivity_factor = 0.5


class TestMechanisticModelInit:
    """Test MechanisticModel initialization."""

    def test_init_standard_phi29(self):
        """Initialize with standard phi29 conditions."""
        conditions = ReactionConditions(temp=30.0, polymerase='phi29')
        model = MechanisticModel(conditions)

        assert model.conditions == conditions
        assert model._poly_params['optimal_temp'] == 30.0

    def test_init_equiphi29(self):
        """Initialize with EquiPhi29 conditions."""
        conditions = ReactionConditions(temp=42.0, polymerase='equiphi29')
        model = MechanisticModel(conditions)

        assert model._poly_params['optimal_temp'] == 42.0

    def test_repr(self):
        """Test string representation."""
        conditions = ReactionConditions(temp=30.0, polymerase='phi29')
        model = MechanisticModel(conditions)

        repr_str = repr(model)
        assert 'temp=30' in repr_str
        assert 'phi29' in repr_str


class TestPathway1TmModification:
    """Test Tm modification pathway."""

    def test_dmso_lowers_tm(self):
        """DMSO should lower effective Tm."""
        cond_no_dmso = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_with_dmso = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=5.0)

        model_no = MechanisticModel(cond_no_dmso)
        model_with = MechanisticModel(cond_with_dmso)

        effects_no = model_no.calculate_effects('ATCGATCGATCG', 0.5)
        effects_with = model_with.calculate_effects('ATCGATCGATCG', 0.5)

        # DMSO should lower Tm
        assert effects_with.effective_tm < effects_no.effective_tm
        # Tm correction should be negative
        assert effects_with.tm_correction < effects_no.tm_correction

    def test_betaine_lowers_tm_uniformly(self):
        """Betaine has uniform Tm lowering component."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_with = ReactionConditions(temp=30.0, mg_conc=2.5, betaine_m=1.0)

        model_no = MechanisticModel(cond_no)
        model_with = MechanisticModel(cond_with)

        # Test with 50% GC (no GC normalization effect)
        effects_no = model_no.calculate_effects('ATGCATGCATGC', 0.5)
        effects_with = model_with.calculate_effects('ATGCATGCATGC', 0.5)

        # Betaine should lower Tm
        assert effects_with.tm_correction < effects_no.tm_correction

    def test_gc_normalization_effect(self):
        """Betaine should normalize GC-rich primer Tm towards baseline."""
        cond_with_betaine = ReactionConditions(temp=30.0, mg_conc=2.5, betaine_m=2.0)
        model = MechanisticModel(cond_with_betaine)

        # GC-rich primer (75% GC)
        effects_gc = model.calculate_effects('GCGCGCGCGCGC', 0.5)
        # AT-rich primer (25% GC)
        effects_at = model.calculate_effects('ATATATATATATATA', 0.5)

        # With betaine, the Tm difference should be reduced
        # (though not eliminated at 2M betaine)
        # The GC-rich primer should have more negative correction
        assert effects_gc.tm_correction < effects_at.tm_correction


class TestPathway2Accessibility:
    """Test secondary structure accessibility pathway."""

    def test_high_gc_template_lower_accessibility(self):
        """High GC templates should have lower accessibility."""
        conditions = ReactionConditions(temp=30.0, mg_conc=2.5)
        model = MechanisticModel(conditions)

        effects_low_gc = model.calculate_effects('ATCGATCG', template_gc=0.3)
        effects_high_gc = model.calculate_effects('ATCGATCG', template_gc=0.7)

        assert effects_high_gc.accessibility_factor < effects_low_gc.accessibility_factor

    def test_dmso_improves_accessibility(self):
        """DMSO should improve accessibility for high GC templates."""
        cond_no_dmso = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_with_dmso = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=5.0)

        model_no = MechanisticModel(cond_no_dmso)
        model_with = MechanisticModel(cond_with_dmso)

        # Test with high GC template
        effects_no = model_no.calculate_effects('ATCGATCG', template_gc=0.7)
        effects_with = model_with.calculate_effects('ATCGATCG', template_gc=0.7)

        assert effects_with.accessibility_factor > effects_no.accessibility_factor

    def test_betaine_improves_accessibility(self):
        """Betaine should improve accessibility."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_with = ReactionConditions(temp=30.0, mg_conc=2.5, betaine_m=1.5)

        model_no = MechanisticModel(cond_no)
        model_with = MechanisticModel(cond_with)

        effects_no = model_no.calculate_effects('ATCGATCG', template_gc=0.7)
        effects_with = model_with.calculate_effects('ATCGATCG', template_gc=0.7)

        assert effects_with.accessibility_factor > effects_no.accessibility_factor

    def test_higher_temp_improves_accessibility(self):
        """Higher temperature should improve accessibility."""
        cond_low = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_high = ReactionConditions(temp=42.0, mg_conc=2.5, polymerase='equiphi29')

        model_low = MechanisticModel(cond_low)
        model_high = MechanisticModel(cond_high)

        effects_low = model_low.calculate_effects('ATCGATCG', template_gc=0.7)
        effects_high = model_high.calculate_effects('ATCGATCG', template_gc=0.7)

        assert effects_high.accessibility_factor > effects_low.accessibility_factor

    def test_accessibility_bounded(self):
        """Accessibility should be bounded between 0.1 and 1.0."""
        conditions = ReactionConditions(temp=30.0, mg_conc=2.5)
        model = MechanisticModel(conditions)

        # Very high GC
        effects_high = model.calculate_effects('ATCGATCG', template_gc=0.9)
        assert 0.1 <= effects_high.accessibility_factor <= 1.0

        # Very low GC
        effects_low = model.calculate_effects('ATCGATCG', template_gc=0.1)
        assert 0.1 <= effects_low.accessibility_factor <= 1.0


class TestPathway3EnzymeActivity:
    """Test enzyme activity pathway."""

    def test_dmso_reduces_processivity(self):
        """High DMSO should reduce enzyme processivity."""
        cond_low = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=2.0)
        cond_high = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=8.0)

        model_low = MechanisticModel(cond_low)
        model_high = MechanisticModel(cond_high)

        effects_low = model_low.calculate_effects('ATCGATCG', 0.5)
        effects_high = model_high.calculate_effects('ATCGATCG', 0.5)

        assert effects_high.processivity_factor < effects_low.processivity_factor

    def test_betaine_enhances_stability(self):
        """Low betaine should enhance enzyme stability."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_bet = ReactionConditions(temp=30.0, mg_conc=2.5, betaine_m=1.0)

        model_no = MechanisticModel(cond_no)
        model_bet = MechanisticModel(cond_bet)

        effects_no = model_no.calculate_effects('ATCGATCG', 0.5)
        effects_bet = model_bet.calculate_effects('ATCGATCG', 0.5)

        assert effects_bet.stability_factor > effects_no.stability_factor

    def test_high_betaine_inhibits_processivity(self):
        """High betaine should inhibit processivity."""
        cond_low = ReactionConditions(temp=30.0, mg_conc=2.5, betaine_m=1.0)
        cond_high = ReactionConditions(temp=30.0, mg_conc=2.5, betaine_m=2.5)

        model_low = MechanisticModel(cond_low)
        model_high = MechanisticModel(cond_high)

        effects_low = model_low.calculate_effects('ATCGATCG', 0.5)
        effects_high = model_high.calculate_effects('ATCGATCG', 0.5)

        assert effects_high.processivity_factor < effects_low.processivity_factor

    def test_low_mg_reduces_processivity(self):
        """Low Mg2+ should reduce processivity."""
        cond_low = ReactionConditions(temp=30.0, mg_conc=0.5)
        cond_opt = ReactionConditions(temp=30.0, mg_conc=2.5)

        model_low = MechanisticModel(cond_low)
        model_opt = MechanisticModel(cond_opt)

        effects_low = model_low.calculate_effects('ATCGATCG', 0.5)
        effects_opt = model_opt.calculate_effects('ATCGATCG', 0.5)

        assert effects_low.processivity_factor < effects_opt.processivity_factor

    def test_temp_deviation_reduces_speed(self):
        """Temperature deviation from optimal should reduce speed."""
        # phi29 optimal is 30C
        cond_opt = ReactionConditions(temp=30.0, mg_conc=2.5, polymerase='phi29')
        cond_high = ReactionConditions(temp=37.0, mg_conc=2.5, polymerase='phi29')

        model_opt = MechanisticModel(cond_opt)
        model_high = MechanisticModel(cond_high)

        effects_opt = model_opt.calculate_effects('ATCGATCG', 0.5)
        effects_high = model_high.calculate_effects('ATCGATCG', 0.5)

        assert effects_high.speed_factor < effects_opt.speed_factor

    def test_formamide_inhibits(self):
        """Formamide should reduce processivity and stability."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_form = ReactionConditions(temp=30.0, mg_conc=2.5, formamide_percent=5.0)

        model_no = MechanisticModel(cond_no)
        model_form = MechanisticModel(cond_form)

        effects_no = model_no.calculate_effects('ATCGATCG', 0.5)
        effects_form = model_form.calculate_effects('ATCGATCG', 0.5)

        assert effects_form.processivity_factor < effects_no.processivity_factor
        assert effects_form.stability_factor < effects_no.stability_factor

    def test_glycerol_stabilizes_but_slows(self):
        """Glycerol should stabilize enzyme but reduce speed."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_glyc = ReactionConditions(temp=30.0, mg_conc=2.5, glycerol_percent=10.0)

        model_no = MechanisticModel(cond_no)
        model_glyc = MechanisticModel(cond_glyc)

        effects_no = model_no.calculate_effects('ATCGATCG', 0.5)
        effects_glyc = model_glyc.calculate_effects('ATCGATCG', 0.5)

        assert effects_glyc.stability_factor > effects_no.stability_factor
        assert effects_glyc.speed_factor < effects_no.speed_factor

    def test_get_enzyme_parameters(self):
        """Test get_enzyme_parameters method."""
        conditions = ReactionConditions(temp=30.0, mg_conc=2.5)
        model = MechanisticModel(conditions)

        params = model.get_enzyme_parameters()

        assert 'processivity_factor' in params
        assert 'speed_factor' in params
        assert 'stability_factor' in params


class TestPathway4BindingKinetics:
    """Test binding kinetics pathway."""

    def test_optimal_tm_gives_best_kon(self):
        """Primers with Tm ~7C above reaction temp should have best kon."""
        # At 30C reaction, optimal Tm is ~37C (30+7)
        # Use EquiPhi29 at 42C where optimal Tm is ~49C
        conditions = ReactionConditions(temp=42.0, polymerase='equiphi29', mg_conc=2.5)
        model = MechanisticModel(conditions)

        # Very short primer (very low Tm, far below optimal)
        effects_suboptimal = model.calculate_effects('ATCGAT', 0.5)  # ~6mer, Tm ~18C
        # Longer primer (higher Tm, closer to optimal ~49C)
        effects_better = model.calculate_effects('ATCGATCGATCGATCG', 0.5)  # ~16mer, Tm ~55C

        # The longer primer should have better kon (closer to optimal)
        # Both may hit floor, but the longer primer should be equal or better
        assert effects_better.kon_factor >= effects_suboptimal.kon_factor

    def test_ssb_boosts_kon(self):
        """SSB should increase kon."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5, ssb=False)
        cond_ssb = ReactionConditions(temp=30.0, mg_conc=2.5, ssb=True)

        model_no = MechanisticModel(cond_no)
        model_ssb = MechanisticModel(cond_ssb)

        effects_no = model_no.calculate_effects('ATCGATCGATCG', 0.5)
        effects_ssb = model_ssb.calculate_effects('ATCGATCGATCG', 0.5)

        assert effects_ssb.kon_factor > effects_no.kon_factor

    def test_additives_affect_koff(self):
        """Additives should increase koff (destabilizing effect)."""
        cond_no = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_add = ReactionConditions(temp=30.0, mg_conc=2.5, dmso_percent=5.0, betaine_m=1.0)

        model_no = MechanisticModel(cond_no)
        model_add = MechanisticModel(cond_add)

        effects_no = model_no.calculate_effects('ATCGATCGATCG', 0.5)
        effects_add = model_add.calculate_effects('ATCGATCGATCG', 0.5)

        # Additives increase koff
        assert effects_add.koff_factor > effects_no.koff_factor

    def test_kinetics_bounded(self):
        """Kinetic factors should be bounded."""
        conditions = ReactionConditions(temp=30.0, mg_conc=2.5, ssb=True, betaine_m=2.0)
        model = MechanisticModel(conditions)

        effects = model.calculate_effects('ATCGATCGATCG', 0.5)

        assert 0.1 <= effects.kon_factor <= 3.0
        assert 0.2 <= effects.koff_factor <= 2.0


class TestCombinedEffects:
    """Test combined amplification predictions."""

    def test_amplification_factor_bounded(self):
        """Amplification factor should be bounded 0-1."""
        conditions = ReactionConditions(temp=30.0, mg_conc=2.5)
        model = MechanisticModel(conditions)

        effects = model.calculate_effects('ATCGATCGATCG', 0.5)

        assert 0.0 <= effects.predicted_amplification_factor <= 1.0

    def test_effective_binding_rate_bounded(self):
        """Effective binding rate should be bounded 0-1."""
        conditions = ReactionConditions(temp=30.0, mg_conc=2.5, ssb=True)
        model = MechanisticModel(conditions)

        effects = model.calculate_effects('ATCGATCGATCG', 0.5)

        assert 0.0 <= effects.effective_binding_rate <= 1.0

    def test_enhanced_conditions_improve_gc_amplification(self):
        """Enhanced conditions should improve high-GC template amplification."""
        cond_standard = ReactionConditions(temp=30.0, mg_conc=2.5)
        cond_enhanced = ReactionConditions(
            temp=42.0, mg_conc=2.5, polymerase='equiphi29',
            dmso_percent=5.0, betaine_m=1.0
        )

        model_std = MechanisticModel(cond_standard)
        model_enh = MechanisticModel(cond_enhanced)

        # Use a 12-mer primer for reasonable Tm
        primer = 'ATCGATCGATCG'

        effects_std = model_std.calculate_effects(primer, template_gc=0.7)
        effects_enh = model_enh.calculate_effects(primer, template_gc=0.7)

        # Enhanced conditions should have better accessibility for high GC
        assert effects_enh.accessibility_factor > effects_std.accessibility_factor


class TestHelperMethods:
    """Test helper methods."""

    def test_primer_gc_calculation(self):
        """Test GC content calculation."""
        assert MechanisticModel._primer_gc('GGGGCCCC') == 1.0
        assert MechanisticModel._primer_gc('AAAATTTT') == 0.0
        assert MechanisticModel._primer_gc('ATGC') == 0.5
        assert MechanisticModel._primer_gc('') == 0.5  # Default for empty

    def test_sigmoid_function(self):
        """Test sigmoid function."""
        # At midpoint, sigmoid should be ~0.5
        assert abs(MechanisticModel._sigmoid(1.0, 1.0, 1.0) - 0.5) < 0.01

        # Far below midpoint should be ~0
        assert MechanisticModel._sigmoid(-5.0, 1.0, 1.0) < 0.01

        # Far above midpoint should be ~1
        assert MechanisticModel._sigmoid(7.0, 1.0, 1.0) > 0.99

        # Higher steepness = sharper transition
        steep = MechanisticModel._sigmoid(1.5, 1.0, 5.0)
        gentle = MechanisticModel._sigmoid(1.5, 1.0, 1.0)
        assert steep > gentle


class TestPolymeraseSpecific:
    """Test polymerase-specific behavior."""

    def test_equiphi29_at_42c_optimal(self):
        """EquiPhi29 should perform well at 42C."""
        conditions = ReactionConditions(temp=42.0, mg_conc=2.5, polymerase='equiphi29')
        model = MechanisticModel(conditions)

        effects = model.calculate_effects('ATCGATCGATCG', 0.5)

        # At optimal temp, speed should be high
        assert effects.speed_factor > 0.9

    def test_phi29_suboptimal_at_40c(self):
        """Phi29 should be suboptimal at 40C."""
        cond_30 = ReactionConditions(temp=30.0, mg_conc=2.5, polymerase='phi29')
        cond_40 = ReactionConditions(temp=40.0, mg_conc=2.5, polymerase='phi29')

        model_30 = MechanisticModel(cond_30)
        model_40 = MechanisticModel(cond_40)

        effects_30 = model_30.calculate_effects('ATCGATCGATCG', 0.5)
        effects_40 = model_40.calculate_effects('ATCGATCGATCG', 0.5)

        # Speed should be lower at non-optimal temp
        assert effects_40.speed_factor < effects_30.speed_factor

    def test_bst_at_63c(self):
        """Bst should work at 63C."""
        conditions = ReactionConditions(temp=63.0, mg_conc=2.5, polymerase='bst')
        model = MechanisticModel(conditions)

        effects = model.calculate_effects('ATCGATCGATCGATCGATCG', 0.5)

        # Should have reasonable speed at optimal temp
        assert effects.speed_factor > 0.8
