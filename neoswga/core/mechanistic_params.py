"""
Literature-based parameters for mechanistic SWGA model.

All parameters are derived from published literature with references.
This module provides the constants used by MechanisticModel to calculate
the effects of reaction conditions on SWGA performance.

Four-pathway model:
1. Tm modification - how additives affect primer-template stability
2. Secondary structure accessibility - template melting effects
3. Enzyme activity - polymerase performance modifiers
4. Binding kinetics - association/dissociation rate effects

References:
    - Varadaraj & Skinner (1994) Gene 140:1-5 (DMSO)
    - Blake & Delcourt (1996) NAR 24:2095-2103 (formamide)
    - Spiess et al. (2004) Biotechniques 36:732-736 (trehalose)
    - Cheng et al. (1994) PNAS 91:5695-5699 (ethanol)
    - Rees et al. (1993) Biochemistry 32:137-144 (betaine)
    - Hutton (1977) NAR 4:3537-3555 (urea)
    - Melchior & von Hippel (1973) PNAS 70:298-302 (TMAC)
    - Blanco et al. (1989) JBC 264:8935-8940 (phi29 processivity)
    - Rahman et al. (2014) PLoS One 9:e112515 (Mg optimization)
"""

from typing import Dict, Any, Tuple


# =============================================================================
# Additive Tm Correction Parameters (Arrhenius-based)
# =============================================================================
# Literature-based parameters for temperature-dependent Tm corrections.
# Uses Arrhenius kinetics: dTm(T) = dTm_ref * exp(-Ea/R * (1/T - 1/T_ref))
#
# References:
#   - Chester & Marshak (1993) Anal Biochem 209:284-290 (DMSO multi-temp)
#   - Rees et al. (1993) Biochemistry 32:137-144 (betaine)
#   - Henke et al. (1997) NAR 25:3957-3958 (betaine PCR)
#   - Blake & Delcourt (1996) NAR 24:2095-2103 (formamide)
#   - McConaughy et al. (1969) Biochemistry 8:3289-3295 (formamide)
#   - Spiess et al. (2004) Biotechniques 36:732-736 (trehalose)
#   - Lesnick & Bhalla (1995) NAR 23:4665-4666 (urea)
#   - Hutton (1977) NAR 4:3537-3555 (urea)
#   - Melchior & von Hippel (1973) PNAS 70:298-302 (TMAC)
#   - Cheng et al. (1994) PNAS 91:5695-5699 (ethanol)
# =============================================================================

ADDITIVE_TM_PARAMS: Dict[str, Dict[str, Any]] = {
    'dmso': {
        'ref_coef': -0.55,            # C per % at T_ref (Chester 1993)
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 2500.0,  # J/mol (estimated from Chester 1993)
        'max_concentration': 10.0,    # %
        'gc_dependent': False,
        'description': 'Destabilizes AT base pairs, reduces secondary structure',
    },
    'betaine': {
        # Published value: -1.0 C/M (Rees 1993), -1.3 C/M (Henke 1997 PCR context).
        # Using -1.2 C/M as a weighted average of the two studies. The Rees
        # (1993) value was measured on long DNA, while Henke (1997) measured on
        # short PCR amplicons closer to SWGA primer lengths. Neither study was
        # performed at SWGA-relevant concentrations (0.5-2.5 M), so the value
        # represents a literature consensus rather than a single source.
        'ref_coef': -1.2,             # C per M at T_ref
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 1800.0,  # J/mol (estimated from Rees 1993 multi-temp data)
        'max_concentration': 2.5,     # M
        'gc_dependent': True,
        'gc_equalization_conc': 5.2,  # M for full GC independence (Rees 1993)
        'description': 'Equalizes AT/GC stability, enables longer primers',
    },
    'formamide': {
        'ref_coef': -0.65,            # C per % (Blake 1996)
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 3000.0,  # J/mol (estimated from McConaughy 1969)
        'max_concentration': 10.0,    # %
        'gc_dependent': False,
        'description': 'Destabilizes hydrogen bonding',
    },
    'trehalose': {
        # Published value: Spiess et al. (2004) reported Tm depression of
        # ~2-4 C/M depending on sequence context. We use -3.0 C/M as the
        # midpoint of their observed range. Their measurements used real-time
        # PCR with short amplicons at 0.2-0.8 M trehalose, which is within
        # the practical SWGA concentration range.
        'ref_coef': -3.0,             # C per M (midpoint of Spiess 2004 range)
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 1500.0,  # J/mol (estimated; no multi-temp data available)
        'max_concentration': 1.0,     # M
        'gc_dependent': False,
        'description': 'Stabilizes proteins, modifies water structure',
    },
    'urea': {
        # Published values: Hutton (1977) reported -5.0 C/M for long genomic
        # DNA at high urea concentrations (4-8 M). Lesnick & Bhalla (1995)
        # measured -2.0 to -3.0 C/M for short oligonucleotides at 0.5-2 M,
        # which is the concentration range relevant to SWGA. The discrepancy
        # is likely due to cooperative denaturation effects in long DNA that
        # do not apply to short primer-template duplexes. We use -2.5 C/M
        # as the midpoint of the Lesnick (1995) range.
        'ref_coef': -2.5,             # C per M (Lesnick 1995, short oligos)
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 2000.0,  # J/mol (estimated from Hutton 1977 multi-temp data)
        'max_concentration': 2.0,     # M
        'gc_dependent': True,
        'gc_preference': 1.3,         # 30% stronger effect on GC-rich sequences
        'description': 'Preferentially destabilizes GC base pairs',
    },
    'tmac': {
        'ref_coef': -0.5,             # C per M uniform component (minimal)
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 1000.0,  # J/mol (estimated)
        'max_concentration': 0.1,     # M (practical SWGA range)
        'gc_dependent': True,
        'gc_equalization_conc': 3.0,  # M for full GC independence (Melchior 1973)
        'description': 'Equalizes AT/GC Tm, isostabilizing agent',
    },
    'ethanol': {
        'ref_coef': -0.4,             # C per % (Cheng 1994)
        'ref_temp': 310.15,           # K (37C)
        'activation_energy': 2200.0,  # J/mol (estimated)
        'max_concentration': 5.0,     # %
        'gc_dependent': False,
        'description': 'Reduces secondary structure formation',
    },
}


MECHANISTIC_MODEL_PARAMS: Dict[str, Dict[str, Any]] = {
    # Pathway 1: Tm modification
    # How additives affect primer melting temperature
    # NOTE: Use ADDITIVE_TM_PARAMS for Arrhenius-based corrections
    'tm': {
        # Uniform Tm corrections (C per unit concentration).
        # See ADDITIVE_TM_PARAMS above for detailed provenance of each value.
        'dmso_coef': -0.55,           # C per %, Chester 1993, Varadaraj 1994
        'formamide_coef': -0.65,      # C per %, Blake 1996
        'trehalose_coef': -3.0,       # C per M, midpoint of Spiess 2004 range
        'ethanol_coef': -0.4,         # C per %, Cheng 1994
        'betaine_uniform_coef': -1.2, # C per M, avg of Rees 1993 / Henke 1997
        'urea_coef': -2.5,            # C per M, Lesnick 1995 (short oligos)
        'tmac_uniform_coef': -0.5,    # C per M, Melchior 1973

        # GC normalization sigmoid parameters
        # Betaine reaches full GC equalization at ~5.2M (Rees 1993)
        # Using sigmoid for more accurate dose-response
        'betaine_gc_midpoint': 1.5,   # M for sigmoid center
        'betaine_gc_steepness': 1.5,  # sigmoid steepness

        # TMAC reaches full equalization at ~3M (Melchior 1973)
        'tmac_gc_midpoint': 1.0,      # M for sigmoid center
        'tmac_gc_steepness': 2.0,     # sigmoid steepness
    },

    # Pathway 2: Secondary structure accessibility
    # How additives affect template secondary structure
    'structure': {
        # DMSO melts secondary structure (saturating effect)
        'dmso_melt_rate': 0.5,        # rate constant for 1-exp(-rate*conc)
        'dmso_max_effect': 0.6,       # maximum structure reduction

        # Betaine also helps melt structure
        'betaine_melt_rate': 0.3,
        'betaine_max_effect': 0.4,

        # Temperature effect on structure
        'temp_structure_coef': 0.02,  # per degree above 25C

        # GC content effect on baseline structure
        'gc_structure_coef': 1.0,     # penalty for GC > 50%
    },

    # Pathway 3: Enzyme activity modifiers
    # How conditions affect polymerase performance
    'enzyme': {
        # Phi29 polymerase characteristics
        'phi29': {
            'optimal_temp': 30.0,     # C
            'processivity': 70000,    # bp, Blanco 1989
            'processivity_step': 0.99857,  # Per 100bp step for stochastic model
            'extension_rate': 150,    # nt/s at optimal temp (median of 100-170 range)
            'dmso_threshold': 5.0,    # % before steep inhibition
            'dmso_mild_coef': 0.02,   # activity reduction per % below threshold
            'dmso_steep_coef': 0.12,  # activity reduction per % above threshold
        },
        # EquiPhi29 (thermostable variant)
        'equiphi29': {
            'optimal_temp': 42.0,
            'processivity': 80000,
            'processivity_step': 0.99875,  # Per 100bp step for stochastic model
            'extension_rate': 200,    # nt/s (thermally enhanced)
            'dmso_threshold': 4.0,    # slightly less tolerant
            'dmso_mild_coef': 0.025,
            'dmso_steep_coef': 0.15,
        },
        # Bst polymerase
        'bst': {
            'optimal_temp': 63.0,
            'processivity': 2000,
            'processivity_step': 0.9512,   # Per 100bp step for stochastic model
            'extension_rate': 100,
            'dmso_threshold': 3.0,
            'dmso_mild_coef': 0.03,
            'dmso_steep_coef': 0.18,
        },
        # Klenow fragment
        'klenow': {
            'optimal_temp': 37.0,
            'processivity': 10000,
            'processivity_step': 0.99005,  # Per 100bp step for stochastic model
            'extension_rate': 50,
            'dmso_threshold': 5.0,
            'dmso_mild_coef': 0.02,
            'dmso_steep_coef': 0.10,
        },

        # Betaine effect: enhancement at low conc, inhibition at high
        'betaine_peak': 1.0,          # M for maximum enhancement
        'betaine_enhancement': 0.12,  # max stability enhancement
        'betaine_inhibition_start': 1.5,  # M where inhibition begins
        'betaine_inhibition_coef': 0.15,  # processivity reduction per M excess

        # Formamide is always inhibitory
        'formamide_coef': 0.06,       # activity reduction per %

        # Mg2+ optimum
        'mg_optimal': 2.5,            # mM optimal concentration
        'mg_low_threshold': 1.0,      # mM below which activity drops
        'mg_high_threshold': 6.0,     # mM above which activity drops

        # Glycerol: stabilizes enzyme but slows it
        'glycerol_stability': 0.02,   # stability increase per %
        'glycerol_speed_penalty': 0.02,  # speed reduction per %

        # Temperature activity coefficient
        'temp_activity_coef': 0.03,   # activity reduction per C from optimal
    },

    # Pathway 4: Binding kinetics
    # How conditions affect primer binding rates
    'kinetics': {
        # Optimal Tm is above reaction temp for stable binding
        'optimal_delta_t': 7.0,       # C above reaction temp
        'delta_t_width': 5.0,         # Gaussian width for Tm effect

        # Betaine effects on kinetics
        'betaine_kon_boost': 0.10,    # kon increase per M
        'betaine_koff_penalty': 0.05, # koff increase per M

        # DMSO effects on kinetics
        'dmso_kon_boost': 0.02,       # kon increase per %
        'dmso_koff_penalty': 0.03,    # koff increase per %

        # SSB (single-stranded binding protein) effect
        'ssb_kon_multiplier': 2.0,    # SSB doubles kon

        # Temperature effect on koff
        'koff_temp_coef': 0.1,        # koff increase per C above optimal
    },

    # Additive interactions (synergistic/antagonistic)
    'interactions': {
        # DMSO chelates Mg2+ at high concentrations
        'dmso_mg_chelation': 0.08,    # mM Mg per % DMSO
        'dmso_mg_threshold': 3.0,     # % DMSO where chelation starts

        # Betaine + trehalose synergy for enzyme stability
        'betaine_trehalose_synergy': 0.15,

        # DMSO + formamide antagonism (both destabilize)
        'dmso_formamide_antagonism': 0.03,
    },
}


# Application profiles for automatic set size optimization
#
# The coverage vs fg/bg ratio tradeoff is NOT a simple inverse relationship.
# Good primer selection can improve BOTH coverage AND specificity by choosing
# primers with high fg/bg ratios that also bind unique target regions.
#
# The Pareto frontier approach generates sets of different sizes and shows
# the achievable (coverage, fg_bg_ratio) combinations. Application profiles
# determine where to select on this frontier:
#
# - priority='coverage': maximize coverage while meeting fg_bg_ratio constraint
# - priority='specificity': maximize fg_bg_ratio while meeting coverage needs
# - priority='balanced': find the knee of the curve (best tradeoff)
#
# Note: 'min_specificity' is DEPRECATED, use 'default_min_fg_bg_ratio' instead.
# The fg/bg ratio is a more meaningful metric than abstract "specificity".
APPLICATION_PROFILES: Dict[str, Dict[str, Any]] = {
    'discovery': {
        'priority': 'coverage',
        'default_target_coverage': 0.90,
        'default_min_fg_bg_ratio': 2.0,
        'typical_size': (10, 15),
        'description': 'Pathogen discovery - maximize sensitivity',
        # Legacy fields for backward compatibility
        'target_coverage': 0.90,
        'min_specificity': 0.60,
    },
    'clinical': {
        'priority': 'specificity',
        'default_target_coverage': 0.70,
        'default_min_fg_bg_ratio': 10.0,
        'typical_size': (6, 10),
        'description': 'Clinical diagnostics - minimize false positives',
        # Legacy fields for backward compatibility
        'target_coverage': 0.70,
        'min_specificity': 0.90,
    },
    'enrichment': {
        'priority': 'balanced',
        'default_target_coverage': 0.80,
        'default_min_fg_bg_ratio': 5.0,
        'typical_size': (8, 12),
        'description': 'Sequencing enrichment - balanced approach',
        # Legacy fields for backward compatibility
        'target_coverage': 0.80,
        'min_specificity': 0.75,
    },
    'metagenomics': {
        'priority': 'coverage',
        'default_target_coverage': 0.95,
        'default_min_fg_bg_ratio': 1.5,
        'typical_size': (15, 20),
        'description': 'Metagenomics - capture diversity',
        # Legacy fields for backward compatibility
        'target_coverage': 0.95,
        'min_specificity': 0.50,
    },
}


def get_polymerase_params(polymerase: str) -> Dict[str, Any]:
    """
    Get parameters for a specific polymerase.

    Args:
        polymerase: Polymerase name ('phi29', 'equiphi29', 'bst', 'klenow')

    Returns:
        Dictionary of polymerase-specific parameters

    Raises:
        ValueError: If polymerase is not recognized
    """
    poly = polymerase.lower()
    enzyme_params = MECHANISTIC_MODEL_PARAMS['enzyme']

    if poly not in enzyme_params:
        available = [k for k in enzyme_params.keys()
                     if isinstance(enzyme_params[k], dict) and 'optimal_temp' in enzyme_params[k]]
        raise ValueError(
            f"Unknown polymerase '{polymerase}'. "
            f"Available: {', '.join(available)}"
        )

    return enzyme_params[poly]


def get_application_profile(application: str) -> Dict[str, Any]:
    """
    Get application profile for set size optimization.

    Args:
        application: Application name ('discovery', 'clinical', 'enrichment', 'metagenomics')

    Returns:
        Dictionary with target_coverage, min_specificity, typical_size, description

    Raises:
        ValueError: If application is not recognized
    """
    if application not in APPLICATION_PROFILES:
        raise ValueError(
            f"Unknown application '{application}'. "
            f"Available: {', '.join(APPLICATION_PROFILES.keys())}"
        )

    return APPLICATION_PROFILES[application]


def list_applications() -> Dict[str, str]:
    """
    List available application profiles with descriptions.

    Returns:
        Dictionary mapping application name to description
    """
    return {name: profile['description'] for name, profile in APPLICATION_PROFILES.items()}


def list_polymerases() -> Dict[str, Tuple[float, int]]:
    """
    List available polymerases with their optimal temp and processivity.

    Returns:
        Dictionary mapping polymerase name to (optimal_temp, processivity)
    """
    enzyme_params = MECHANISTIC_MODEL_PARAMS['enzyme']
    result = {}

    for name, params in enzyme_params.items():
        if isinstance(params, dict) and 'optimal_temp' in params:
            result[name] = (params['optimal_temp'], params['processivity'])

    return result


def get_additive_tm_params(additive: str) -> Dict[str, Any]:
    """
    Get Arrhenius-based Tm correction parameters for an additive.

    Args:
        additive: Additive name (dmso, betaine, formamide, trehalose,
                  urea, tmac, ethanol)

    Returns:
        Dictionary with ref_coef, ref_temp, activation_energy, max_concentration,
        gc_dependent, and additive-specific parameters

    Raises:
        ValueError: If additive is not recognized
    """
    additive_lower = additive.lower()
    if additive_lower not in ADDITIVE_TM_PARAMS:
        available = list(ADDITIVE_TM_PARAMS.keys())
        raise ValueError(
            f"Unknown additive '{additive}'. "
            f"Available: {', '.join(available)}"
        )
    return ADDITIVE_TM_PARAMS[additive_lower]


def list_additives() -> Dict[str, str]:
    """
    List available additives with descriptions.

    Returns:
        Dictionary mapping additive name to description
    """
    return {name: params['description'] for name, params in ADDITIVE_TM_PARAMS.items()}
