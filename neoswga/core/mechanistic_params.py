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


MECHANISTIC_MODEL_PARAMS: Dict[str, Dict[str, Any]] = {
    # Pathway 1: Tm modification
    # How additives affect primer melting temperature
    'tm': {
        # Uniform Tm corrections (C per unit concentration)
        'dmso_coef': -0.6,            # C per %, Varadaraj 1994
        'formamide_coef': -0.65,      # C per %, Blake 1996
        'trehalose_coef': -5.0,       # C per M, Spiess 2004
        'ethanol_coef': -0.5,         # C per %, Cheng 1994
        'betaine_uniform_coef': -0.5, # C per M, Rees 1993
        'urea_coef': -5.0,            # C per M, Hutton 1977
        'tmac_uniform_coef': -10.0,   # C per M, Melchior 1973

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
            'extension_rate': 167,    # nt/s at optimal temp
            'dmso_threshold': 5.0,    # % before steep inhibition
            'dmso_mild_coef': 0.02,   # activity reduction per % below threshold
            'dmso_steep_coef': 0.12,  # activity reduction per % above threshold
        },
        # EquiPhi29 (thermostable variant)
        'equiphi29': {
            'optimal_temp': 42.0,
            'processivity': 80000,
            'extension_rate': 200,
            'dmso_threshold': 4.0,    # slightly less tolerant
            'dmso_mild_coef': 0.025,
            'dmso_steep_coef': 0.15,
        },
        # Bst polymerase
        'bst': {
            'optimal_temp': 63.0,
            'processivity': 2000,
            'extension_rate': 100,
            'dmso_threshold': 3.0,
            'dmso_mild_coef': 0.03,
            'dmso_steep_coef': 0.18,
        },
        # Klenow fragment
        'klenow': {
            'optimal_temp': 37.0,
            'processivity': 10000,
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
APPLICATION_PROFILES: Dict[str, Dict[str, Any]] = {
    'discovery': {
        'target_coverage': 0.90,
        'min_specificity': 0.60,
        'typical_size': (10, 15),
        'description': 'Pathogen discovery - maximize sensitivity',
    },
    'clinical': {
        'target_coverage': 0.70,
        'min_specificity': 0.90,
        'typical_size': (6, 10),
        'description': 'Clinical diagnostics - minimize false positives',
    },
    'enrichment': {
        'target_coverage': 0.80,
        'min_specificity': 0.75,
        'typical_size': (8, 12),
        'description': 'Sequencing enrichment - balanced approach',
    },
    'metagenomics': {
        'target_coverage': 0.95,
        'min_specificity': 0.50,
        'typical_size': (15, 20),
        'description': 'Metagenomics - capture diversity',
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
