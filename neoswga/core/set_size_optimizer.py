"""
Automatic primer set size optimization.

Determines optimal number of primers based on:
- Application profile (discovery, clinical, enrichment, metagenomics)
- Reaction conditions (affects primer effectiveness)
- Coverage/specificity tradeoff

The set size recommendation is based on first-principles modeling of
SWGA coverage combined with empirical application profiles.

Usage:
    from neoswga.core.set_size_optimizer import recommend_set_size
    from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
    from neoswga.core.reaction_conditions import ReactionConditions

    # Get mechanistic effects for a representative primer
    conditions = ReactionConditions(temp=42.0, polymerase='equiphi29', mg_conc=2.5)
    model = MechanisticModel(conditions)
    effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)

    # Get recommendation
    recommendation = recommend_set_size(
        application='clinical',
        genome_length=1_000_000,
        primer_length=12,
        mech_effects=effects,
    )
    print(f"Recommended size: {recommendation['recommended_size']}")
"""

import math
import logging
from typing import List, Tuple, Dict, Optional, Any, TYPE_CHECKING

import numpy as np

from neoswga.core.mechanistic_params import APPLICATION_PROFILES, get_application_profile

if TYPE_CHECKING:
    from neoswga.core.mechanistic_model import MechanisticEffects

logger = logging.getLogger(__name__)


def estimate_optimal_set_size(
    genome_length: int,
    primer_length: int,
    target_coverage: float,
    processivity: int,
    mech_effects: 'MechanisticEffects'
) -> int:
    """
    Estimate optimal set size from first principles.

    Uses a probabilistic model of genome coverage:
    Coverage = 1 - exp(-N * sites_per_primer * effective_processivity / genome_length)

    Args:
        genome_length: Target genome length in bp
        primer_length: Primer length in bp (k-mer size)
        target_coverage: Desired genome coverage fraction (0-1)
        processivity: Base polymerase processivity in bp
        mech_effects: MechanisticEffects from the model

    Returns:
        Estimated optimal number of primers

    Example:
        >>> from neoswga.core.mechanistic_model import MechanisticEffects
        >>> effects = MechanisticEffects(
        ...     tm_correction=0, effective_tm=35.0,
        ...     accessibility_factor=1.0, processivity_factor=1.0,
        ...     speed_factor=1.0, stability_factor=1.0,
        ...     kon_factor=1.0, koff_factor=1.0,
        ...     effective_binding_rate=0.8, effective_extension_rate=1.0,
        ...     predicted_amplification_factor=0.8,
        ... )
        >>> estimate_optimal_set_size(1_000_000, 10, 0.8, 70000, effects)
        8
    """
    # Expected binding sites per primer (random sequence model)
    # For a k-mer in a random sequence of length L, expected count is:
    # E[sites] = (L - k + 1) * (1/4)^k * 2 (both strands)
    # Simplified: ~ L * 4^(-k) * 2
    expected_sites = (4 ** (-primer_length)) * genome_length * 2

    # Effective processivity accounting for reaction conditions
    effective_proc = processivity * mech_effects.processivity_factor

    # Coverage contribution per primer
    # Each primer contributes coverage proportional to:
    # sites * effective_processivity * binding_rate / genome_length
    coverage_per_primer = (
        expected_sites * effective_proc * mech_effects.effective_binding_rate
        / genome_length
    )

    if coverage_per_primer <= 0:
        logger.warning("Coverage per primer is zero or negative, using default size")
        return 20

    # Solve for N in: target_coverage = 1 - exp(-N * coverage_per_primer)
    # => N = -ln(1 - target_coverage) / coverage_per_primer
    try:
        optimal_n = -math.log(1 - target_coverage) / coverage_per_primer
    except (ValueError, ZeroDivisionError):
        logger.warning("Could not compute optimal N, using default size")
        return 15

    # Apply overlap correction factor
    # Primers don't cover genome independently; there's overlap
    # Empirical correction: multiply by ~1.3
    optimal_n *= 1.3

    # Bound to reasonable range
    return max(4, min(20, int(optimal_n + 0.5)))


def recommend_set_size(
    application: str,
    genome_length: int,
    primer_length: int,
    mech_effects: 'MechanisticEffects',
    processivity: int = 70000
) -> Dict[str, Any]:
    """
    Recommend primer set size based on application profile.

    Combines first-principles estimation with application-specific
    requirements for coverage and specificity.

    Args:
        application: Application profile name
            - 'discovery': Maximize sensitivity for pathogen discovery
            - 'clinical': Minimize false positives for diagnostics
            - 'enrichment': Balanced for sequencing enrichment
            - 'metagenomics': Capture diversity
        genome_length: Target genome length in bp
        primer_length: Primer length in bp
        mech_effects: MechanisticEffects from mechanistic model
        processivity: Base polymerase processivity (default 70000 for phi29)

    Returns:
        Dictionary with:
            - recommended_size: Recommended number of primers
            - size_range: (min, max) typical range for this application
            - target_coverage: Target coverage fraction
            - min_specificity: Minimum specificity requirement
            - rationale: Description of the application
            - base_estimate: Raw estimate before application adjustment

    Example:
        >>> recommendation = recommend_set_size('clinical', 1_000_000, 12, effects)
        >>> print(f"Use {recommendation['recommended_size']} primers")
    """
    # Get application profile
    try:
        profile = get_application_profile(application)
    except ValueError:
        logger.warning(f"Unknown application '{application}', using 'enrichment'")
        application = 'enrichment'
        profile = get_application_profile(application)

    # Get base estimate from first principles
    base_estimate = estimate_optimal_set_size(
        genome_length=genome_length,
        primer_length=primer_length,
        target_coverage=profile['target_coverage'],
        processivity=processivity,
        mech_effects=mech_effects,
    )

    # Adjust for specificity requirement
    # Higher specificity -> fewer primers (more selective)
    # Lower specificity -> more primers (broader coverage)
    if profile['min_specificity'] > 0.80:
        # High specificity applications need fewer, more selective primers
        adjusted = int(base_estimate * 0.8)
    elif profile['min_specificity'] < 0.60:
        # Low specificity applications can use more primers
        adjusted = int(base_estimate * 1.2)
    else:
        adjusted = base_estimate

    # Constrain to typical range for this application
    min_size, max_size = profile['typical_size']
    recommended = max(min_size, min(max_size, adjusted))

    return {
        'recommended_size': recommended,
        'size_range': profile['typical_size'],
        'target_coverage': profile['target_coverage'],
        'min_specificity': profile['min_specificity'],
        'rationale': profile['description'],
        'base_estimate': base_estimate,
        'application': application,
    }


def find_optimal_size_by_elbow(
    optimizer,
    candidates: List[str],
    min_size: int = 4,
    max_size: int = 15,
    verbose: bool = True
) -> Tuple[int, List[Dict[str, Any]]]:
    """
    Find optimal set size using elbow method on coverage curve.

    Runs optimization at multiple set sizes and finds the "elbow"
    where adding more primers yields diminishing returns.

    Args:
        optimizer: NetworkOptimizer instance (or compatible optimizer)
        candidates: List of candidate primer sequences
        min_size: Minimum set size to evaluate
        max_size: Maximum set size to evaluate
        verbose: Whether to log progress

    Returns:
        Tuple of (optimal_size, results_list)
        - optimal_size: Recommended set size at elbow
        - results_list: List of dicts with size, primers, coverage, score

    Note:
        This method is computationally expensive as it runs full
        optimization at each size. Use recommend_set_size() for
        a faster estimate.
    """
    results = []

    for size in range(min_size, max_size + 1):
        if verbose:
            logger.info(f"Evaluating set size {size}...")

        # Run optimization at this size
        primers = optimizer.optimize_greedy(candidates, num_primers=size)
        score = optimizer.score_primer_set(primers)

        result = {
            'size': size,
            'primers': primers,
            'coverage': score.get('target_coverage', score.get('enrichment', 0)),
            'enrichment': score.get('enrichment', 0),
            'score': score,
        }
        results.append(result)

        if verbose:
            logger.info(f"  Size {size}: coverage={result['coverage']:.1%}, "
                       f"enrichment={result['enrichment']:.1f}x")

    # Find elbow using curvature analysis
    if len(results) < 3:
        return results[len(results) // 2]['size'], results

    sizes = np.array([r['size'] for r in results])
    coverages = np.array([r['coverage'] for r in results])

    # Normalize to [0, 1] range
    size_range = sizes.max() - sizes.min()
    coverage_range = coverages.max() - coverages.min()

    if size_range == 0 or coverage_range < 1e-10:
        # No variation, return middle
        return results[len(results) // 2]['size'], results

    sizes_norm = (sizes - sizes.min()) / size_range
    coverages_norm = (coverages - coverages.min()) / coverage_range

    # Compute curvature using finite differences
    # First derivative
    d1 = np.gradient(coverages_norm, sizes_norm)
    # Second derivative
    d2 = np.gradient(d1, sizes_norm)

    # Curvature formula: |d2| / (1 + d1^2)^1.5
    curvature = np.abs(d2) / np.power(1 + d1**2, 1.5)

    # Find maximum curvature (excluding endpoints which can be noisy)
    if len(curvature) > 2:
        # Exclude first and last points
        interior_curvature = curvature[1:-1]
        elbow_idx = np.argmax(interior_curvature) + 1
    else:
        elbow_idx = len(results) // 2

    optimal_size = results[elbow_idx]['size']

    if verbose:
        logger.info(f"Elbow detected at size {optimal_size}")

    return optimal_size, results


def get_size_recommendation_summary(
    application: str,
    genome_length: int,
    primer_length: int,
    mech_effects: 'MechanisticEffects',
    processivity: int = 70000
) -> str:
    """
    Get a human-readable summary of set size recommendation.

    Args:
        application: Application profile name
        genome_length: Target genome length in bp
        primer_length: Primer length in bp
        mech_effects: MechanisticEffects from mechanistic model
        processivity: Base polymerase processivity

    Returns:
        Formatted string with recommendation details
    """
    rec = recommend_set_size(
        application=application,
        genome_length=genome_length,
        primer_length=primer_length,
        mech_effects=mech_effects,
        processivity=processivity,
    )

    lines = [
        f"Set Size Recommendation for '{application}' Application",
        "=" * 50,
        f"",
        f"Recommended primers: {rec['recommended_size']}",
        f"Typical range: {rec['size_range'][0]}-{rec['size_range'][1]}",
        f"",
        f"Application: {rec['rationale']}",
        f"Target coverage: {rec['target_coverage']:.0%}",
        f"Min specificity: {rec['min_specificity']:.0%}",
        f"",
        f"Base estimate (first principles): {rec['base_estimate']}",
        f"",
        f"Input parameters:",
        f"  Genome length: {genome_length:,} bp",
        f"  Primer length: {primer_length} bp",
        f"  Processivity: {processivity:,} bp",
        f"  Effective binding: {mech_effects.effective_binding_rate:.2f}",
        f"  Processivity factor: {mech_effects.processivity_factor:.2f}",
    ]

    return "\n".join(lines)


def create_baseline_effects() -> 'MechanisticEffects':
    """
    Create baseline MechanisticEffects for quick estimates.

    Returns MechanisticEffects with all factors at 1.0 (optimal conditions).
    Use this when you don't have specific reaction conditions.

    Returns:
        MechanisticEffects with baseline values
    """
    from neoswga.core.mechanistic_model import MechanisticEffects

    return MechanisticEffects(
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


def quick_size_estimate(
    application: str,
    genome_length: int,
    primer_length: int = 10,
    processivity: int = 70000
) -> int:
    """
    Quick set size estimate without full mechanistic model.

    Convenience function for rapid estimates when you don't need
    detailed reaction condition modeling.

    Args:
        application: Application profile name
        genome_length: Target genome length in bp
        primer_length: Primer length (default 10)
        processivity: Polymerase processivity (default 70000 for phi29)

    Returns:
        Recommended number of primers

    Example:
        >>> quick_size_estimate('clinical', 1_000_000)
        8
        >>> quick_size_estimate('discovery', 5_000_000)
        12
    """
    effects = create_baseline_effects()
    rec = recommend_set_size(
        application=application,
        genome_length=genome_length,
        primer_length=primer_length,
        mech_effects=effects,
        processivity=processivity,
    )
    return rec['recommended_size']
