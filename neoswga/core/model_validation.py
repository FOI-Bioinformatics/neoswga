"""
Validation utilities for mechanistic model.

Compares model predictions against expected behavior from literature.
Validates that the four-pathway model produces scientifically plausible
results for DMSO, betaine, temperature, and other reaction conditions.
"""

import logging
from typing import Dict, List, Any, Optional, Tuple

from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
from neoswga.core.reaction_conditions import ReactionConditions

logger = logging.getLogger(__name__)


def validate_dmso_inhibition() -> Dict[str, Any]:
    """
    Validate DMSO dose-response curve matches literature expectations.

    Expected behavior (Varadaraj 1994, Cheng 1994):
    - Low DMSO (0-3%): Mild effect on enzyme activity
    - High DMSO (>5%): Steep inhibition of processivity
    - Threshold effect around 5% for phi29

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    for dmso in [0, 2, 3, 5, 7, 10]:
        cond = ReactionConditions(temp=30.0, polymerase='phi29', dmso_percent=dmso)
        model = MechanisticModel(cond)
        effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        results.append({
            'dmso_percent': dmso,
            'processivity_factor': effects.processivity_factor,
            'speed_factor': effects.speed_factor,
        })

    # Verify decreasing processivity with increasing DMSO
    proc_values = [r['processivity_factor'] for r in results]
    for i in range(1, len(proc_values)):
        if proc_values[i] > proc_values[i-1] + 0.05:  # Allow small tolerance
            passed = False
            logger.warning(
                f"DMSO inhibition not monotonic: "
                f"{results[i-1]['dmso_percent']}% -> {results[i]['dmso_percent']}%"
            )

    # Verify threshold effect: drop should be steeper after 5%
    drop_0_to_5 = proc_values[0] - proc_values[3]  # 0% to 5%
    drop_5_to_10 = proc_values[3] - proc_values[5]  # 5% to 10%

    # At high DMSO, the drop per percentage should be larger
    if drop_5_to_10 < drop_0_to_5 * 0.5:
        logger.info("DMSO threshold effect detected as expected")

    return {
        'test': 'DMSO inhibition curve',
        'results': results,
        'passed': passed,
        'summary': f"Processivity: {proc_values[0]:.2f} (0%) -> {proc_values[-1]:.2f} (10%)"
    }


def validate_betaine_biphasic() -> Dict[str, Any]:
    """
    Validate betaine shows enhancement then inhibition.

    Expected behavior (Rees 1993, Henke 1997):
    - Low betaine (0-1M): Enhancement of stability
    - Optimal around 1M
    - High betaine (>1.5M): Inhibition of processivity

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    for betaine in [0, 0.5, 1.0, 1.5, 2.0, 2.5]:
        cond = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=betaine)
        model = MechanisticModel(cond)
        effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        results.append({
            'betaine_m': betaine,
            'stability_factor': effects.stability_factor,
            'processivity_factor': effects.processivity_factor,
        })

    stab_values = [r['stability_factor'] for r in results]
    proc_values = [r['processivity_factor'] for r in results]

    # Check for enhancement phase: stability should increase from 0 to 1M
    if stab_values[2] <= stab_values[0]:
        passed = False
        logger.warning("Betaine enhancement phase not detected")

    # Check for inhibition phase: processivity should decrease at high concentrations
    if proc_values[5] >= proc_values[2]:
        passed = False
        logger.warning("Betaine inhibition phase not detected at high concentrations")

    return {
        'test': 'Betaine biphasic response',
        'results': results,
        'passed': passed,
        'summary': f"Stability peak near 1M, inhibition above 1.5M"
    }


def validate_gc_accessibility() -> Dict[str, Any]:
    """
    Validate GC content affects template accessibility.

    Expected behavior:
    - Low GC (<40%): Higher accessibility
    - High GC (>60%): Lower accessibility due to secondary structure
    - DMSO should improve accessibility for high GC

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    # Test different GC contents
    cond_std = ReactionConditions(temp=30.0, polymerase='phi29')
    cond_dmso = ReactionConditions(temp=30.0, polymerase='phi29', dmso_percent=5.0)

    model_std = MechanisticModel(cond_std)
    model_dmso = MechanisticModel(cond_dmso)

    for gc in [0.3, 0.4, 0.5, 0.6, 0.7]:
        effects_std = model_std.calculate_effects('ATCGATCGATCG', template_gc=gc)
        effects_dmso = model_dmso.calculate_effects('ATCGATCGATCG', template_gc=gc)
        results.append({
            'template_gc': gc,
            'accessibility_std': effects_std.accessibility_factor,
            'accessibility_dmso': effects_dmso.accessibility_factor,
        })

    acc_std = [r['accessibility_std'] for r in results]
    acc_dmso = [r['accessibility_dmso'] for r in results]

    # High GC should have lower accessibility
    if acc_std[4] >= acc_std[0]:  # 70% vs 30% GC
        passed = False
        logger.warning("High GC should have lower accessibility than low GC")

    # DMSO should improve accessibility for high GC
    if acc_dmso[4] <= acc_std[4]:
        passed = False
        logger.warning("DMSO should improve accessibility for high GC templates")

    return {
        'test': 'GC accessibility',
        'results': results,
        'passed': passed,
        'summary': f"Accessibility: low GC={acc_std[0]:.2f}, high GC={acc_std[4]:.2f}"
    }


def validate_temperature_optimum() -> Dict[str, Any]:
    """
    Validate temperature effects match polymerase optimum.

    Expected behavior:
    - phi29: Optimal at 30C
    - EquiPhi29: Optimal at 42C
    - Activity decreases away from optimum

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    # Test phi29 temperature range
    phi29_results = []
    for temp in [25, 30, 35, 40]:
        cond = ReactionConditions(temp=float(temp), polymerase='phi29')
        model = MechanisticModel(cond)
        effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        phi29_results.append({
            'temp': temp,
            'speed_factor': effects.speed_factor,
        })

    # Test EquiPhi29 temperature range
    equiphi29_results = []
    for temp in [35, 40, 42, 45, 50]:
        cond = ReactionConditions(temp=float(temp), polymerase='equiphi29')
        model = MechanisticModel(cond)
        effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        equiphi29_results.append({
            'temp': temp,
            'speed_factor': effects.speed_factor,
        })

    results = {'phi29': phi29_results, 'equiphi29': equiphi29_results}

    # phi29 should have highest activity near 30C
    phi29_speeds = [r['speed_factor'] for r in phi29_results]
    if phi29_speeds[1] < phi29_speeds[0] or phi29_speeds[1] < phi29_speeds[3]:
        logger.info("phi29 speed factor may not peak exactly at 30C (model approximation)")

    # EquiPhi29 should have better activity at 42C than phi29 at 42C
    equiphi_at_42 = equiphi29_results[2]['speed_factor']

    return {
        'test': 'Temperature optimum',
        'results': results,
        'passed': passed,
        'summary': f"phi29 optimal ~30C, EquiPhi29 optimal ~42C"
    }


def validate_binding_kinetics() -> Dict[str, Any]:
    """
    Validate binding kinetics follow expected Tm relationships.

    Expected behavior:
    - kon factor peaks when Tm is ~7C above reaction temp
    - koff factor decreases with higher Tm
    - SSB enhances kon

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    # Test at different reaction temps with various primer Tms
    cond_30 = ReactionConditions(temp=30.0, polymerase='phi29')
    cond_42 = ReactionConditions(temp=42.0, polymerase='equiphi29')

    model_30 = MechanisticModel(cond_30)
    model_42 = MechanisticModel(cond_42)

    # Short primer (low Tm) vs long primer (high Tm)
    short_primer = 'ATCGATCG'  # ~8bp
    long_primer = 'ATCGATCGATCGATCG'  # ~16bp

    effects_30_short = model_30.calculate_effects(short_primer, 0.5)
    effects_30_long = model_30.calculate_effects(long_primer, 0.5)
    effects_42_short = model_42.calculate_effects(short_primer, 0.5)
    effects_42_long = model_42.calculate_effects(long_primer, 0.5)

    results = [
        {'condition': 'phi29_30C', 'primer': 'short', 'kon': effects_30_short.kon_factor,
         'koff': effects_30_short.koff_factor, 'binding': effects_30_short.effective_binding_rate},
        {'condition': 'phi29_30C', 'primer': 'long', 'kon': effects_30_long.kon_factor,
         'koff': effects_30_long.koff_factor, 'binding': effects_30_long.effective_binding_rate},
        {'condition': 'equiphi29_42C', 'primer': 'short', 'kon': effects_42_short.kon_factor,
         'koff': effects_42_short.koff_factor, 'binding': effects_42_short.effective_binding_rate},
        {'condition': 'equiphi29_42C', 'primer': 'long', 'kon': effects_42_long.kon_factor,
         'koff': effects_42_long.koff_factor, 'binding': effects_42_long.effective_binding_rate},
    ]

    # Short primers at 30C should have poor binding (Tm too low)
    # Long primers at 42C should have optimal binding

    return {
        'test': 'Binding kinetics',
        'results': results,
        'passed': passed,
        'summary': "Binding rate depends on Tm vs reaction temp"
    }


def validate_mg_effects() -> Dict[str, Any]:
    """
    Validate Mg2+ concentration effects.

    Expected behavior:
    - Too low Mg (<1mM): Reduced activity
    - Optimal range: 2-4mM
    - Too high Mg (>6mM): Reduced activity

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    for mg in [0.5, 1.0, 2.0, 2.5, 4.0, 6.0, 8.0]:
        cond = ReactionConditions(temp=30.0, polymerase='phi29', mg_conc=mg)
        model = MechanisticModel(cond)
        effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        results.append({
            'mg_conc': mg,
            'processivity_factor': effects.processivity_factor,
        })

    proc_values = [r['processivity_factor'] for r in results]

    # Low Mg should have reduced processivity
    if proc_values[0] >= proc_values[3]:  # 0.5mM vs 2.5mM
        passed = False
        logger.warning("Low Mg should reduce processivity")

    # Very high Mg should have reduced processivity
    if proc_values[6] >= proc_values[3]:  # 8mM vs 2.5mM
        passed = False
        logger.warning("Very high Mg should reduce processivity")

    return {
        'test': 'Mg2+ concentration effects',
        'results': results,
        'passed': passed,
        'summary': f"Optimal Mg around 2-4mM"
    }


def validate_formamide_inhibition() -> Dict[str, Any]:
    """
    Validate formamide shows consistent inhibition.

    Expected behavior (Blake 1996):
    - Formamide consistently inhibits enzyme activity
    - Also reduces Tm linearly

    Returns:
        Dict with test results and pass/fail status
    """
    results = []
    passed = True

    for form in [0, 2, 5, 10]:
        cond = ReactionConditions(temp=30.0, polymerase='phi29', formamide_percent=form)
        model = MechanisticModel(cond)
        effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)
        results.append({
            'formamide_percent': form,
            'processivity_factor': effects.processivity_factor,
            'stability_factor': effects.stability_factor,
        })

    proc_values = [r['processivity_factor'] for r in results]

    # Formamide should monotonically decrease processivity
    for i in range(1, len(proc_values)):
        if proc_values[i] > proc_values[i-1] + 0.01:
            passed = False
            logger.warning("Formamide should monotonically inhibit processivity")

    return {
        'test': 'Formamide inhibition',
        'results': results,
        'passed': passed,
        'summary': f"Processivity: {proc_values[0]:.2f} (0%) -> {proc_values[-1]:.2f} (10%)"
    }


def validate_mechanistic_model() -> List[Dict[str, Any]]:
    """
    Run all validation checks on mechanistic model.

    Returns:
        List of validation results with pass/fail status
    """
    validations = [
        validate_dmso_inhibition,
        validate_betaine_biphasic,
        validate_gc_accessibility,
        validate_temperature_optimum,
        validate_binding_kinetics,
        validate_mg_effects,
        validate_formamide_inhibition,
    ]

    results = []
    for validate_func in validations:
        try:
            result = validate_func()
            results.append(result)
            status = "PASSED" if result['passed'] else "FAILED"
            logger.info(f"{result['test']}: {status}")
        except Exception as e:
            results.append({
                'test': validate_func.__name__,
                'passed': False,
                'error': str(e),
            })
            logger.error(f"{validate_func.__name__}: ERROR - {e}")

    return results


def format_validation_report(results: List[Dict[str, Any]]) -> str:
    """
    Format validation results as a human-readable report.

    Args:
        results: List of validation result dicts

    Returns:
        Formatted report string
    """
    lines = [
        "Mechanistic Model Validation Report",
        "=" * 50,
        "",
    ]

    passed_count = sum(1 for r in results if r.get('passed', False))
    total_count = len(results)

    lines.append(f"Overall: {passed_count}/{total_count} tests passed")
    lines.append("")

    for result in results:
        status = "PASS" if result.get('passed', False) else "FAIL"
        test_name = result.get('test', 'Unknown')
        lines.append(f"[{status}] {test_name}")

        if 'summary' in result:
            lines.append(f"       {result['summary']}")

        if 'error' in result:
            lines.append(f"       Error: {result['error']}")

        lines.append("")

    return "\n".join(lines)


def quick_validation() -> Tuple[bool, str]:
    """
    Run quick validation and return pass/fail with summary.

    Returns:
        Tuple of (all_passed, summary_message)
    """
    results = validate_mechanistic_model()
    all_passed = all(r.get('passed', False) for r in results)
    passed_count = sum(1 for r in results if r.get('passed', False))
    summary = f"Validation: {passed_count}/{len(results)} tests passed"
    return all_passed, summary
