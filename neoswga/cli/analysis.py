"""Analysis CLI handlers: presets display, condition optimisation, and the
``analyze-*`` primer/genome/dimer/stability commands.

Extracted from cli_unified.py. Shared helpers come from neoswga.cli._common.
"""

import json
import logging
import os
import sys

from neoswga.cli._common import PRESETS, validate_primer_sequence

logger = logging.getLogger(__name__)


def show_presets():
    """Show available reaction condition presets"""
    print("\n" + "=" * 70)
    print("Available Reaction Condition Presets")
    print("=" * 70 + "\n")

    # Show presets from the PRESETS dict
    for name, config in PRESETS.items():
        print(f"{name}:")
        print(f"  Temperature: {config['temperature']}C")
        print(f"  Polymerase: {config['polymerase']}")
        print(f"  DMSO: {config['dmso_percent']}%")
        print(f"  Betaine: {config['betaine_m']} M")
        print(f"  SSB: {config['ssb']}")
        print(f"  Optimization: {config['optimization_method']}")
        print()

    # Show additional presets available via --preset but not in PRESETS dict
    additional = [
        "q_solution",
        "gc_melt",
        "crude_sample",
        "low_temp",
        "bst",
        "klenow",
        "extreme_gc",
    ]
    try:
        from neoswga.core import reaction_conditions as rc

        preset_funcs = {
            "q_solution": ("Q-solution equivalent", rc.get_q_solution_equivalent),
            "gc_melt": ("GC-melt conditions", rc.get_gc_melt_conditions),
            "crude_sample": ("Crude sample conditions", rc.get_crude_sample_conditions),
            "low_temp": ("Low temperature conditions", rc.get_low_temp_conditions),
            "bst": ("Bst polymerase conditions", rc.get_bst_conditions),
            "klenow": ("Klenow polymerase conditions", rc.get_klenow_conditions),
            "extreme_gc": ("Extreme GC content conditions", rc.get_extreme_gc_conditions),
        }
        for name in additional:
            desc, func = preset_funcs[name]
            try:
                cond = func()
                print(f"{name}:")
                print(f"  {desc}")
                print(
                    f"  Temperature: {cond.temp}C, Na+: {cond.na_conc} mM, Mg2+: {cond.mg_conc} mM"
                )
                if hasattr(cond, "dmso_percent") and cond.dmso_percent:
                    print(f"  DMSO: {cond.dmso_percent}%")
                if hasattr(cond, "betaine_m") and cond.betaine_m:
                    print(f"  Betaine: {cond.betaine_m} M")
                print()
            except Exception:
                print(f"{name}: {desc}")
                print()
    except ImportError:
        for name in additional:
            print(f"{name}: (use --preset {name} with filter command)")
        print()


def optimize_conditions(args):
    """Optimize reaction conditions for a genome"""
    from Bio import SeqIO

    from neoswga.core import reaction_conditions as rc

    logger.info("Analyzing genome and recommending optimal conditions...")

    # Load genome
    records = list(SeqIO.parse(args.fg, "fasta"))
    genome_seq = str(records[0].seq).upper()
    genome_length = len(genome_seq)

    # Calculate GC content
    gc_content = (genome_seq.count("G") + genome_seq.count("C")) / genome_length

    print(f"\nGenome properties:")
    print(f"  Length: {genome_length:,} bp")
    print(f"  GC content: {gc_content:.1%}")

    # Recommend conditions
    recommendations = rc.recommend_conditions(genome_seq, target_k=args.target_k)

    print(f"\nRecommended conditions:")
    print(f"  Optimal k-mer length: {recommendations['optimal_k']}")
    print(f"  Recommended temperature: {recommendations['temperature']:.1f}C")
    print(f"  Polymerase: {recommendations['polymerase']}")
    print(f"  DMSO: {recommendations['dmso_percent']:.1f}%")
    print(f"  Betaine: {recommendations['betaine_m']:.1f} M")
    if recommendations["ssb"]:
        print(f"  SSB: Recommended")

    # Save recommendations
    os.makedirs(args.output, exist_ok=True)
    rec_file = os.path.join(args.output, "recommended_conditions.json")
    with open(rec_file, "w") as f:
        json.dump(recommendations, f, indent=2)

    print(f"\nRecommendations saved to: {rec_file}")


def analyze_primer_set(args):
    """Analyze an existing primer set"""
    from neoswga.core import reaction_conditions as rc
    from neoswga.core import secondary_structure as ss
    from neoswga.core import thermodynamics as thermo

    logger.info("Analyzing primer set...")

    # Load preset config
    config_dict = PRESETS[args.preset]
    conditions = rc.ReactionConditions(
        temp=config_dict["temperature"],
        dmso_percent=config_dict["dmso_percent"],
        betaine_m=config_dict["betaine_m"],
        polymerase=config_dict["polymerase"],
    )

    primers = args.primers

    print(f"\nAnalyzing {len(primers)} primers:")
    for primer in primers:
        print(f"  {primer}")
    print()

    # Thermodynamic analysis
    print("Thermodynamic Properties:")
    for primer in primers:
        tm = conditions.calculate_effective_tm(primer)
        dg = thermo.calculate_free_energy(primer, conditions.temp)
        gc = thermo.gc_content(primer)
        print(f"  {primer}:")
        print(f"    Tm: {tm:.1f}C")
        print(f"    dG: {dg:.2f} kcal/mol")
        print(f"    GC: {gc:.1%}")

    # Secondary structure analysis
    print("\nSecondary Structure Analysis:")
    for primer in primers:
        hairpins = ss.check_hairpins(primer, conditions)
        homodimer = ss.check_homodimer(primer, conditions)
        print(f"  {primer}:")
        if hairpins:
            worst = min(hairpins, key=lambda h: h.energy)
            print(f"    Hairpin: dG={worst.energy:.2f} kcal/mol, stem={worst.stem_length}bp")
        else:
            print(f"    Hairpin: none detected")
        print(f"    Homodimer: dG={homodimer.energy:.2f} kcal/mol")

    # Check heterodimers
    print("\nHeterodimer Analysis:")
    for i, p1 in enumerate(primers):
        for j, p2 in enumerate(primers):
            if i < j:
                result = ss.check_heterodimer(p1, p2, conditions)
                if result.energy < -6.0:
                    print(f"  {p1} x {p2}: dG={result.energy:.2f} kcal/mol (warning)")

    print("\nAnalysis complete!")


# =========================================================================
# CATEGORY 1: Handler functions for orphaned features
# =========================================================================


def run_analyze_genome(args):
    """Analyze genome suitability for SWGA"""
    from neoswga.core import genome_analysis

    logger.info(f"Analyzing genome: {args.genome}")

    try:
        results = genome_analysis.analyze_genome(
            genome_path=args.genome, window_size=args.window_size, output_dir=args.output
        )

        logger.info(f"Analysis complete! Report saved to: {args.output}")

    except Exception as e:
        logger.error(f"Genome analysis failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def run_analyze_dimers(args):
    """Analyze primer dimer interaction network"""
    from neoswga.core.dimer_network_analyzer import DimerNetworkAnalyzer

    # Validate primer sequences
    for primer in args.primers:
        validate_primer_sequence(primer, name="--primers")

    logger.info(f"Analyzing dimer network for {len(args.primers)} primers")

    try:
        analyzer = DimerNetworkAnalyzer(severity_threshold=args.threshold)
        metrics, profiles, matrix = analyzer.analyze_primer_set(args.primers)

        os.makedirs(args.output, exist_ok=True)

        # Write summary
        summary_path = os.path.join(args.output, "dimer_summary.txt")
        with open(summary_path, "w") as f:
            f.write(str(metrics) + "\n")
            for primer, profile in profiles.items():
                f.write(f"\n{primer}: {profile}\n")

        logger.info(f"Dimer analysis complete! Results saved to: {args.output}")

        if getattr(metrics, "num_hub_primers", 0) > 0:
            logger.warning(f"Found {metrics.num_hub_primers} hub primers with many interactions")

    except Exception as e:
        logger.error(f"Dimer analysis failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


def run_analyze_stability(args):
    """Analyze 3' end stability and specificity"""
    from neoswga.core.three_prime_stability import ThreePrimeStabilityAnalyzer

    # Validate primer sequences
    for primer in args.primers:
        validate_primer_sequence(primer, name="--primers")

    logger.info(f"Analyzing 3' stability for {len(args.primers)} primers")

    try:
        analyzer = ThreePrimeStabilityAnalyzer()
        results = []
        for primer in args.primers:
            stability = analyzer.analyze_primer(primer)
            results.append(stability)

        # Ensure output directory exists
        output_dir = os.path.dirname(args.output)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        # Write report
        with open(args.output, "w") as f:
            for primer, stability in zip(args.primers, results):
                f.write(f"{primer}\t{stability}\n")

        logger.info(f"Stability analysis complete! Report saved to: {args.output}")

    except Exception as e:
        logger.error(f"Stability analysis failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
