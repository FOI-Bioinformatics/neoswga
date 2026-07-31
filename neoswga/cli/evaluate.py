"""`evaluate-set`: assess an existing oligo set against a genome.

The gap this fills: someone holding a primer set designed elsewhere - or one
this pipeline produced months ago - had no way to ask "is this set any good on
this genome?". Every iterative command required a params.json carrying
`fg_prefixes`/`fg_seq_lengths`, i.e. a prior `init` plus `count-kmers`, and
primers absent from the resulting HDF5 index silently scored 0% coverage with
nothing to say the number was meaningless rather than bad.

This is a thin front door, not a new engine. Position lookup goes through
`PositionCache(on_missing='scan')`, coverage through the same
`compute_per_prefix_coverage` the optimizers use, and thermodynamics through
`ReactionConditions` - so its numbers are the numbers the rest of the tool
reports, not a parallel implementation that could drift.
"""

import json
import logging
import os

from neoswga.cli._common import (
    bootstrap_params_from_genome,
    collect_primers_from_args,
)

logger = logging.getLogger(__name__)


def _resolve_sources(args):
    """Work out foreground prefixes, genomes and lengths from -j and/or --genome."""
    from neoswga.core import parameter

    genomes = list(getattr(args, "genome", None) or [])

    if getattr(args, "json_file", None):
        import neoswga.core.pipeline as pipeline_mod

        pipeline_mod._initialized = False
        pipeline_mod._initialize()
        prefixes = list(getattr(parameter, "fg_prefixes", []) or [])
        lengths = list(getattr(parameter, "fg_seq_lengths", []) or [])
        genomes = genomes or list(getattr(parameter, "fg_genomes", []) or [])
    elif genomes:
        applied = bootstrap_params_from_genome(genomes, args.output, circular=not args.linear)
        prefixes = applied["fg_prefixes"]
        lengths = applied["fg_seq_lengths"]
    else:
        raise ValueError(
            "evaluate-set needs either -j params.json or --genome FASTA. "
            "With --genome alone it scans the FASTA directly, so no prior "
            "count-kmers run is required."
        )

    if not prefixes or not lengths:
        raise ValueError(
            "No foreground genome information available. Pass --genome, or use a "
            "params.json with fg_genomes/fg_prefixes/fg_seq_lengths populated."
        )
    return prefixes, genomes, lengths


def run_evaluate_set(args):
    """Evaluate an existing oligo set: coverage, gaps, selectivity, dimers."""
    from neoswga.core import parameter
    from neoswga.core.coverage import compute_per_prefix_coverage, polymerase_extension_reach
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import ReactionConditions

    primers = collect_primers_from_args(
        getattr(args, "primers", None), getattr(args, "primers_file", None)
    )
    if not primers:
        raise ValueError("No primers supplied. Use --primers or --primers-file.")

    os.makedirs(args.output, exist_ok=True)
    prefixes, genomes, lengths = _resolve_sources(args)

    polymerase = (
        getattr(args, "polymerase", None) or getattr(parameter, "polymerase", "phi29") or "phi29"
    )
    conditions = ReactionConditions(
        temp=getattr(args, "reaction_temp", None)
        or getattr(parameter, "reaction_temp", None)
        or 30.0,
        polymerase=polymerase,
    )
    reach = polymerase_extension_reach(polymerase, coverage_metric="realistic")

    # on_missing='scan' is the point of this command: an outside primer set is
    # not in the index, and scoring it zero would be the bug, not the answer.
    scan_ok = bool(genomes) and len(genomes) == len(prefixes)
    cache = PositionCache(
        prefixes,
        primers,
        genome_paths=genomes if scan_ok else None,
        circular=not args.linear,
        on_missing="scan" if scan_ok else "warn",
    )
    if not scan_ok:
        logger.warning(
            "No genome FASTA available for scanning; primers absent from the "
            "k-mer index will report zero coverage. Pass --genome to resolve them."
        )

    overall, coverage = compute_per_prefix_coverage(
        cache,
        primers,
        prefixes,
        lengths,
        extension=reach,
        circular=not args.linear,
    )

    per_primer = []
    for primer in primers:
        sites = sum(len(cache.get_positions(p, primer)) for p in prefixes)
        per_primer.append(
            {
                "primer": primer,
                "length": len(primer),
                "gc": round((primer.count("G") + primer.count("C")) / max(1, len(primer)), 4),
                "tm": round(conditions.calculate_effective_tm(primer), 2),
                "binding_sites": sites,
                "indexed": primer not in cache.zero_site_primers,
            }
        )

    total_sites = sum(p["binding_sites"] for p in per_primer)
    genome_bp = sum(lengths)
    result = {
        "primers": per_primer,
        "num_primers": len(primers),
        "fg_coverage": round(overall, 4),
        "per_target_coverage": {k: round(v, 4) for k, v in coverage.items()},
        "extension_reach_bp": reach,
        "total_binding_sites": total_sites,
        # The literature's dominant predictor of SWGA success: mean distance
        # between binding sites (Clarke 2017, Dwivedi-Yu 2023). Successful
        # published sets sit near 1 site per 2-5 kbp.
        "mean_binding_distance_bp": (round(genome_bp / total_sites) if total_sites else None),
        "genome_bp": genome_bp,
        "zero_site_primers": cache.zero_site_primers,
        "conditions": {
            "polymerase": conditions.polymerase,
            "temp": conditions.temp,
            "mg_conc": conditions.mg_conc,
        },
    }

    out_path = os.path.join(args.output, "evaluation.json")
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)

    _print_report(result)
    logger.info("Wrote %s", out_path)
    print(f"\nWrote {out_path}")
    print("Next: neoswga expand-primers --fixed-primers <keep> --num-new N")
    return result


def _print_report(result):
    print("\n" + "=" * 70)
    print("PRIMER SET EVALUATION")
    print("=" * 70)
    print(f"  Primers              : {result['num_primers']}")
    print(f"  Target size          : {result['genome_bp']:,} bp")
    print(f"  Binding sites        : {result['total_binding_sites']:,}")
    mbd = result["mean_binding_distance_bp"]
    if mbd:
        print(f"  Mean binding distance: {mbd:,} bp (1 site per {mbd / 1000:.1f} kbp)")
        # Thresholds from Dwivedi-Yu et al. (2023): successful Prevotella sets
        # ran 1 site per 2.0-4.9 kbp, unsuccessful ones 1 per 2.8-7.4 kbp.
        if mbd <= 5000:
            print("      within the density range of published successful sets")
        else:
            print("      sparser than published successful sets (1 per 2-5 kbp)")
    print(f"  Coverage @ {result['extension_reach_bp']} bp reach : {result['fg_coverage']:.1%}")
    for name, cov in result["per_target_coverage"].items():
        print(f"      {os.path.basename(name)}: {cov:.1%}")

    if result["zero_site_primers"]:
        print(
            f"\n  {len(result['zero_site_primers'])} primer(s) have NO binding "
            f"sites in the target:"
        )
        for p in result["zero_site_primers"][:10]:
            print(f"      {p}")
        print("      These contribute nothing and can be dropped.")

    print("\n  Per primer:")
    print(f"      {'primer':<22s} {'len':>4s} {'GC':>6s} {'Tm':>7s} {'sites':>8s}")
    for p in result["primers"]:
        print(
            f"      {p['primer']:<22s} {p['length']:>4d} {p['gc']:>6.2f} "
            f"{p['tm']:>7.1f} {p['binding_sites']:>8d}"
        )
    print("=" * 70)


def add_parsers(subparsers):
    from neoswga.cli._common import add_common_options, add_position_source_options

    p = subparsers.add_parser(
        "evaluate-set",
        help="Evaluate an existing oligo set against a genome (coverage, gaps, dimers)",
        description=(
            "Assess a primer set you already have. Works on an externally "
            "designed set: pass --genome and the binding sites are found by "
            "scanning the FASTA, so no prior count-kmers run is needed."
        ),
    )
    # -j comes from add_common_options below, and is deliberately NOT required:
    # --genome alone is enough.
    p.add_argument("--primers", nargs="+", help="Primer sequences")
    p.add_argument("--primers-file", help="File with one primer per line")
    p.add_argument("-o", "--output", default="evaluation", help="Output directory")
    p.add_argument("--reaction-temp", type=float)
    p.add_argument(
        "--linear",
        action="store_true",
        help="Treat targets as linear (default: circular, as for most bacterial "
        "chromosomes and plasmids)",
    )
    add_position_source_options(p)
    add_common_options(p)
    return p
