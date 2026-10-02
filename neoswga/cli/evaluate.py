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

The per-reference blocks added on 2026-10-02 keep that contract: they are
composed by `core/reference_panel_evaluation.py`, which this handler calls once.
What they add is the question no command could ask before - how does this set do
on EACH target and EACH host, including a host that was never in the design's
params.json - and every pooled field keeps the name and the arithmetic it had.
"""

import json
import logging
import os

from neoswga.cli._common import (
    background_references_from_genomes,
    bootstrap_params_from_genome,
    collect_primers_from_args,
    merge_args_to_parameter,
    params_command,
)

logger = logging.getLogger(__name__)


def _resolve_sources(args):
    """Work out foreground prefixes, genomes and lengths from -j and/or --genome."""
    from neoswga.core import parameter

    genomes = list(getattr(args, "genome", None) or [])

    if getattr(args, "json_file", None):
        import neoswga.core.pipeline as pipeline_mod

        # `_initialize()` reads `parameter.json_file`, not an argument, so the
        # path has to be on the module global before the call. Setting it here
        # as well as in the @params_command decorator is deliberate: this
        # function is only correct if the assignment happened, and depending on
        # a caller to have done it is how it came to be missing. The merge
        # skips None, so it is idempotent.
        merge_args_to_parameter(args, parameter, ["json_file"])

        pipeline_mod._initialized = False
        try:
            pipeline_mod._initialize()
        except KeyError as e:
            # `_initialize` is shared with the pipeline steps, where a missing
            # foreground key is fatal and a bare KeyError is fine. Here it is
            # not: this command can work from a FASTA alone, so it can offer a
            # way forward that the generic handler cannot.
            raise ValueError(
                f"params.json is missing {e} needed to locate the foreground "
                f"genome. Pass --genome FASTA to scan it directly instead, or "
                f"run count-kmers to populate the k-mer index."
            ) from e
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


def _resolve_primers(args):
    """The oligo set: either named on the command line, or a delivered set.

    `--from-results DIR` reads ONE set through `core/delivered_set.py`, which is
    the same reader `export`, `interpret`, `report` and `simulate` use. The
    alternatives in `step4_improved_df.csv` are separate answers to one design,
    so evaluating their union would describe a tube nobody would set up and a
    dimer screen has never compared.
    """
    from_results = getattr(args, "from_results", None)
    explicit = getattr(args, "primers", None) or getattr(args, "primers_file", None)
    if from_results and explicit:
        raise ValueError(
            "Pass either --from-results or --primers/--primers-file, not both: "
            "a delivered set and a hand-written list are two different panels."
        )
    if not from_results:
        return collect_primers_from_args(
            getattr(args, "primers", None), getattr(args, "primers_file", None)
        )

    from neoswga.core.delivered_set import read_delivered_set

    path = from_results
    if os.path.isdir(path):
        path = os.path.join(path, "step4_improved_df.csv")
    delivered = read_delivered_set(path, getattr(args, "set_index", 0))
    primers = list(delivered.primers)
    if not primers:
        raise ValueError(f"No primers in {path}")
    logger.info("Evaluating %s", delivered.describe())
    return primers


def _reference_specs(args, prefixes, genomes, lengths, bg_prefixes, bg_genomes, bg_lengths):
    """One `ReferenceSpec` per genome, targets and hosts together.

    Prefix, genome and length travel together on each spec. The hosts given by
    `--background` are scanned only when `--scan-background` says so; otherwise
    their sites come from k-mer counts and their coverage is reported as
    unavailable with that as the reason, which is the honest answer for a
    host-sized genome rather than a zero.
    """
    from neoswga.core import parameter
    from neoswga.core.reference_panel_evaluation import ROLE_HOST, ROLE_TARGET, ReferenceSpec

    circular_fg = not args.linear
    scan_bg = bool(getattr(args, "scan_background", False))
    # The same key the existing background block reads, so a host configured in
    # params.json is treated identically in both places.
    bg_circular = bool(getattr(parameter, "bg_circular", False))

    specs = []
    aligned_genomes = genomes if len(genomes) == len(prefixes) else [None] * len(prefixes)
    for prefix, genome, length in zip(prefixes, aligned_genomes, lengths, strict=True):
        specs.append(
            ReferenceSpec(
                prefix=prefix,
                genome=genome,
                length=int(length),
                role=ROLE_TARGET,
                circular=circular_fg,
                scan=True,
            )
        )

    host_genomes = bg_genomes if len(bg_genomes) == len(bg_prefixes) else [None] * len(bg_prefixes)
    for prefix, genome, length in zip(bg_prefixes, host_genomes, bg_lengths, strict=True):
        specs.append(
            ReferenceSpec(
                prefix=prefix,
                genome=genome,
                length=int(length),
                role=ROLE_HOST,
                circular=bg_circular,
                scan=scan_bg,
            )
        )
    return specs


def _per_reference_blocks(args, primers, conditions, reach, foreground, background, totals):
    """The per-reference blocks, and the pooled host fields they may fill.

    Per reference, separately: the question the pooled figures cannot answer. A
    primer that binds only one of two targets is invisible in `fg_coverage`, and
    a host that is a small share of the pooled bases is invisible in
    `selectivity_ratio` -- Known Issue 6 in a second place.

    Hosts given by `--background` are kept apart from the configured ones all the
    way through. `totals` comes back unchanged unless the ONLY hosts are ones
    given by path: with a host configured in params.json the pooled fields keep
    the arithmetic they have always had, so no existing output moves, and with
    none they were simply None and a measurement now exists.
    """
    from neoswga.core.reference_panel_evaluation import evaluate_reference_panel

    prefixes, genomes, lengths = foreground
    bg_prefixes, bg_genomes, bg_lengths = background
    total_sites, total_bg_sites, selectivity = totals

    extra_hosts = []
    if getattr(args, "background", None):
        extra_hosts = background_references_from_genomes(args.background, args.output)

    panel = evaluate_reference_panel(
        primers,
        _reference_specs(
            args,
            prefixes,
            genomes,
            lengths,
            list(bg_prefixes) + [host["prefix"] for host in extra_hosts],
            list(bg_genomes) + [host["genome"] for host in extra_hosts],
            list(bg_lengths) + [host["length"] for host in extra_hosts],
        ),
        conditions,
        extension=reach,
    )
    if not bg_prefixes and extra_hosts:
        total_bg_sites, selectivity = _pooled_host_sites(panel, total_sites)
    return panel.as_dict(), total_bg_sites, selectivity


def _pooled_host_sites(panel, total_sites):
    """The pooled host site count and ratio, for hosts given by `--background`.

    None when ANY of those hosts could not be measured. A sum over the hosts
    that happened to answer is not the panel's host binding, and it reads as
    though it were -- the same shape as the reductions inside
    `reference_panel_evaluation`. The per-host rows sit beside it in the JSON,
    so a reader never has only the pooled figure for hosts of different sizes.
    """
    hosts = panel.hosts
    if not hosts or any(host.sites.value is None for host in hosts):
        return None, None
    pooled = int(sum(host.sites.value for host in hosts))
    ratio = round(total_sites / pooled, 4) if pooled else None
    return pooled, ratio


def _gap_statistics(cache, primers, prefixes, lengths, circular):
    """Mean, max and Gini of the distances between binding sites.

    The literature disagrees about which of these predicts SWGA success --
    see docs/validation/published_primer_sets.md -- so all three are
    reported rather than one being chosen.
    """
    import numpy as _np

    gaps = []
    for prefix, length in zip(prefixes, lengths, strict=True):
        pos = sorted({int(x) for primer in primers for x in cache.get_positions(prefix, primer)})
        if len(pos) < 2:
            continue
        d = _np.diff(pos).tolist()
        if circular:
            d.append(length - pos[-1] + pos[0])  # wrap-around gap
        gaps.extend(d)

    if not gaps:
        return 0.0, 0.0, 0.0

    arr = _np.array(sorted(gaps), dtype=float)
    n = len(arr)
    gini = (
        float(2.0 * _np.sum((_np.arange(1, n + 1)) * arr) / (n * arr.sum()) - (n + 1) / n)
        if arr.sum() > 0
        else 0.0
    )
    return float(arr.mean()), float(arr.max()), gini


def _occupancy_selectivity(primers, fg_prefixes, bg_prefixes, conditions):
    """Occupancy-weighted selectivity, matching what the optimizers report.

    Kept beside the exact-count ratio rather than replacing it: this command
    exists to evaluate an outside oligo set, and a reader comparing against a
    published fg/bg figure needs the count as well as the weighted value. On a
    set with no exact background matches the two disagree completely -- the
    count says "no background binding" while this finds the near-match load
    stringency can actually act on.

    Returns None when the jellyfish count files are absent, which is the case
    with `--genome` alone: positions are scanned straight from a FASTA and
    there are no counts. Absent rather than silently substituted by the count
    ratio, since they are different numbers.
    """
    from neoswga.core.occupancy import weighted_site_load

    if not bg_prefixes:
        return None

    try:
        fg_load = weighted_site_load(primers, fg_prefixes, conditions, 1)
        bg_load = weighted_site_load(primers, bg_prefixes, conditions, 1)
    except (FileNotFoundError, OSError):
        logger.info("No k-mer count files available; reporting exact-match selectivity only.")
        return None

    return round(fg_load / bg_load, 4) if bg_load > 0 else None


# Every other params-taking command carries this decorator; evaluate-set was
# the one that did not. It validates the path and, critically, merges
# `args.json_file` onto `parameter.json_file` -- which is what
# `pipeline._initialize()` reads. Without it `-j` was accepted and ignored:
# the file was never opened, and the command then reported the empty globals
# it found as though they were the user's configuration ("Missing required
# parameter: 'fg_prefixes'" against a file that lists one). The decorator
# skips None, so the --genome-only path is unaffected.
@params_command
def run_evaluate_set(args):
    """Evaluate an existing oligo set: coverage, gaps, selectivity, dimers."""
    from neoswga.core import parameter
    from neoswga.core.coverage import (
        compute_per_prefix_coverage,
        gap_regime_note,
        interpret_gap_metrics,
        polymerase_extension_reach,
    )
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import build_reaction_conditions

    primers = _resolve_primers(args)
    os.makedirs(args.output, exist_ok=True)
    prefixes, genomes, lengths = _resolve_sources(args)
    bg_prefixes = list(getattr(parameter, "bg_prefixes", []) or [])
    bg_genomes = list(getattr(parameter, "bg_genomes", []) or [])
    bg_lengths = list(getattr(parameter, "bg_seq_lengths", []) or [])

    polymerase = (
        getattr(args, "polymerase", None) or getattr(parameter, "polymerase", "phi29") or "phi29"
    )
    # This passed two of twenty parameters, so every additive and buffer
    # species in params.json was absent from the Tm and dimer energies this
    # command reports -- the point of the command being to evaluate an oligo
    # set under the user's chemistry.
    conditions = build_reaction_conditions(args, polymerase=polymerase)
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

    # Three outcomes, and conflating them is precisely the bug this command
    # exists to prevent:
    #   ok           - looked up, binds somewhere
    #   no_sites     - looked up, binds nowhere. A real result.
    #   not_indexed  - never resolved. Its zero is meaningless, not low.
    # `missing_primers` is populated only when nothing could be scanned, and is
    # cleared once a scan resolves everything, so the two lists are disjoint.
    unresolved = sorted({p for _, p in cache.missing_primers})
    zero_site = set(cache.zero_site_primers)

    per_primer = []
    for primer in primers:
        sites = sum(len(cache.get_positions(p, primer)) for p in prefixes)
        if primer in unresolved:
            status = "not_indexed"
        elif sites == 0 or primer in zero_site:
            status = "no_sites"
        else:
            status = "ok"
        per_primer.append(
            {
                "primer": primer,
                "length": len(primer),
                "gc": round((primer.count("G") + primer.count("C")) / max(1, len(primer)), 4),
                "tm": round(conditions.calculate_effective_tm(primer), 2),
                "binding_sites": sites,
                "status": status,
            }
        )

    # Background sites, and the selectivity they give.
    #
    # This command is documented as reporting "coverage, gaps, selectivity,
    # dimers" and had no background handling at all: params.json carries
    # bg_prefixes, every other command counts them, and this one ignored them.
    # A user with a host genome configured got a report with no specificity in
    # it -- for SWGA, the half that decides whether a design is usable.
    # `--scan-background` was a flag for the missing capability.
    #
    # Scanning is opt-in for the reason its help gives: a host-sized genome
    # held in memory is expensive, and background specificity is a frequency
    # question the k-mer counts already answer.
    bg_per_primer = {}
    occupancy_selectivity = None
    if bg_prefixes:
        scan_bg = bool(getattr(args, "scan_background", False)) and len(bg_genomes) == len(
            bg_prefixes
        )
        bg_cache = PositionCache(
            bg_prefixes,
            primers,
            genome_paths=bg_genomes if scan_bg else None,
            circular=bool(getattr(parameter, "bg_circular", False)),
            on_missing="scan" if scan_bg else "warn",
        )
        for primer in primers:
            bg_per_primer[primer] = sum(
                len(bg_cache.get_positions(prefix, primer)) for prefix in bg_prefixes
            )
        for record in per_primer:
            record["background_sites"] = bg_per_primer.get(record["primer"], 0)

    total_sites = sum(p["binding_sites"] for p in per_primer)
    total_bg_sites = sum(bg_per_primer.values()) if bg_prefixes else None
    # None rather than 0 where no background is configured: an absent
    # background is "not measured", and reporting 0 would read as perfect
    # specificity -- a confident claim from no evidence.
    selectivity = None
    if bg_prefixes:
        selectivity = round(total_sites / total_bg_sites, 4) if total_bg_sites else None

        occupancy_selectivity = _occupancy_selectivity(primers, prefixes, bg_prefixes, conditions)
    genome_bp = sum(lengths)

    # Gap statistics from the union of binding positions. Reported alongside
    # mean spacing because the two published benchmarks disagree about which of
    # the three predicts success - see docs/validation/.
    mean_gap, max_gap, gap_gini = _gap_statistics(
        cache, primers, prefixes, lengths, circular=not args.linear
    )

    panel_block, total_bg_sites, selectivity = _per_reference_blocks(
        args,
        primers,
        conditions,
        reach,
        (prefixes, genomes, lengths),
        (bg_prefixes, bg_genomes, bg_lengths),
        (total_sites, total_bg_sites, selectivity),
    )

    result = {
        "primers": per_primer,
        "num_primers": len(primers),
        "fg_coverage": round(overall, 4),
        "per_target_coverage": {k: round(v, 4) for k, v in coverage.items()},
        "extension_reach_bp": reach,
        "total_binding_sites": total_sites,
        "total_background_sites": total_bg_sites,
        "selectivity_ratio": selectivity,
        "occupancy_selectivity_ratio": occupancy_selectivity,
        # The literature's dominant predictor of SWGA success: mean distance
        # between binding sites (Clarke 2017, Dwivedi-Yu 2023). Successful
        # published sets sit near 1 site per 2-5 kbp.
        "mean_binding_distance_bp": (round(genome_bp / total_sites) if total_sites else None),
        "genome_bp": genome_bp,
        "mean_gap_bp": round(mean_gap),
        "max_gap_bp": round(max_gap),
        "gap_gini": round(gap_gini, 4),
        "gap_interpretation": interpret_gap_metrics(
            mean_gap, max_gap, gap_gini, extension_reach=reach
        ),
        "zero_site_primers": cache.zero_site_primers,
        # Primers whose coverage could not be determined at all. Reported
        # separately because their zero is an absence of evidence, not
        # evidence of absence, and a caller reading only the JSON would
        # otherwise see a confident-looking 0% with nothing to flag it.
        "unresolved_primers": unresolved,
        "_regime_note": gap_regime_note(gap_gini),
        "conditions": {
            "polymerase": conditions.polymerase,
            "temp": conditions.temp,
            "mg_conc": conditions.mg_conc,
        },
        # Per reference. `per_target_coverage` above keeps its shape and its
        # meaning; these blocks carry the site counts, densities, gaps and
        # weighted loads per genome, each measured against its OWN length, plus
        # the target-against-host table and the two reductions over it.
        "per_target": panel_block["per_target"],
        "per_host": panel_block["per_host"],
        "target_host_pairs": panel_block["target_host_pairs"],
        "worst_target_coverage": panel_block["worst_target_coverage"],
        "worst_host_selectivity_density": panel_block["worst_host_selectivity_density"],
        "per_reference_notes": panel_block["notes"],
    }

    out_path = os.path.join(args.output, "evaluation.json")
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)

    _print_report(result)
    logger.info("Wrote %s", out_path)
    print(f"\nWrote {out_path}")
    print("Next: neoswga expand-primers --fixed-primers <keep> --num-new N")
    return result


def _measured(entry):
    """A measurement's value, or None. One reader for the JSON shape."""
    return entry.get("value") if isinstance(entry, dict) else None


def _format_measurement(entry, fmt="{:.4g}"):
    value = _measured(entry)
    if value is None:
        reason = (entry or {}).get("unavailable") or "not measured"
        return f"unavailable ({reason})"
    return fmt.format(value)


def _print_per_reference(result):
    """The per-reference rows, and the two reductions over them.

    Printed with the reductions LAST, because they are the figures a reader acts
    on and an unavailable one has to be read as unknown rather than skipped.
    """
    targets = result.get("per_target") or {}
    hosts = result.get("per_host") or {}
    if not targets and not hosts:
        return

    print("\n  Per reference (each against its own length):")
    for name, record in list(targets.items()) + list(hosts.items()):
        label = os.path.basename(name)
        sites = _format_measurement(record.get("sites"), "{:.0f}")
        density = _format_measurement(record.get("sites_per_mb"), "{:.1f}")
        coverage = record.get("coverage") or {}
        covered = _measured(coverage)
        coverage_text = f"{covered:.1%}" if covered is not None else "coverage unavailable"
        print(
            f"      {record.get('role', '?'):<6s} {label:<22s} sites {sites:>8s}  "
            f"{density:>8s}/Mb  {coverage_text}"
        )
        if covered is None and coverage.get("unavailable"):
            print(f"          why: {coverage['unavailable']}")

    pairs = result.get("target_host_pairs") or []
    if pairs:
        print("\n  Target against host (density is the comparable figure;")
        print("  the ratio moves with host size -- Known Issue 6):")
        for pair in pairs:
            print(
                f"      {os.path.basename(pair['target']):<18s} vs "
                f"{os.path.basename(pair['host']):<18s} "
                f"density {_format_measurement(pair.get('selectivity_density'))}  "
                f"ratio {_format_measurement(pair.get('selectivity_ratio'))}"
            )

    worst_coverage = result.get("worst_target_coverage")
    worst_density = result.get("worst_host_selectivity_density")
    if worst_coverage or worst_density:
        print("\n  Worst case over the panel:")
        print(f"      worst target coverage      : {_format_measurement(worst_coverage, '{:.1%}')}")
        print(f"      worst host density         : {_format_measurement(worst_density)}")


def _print_report(result):
    print("\n" + "=" * 70)
    print("PRIMER SET EVALUATION")
    print("=" * 70)
    print(f"  Primers              : {result['num_primers']}")
    print(f"  Target size          : {result['genome_bp']:,} bp")
    print(f"  Binding sites        : {result['total_binding_sites']:,}")
    if result.get("total_background_sites") is not None:
        print(f"  Background sites     : {result['total_background_sites']:,}")
        ratio = result.get("selectivity_ratio")
        # "no background sites at all" is the best possible outcome, not a
        # missing measurement, so it gets said rather than left blank.
        print(
            f"  Selectivity (exact)  : "
            + (f"{ratio:.2f}x" if ratio is not None else "no exact background matches")
        )
    occupancy_ratio = result.get("occupancy_selectivity_ratio")
    if occupancy_ratio is not None:
        # Shown separately because the two answer different questions, and on
        # a set with no exact background matches they disagree completely: the
        # count says "no background binding" while the weighted figure finds
        # real near-match load that stringency can act on.
        print(
            f"  Selectivity (weighted): {occupancy_ratio:.2f}x  (near-matches, at reaction conditions)"
        )
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

    _print_per_reference(result)

    if result.get("gap_interpretation"):
        print("\n  Gap analysis (all three shown: the published benchmarks")
        print("  disagree about which one predicts success):")
        for r in result["gap_interpretation"]:
            print(f"      {r['label']:<22s} {r['value']:>9s}  {r['verdict']}")
        print(f"      -> {result['_regime_note']}")

    if result.get("unresolved_primers"):
        print(
            f"\n  {len(result['unresolved_primers'])} primer(s) could NOT be "
            f"looked up; their coverage is unknown, not zero:"
        )
        for p in result["unresolved_primers"][:10]:
            print(f"      {p}")
        print("      Pass --genome so they can be scanned directly.")

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
    p.add_argument(
        "--from-results",
        metavar="DIR",
        help="Read the oligo set from a finished run's step4_improved_df.csv "
        "instead of --primers. One set, not the union of the alternatives.",
    )
    p.add_argument(
        "--set",
        dest="set_index",
        type=int,
        default=0,
        help="Which of the alternative primer sets --from-results reads "
        "(default: 0, the set the run summary describes). Alternatives are "
        "separate answers to the same design, not additions to it.",
    )
    p.add_argument(
        "--background",
        nargs="+",
        metavar="FASTA",
        help="Host genome FASTA(s) to score against, needing no params.json "
        "entry. Sites come from k-mer counts, which is a frequency question a "
        "table or a single pass over the reference answers; pass "
        "--scan-background as well to locate them and measure host coverage, "
        "which holds the host in memory.",
    )
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
