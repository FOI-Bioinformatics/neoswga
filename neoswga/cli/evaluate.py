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


def _resolve_reach(args, polymerase):
    """The per-primer reach coverage is measured at, and where it came from.

    `coverage_reach` in params.json is honoured, through
    `coverage.resolve_coverage_reach`, the reader `optimize`, `plan-pool`,
    `expand-primers` and `improve-set` resolve it with. Until 2026-10-02 this
    command read only the polymerase, so a design made at a configured reach
    was evaluated at another one and the key was inert here (Known Issue 8):
    with `coverage_reach: 800` it reported 3000.

    The key is read only when a params file was given. With `--genome` alone
    the `parameter` module holds whatever an earlier call in this process left
    there, which is not this run's configuration.

    No params file this command accepted is refused by this: a value the reader
    would reject (zero, negative, not an integer) already fails the schema
    check in `validate_params_json_file`, before this function runs.
    """
    from neoswga.core import parameter
    from neoswga.core.coverage import resolve_coverage_reach

    configured = None
    if getattr(args, "json_file", None):
        configured = getattr(parameter, "coverage_reach", None)
    reach = resolve_coverage_reach(polymerase, override=configured)
    source = "coverage_reach in params.json" if configured is not None else "polymerase default"
    return reach, source


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


def _configured_background():
    """The hosts params.json configures, as (prefixes, genomes, lengths).

    Read in one place and carried as one bundle, so prefix, genome and length
    travel together: `tests/test_pairing_a_prefix_with_its_length_cannot_
    truncate.py` and `tests/test_expansion_counts_background.py` both exist
    because these three lists came apart somewhere.
    """
    from neoswga.core import parameter

    return (
        list(getattr(parameter, "bg_prefixes", []) or []),
        list(getattr(parameter, "bg_genomes", []) or []),
        list(getattr(parameter, "bg_seq_lengths", []) or []),
    )


def _background_totals(args, primers, per_primer, total_sites, prefixes, background, conditions):
    """Host sites per primer, the pooled total, and the two selectivity figures.

    This command is documented as reporting "coverage, gaps, selectivity,
    dimers" and had no background handling at all: params.json carries
    bg_prefixes, every other command counts them, and this one ignored them. A
    user with a host genome configured got a report with no specificity in it
    -- for SWGA, the half that decides whether a design is usable.

    Scanning is opt-in for the reason its help gives: a host-sized genome held
    in memory is expensive, and background specificity is a frequency question
    the k-mer counts already answer.

    `total_bg_sites` and the ratio are None where no background is configured.
    An absent background is "not measured", and 0 would read as perfect
    specificity -- a confident claim from no evidence.
    """
    from neoswga.core import parameter
    from neoswga.core.position_cache import PositionCache

    bg_prefixes, bg_genomes, _bg_lengths = background
    if not bg_prefixes:
        return None, None, None

    scan_bg = bool(getattr(args, "scan_background", False)) and len(bg_genomes) == len(bg_prefixes)
    bg_cache = PositionCache(
        bg_prefixes,
        primers,
        genome_paths=bg_genomes if scan_bg else None,
        circular=bool(getattr(parameter, "bg_circular", False)),
        on_missing="scan" if scan_bg else "warn",
    )
    bg_per_primer = {
        primer: sum(len(bg_cache.get_positions(prefix, primer)) for prefix in bg_prefixes)
        for primer in primers
    }
    for record in per_primer:
        record["background_sites"] = bg_per_primer.get(record["primer"], 0)

    total_bg_sites = sum(bg_per_primer.values())
    selectivity = round(total_sites / total_bg_sites, 4) if total_bg_sites else None
    occupancy = _occupancy_selectivity(primers, prefixes, bg_prefixes, conditions)
    return total_bg_sites, selectivity, occupancy


def _same_file(left, right):
    try:
        return os.path.samefile(left, right)
    except OSError:
        return os.path.realpath(left) == os.path.realpath(right)


def _variant_reference(prefixes, genomes, lengths, named=None):
    """The one target a variant table is placed against, or a refusal.

    A variant table states positions against ONE assembly. When the run has one
    foreground genome that is it. When it has more than one -- readable or not
    -- this refuses unless `--variants-reference FASTA` names one of them,
    because choosing the only readable one, or the first, would place the
    table in a coordinate space nobody chose, and the REF check cannot always
    tell near-identical strains apart. Each refusal names an action that
    works with the inputs the user already has, so following it cannot lead
    into another refusal.
    """
    from neoswga.core.exceptions import ReferenceDataError

    if len(genomes) != len(prefixes):
        raise ReferenceDataError(
            "variant route",
            f"the run has {len(prefixes)} foreground prefix(es) and "
            f"{len(genomes)} foreground FASTA path(s), so no prefix can be "
            f"paired with the FASTA its coordinates come from",
            "with -j, list one FASTA per fg_prefixes entry under fg_genomes in "
            "params.json and leave --genome off (it replaces that list); "
            "without -j, pass --genome FASTA",
        )
    if named:
        chosen = [i for i, genome in enumerate(genomes) if genome and _same_file(genome, named)]
        if not chosen:
            raise ReferenceDataError(
                "variant route",
                f"--variants-reference {named} is not one of this run's foreground "
                f"genomes ({', '.join(str(g) for g in genomes)})",
                "name one of the foreground genomes listed, exactly as the run configures it",
            )
        index = chosen[0]
    elif len(prefixes) > 1:
        raise ReferenceDataError(
            "variant route",
            f"a variant table is stated against one assembly, and this run "
            f"configures {len(prefixes)} foreground genomes, so which one the "
            f"table was made against cannot be inferred",
            "keep the run as it is and add --variants-reference FASTA naming "
            "one of: " + ", ".join(str(g) for g in genomes),
        )
    else:
        index = 0
    genome = genomes[index]
    if not genome or not os.path.isfile(genome):
        raise ReferenceDataError(
            "variant route",
            f"the foreground genome {genome!r} cannot be read, and it is the "
            f"reference the variant table would be placed against",
            "check the path in fg_genomes or --genome",
        )
    return prefixes[index], genome, lengths[index]


def _variant_blocks(args, primers, per_primer, cache, foreground, reach):
    """The `per_strain` block, or None when `--variants` was not given.

    Every figure in it comes from the functions the per-reference blocks use,
    over the intact sites only, so a strain carrying no variants reproduces the
    reference figures exactly. Nothing here is a score: an affected site is
    counted as lost, never weighted and never called tolerated, because the
    shipped mismatch model is uniform and its evidence status is `assumed`.
    """
    path = getattr(args, "variants", None)
    named = getattr(args, "variants_reference", None)
    if not path:
        if named:
            raise ValueError(
                "--variants-reference names the reference a variant table is "
                "placed against, and no --variants table was given"
            )
        return None

    from neoswga.core.exceptions import ReferenceDataError
    from neoswga.core.reference_layout import LayoutMismatch, verify_layout
    from neoswga.core.variant_sites import evaluate_strain_panel
    from neoswga.core.variant_table import open_variants

    prefixes, genomes, lengths = foreground
    prefix, genome, length = _variant_reference(prefixes, genomes, lengths, named)
    table = open_variants(path, genome, prefix=prefix)
    try:
        # The record starts are the ones the position index stored when
        # `filter` built it; a FASTA edited since then puts every variant at an
        # offset the positions no longer use. When this run scanned the FASTA
        # itself (no index for this prefix) the list is empty and
        # `verify_layout` records the starts check as not run, which is right:
        # positions scanned now come from this very file. The length half still
        # applies, since the configured length is every coverage denominator.
        verify_layout(table.layout, cache.get_record_starts(prefix), int(length))
    except LayoutMismatch as exc:
        raise ReferenceDataError(
            f"reference {genome}",
            str(exc),
            "re-run `neoswga filter` so the index describes this FASTA, or "
            "place the table against the FASTA the index was built from",
        ) from exc
    logger.info("Variant table %s placed against %s (prefix %s)", path, genome, prefix)
    block = evaluate_strain_panel(
        cache,
        prefix,
        primers,
        int(length),
        table,
        extension=int(reach),
        circular=not args.linear,
    )
    _attach_per_primer_intact(per_primer, block)
    return block


def _attach_per_primer_intact(per_primer, block):
    """Per primer, the fraction of its sites intact in every strain.

    None where any strain's answer is unknown, and None for a primer with no
    site on the reference. Zero of zero sites intact is not 0.0, and a 0.0 here
    would read as a primer whose sites the variants destroyed.
    """
    lowest = block.get("per_primer_intact_fraction_every_strain") or {}
    by_strain = {
        name: (record.get("per_primer_intact_fraction") or {})
        for name, record in (block.get("per_strain") or {}).items()
    }
    for record in per_primer:
        primer = record["primer"]
        value = lowest.get(primer)
        record["intact_fraction_every_strain"] = (
            round(float(value), 4) if value is not None else None
        )
        record["intact_fraction_by_strain"] = {
            name: (
                round(float(fractions[primer]), 4) if fractions.get(primer) is not None else None
            )
            for name, fractions in by_strain.items()
        }


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
    )
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import build_reaction_conditions

    primers = _resolve_primers(args)
    os.makedirs(args.output, exist_ok=True)
    prefixes, genomes, lengths = _resolve_sources(args)
    background = _configured_background()

    polymerase = (
        getattr(args, "polymerase", None) or getattr(parameter, "polymerase", "phi29") or "phi29"
    )
    # This passed two of twenty parameters, so every additive and buffer
    # species in params.json was absent from the Tm and dimer energies this
    # command reports -- the point of the command being to evaluate an oligo
    # set under the user's chemistry.
    conditions = build_reaction_conditions(args, polymerase=polymerase)
    reach, reach_source = _resolve_reach(args, polymerase)

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

    total_sites = sum(p["binding_sites"] for p in per_primer)
    total_bg_sites, selectivity, occupancy_selectivity = _background_totals(
        args, primers, per_primer, total_sites, prefixes, background, conditions
    )
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
        background,
        (total_sites, total_bg_sites, selectivity),
    )
    # Per strain, from a variant table against this reference. Added last so
    # the per-primer records it annotates are already built.
    variant_block = _variant_blocks(
        args, primers, per_primer, cache, (prefixes, genomes, lengths), reach
    )

    result = {
        "primers": per_primer,
        "num_primers": len(primers),
        "fg_coverage": round(overall, 4),
        "per_target_coverage": {k: round(v, 4) for k, v in coverage.items()},
        "extension_reach_bp": reach,
        "extension_reach_source": reach_source,
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
    if variant_block is not None:
        result["per_strain"] = variant_block["per_strain"]
        result["variant_route"] = {
            key: value for key, value in variant_block.items() if key != "per_strain"
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


def _print_per_strain(result):
    """The per-strain rows, the reductions, and the three limits.

    The limits are printed with the figures rather than left to the JSON: the
    route cannot see a site a variant creates, treats an indel as affecting
    every site it overlaps, and says nothing about sequence absent from the
    reference. A reader acting on the intact fraction needs all three.
    """
    strains = result.get("per_strain") or {}
    route = result.get("variant_route") or {}
    if not strains:
        return

    window = route.get("three_prime_window_nt")
    geometry = route.get("geometry") or {}
    print(
        f"\n  Variants placed against {route.get('reference_fasta')} "
        f"(prefix {os.path.basename(str(route.get('reference_prefix')))}, "
        f"{geometry.get('length_bp')} bp, "
        f"{'circular' if geometry.get('circular') else 'linear'})"
    )
    print(f"  Per strain, over {route.get('reference_sites', 0)} reference site(s):")
    for name, record in strains.items():
        if record.get("status") != "measured":
            print(f"      {name:<20s} unavailable: {record.get('unavailable')}")
            continue
        intact = _format_measurement(record.get("intact_sites"), "{:.0f}")
        affected = _format_measurement(record.get("affected_sites"), "{:.0f}")
        not_assessed = _measured(record.get("not_assessed_sites"))
        if not_assessed:
            print(f"      {name:<20s} {not_assessed:.0f} site(s) not assessed")
        fraction = _format_measurement(record.get("intact_site_fraction"), "{:.1%}")
        coverage = _format_measurement(record.get("coverage_on_intact_sites"), "{:.1%}")
        print(
            f"      {name:<20s} intact {intact:>7s}  affected {affected:>7s}  "
            f"intact fraction {fraction}  coverage {coverage}"
        )
        print(
            f"          affected within {window} base(s) of the 3' end: "
            f"{_format_measurement(record.get('affected_three_prime_proximal'), '{:.0f}')}"
            f", further out: "
            f"{_format_measurement(record.get('affected_distal'), '{:.0f}')}"
        )

    print("\n  Worst case over the strains:")
    print(
        f"      worst strain coverage      : "
        f"{_format_measurement(route.get('worst_strain_coverage'), '{:.1%}')}"
    )
    print(
        f"      worst intact site fraction : "
        f"{_format_measurement(route.get('worst_strain_intact_fraction'), '{:.1%}')}"
    )
    print("\n  What this route cannot see:")
    for limit in route.get("limits") or []:
        print(f"      - {limit}")


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
    if result.get("extension_reach_source"):
        print(f"      reach from: {result['extension_reach_source']}")
    for name, cov in result["per_target_coverage"].items():
        print(f"      {os.path.basename(name)}: {cov:.1%}")

    _print_per_reference(result)
    _print_per_strain(result)

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
            f"\n  {len(result['zero_site_primers'])} primer(s) have NO binding sites in the target:"
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
    p.add_argument(
        "--variants",
        metavar="FILE",
        help="Variants of the target strains against this reference: a VCF/BCF, "
        "or a TSV with a header of chrom/pos/ref/alt and one optional 0/1 "
        "column per strain (pos is 1-based, as in VCF). Reports which binding "
        "sites survive in which strain. A site is counted as lost when any "
        "variant falls inside it; no mismatch is scored, and sites a variant "
        "creates are invisible to this route.",
    )
    p.add_argument(
        "--variants-reference",
        metavar="FASTA",
        help="With --variants and more than one foreground genome: the one the "
        "variant table was made against. Must be one of the run's foreground "
        "genomes. Not needed, and not guessed, when there is only one.",
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
