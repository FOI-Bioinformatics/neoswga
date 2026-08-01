"""Iterative set-editing CLI handlers: ``expand-primers``, ``swap-primer``,
``contract-set``, ``rescore-set`` (and their private primer/candidate loaders).

Extracted from cli_unified.py. Shared helpers come from neoswga.cli._common.
"""

import json
import logging
import os
import sys

from neoswga.cli._common import (
    _apply_seed,
    _record_run_manifest,
    collect_primers_from_args,
    merge_args_to_parameter,
    params_command,
    validate_params_json_file,
)
from neoswga.core.position_cache import MissingPositionsError
from neoswga.core.registry import polymerase_names as _polymerase_names

logger = logging.getLogger(__name__)


def _resolve_expansion_candidates(candidates_file, data_dir, quiet=False):
    """Load the candidate pool for expand-primers.

    Prefers an explicit --candidates-file, falling back to the scored
    ``step3_df.csv``. Before --candidates-file existed, step3 was a hard
    precondition, so expanding a set required re-running the whole pipeline
    even when the caller already had a candidate list.
    """
    import pandas as _pd

    from neoswga.core.io_utils import primer_column

    if candidates_file:
        if not os.path.exists(candidates_file):
            logger.error(f"Candidates file not found: {candidates_file}")
            sys.exit(1)
        source = candidates_file
        logger.info(f"Using candidate pool from {candidates_file}")
    else:
        source = os.path.join(data_dir, "step3_df.csv")
        if not os.path.exists(source):
            logger.error(f"Step 3 output not found: {source}")
            logger.error(
                "Run 'neoswga score' first, or pass --candidates-file with your "
                "own candidate list."
            )
            sys.exit(1)

    if not quiet:
        logger.info(f"Loading candidates from {source}...")

    df = _pd.read_csv(source)
    try:
        primer_col = primer_column(df)
    except KeyError:
        logger.error(
            f"Invalid candidates file: missing 'primer' or 'seq' column. "
            f"Found: {list(df.columns)}"
        )
        sys.exit(1)
    return df[primer_col].tolist()


def run_expand_primers(args):
    """
    Expand existing primer set with additional primers.

    Useful for iterative design where some primers have been validated
    experimentally and user wants to add more to fill coverage gaps.
    """
    import json as json_module

    import pandas as pd

    from neoswga.core import parameter
    from neoswga.core import pipeline as core_pipeline
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.primer_expansion import PrimerExpander

    quiet = getattr(args, "quiet", False)
    _apply_seed(args)

    if not quiet:
        logger.info("Primer Set Expansion")
        logger.info("=" * 60)

    # Load and validate primers using consolidated helper
    fixed_primers = collect_primers_from_args(
        cli_primers=args.fixed_primers,
        primers_file=args.fixed_primers_file,
        name="fixed primer",
        allow_empty=True,  # Fixed primers can be empty for initial design
    )

    failed_primers = collect_primers_from_args(
        cli_primers=args.failed_primers,
        primers_file=args.failed_primers_file,
        name="failed primer",
        allow_empty=True,  # Failed primers are optional
    )

    if not quiet:
        logger.info(f"Fixed primers: {len(fixed_primers)}")
        logger.info(f"Failed primers (excluded): {len(failed_primers)}")
        logger.info(f"Target new primers: {args.num_new}")

    # Set json_file for pipeline initialization
    merge_args_to_parameter(args, parameter, ["json_file"])

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    try:
        # Initialize pipeline to get genome info
        if not quiet:
            logger.info("Initializing pipeline...")
        core_pipeline._initialize()

        fg_prefixes = core_pipeline.fg_prefixes
        bg_prefixes = core_pipeline.bg_prefixes
        fg_seq_lengths = core_pipeline.fg_seq_lengths
        bg_seq_lengths = core_pipeline.bg_seq_lengths

        candidates = _resolve_expansion_candidates(
            getattr(args, "candidates_file", None),
            getattr(parameter, "data_dir", "."),
            quiet=quiet,
        )

        if not quiet:
            logger.info(f"Candidate pool: {len(candidates)} primers")

        # Initialize position cache with candidates and fixed primers
        all_primers = list(set(candidates + fixed_primers))
        if not quiet:
            logger.info("Loading position data...")

        cache = PositionCache(fg_prefixes, all_primers)

        # Create expander
        expander = PrimerExpander(
            position_cache=cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
        )

        # Build target gaps: in-silico gaps from the fixed set, optionally
        # merged with low sequencing-depth regions from a mapped BAM.
        bam_gaps_list = None
        if getattr(args, "bam", None):
            try:
                from neoswga.core.bam_coverage import bam_gaps
            except Exception as e:
                logger.error(str(e))
                sys.exit(1)
            aliases = {}
            for item in getattr(args, "contig_alias", None) or []:
                if "=" in item:
                    k, v = item.split("=", 1)
                    aliases[k] = v
            fg_circular = bool(getattr(parameter, "fg_circular", False))
            try:
                bam_gaps_list = bam_gaps(
                    args.bam,
                    fg_prefixes,
                    fg_seq_lengths,
                    min_depth=args.min_depth,
                    min_gap_size=args.min_gap_size,
                    circular=fg_circular,
                    contig_aliases=aliases or None,
                )
            except RuntimeError as e:  # pysam missing
                logger.error(str(e))
                sys.exit(1)
            if not quiet:
                logger.info(f"BAM low-depth gaps: {len(bam_gaps_list)}")
            # Actionable hint for the most common BAM failure: contig names in
            # the BAM header (e.g. 'chr1') not matching the fg prefixes.
            if not bam_gaps_list:
                try:
                    import pysam

                    with pysam.AlignmentFile(args.bam, "rb") as _bam:
                        _bam_contigs = list(_bam.references)
                    logger.warning(
                        "No BAM gaps produced. If this is unexpected, the BAM "
                        "contig names may not match the foreground prefixes.\n"
                        f"  BAM contigs: {_bam_contigs[:10]}\n"
                        f"  fg prefixes: {[os.path.basename(p) for p in fg_prefixes]}\n"
                        "  Map them with --contig-alias FG=BAMCONTIG (repeatable)."
                    )
                except Exception:
                    pass

        target_gaps = expander.identify_gaps(
            fixed_primers,
            min_gap_size=args.min_gap_size,
            extra_gaps=bam_gaps_list,
            merge=True,
        )

        os.makedirs(args.output, exist_ok=True)
        if target_gaps:
            from neoswga.core.export import export_gaps_to_bed

            export_gaps_to_bed(target_gaps, os.path.join(args.output, "merged_gaps.bed"))

        # Run expansion (focus candidates on the gaps when we have any)
        result = expander.expand(
            candidates=candidates,
            fixed_primers=fixed_primers,
            failed_primers=failed_primers,
            target_new=args.num_new,
            optimization_method=args.optimization_method,
            verbose=not quiet,
            target_gaps=target_gaps or None,
        )

        # Save results
        result_file = os.path.join(args.output, "expansion_result.json")
        with open(result_file, "w") as f:
            json_module.dump(result.to_dict(), f, indent=2)

        # Save primers as CSV
        primers_file = os.path.join(args.output, "expanded_primers.csv")
        pd.DataFrame(
            {
                "primer": result.combined_set,
                "type": ["fixed"] * result.n_fixed + ["new"] * result.n_new,
            }
        ).to_csv(primers_file, index=False)

        # Print summary
        if not quiet:
            logger.info("")
            logger.info("=" * 60)
            logger.info("EXPANSION COMPLETE")
            logger.info("=" * 60)
            logger.info(f"Fixed primers: {result.n_fixed}")
            logger.info(f"New primers: {result.n_new}")
            logger.info(f"Total set: {result.n_total}")
            logger.info("")
            logger.info(f"Coverage: {result.coverage_before:.1%} -> {result.coverage_after:.1%}")
            logger.info(f"Improvement: {result.predicted_improvement:.2f}x")
            logger.info(f"Gaps remaining: {result.gaps_remaining}")
            logger.info("")
            logger.info("New primers:")
            for primer in result.new_primers:
                logger.info(f"  {primer}")
            logger.info("")
            logger.info(f"Results saved to: {args.output}")
            logger.info("")
            logger.info("Next steps:")
            logger.info(f"  neoswga interpret -d {args.output}    # assess the expanded set")
            logger.info(f"  neoswga report -d {args.output} --level full   # full report")
            logger.info(f"  neoswga export -d {args.output} -o order/      # order for synthesis")

        _record_run_manifest(
            "expand-primers",
            args,
            parameter,
            input_files=[os.path.join(args.output, "expanded_primers.csv")],
        )

    except SystemExit:
        raise
    except (FileNotFoundError, KeyError) as e:
        logger.error(f"Primer expansion failed: {e}")
        logger.error("Check that 'neoswga score' has run and params.json is valid.")
        sys.exit(1)
    except RuntimeError as e:
        logger.error(f"Primer expansion failed: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Primer expansion failed: {e}")
        if logger.level <= logging.DEBUG:
            import traceback

            traceback.print_exc()
        else:
            logger.error("Run with --verbose for the full traceback.")
        sys.exit(1)


def _load_primer_list(cli_primers, primers_file, name="primer"):
    """Collect primers from --primers or --primers-file."""
    primers = list(cli_primers or [])
    if primers_file:
        with open(primers_file, "r") as fh:
            for line in fh:
                p = line.strip().split()[0] if line.strip() else ""
                if p and not p.startswith("#"):
                    primers.append(p)
    primers = [p.upper().strip() for p in primers if p.strip()]
    if not primers:
        logger.error(f"No {name} sequences provided")
        sys.exit(1)
    return primers


def _load_candidate_pool(candidates_file, data_dir):
    """Load a candidate-pool primer list from --candidates-file or data_dir/step2_df.csv."""
    import pandas as pd

    candidates = []
    if candidates_file and os.path.exists(candidates_file):
        with open(candidates_file, "r") as fh:
            for line in fh:
                p = line.strip().split()[0] if line.strip() else ""
                if p and not p.startswith("#"):
                    candidates.append(p)
    else:
        step2 = os.path.join(data_dir or ".", "step2_df.csv")
        if os.path.isfile(step2):
            df = pd.read_csv(step2)
            if "primer" in df.columns:
                candidates = df["primer"].astype(str).tolist()
    return [c.upper().strip() for c in candidates if c.strip()]


@params_command(seed=True)
def run_swap_primer(args):
    """Swap under-performing primers in a set with better candidates."""
    import json as _json

    # Ensure parameter globals are populated from params.json
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter
    from neoswga.core.dimer_network_analyzer import create_dimer_network_analyzer
    from neoswga.core.reaction_conditions import ReactionConditions

    pipeline_mod._initialized = False
    pipeline_mod._initialize()

    current = _load_primer_list(args.primers, args.primers_file, name="current primer")
    data_dir = getattr(parameter, "data_dir", ".")
    candidates = _load_candidate_pool(args.candidates_file, data_dir)
    if not candidates:
        logger.error(
            "No candidate pool available. Provide --candidates-file or "
            "run 'neoswga filter' first to generate step2_df.csv."
        )
        sys.exit(1)

    # Blacklist guard (review I4): when params.json carries bl_genomes,
    # double-check the candidate pool has no primer that k-mer-hits the
    # blacklist. Prevents swap-primer from re-injecting primers a user
    # loads from an external --candidates-file that was not passed through
    # the neoswga filter step.
    bl_prefixes = list(getattr(parameter, "bl_prefixes", []) or [])
    bl_lengths = list(getattr(parameter, "bl_seq_lengths", []) or [])
    forbidden_primers: set = set()
    if bl_prefixes and bl_lengths:
        try:
            from neoswga.core.pipeline import _filter_blacklist_penalty

            _, bl_freqs = _filter_blacklist_penalty(
                candidates,
                bl_prefixes,
                bl_lengths,
                max_bl_freq=getattr(parameter, "max_bl_freq", 0.0) or 0.0,
            )
            forbidden_primers = {
                p
                for p, f in zip(candidates, bl_freqs)
                if f > (getattr(parameter, "max_bl_freq", 0.0) or 0.0)
            }
            if forbidden_primers:
                logger.warning(
                    f"Excluded {len(forbidden_primers)} blacklist-hitting "
                    f"primer(s) from the swap-primer candidate pool. "
                    f"Pass --candidates-file with a filtered pool to avoid "
                    f"this at load time."
                )
                candidates = [c for c in candidates if c not in forbidden_primers]
        except Exception as e:
            logger.debug(f"Blacklist guard skipped: {e}")

    conditions = ReactionConditions(
        temp=getattr(parameter, "reaction_temp", 30.0) or 30.0,
        polymerase=getattr(parameter, "polymerase", "phi29"),
        na_conc=getattr(parameter, "na_conc", 50.0),
        mg_conc=getattr(parameter, "mg_conc", 10.0),
        dmso_percent=getattr(parameter, "dmso_percent", 0.0),
        betaine_m=getattr(parameter, "betaine_m", 0.0),
        trehalose_m=getattr(parameter, "trehalose_m", 0.0),
        formamide_percent=getattr(parameter, "formamide_percent", 0.0),
        ethanol_percent=getattr(parameter, "ethanol_percent", 0.0),
        urea_m=getattr(parameter, "urea_m", 0.0),
        tmac_m=getattr(parameter, "tmac_m", 0.0),
    )
    analyzer = create_dimer_network_analyzer(conditions)

    # Snapshot "before" state: analyze_primer_set returns (metrics, profiles, matrix)
    before_metrics, _before_profiles, _before_matrix = analyzer.analyze_primer_set(current)
    # Greedy swap
    improved, after_metrics = analyzer.optimize_set_greedy(
        list(current),
        candidates,
        max_iterations=args.max_swaps,
    )

    swapped_in = [p for p in improved if p not in current]
    swapped_out = [p for p in current if p not in improved]
    result = {
        "before_set": list(current),
        "after_set": list(improved),
        "swapped_in": swapped_in,
        "swapped_out": swapped_out,
        "num_swaps": len(swapped_in),
        "metrics_before": {
            "total_interactions": before_metrics.total_interactions,
            "mean_severity": before_metrics.mean_severity,
            "max_severity": before_metrics.max_severity,
            "num_hubs": before_metrics.num_hub_primers,
        },
        "metrics_after": {
            "total_interactions": after_metrics.total_interactions,
            "mean_severity": after_metrics.mean_severity,
            "max_severity": after_metrics.max_severity,
            "num_hubs": after_metrics.num_hub_primers,
        },
    }

    output_json = _json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(output_json + "\n")
        if not args.quiet:
            logger.info(f"swap-primer results written to {args.output}")
    else:
        print(output_json)


def _realistic_reach():
    """Per-primer coverage reach for the configured polymerase.

    Coverage is scored on the REALISTIC reach (~3 kb for phi29), not
    single-molecule processivity (~70 kb). contract-set and rescore-set both
    used processivity, which on a multi-Mb target let one binding site cover a
    large fraction of the genome -- turning --min-coverage into an on/off switch
    rather than a threshold.
    """
    from neoswga.core import parameter
    from neoswga.core.coverage import polymerase_extension_reach

    return polymerase_extension_reach(
        getattr(parameter, "polymerase", "phi29") or "phi29",
        coverage_metric="realistic",
    )


def _prefix_lengths(parameter, kind):
    """Resolve (prefixes, sequence lengths) for 'fg' or 'bg' from module state.

    `get_params` publishes `<kind>_seq_lengths` as a module global, deriving it
    from the FASTAs when params.json omits the key. The `_json_data` fallback is
    kept only for a caller that populated the raw JSON without running
    get_params.

    That fallback used to be the PRIMARY path, which was the bug: the derived
    value never lands in `_json_data`, so a params.json without an explicit
    `fg_seq_lengths` -- which is every config `neoswga init` writes -- made
    contract-set and rescore-set exit with "fg_seq_lengths missing from params;
    cannot compute coverage" while the pipeline itself ran fine. Only the two
    shipped examples declare the key by hand, which is why it went unseen.
    """
    prefixes = list(getattr(parameter, f"{kind}_prefixes", []) or [])
    lengths = list(getattr(parameter, f"{kind}_seq_lengths", []) or [])
    if not lengths:
        raw = getattr(parameter, "_json_data", {}) or {}
        lengths = list(raw.get(f"{kind}_seq_lengths", []) or [])
    return prefixes, lengths


@params_command(seed=True)
def run_contract_set(args):
    """Greedy leave-one-out contraction that keeps coverage above threshold."""
    import json as _json

    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter
    from neoswga.core.position_cache import PositionCache

    pipeline_mod._initialized = False
    pipeline_mod._initialize()

    current = _load_primer_list(args.primers, args.primers_file, name="primer")
    fg_prefixes, fg_lengths = _prefix_lengths(parameter, "fg")
    if not fg_prefixes:
        logger.error("No fg_prefixes available; cannot compute coverage.")
        sys.exit(1)
    if not fg_lengths:
        logger.error("fg_seq_lengths missing from params; cannot compute coverage.")
        sys.exit(1)

    cache = PositionCache(fg_prefixes, current)

    def _coverage(primers):
        # Union of binding positions across all fg prefixes, binned at 1 bp granularity.
        total_genome = sum(fg_lengths) if fg_lengths else 0
        if not total_genome:
            return 0.0
        covered = 0
        extension = _realistic_reach()
        import numpy as _np

        for prefix, length in zip(fg_prefixes, fg_lengths or [0] * len(fg_prefixes)):
            if length <= 0:
                continue
            occupied = _np.zeros(length, dtype=bool)
            for primer in primers:
                positions = cache.get_positions(prefix, primer, "both")
                for pos in positions:
                    start = max(0, int(pos) - extension)
                    end = min(length, int(pos) + extension)
                    occupied[start:end] = True
            covered += int(occupied.sum())
        return covered / total_genome if total_genome else 0.0

    baseline = _coverage(current)
    if not args.quiet:
        logger.info(f"Baseline coverage with {len(current)} primers: {baseline:.1%}")

    # Quality-weighted removal ranking (Phase 12C).
    # When multiple primers can be dropped without falling below the coverage
    # threshold, prefer to drop the one with the worst quality: highest
    # dimer interactions with the rest of the set, worst Tm fit under the
    # reaction conditions, lowest per-primer strand alternation. Weights:
    #   w_cov   = 0.40  (penalize loss of coverage contribution)
    #   w_dimer = 0.25  (prefer removing dimer-prone primers)
    #   w_tm    = 0.20  (prefer removing primers with poor Tm fit)
    #   w_bg    = 0.15  (prefer removing primers with high bg frequency)
    from neoswga.core.dimer_network_analyzer import create_dimer_network_analyzer
    from neoswga.core.integrated_quality_scorer import IntegratedQualityScorer
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(
        temp=getattr(parameter, "reaction_temp", 30.0) or 30.0,
        polymerase=getattr(parameter, "polymerase", "phi29"),
        na_conc=getattr(parameter, "na_conc", 50.0),
        mg_conc=getattr(parameter, "mg_conc", 10.0),
        dmso_percent=getattr(parameter, "dmso_percent", 0.0),
        betaine_m=getattr(parameter, "betaine_m", 0.0),
        trehalose_m=getattr(parameter, "trehalose_m", 0.0),
        formamide_percent=getattr(parameter, "formamide_percent", 0.0),
        ethanol_percent=getattr(parameter, "ethanol_percent", 0.0),
        urea_m=getattr(parameter, "urea_m", 0.0),
        tmac_m=getattr(parameter, "tmac_m", 0.0),
    )
    scorer = IntegratedQualityScorer(conditions=conditions)
    dimer_analyzer = create_dimer_network_analyzer(conditions=conditions)

    W_COV, W_DIMER, W_TM, W_BG = 0.40, 0.25, 0.20, 0.15
    min_tm = getattr(parameter, "min_tm", 15.0) or 15.0
    max_tm = getattr(parameter, "max_tm", 45.0) or 45.0
    target_tm = (min_tm + max_tm) / 2.0
    tm_window = max(1.0, (max_tm - min_tm) / 2.0)

    def _deficit_score(primer: str, primers_in_set: list, set_cov: float) -> float:
        # 1. Coverage contribution: how much does removing this primer cost?
        without = [p for p in primers_in_set if p != primer]
        cov_without = _coverage(without) if without else 0.0
        cov_contribution = max(0.0, set_cov - cov_without)

        # 2. Dimer risk: pairwise severity with rest of set
        if len(primers_in_set) > 1:
            try:
                bmetrics, profiles, _ = dimer_analyzer.analyze_primer_set(
                    primers_in_set,
                    verbose=False,
                )
                prof = profiles.get(primer)
                dimer_risk = float(prof.mean_severity) if prof else 0.0
            except Exception:
                dimer_risk = 0.0
        else:
            dimer_risk = 0.0

        # 3. Tm fit under conditions
        try:
            gc = sum(1 for b in primer if b in "GC") / max(len(primer), 1)
            tm = scorer._primer_tm(primer, gc, len(primer))
            tm_deficit = abs(tm - target_tm) / tm_window
        except Exception:
            tm_deficit = 0.0

        # 4. Background frequency (approximate from position cache counts
        # across bg_prefixes if available). Zero when no background
        # data present.
        bg_prefixes = list(getattr(parameter, "bg_prefixes", []) or [])
        bg_lengths = list(getattr(parameter, "bg_seq_lengths", []) or [])
        if not bg_lengths and getattr(parameter, "_json_data", None):
            bg_lengths = list(parameter._json_data.get("bg_seq_lengths", []) or [])
        bg_freq = 0.0
        if bg_prefixes and bg_lengths:
            try:
                bg_cache = PositionCache(bg_prefixes, [primer])
                hits = 0
                total = 0
                for p, L in zip(bg_prefixes, bg_lengths):
                    hits += len(bg_cache.get_positions(p, primer, "both"))
                    total += L
                bg_freq = hits / total if total else 0.0
            except Exception:
                bg_freq = 0.0

        # Higher deficit = more reason to remove.
        # We INVERT cov_contribution — losing a lot of coverage means low
        # removal preference.
        deficit = (
            -W_COV * cov_contribution
            + W_DIMER * dimer_risk
            + W_TM * min(tm_deficit, 1.0)
            + W_BG * min(bg_freq * 1000, 1.0)  # scale bg_freq for visibility
        )
        return deficit

    kept = list(current)
    removed = []
    while True:
        set_cov = _coverage(kept)
        # Rank candidates by deficit score (higher = remove first)
        ranked = sorted(
            kept,
            key=lambda p: _deficit_score(p, kept, set_cov),
            reverse=True,
        )

        dropped_this_round = False
        for candidate in ranked:
            trial = [p for p in kept if p != candidate]
            if not trial:
                break
            trial_cov = _coverage(trial)
            if trial_cov >= args.min_coverage:
                kept = trial
                removed.append(
                    {
                        "primer": candidate,
                        "resulting_coverage": trial_cov,
                    }
                )
                if not args.quiet:
                    logger.info(
                        f"Removed {candidate} (quality-weighted) -> coverage {trial_cov:.1%}"
                    )
                dropped_this_round = True
                break
        if not dropped_this_round:
            break

    result = {
        "original_set": list(current),
        "contracted_set": kept,
        "removed_primers": [r["primer"] for r in removed],
        "removal_trace": removed,
        "baseline_coverage": baseline,
        "final_coverage": _coverage(kept),
        "min_coverage_threshold": args.min_coverage,
    }

    output_json = _json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(output_json + "\n")
        if not args.quiet:
            logger.info(f"contract-set results written to {args.output}")
    else:
        print(output_json)


@params_command()
def run_rescore_set(args):
    """Rescore an existing primer set under arbitrary reaction conditions."""
    import json as _json

    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter
    from neoswga.core.integrated_quality_scorer import IntegratedQualityScorer
    from neoswga.core.reaction_conditions import ReactionConditions

    pipeline_mod._initialized = False
    pipeline_mod._initialize()

    primers = _load_primer_list(args.primers, args.primers_file, name="primer")

    def _pick(attr, fallback):
        cli_val = getattr(args, attr, None)
        if cli_val is not None:
            return cli_val
        pv = getattr(parameter, attr, None)
        return pv if pv is not None else fallback

    # A temperature configured for one polymerase is not meaningful for another:
    # params.json carrying phi29's 30 C makes `--polymerase bst` fail validation
    # outright (bst runs 50-72 C), which is the most obvious use of a command
    # whose purpose is rescoring under a different chemistry. When the
    # polymerase is overridden on the command line and the temperature is not,
    # take the new polymerase's optimum rather than the old one's setting.
    polymerase = _pick("polymerase", "phi29")
    if getattr(args, "polymerase", None) and getattr(args, "reaction_temp", None) is None:
        reaction_temp = parameter.default_reaction_temp(polymerase)
        logger.info(
            "Using %s's optimal temperature %.1f C (params.json specifies %s C, "
            "which belongs to the previous polymerase). Pass --reaction-temp to override.",
            polymerase,
            reaction_temp,
            getattr(parameter, "reaction_temp", None),
        )
    else:
        reaction_temp = _pick("reaction_temp", 30.0) or 30.0

    conditions = ReactionConditions(
        temp=reaction_temp,
        polymerase=polymerase,
        na_conc=_pick("na_conc", 50.0),
        mg_conc=_pick("mg_conc", 10.0),
        dmso_percent=_pick("dmso_percent", 0.0),
        betaine_m=_pick("betaine_m", 0.0),
        trehalose_m=_pick("trehalose_m", 0.0),
        formamide_percent=_pick("formamide_percent", 0.0),
        ethanol_percent=_pick("ethanol_percent", 0.0),
        urea_m=_pick("urea_m", 0.0),
        tmac_m=_pick("tmac_m", 0.0),
    )

    scorer = IntegratedQualityScorer(conditions=conditions)
    per_primer = []
    for p in primers:
        q = scorer.score_primer(p)
        per_primer.append(
            {
                "primer": p,
                "overall_score": float(q.overall_score),
                "strand_bias_score": float(q.strand_bias_score),
                "dimer_score": float(q.dimer_score),
                "three_prime_score": float(q.three_prime_score),
                "complexity_score": float(q.complexity_score),
                "thermo_score": float(q.thermo_score),
                "passes_all": bool(q.passes_all),
                "failure_reasons": [str(r) for r in q.failure_reasons],
                "tm_effective_c": float(
                    scorer._primer_tm(
                        p,
                        sum(1 for b in p if b in "GC") / max(len(p), 1),
                        len(p),
                    )
                ),
            }
        )

    _, set_score = scorer.analyze_primer_set(primers)

    # Coverage + background computation under the chosen reaction conditions
    # (Phase 12D). Extension reach follows polymerase processivity, so the
    # same set can look 90% covered under phi29 (70 kb reach) and 78% under
    # klenow (10 kb reach). Per-target coverage is reported separately so a
    # multi-target user can see uneven distribution.
    coverage_block: dict = {}
    try:
        from neoswga.core.coverage import polymerase_extension_reach
        from neoswga.core.position_cache import PositionCache

        # Realistic per-primer reach, matching how the optimizer scores
        # fg_coverage. This used single-molecule processivity (~70 kb), which
        # overstated coverage by roughly an order of magnitude.
        extension = polymerase_extension_reach(conditions.polymerase, coverage_metric="realistic")

        fg_prefixes, fg_lengths = _prefix_lengths(parameter, "fg")

        bg_prefixes, bg_lengths = _prefix_lengths(parameter, "bg")

        all_prefixes = fg_prefixes + bg_prefixes
        cache = PositionCache(all_prefixes, primers) if all_prefixes else None

        def _coverage_across(prefix_list, length_list):
            if not prefix_list or not length_list or cache is None:
                return 0.0, {}
            import numpy as _np

            per_prefix: dict = {}
            total_cov = 0
            total_len = 0
            for prefix, length in zip(prefix_list, length_list):
                if length <= 0:
                    per_prefix[prefix] = 0.0
                    continue
                occupied = _np.zeros(length, dtype=bool)
                for primer in primers:
                    try:
                        positions = cache.get_positions(prefix, primer, "both")
                    except Exception:
                        continue
                    for pos in positions:
                        start = max(0, int(pos) - extension)
                        end = min(length, int(pos) + extension)
                        occupied[start:end] = True
                covered = int(occupied.sum())
                per_prefix[prefix] = covered / length if length else 0.0
                total_cov += covered
                total_len += length
            agg = total_cov / total_len if total_len else 0.0
            return agg, per_prefix

        fg_cov, fg_per_target = _coverage_across(fg_prefixes, fg_lengths)
        bg_cov, bg_per_prefix = _coverage_across(bg_prefixes, bg_lengths)

        coverage_block = {
            "fg_coverage": float(fg_cov),
            "per_target_coverage": {k: float(v) for k, v in fg_per_target.items()},
            "bg_coverage": float(bg_cov),
            "per_background_coverage": {k: float(v) for k, v in bg_per_prefix.items()},
            "selectivity_ratio": float((fg_cov / bg_cov) if bg_cov > 1e-12 else (fg_cov * 1000.0)),
            "extension_reach_bp": int(extension),
        }
    except MissingPositionsError as e:
        # Never swallow this one. It means the coverage number would be a
        # meaningless zero rather than a low result, which is exactly the
        # failure this command is used to diagnose.
        logger.error(
            "Cannot compute coverage: %s\n"
            "Pass --genome to scan for these primers directly, or run "
            "'neoswga count-kmers' so they are indexed.",
            e,
        )
        raise
    except Exception as e:
        logger.error(
            "rescore-set coverage block failed: %s. Reported coverage is " "incomplete.",
            e,
            exc_info=True,
        )
        coverage_block = {"error": str(e)}

    result = {
        "conditions": {
            "temp": float(conditions.temp),
            "polymerase": conditions.polymerase,
            "na_conc": float(conditions.na_conc),
            "mg_conc": float(conditions.mg_conc),
            "dmso_percent": float(conditions.dmso_percent),
            "betaine_m": float(conditions.betaine_m),
            "trehalose_m": float(conditions.trehalose_m),
            "formamide_percent": float(conditions.formamide_percent),
            "ethanol_percent": float(conditions.ethanol_percent),
            "urea_m": float(conditions.urea_m),
            "tmac_m": float(conditions.tmac_m),
        },
        "primers": per_primer,
        "set_summary": {
            "mean_overall_score": float(set_score.mean_overall_score),
            "min_overall_score": float(set_score.min_overall_score),
            "set_dimer_score": float(set_score.set_dimer_score),
            "passes": bool(set_score.passes),
        },
        "coverage": coverage_block,
    }

    output_json = _json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(output_json + "\n")
        if not args.quiet:
            logger.info(f"rescore-set results written to {args.output}")
    else:
        print(output_json)


def add_parsers(subparsers):
    """Register this group's subcommands on the shared subparsers object.

    Called by neoswga.cli_unified.create_parser(). Extracted from the former
    monolithic create_parser() so each command group owns its argparse setup
    next to its handlers.
    """
    import argparse  # noqa: F401  (used by some command blocks)

    from neoswga.cli._common import add_common_options  # noqa: F401

    expand_parser = subparsers.add_parser(
        "expand-primers",
        help="Expand existing primer set with additional primers to fill coverage gaps",
    )
    expand_parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    expand_parser.add_argument(
        "--fixed-primers",
        nargs="+",
        required=True,
        help="Primer sequences to keep (already validated)",
    )
    expand_parser.add_argument(
        "--fixed-primers-file", help="File with fixed primers (one per line)"
    )
    expand_parser.add_argument(
        "--failed-primers", nargs="+", help="Primer sequences to exclude (failed in wet lab)"
    )
    expand_parser.add_argument(
        "--candidates-file",
        help="CSV of candidate primers to draw new ones from. Without this, "
        "expand-primers requires data_dir/step3_df.csv from a prior 'score' run; "
        "with it, that becomes a fallback rather than a precondition.",
    )
    expand_parser.add_argument(
        "--failed-primers-file", help="File with failed primers to exclude (one per line)"
    )
    expand_parser.add_argument(
        "--num-new", type=int, default=6, help="Number of new primers to add (default: 6)"
    )
    expand_parser.add_argument(
        "--optimization-method",
        default="hybrid",
        choices=["hybrid", "dominating-set", "network", "background-aware", "ensemble"],
        help="Optimization method (default: hybrid)",
    )
    expand_parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for reproducible expansion.",
    )
    expand_parser.add_argument(
        "--bam",
        help="Mapped BAM of real sequencing reads against the "
        "target genome. Low-depth regions are merged with "
        "in-silico coverage gaps and the candidate pool is "
        "focused on primers that bind inside the gaps. "
        "Requires the [bam] extra (pip install 'neoswga[bam]').",
    )
    expand_parser.add_argument(
        "--min-depth",
        type=int,
        default=5,
        help="Sequencing depth below which a base counts as a " "BAM coverage gap (default: 5).",
    )
    expand_parser.add_argument(
        "--min-gap-size",
        type=int,
        default=10000,
        help="Minimum gap length to act on, bp (default: 10000).",
    )
    expand_parser.add_argument(
        "--contig-alias",
        action="append",
        default=None,
        metavar="FG=BAMCONTIG",
        help="Map a foreground prefix/basename to a BAM contig "
        "name when they differ. Repeatable.",
    )
    expand_parser.add_argument(
        "--output", "-o", required=True, help="Output directory for expanded primer set"
    )
    expand_parser.add_argument(
        "--quiet", "-q", action="store_true", help="Suppress progress output"
    )

    # Analyze coverage - read-only inspection of in-silico + BAM gaps

    swap_parser = subparsers.add_parser(
        "swap-primer",
        help="Replace under-performing primers in an existing set with better "
        "candidates from a pool (greedy dimer-aware optimization).",
    )
    swap_parser.add_argument("-j", "--json-file", required=True, help="Parameters JSON file")
    swap_parser.add_argument("--primers", nargs="+", required=True, help="Current primer set")
    swap_parser.add_argument("--primers-file", help="File with current primers (one per line)")
    swap_parser.add_argument(
        "--candidates-file",
        help="File with candidate pool for replacement, "
        "one primer per line. Defaults to step2_df.csv "
        "primers if present in data_dir.",
    )
    swap_parser.add_argument(
        "--max-swaps", type=int, default=3, help="Maximum replacements to attempt (default: 3)"
    )
    swap_parser.add_argument(
        "--output",
        "-o",
        help="Output JSON with before/after set and per-swap " "deltas (default: stdout)",
    )
    swap_parser.add_argument("--quiet", "-q", action="store_true")
    swap_parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducible swaps."
    )

    # Contract set: remove redundant primers while keeping coverage above threshold

    contract_parser = subparsers.add_parser(
        "contract-set",
        help="Remove redundant primers from a set while keeping coverage above "
        "--min-coverage. Greedy leave-one-out.",
    )
    contract_parser.add_argument("-j", "--json-file", required=True)
    contract_parser.add_argument("--primers", nargs="+", required=True)
    contract_parser.add_argument("--primers-file")
    contract_parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.70,
        help="Minimum foreground coverage to keep (default: 0.70)",
    )
    contract_parser.add_argument("--output", "-o")
    contract_parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducible contraction."
    )
    contract_parser.add_argument("--quiet", "-q", action="store_true")

    # Rescore set under arbitrary conditions

    rescore_parser = subparsers.add_parser(
        "rescore-set",
        help="Rescore an existing primer set under arbitrary reaction "
        "conditions (polymerase + additives) without re-optimizing.",
    )
    rescore_parser.add_argument("-j", "--json-file", required=True)
    rescore_parser.add_argument("--primers", nargs="+", required=True)
    rescore_parser.add_argument("--primers-file")
    rescore_parser.add_argument("--polymerase", choices=_polymerase_names())
    rescore_parser.add_argument("--reaction-temp", type=float)
    rescore_parser.add_argument("--mg-conc", type=float)
    rescore_parser.add_argument("--na-conc", type=float)
    rescore_parser.add_argument("--dmso-percent", type=float)
    rescore_parser.add_argument("--betaine-m", type=float)
    rescore_parser.add_argument("--trehalose-m", type=float)
    rescore_parser.add_argument("--formamide-percent", type=float)
    rescore_parser.add_argument("--ethanol-percent", type=float)
    rescore_parser.add_argument("--urea-m", type=float)
    rescore_parser.add_argument("--tmac-m", type=float)
    rescore_parser.add_argument("--output", "-o")
    rescore_parser.add_argument("--quiet", "-q", action="store_true")

    # Predict efficiency - unified confidence score for primer set
