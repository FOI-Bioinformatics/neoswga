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
    validate_params_json_file,
)

logger = logging.getLogger(__name__)


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

        # Load candidate primers from scored step3 output
        data_dir = getattr(parameter, "data_dir", ".")
        step3_file = os.path.join(data_dir, "step3_df.csv")

        if not os.path.exists(step3_file):
            logger.error(f"Step 3 output not found: {step3_file}")
            logger.error("Run 'neoswga score' first to generate candidate primers.")
            sys.exit(1)

        if not quiet:
            logger.info(f"Loading candidates from {step3_file}...")

        df = pd.read_csv(step3_file)
        from neoswga.core.io_utils import primer_column

        try:
            primer_col = primer_column(df)
        except KeyError:
            logger.error(
                f"Invalid step3 file: missing 'primer' or 'seq' column. Found: {list(df.columns)}"
            )
            sys.exit(1)
        candidates = df[primer_col].tolist()

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


def run_swap_primer(args):
    """Swap under-performing primers in a set with better candidates."""
    import json as _json

    from neoswga.core import parameter
    from neoswga.core.dimer_network_analyzer import create_dimer_network_analyzer
    from neoswga.core.reaction_conditions import ReactionConditions

    _apply_seed(args)
    validate_params_json_file(args.json_file)
    merge_args_to_parameter(args, parameter, ["json_file"])
    # Ensure parameter globals are populated from params.json
    import neoswga.core.pipeline as pipeline_mod

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


def run_contract_set(args):
    """Greedy leave-one-out contraction that keeps coverage above threshold."""
    import json as _json

    from neoswga.core import parameter
    from neoswga.core.position_cache import PositionCache

    _apply_seed(args)
    validate_params_json_file(args.json_file)
    merge_args_to_parameter(args, parameter, ["json_file"])
    import neoswga.core.pipeline as pipeline_mod

    pipeline_mod._initialized = False
    pipeline_mod._initialize()

    current = _load_primer_list(args.primers, args.primers_file, name="primer")
    fg_prefixes = list(getattr(parameter, "fg_prefixes", []) or [])
    # fg_seq_lengths flows through _json_data / data dict rather than the module
    # globals, so fall back to _json_data if the global is empty.
    fg_lengths = list(getattr(parameter, "fg_seq_lengths", []) or [])
    if not fg_lengths:
        raw = getattr(parameter, "_json_data", {}) or {}
        fg_lengths = list(raw.get("fg_seq_lengths", []) or [])
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
        extension = 70000
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


def run_rescore_set(args):
    """Rescore an existing primer set under arbitrary reaction conditions."""
    import json as _json

    from neoswga.core import parameter
    from neoswga.core.integrated_quality_scorer import IntegratedQualityScorer
    from neoswga.core.reaction_conditions import ReactionConditions

    validate_params_json_file(args.json_file)
    merge_args_to_parameter(args, parameter, ["json_file"])
    import neoswga.core.pipeline as pipeline_mod

    pipeline_mod._initialized = False
    pipeline_mod._initialize()

    primers = _load_primer_list(args.primers, args.primers_file, name="primer")

    def _pick(attr, fallback):
        cli_val = getattr(args, attr, None)
        if cli_val is not None:
            return cli_val
        pv = getattr(parameter, attr, None)
        return pv if pv is not None else fallback

    conditions = ReactionConditions(
        temp=_pick("reaction_temp", 30.0) or 30.0,
        polymerase=_pick("polymerase", "phi29"),
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
        from neoswga.core.position_cache import PositionCache
        from neoswga.core.reaction_conditions import get_polymerase_processivity

        try:
            extension = int(get_polymerase_processivity(conditions.polymerase))
        except Exception:
            extension = 70_000

        fg_prefixes = list(getattr(parameter, "fg_prefixes", []) or [])
        fg_lengths = list(getattr(parameter, "fg_seq_lengths", []) or [])
        if not fg_lengths:
            raw = getattr(parameter, "_json_data", {}) or {}
            fg_lengths = list(raw.get("fg_seq_lengths", []) or [])

        bg_prefixes = list(getattr(parameter, "bg_prefixes", []) or [])
        bg_lengths = list(getattr(parameter, "bg_seq_lengths", []) or [])
        if not bg_lengths:
            raw = getattr(parameter, "_json_data", {}) or {}
            bg_lengths = list(raw.get("bg_seq_lengths", []) or [])

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
    except Exception as e:
        logger.debug(f"rescore-set coverage block skipped: {e}")
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
