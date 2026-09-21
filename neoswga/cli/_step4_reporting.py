"""Reporting that runs after `optimize` has produced a set.

Split out of ``cli/pipeline.py``, which was over its module size budget.
Neither function selects or scores anything: they describe a set that has
already been chosen, so they belong together and away from the step body.
"""

import logging
import os

logger = logging.getLogger(__name__)


def _report_occupancy_by_length(parameter, primers):
    """What a mixed-length delivered panel asks of one reaction.

    A panel size counts a 7-mer bound 4% of the time and an 11-mer bound 99%
    of the time as one member each. Silent on a single-length panel, which is
    every panel this project delivered before 2026-09-21.

    Best-effort: a diagnostic must not fail the step that produced the panel.
    """
    try:
        from neoswga.core.length_occupancy import log_occupancy_spread
        from neoswga.core.reaction_conditions import build_reaction_conditions

        log_occupancy_spread(
            [str(primer) for primer in primers],
            build_reaction_conditions(parameter),
            label="delivered panel",
        )
    except Exception as exc:  # pragma: no cover - diagnostic only
        logger.debug(f"Could not report occupancy by length: {exc}")


def _report_marginal_coverage(parameter, cache, primers):
    """Print cumulative foreground coverage as the delivered primers are added.

    Set size is the most consequential choice a user makes and nothing showed
    the shape of the return curve: `--auto-size` estimates a size without
    reading the candidate pool and `--show-frontier` stops at 20 primers.

    Measured on THIS set in its delivered order, at the same reach the run was
    selected and scored on, so the last row equals the reported fg_coverage. A
    prefix of this set is not the set an optimizer would return at that size,
    so each row is a lower bound.
    """
    if cache is None or not primers:
        return
    try:
        from neoswga.core.coverage import marginal_coverage_curve, resolve_coverage_reach

        reach = resolve_coverage_reach(
            getattr(parameter, "polymerase", "phi29") or "phi29",
            override=getattr(parameter, "coverage_reach", None),
        )
        curve = marginal_coverage_curve(
            cache=cache,
            primers=list(primers),
            prefixes=getattr(parameter, "fg_prefixes", []) or [],
            seq_lengths=getattr(parameter, "fg_seq_lengths", []) or [],
            extension=reach,
            circular=bool(getattr(parameter, "fg_circular", False)),
        )
    except Exception as e:  # pragma: no cover - reporting only
        logger.debug(f"Marginal coverage curve skipped: {e}")
        return

    if not curve:
        return

    # Geometric checkpoints plus the last row: a 160-primer set must not print
    # 160 lines.
    marks = {1, 2, 4, 8, 16, 32, 64, 96, 128, 160, len(curve)}
    rows = [r for r in curve if r["n"] in marks]
    logger.info("")
    logger.info("=" * 60)
    logger.info(f"Marginal foreground coverage (reach {reach} bp)")
    logger.info("=" * 60)
    logger.info(f"{'n':>5}  {'coverage':>9}  {'pp/primer':>10}")
    logger.info("-" * 28)
    prev_n, prev_cov = 0, 0.0
    for row in rows:
        span = row["n"] - prev_n
        rate = (row["coverage"] - prev_cov) * 100.0 / span if span else 0.0
        logger.info(f"{row['n']:>5}  {row['coverage']:>9.3f}  {rate:>10.2f}")
        prev_n, prev_cov = row["n"], row["coverage"]
    logger.info(
        "Measured on this set in its delivered order, so each row is a lower "
        "bound on re-optimizing at that size. A flat pp/primer column means "
        "more primers buy little coverage; it says nothing about specificity, "
        "which peaks and then decays. Use --show-frontier for that trade-off, "
        "up to 20 primers."
    )


def _report_pareto_frontier(
    args,
    parameter,
    results,
    cache,
    show_frontier,
    application,
    quick_estimate,
    min_fg_bg_ratio,
):
    """Build and report the coverage / specificity frontier.

    Extracted from ``run_step4`` unchanged apart from the parameters it now
    takes explicitly, to keep that function inside its length budget.
    """
    if show_frontier and results and cache is not None:
        try:
            import pandas as pd

            from neoswga.core.pareto_frontier import (
                generate_frontier_report,
                plot_frontier,
                summarize_frontier_for_cli,
            )
            from neoswga.core.set_size_optimizer import (
                ParetoFrontierGenerator,
                select_from_frontier,
            )

            logger.info("")
            logger.info("=" * 60)
            logger.info("Pareto Frontier Analysis")
            logger.info("=" * 60)

            # Load step2 or step3 DataFrame for primer pool
            data_dir = parameter.data_dir
            step3_file = os.path.join(data_dir, "step3_df.csv")
            step2_file = os.path.join(data_dir, "step2_df.csv")

            if os.path.exists(step3_file):
                primer_pool = pd.read_csv(step3_file)
            elif os.path.exists(step2_file):
                primer_pool = pd.read_csv(step2_file)
            else:
                raise FileNotFoundError("No primer pool CSV found")

            # Get genome lengths
            # `parameter.fg_lengths` is not a module attribute -- it is a
            # local inside `get_params`. The global is `fg_seq_lengths`, so
            # this read always fell through to its default and the frontier
            # scored every design against a 1 Mb genome and no background.
            fg_lengths = getattr(parameter, "fg_seq_lengths", []) or []
            bg_lengths = getattr(parameter, "bg_seq_lengths", []) or []
            if not fg_lengths:
                raise ValueError(
                    "No foreground genome lengths available; cannot build a "
                    "coverage frontier. Run count-kmers first."
                )
            fg_prefixes = parameter.fg_prefixes
            bg_prefixes = getattr(parameter, "bg_prefixes", [])

            # Create frontier generator. The generator uses `processivity`
            # as the coverage read-length, so pass the REALISTIC per-primer
            # reach (phi29 ~3 kb), not single-molecule processivity (70 kb),
            # which would inflate the frontier's coverage estimates.
            from neoswga.core.coverage import polymerase_extension_reach

            _frontier_reach = polymerase_extension_reach(
                getattr(parameter, "polymerase", "phi29") or "phi29",
                coverage_metric="realistic",
            )
            generator = ParetoFrontierGenerator(
                primer_pool=primer_pool,
                position_cache=cache,
                fg_prefixes=fg_prefixes,
                bg_prefixes=bg_prefixes,
                fg_seq_lengths=fg_lengths,
                bg_seq_lengths=bg_lengths,
                processivity=_frontier_reach,
            )

            # Generate frontier (quick estimation only if requested)
            frontier_result = generator.generate_frontier(
                min_size=4,
                max_size=min(20, len(primer_pool)),
                quick_only=quick_estimate,
                verbose=not args.quiet,
            )

            # Select from frontier based on application
            selected, explanation = select_from_frontier(
                frontier_result.pareto_points,
                application=application,
                min_fg_bg_ratio=min_fg_bg_ratio,
            )
            frontier_result.selected_point = selected
            frontier_result.selection_explanation = explanation

            # Display summary
            logger.info(summarize_frontier_for_cli(frontier_result, application))
            logger.info("")
            logger.info(explanation)

            # Try to save plot
            try:
                fig = plot_frontier(frontier_result, application=application)
                plot_path = os.path.join(data_dir, "pareto_frontier.png")
                fig.savefig(plot_path, dpi=150, bbox_inches="tight")
                logger.info(f"Pareto frontier plot saved to: {plot_path}")
                import matplotlib.pyplot as plt

                plt.close(fig)
            except Exception as e:
                logger.debug(f"Could not save frontier plot: {e}")

            logger.info("=" * 60)

        except Exception as e:
            logger.warning(f"Pareto frontier analysis failed: {e}")
            import traceback

            logger.debug(traceback.format_exc())
