"""Plan the number of oligos needed for coverage and specificity targets."""

import hashlib
import inspect
import json
import logging
from pathlib import Path
from types import SimpleNamespace

logger = logging.getLogger(__name__)


def run_report_pool(args):
    """Render saved design results without genome indexes or optimization."""
    from neoswga.core.pool_plan_report import write_pool_plan

    source = Path(args.input)
    if source.is_dir():
        source = source / "pool_plan.json"
    # The shared handler answers a bare FileNotFoundError with advice about
    # genome files and data paths. This command reads neither, so say what is
    # actually missing and where to get it.
    if not source.is_file():
        raise FileNotFoundError(
            f"No saved plan at {source.resolve()}. Pass a pool_plan.json written by "
            f"'neoswga plan-pool', or the directory holding one."
        )
    plan = json.loads(source.read_text())
    if args.title:
        plan["title"] = args.title
    report = write_pool_plan(plan, args.output)
    print(f"Report: {report.resolve()}")


def load_grid_file(path, baseline):
    """Read a design grid and resolve it against the run's own chemistry.

    A grid entry names what CHANGES. Each is applied on top of `baseline`, the
    reaction this run already resolved, so a grid varying only DMSO keeps the
    user's buffer, salts and per-oligo concentration. Rebuilding from the
    overrides alone would compare designs against library defaults rather than
    against the reaction being run, and that comparison would read as a
    chemistry result.
    """
    from neoswga.core.pool_design_sweep import load_design_grid

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(
            f"No design grid at {source.resolve()}. It is a JSON object with "
            f'"lengths" and "conditions", where each condition names only the '
            f"fields it changes."
        )
    return load_design_grid(json.loads(source.read_text()), baseline)


def _open_source(data_dir, conditions, args, fallback):
    """The inventory when this directory has one, the CSV list otherwise.

    An explicit `--candidates` file always wins: the user named the pool, and
    quietly designing over a different one would be worse than useless. A
    directory written before the inventory existed falls back to the list it
    was going to use anyway, with a line saying which it took, because "which
    pool did this run search" should not have to be inferred.
    """
    from neoswga.core.candidate_source import open_candidate_source

    if args.candidates:
        return open_candidate_source(data_dir, "", [args.primer_length], candidates=fallback)
    try:
        source = open_candidate_source(
            data_dir,
            conditions.fingerprint(),
            [args.primer_length],
            frontier=len(fallback),
        )
    except ValueError as exc:
        logger.info("Designing from the candidate list: %s", exc)
        return open_candidate_source(data_dir, "", [args.primer_length], candidates=fallback)
    described = source.describe()
    # No examined count here: the frontier is drawn inside `plan_pool`, so at
    # this point nothing has been looked at and `unexamined` is the whole
    # universe. Quoting it as the shortfall read as "this run examined none of
    # it". The figure the reader wants is in `candidate_source` in the saved
    # plan, after the search.
    frontier = described["frontier"] or described["universe"]
    logger.info(
        "Candidate universe: %d eligible in the inventory; this run will search a "
        "frontier of %d, leaving %d unexamined.",
        described["universe"],
        frontier,
        max(0, described["universe"] - frontier),
    )
    return source


def _sweep_over_grid(
    args, params, base, fg, bg, config, method, grid_path, out, resolved_conditions
):
    """Design once per condition and length in the grid, then compare.

    The production caller `design_sweep` was written for and never had. It was
    built, tested, merged and unreachable, and `--design-grid` was parsed,
    documented and read by nothing: audit finding F4, and one of the six
    capabilities `tests/test_no_capability_is_unreachable.py` was created to
    catch.

    A cache and an optimizer are rebuilt per (condition, length) rather than
    shared. The index is per length, and the chemistry is what the grid varies,
    so sharing either would compare designs that were not designed under the
    conditions they are labelled with.
    """
    from neoswga.core.candidate_inventory import STAGE2_INVENTORY_NAME, CandidateInventory
    from neoswga.core.candidate_source import InventoryCandidateSource
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.pool_design_sweep import design_sweep
    from neoswga.core.pool_objective import PoolConstraints
    from neoswga.core.pool_planner import plan_pool
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import ReactionConditions

    data_dir = base / params.get("data_dir", ".")
    # The baseline is the RESOLVED reaction, not params.json: a grid varying
    # only DMSO must keep this run's buffer, salts and oligo concentration, and
    # rebuilding from the overrides alone would compare designs against library
    # defaults. Filtered to the constructor's own parameters, because params.json
    # carries genome prefixes and lengths that are not reaction fields.
    accepted = set(inspect.signature(ReactionConditions.__init__).parameters) - {"self"}
    baseline = {
        name: getattr(resolved_conditions, name)
        for name in accepted
        if hasattr(resolved_conditions, name)
    }
    lengths, conditions = load_grid_file(grid_path, baseline)
    inventory_path = data_dir / STAGE2_INVENTORY_NAME
    if not inventory_path.is_file():
        raise FileNotFoundError(
            f"A design grid needs the candidate inventory at {inventory_path}, "
            "which `neoswga filter` writes. Without it there is nothing to "
            "design each condition from."
        )

    from dataclasses import replace
    from types import SimpleNamespace

    from neoswga.core.panel_acceptance import constraints_from_parameter

    constraints = replace(
        constraints_from_parameter(SimpleNamespace(**params)) or PoolConstraints(),
        coverage_metric=args.coverage_metric,
        min_selectivity_density=(
            args.min_selectivity_density
            if args.min_selectivity_density is not None
            else params.get("min_selectivity_density")
        ),
        max_background_sites=(
            args.max_background_sites
            if args.max_background_sites is not None
            else params.get("max_background_sites")
        ),
    )

    def run_design(condition, length, sizes, coverage_targets, constraints, provider):
        candidates = provider._eligible_in_search_order()
        cache = PositionCache(fg + bg, candidates, on_missing="error")
        cache.require_record_metadata(fg + bg)
        optimizer = OptimizerFactory.create(
            method,
            cache,
            fg,
            params["fg_seq_lengths"],
            bg,
            params.get("bg_seq_lengths", []) if bg else [],
            config=config,
            conditions=condition,
            polymerase=params.get("polymerase", "phi29"),
        )
        source = InventoryCandidateSource(provider, frontier=len(candidates))
        return plan_pool(
            optimizer,
            source,
            list(sizes),
            list(coverage_targets),
            primer_length=length,
            min_selectivity_density=constraints.min_selectivity_density,
            max_background_sites=constraints.max_background_sites,
            coverage_metric=constraints.coverage_metric,
            repair=not args.no_repair,
            panel_constraints=constraints,
        )

    with CandidateInventory(inventory_path) as inventory:
        sweep = design_sweep(
            inventory,
            conditions,
            lengths,
            list(range(args.min_size, args.max_size + 1)),
            args.coverage_targets,
            constraints,
            run_design,
        )

    out.mkdir(parents=True, exist_ok=True)
    written = out / "design_sweep.json"
    written.write_text(json.dumps(sweep, indent=2, default=str))
    print(f"Designed {len(sweep['designs'])} condition/length combinations")
    for design in sweep["designs"]:
        if design.get("result") is None:
            print(f"  {design['condition']} at {design['length']}mer: {design.get('note')}")
            continue
        sizes = [r["size"] for r in design["result"]["rows"] if r.get("eligible")]
        print(
            f"  {design['condition']} at {design['length']}mer: "
            f"{design['eligible_candidates']} eligible, "
            f"{len(sizes)} qualifying panel(s)"
        )
    print(f"Sweep: {written.resolve()}")
    return sweep


def run_plan_pool(args):
    import pandas as pd

    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.coverage import resolve_coverage_reach
    from neoswga.core.design_context import design_context_from_params
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.pool_plan_report import write_pool_plan
    from neoswga.core.pool_planner import plan_pool
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.unified_optimizer import _ensure_optimizers_registered

    source = Path(args.json_file).resolve()
    params = json.loads(source.read_text())
    base = source.parent
    fg = [str((base / p).resolve()) for p in params["fg_prefixes"]]
    bg = (
        []
        if args.no_background
        else [str((base / p).resolve()) for p in params.get("bg_prefixes", [])]
    )
    if not bg and not args.no_background:
        raise ValueError("Supply background indexes or explicitly use --no-background")
    pool_path = (
        Path(args.candidates).resolve()
        if args.candidates
        else base / params.get("data_dir", ".") / "step3_df.csv"
    )
    candidates = pd.read_csv(pool_path)["primer"].astype(str).str.upper().tolist()
    candidates = list(dict.fromkeys(p for p in candidates if len(p) == args.primer_length))
    if not candidates:
        raise ValueError(
            f"No {args.primer_length}-mers in {pool_path}; run count-kmers/filter/score at this length first"
        )
    data_dir = base / params.get("data_dir", ".")
    for prefix in fg + bg:
        path = Path(f"{prefix}_{args.primer_length}mer_positions.h5")
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing position index {path}; run count-kmers and filter first"
            )
    if args.min_size < 1 or args.max_size < args.min_size:
        raise ValueError("Require 1 <= min-size <= max-size")
    out = Path(args.output)
    # Checked here as well as in write_pool_plan, so a busy output path fails
    # before the indexes are loaded rather than after the search has run.
    if out.exists() and not out.is_dir():
        raise ValueError(f"{out} is not a directory; pass a new or empty directory")
    if out.exists() and any(out.iterdir()):
        raise ValueError(
            "Output directory is not empty; use a new directory to preserve earlier designs"
        )
    cache = PositionCache(fg + bg, candidates, on_missing="error")
    # Before any panel is evaluated, and for both references. An index written
    # before record-aware geometry looks complete and silently lets coverage
    # windows cross contig boundaries; a new design must not run on one.
    cache.require_record_metadata(fg + bg)
    # One resolution, shared with `expand-primers`. Listing the fields here as
    # well is how `max_dimer_bp` came to be 3 in one place and 4 in another.
    context = design_context_from_params(params, coverage_reach_override=args.coverage_reach)
    conditions = context.conditions
    reach = context.coverage_reach
    # An absent flag takes OptimizerConfig's own default rather than a second,
    # different one declared here; the two used to disagree, 100000 against 10000.
    swap_max_evaluations = (
        args.swap_max_evaluations
        if args.swap_max_evaluations is not None
        else params.get("swap_max_evaluations", OptimizerConfig.swap_max_evaluations)
    )
    config = context.optimizer_config(
        refinement_method="swap",
        swap_max_evaluations=swap_max_evaluations,
        allow_dimer_relaxation=False,
        verbose=False,
    )
    _ensure_optimizers_registered()
    method = args.method or ("background-aware" if bg else "hybrid")
    optimizer = OptimizerFactory.create(
        method,
        cache,
        fg,
        params["fg_seq_lengths"],
        bg,
        params.get("bg_seq_lengths", []) if bg else [],
        config=config,
        conditions=conditions,
        polymerase=params.get("polymerase", "phi29"),
    )
    if args.design_grid:
        # A grid asks for several designs, so it does not produce the single
        # `pool_plan` report the rest of this function writes.
        _sweep_over_grid(
            args, params, base, fg, bg, config, method, args.design_grid, out, conditions
        )
        return

    # Read the candidates through the shared source. The inventory holds every
    # candidate that cleared hard QC, and until now nothing in production read
    # it: the design searched the CSV shortlist and the rest were stored,
    # indexed, and unable to affect any panel.
    #
    # The frontier starts at exactly that shortlist, so this changes no
    # delivered panel. What it adds is the seam expansion will use, and the
    # counts that let a reader tell a shortlist from a universe.
    candidate_source = _open_source(data_dir, conditions, args, candidates)
    plan = plan_pool(
        optimizer,
        candidate_source,
        range(args.min_size, args.max_size + 1),
        args.coverage_targets,
        primer_length=args.primer_length,
        min_selectivity_density=args.min_selectivity_density,
        max_background_sites=args.max_background_sites,
        coverage_metric=args.coverage_metric,
        repair=not args.no_repair,
        progress=lambda n: logger.info("Evaluating a pool of up to %d oligos", n),
        panel_constraints=context.constraints,
    )
    plan["inputs"] = dict(
        params=str(source),
        candidates=str(pool_path),
        params_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        candidates_sha256=hashlib.sha256(pool_path.read_bytes()).hexdigest(),
        foreground_prefixes=fg,
        background_prefixes=bg,
    )
    plan["title"] = args.title or "Oligo pool design"
    plan["design_parameters"] = params
    plan["search"] = dict(
        method=method,
        min_size=args.min_size,
        max_size=args.max_size,
        refinement_method="swap",
        swap_max_evaluations=swap_max_evaluations,
        swap_max_seconds=config.swap_max_seconds,
    )
    report = write_pool_plan(plan, out)
    for rec in plan["recommendations"]:
        count = (
            str(rec["size"]) + " oligos" if rec["size"] is not None else "no qualifying panel found"
        )
        print(f"{rec['target_coverage']:.0%} estimated coverage: {count}")
    print(f"Report: {report.resolve()}")
    print(plan["interpretation"])


def add_parsers(subparsers):
    report = subparsers.add_parser(
        "report-pool", help="Generate a pool design report from saved plan-pool results"
    )
    report.add_argument("--input", required=True, help="pool_plan.json or its directory")
    report.add_argument("--title", help="Override the report title")
    report.add_argument("-o", "--output", required=True, help="New or empty output directory")
    p = subparsers.add_parser(
        "plan-pool", help="Find small oligo pools meeting coverage and specificity targets"
    )
    p.add_argument("-j", "--json-file", required=True)
    p.add_argument("--candidates", help="Candidate CSV (default: data_dir/step3_df.csv)")
    p.add_argument("--title", help="Title for the saved report")
    p.add_argument("--primer-length", type=int, default=12)
    p.add_argument("--min-size", type=int, default=1)
    p.add_argument("--max-size", type=int, default=24)
    p.add_argument("--coverage-targets", type=float, nargs="+", default=[0.9, 0.95])
    p.add_argument("--coverage-metric", choices=["raw", "effective"], default="effective")
    p.add_argument("--coverage-reach", type=int)
    p.add_argument(
        "--no-repair",
        action="store_true",
        help=(
            "Do not attempt a bounded second pass on a panel that misses a "
            "constraint; report it as the optimizer returned it"
        ),
    )
    p.add_argument(
        "--min-selectivity-density",
        type=float,
        help="Minimum target/background binding-site density ratio",
    )
    p.add_argument(
        "--max-background-sites",
        type=int,
        help="Maximum exact background sites in a delivered panel",
    )
    p.add_argument(
        "--no-background", action="store_true", help="Explicitly omit specificity assessment"
    )
    p.add_argument("--method", choices=["dominating-set", "hybrid", "background-aware", "clique"])
    p.add_argument(
        "--swap-max-evaluations",
        type=int,
        default=None,
        help="Swap evaluations (default: the OptimizerConfig default)",
    )
    p.add_argument(
        "--design-grid",
        default=None,
        help="JSON file with 'lengths' and 'conditions'; each condition names "
        "only the fields it changes, applied over this run's resolved chemistry. "
        "Designs one pool per combination and reports the trade-off frontier.",
    )
    p.add_argument("-o", "--output", default="pool_plan")
    return p
