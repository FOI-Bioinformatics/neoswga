"""Plan the number of oligos needed for coverage and specificity targets."""

import hashlib
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
    logger.info(
        "Candidate universe: %d eligible in the inventory, frontier %d. "
        "%d were not examined by this run.",
        described["universe"],
        described["frontier"] or described["universe"],
        described.get("unexamined", 0),
    )
    return source


def run_plan_pool(args):
    import pandas as pd

    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.coverage import resolve_coverage_reach
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core.pool_plan_report import write_pool_plan
    from neoswga.core.pool_planner import plan_pool
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import build_reaction_conditions
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
    conditions = build_reaction_conditions(SimpleNamespace(**params))
    reach = resolve_coverage_reach(
        params.get("polymerase", "phi29"),
        override=args.coverage_reach or params.get("coverage_reach"),
    )
    # An absent flag takes OptimizerConfig's own default rather than a second,
    # different one declared here; the two used to disagree, 100000 against 10000.
    swap_max_evaluations = (
        args.swap_max_evaluations
        if args.swap_max_evaluations is not None
        else OptimizerConfig.swap_max_evaluations
    )
    config = OptimizerConfig(
        max_dimer_bp=params.get("max_dimer_bp", 3),
        max_self_dimer_bp=params.get("max_self_dimer_bp", 4),
        min_tm=params.get("min_tm", 20),
        max_tm=params.get("max_tm", 50),
        extension_reach=reach,
        fg_circular=params.get("fg_circular", False),
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
