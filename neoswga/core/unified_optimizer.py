"""
Unified optimizer entry point for CLI and programmatic use.

This module provides a clean interface for running primer set optimization
using any registered optimizer. It replaces the complex conditional logic
in cli_unified.py with a simple factory-based approach.

Usage:
    # From CLI
    result = run_optimization('greedy', candidates, config)

    # Programmatic with custom config
    result = run_optimization(
        method='dominating-set',
        candidates=candidates,
        fg_prefixes=['data/target'],
        fg_seq_lengths=[1000000],
        target_size=10,
    )
"""

import json
import logging
import os
import threading
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from . import parameter
from .base_optimizer import OptimizationResult, OptimizationStatus, OptimizerConfig
from .candidate_source import describe_reach, order_candidates_by_background
from .design_result import panel_validation_is_ok
from .dimer import dimer_validation_issue, worst_heterodimer
from .ensemble_comparison import _ensemble_error_row, _select_ensemble_winner
from .exceptions import DesignError, ModelEvaluationError, ReferenceDataError
from .optimization_reporting import log_optimization_outcome
from .optimizer_factory import OptimizerFactory, OptimizerRegistry
from .panel_acceptance import (
    apply_configured_limits,
    check_per_target_coverage,
    constraints_from_parameter,
    per_target_floor,
    report_per_target,
)
from .panel_regime import assess_from_parameter, log_regime
from .position_cache import PositionCache, StreamingPositionCache
from .progress import progress_context
from .search_control import SearchBudgetExhausted, collect_alternative_sets  # noqa: F401
from .step4_output import _write_validation_report, save_results
from .validation_record import limit_violation_issue, panel_assessment

logger = logging.getLogger(__name__)

# Thread-safe registration tracking
_optimizers_registered = False
_registration_lock = threading.Lock()

# Last raw OptimizationResult from optimize_step4 / run_optimization.
# Exposed so the CLI can surface pareto_front / validation without changing
# the legacy (primer_sets, scores, cache) return contract. Access via
# unified_optimizer._LAST_RESULT; do not rely on this for library use.
_LAST_RESULT: Optional["OptimizationResult"] = None

# Alternative primer sets from the most recent run, best first. Populated
# when `max_sets` asks for more than one; always at least the primary set.
_LAST_PRIMER_SETS: List[Tuple[str, ...]] = []

#: How much of the available candidate pool the last run could reach.
#: Stashed rather than returned because `optimize_step4` writes the summary
#: and `run_optimization` is where the frontier is in scope -- the same
#: reason `_LAST_PRIMER_SETS` exists. None when the caller passed a plain
#: list, which is absence rather than a reach of zero.
_LAST_CANDIDATE_REACH: Optional[Dict[str, Any]] = None


@dataclass
class OptimizationConfig:
    """
    Unified configuration for optimization runs.

    Combines all parameters needed for optimization in one place.
    """

    # Optimizer selection
    method: str = "hybrid"

    # Target parameters
    target_set_size: int = 6
    max_iterations: int = 100

    # File paths
    data_dir: str = "."
    step3_file: str = "step3_df.csv"
    output_file: str = "step4_improved_df.csv"

    # Genome info (loaded from pipeline if not provided)
    fg_prefixes: Optional[List[str]] = None
    fg_seq_lengths: Optional[List[int]] = None
    bg_prefixes: Optional[List[str]] = None
    bg_seq_lengths: Optional[List[int]] = None

    # Performance options
    use_position_cache: bool = True
    use_background_filter: bool = True

    # Output options
    verbose: bool = True
    quiet: bool = False

    # Advanced options
    allow_dimer_relaxation: bool = False
    refinement_method: str = "network"
    swap_max_evaluations: int = 10000
    stage1_objective_width: Optional[int] = None
    swap_max_seconds: float = 10.0
    uniformity_weight: float = 0.0
    minimize_primers: bool = False
    target_coverage: float = 0.70


def list_available_optimizers() -> Dict[str, str]:
    """
    List all registered optimizers with descriptions.

    Returns:
        Dict mapping optimizer names to descriptions
    """
    # Ensure all optimizers are registered
    _ensure_optimizers_registered()

    return OptimizerFactory.list_optimizers()


def _ensure_optimizers_registered():
    """
    Ensure all optimizer modules are imported and registered.

    Thread-safe and idempotent - only performs registration once.
    """
    global _optimizers_registered

    # Fast path: already registered
    if _optimizers_registered:
        return

    with _registration_lock:
        # Double-checked locking
        if _optimizers_registered:
            return

        # Import optimizer modules to trigger factory registration
        try:
            from . import background_aware_optimizer  # Clinical/background-aware
            from . import clique_optimizer  # Structurally dimer-free sets
            from . import dominating_set_adapter  # Graph-based set cover
            from . import hybrid_optimizer  # Two-stage hybrid
            from . import network_optimizer  # Network connectivity
        except ImportError as e:
            logger.warning(f"Some optimizer modules not available: {e}")

        _optimizers_registered = True
        logger.debug("Optimizer registration complete")


def _reseed(seed) -> None:
    """Re-seed all RNGs so ensemble methods are mutually reproducible and
    order-independent. No-op when seed is None."""
    if seed is None:
        return
    import random

    import numpy as np

    random.seed(seed)
    np.random.seed(seed)
    try:
        from .rf_preprocessing import set_kmer_sampling_seed
    except ImportError:
        # The amplification model was retired from the default path, so its
        # sampling RNG may not be importable at all. Nothing consumed it, so
        # there is nothing to make reproducible.
        return
    try:
        set_kmer_sampling_seed(seed)
    except Exception as exc:
        # A seed that was asked for and not applied is a broken promise, and
        # the caller logs "set for reproducibility" immediately afterwards.
        raise ModelEvaluationError("random seed", seed, str(exc)) from exc


def _resolve_target_coverage(kwargs, default=0.70):
    """The coverage floor `--minimize-primers` trims down to.

    Pops the key, so it is not forwarded twice.

    This was `float(kwargs.pop("target_coverage", 0.70) or 0.70)`, and `or`
    treats an explicit 0 as absent -- so `--target-coverage 0`, a legitimate
    request meaning "trim as far as you can", silently became 0.70. Only a
    MISSING value may fall back.
    """
    value = kwargs.pop("target_coverage", None)
    return float(default) if value is None else float(value)


def _run_ensemble(
    methods: List[str],
    cache,
    candidates: List[str],
    fg_prefixes: List[str],
    fg_seq_lengths: List[int],
    bg_prefixes: Optional[List[str]],
    bg_seq_lengths: Optional[List[int]],
    target_size: int,
    config,
    conditions,
    application: str = "balanced",
    seed=None,
    verbose: bool = True,
    combine: str = "best",
    **kwargs,
) -> OptimizationResult:
    """Compare proposals under one evaluator and attach per-method metrics.

    Method failures are reported and skipped. Normalized application scores
    remain descriptive; an available shared objective decides the winner.
    """
    from dataclasses import replace as _dc_replace

    rows: List[Dict[str, Any]] = []
    results: Dict[str, OptimizationResult] = {}
    from .search_control import SearchBudget

    shared_objective = None
    search_budget = kwargs.pop("search_budget", None) or SearchBudget.from_config(config)
    candidate_source = kwargs.pop("candidate_source", None)

    for m in methods:
        _reseed(seed)  # order-independent reproducibility
        try:
            optimizer = OptimizerFactory.create(
                name=m,
                position_cache=cache,
                fg_prefixes=fg_prefixes,
                fg_seq_lengths=fg_seq_lengths,
                bg_prefixes=bg_prefixes,
                bg_seq_lengths=bg_seq_lengths,
                config=config,
                conditions=conditions,
                **kwargs,
            )
            from .optimization_service import OptimizationRequest, run_panel_search
            from .panel_refinement import objective_for_optimizer

            if shared_objective is None:
                shared_objective = objective_for_optimizer(
                    optimizer, constraints_from_parameter(parameter)
                )
            else:
                from .swap_refinement import attach_search_config

                attach_search_config(optimizer, "pool_objective", shared_objective)
            with progress_context(f"  ensemble: {optimizer.name}", disable=not verbose):
                res = run_panel_search(
                    OptimizationRequest(
                        optimizer,
                        tuple(candidates),
                        target_size,
                        budget=search_budget,
                        candidate_source=candidate_source,
                    )
                )
        except DesignError:
            # A failed calculation, a missing reference answer or an
            # unsupported model is not a property of this ensemble member: the
            # next member would compute the same quantity from the same data.
            # Recording it as one method's error and letting another method
            # win returns a panel built on the same unmeasured ground, with
            # nothing in the comparison table to say so.
            raise
        except Exception as e:
            # A member that fails for its own reasons -- an algorithm that
            # cannot run on this pool -- stays a visible row in
            # `ensemble_comparison` rather than failing the command. The
            # ensemble was asked for several methods precisely because any one
            # of them may not apply.
            logger.warning(f"  ensemble method '{m}' failed: {e}")
            rows.append(_ensemble_error_row(m))
            continue

        norm = res.metrics.normalized_score(application=application)
        # Only a method that actually produced a set is a candidate to win.
        #
        # Hybrid, network and background-aware all return NO_CONVERGENCE with an
        # empty primer tuple rather than raising, and an empty
        # `PrimerSetMetrics` does NOT score zero -- a (0.0, 0.0) Tm range reads
        # as perfectly tight and collects the whole Tm term, so it scores 0.10
        # under balanced weights and 0.05 under clinical. Such a result was
        # therefore rankable, and when every method came back empty the ensemble
        # returned one of them as the "winner", after which `optimize_step4`
        # logged "Optimization failed: " with an empty message. The clear
        # all-methods-failed error below already existed and was unreachable.
        if res.primers:
            results[m] = res
        rows.append(
            {
                "method": m,
                "normalized_score": round(float(norm), 4),
                "score": float(res.score),
                "n_primers": res.num_primers,
                "fg_coverage": round(float(res.metrics.fg_coverage), 4),
                "bg_coverage": round(float(res.metrics.bg_coverage), 4),
                "status": res.status.value,
                "selected": False,
            }
        )

    if not results:
        return OptimizationResult.failure("ensemble", f"All ensemble methods failed: {methods}")

    best_method = _select_ensemble_winner(results, application, shared_objective)
    for row in rows:
        row["selected"] = row["method"] == best_method and row["status"] != "error"

    if verbose:
        logger.info(
            f"Ensemble winner: '{best_method}' "
            f"(normalized_score={results[best_method].metrics.normalized_score(application=application):.4f}, "
            f"application='{application}')"
        )

    winner = results[best_method]

    # Union-combine: pool the (diverse) primers every method selected and
    # re-optimize the winning method on that small pool. Other methods often
    # found primers covering regions the winner missed; re-optimizing over the
    # union can beat any single method. Guarded to NEVER worsen the winner.
    if combine == "union" and len(results) > 1:
        union_pool = list(dict.fromkeys(p for r in results.values() for p in r.primers))
        if len(union_pool) > target_size:
            _reseed(seed)
            try:
                combiner = OptimizerFactory.create(
                    name=best_method,
                    position_cache=cache,
                    fg_prefixes=fg_prefixes,
                    fg_seq_lengths=fg_seq_lengths,
                    bg_prefixes=bg_prefixes,
                    bg_seq_lengths=bg_seq_lengths,
                    config=config,
                    conditions=conditions,
                    **kwargs,
                )
                if shared_objective is not None:
                    from .swap_refinement import attach_search_config

                    attach_search_config(combiner, "pool_objective", shared_objective)
                else:
                    objective_for_optimizer(combiner, constraints_from_parameter(parameter))
                combined = run_panel_search(
                    OptimizationRequest(
                        combiner, tuple(union_pool), target_size, budget=search_budget
                    )
                )
                c_norm = combined.metrics.normalized_score(application=application)
                w_norm = winner.metrics.normalized_score(application=application)
                from .optimization_service import panel_rank

                improved = (
                    panel_rank(shared_objective, combined.primers)
                    >= panel_rank(shared_objective, winner.primers)
                    if shared_objective is not None and combined.primers
                    else c_norm >= w_norm
                )
                from .optimization_service import panel_violations

                if (
                    combined.primers
                    and improved
                    and not panel_violations(combiner, combined.primers)
                ):
                    winner = combined
                    for row in rows:
                        row["selected"] = False
                    rows.append(
                        {
                            "method": f"union({best_method})",
                            "normalized_score": round(float(c_norm), 4),
                            "score": float(combined.score),
                            "n_primers": combined.num_primers,
                            "fg_coverage": round(float(combined.metrics.fg_coverage), 4),
                            "bg_coverage": round(float(combined.metrics.bg_coverage), 4),
                            "status": combined.status.value,
                            "selected": True,
                        }
                    )
                    if verbose:
                        logger.info(
                            f"Ensemble union-combine improved the set "
                            f"(normalized_score {w_norm:.4f} -> {c_norm:.4f})"
                        )
            except Exception as e:
                logger.warning(f"  ensemble union-combine skipped: {e}")

    return _dc_replace(winner, ensemble_comparison=tuple(rows))


def _vet_the_index(cache, fg_prefixes, bg_prefixes, candidate_count, pipeline_path):
    """Refuse an index a new design cannot honestly be scored against.

    Two checks with two different subjects. The first asks whether the index is
    of a shape a current design can use: record geometry, so a coverage window
    stops at a contig edge, and an on-disk format this version understands. The
    second asks whether it covers the pool, rather than letting selection run
    over whichever part happens to be present.

    Reference IDENTITY is deliberately NOT checked here. It is a relation
    between a prefix and a genome, and only the resolved request names both;
    `cli/pipeline.run_step4` does it there. Reading `parameter.fg_genomes` at
    this point pairs the prefixes this call was GIVEN with whatever genomes the
    module currently holds, and under `pytest -n 8` that paired a test's own
    prefix with another test's FASTA. A mutable global is not a manifest.

    Both are skipped for a caller that supplied its own pool, which is the
    library path: it has not asked step 4 to read step 3 and is not subject to
    step 4's prerequisites.
    """
    from .pipeline import validate_index_covers_candidates

    if pipeline_path:
        cache.require_record_metadata(list(fg_prefixes) + list(bg_prefixes or []))
    return validate_index_covers_candidates(
        cache, fg_prefixes, candidate_count, refuse=pipeline_path
    )


def _make_position_cache(prefixes, primers, use_cache=True):
    """The position cache a run should use.

    `--no-position-cache` sets `use_cache=False`. It used to change nothing:
    `optimize_step4` accepted the value and never referenced it, while the CLI
    logged "Position cache: False" beforehand, which reads as confirmation the
    request was honoured.

    False now selects `StreamingPositionCache`, the memory-mapped
    implementation, which is what "disable the in-memory cache (slower)" should
    mean on a background too large to hold.
    """
    if use_cache:
        return PositionCache(prefixes, primers)

    logger.info("Position cache disabled; using the memory-mapped streaming cache")
    return StreamingPositionCache(prefixes, primers)


def _build_optimizer_config(
    target_size, verbose, extension_reach, fg_circular, kwargs, method=None
):
    """Populate every configurable field, not just the ones a caller passes.

    This used to set five fields and leave the rest at their dataclass
    defaults, so `max_dimer_bp`, `max_self_dimer_bp`, `min_tm` and `max_tm`
    never arrived however params.json configured them -- every optimizer read
    the default instead.

    It mattered most on the clique optimizer, where the effect was a broken
    guarantee rather than a shifted score. A params.json asking for
    max_dimer_bp=3 got the default 4, which is LOOSER, so the compatibility
    graph gained edges and the returned "dimer-free" clique could contain a
    pair that dimerises at 3. The opposite direction would have been harmless
    -- a clique in a sparser graph is still a clique in a denser one -- which
    is why this only ever failed one way, and why unit tests that build the
    config themselves never saw it.

    Precedence: explicit kwarg, then the `parameter` global, then the default.

    When `method` names an optimizer with its own `OptimizerConfig` subclass,
    that subclass is built instead of the base, and its extra fields are
    resolved through the same precedence. Naming each one here by hand is what
    produced this defect class in the first place, so they are enumerated from
    the dataclass rather than listed.
    """

    def pick(name, default):
        value = kwargs.get(name)
        if value is None:
            value = getattr(parameter, name, None)
        return default if value is None else value

    config_class = OptimizerConfig
    if method:
        from .optimizer_factory import OptimizerRegistry

        config_class = OptimizerRegistry.config_class(method)

    # Fields the subclass adds on top of the base. Resolved generically so a
    # newly added field is reachable from params.json the day it is declared.
    extra = {}
    if config_class is not OptimizerConfig:
        import dataclasses

        base_names = {f.name for f in dataclasses.fields(OptimizerConfig)}
        for field in dataclasses.fields(config_class):
            if field.name in base_names:
                continue
            default = (
                field.default
                if field.default is not dataclasses.MISSING
                else field.default_factory()  # type: ignore[misc]
            )
            extra[field.name] = pick(field.name, default)

    from .search_control import SEARCH_CONTROL_DEFAULTS

    return config_class(
        **{key: pick(key, default) for key, default in SEARCH_CONTROL_DEFAULTS.items()},
        **extra,
        target_set_size=target_size,
        # params.json calls this `iterations` and `get_params` assigns that
        # name; `OptimizerConfig` calls it `max_iterations`. A bare
        # `pick("max_iterations", ...)` looks for a global that does not exist
        # and always took the default, so the key was inert whichever way it was
        # set. Both names are tried, the config's first.
        max_iterations=pick("max_iterations", None) or pick("iterations", 100),
        verbose=verbose,
        extension_reach=extension_reach,
        fg_circular=fg_circular,
        max_dimer_bp=pick("max_dimer_bp", 4),
        max_dimer_dg=pick("max_dimer_dg", None),
        allow_dimer_relaxation=pick("allow_dimer_relaxation", False),
        refinement_method=pick("refinement_method", "network"),
        swap_max_evaluations=pick("swap_max_evaluations", 10000),
        stage1_objective_width=pick("stage1_objective_width", None),
        swap_max_seconds=pick("swap_max_seconds", 10.0),
        max_self_dimer_bp=pick("max_self_dimer_bp", 5),
        # `parameter.max_mismatches` is assigned by _apply_params_only_keys and
        # `base_optimizer._weighted_loads` reads `self.config.max_mismatches`,
        # but nothing connected the two, so the modelled site load ran at 1
        # mismatch whatever params.json asked for. Tenth instance of the
        # documented-accepted-and-read-by-nothing class, and a missing line.
        max_mismatches=pick("max_mismatches", 1),
        min_tm=pick("min_tm", 20.0),
        max_tm=pick("max_tm", 50.0),
    )


_SELECTION_WEIGHT_KEYS = ("tm_weight", "uniformity_weight", "dimer_penalty")


def _resolve_selection_weights(kwargs: dict, application, verbose: bool) -> None:
    """Fill the selection weights from the `--application` profile in place.

    These three knobs shape the greedy selection inside NetworkOptimizer, so
    `--application` changes which primers are returned rather than only how
    they are displayed.

    None means "the caller expressed no preference"; 0.0 is a preference and
    is kept. The two cannot be told apart once a caller has substituted a
    concrete default, which is how the profile's `uniformity_weight` came to
    be unreachable: `optimize_step4` and the CLI forwarded 0.0 on every call,
    so the key was always present and the profile could never be installed.
    The sentinel therefore has to survive from argparse down to here.
    """
    for key in _SELECTION_WEIGHT_KEYS:
        if key in kwargs and kwargs[key] is None:
            del kwargs[key]

    if not application:
        return

    try:
        from .base_optimizer import OPTIMIZER_APPLICATION_WEIGHTS

        profile = OPTIMIZER_APPLICATION_WEIGHTS.get(
            str(application).lower(), OPTIMIZER_APPLICATION_WEIGHTS["balanced"]
        )
        for key, value in profile.items():
            kwargs.setdefault(key, value)
        if verbose:
            logger.info(
                f"application='{application}' → selection weights applied: "
                f"tm={kwargs['tm_weight']:.2f}, "
                f"uniformity={kwargs['uniformity_weight']:.2f}, "
                f"dimer_penalty={kwargs['dimer_penalty']:.2f}"
            )
    except KeyError as exc:
        # `--application clinical` that quietly has no effect optimizes under
        # weights the user did not ask for and reports normally. The profile
        # table is a package constant, so a miss here is a programming error
        # rather than a user one, and it must not be absorbed.
        raise ModelEvaluationError("application weights", application, str(exc)) from exc


def _collect_forbidden_primers(candidates, verbose: bool) -> list:
    """Candidates that exceed the blacklist frequency ceiling.

    Phase 15C library-level guard. The validator used to receive
    forbidden_primers=None, so library callers who bypassed the swap-primer CLI
    could silently inject blacklist candidates. Cross-checks the pool against
    params.bl_prefixes with the same helper the filter step uses, and hands the
    hits to validate() so any that reached the selected set are flagged.
    """
    bl_prefixes_list = list(getattr(parameter, "bl_prefixes", []) or [])
    bl_lengths_list = list(getattr(parameter, "bl_seq_lengths", []) or [])
    if not (bl_prefixes_list and bl_lengths_list and candidates):
        return []

    try:
        from .pipeline import _filter_blacklist_penalty

        max_bl = getattr(parameter, "max_bl_freq", 0.0) or 0.0
        _mask, freqs = _filter_blacklist_penalty(
            list(candidates), bl_prefixes_list, bl_lengths_list, max_bl_freq=max_bl
        )
        forbidden = [p for p, f in zip(candidates, freqs, strict=True) if f > max_bl]
        if forbidden and verbose:
            logger.warning(
                f"Library blacklist guard: {len(forbidden)} candidate(s) exceed "
                f"max_bl_freq={max_bl}; validator will flag any that reached the set."
            )
        return forbidden
    except Exception as exc:
        # An empty list here is indistinguishable from "no candidate is
        # blacklisted", so the validator then receives `forbidden_primers=None`
        # and a blacklisted primer in the delivered set goes unflagged. The
        # blacklist is configured, so failing to apply it is a design failure.
        raise ReferenceDataError(
            "blacklist k-mer counts",
            str(exc),
            "Re-run `neoswga count-kmers` for the blacklist genomes, or unset bl_prefixes.",
        ) from exc


def _dispatch_optimizer(
    method,
    cache,
    candidates,
    fg_prefixes,
    fg_seq_lengths,
    bg_prefixes,
    bg_seq_lengths,
    target_size,
    config,
    conditions,
    seed,
    verbose,
    _application,
    kwargs,
):
    """Run `method` and return `(result, optimizer)`.

    The optimizer comes back beside the result because the post-processing in
    `run_optimization` measures coverage with it.
    """
    optimizer = None
    search_budget = kwargs.pop("search_budget", None)
    candidate_source = kwargs.pop("candidate_source", None)

    if method in ("ensemble", "auto", "all"):
        ensemble_methods = kwargs.pop("ensemble_methods", None) or [
            "hybrid",
            "dominating-set",
            "network",
            "background-aware",
        ]
        # `seed` is captured separately and passed explicitly; remove it from
        # kwargs so it is not forwarded twice into _run_ensemble.
        kwargs.pop("seed", None)
        ensemble_combine = kwargs.pop("ensemble_combine", "best")
        result = _run_ensemble(
            methods=ensemble_methods,
            cache=cache,
            candidates=candidates,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
            target_size=target_size,
            config=config,
            conditions=conditions,
            application=(_application or "balanced"),
            seed=seed,
            verbose=verbose,
            combine=ensemble_combine,
            search_budget=search_budget,
            candidate_source=candidate_source,
            **kwargs,
        )

        # The post-processing below measures coverage with
        # `optimizer.compute_metrics`, which is `BaseOptimizer`'s and does not
        # depend on which algorithm selected the set. Build one for the winning
        # method so `--minimize-primers` works on an ensemble run rather than
        # silently doing nothing: handing the minimiser None would put an
        # AttributeError inside its `try` and return the set untrimmed, which
        # is the same defect this class of fix exists to remove.
        try:
            optimizer = OptimizerFactory.create(
                name=getattr(result, "optimizer_name", None) or ensemble_methods[0],
                position_cache=cache,
                fg_prefixes=fg_prefixes,
                fg_seq_lengths=fg_seq_lengths,
                bg_prefixes=bg_prefixes,
                bg_seq_lengths=bg_seq_lengths,
                config=config,
                conditions=conditions,
                **kwargs,
            )
        except Exception as exc:
            # Without this optimizer `_hold_to_configured_limits` returns the
            # panel unrepaired, so every configured panel limit is silently
            # not enforced and the run reports a clean result. A limit the
            # user set is a hard constraint, not a preference.
            raise ModelEvaluationError(
                "ensemble winner evaluator",
                getattr(result, "optimizer_name", None) or ensemble_methods[0],
                str(exc),
            ) from exc
    else:
        # Create optimizer via factory
        try:
            optimizer = OptimizerFactory.create(
                name=method,
                position_cache=cache,
                fg_prefixes=fg_prefixes,
                fg_seq_lengths=fg_seq_lengths,
                bg_prefixes=bg_prefixes,
                bg_seq_lengths=bg_seq_lengths,
                config=config,
                conditions=conditions,
                **kwargs,
            )
        except Exception as e:
            logger.error(f"Failed to create optimizer '{method}': {e}")
            return OptimizationResult.failure(method, str(e)), None

        from .optimization_service import OptimizationRequest, run_panel_search
        from .panel_refinement import objective_for_optimizer

        objective_for_optimizer(optimizer, constraints_from_parameter(parameter))
        # Run optimization
        with progress_context(f"Running {optimizer.name} optimizer", disable=not verbose):
            result = run_panel_search(
                OptimizationRequest(
                    optimizer,
                    tuple(candidates),
                    target_size,
                    budget=search_budget,
                    candidate_source=candidate_source,
                )
            )

    if optimizer is not None:
        from .panel_refinement import objective_for_optimizer

        objective_for_optimizer(optimizer, constraints_from_parameter(parameter))
    return result, optimizer


def _minimize_primer_count(result, optimizer, target_coverage: float, verbose: bool):
    """Remove redundant primers while preserving the shared panel contract.

    Each deletion must meet the requested coverage under the design's selected
    metric and all configured panel limits. Condition-free library evaluators
    retain geometric coverage; configured designs use effective coverage.
    """
    from .optimization_service import reduce_result

    return reduce_result(result, optimizer, target_coverage)


def _pool_for_this_run(fg_prefixes, conditions):
    """The candidates `optimize` may select, when the caller named none.

    Through the shared source rather than straight out of `step3_df.csv`. The
    inventory holds every candidate that cleared hard QC and the CSV holds the
    `max_primer` shortlist, so reading the CSV here made the rest unreachable
    however large the inventory was: audit finding F1 named this call site,
    `expand-primers` and `plan-pool`, and this is the last of the three.

    The frontier opens at the shortlist's own size, so this changes no delivered
    panel. What it adds is the counts and the seam increment 5's refill reaches
    past.

    The step-4 prerequisites are checked first, because without them the
    position cache falls back to `on_missing="warn"`, which calls the coverage
    number meaningless and then lets the run report one anyway.
    """
    from .candidate_source import open_design_source
    from .pipeline import StepPrerequisiteError, validate_step4_prerequisites

    validation = validate_step4_prerequisites(parameter.data_dir, list(fg_prefixes or []))
    if not validation.valid:
        raise StepPrerequisiteError(4, validation)

    step3_path = os.path.join(parameter.data_dir, "step3_df.csv")
    shortlist = pd.read_csv(step3_path)["primer"].astype(str).tolist()
    source = open_design_source(
        parameter.data_dir,
        conditions.fingerprint() if hasattr(conditions, "fingerprint") else "",
        sorted({len(primer) for primer in shortlist}),
        fallback=shortlist,
    )
    from .candidate_source import CandidateFrontier

    return CandidateFrontier(source)


def _report_candidate_reach(candidates, verbose: bool) -> Optional[Dict[str, Any]]:
    """What the search could have examined, against what it did.

    A run that qualifies on its opening frontier never widens, so this is
    routinely a small fraction of the inventory and nothing said so. Both
    numbers were known when the pool was opened.

    Returns None for a plain candidate list, where the question has no answer.
    """
    reach = describe_reach(candidates)
    if verbose and reach and not reach["complete"]:
        logger.info(
            f"Searched {reach['examined']:,} of {reach['universe']:,} available "
            f"candidates ({reach['fraction']:.1%}). The search stops at the "
            f"first qualifying frontier; looking further was measured to cost "
            f"specificity -- see docs/validation/"
            f"looking_further_costs_specificity_2026-09-22.md"
        )
    return reach


def run_optimization(
    method: str = "hybrid",
    candidates: Optional[List[str]] = None,
    fg_prefixes: Optional[List[str]] = None,
    fg_seq_lengths: Optional[List[int]] = None,
    bg_prefixes: Optional[List[str]] = None,
    bg_seq_lengths: Optional[List[int]] = None,
    target_size: int = 6,
    verbose: bool = True,
    design_request=None,
    **kwargs,
) -> OptimizationResult:
    """
    Run primer set optimization using specified method.

    The main entry point for optimization: loads candidates if they were not
    provided, builds the position cache, instantiates the optimizer via the
    factory, runs it, and returns a typed result.

    Args:
        method: Optimizer method name ('greedy', 'dominating-set', etc.)
        candidates: List of candidate primers (loaded from step3 if None)
        fg_prefixes: Foreground genome HDF5 prefixes
        fg_seq_lengths: Foreground genome lengths
        bg_prefixes: Background genome HDF5 prefixes (optional)
        bg_seq_lengths: Background genome lengths (optional)
        target_size: Desired primer set size
        verbose: Print progress information
        design_request: The resolved `DesignRequest`, or None when a caller
            drives the optimizer programmatically without one. It carries
            `request_hash` and is what `panel_evaluation.evaluate_panel` takes.
            `tests/test_the_request_reaches_the_optimizer.py` says why it was
            resolved before the search and did not arrive here.
        **kwargs: Additional optimizer-specific parameters

    Returns:
        OptimizationResult with selected primers and metrics

    Raises:
        StepPrerequisiteError: only when `candidates` is None, so the pool is
            read from step3_df.csv here. Raised when that file is missing or
            empty, when the position files are absent, or when the index covers
            only part of the pool. A caller supplying its own candidates is
            never blocked and gets a failure OptimizationResult instead.

    Example:
        result = run_optimization(
            method='dominating-set',
            candidates=['ATCGATCG', 'GCTAGCTA', ...],
            fg_prefixes=['data/target'],
            fg_seq_lengths=[1000000],
            target_size=10,
        )

        if result.is_success:
            print(f"Selected: {result.primers}")
            print(f"Coverage: {result.metrics.fg_coverage:.1%} (measured)")
    """
    _ensure_optimizers_registered()

    # Reset the last-result stash at entry so a second call that fails
    # before the bottom-of-function assignment cannot leak the previous
    # successful result to CLI consumers (--show-frontier, audit tooling).
    # See Phase 16 critical gap #1.
    global _LAST_RESULT
    _LAST_RESULT = None

    # Load parameters from pipeline if needed
    if fg_prefixes is None or fg_seq_lengths is None:
        from . import pipeline as core_pipeline

        core_pipeline._initialize()
        fg_prefixes = core_pipeline.fg_prefixes
        fg_seq_lengths = core_pipeline.fg_seq_lengths
        bg_prefixes = bg_prefixes or getattr(core_pipeline, "bg_prefixes", [])
        bg_seq_lengths = bg_seq_lengths or getattr(core_pipeline, "bg_seq_lengths", [])

    # Handle no-background (host-free) mode
    if kwargs.get("no_background", False):
        bg_prefixes = []
        bg_seq_lengths = []
        if verbose:
            logger.info("  Host-free mode: no background genome data used")

    # Invalid configured chemistry must stop the design before candidate loading.
    conditions = kwargs.pop("conditions", None)
    if conditions is None:
        from .reaction_conditions import build_reaction_conditions

        conditions = build_reaction_conditions()

    # Load candidates from step3, after checking that the index step 4 will
    # score against actually exists. Without this the position cache falls back
    # to on_missing="warn", which says the coverage number is "meaningless, not
    # low" and then lets the run select a set and report it.
    _pool_read_from_step3 = candidates is None
    if candidates is None:
        candidates = _pool_for_this_run(fg_prefixes, conditions)
    _candidate_source = getattr(candidates, "source", None)

    # Empty candidate pool guard: downstream optimizers behave inconsistently
    # on an empty pool (some raise, some return an empty result without saying
    # why). A pipeline run cannot reach this branch -- with candidates=None the
    # step-4 validator above rejects a missing or empty step3_df.csv first, and
    # names the file -- so the only caller here is one that passed [] itself.
    if not candidates:
        msg = (
            "No candidate primers available for optimization: the caller "
            "passed an empty candidate list. Pass a non-empty list of primer "
            "sequences, or pass candidates=None to load the pool from "
            "step3_df.csv in data_dir. Status: "
            f"data_dir={getattr(parameter, 'data_dir', '?')}, "
            f"method={method}."
        )
        logger.error(msg)
        return OptimizationResult.failure(method, msg)

    # Set global random seed if provided (ensures reproducibility across all
    # optimizer components, not just those that accept a seed parameter)
    seed = kwargs.get("seed")
    if seed is not None:
        import random

        import numpy as np

        random.seed(seed)
        np.random.seed(seed)
        # Also seed the RF k-mer sampling RNG so any re-scoring during
        # optimization is reproducible (k-mer sampling is on by default).
        # This used to log "set for reproducibility" whether or not the seed
        # had actually been applied. `_reseed` now raises instead of passing,
        # so the line below is only reached when the claim is true.
        _reseed(seed)
        if verbose:
            logger.info(f"Random seed set to {seed} for reproducibility")

    if verbose:
        logger.info(f"Running {method} optimization")
        logger.info(f"  Candidates: {len(candidates)}")
        logger.info(f"  Target size: {target_size}")

    # Create position cache
    with progress_context("Loading position data", disable=not verbose):
        cache = _make_position_cache(
            fg_prefixes + (bg_prefixes or []), candidates, kwargs.get("use_cache", True)
        )

    # Refuse an index a new design cannot be scored against: missing record
    # geometry, an older on-disk format, or one built from another reference.
    # Only `plan-pool` used to check, so the two commands most people run
    # scored against whatever the directory happened to hold. Skipped when the
    # caller supplied its own pool AND no params file named the genomes, which
    # is the library path; `_pool_read_from_step3` marks the pipeline path.
    _unindexed = _vet_the_index(
        cache, fg_prefixes, bg_prefixes, len(candidates), _pool_read_from_step3
    )

    # Orders by fg/bg ratio; deletes nothing. See order_candidates_by_background
    # for why it used to delete and why `bg_max_removal` went with that.
    if bg_prefixes and kwargs.get("bg_prefilter", True):
        candidates, _deprioritised = order_candidates_by_background(
            cache,
            candidates,
            fg_prefixes,
            bg_prefixes,
            # `or` not `.get` default: the CLI forwards None for an unset flag.
            min_ratio=kwargs.get("bg_min_ratio") or 1.0,
            verbose=verbose,
        )

    from .coverage import resolve_extension_reach

    extension_reach = resolve_extension_reach(kwargs, verbose)

    # Build optimizer config
    fg_circular = kwargs.get("fg_circular")
    if fg_circular is None:
        fg_circular = bool(getattr(parameter, "fg_circular", False))

    config = _build_optimizer_config(
        target_size=target_size,
        verbose=verbose,
        extension_reach=extension_reach,
        fg_circular=fg_circular,
        kwargs=kwargs,
        # So an optimizer with its own config subclass gets its own fields
        # resolved from params.json. The ensemble path passes its per-method
        # name too, via the same helper.
        method=method,
    )

    # Popped rather than read off `config`: these arrive as kwargs, and the
    # `OptimizerConfig` built here has no field for them -- `minimize_primers`
    # lives on the separate `OptimizationConfig` used by
    # run_optimization_from_config. Reading them from the wrong object is how
    # they came to be set, forwarded, and never acted on. Popping also keeps
    # them out of the optimizer constructors, which have no use for them.
    _minimize_primers = bool(kwargs.pop("minimize_primers", False))
    _target_coverage = _resolve_target_coverage(kwargs)

    _application = kwargs.pop("application", None)
    _resolve_selection_weights(kwargs, _application, verbose)

    # Ensemble path: run several methods on the SAME position cache and keep
    # the best by normalized_score. Falls through to the shared
    # post-processing below so the winner gets per-target coverage, saturation
    # checks, and the validation report like any single-method run.
    # Bound before the branch. The post-processing below reads `optimizer`, and
    # only the single-method branch used to assign it, so `--minimize-primers`
    # with `--optimization-method ensemble` raised UnboundLocalError while the
    # call arguments were being evaluated -- before the `try` inside
    # `_minimize_primer_count` could see it. The step then exited 1 having run
    # every ensemble method to completion and thrown all of it away: no CSV, no
    # summary, no validation report. The `auto` and `all` aliases hit the same
    # branch.
    from .search_control import SearchBudget

    search_budget = SearchBudget.from_config(config)
    kwargs["search_budget"] = search_budget
    kwargs["candidate_source"] = _candidate_source
    result, optimizer = _dispatch_optimizer(
        method,
        cache,
        candidates,
        fg_prefixes,
        fg_seq_lengths,
        bg_prefixes,
        bg_seq_lengths,
        target_size,
        config,
        conditions,
        seed,
        verbose,
        _application,
        kwargs,
    )
    if _candidate_source is not None:
        candidates = list(_candidate_source.frontier())
    if result is not None and getattr(result, "status", None) is OptimizationStatus.ERROR:
        return result

    # --minimize-primers / --target-coverage post-processing.
    #
    # Both flags were reaching this function and going no further: they were
    # stored on OptimizationConfig, forwarded through **kwargs, and read by no
    # optimizer. `MinimalPrimerSelector` -- which exists precisely to trim a set
    # to the smallest one still meeting a coverage target -- had no callers
    # anywhere outside its own module. So `--minimize-primers` did nothing,
    # which also made the default path look like it was minimising when the
    # real cause was the coverage-binning bug in dominating_set_optimizer.
    if _minimize_primers and result.primers and optimizer is None:
        logger.warning(
            "--minimize-primers was requested but no optimizer is available to "
            "measure coverage with; the set is returned untrimmed."
        )
    elif _minimize_primers and result.primers:
        from .optimization_service import OptimizationRequest, run_panel_search

        result = run_panel_search(
            OptimizationRequest(
                optimizer,
                tuple(candidates),
                target_size,
                minimize=True,
                target_coverage=_target_coverage,
                budget=search_budget,
            ),
            initial_result=result,
        )

    # Resolve any remaining repair before derived coverage and validation reports.
    result, acceptance = _hold_to_configured_limits(
        result, optimizer, candidates, config, verbose, budget=search_budget
    )

    # Alternative sets, when `max_sets` asks for more than one. Done here, on
    # the finished primary result, so an alternative is a genuinely different
    # set of oligos rather than a reordering.
    global _LAST_PRIMER_SETS, _LAST_CANDIDATE_REACH
    _LAST_PRIMER_SETS = [tuple(result.primers)] if result.primers else []

    _LAST_CANDIDATE_REACH = _report_candidate_reach(candidates, verbose)
    _max_sets = int(kwargs.get("max_sets") or getattr(parameter, "max_sets", 1) or 1)
    if _max_sets > 1 and result.primers and optimizer is not None:
        _LAST_PRIMER_SETS = collect_alternative_sets(
            primary=result,
            optimizer=optimizer,
            candidates=candidates,
            target_size=target_size,
            max_sets=_max_sets,
            max_iterations=getattr(config, "max_iterations", 8),
        )
        if verbose and len(_LAST_PRIMER_SETS) > 1:
            logger.info(f"max_sets={_max_sets}: offering {len(_LAST_PRIMER_SETS)} distinct sets")

    # Carry the unindexed-candidate count onto the result so it reaches
    # step4_improved_df_summary.json even on the paths that do not refuse.
    if _unindexed:
        from dataclasses import replace as _dc_replace

        result = _dc_replace(result, unindexed_candidates=_unindexed)

    # Phase 15A: populate per_target_coverage on the result so multi-
    # target runs surface "target A 95% / target B 40%" instead of one
    # aggregate. Computed in the caller (here) rather than in every
    # optimizer so all 16 methods get the same treatment uniformly. The
    # `dataclasses.replace` preserves frozen-dataclass immutability.
    _saturation_warnings: list = []  # collected here, appended to validation dict below
    _coverage_extension: int = extension_reach
    try:
        if result.primers and fg_prefixes:
            from dataclasses import replace as _dc_replace

            from .coverage import compute_per_prefix_coverage

            # The reach resolved above, which honours `coverage_reach` from
            # params.json or `--coverage-reach`. This used to call
            # `polymerase_extension_reach` again and discard the override, so a
            # run selected at 8 kb and reported at phi29's 3 kb default -- about
            # a third of the coverage the set was actually chosen for, with
            # nothing saying the two figures came from different assumptions.
            _ext = extension_reach
            _coverage_extension = _ext
            _, per_target = compute_per_prefix_coverage(
                cache=cache,
                primers=list(result.primers),
                prefixes=fg_prefixes,
                seq_lengths=fg_seq_lengths,
                extension=_ext,
                circular=fg_circular,
            )
            if per_target:
                new_metrics = _dc_replace(result.metrics, per_target_coverage=per_target)
                result = _dc_replace(result, metrics=new_metrics)

                # Small-genome coverage saturation check. When
                # (num_primers * 2 * extension) >= genome_len, a primer
                # set can trivially cover the genome regardless of design
                # quality — the coverage number is saturation-bounded,
                # not a real property. Warn so the user does not mistake
                # plasmid-scale scenarios for "perfect design".
                n_primers = len(result.primers)
                if n_primers > 0 and _ext > 0:
                    for prefix, length in zip(fg_prefixes, fg_seq_lengths, strict=True):
                        if length <= 0:
                            continue
                        expected_window = n_primers * 2 * _ext
                        saturation = length / expected_window if expected_window else 1.0
                        cov = per_target.get(prefix, 0.0)
                        if saturation < 1.0 and cov >= 0.95:
                            _saturation_warnings.append(
                                {
                                    "level": "warning",
                                    "code": "coverage_saturated_on_small_genome",
                                    "detail": (
                                        f"{prefix}: genome={length} bp, "
                                        f"{n_primers} primers x 2x {_ext} bp "
                                        f"reach = {expected_window} bp; "
                                        f"coverage={cov:.1%} is saturation-"
                                        f"bounded (metric unreliable)"
                                    ),
                                }
                            )
    except DesignError:
        raise
    except Exception as exc:
        # `base_optimizer` gates the per-target floor on a non-empty dict, so
        # an empty one makes a REQUESTED floor pass vacuously: the panel that
        # starves one target is reported as satisfying the limit set to stop
        # exactly that. Aggregate coverage cannot reveal it either, which is
        # why the floor exists.
        raise ModelEvaluationError("per_target_coverage", sorted(fg_prefixes), str(exc)) from exc

    # Post-optimization sanity validation. Catches duplicates, size drift,
    # zero coverage, and accidental blacklist re-injection. The result is
    # written to data_dir/step4_improved_df_validation.json and logged here
    # so downstream CSV consumers can trust the output shape.
    try:
        forbidden = _collect_forbidden_primers(candidates, verbose)

        validation = result.validate(
            target_size=target_size,
            min_per_target_coverage=per_target_floor(kwargs, parameter) or 0.0,
            forbidden_primers=forbidden or None,
        )
    except DesignError:
        raise
    except Exception as exc:
        # Skipping the validator disarms every post-optimization guard at once
        # -- duplicates, size drift, zero coverage, blacklist re-injection and
        # the delivered-panel dimer check -- and writes no validation file, so
        # `neoswga export` finds nothing to block on and prints "ready for
        # ordering".
        raise ModelEvaluationError("post-optimization validation", method, str(exc)) from exc

    # Attach saturation warnings so the HTML report surfaces them via
    # the validator-banner path. Warnings only; the `ok` flag is
    # unchanged because saturation is not a correctness bug.
    if validation is not None and _saturation_warnings:
        validation.setdefault("issues", []).extend(_saturation_warnings)
        if verbose:
            for w in _saturation_warnings:
                logger.warning(f"{w['code']}: {w['detail']}")

    # Screen the pool that is actually being delivered against the threshold
    # the user configured. The optimizers penalise dimers (except clique, which
    # constrains them), but a penalty can be outweighed, so the only way to know
    # what came back is to measure it.
    if validation is not None:
        validation["dimer_policy"] = (
            "allow-relaxation" if config.allow_dimer_relaxation else "strict"
        )
        validation["refinement_method"] = config.refinement_method
        validation["swap_max_evaluations"] = config.swap_max_evaluations
        validation["swap_max_seconds"] = config.swap_max_seconds
        validation["requested_primers"] = target_size
        validation["delivered_primers"] = len(result.primers)
    if validation is not None and result.primers:
        _dimer_issue = dimer_validation_issue(
            list(result.primers), getattr(config, "max_dimer_bp", None)
        )
        if _dimer_issue is not None:
            validation.setdefault("issues", []).append(_dimer_issue)
            if verbose:
                logger.warning(f"{_dimer_issue['code']}: {_dimer_issue['detail']}")

    # Reported, not repaired; `panel_acceptance` says why.
    if verbose:
        report_per_target(
            check_per_target_coverage(result.metrics, per_target_floor(kwargs, parameter))
        )

    if validation is not None and design_request is not None and result.primers:
        validation["assessment"] = panel_assessment(design_request, result, optimizer)

    _limit_issue = limit_violation_issue(acceptance)
    if validation is not None and _limit_issue is not None:
        validation.setdefault("issues", []).append(_limit_issue)

    if validation is not None:
        # `result.validate` folded `ok` out of the issues it had built, and
        # the saturation warnings and the dimer finding are appended above.
        # Recompute from the final list or a blocking finding cannot move it.
        validation["ok"] = panel_validation_is_ok(validation.get("issues") or [])
        _write_validation_report(validation)

    if verbose:
        log_optimization_outcome(result, target_size, method, validation)

    return result


def _hold_to_configured_limits(result, optimizer, candidates, config, verbose, budget=None):
    """Report configured panel limits and resolve any remaining bounded repair.

    Returns `(result, report)`. The report is None when this run configured no
    limit, which is absence rather than compliance: with none configured
    `constraints_from_parameter` returns None, no objective is built and
    nothing here is evaluated.
    """
    constraints = constraints_from_parameter(parameter)
    if constraints is None or not result.primers or optimizer is None:
        return result, None
    replacement, report = apply_configured_limits(
        result,
        optimizer,
        candidates=candidates,
        config=config,
        constraints=constraints,
        verbose=verbose,
        background_available=bool(getattr(parameter, "bg_prefixes", None)),
        budget=budget,
    )
    return (replacement if replacement is not None else result), report


def run_optimization_from_config(config: OptimizationConfig) -> OptimizationResult:
    """
    Run optimization using a config object.

    Args:
        config: OptimizationConfig with all settings

    Returns:
        OptimizationResult with selected primers
    """
    return run_optimization(
        method=config.method,
        allow_dimer_relaxation=config.allow_dimer_relaxation,
        refinement_method=config.refinement_method,
        swap_max_evaluations=config.swap_max_evaluations,
        swap_max_seconds=config.swap_max_seconds,
        stage1_objective_width=config.stage1_objective_width,
        fg_prefixes=config.fg_prefixes,
        fg_seq_lengths=config.fg_seq_lengths,
        bg_prefixes=config.bg_prefixes,
        bg_seq_lengths=config.bg_seq_lengths,
        target_size=config.target_set_size,
        verbose=config.verbose,
        max_iterations=config.max_iterations,
        # Nine of the eighteen fields this dataclass declares were dropped here,
        # so a library caller writing
        #     OptimizationConfig(minimize_primers=True, target_coverage=0.9)
        # got neither. That also contradicted the comment justifying how
        # `minimize_primers` is handled in `run_optimization`, which says the
        # flag "lives on the separate OptimizationConfig used by
        # run_optimization_from_config" -- it does live there, and this did not
        # forward it.
        minimize_primers=config.minimize_primers,
        target_coverage=config.target_coverage,
        uniformity_weight=config.uniformity_weight,
        use_cache=config.use_position_cache,
        use_background_filter=config.use_background_filter,
    )


def _simulation_rescore(
    result: OptimizationResult,
    fg_prefixes: List[str],
    fg_seq_lengths: List[int],
    simulation_time: float = 1800.0,
    verbose: bool = True,
) -> Optional[Dict[str, Any]]:
    """Re-score optimization result using simulation-based fitness.

    Uses the agent-based replication simulator to predict actual
    amplification performance of the selected primer set.

    Args:
        result: OptimizationResult from an optimizer.
        fg_prefixes: Foreground genome prefixes.
        fg_seq_lengths: Foreground genome lengths.
        simulation_time: Simulation duration in seconds.
        verbose: Print progress.

    Returns:
        Dictionary with simulation fitness metrics, or None on failure.
    """
    if not result.is_success or not result.primers:
        return None

    try:
        from .genome_io import GenomeLoader
        from .position_cache import PositionCache
        from .simulation_fitness import SimulationBasedEvaluator
    except ImportError as e:
        if verbose:
            logger.warning(f"Simulation modules not available: {e}")
        return None

    try:
        # Load genome sequence for simulation
        genome_seq = None
        genome_length = sum(fg_seq_lengths)

        fg_genomes = getattr(parameter, "fg_genomes", [])
        if fg_genomes:
            try:
                loader = GenomeLoader()
                genome_seq = loader.load_genome(fg_genomes[0], return_stats=False)
                genome_length = len(genome_seq)
            except Exception:
                pass

        if genome_seq is None:
            if verbose:
                logger.info("  Genome sequence not available for simulation")
            return None

        # Build position cache for simulation
        primers = list(result.primers)
        cache = PositionCache(fg_prefixes, primers)

        evaluator = SimulationBasedEvaluator(
            genome_sequence=genome_seq,
            genome_length=genome_length,
            position_cache=cache,
            n_replicates=2,
            simulation_duration=simulation_time,
        )

        fitness = evaluator.evaluate(primers, verbose=verbose)

        return {
            "simulation_coverage": fitness.mean_coverage,
            "simulation_uniformity": fitness.coverage_uniformity,
            "simulation_fitness": fitness.fitness_score,
            "simulation_forks": fitness.mean_forks_created,
            "simulation_time": fitness.simulation_time,
        }
    except Exception as e:
        if verbose:
            logger.warning(f"Simulation re-scoring failed: {e}")
        return None


def optimize_step4(
    use_cache: bool = True,
    use_background_filter: bool = True,
    optimization_method: str = "hybrid",
    verbose: bool = True,
    uniformity_weight: Optional[float] = None,
    minimize_primers: bool = False,
    target_coverage: float = 0.70,
    design_request=None,
    **kwargs,
) -> Tuple[List[List[str]], List[float], Any]:
    """
    Drop-in replacement for improved_step4 from pipeline_integration.

    Maintains backward compatibility with existing CLI while using
    the new optimizer framework internally.

    Args:
        use_cache: True builds the in-memory `PositionCache`; False builds the
            memory-mapped `StreamingPositionCache` instead, for a background
            too large to hold. This said "always True with new system" and was
            never read, so `--no-position-cache` changed nothing while the CLI
            logged "Position cache: False".
        use_background_filter: Use background filtering
        optimization_method: Optimizer method name
        verbose: Print progress
        uniformity_weight: Weight for coverage uniformity. None leaves it to
            the `--application` profile (see _resolve_selection_weights);
            0.0 is an explicit request for no uniformity term.
        minimize_primers: Minimize primer count
        target_coverage: Target coverage for minimization
        **kwargs: Additional parameters

    Returns:
        Tuple of (primer_sets, scores, cache) for CLI compatibility
    """
    _ensure_optimizers_registered()

    # Reset the last-result stash at entry. run_optimization also resets,
    # but optimize_step4 may early-return via no-candidates path before the
    # run_optimization call so we clear explicitly here too.
    global _LAST_RESULT
    _LAST_RESULT = None

    # Get target size: prefer explicit kwarg, then parameter module, then default
    target_size = kwargs.pop("target_size", None)
    if target_size is None:
        target_size = getattr(parameter, "num_primers", 6)
        target_size = getattr(parameter, "target_set_size", target_size)

    if verbose:
        logger.info(f"Unified optimizer: method={optimization_method}")

    # Run optimization
    result = run_optimization(
        method=optimization_method,
        target_size=target_size,
        verbose=verbose,
        uniformity_weight=uniformity_weight,
        minimize_primers=minimize_primers,
        target_coverage=target_coverage,
        # Named on this function's signature, so it is NOT in **kwargs and has
        # to be forwarded explicitly. Being absent here is why the flag stayed
        # inert even once `run_optimization` knew what to do with it.
        use_cache=use_cache,
        design_request=design_request,
        **kwargs,
    )

    # Optional simulation-based re-scoring
    if kwargs.get("validate_with_simulation", False) and result.is_success:
        if verbose:
            logger.info("")
            logger.info("=" * 60)
            logger.info("Simulation-Based Validation")
            logger.info("=" * 60)

        # Resolve genome prefixes for simulation
        sim_fg_prefixes = getattr(parameter, "fg_prefixes", [])
        sim_fg_seq_lengths = getattr(parameter, "fg_seq_lengths", [])

        sim_time = kwargs.get("simulation_time", 1800.0)
        sim_results = _simulation_rescore(
            result,
            sim_fg_prefixes,
            sim_fg_seq_lengths,
            simulation_time=sim_time,
            verbose=verbose,
        )

        if sim_results:
            logger.info(f"  Simulated coverage: {sim_results['simulation_coverage']:.1%}")
            logger.info(f"  Simulated uniformity: {sim_results['simulation_uniformity']:.2f}")
            logger.info(f"  Simulation fitness: {sim_results['simulation_fitness']:.3f}")

    # Stash for the CLI / pareto-front rendering to consult.
    # (The `global _LAST_RESULT` declaration is at the top of optimize_step4;
    # it is not re-declared here because Python does not allow re-declaration
    # of a global within the same function.)
    _LAST_RESULT = result

    # Convert to legacy format
    if result.is_success or (result.status == OptimizationStatus.PARTIAL and result.primers):
        # Every set the run produced, best first. This was always a
        # single-entry list, which is why `max_sets` had nothing to change.
        primer_sets = [list(s) for s in (_LAST_PRIMER_SETS or [result.primers])]
        scores = [result.score]

        # Save results
        output_path = os.path.join(parameter.data_dir, "step4_improved_df.csv")
        # The profile this run was scored under, so the CSV's
        # `normalized_score` column agrees with the ensemble comparison
        # table rather than silently using the balanced weights.
        # Which criterion limits the delivered panel, printed and recorded.
        # The two wet-lab benchmarks disagree about which spacing property
        # predicts success, so the useful thing to report is which one is
        # currently binding rather than a weighted blend of all of them.
        # `None` when a reference could not be resolved, and then nothing is
        # printed rather than a limit nothing established.
        regime = assess_from_parameter(
            result.metrics,
            parameter,
            delivered_size=len(result.primers),
            application=kwargs.get("application"),
        )
        log_regime(regime)

        save_results(
            result,
            output_path,
            application=kwargs.get("application"),
            primer_sets=_LAST_PRIMER_SETS or None,
            regime=regime,
            candidate_reach=_LAST_CANDIDATE_REACH,
        )

        # Return cache placeholder (for compatibility)
        from .position_cache import PositionCache

        fg_prefixes = getattr(parameter, "fg_prefixes", [])
        cache = PositionCache(fg_prefixes, list(result.primers)) if fg_prefixes else None

        return primer_sets, scores, cache
    else:
        logger.error(f"Optimization failed: {result.message}")
        return [], [], None


# CLI entry point
def main():
    """CLI entry point for standalone optimization."""
    import argparse

    parser = argparse.ArgumentParser(description="NeoSWGA Primer Set Optimization")
    parser.add_argument("-m", "--method", default="hybrid", help="Optimization method")
    parser.add_argument("-j", "--json-file", help="Parameters JSON file")
    parser.add_argument("-n", "--num-primers", type=int, default=6, help="Target number of primers")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument(
        "--list-methods", action="store_true", help="List available optimization methods"
    )

    args = parser.parse_args()

    if args.list_methods:
        print("Available optimization methods:")
        for name, desc in list_available_optimizers().items():
            print(f"  {name}: {desc}")
        return

    # Load parameters
    if args.json_file:
        parameter.json_file = args.json_file
        from . import pipeline as core_pipeline

        core_pipeline._initialize()

    # Run optimization
    result = run_optimization(
        method=args.method,
        target_size=args.num_primers,
        verbose=args.verbose,
    )

    if result.is_success:
        print(f"\nSelected {result.num_primers} primers:")
        for primer in result.primers:
            print(f"  {primer}")
        print(f"\nScore: {result.score:.4f}")
        print(f"Coverage: {result.metrics.fg_coverage:.1%} (measured)")
    else:
        print(f"Optimization failed: {result.message}")
        exit(1)


if __name__ == "__main__":
    main()
