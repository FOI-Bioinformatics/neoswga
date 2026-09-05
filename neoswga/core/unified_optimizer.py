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
from .dimer import dimer_validation_issue, worst_heterodimer
from .optimizer_factory import OptimizerFactory, OptimizerRegistry
from .position_cache import PositionCache, StreamingPositionCache
from .progress import progress_context
from .step4_output import save_results

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


def _prefilter_by_background(
    cache,
    candidates,
    fg_prefixes,
    bg_prefixes,
    min_ratio=1.0,
    max_removal_fraction=0.20,
    verbose=False,
):
    """Remove candidates with poor foreground/background binding ratio.

    Keeps all primers with ratio >= min_ratio. If that would remove more
    than max_removal_fraction of candidates, keeps the top
    (1 - max_removal_fraction) by ratio instead.

    Args:
        cache: PositionCache with loaded positions
        candidates: List of candidate primer sequences
        fg_prefixes: Foreground genome prefixes
        bg_prefixes: Background genome prefixes
        min_ratio: Minimum fg/bg ratio to keep (default 1.0)
        max_removal_fraction: Maximum fraction of candidates to remove (default 0.20)
        verbose: Log filtering details

    Returns:
        Filtered list of candidates (never empty)
    """
    ratios = []
    for primer in candidates:
        fg_count = sum(len(cache.get_positions(p, primer, "both")) for p in fg_prefixes)
        bg_count = sum(len(cache.get_positions(p, primer, "both")) for p in bg_prefixes)
        ratios.append((fg_count / (bg_count + 1), primer))

    # Keep primers above threshold
    above = [(r, p) for r, p in ratios if r >= min_ratio]

    # Ensure we don't remove too many
    min_keep = int(len(candidates) * (1 - max_removal_fraction))
    if len(above) < min_keep:
        ratios.sort(reverse=True)
        above = ratios[:min_keep]

    filtered = [p for _, p in above]
    removed = len(candidates) - len(filtered)
    if removed > 0 and verbose:
        logger.info(f"  Background pre-filter: removed {removed} candidates with poor fg/bg ratio")

    return filtered if filtered else candidates  # never return empty


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

        set_kmer_sampling_seed(seed)
    except Exception:
        pass


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


def collect_alternative_sets(
    primary, optimizer, candidates, target_size, max_sets=1, max_iterations=8
):
    """Up to `max_sets` distinct primer sets, best first.

    Alternatives are found by removing the primers already chosen from the
    candidate pool and running selection again, so each one is a genuinely
    different set rather than a reordering of the same oligos. `max_iterations`
    caps how many such attempts are made, which matters because a pool can run
    out of usable candidates long before `max_sets` is reached.

    Fewer than `max_sets` is a normal outcome, not a failure: a small pool
    simply cannot yield many disjoint sets. The primary result is always first.

    Both parameters were documented and inert -- `max_sets` reached only
    `search_context.BFSConfig`, which has no callers, and no optimizer reads
    `config.max_iterations`. The output format already anticipated this: the
    `set_index` column of step4_improved_df.csv was hardcoded to 0.

    Note `max_iterations` bounds the search for ALTERNATIVES only. Bounding the
    primary selection with it would cap how many primers a run can choose, so
    `iterations: 8` would quietly truncate a 96-oligo panel.
    """
    sets = [tuple(primary.primers)]
    if not primary.primers or max_sets <= 1:
        return sets

    seen = {frozenset(primary.primers)}
    remaining = [p for p in candidates if p not in set(primary.primers)]
    attempts = 0

    while len(sets) < max_sets and attempts < max(1, int(max_iterations)):
        attempts += 1
        if len(remaining) < target_size:
            break
        try:
            result = optimizer.optimize(remaining, target_size)
        except Exception as exc:
            logger.debug(f"Alternative set search stopped: {exc}")
            break

        chosen = tuple(result.primers)
        if not chosen or frozenset(chosen) in seen:
            break

        sets.append(chosen)
        seen.add(frozenset(chosen))
        remaining = [p for p in remaining if p not in set(chosen)]

    return sets


def _select_ensemble_winner(results, application):
    """Which method's set to deliver, decided by a stated rule.

    Ranked on `normalized_score` (application-weighted), then -- and this is the
    part that was missing -- on explicit tie-breaks rather than on whichever
    method happened to be listed first.

    `max()` over the results dict returns the FIRST key attaining the maximum,
    and that dict is built in `--ensemble-methods` order, so flag order decided
    every tie. Ties are the common case: on the plasmid example hybrid,
    dominating-set and background-aware return the identical set and score
    0.9003 alike, and reordering the flag moved the winner each time. Here the
    tied methods returned the same primers so only provenance moved, but the tie
    is on the composite score rather than on the set, so two methods returning
    different pools can tie and then flag order picks what gets ordered.

    Tie-breaks, in order:

    1. the higher `normalized_score`;
    2. the SMALLER set -- fewer oligos at equal quality is cheaper to
       synthesise and simpler to run;
    3. the method name, alphabetically, purely so the answer is stable.

    Rule 3 is arbitrary on purpose. Something has to settle a true tie, and an
    arbitrary rule that is written down and reproducible is better than an
    arbitrary rule that depends on argument order and is documented as
    order-independent.
    """

    def rank(method):
        result = results[method]
        return (
            result.metrics.normalized_score(application=application),
            -len(result.primers),
            # `max` takes the largest, so invert the name to prefer the earliest.
            tuple(-ord(c) for c in method),
        )

    return max(results, key=rank)


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
    """Run several optimizers on one shared cache; keep the best.

    Methods are compared by ``normalized_score`` (a [0,1] value comparable
    across optimizer types, weighted by ``application``) — NOT raw ``score``,
    which is on per-optimizer scales. A per-method comparison table is attached
    to the returned result as ``ensemble_comparison``. Per-method failures are
    logged and skipped (graceful degradation); if every method fails the result
    is an ERROR.

    Building the PositionCache once (passed in as ``cache``) is the key speed
    win over calling run_optimization N times.
    """
    from dataclasses import replace as _dc_replace

    rows: List[Dict[str, Any]] = []
    results: Dict[str, OptimizationResult] = {}

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
            with progress_context(f"  ensemble: {optimizer.name}", disable=not verbose):
                res = optimizer.optimize(candidates, target_size)
        except Exception as e:
            logger.warning(f"  ensemble method '{m}' failed: {e}")
            rows.append(
                {
                    "method": m,
                    "normalized_score": 0.0,
                    "score": float("-inf"),
                    "n_primers": 0,
                    "fg_coverage": 0.0,
                    "bg_coverage": 0.0,
                    "status": "error",
                    "selected": False,
                }
            )
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

    best_method = _select_ensemble_winner(results, application)
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
                combined = combiner.optimize(union_pool, target_size)
                c_norm = combined.metrics.normalized_score(application=application)
                w_norm = winner.metrics.normalized_score(application=application)
                if combined.primers and c_norm >= w_norm:
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

    return config_class(
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
    except Exception as e:
        logger.debug(f"application weight lookup skipped: {e}")


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
        forbidden = [p for p, f in zip(candidates, freqs) if f > max_bl]
        if forbidden and verbose:
            logger.warning(
                f"Library blacklist guard: {len(forbidden)} candidate(s) exceed "
                f"max_bl_freq={max_bl}; validator will flag any that reached the set."
            )
        return forbidden
    except Exception as e:
        logger.debug(f"library blacklist guard skipped ({e})")
        return []


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
            logger.debug(f"No measuring optimizer for the ensemble winner: {exc}")
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

        # Run optimization
        with progress_context(f"Running {optimizer.name} optimizer", disable=not verbose):
            result = optimizer.optimize(candidates, target_size)

    return result, optimizer


def _minimize_primer_count(result, optimizer, target_coverage: float, verbose: bool):
    """Trim a selected set to the fewest primers still meeting a coverage target.

    Repeatedly drops whichever primer costs the least coverage, for as long as
    what remains still clears `target_coverage`. This is deliberately a
    post-process rather than a selection criterion: the optimizer picks the best
    set it can at the requested size, and this asks separately whether a smaller
    subset would still do the job.

    Coverage is measured with `optimizer.compute_metrics`, the same base-level
    figure written to `step4_improved_df_summary.json`. That matters more than
    it looks. `MinimalPrimerSelector`, the module apparently written for this
    job and never called by anything, counts `covered_positions` as the set of
    binding-site COORDINATES -- so a primer covers as many bases as it has
    sites, and extension is ignored entirely. On a 30 kb genome, 30 sites reads
    as 0.1% coverage and no target is ever reachable. Using it here would have
    given the flag a criterion that cannot fire. One coverage semantics, and it
    is the reported one.

    Guarded in both directions: a trimmed set is accepted only if it is
    genuinely smaller and still clears the target, so the flag can shrink a set
    but never silently degrade one.
    """
    from dataclasses import replace as _replace

    try:
        current = list(result.primers)
        if len(current) < 2:
            return result

        best_metrics = None
        while len(current) > 1:
            candidate_drops = []
            for primer in current:
                remaining = [p for p in current if p != primer]
                metrics = optimizer.compute_metrics(remaining)
                if metrics.fg_coverage >= target_coverage:
                    candidate_drops.append((metrics.fg_coverage, primer, remaining, metrics))

            if not candidate_drops:
                break

            # Drop the primer whose removal leaves the most coverage standing.
            _cov, _primer, remaining, metrics = max(candidate_drops, key=lambda d: d[0])
            current, best_metrics = remaining, metrics

        if len(current) >= len(result.primers):
            if verbose:
                logger.info(
                    f"--minimize-primers: no smaller subset holds "
                    f"{target_coverage:.0%} coverage; keeping {len(result.primers)}."
                )
            return result

        if verbose:
            logger.info(
                f"--minimize-primers: {len(result.primers)} -> {len(current)} primers "
                f"at {best_metrics.fg_coverage:.1%} coverage (target {target_coverage:.0%})."
            )
        # `score` has to move with the set. Replacing primers and metrics and
        # leaving `score` alone reported the FOUR-primer score for a set trimmed
        # to one, in both the CSV `score` column and the summary JSON -- a
        # number describing a set that no longer exists.
        #
        # `normalized_score` is the comparable [0,1] composite and is derived
        # from the metrics, which have been recomputed, so it is the honest
        # value to carry. The per-optimizer raw score cannot be recomputed here
        # without re-running the optimizer on the trimmed set.
        return _replace(
            result,
            primers=tuple(current),
            metrics=best_metrics,
            score=float(best_metrics.normalized_score()),
        )

    except Exception as e:
        logger.warning(f"--minimize-primers skipped ({e}); keeping the full set.")
        return result


def run_optimization(
    method: str = "hybrid",
    candidates: Optional[List[str]] = None,
    fg_prefixes: Optional[List[str]] = None,
    fg_seq_lengths: Optional[List[int]] = None,
    bg_prefixes: Optional[List[str]] = None,
    bg_seq_lengths: Optional[List[int]] = None,
    target_size: int = 6,
    verbose: bool = True,
    **kwargs,
) -> OptimizationResult:
    """
    Run primer set optimization using specified method.

    This is the main entry point for optimization. It:
    1. Loads candidates if not provided
    2. Creates position cache
    3. Instantiates optimizer via factory
    4. Runs optimization
    5. Returns typed result

    Args:
        method: Optimizer method name ('greedy', 'dominating-set', etc.)
        candidates: List of candidate primers (loaded from step3 if None)
        fg_prefixes: Foreground genome HDF5 prefixes
        fg_seq_lengths: Foreground genome lengths
        bg_prefixes: Background genome HDF5 prefixes (optional)
        bg_seq_lengths: Background genome lengths (optional)
        target_size: Desired primer set size
        verbose: Print progress information
        **kwargs: Additional optimizer-specific parameters

    Returns:
        OptimizationResult with selected primers and metrics

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

    # Validate user-supplied reaction conditions BEFORE loading candidates
    # so a params.json bounds error (e.g. formamide>10%) is reported as
    # the specific config problem it is — not masked by a downstream
    # "empty candidate pool" message. The failure is non-fatal: the
    # optimizer still runs (additive-blind), but the user sees it in both
    # the log and the validator-banner.
    conditions = kwargs.pop("conditions", None)
    _conditions_init_error: str = ""
    if conditions is None:
        try:
            from .reaction_conditions import build_reaction_conditions

            conditions = build_reaction_conditions()
        except (ValueError, TypeError, KeyError) as e:
            _conditions_init_error = str(e)
            logger.error(
                f"ReactionConditions construction failed: {e}. "
                f"Optimizer will run additive-BLIND — Tm calculations will "
                f"fall back to the Wallace formula and your DMSO / betaine / "
                f"formamide / etc. concentrations will be ignored. Check "
                f"that your params.json values are within the bounds listed "
                f"in `neoswga suggest --help` (e.g. DMSO 0-10%, betaine "
                f"0-2.5M, formamide 0-10%)."
            )
            conditions = None

    # Load candidates from step3 if not provided
    if candidates is None:
        step3_path = os.path.join(parameter.data_dir, "step3_df.csv")
        if not os.path.exists(step3_path):
            return OptimizationResult.failure(method, f"Step 3 output not found: {step3_path}")
        step3_df = pd.read_csv(step3_path)
        candidates = step3_df["primer"].tolist()

    # Empty candidate pool guard. If filter/score removed every primer,
    # downstream optimizers behave inconsistently (some crash, some
    # return empty results silently). Fail with an actionable message
    # instead of dispatching into the factory with zero candidates.
    if not candidates:
        msg = (
            "No candidate primers available for optimization. "
            "Common causes: (1) `neoswga filter` produced an empty "
            "step2_df.csv (try relaxing max_bg_freq / max_gini or "
            "widening min_k-max_k); (2) `neoswga score` filtered all "
            "candidates below min_amp_pred (try lowering min_amp_pred); "
            "(3) caller passed an empty candidate list. Status: "
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
        try:
            from .rf_preprocessing import set_kmer_sampling_seed

            set_kmer_sampling_seed(seed)
        except Exception as e:
            logger.debug(f"Could not seed k-mer sampling RNG: {e}")
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

    # Optional: pre-filter candidates by fg/bg binding ratio
    if bg_prefixes and kwargs.get("bg_prefilter", True):
        candidates = _prefilter_by_background(
            cache,
            candidates,
            fg_prefixes,
            bg_prefixes,
            # `or` rather than a `.get` default: the CLI now forwards these,
            # and forwards None when the flag was not given, which a bare
            # `.get(name, default)` would hand straight through as None.
            min_ratio=kwargs.get("bg_min_ratio") or 1.0,
            max_removal_fraction=kwargs.get("bg_max_removal") or 0.20,
            verbose=verbose,
        )

    # Determine polymerase extension reach for coverage computation.
    # Uses realistic per-primer reach (phi29 ~3 kb, equiphi29 ~4 kb), not
    # single-molecule processivity. See coverage.polymerase_extension_reach
    # for the Clarke 2017 / Dwivedi-Yu 2023 rationale.
    # An explicit `coverage_reach` (params.json or --coverage-reach) overrides
    # it. The right value is an empirical property of the reaction, and
    # `neoswga calibrate-reach --bam` estimates it from sequencing depth.
    polymerase = kwargs.get("polymerase") or getattr(parameter, "polymerase", "phi29")
    reach_override = kwargs.get("coverage_reach")
    if reach_override is None:
        reach_override = getattr(parameter, "coverage_reach", None)
    try:
        from .coverage import resolve_coverage_reach

        extension_reach = resolve_coverage_reach(polymerase, override=reach_override)
    except ImportError:
        extension_reach = 3000  # phi29 realistic per-primer reach
    if reach_override is not None and verbose:
        logger.info(
            f"Coverage reach: {extension_reach:,} bp (explicit; polymerase default "
            f"would be used otherwise). Coverage figures are not comparable across "
            f"different reaches."
        )

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

    # By this point `conditions` is either a valid object or None (with
    # `_conditions_init_error` populated for validator-banner reporting
    # further down).

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
        result = _minimize_primer_count(
            result=result,
            optimizer=optimizer,
            target_coverage=_target_coverage,
            verbose=verbose,
        )

    # Alternative sets, when `max_sets` asks for more than one. Done here, on
    # the finished primary result, so an alternative is a genuinely different
    # set of oligos rather than a reordering.
    global _LAST_PRIMER_SETS
    _LAST_PRIMER_SETS = [tuple(result.primers)] if result.primers else []
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
                    for prefix, length in zip(fg_prefixes, fg_seq_lengths):
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
    except Exception as e:
        logger.debug(f"per_target_coverage population skipped ({e})")

    # Post-optimization sanity validation. Catches duplicates, size drift,
    # zero coverage, and accidental blacklist re-injection. The result is
    # written to data_dir/step4_improved_df_validation.json and logged here
    # so downstream CSV consumers can trust the output shape.
    try:
        forbidden = _collect_forbidden_primers(candidates, verbose)

        validation = result.validate(
            target_size=target_size,
            min_coverage=0.0,  # soft by default; caller can tighten
            min_per_target_coverage=kwargs.get("min_per_target_coverage", 0.0),
            forbidden_primers=forbidden or None,
        )
    except Exception as e:
        logger.debug(f"Post-optimization validator crashed ({e}); skipping")
        validation = None

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
    if validation is not None and result.primers:
        _dimer_issue = dimer_validation_issue(
            list(result.primers), getattr(config, "max_dimer_bp", None)
        )
        if _dimer_issue is not None:
            validation.setdefault("issues", []).append(_dimer_issue)
            if verbose:
                logger.warning(f"{_dimer_issue['code']}: {_dimer_issue['detail']}")

    # Surface the ReactionConditions init failure in the banner too.
    # Marked as a warning (not error) because the optimizer still
    # produced a primer set — just with additive-blind Tm.
    if validation is not None and _conditions_init_error:
        validation.setdefault("issues", []).append(
            {
                "level": "warning",
                "code": "reaction_conditions_init_failed",
                "detail": (
                    f"{_conditions_init_error}. Optimizer ran additive-blind; "
                    f"Tm calculations fell back to Wallace. Re-check params.json "
                    f"against the bounds documented in `neoswga suggest --help`."
                ),
            }
        )

    if validation is not None:
        try:
            import json as _json

            data_dir = getattr(parameter, "data_dir", None) or os.getcwd()
            out_path = os.path.join(data_dir, "step4_improved_df_validation.json")
            with open(out_path, "w") as fh:
                _json.dump(validation, fh, indent=2)
        except Exception as e:
            logger.debug(f"Could not write validation report ({e}); continuing")

    if verbose:
        if result.is_success:
            logger.info(f"Selected {result.num_primers} primers")
            logger.info(f"Coverage: {result.metrics.fg_coverage:.1%} (measured)")
            logger.info(f"Score: {result.score:.4f}")
        elif result.status == OptimizationStatus.PARTIAL:
            if result.num_primers < target_size:
                logger.warning(
                    f"PARTIAL result: found {result.num_primers} primers "
                    f"but target was {target_size}"
                )
            else:
                logger.warning(
                    f"PARTIAL result: found {result.num_primers} primers "
                    f"but genome coverage is below threshold"
                )
            logger.warning("Suggestions to improve results:")
            logger.warning("  - Relax filtering thresholds (increase max_bg_freq or max_gini)")
            logger.warning("  - Widen the k-mer range (decrease min_k or increase max_k)")
            logger.warning("  - Increase the candidate pool (raise max_primer in filter step)")
            logger.warning(f"  - Try a different optimizer (current: {method})")
        else:
            logger.warning(f"Optimization failed: {result.message}")

        if validation is not None and validation.get("issues"):
            for issue in validation["issues"]:
                lvl = issue.get("level", "warning")
                msg = f"Post-opt {lvl}: {issue.get('code')} — {issue.get('detail')}"
                (logger.error if lvl == "error" else logger.warning)(msg)

    return result


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
        save_results(
            result,
            output_path,
            application=kwargs.get("application"),
            primer_sets=_LAST_PRIMER_SETS or None,
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
