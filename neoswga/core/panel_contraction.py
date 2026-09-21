"""Contract an existing panel under the shared design objective."""

import numpy as np

from .base_optimizer import OptimizationResult, OptimizationStatus
from .optimization_service import OptimizationRequest, panel_violations, run_panel_search
from .panel_refinement import objective_for_optimizer


def scan_panel_positions(primers, references):
    """Read-only exact scans for supplied oligos, including external primers.

    Each reference supplies (prefix, FASTA, length, circularity). Scanning the
    small supplied panel avoids treating an unindexed primer as a zero-site
    primer and does not modify existing pipeline indexes.
    """
    from . import string_search
    from .position_cache import POSITION_DTYPE, PositionCache
    from .thermodynamics import reverse_complement

    cache = PositionCache([r[0] for r in references], [])
    cache.primers = set(primers)
    queries = set(primers) | {reverse_complement(p) for p in primers}
    by_length = {k: sorted(p for p in queries if len(p) == k) for k in {len(p) for p in queries}}
    for prefix, fasta, length, circular in references:
        sequence = string_search.get_cached_genome_sequence(fasta)
        if len(sequence) != length:
            raise ValueError(f"Reference length differs from configured length: {prefix}")
        if string_search.AHOCORASICK_AVAILABLE:
            positions = string_search.get_all_positions_multi_k(by_length, fasta, circular)
        else:
            positions = {}
            for query in by_length.values():
                positions.update(string_search.get_all_positions_per_k(query, fasta, circular))
        cache.record_starts[prefix] = [0] + string_search.get_cached_record_boundaries(fasta)
        for primer in primers:
            for strand, query in (("forward", primer), ("reverse", reverse_complement(primer))):
                cache.cache[(prefix, primer, strand)] = np.asarray(
                    positions[query], dtype=POSITION_DTYPE
                )
    return cache


def contract_panel(
    optimizer, primers, minimum_coverage, constraints=None, positions_source="provided_cache"
):
    """Keep the best qualifying subset found, reporting unmet targets explicitly."""
    objective = objective_for_optimizer(optimizer, constraints)
    metrics = optimizer.compute_metrics(primers)
    initial = OptimizationResult(
        tuple(primers),
        metrics.normalized_score(),
        OptimizationStatus.SUCCESS,
        metrics,
        0,
        "contract-set",
    )
    result = run_panel_search(
        OptimizationRequest(
            optimizer,
            tuple(primers),
            len(primers),
            constraints=constraints,
            minimize=True,
            target_coverage=minimum_coverage,
            refine=False,
            repair=False,
        ),
        initial_result=initial,
    )
    coverage = objective.coverage(result.primers)
    violations = panel_violations(optimizer, result.primers)
    removals = [removal for stage in result.stage_history for removal in stage.get("removals", [])]
    return dict(
        positions_source=positions_source,
        original_set=list(primers),
        contracted_set=list(result.primers),
        removed_primers=[r["primer"] for r in removals],
        removal_trace=removals,
        stage_history=list(result.stage_history),
        baseline_coverage=objective.coverage(primers),
        final_coverage=coverage,
        baseline_raw_coverage=metrics.fg_coverage,
        final_raw_coverage=result.metrics.fg_coverage,
        coverage_metric=objective.constraints.coverage_metric,
        min_coverage_threshold=minimum_coverage,
        meets_target=coverage is not None and coverage >= minimum_coverage and not violations,
        failed_constraints=list(violations),
        conditions_fingerprint=optimizer.conditions.fingerprint() if optimizer.conditions else None,
    )
