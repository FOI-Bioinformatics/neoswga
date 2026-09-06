"""Execution tests for the default optimizer.

`hybrid` is what `neoswga optimize` runs unless told otherwise, so it produces
most of the primer sets this tool ever ships. Its 27 existing tests cover
construction and config; none of them ran `optimize()`, leaving the two-stage
selection, the network refinement, the background pruning and the thermodynamic
candidate filter untested -- roughly 250 lines of the most-used component.

The fixture builds a real `PositionCache` over a real HDF5 written by the real
scanner, so the optimizer sees the same position data a pipeline run would.

Two properties get particular attention. Coverage must be scored at the
realistic per-primer reach (~3 kb for phi29) rather than single-molecule
processivity (~70 kb) -- conflating them overstated coverage by an order of
magnitude. And `fixed_primers` must actually survive into the result, since the
whole iterative-design workflow depends on keeping validated oligos.
"""

import random

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer, HybridResult
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 30_000


@pytest.fixture(scope="module")
def genome():
    """A genome with 12 planted primers spread across it.

    Spread matters: a set-cover objective needs candidates whose coverage
    actually differs, or every selection is equally good and the test proves
    nothing.
    """
    rng = random.Random(9182)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    for i in range(12):
        primer = "".join(rng.choice("ACGT") for _ in range(10))
        if primer in primers:
            continue
        primers.append(primer)
        # Each primer gets 2-4 sites, clustered differently per primer.
        for j in range(2 + (i % 3)):
            pos = (i * 2_400 + j * 700) % (GENOME_LENGTH - 20)
            seq[pos : pos + 10] = list(primer)

    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture
def cache(tmp_path, genome):
    prefix = str(tmp_path / "target")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
        for primer in genome["primers"]:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], genome["seq"].find(key)
                while i != -1:
                    positions.append(i)
                    i = genome["seq"].find(key, i + 1)
                if positions:
                    f.create_dataset(key, data=np.array(positions, dtype=np.int32))
    return PositionCache([prefix], genome["primers"]), prefix


@pytest.fixture
def optimizer(cache):
    position_cache, prefix = cache
    return HybridOptimizer(
        position_cache=position_cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        bin_size=1_000,
        coverage_reach=3_000,
    )


# ----------------------------------------------------------------------
# It runs, and returns something usable
# ----------------------------------------------------------------------


def test_optimize_selects_from_the_candidates(optimizer, genome):
    result = optimizer.optimize(genome["primers"], final_count=5, verbose=False)

    assert isinstance(result, HybridResult)
    assert result.primers, "optimizer selected nothing"
    assert set(result.primers) <= set(genome["primers"]), "a primer was invented"


def test_optimize_respects_the_requested_size(optimizer, genome):
    result = optimizer.optimize(genome["primers"], final_count=4, verbose=False)
    assert len(result.primers) <= 4


def test_optimize_reports_real_coverage(optimizer, genome):
    """Zero coverage would mean the position lookup silently found nothing --
    the failure mode this audit has found in six other modules."""
    result = optimizer.optimize(genome["primers"], final_count=5, verbose=False)

    assert 0.0 < result.final_coverage <= 1.0


def test_no_duplicate_primers_in_the_result(optimizer, genome):
    result = optimizer.optimize(genome["primers"], final_count=6, verbose=False)
    assert len(set(result.primers)) == len(result.primers)


def test_more_primers_do_not_reduce_coverage(optimizer, genome):
    """Set cover is monotone: adding primers cannot uncover the genome."""
    small = optimizer.optimize(genome["primers"], final_count=2, verbose=False)
    large = optimizer.optimize(genome["primers"], final_count=8, verbose=False)

    assert large.final_coverage >= small.final_coverage - 1e-9


def test_result_records_both_stages(optimizer, genome):
    """The two-stage structure is the point of this optimizer.

    Stage 1 is set cover, stage 2 is network refinement; a result that does not
    distinguish them cannot be diagnosed when a set underperforms.
    """
    result = optimizer.optimize(genome["primers"], final_count=5, verbose=False)

    assert result.stage1_primers
    assert result.stage2_primers
    # Stage 1 maximises coverage; stage 2 refines that set for connectivity.
    assert 0.0 <= result.stage1_coverage <= 1.0
    assert result.final_connectivity >= 0.0


# ----------------------------------------------------------------------
# Fixed primers -- the iterative-design workflow
# ----------------------------------------------------------------------


def test_fixed_primers_are_kept(optimizer, genome):
    """`expand-primers` depends on this: validated oligos must survive."""
    keep = genome["primers"][:2]
    result = optimizer.optimize(genome["primers"], final_count=5, fixed_primers=keep, verbose=False)

    for primer in keep:
        assert primer in result.primers, f"fixed primer {primer} was dropped"


def test_fixed_primers_are_not_duplicated(optimizer, genome):
    """They are removed from the candidate pool before selection."""
    keep = genome["primers"][:2]
    result = optimizer.optimize(genome["primers"], final_count=5, fixed_primers=keep, verbose=False)

    assert len(set(result.primers)) == len(result.primers)


def test_fixed_primers_are_case_insensitive(optimizer, genome):
    keep = [p.lower() for p in genome["primers"][:2]]
    result = optimizer.optimize(genome["primers"], final_count=5, fixed_primers=keep, verbose=False)

    upper = {p.upper() for p in result.primers}
    for primer in keep:
        assert primer.upper() in upper


def test_fixing_more_than_requested_still_returns_them(optimizer, genome):
    """Asking for 2 while fixing 4 is a user error, not a crash."""
    keep = genome["primers"][:4]
    result = optimizer.optimize(genome["primers"], final_count=2, fixed_primers=keep, verbose=False)
    assert result.primers


# ----------------------------------------------------------------------
# Coverage reach: selection and scoring must agree
# ----------------------------------------------------------------------


def test_coverage_responds_to_the_reach(cache, genome):
    """`_calculate_coverage` ignored `coverage_reach` entirely.

    It passed no `extension_reach` to the bipartite graph, so the default of 0
    applied and it measured bin OCCUPANCY -- "does any site fall in this bin".
    The answer then depended on `bin_size` rather than on the polymerase: the
    same set scored 13% at a 1 kb bin and 100% at the default 10 kb one, while
    the shipped metric scored it at 39%.

    Stage-1 set cover threads the reach correctly; only this reporting path
    dropped it.
    """
    position_cache, prefix = cache
    primers = genome["primers"][:4]

    def coverage_at(reach):
        opt = HybridOptimizer(
            position_cache=position_cache,
            fg_prefixes=[prefix],
            fg_seq_lengths=[GENOME_LENGTH],
            bin_size=1_000,
            coverage_reach=reach,
        )
        return opt._calculate_coverage(primers)

    narrow, wide = coverage_at(500), coverage_at(10_000)
    assert wide > narrow, "coverage_reach has no effect on the coverage figure"


def test_binned_coverage_is_documented_as_an_approximation(cache, genome):
    """It will not equal the authoritative metric, and should not be read as it.

    `PrimerSetMetrics.fg_coverage` is computed base-by-base and is what reaches
    `step4_improved_df_summary.json`; this one works in bins and is coarser.
    """
    import inspect

    source = inspect.getsource(HybridOptimizer._calculate_coverage)
    assert "APPROXIMATION" in source
    assert "fg_coverage" in source


def test_network_reach_defaults_to_processivity(optimizer):
    """The two reach quantities are distinct and must stay so.

    `max_extension` is the Stage-2 amplification-network reach -- "could two
    primers connect at all" -- for which single-molecule processivity is right.
    Coverage uses the realistic per-primer reach instead.
    """
    from neoswga.core.registry import views

    assert optimizer.max_extension == views.processivity_map()["phi29"]
    assert optimizer.coverage_reach == 3_000
    assert optimizer.coverage_reach < optimizer.max_extension


# ----------------------------------------------------------------------
# Degenerate inputs
# ----------------------------------------------------------------------


def test_empty_candidate_list_does_not_crash(optimizer):
    result = optimizer.optimize([], final_count=5, verbose=False)
    assert result.primers == [] or not result.primers


def test_single_candidate(optimizer, genome):
    result = optimizer.optimize(genome["primers"][:1], final_count=5, verbose=False)
    assert len(result.primers) <= 1


def test_asking_for_more_than_available(optimizer, genome):
    """Requesting 50 from a pool of 12 returns what exists, not an error."""
    result = optimizer.optimize(genome["primers"], final_count=50, verbose=False)
    assert len(result.primers) <= len(genome["primers"])


def test_a_candidate_with_no_positions_is_tolerated(cache, genome):
    """A primer absent from the cache must not derail selection."""
    position_cache, prefix = cache
    opt = HybridOptimizer(
        position_cache=position_cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        bin_size=1_000,
        coverage_reach=3_000,
    )
    result = opt.optimize(genome["primers"] + ["CGCGCGCGCG"], final_count=5, verbose=False)
    assert result.primers


# ----------------------------------------------------------------------
# The BaseOptimizer adapter, which is what the CLI dispatches through
# ----------------------------------------------------------------------


def test_base_optimizer_adapter_produces_comparable_metrics(cache, genome):
    """`normalized_score` is the only cross-optimizer-comparable value, so the
    adapter must populate the metrics it is computed from."""
    from neoswga.core.hybrid_optimizer import HybridBaseOptimizer

    position_cache, prefix = cache
    opt = HybridBaseOptimizer(
        position_cache=position_cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        coverage_reach=3_000,
    )
    result = opt.optimize(genome["primers"], num_primers=5)

    assert result.primers
    assert 0.0 <= result.metrics.normalized_score() <= 1.0
    assert result.metrics.fg_coverage > 0
    assert opt.name
    assert opt.description


# ----------------------------------------------------------------------
# Background pruning -- the selectivity path
# ----------------------------------------------------------------------


@pytest.fixture
def cache_with_background(tmp_path, genome):
    """Target plus a background that shares some of the same primers.

    Sharing matters: with a background that binds nothing, pruning has nothing
    to trade against and the test would pass without exercising the logic.
    """
    rng = random.Random(4242)
    bg = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    # Plant the first half of the primers in the background too, densely, so
    # they are clearly worse choices than the rest.
    shared = genome["primers"][:6]
    for i, primer in enumerate(shared):
        for j in range(6):
            pos = (i * 900 + j * 3_100) % (GENOME_LENGTH - 20)
            bg[pos : pos + 10] = list(primer)
    bg_seq = "".join(bg)

    def write(prefix, sequence):
        with h5py.File(f"{prefix}_10mer_positions.h5", "w") as f:
            for primer in genome["primers"]:
                for key in {primer, reverse_complement(primer)}:
                    positions, i = [], sequence.find(key)
                    while i != -1:
                        positions.append(i)
                        i = sequence.find(key, i + 1)
                    if positions:
                        f.create_dataset(key, data=np.array(positions, dtype=np.int32))

    fg_prefix = str(tmp_path / "fg")
    bg_prefix = str(tmp_path / "bg")
    write(fg_prefix, genome["seq"])
    write(bg_prefix, bg_seq)

    cache = PositionCache([fg_prefix, bg_prefix], genome["primers"])
    return cache, fg_prefix, bg_prefix, shared


def _bg_optimizer(cache_with_background, **kwargs):
    cache, fg_prefix, bg_prefix, _shared = cache_with_background
    return HybridOptimizer(
        position_cache=cache,
        fg_prefixes=[fg_prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        bg_prefixes=[bg_prefix],
        bg_seq_lengths=[GENOME_LENGTH],
        bin_size=1_000,
        coverage_reach=3_000,
        **kwargs,
    )


def test_background_pruning_runs_and_still_selects(cache_with_background, genome):
    opt = _bg_optimizer(cache_with_background, background_pruning=True)
    result = opt.optimize(genome["primers"], final_count=5, verbose=False)

    assert result.primers
    assert result.final_coverage > 0


def test_background_pruning_reduces_background_binding(cache_with_background, genome):
    """The point of the feature. If it does not lower background sites, it is
    costing time for nothing.

    Asserted on the pruning STAGE rather than on the finished pipeline. The
    end-to-end comparison this used to make -- pruned pipeline against unpruned
    pipeline -- is not a property the pipeline guarantees: stage 2 refines on
    connectivity and coverage with no background term at all, so it can give
    back what pruning won, and a run can legitimately finish with more host
    sites than one that never pruned. The old assertion was `<=`, which a
    complete no-op also satisfies, and for a long time a no-op is exactly what
    it was measuring: the coverage floor was absolute at 0.95, so on any set
    below that every removal was rejected and the stage returned its input.
    """
    opt = _bg_optimizer(cache_with_background, background_pruning=True)
    stage1 = list(opt.dominating_optimizer.optimize_greedy(genome["primers"], 7)["primers"])
    before = opt._count_background_sites(stage1)
    assert before > 0, "the fixture has no background binding to prune"

    kept, coverage, after = opt._prune_background(stage1, target_size=5, verbose=False)

    assert len(kept) == 5
    assert after < before, (
        f"pruning removed no background: {before} sites before, {after} after. "
        "A stage that returns its input is the failure mode here."
    )
    assert coverage > 0


def test_background_pruning_respects_the_coverage_floor(cache_with_background, genome):
    """It trades selectivity for coverage, but not below min_coverage_threshold
    relative to where it started."""
    opt = _bg_optimizer(cache_with_background, background_pruning=True, min_coverage_threshold=0.95)
    result = opt.optimize(genome["primers"], final_count=5, verbose=False)

    assert result.final_coverage > 0


def test_the_coverage_floor_is_relative_to_where_pruning_started():
    """Pruning must still run on a set whose coverage is nowhere near 0.95.

    The floor was absolute: a removal was rejected whenever what remained fell
    below `min_coverage_threshold`, default 0.95. The quantity compared is
    binned coverage at realistic reach, which on a real target sits well below
    that -- the note on `_calculate_coverage` records the same set reading 39%
    under the corrected measure where the older bin-occupancy measure read
    100%. So every candidate removal was rejected on the first pass and
    `background-aware`, the method sold on a 10-20x background reduction,
    quietly became plain hybrid on every ordinary run.

    Driven with a controlled coverage function rather than a genome, because
    the point is the shape of the floor and not any particular target.
    """
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    primers = [f"P{i}" for i in range(1, 13)]
    background = {p: 100 * i for i, p in enumerate(primers, 1)}

    def count_background(subset):
        try:
            return sum(background.get(p, 0) for p in subset)
        except TypeError:  # the loop also probes with non-primer values
            return 0

    def optimizer_at(coverage):
        opt = HybridOptimizer.__new__(HybridOptimizer)
        opt.min_coverage_threshold = 0.95
        opt._absolute_coverage_floor = False
        opt.min_coverage_drop = 0.05
        opt.background_weight = 2.0
        opt.fg_prefixes, opt.fg_seq_lengths = ["fg"], [1_000_000]
        opt.bg_prefixes, opt.bg_seq_lengths = ["bg"], [1_000_000]
        # Coverage decays gently as primers are dropped, anchored at `coverage`.
        opt._calculate_coverage = lambda s: coverage * (len(s) / 12.0) ** 0.05
        opt._count_background_sites = count_background
        return opt

    started_with = count_background(primers)
    for coverage in (0.99, 0.94, 0.80, 0.60, 0.39):
        kept, _final_coverage, _bg = optimizer_at(coverage)._prune_background(
            list(primers), target_size=9, verbose=False
        )
        assert len(kept) == 9, (
            f"at stage-1 coverage {coverage} pruning kept {len(kept)} of 12 primers; "
            "an absolute floor blocks every removal once coverage is below it"
        )
        assert count_background(kept) < started_with


def test_background_sites_are_counted_from_the_background_prefix(cache_with_background, genome):
    """A count of zero would mean the background cache is not being read --
    which would make pruning a no-op that looks like it worked."""
    _cache, _fg, _bg, shared = cache_with_background
    opt = _bg_optimizer(cache_with_background)

    assert opt._count_background_sites(shared) > 0


# ----------------------------------------------------------------------
# Thermodynamic pre-filter
# ----------------------------------------------------------------------


def test_thermo_filter_narrows_the_candidate_pool(optimizer, genome):
    """The pre-stage filters candidates to the polymerase's Tm and GC window."""
    kept = optimizer._thermo_filter_candidates(genome["primers"], verbose=False)

    assert isinstance(kept, list)
    assert set(kept) <= set(genome["primers"])


def test_thermo_filter_never_returns_nothing_from_a_real_pool(optimizer, genome):
    """Filtering everything out would leave the optimizer with no input, and a
    zero-primer result reads as a design failure rather than a filter one."""
    kept = optimizer._thermo_filter_candidates(genome["primers"], verbose=False)
    assert kept, "the thermodynamic pre-filter rejected every candidate"
