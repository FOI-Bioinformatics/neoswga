"""The background scan runs after the cut, over the pool `candidate_retention` keeps.

Audit finding B1: `step2` scanned background positions for every candidate that
passed the frequency and sequence-quality gates, and cut to `max_primer` only
afterwards. Nothing between the two points needs background positions --
`filter.get_gini` reads the foreground only, and `_rank_and_cut_candidates`
reaches `occupancy.weighted_site_load`, which reads the jellyfish count tables.
So the scan moved behind the cut with no approximation. That ORDERING is what
`test_the_background_scan_runs_after_the_cut` still pins.

**The SIZE it is handed is a separate contract, and it changed on 2026-09-16.**
This file used to assert the scan received exactly the `max_primer` shortlist.
`candidate_retention` deliberately replaced that: indexing only the shortlist
left every other hard-QC candidate with no background index, so a design
reaching one scored it against an empty background and read as perfectly
specific -- the silent-zero shape of Known Issues 5, 6 and 13. The removed
`legacy` mode is exactly the behaviour the old assertion demanded.

So the size is now whatever `candidate_inventory.background_scan_pool` returns
for the configured mode, and both modes are pinned below. Measured on this
fixture with `max_primer` 100: `all_qc` scans 5,190, `post_gini` scans 37, and
the shortlist is 37 either way.

The old assertion failed only once `examples/plasmid_example` was primed, which
needs jellyfish, so it was invisible on an unprimed checkout and on CI.

Observed by spying on `string_search.get_positions` and recording how many
primers each call was handed. A timing assertion would not separate the two
cases on a fixture small enough to run in the unit suite.
"""

import json
import os
import shutil
from pathlib import Path

import pytest

import neoswga.core.pipeline as pipeline
from neoswga.core import parameter, string_search
from tests.conftest import plasmid_example_ready

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = ROOT / "examples" / "plasmid_example"


def _reset_pipeline_state(params_file):
    pipeline._initialized = False
    pipeline.fg_prefixes = None
    pipeline.bg_prefixes = None
    pipeline.fg_genomes = None
    pipeline.bg_genomes = None
    pipeline.fg_seq_lengths = None
    pipeline.bg_seq_lengths = None
    pipeline.fg_circular = None
    pipeline.bg_circular = None
    parameter.json_file = str(params_file)


@pytest.fixture
def plasmid(tmp_path, monkeypatch, request):
    """The plasmid example, copied out and cut hard enough to be observable.

    Indirectly parameterisable with a `candidate_retention` mode; the default
    is whatever params.json carries, which is what a plain run uses.
    """
    if not plasmid_example_ready():
        pytest.skip("plasmid example not available")

    for name in os.listdir(EXAMPLE_DIR):
        src = EXAMPLE_DIR / name
        if src.is_file():
            shutil.copy2(src, tmp_path / name)

    # Scan from scratch. The checked-in position files would otherwise be
    # reused, which would leave this test measuring nothing.
    for stale in tmp_path.glob("*_positions.h5"):
        stale.unlink()

    params_path = tmp_path / "params.json"
    params = json.loads(params_path.read_text())
    params["max_primer"] = 100
    retention = getattr(request, "param", None)
    if retention is not None:
        params["candidate_retention"] = retention
    params_path.write_text(json.dumps(params, indent=2))

    before = getattr(parameter, "json_file", None)
    monkeypatch.chdir(tmp_path)
    _reset_pipeline_state(params_path)
    string_search.clear_genome_cache()
    yield tmp_path
    _reset_pipeline_state(params_path)
    parameter.json_file = before
    pipeline._initialized = False
    string_search.clear_genome_cache()


def _spy_on_scans(monkeypatch):
    """Record the prefixes and pool size of every position scan, in order."""
    calls = []
    real = string_search.get_positions

    def spy(primer_list, fname_prefixes, fname_genomes, circular, **kwargs):
        primers = list(primer_list)
        calls.append((tuple(fname_prefixes), len(primers)))
        return real(primers, fname_prefixes, fname_genomes, circular, **kwargs)

    monkeypatch.setattr(string_search, "get_positions", spy)
    return calls


def test_the_background_scan_runs_after_the_cut(plasmid, monkeypatch):
    """Audit finding B1's ordering property, which is independent of size.

    One foreground scan, then one background scan. If the background scan ran
    in the same pass as the foreground one it could not be handed a different
    pool at all, whatever `candidate_retention` says.
    """
    calls = _spy_on_scans(monkeypatch)

    pipeline.step2()

    assert [prefixes for prefixes, _ in calls] == [
        ("pcDNA",),
        ("pLTR",),
    ], f"expected one foreground scan then one background scan, got {calls}"


@pytest.mark.parametrize(
    "plasmid,expected",
    [("all_qc", "every hard-QC candidate"), ("post_gini", "the post-Gini pool")],
    indirect=["plasmid"],
)
def test_the_background_scan_gets_the_pool_retention_asks_for(plasmid, monkeypatch, expected):
    """The size contract, which `candidate_retention` owns.

    Under `all_qc` the background scan must receive exactly what the foreground
    scan received: every candidate clearing the hard gates. Anything smaller
    leaves a retained candidate with no background index, which reads as
    perfect specificity rather than as a missing measurement.

    Under `post_gini` it must receive fewer, because the evenness gate is an
    admission rule there.
    """
    calls = _spy_on_scans(monkeypatch)

    result = pipeline.step2()

    foreground_count, background_count = calls[0][1], calls[1][1]
    retention = json.loads((plasmid / "params.json").read_text())["candidate_retention"]

    if retention == "all_qc":
        assert background_count == foreground_count, (
            f"{expected} must be indexed, but the background scan got "
            f"{background_count} of {foreground_count}. A retained candidate "
            "with no background index scores as perfectly specific."
        )
    else:
        assert background_count < foreground_count, (
            f"{expected} must be smaller than the hard-QC pool, but the "
            f"background scan got {background_count} of {foreground_count}"
        )
    assert background_count >= len(result), (
        f"the background scan got {background_count} primers, fewer than the "
        f"{len(result)} the shortlist delivers, so a delivered primer has no "
        "background index at all"
    )


def test_the_foreground_is_still_scanned_before_the_gini_gate(plasmid, monkeypatch):
    """Gini reads foreground positions, so that scan must not move."""
    calls = []
    real = string_search.get_positions

    def spy(primer_list, fname_prefixes, fname_genomes, circular, **kwargs):
        primers = list(primer_list)
        calls.append((tuple(fname_prefixes), len(primers)))
        return real(primers, fname_prefixes, fname_genomes, circular, **kwargs)

    monkeypatch.setattr(string_search, "get_positions", spy)

    result = pipeline.step2()

    assert calls[0][0] == ("pcDNA",)
    assert calls[0][1] > len(result), (
        "the foreground scan was handed the cut pool; the Gini gate runs "
        "before the cut and would have no positions to measure"
    )
    assert result["gini"].notna().all(), "Gini came back unmeasurable"


def test_the_genome_cache_is_released_when_the_scans_are_done(plasmid):
    """Nothing past step 2's scans reads the genome strings.

    `_genome_cache` is a module-level dict, so without a release the foreground
    and background genomes stay resident for the rest of the process. For hg38
    that is 3.3 GB held alongside the 1.24 GB k-mer count table
    `mismatch_counts` caches for the occupancy ranking.

    `clear_genome_cache` has existed since the cache did and had no caller
    outside the tests.
    """
    assert string_search.get_genome_cache_stats()["num_genomes"] == 0

    pipeline.step2()

    stats = string_search.get_genome_cache_stats()
    assert stats["num_genomes"] == 0, (
        f"step 2 left {stats['num_genomes']} genome(s) and {stats['total_bp']:,} bp "
        "resident after both scans finished"
    )


def test_the_release_happens_after_the_background_scan(plasmid, monkeypatch):
    """Releasing too early would force the background genome to be re-parsed."""
    seen = []
    real_scan = string_search.get_positions
    real_clear = string_search.clear_genome_cache

    def scan_spy(primer_list, fname_prefixes, fname_genomes, circular, **kwargs):
        seen.append(("scan", tuple(fname_prefixes)))
        return real_scan(list(primer_list), fname_prefixes, fname_genomes, circular, **kwargs)

    def clear_spy():
        seen.append(("clear", ()))
        return real_clear()

    monkeypatch.setattr(string_search, "get_positions", scan_spy)
    monkeypatch.setattr(string_search, "clear_genome_cache", clear_spy)

    pipeline.step2()

    assert seen == [
        ("scan", ("pcDNA",)),
        ("scan", ("pLTR",)),
        ("clear", ()),
    ], f"expected both scans then one release, got {seen}"


@pytest.fixture
def plasmid_without_background(plasmid):
    """The same fixture with the background removed.

    Every other test here carries `bg_prefixes: ["pLTR"]`, so a release placed
    inside `if len(bg_prefixes) > 0` passes all of them while leaving the
    foreground genome resident on a run that configures no background --
    phi29_baseline is exactly that case.
    """
    params_path = plasmid / "params.json"
    params = json.loads(params_path.read_text())
    params["bg_genomes"] = []
    params["bg_prefixes"] = []
    params["bg_seq_lengths"] = []
    params_path.write_text(json.dumps(params, indent=2))
    _reset_pipeline_state(params_path)
    string_search.clear_genome_cache()
    return plasmid


def test_the_genome_cache_is_released_when_no_background_is_configured(
    plasmid_without_background,
):
    """The release is not conditional on there being a background to scan."""
    assert string_search.get_genome_cache_stats()["num_genomes"] == 0

    result = pipeline.step2()
    assert len(result) > 0, "the fixture produced no candidates, so it measures nothing"

    stats = string_search.get_genome_cache_stats()
    assert stats["num_genomes"] == 0, (
        f"a run with no background left {stats['num_genomes']} genome(s) and "
        f"{stats['total_bp']:,} bp resident; the release is inside the "
        "background branch"
    )


def test_a_second_step2_reuses_the_position_files(plasmid, monkeypatch):
    """Reuse through the pipeline, not just through `get_positions`.

    Every other reuse test calls `string_search.get_positions` directly. The
    failure the reuse could cause is a step-2 outcome: a primer left out of the
    returned map reaches the Gini gate as a primer that binds nowhere, scores
    NaN and is dropped, which on a second run would empty the whole pool while
    the step logged success.
    """
    first = pipeline.step2()
    assert len(first) > 0

    scanned = []
    real = string_search.get_all_positions_multi_k

    def spy(primer_lists_by_k, seq_fname, circular, chunk_size=None):
        scanned.append(sum(len(v) for v in primer_lists_by_k.values()))
        return real(primer_lists_by_k, seq_fname, circular, chunk_size)

    monkeypatch.setattr(string_search, "get_all_positions_multi_k", spy)
    string_search.clear_genome_cache()
    _reset_pipeline_state(plasmid / "params.json")

    second = pipeline.step2()

    assert scanned == [], (
        "the second step2 rescanned; pattern counts per automaton build were " f"{scanned}"
    )
    assert len(second) == len(first), (
        f"the second run kept {len(second)} candidates where the first kept "
        f"{len(first)}; a short position map empties the pool through the Gini gate"
    )
    assert second["gini"].notna().all(), "Gini came back unmeasurable on the reused run"
    assert list(second["primer"]) == list(first["primer"])
