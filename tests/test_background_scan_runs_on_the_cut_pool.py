"""The background is scanned for the pool that survives, not the pool that does not.

`step2` scanned background positions for every candidate that passed the
frequency and sequence-quality gates, and cut the pool to `max_primer` only
afterwards. On the run that found this, 20,301 primers were scanned and 3,000
kept.

Nothing between the two points needs background positions. `filter.get_gini`
reads the foreground only, and `_rank_and_cut_candidates` reaches
`occupancy.weighted_site_load`, which reads the jellyfish count tables --
neither `occupancy.py` nor `mismatch_counts.py` imports h5py. So the scan can
move behind the cut with no approximation.

Aho-Corasick cost grows with the pattern count, so on a whole-genome host this
is the difference between about 15 minutes and about 1 minute, and it accounts
for very nearly all of the 924 s the project has recorded for a filter step
against a large host.

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
def plasmid(tmp_path, monkeypatch):
    """The plasmid example, copied out and cut hard enough to be observable."""
    if not EXAMPLE_DIR.is_dir():
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


def test_the_background_is_scanned_for_the_cut_pool(plasmid, monkeypatch):
    """The background call must receive the kept pool, not the pre-cut pool."""
    calls = []
    real = string_search.get_positions

    def spy(primer_list, fname_prefixes, fname_genomes, circular, **kwargs):
        primers = list(primer_list)
        calls.append((tuple(fname_prefixes), len(primers)))
        return real(primers, fname_prefixes, fname_genomes, circular, **kwargs)

    monkeypatch.setattr(string_search, "get_positions", spy)

    result = pipeline.step2()

    assert [prefixes for prefixes, _ in calls] == [
        ("pcDNA",),
        ("pLTR",),
    ], f"expected one foreground scan then one background scan, got {calls}"
    foreground_count, background_count = calls[0][1], calls[1][1]

    assert background_count == len(result), (
        f"the background scan was handed {background_count} primers but only "
        f"{len(result)} survive the cut; the scan is still running ahead of "
        "_rank_and_cut_candidates"
    )
    assert background_count < foreground_count, (
        "the cut removed nothing, so this fixture cannot tell the two orders "
        f"apart (foreground {foreground_count}, background {background_count})"
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
