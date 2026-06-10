"""Tests for BAM coverage-gap detection and gap-targeted expansion."""

import numpy as np
import pytest

from neoswga.core.bam_coverage import (
    find_low_depth_gaps,
    match_contigs,
    bam_gaps,
)
from neoswga.core.primer_expansion import (
    CoverageGap,
    merge_gap_intervals,
    PrimerExpander,
)


# ---------------------------------------------------------------------------
# find_low_depth_gaps
# ---------------------------------------------------------------------------

def test_find_low_depth_gaps_linear():
    # depth: high everywhere except a 0-depth window [30, 60)
    depth = np.full(100, 10, dtype=np.int32)
    depth[30:60] = 0
    gaps = find_low_depth_gaps(depth, "chr", min_depth=5, min_gap_size=10)
    assert len(gaps) == 1
    g = gaps[0]
    assert (g.chromosome, g.start, g.end, g.size) == ("chr", 30, 60, 30)


def test_find_low_depth_gaps_respects_min_gap_size():
    depth = np.full(100, 10, dtype=np.int32)
    depth[30:35] = 0   # only 5 bp low
    gaps = find_low_depth_gaps(depth, "chr", min_depth=5, min_gap_size=10)
    assert gaps == []


def test_find_low_depth_gaps_circular_wrap_merges_origin():
    # low at the very end and very start -> one origin-spanning gap.
    depth = np.full(100, 10, dtype=np.int32)
    depth[90:100] = 0
    depth[0:10] = 0
    gaps = find_low_depth_gaps(depth, "chr", min_depth=5, min_gap_size=5,
                               circular=True)
    assert len(gaps) == 1
    g = gaps[0]
    assert g.start == 90 and g.end == 110 and g.size == 20  # wrap encoded


def test_find_low_depth_gaps_no_low_returns_empty():
    depth = np.full(50, 8, dtype=np.int32)
    assert find_low_depth_gaps(depth, "chr", min_depth=5, min_gap_size=1) == []


# ---------------------------------------------------------------------------
# match_contigs
# ---------------------------------------------------------------------------

def test_match_contigs_exact_and_basename():
    m = match_contigs(["genomeA"], [1000], ["data/genomeA"], [1000])
    assert m == {"data/genomeA": "genomeA"}


def test_match_contigs_chr_prefix():
    m = match_contigs(["I"], [1000], ["chrI"], [1000])
    assert m == {"chrI": "I"}


def test_match_contigs_length_fallback():
    m = match_contigs(["weird_name"], [12345], ["fg"], [12345])
    assert m == {"fg": "weird_name"}


def test_match_contigs_alias_and_unmatched():
    m = match_contigs(
        ["bamctg"], [999], ["fg"], [1000],   # length differs -> no fallback
        aliases={"fg": "bamctg"},
    )
    assert m == {"fg": "bamctg"}
    # No alias, no match -> empty
    m2 = match_contigs(["bamctg"], [999], ["fg"], [1000])
    assert m2 == {}


# ---------------------------------------------------------------------------
# merge_gap_intervals
# ---------------------------------------------------------------------------

def test_merge_gap_intervals_overlap_adjacent_disjoint():
    gaps = [
        CoverageGap("c", 0, 100, 100),
        CoverageGap("c", 90, 150, 60),    # overlaps previous
        CoverageGap("c", 150, 200, 50),   # touches previous
        CoverageGap("c", 500, 600, 100),  # disjoint
        CoverageGap("d", 10, 20, 10),     # other chromosome
    ]
    merged = merge_gap_intervals(gaps)
    intervals = sorted((g.chromosome, g.start, g.end) for g in merged)
    assert intervals == [("c", 0, 200), ("c", 500, 600), ("d", 10, 20)]


# ---------------------------------------------------------------------------
# Gap-restricted candidate pre-filter
# ---------------------------------------------------------------------------

class _Cache:
    def __init__(self, mapping):
        self._m = mapping  # (prefix, primer) -> np.array

    def get_positions(self, prefix, primer, strand="both"):
        return self._m.get((prefix, primer), np.array([], dtype=np.int64))


def test_filter_candidates_to_gaps_keeps_only_in_gap():
    cache = _Cache({
        ("fg", "IN"): np.array([55]),       # inside gap [30,60)
        ("fg", "OUT"): np.array([5, 80]),   # outside the gap
    })
    expander = PrimerExpander(cache, ["fg"], [100])
    gaps = [CoverageGap("fg", 30, 60, 30)]
    kept = expander._filter_candidates_to_gaps(["IN", "OUT"], gaps)
    assert kept == ["IN"]


def test_filter_candidates_to_gaps_handles_wrap_gap():
    # wrap gap encoded as end>length: covers [90,100) and [0,10)
    cache = _Cache({
        ("fg", "WRAP"): np.array([5]),      # inside [0,10)
        ("fg", "MID"): np.array([50]),      # not in wrap
    })
    expander = PrimerExpander(cache, ["fg"], [100])
    gaps = [CoverageGap("fg", 90, 110, 20)]
    kept = expander._filter_candidates_to_gaps(["WRAP", "MID"], gaps)
    assert kept == ["WRAP"]


# ---------------------------------------------------------------------------
# identify_gaps merges extra (BAM) gaps
# ---------------------------------------------------------------------------

def test_identify_gaps_merges_extra_gaps():
    # No primers -> whole genome is a gap; merging extra is a no-op union.
    cache = _Cache({})
    expander = PrimerExpander(cache, ["fg"], [1000])
    extra = [CoverageGap("fg", 200, 400, 200)]
    gaps = expander.identify_gaps([], min_gap_size=10, extra_gaps=extra)
    # whole-genome [0,1000) already subsumes [200,400)
    assert len(gaps) == 1
    assert gaps[0].start == 0 and gaps[0].end == 1000


# ---------------------------------------------------------------------------
# pysam-gated end-to-end
# ---------------------------------------------------------------------------

def test_require_pysam_friendly_error(monkeypatch):
    import builtins
    import neoswga.core.bam_coverage as bc

    real_import = builtins.__import__

    def fake_import(name, *a, **k):
        if name == "pysam":
            raise ImportError("no pysam")
        return real_import(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(RuntimeError, match=r"neoswga\[bam\]"):
        bc._require_pysam()


def test_bam_gaps_end_to_end(tmp_path):
    pysam = pytest.importorskip("pysam")

    # Build a tiny reference + BAM with reads only over part of the contig,
    # leaving a low-depth gap in the middle.
    length = 300
    ref_name = "fg"
    header = {"HD": {"VN": "1.0"},
              "SQ": [{"SN": ref_name, "LN": length}]}
    bam_path = tmp_path / "reads.bam"
    seq = "A" * 50
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out:
        # Cover [0,100) and [200,300); leave [100,200) empty.
        for start in list(range(0, 100, 10)) + list(range(200, 300, 10)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{start}"
            a.query_sequence = seq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigar = [(0, 50)]  # 50M
            out.write(a)
    pysam.index(str(bam_path))

    gaps = bam_gaps(str(bam_path), [ref_name], [length],
                    min_depth=1, min_gap_size=20)
    # Reads start at 0..90 (50bp) so coverage reaches base 140; next read at
    # 200. The uncovered window is therefore [140, 200).
    assert len(gaps) == 1
    g = gaps[0]
    assert g.chromosome == ref_name
    assert 130 <= g.start <= 145
    assert 195 <= g.end <= 205
