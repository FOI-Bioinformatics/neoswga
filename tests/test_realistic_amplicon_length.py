"""Phase 16 critical gap #2 + Phase 16.5 semantics fix: coverage uses the
effective per-primer reach in a dense SWGA design by default, not
polymerase processivity.

Processivity is a single-molecule theoretical maximum (phi29 ~70 kb,
Blanco 1989). In a dense SWGA reaction, extension from any one primer
is truncated by neighbouring primers' strand-displacement products
after ~2-5 kb — Clarke et al. (2017) post-hoc filter M. tuberculosis
sets to <5 kb mean inter-primer-site spacing; Dwivedi-Yu et al. (2023)
PLOS Comput Biol 19:e1010137 report successful Prevotella sets at
1/2-5 kbp site densities. Reporting coverage with a 70 kb window
therefore over-states the reach of a primer set by 10-30x for phi29.

NB: this is NOT gel fragment size (unselective phi29 MDA peaks near
~10 kb per Dean 2002 / Picher 2016); neoswga stores per-primer reach
under `typical_amplicon_length` specifically for the coverage metric.

These tests lock in:
- every registered polymerase exposes a `typical_amplicon_length` field;
- `get_typical_amplicon_length` returns the per-primer-reach values
  (phi29 3 kb, equiphi29 4 kb, bst 1 kb, klenow 1.5 kb);
- `compute_per_prefix_coverage`'s default extension is per-primer reach,
  not processivity;
- `polymerase_extension_reach(..., coverage_metric='processivity')`
  still returns the legacy values for callers that need graph reach.
"""

import pytest


@pytest.mark.parametrize("polymerase,expected", [
    ("phi29", 3000),
    ("equiphi29", 4000),
    ("bst", 1000),
    ("klenow", 1500),
])
def test_typical_amplicon_length_matches_literature(polymerase, expected):
    from neoswga.core.reaction_conditions import (
        POLYMERASE_CHARACTERISTICS, get_typical_amplicon_length,
    )
    assert POLYMERASE_CHARACTERISTICS[polymerase]["typical_amplicon_length"] == expected
    assert get_typical_amplicon_length(polymerase) == expected


def test_typical_amplicon_strictly_shorter_than_processivity():
    """Sanity: per-primer reach < single-event processivity for every
    polymerase."""
    from neoswga.core.reaction_conditions import (
        get_typical_amplicon_length, get_polymerase_processivity,
    )
    for poly in ("phi29", "equiphi29", "bst", "klenow"):
        reach = get_typical_amplicon_length(poly)
        processivity = get_polymerase_processivity(poly)
        assert reach <= processivity, (
            f"{poly}: per-primer reach ({reach}) must not exceed "
            f"single-event processivity ({processivity})"
        )
        # phi29 / equiphi29: reach is spacing-bounded (Clarke 2017,
        # Dwivedi-Yu 2023) so well below processivity. bst / klenow:
        # reach is already processivity-limited.
        if poly in ("phi29", "equiphi29"):
            assert reach * 5 <= processivity, (
                f"{poly}: per-primer reach ({reach}) should be at most "
                f"1/5 of processivity ({processivity}) in a dense SWGA "
                f"design"
            )


def test_polymerase_extension_reach_defaults_to_realistic():
    """Phase 16: the default is per-primer reach, not processivity."""
    from neoswga.core.coverage import polymerase_extension_reach
    assert polymerase_extension_reach("phi29") == 3000
    assert polymerase_extension_reach("equiphi29") == 4000


def test_polymerase_extension_reach_processivity_opt_in():
    """Callers that explicitly ask for processivity still get it."""
    from neoswga.core.coverage import polymerase_extension_reach
    assert polymerase_extension_reach("phi29", coverage_metric="processivity") == 70000
    assert polymerase_extension_reach("equiphi29", coverage_metric="processivity") == 80000
    assert polymerase_extension_reach("bst", coverage_metric="processivity") == 2000


def test_compute_per_prefix_coverage_default_extension_is_realistic():
    """Directly exercise the default: with a single primer at mid-genome
    and the new 3 kb default, the covered window is 6 kb / 20 kb = 30%
    instead of the old 40 kb (clamped at genome length = 100%)."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    class FakeCache:
        def get_positions(self, prefix, primer, strand):
            return [10000]  # single mid-genome position

    agg, per = compute_per_prefix_coverage(
        cache=FakeCache(), primers=["ACGT"],
        prefixes=["x"], seq_lengths=[20000],
        # NB: no `extension=` — exercising the default
    )
    # Default extension = 3000 → window [7000, 13000) = 6000 bp covered
    assert abs(agg - 0.30) < 1e-9, (
        f"Default per-primer-reach coverage should yield 30% on 20 kb "
        f"genome with one mid-position primer; got {agg:.3f}. If this "
        f"fails, the default likely changed from 3000 bp."
    )


def test_compute_per_prefix_coverage_explicit_extension_honoured():
    """Callers that explicitly pass extension override the default."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    class FakeCache:
        def get_positions(self, prefix, primer, strand):
            return [10000]

    # 70 kb extension on a 20 kb genome → full coverage
    agg, _ = compute_per_prefix_coverage(
        cache=FakeCache(), primers=["ACGT"],
        prefixes=["x"], seq_lengths=[20000], extension=70000,
    )
    assert agg == 1.0


def test_coverage_change_documented_in_citations():
    """Phase 16/16.5: SCIENCE_CITATIONS.md must surface the per-primer
    reach numbers so scientific reviewers can find them. Accept either
    unformatted (3000) or comma-separated (3,000) rendering."""
    from pathlib import Path
    doc = Path(__file__).resolve().parent.parent / "docs" / "SCIENCE_CITATIONS.md"
    if not doc.is_file():
        pytest.skip("SCIENCE_CITATIONS.md not present")
    text = doc.read_text()
    # The field name must appear verbatim.
    assert "typical_amplicon_length" in text, (
        "SCIENCE_CITATIONS.md should mention 'typical_amplicon_length'"
    )
    # Numbers may be rendered with or without thousands separator.
    for value, comma_form in (("3000", "3,000"), ("4000", "4,000")):
        assert value in text or comma_form in text, (
            f"SCIENCE_CITATIONS.md should document the per-primer reach "
            f"{value!r} (or {comma_form!r})"
        )
