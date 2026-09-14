"""An extension window must not reach out of one record and into the next.

Audit finding F3, second half. The first half stopped a k-mer being MATCHED
across the join between two FASTA records. This stops a coverage window
anchored near the end of one record from marking bases in the next, which is a
different error with the same cause: the concatenation is the coordinate system
and nothing in it knew where one record ended.

Measured on the shipped Wolbachia/Drosophila pair the effect is small -- 3 of
the delivered panel's 148 background sites sit within one reach of a record
edge, and the overcount is 11,289 bases, 0.008% of the assembly -- because that
background is six chromosome-scale records holding 94.8% of the sequence.

The exposure is in draft assemblies, which SWGA is routinely pointed at. 92.3%
of that same reference's 1870 records are shorter than twice the 3 kb reach, and
on an assembly made only of contigs that size every window would cross an edge
and coverage would be inflated throughout.
"""

import numpy as np
import pytest

from neoswga.core.coverage import _mark_window

LENGTH = 1_000
# Three records: [0, 400), [400, 700), [700, 1000).
RECORD_STARTS = [0, 400, 700]


def _marked(pos, extension, circular=False, record_starts=RECORD_STARTS):
    occupied = np.zeros(LENGTH, dtype=bool)
    _mark_window(occupied, pos, extension, LENGTH, circular, record_starts=record_starts)
    return occupied


def test_a_window_is_clipped_to_the_record_holding_its_site():
    """A site at 390 with a 100 bp reach must not touch the record at 400."""
    occupied = _marked(390, 100)

    assert occupied[300:400].all()
    assert not occupied[400:].any(), "the window reached into the next record"


def test_a_window_does_not_reach_backwards_into_the_previous_record():
    occupied = _marked(410, 100)

    assert occupied[400:510].all()
    assert not occupied[:400].any()


def test_a_site_in_the_interior_is_unaffected():
    occupied = _marked(550, 100)

    assert occupied[450:650].all()
    assert occupied.sum() == 200


def test_without_record_starts_the_behaviour_is_unchanged():
    """Single-record references, and every existing caller, must not move."""
    occupied = _marked(390, 100, record_starts=None)

    assert occupied[290:490].all()
    assert occupied.sum() == 200


def test_a_circular_window_wraps_within_its_own_record_only():
    """A multi-record file is not one circle.

    Wrapping a window from the end of the last record to the start of the first
    would join sequences that are not adjacent, which is the same error the
    clipping exists to prevent.
    """
    occupied = _marked(950, 100, circular=True)

    assert occupied[850:1000].all()
    assert not occupied[:700].any(), "wrapped into an unrelated record"


@pytest.mark.parametrize("pos", [0, 399, 400, 699, 700, 999])
def test_every_site_marks_only_its_own_record(pos):
    occupied = _marked(pos, 5_000)
    lo = max(s for s in RECORD_STARTS if s <= pos)
    hi = min([s for s in RECORD_STARTS if s > pos] + [LENGTH])

    assert occupied[lo:hi].all()
    assert occupied.sum() == hi - lo


def test_the_clipping_is_actually_wired_to_the_position_index(tmp_path):
    """The whole path: scan a two-record FASTA, then measure coverage.

    A capability nothing calls is the failure this project already recorded
    once, when a background genome was read and never acted on. So this drives
    the real chain -- scan writes record starts into the index, `PositionCache`
    reads them, `compute_per_prefix_coverage` passes them to `_mark_window` --
    rather than calling the helper directly as the tests above do.
    """
    h5py = pytest.importorskip("h5py")

    from neoswga.core import string_search
    from neoswga.core.coverage import compute_per_prefix_coverage
    from neoswga.core.position_cache import PositionCache

    # Two 60 bp records; the primer sits at the very end of the first.
    primer = "ACGTACGTACGT"
    first = "T" * 48 + primer
    second = "G" * 60
    fasta = tmp_path / "two.fna"
    fasta.write_text(f">one\n{first}\n>two\n{second}\n")
    string_search.clear_genome_cache()

    prefix = str(tmp_path / "idx")
    found = string_search.get_all_positions_multi_k({12: [primer]}, str(fasta), circular=False)
    with h5py.File(f"{prefix}_12mer_positions.h5", "w"):
        pass
    string_search.write_to_h5py(
        {primer: found[primer]},
        prefix,
        record_starts=[0] + string_search.get_cached_record_boundaries(str(fasta)),
    )

    cache = PositionCache([prefix], [primer])
    assert cache.get_record_starts(prefix) == [0, 60], "record starts did not survive the index"

    aggregate, _ = compute_per_prefix_coverage(
        cache=cache,
        primers=[primer],
        prefixes=[prefix],
        seq_lengths=[120],
        extension=30,
        circular=False,
    )

    # The site is at 48. A 30 bp reach would run to 78, twelve bases into the
    # second record; confined to record one it stops at 60.
    assert aggregate == pytest.approx((60 - 18) / 120)
