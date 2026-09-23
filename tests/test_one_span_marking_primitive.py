"""`_mark_window` marks a span; it should say so in its arguments.

The body at `coverage.py:262-293` works entirely on `[start, end)`. The site
position is used for exactly one thing -- a `bisect` to decide which record the
site sits in -- and everything else is span arithmetic: clip to the record,
wrap on a circular reference, clip to the genome otherwise.

That is a span primitive with a symmetric window baked into its signature, and
it is why a directional span has nowhere to go. `mark_span` is the same body
with the span passed in and the site kept only as the record anchor;
`_mark_window` becomes the symmetric caller of it.

**Nothing about the symmetric path may change**, and this file's job is to
prove it rather than assert it. `test_coverage_independent_oracle.py` already
checks `_mark_window` against an oracle that shares none of its code, so the
refactor is verifiable against something that does not share its assumptions;
the tests here pin the seam itself, including the two cases the delegation is
most likely to get wrong -- a window that wraps past both ends of a circular
reference, and a window whose site sits in one record while its span would
reach into the next.
"""

import numpy as np
import pytest

from neoswga.core.coverage import _mark_window, mark_span


def marked(fn, size, **kwargs):
    """`size` rather than `length`, which is also one of the marked functions'
    own arguments and would bind twice."""
    occupied = np.zeros(size, dtype=bool)
    fn(occupied, **kwargs)
    return occupied


# ---------------------------------------------------------------------------
# The seam: the symmetric caller and the primitive agree
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("circular", [False, True])
@pytest.mark.parametrize("pos,extension", [(50, 10), (0, 10), (99, 10), (50, 200), (5, 5)])
def test_the_window_is_the_span_it_names(pos, extension, circular):
    """`_mark_window(pos, ext)` must mark exactly `mark_span(pos-ext, pos+ext)`."""
    window = marked(_mark_window, 100, pos=pos, extension=extension, length=100, circular=circular)
    span = marked(
        mark_span,
        100,
        start=pos - extension,
        end=pos + extension,
        length=100,
        circular=circular,
        anchor=pos,
    )

    assert np.array_equal(window, span), f"pos={pos} ext={extension} circular={circular}"


@pytest.mark.parametrize("pos,extension", [(15, 10), (25, 10), (0, 3), (39, 3)])
def test_they_agree_with_record_starts_too(pos, extension):
    """Record confinement is the half most likely to break: it is the only use
    of the site position, and the delegation moves it from an implicit
    parameter to an explicit one."""
    starts = [0, 20, 30]
    window = marked(
        _mark_window,
        40,
        pos=pos,
        extension=extension,
        length=40,
        circular=False,
        record_starts=starts,
    )
    span = marked(
        mark_span,
        40,
        start=pos - extension,
        end=pos + extension,
        length=40,
        circular=False,
        record_starts=starts,
        anchor=pos,
    )

    assert np.array_equal(window, span)


# ---------------------------------------------------------------------------
# What the primitive does on its own, for a span that is not symmetric
# ---------------------------------------------------------------------------


def test_a_one_sided_span_marks_only_its_own_side():
    """The reason this exists. A directional span reaches one way from its
    site, and no symmetric caller can express that."""
    occupied = marked(mark_span, 100, start=50, end=70, length=100, circular=False, anchor=50)

    assert occupied[50:70].all()
    assert not occupied[:50].any()
    assert not occupied[70:].any()


def test_a_one_sided_span_is_confined_by_the_record_holding_its_anchor():
    """A polymerase extending from a site near the end of one contig does not
    continue into the next, and the anchor is what says which contig it is in.
    Passing the span's start instead would put a leftward span in the previous
    record."""
    occupied = marked(
        mark_span,
        40,
        start=5,
        end=25,
        length=40,
        circular=False,
        record_starts=[0, 20, 30],
        anchor=19,
    )

    assert occupied[5:20].all()
    assert not occupied[20:].any(), "the span crossed into the next record"


def test_an_anchorless_span_is_not_record_confined():
    """Absence of an anchor is not a record of zero. A caller with no site to
    name -- a gap interval, say -- gets the unconfined span rather than being
    silently placed in the first record."""
    occupied = marked(
        mark_span, 40, start=15, end=25, length=40, circular=False, record_starts=[0, 20, 30]
    )

    assert occupied[15:25].all()


def test_a_circular_span_wraps_past_the_origin():
    occupied = marked(mark_span, 100, start=-5, end=10, length=100, circular=True, anchor=0)

    assert occupied[0:10].all()
    assert occupied[95:100].all()
    assert not occupied[10:95].any()


def test_a_linear_span_is_clipped_rather_than_wrapped():
    occupied = marked(mark_span, 100, start=-5, end=10, length=100, circular=False, anchor=0)

    assert occupied[0:10].all()
    assert not occupied[10:].any()


def test_an_empty_span_marks_nothing():
    """A zero-width or inverted span is not an error and must not mark. A reach
    of zero is a legitimate request meaning "sites only"."""
    assert not marked(mark_span, 50, start=20, end=20, length=50, circular=False).any()
    assert not marked(mark_span, 50, start=30, end=20, length=50, circular=False).any()
