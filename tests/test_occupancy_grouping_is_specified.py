"""The occupancy approximation groups by primer, not by site.

Audit finding F5. The docstring of `_compute_effective_coverage` stated

    P(covered at x) = 1 - PRODUCT over sites s reaching x of (1 - theta)

while the implementation forms the union of one primer's windows, applies that
primer's occupancy once, and only then combines distinct primers. The two are
different models and give different numbers, and nothing said which was meant.

Neither has been validated against a measured reaction. The grouping is the one
that ships and it is deliberate: a primer's own overlapping windows do not stack
with each other, because the primer does not compete with itself for its own
site. Per-site independence would multiply (1 - theta) in once per overlapping
window and report a larger number, and the audit is explicit that swapping to it
and presenting the larger figure as validated is not an improvement.

So this test does not argue for either model. It pins the one in use, with the
audit's own arithmetic, so a silent change to the other fails here.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import BaseOptimizer

LENGTH = 1_000
REACH = 100
THETA = 0.5  # what site_occupancy returns at T = Tm


@pytest.fixture
def fixed_occupancy(monkeypatch):
    """Impose theta = 0.5 so the arithmetic below is exact."""
    monkeypatch.setattr("neoswga.core.occupancy.site_occupancy", lambda *a, **k: THETA)


def _coverage(positions_by_primer):
    holder = SimpleNamespace(
        conditions=SimpleNamespace(temp=30.0, calculate_effective_tm=lambda seq, **k: 30.0),
        config=SimpleNamespace(extension_reach=REACH, fg_circular=False),
    )
    return BaseOptimizer._compute_effective_coverage(holder, positions_by_primer, LENGTH)


def test_one_primer_with_two_overlapping_sites_is_weighted_once(fixed_occupancy):
    """Union of [400,610) is 210 bases, carried at theta once: 105 of 1000."""
    coverage = _coverage({"AAAAAAAAAAAA": [500, 510]})

    assert coverage == pytest.approx(0.105, abs=5e-4)


def test_two_distinct_primers_at_the_same_sites_combine(fixed_occupancy):
    """Their windows DO stack where they overlap: 152.5 of 1000.

    190 overlapping bases at 1 - 0.25, plus 20 bases reached by one primer
    only at 0.5.
    """
    coverage = _coverage({"AAAAAAAAAAAA": [500], "CCCCCCCCCCCC": [510]})

    assert coverage == pytest.approx(0.1525, abs=5e-4)


def test_the_grouping_is_what_makes_the_two_differ(fixed_occupancy):
    """Stated as a relation, so it survives a change of constants."""
    one = _coverage({"AAAAAAAAAAAA": [500, 510]})
    two = _coverage({"AAAAAAAAAAAA": [500], "CCCCCCCCCCCC": [510]})

    assert two > one, "distinct primers must stack where one primer's own windows do not"
