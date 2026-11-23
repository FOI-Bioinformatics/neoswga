"""`rescore-set` scores at the configured reach and invents no selectivity.

Two defects (2026-11-23 sweep, tier 2). The reach was the polymerase's
realistic reach whatever `coverage_reach` said, so `rescore-set` and `optimize`
scored one set at different reaches. And when host coverage was zero the
selectivity was reported as `fg_coverage * 1000`, a number nothing measured,
including when no host was configured at all.
"""

import pytest

from neoswga.cli.iterate import _rescore_reach, _selectivity_ratio
from neoswga.core import parameter
from neoswga.core.coverage import polymerase_extension_reach
from neoswga.core.selectivity import MAX_SELECTIVITY


def test_the_configured_reach_is_used(monkeypatch):
    monkeypatch.setattr(parameter, "coverage_reach", 300, raising=False)
    assert _rescore_reach("phi29") == 300


def test_without_a_configured_reach_the_polymerase_decides(monkeypatch):
    monkeypatch.setattr(parameter, "coverage_reach", None, raising=False)
    assert _rescore_reach("bst") == polymerase_extension_reach("bst", coverage_metric="realistic")


def test_no_background_means_no_selectivity():
    assert _selectivity_ratio(0.8, 0.0, has_background=False) is None


def test_an_unbound_background_is_the_sentinel_not_a_multiple():
    assert _selectivity_ratio(0.8, 0.0, has_background=True) == MAX_SELECTIVITY


def test_a_measured_background_gives_the_ratio():
    assert _selectivity_ratio(0.8, 0.2, has_background=True) == pytest.approx(4.0)
