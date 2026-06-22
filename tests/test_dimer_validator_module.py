"""Coverage for neoswga.core.dimer_validator.DimerValidator.

DimerValidator is on the optimizer hot path (used by base_optimizer) but was at
0% coverage. Most tests drive a controllable mock of the low-level dimer module
so the validator's own logic (caching, matrix, set compatibility, risk fraction,
filtering) is exercised deterministically; one smoke test uses the real backend.
"""

import numpy as np
import pytest

from neoswga.core.dimer_validator import DimerValidator


@pytest.fixture
def mock_dimer(monkeypatch):
    """Control which pairs are 'dimers' via neoswga.core.dimer.is_dimer_fast."""
    import neoswga.core.dimer as dimer

    state = {"dimer_pairs": set(), "calls": 0}

    def fake_is_dimer_fast(a, b, max_bp):
        state["calls"] += 1
        key = (a, b) if a <= b else (b, a)
        return key in state["dimer_pairs"]

    monkeypatch.setattr(dimer, "is_dimer_fast", fake_is_dimer_fast)
    return state


def test_is_compatible_and_caching(mock_dimer):
    mock_dimer["dimer_pairs"].add(("AAAA", "TTTT"))
    v = DimerValidator()
    assert v.is_compatible("GGGG", "CCCC") is True  # not a dimer
    assert v.is_compatible("AAAA", "TTTT") is False  # is a dimer
    # symmetric + cached: repeated/reversed lookups don't re-call the backend
    calls_before = mock_dimer["calls"]
    assert v.is_compatible("TTTT", "AAAA") is False
    assert mock_dimer["calls"] == calls_before  # served from cache


def test_has_self_dimer(mock_dimer):
    mock_dimer["dimer_pairs"].add(("GCGC", "GCGC"))
    v = DimerValidator()
    assert v.has_self_dimer("GCGC") is True
    assert v.has_self_dimer("ATAT") is False


def test_is_set_compatible(mock_dimer):
    v = DimerValidator()
    assert v.is_set_compatible(["AAAA", "GGGG", "CCCC"]) is True
    mock_dimer["dimer_pairs"].add(("CCCC", "GGGG"))
    v.clear_cache()
    assert v.is_set_compatible(["AAAA", "GGGG", "CCCC"]) is False


def test_dimer_risk_fraction(mock_dimer):
    # 3 primers -> 3 pairs; mark exactly one as a dimer -> risk 1/3.
    mock_dimer["dimer_pairs"].add(("AAAA", "CCCC"))
    v = DimerValidator()
    risk = v.dimer_risk(["AAAA", "CCCC", "GGGG"])
    assert risk == pytest.approx(1 / 3)
    assert v.dimer_risk(["AAAA"]) == 0.0  # <2 primers
    assert v.dimer_risk([]) == 0.0


def test_incompatible_pairs(mock_dimer):
    mock_dimer["dimer_pairs"].add(("AAAA", "TTTT"))
    v = DimerValidator()
    pairs = v.incompatible_pairs(["AAAA", "TTTT", "GGGG"])
    assert pairs == [("AAAA", "TTTT")]


def test_filter_self_dimers(mock_dimer):
    mock_dimer["dimer_pairs"].add(("AAAA", "AAAA"))  # AAAA self-dimers
    v = DimerValidator()
    kept = v.filter_self_dimers(["AAAA", "GCTA"])
    assert kept == ["GCTA"]


def test_clear_cache(mock_dimer):
    v = DimerValidator()
    v.is_compatible("AAAA", "GGGG")
    assert v._pair_cache
    v.clear_cache()
    assert v._pair_cache == {}


def test_build_matrix_fast_path(monkeypatch):
    import neoswga.core.dimer as dimer

    primers = ["AAAA", "GGGG", "CCCC"]
    fake = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=int)
    called = {}

    def fake_fast(p, max_dimer_bp):
        called["fast"] = True
        return fake

    monkeypatch.setattr(dimer, "heterodimer_matrix_fast", fake_fast)
    v = DimerValidator()
    m = v.build_matrix(primers)
    assert called.get("fast") is True
    assert m.shape == (3, 3) and m.dtype == bool
    assert m[0, 1] and not m[0, 2]


def test_build_matrix_parallel_path_for_large_pool(monkeypatch):
    import neoswga.core.dimer as dimer

    n = 201  # above the parallel threshold
    primers = [f"P{i:04d}" for i in range(n)]
    called = {}

    def fake_parallel(p, max_dimer_bp):
        called["parallel"] = True
        return np.zeros((len(p), len(p)), dtype=int)

    monkeypatch.setattr(dimer, "heterodimer_matrix_parallel", fake_parallel)
    v = DimerValidator()
    m = v.build_matrix(primers)
    assert called.get("parallel") is True
    assert m.shape == (n, n)


def test_real_backend_smoke():
    """No mock: the real dimer backend must run and return sane types."""
    v = DimerValidator(max_dimer_bp=4)
    risk = v.dimer_risk(["ATCGATCG", "GCTAGCTA", "TTAACCGG"])
    assert 0.0 <= risk <= 1.0
    assert isinstance(v.is_set_compatible(["ATCGATCG", "GCTAGCTA"]), bool)
