"""Measuring the filter must not allocate a second copy of it.

`memory_usage_mb` pickled the whole filter and took `sys.getsizeof` of the
resulting bytes. That allocates a second copy of a structure sized in tens of
MB for a host background, to produce one log line, and it measures the
SERIALISED form rather than the live allocation. `num_bits` is the allocation.
"""

import pickle

import pytest

pytest.importorskip("pybloom_live")

from neoswga.core.background_filter import BackgroundBloomFilter


def test_memory_is_read_off_the_bit_array(monkeypatch):
    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.01)

    def explode(*args, **kwargs):
        raise AssertionError("memory_usage_mb serialised the filter to measure it")

    monkeypatch.setattr(pickle, "dumps", explode)
    reported = bloom.memory_usage_mb()

    assert reported == pytest.approx(bloom.bloom.num_bits / 8 / 1e6, rel=1e-9)


def test_a_bigger_filter_reports_more():
    small = BackgroundBloomFilter(capacity=100_000, error_rate=0.01)
    large = BackgroundBloomFilter(capacity=1_000_000, error_rate=0.01)
    assert large.memory_usage_mb() > small.memory_usage_mb()


def test_the_figure_is_the_allocation_not_the_payload():
    """An empty filter already occupies its full bit array, so the figure must
    not move with how many k-mers have been added."""
    bloom = BackgroundBloomFilter(capacity=100_000, error_rate=0.01)
    before = bloom.memory_usage_mb()
    for i in range(1000):
        bloom.add(f"ACGTACGT{i:04d}"[:12])
    assert bloom.memory_usage_mb() == before
