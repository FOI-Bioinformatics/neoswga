"""Test edge cases in filter module."""
import pytest


def test_get_bg_rates_via_bloom_empty_list():
    """get_bg_rates_via_bloom should handle an empty primer list without crashing."""
    from neoswga.core.filter import get_bg_rates_via_bloom
    from unittest.mock import MagicMock

    mock_bloom = MagicMock()
    result = get_bg_rates_via_bloom([], mock_bloom)
    assert result == {}
