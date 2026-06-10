"""Test that get_value_or_default logs a warning instead of printing to stdout."""

import logging

from neoswga.core.parameter import get_value_or_default


def test_missing_key_no_stdout(capsys):
    """Missing key must not produce stdout output."""
    result = get_value_or_default(None, {}, "some_missing_key")
    captured = capsys.readouterr()
    assert captured.out == "", f"Unexpected stdout: {captured.out!r}"
    assert result is None


def test_missing_key_logs_warning(caplog):
    """Missing key must emit a warning via the logging module."""
    with caplog.at_level(logging.WARNING, logger="neoswga.core.parameter"):
        result = get_value_or_default(None, {}, "some_missing_key")
    assert result is None
    assert len(caplog.records) == 1
    assert "some_missing_key" in caplog.records[0].message
    assert caplog.records[0].levelno == logging.WARNING
