"""Tests for neoswga.logging_config (pure stdlib logging helpers)."""

import logging

import pytest

from neoswga import logging_config as lc


def test_get_logger_returns_named_logger():
    log = lc.get_logger("neoswga.test.module")
    assert isinstance(log, logging.Logger)
    assert log.name == "neoswga.test.module"


def test_setup_logging_levels_and_handlers():
    lc.setup_logging(level="DEBUG")
    root = logging.getLogger()
    assert root.level == logging.DEBUG
    assert len(root.handlers) >= 1
    # Re-running replaces handlers (no duplication blow-up)
    lc.setup_logging(level="INFO")
    assert logging.getLogger().level == logging.INFO


def test_setup_logging_format_styles():
    for style in ("default", "detailed", "simple"):
        lc.setup_logging(level="INFO", format_style=style)
        assert logging.getLogger().handlers


def test_setup_logging_to_file(tmp_path):
    log_file = tmp_path / "sub" / "neoswga.log"
    lc.setup_logging(level="INFO", log_file=str(log_file))
    logging.getLogger().info("hello-file")
    for h in logging.getLogger().handlers:
        h.flush()
    assert log_file.exists()


def test_setup_logging_module_levels():
    lc.setup_logging(level="WARNING", module_levels={"neoswga.thermodynamics": "DEBUG"})
    assert logging.getLogger("neoswga.thermodynamics").level == logging.DEBUG


def test_temporary_log_level_restores():
    lc.setup_logging(level="INFO")
    root = logging.getLogger()
    assert root.level == logging.INFO
    with lc.temporary_log_level("DEBUG"):
        assert root.level == logging.DEBUG
    assert root.level == logging.INFO


def test_colored_formatter_outputs_message():
    fmt = lc.ColoredFormatter("%(levelname)s: %(message)s")
    record = logging.LogRecord(
        name="x", level=logging.WARNING, pathname=__file__, lineno=1,
        msg="careful", args=(), exc_info=None,
    )
    out = fmt.format(record)
    assert "careful" in out


def test_progress_logger_update_and_finish(caplog):
    log = lc.get_logger("neoswga.test.progress")
    log.setLevel(logging.INFO)
    with caplog.at_level(logging.INFO, logger="neoswga.test.progress"):
        progress = lc.ProgressLogger(log, total=100, step=25, message="Work")
        for _ in range(100):
            progress.update(1)
        progress.finish()
    # Should have logged at least the 25/50/75/100 milestones + completion
    assert any("Complete" in r.message for r in caplog.records)


def test_configure_from_env(monkeypatch):
    monkeypatch.setenv("NeoSWGA_LOG_LEVEL", "WARNING")
    lc.configure_from_env()
    assert logging.getLogger().level == logging.WARNING
    # restore a sane default for the rest of the suite
    lc.setup_logging(level="INFO")
