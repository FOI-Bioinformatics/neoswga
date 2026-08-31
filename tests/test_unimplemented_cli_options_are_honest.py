"""An option that is not implemented must say so.

The sibling file `test_design_options_have_effect.py` pins the options that are
meant to steer a design: vary them and the design has to change. This file
covers the other half of the same failure, the options that steer nothing.

Several flags were accepted, written to the global parameter object, and read
by no module -- while logging a message that read as confirmation. `--use-gpu`
logged "GPU acceleration enabled", `--use-enhanced-features` logged "Using
enhanced 120+ feature model", and neither attribute had a consumer anywhere
under `neoswga/core/`. A user reading the log had no way to tell an option that
worked from one that was decoration, which is worse than an option that is
absent: an absent option fails loudly at the command line.

The flags are kept, because removing a published argument breaks callers'
scripts. What is removed is the false confirmation. Each test below asserts
both halves: the flag is still accepted, and the run says plainly that it does
nothing.
"""

import argparse
import logging
from types import SimpleNamespace

import pytest

from neoswga.cli._common import (
    UNIMPLEMENTED_OPTIONS,
    report_unimplemented_options,
    setup_gpu_acceleration,
)
from neoswga.cli_unified import create_parser


def _parameter_stub():
    return SimpleNamespace()


def all_help_text():
    """Help for every subcommand, not just the top-level one.

    `create_parser().format_help()` lists the subcommand names and none of
    their arguments, so asserting a flag's absence against it passes whether
    or not the flag exists.
    """
    parser = create_parser()
    chunks = [parser.format_help()]
    for action in parser._subparsers._group_actions:
        for sub in action.choices.values():
            chunks.append(sub.format_help())
    return "\n".join(chunks)


# ----------------------------------------------------------------------
# The flags still exist
# ----------------------------------------------------------------------

PUBLISHED_FLAGS = [
    "--use-gpu",
    "--no-gpu",
    "--gpu-device",
    "--use-enhanced-features",
    "--enhanced-model-path",
    "--use-background-filter",
    "--background-bloom-path",
    "--background-sampled-path",
    "--max-optimization-time",
    "--fast-score",
]


@pytest.mark.parametrize("flag", PUBLISHED_FLAGS)
def test_the_flag_is_still_accepted(flag):
    """Making an option honest must not remove it.

    Deleting a published CLI argument turns a working script into a parse
    error, which is a decision for the repository owner and not a side effect
    of correcting a log message.
    """
    assert flag in all_help_text()


# ----------------------------------------------------------------------
# GPU
# ----------------------------------------------------------------------


def test_use_gpu_does_not_claim_acceleration(caplog):
    """--use-gpu logged "GPU acceleration enabled (device 0)" and set
    `parameter.use_gpu`, which no module reads. Nothing was accelerated.
    """
    args = argparse.Namespace(use_gpu=True, no_gpu=False, gpu_device=0)
    parameter = _parameter_stub()

    with caplog.at_level(logging.DEBUG, logger="neoswga.cli._common"):
        setup_gpu_acceleration(args, parameter)

    text = caplog.text
    assert "GPU acceleration enabled" not in text
    assert "not implemented" in text.lower()
    assert "--use-gpu" in text


def test_gpu_autodetection_does_not_claim_acceleration(caplog, monkeypatch):
    """The auto-detect branch made the same claim without the user asking for
    it: install CuPy for any reason and every run announced acceleration it
    was not doing.
    """
    import neoswga.core.gpu_acceleration as gpu

    monkeypatch.setattr(gpu, "GPU_AVAILABLE", True)

    args = argparse.Namespace(use_gpu=False, no_gpu=False, gpu_device=0)

    with caplog.at_level(logging.DEBUG, logger="neoswga.cli._common"):
        setup_gpu_acceleration(args, _parameter_stub())

    assert "auto-detected and enabled" not in caplog.text


def test_no_gpu_is_silent(caplog):
    """--no-gpu asks for no GPU work and gets no GPU work, so there is nothing
    to report. Only a request the tool cannot honour needs a warning.
    """
    args = argparse.Namespace(use_gpu=False, no_gpu=True, gpu_device=0)

    with caplog.at_level(logging.WARNING, logger="neoswga.cli._common"):
        setup_gpu_acceleration(args, _parameter_stub())

    assert caplog.text.strip() == ""


def test_use_gpu_help_does_not_promise_a_speedup():
    """The help text advertised "10-100x speedup for large genomes".

    `gpu_acceleration.batch_calculate_tm` loops in Python and writes
    element-by-element into a CuPy array, which its own docstring describes as
    providing minimal speedup over NumPy. The figure was never measured on the
    path the flag claims to accelerate.
    """
    assert "10-100x" not in all_help_text()


# ----------------------------------------------------------------------
# Enhanced random-forest features
# ----------------------------------------------------------------------


def test_enhanced_features_reports_that_scoring_is_unchanged(caplog):
    """`parameter.use_enhanced_features` had no reader.

    `rf_preprocessing.predict_with_enhanced_features` exists but is called by
    nothing in the scoring path, and the model it loads is a pickle that the
    package does not ship.
    """
    args = argparse.Namespace(
        use_enhanced_features=True,
        enhanced_model_path=None,
    )

    with caplog.at_level(logging.WARNING, logger="neoswga.cli._common"):
        reported = report_unimplemented_options(args)

    assert "--use-enhanced-features" in reported
    assert "Using enhanced 120+ feature model" not in caplog.text
    assert "not implemented" in caplog.text.lower()


def test_enhanced_model_path_reports_that_it_is_not_loaded(caplog):
    args = argparse.Namespace(
        use_enhanced_features=False,
        enhanced_model_path="/tmp/enhanced_rf_model.pkl",
    )

    with caplog.at_level(logging.WARNING, logger="neoswga.cli._common"):
        reported = report_unimplemented_options(args)

    assert "--enhanced-model-path" in reported


# ----------------------------------------------------------------------
# Background filtering in optimize
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "attr,value,flag",
    [
        ("use_background_filter", True, "--use-background-filter"),
        ("background_bloom_path", "/tmp/bg.bloom", "--background-bloom-path"),
        ("background_sampled_path", "/tmp/bg.sampled", "--background-sampled-path"),
    ],
)
def test_optimize_background_filter_options_report_that_they_are_unused(caplog, attr, value, flag):
    """`optimize` builds no Bloom filter and reads no pre-built index.

    `use_background_filter` reached `optimize_step4` as a named argument that
    the function never passed on. The Bloom paths are consumed only by
    `improved_pipeline.ImprovedSWGAPipeline`, which the `optimize` command does
    not use. Candidate pre-filtering against the background is a separate
    mechanism, controlled by --no-bg-prefilter.
    """
    args = argparse.Namespace(**{attr: value})

    with caplog.at_level(logging.WARNING, logger="neoswga.cli._common"):
        reported = report_unimplemented_options(args)

    assert flag in reported
    assert "--no-bg-prefilter" in caplog.text


def test_unset_background_options_are_silent(caplog):
    """Only an option the caller actually set is worth a warning."""
    args = argparse.Namespace(
        use_background_filter=False,
        background_bloom_path=None,
        background_sampled_path=None,
        max_optimization_time=None,
    )

    with caplog.at_level(logging.WARNING, logger="neoswga.cli._common"):
        reported = report_unimplemented_options(args)

    assert reported == []
    assert caplog.text.strip() == ""


# ----------------------------------------------------------------------
# Wall-clock budget
# ----------------------------------------------------------------------


def test_max_optimization_time_reports_that_no_budget_is_enforced(caplog):
    """No optimizer has a wall-clock budget.

    `docs/development/optimizers.md` documented an `OptimizerConfig.
    timeout_seconds` field that does not exist on the dataclass, so the flag
    had nothing to configure and a long run was never interrupted.
    """
    args = argparse.Namespace(max_optimization_time=60)

    with caplog.at_level(logging.WARNING, logger="neoswga.cli._common"):
        reported = report_unimplemented_options(args)

    assert "--max-optimization-time" in reported
    assert "not implemented" in caplog.text.lower()


def test_max_optimization_time_does_not_advertise_a_default():
    """A default of 300 seconds describes an enforcement that does not
    happen. The argument defaults to None so that "unset" and "asked for a
    budget" stay distinguishable.
    """
    args = create_parser().parse_args(["optimize"])

    assert args.max_optimization_time is None


def test_optimizer_config_has_no_timeout_field():
    """Pins the reason --max-optimization-time cannot be honoured, so that
    wiring it becomes a visible change rather than a silent one.
    """
    from neoswga.core.base_optimizer import OptimizerConfig

    assert not hasattr(OptimizerConfig(), "timeout_seconds")


# ----------------------------------------------------------------------
# Every table entry names a real argparse destination
# ----------------------------------------------------------------------


def test_every_unimplemented_option_names_a_real_flag():
    """A stale entry here would warn about a flag the parser does not offer,
    or silently stop warning about one it renamed.
    """
    help_text = all_help_text()

    for flag, _detail in UNIMPLEMENTED_OPTIONS.values():
        assert flag in help_text, flag
