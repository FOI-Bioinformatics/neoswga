"""Phase 17C — empty candidate pool guard.

If step 2/3 filtering removes every candidate primer, the optimizer was
previously invoked with an empty list. Different optimizers then
behaved differently (some crashed, some returned empty results
silently). The guard in run_optimization short-circuits with a
directly actionable message instead of dispatching into the factory.
"""

import os
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from neoswga.core.base_optimizer import OptimizationStatus
from neoswga.core.unified_optimizer import run_optimization


def test_empty_candidate_list_returns_failure_with_actionable_message():
    """Passing candidates=[] explicitly must short-circuit with
    status=error and a message that tells the user where to look."""
    result = run_optimization(
        method="hybrid",
        candidates=[],
        fg_prefixes=["fake"],
        fg_seq_lengths=[100000],
        target_size=6,
        verbose=False,
    )
    assert result.status == OptimizationStatus.ERROR
    msg = result.message.lower()
    # The message should name the likely cause(s) so the user does not
    # have to read source code to debug.
    assert "no candidate primers" in msg
    assert "max_bg_freq" in msg or "filter" in msg
    # And the optimizer name must still be set so CLI error reporting works.
    assert result.optimizer_name == "hybrid"


def test_empty_step3_csv_returns_failure(tmp_path, monkeypatch):
    """If candidates is None, run_optimization loads from step3_df.csv.
    An empty CSV (header only) must trigger the guard, not a crash
    downstream in the optimizer."""
    # Write a step3_df.csv with no rows.
    data_dir = tmp_path
    (data_dir / "step3_df.csv").write_text("primer\n")

    # Stub parameter.data_dir to the tmp path.
    from neoswga.core import parameter as param_mod
    monkeypatch.setattr(param_mod, "data_dir", str(data_dir), raising=False)

    result = run_optimization(
        method="hybrid",
        candidates=None,
        fg_prefixes=["fake"],
        fg_seq_lengths=[100000],
        target_size=6,
        verbose=False,
    )
    assert result.status == OptimizationStatus.ERROR
    assert "no candidate primers" in result.message.lower()


def test_failure_result_still_has_valid_metrics():
    """The failure result must carry a usable (empty) metrics object so
    CLI reporting code does not need to special-case None."""
    result = run_optimization(
        method="greedy",
        candidates=[],
        fg_prefixes=["x"],
        fg_seq_lengths=[1000],
        target_size=3,
        verbose=False,
    )
    # PrimerSetMetrics.empty() shape
    assert result.metrics is not None
    assert result.metrics.fg_coverage == 0.0
    assert result.num_primers == 0


def test_guard_preserves_method_name_in_message():
    """The error message echoes the requested method so downstream
    logging (and the test above) can identify what was attempted."""
    for method in ("hybrid", "greedy", "network", "dominating-set"):
        result = run_optimization(
            method=method,
            candidates=[],
            fg_prefixes=["x"],
            fg_seq_lengths=[1000],
            target_size=3,
            verbose=False,
        )
        assert result.status == OptimizationStatus.ERROR
        assert method in result.message
