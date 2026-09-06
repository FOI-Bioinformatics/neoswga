"""Step 4 refuses to optimize against an index it cannot trust.

Finding E3. validate_step4_prerequisites has existed since the step validators
were written and has never been called. Without it, a missing or truncated
_positions.h5 leaves the position cache warning that the coverage number is
"meaningless, not low", and the run then selects a set and reports it.
"""

import os

import pandas as pd
import pytest

from neoswga.core import unified_optimizer
from neoswga.core.pipeline import StepPrerequisiteError, validate_step4_prerequisites


@pytest.fixture
def data_dir(tmp_path):
    pd.DataFrame({"primer": ["AAACCCGGGTTT", "ACCCGGGTTTAA"],
                  "ratio": [1.0, 2.0], "gini": [0.3, 0.4],
                  "fg_count": [20, 15], "bg_count": [2, 3]}).to_csv(
        tmp_path / "step3_df.csv", index=False
    )
    return str(tmp_path)


def test_validator_rejects_missing_position_files(data_dir, tmp_path):
    result = validate_step4_prerequisites(data_dir, [str(tmp_path / "fg")])
    assert result.valid is False
    assert any("12mer_positions.h5" in f for f in result.missing_files)


def test_validator_accepts_present_position_files(data_dir, tmp_path):
    import h5py

    with h5py.File(str(tmp_path / "fg_12mer_positions.h5"), "w") as fh:
        fh.create_dataset("AAACCCGGGTTT", data=[1, 2, 3])
    result = validate_step4_prerequisites(data_dir, [str(tmp_path / "fg")])
    assert result.valid is True


def test_validator_remediation_does_not_name_the_retired_gate(tmp_path):
    """min_amp_pred was retired on 2026-09-05; the advice must not name it."""
    pd.DataFrame({"primer": []}).to_csv(tmp_path / "step3_df.csv", index=False)
    result = validate_step4_prerequisites(str(tmp_path), [str(tmp_path / "fg")])
    assert result.valid is False
    assert "min-amp-pred" not in result.remediation
    assert "min_amp_pred" not in result.remediation


def test_run_optimization_raises_rather_than_scoring_against_nothing(
    data_dir, tmp_path, monkeypatch
):
    """The guard sits in run_optimization's candidate-loading block, which is
    where step 3 is read. optimize_step4 delegates to it."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", data_dir, raising=False)
    with pytest.raises(StepPrerequisiteError) as excinfo:
        unified_optimizer.run_optimization(
            method="dominating-set",
            fg_prefixes=[str(tmp_path / "fg")],
            fg_seq_lengths=[100000],
            target_size=2,
            verbose=False,
        )
    assert "STEP 4 PREREQUISITE ERROR" in str(excinfo.value)


def test_a_caller_supplying_its_own_candidates_is_not_blocked(data_dir, tmp_path, monkeypatch):
    """expand-primers and the tests pass their own list and have not asked step 4
    to read step 3. The guard must not fire for them."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", data_dir, raising=False)
    result = unified_optimizer.run_optimization(
        method="dominating-set",
        candidates=["AAACCCGGGTTT", "ACCCGGGTTTAA"],
        fg_prefixes=[str(tmp_path / "fg")],
        fg_seq_lengths=[100000],
        target_size=2,
        verbose=False,
    )
    assert result is not None  # it may fail on empty coverage; it must not raise
