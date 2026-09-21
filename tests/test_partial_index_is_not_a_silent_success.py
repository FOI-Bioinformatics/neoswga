"""A partly indexed candidate pool is refused, not quietly optimized over.

The step-4 prerequisite validator checks that the position files EXIST. It does
not check that they cover the pool. With three candidates in step3_df.csv and
only the first present in the HDF5 file, the run returned SUCCESS with one
primer and fg_coverage 1.0: the two unindexed primers cover no region, so the
greedy never selects one, and the coverage reported is correct for the panel
that was delivered. Two thirds of the pool were invisible to selection and
nothing said so.
"""

import h5py
import pandas as pd
import pytest

from neoswga.core import pipeline, unified_optimizer
from neoswga.core.pipeline import StepPrerequisiteError

# Chosen to clear the self-dimer screen, which runs before selection: the
# previous three ("AAACCCGGGTTT" and friends) are self-complementary, so once
# the shared search applied that screen the whole frontier was removed and the
# optimizer received an empty list. That made this file's subject -- what the
# UNINDEXED guard does -- unreachable behind an unrelated rejection.
_CANDIDATES = ["GCTAAAGACAAT", "ACGTCAGCACGA", "CAGTGTGAATCG"]


@pytest.fixture
def partial_index(tmp_path):
    """Three candidates in step 3, one of them present in the position file."""
    pd.DataFrame(
        {
            "primer": _CANDIDATES,
            "ratio": [1.0, 2.0, 3.0],
            "gini": [0.3, 0.4, 0.5],
            "fg_count": [20, 15, 10],
            "bg_count": [2, 3, 4],
        }
    ).to_csv(tmp_path / "step3_df.csv", index=False)

    with h5py.File(str(tmp_path / "fg_12mer_positions.h5"), "w") as fh:
        fh.create_dataset(_CANDIDATES[0], data=[100, 5000, 20000])

    return tmp_path


def test_pipeline_run_refuses_a_partly_indexed_pool(partial_index, monkeypatch):
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", str(partial_index), raising=False)
    with pytest.raises(StepPrerequisiteError) as excinfo:
        unified_optimizer.run_optimization(
            method="dominating-set",
            fg_prefixes=[str(partial_index / "fg")],
            fg_seq_lengths=[100000],
            target_size=3,
            verbose=False,
        )
    message = str(excinfo.value)
    assert "STEP 4 PREREQUISITE ERROR" in message
    assert "2 of 3" in message
    assert "neoswga filter" in message


def test_a_caller_supplying_its_own_candidates_is_not_blocked(partial_index, monkeypatch):
    """Same partial index, but the caller passed the pool. expand-primers and
    the analysis commands do this and have not asked step 4 to read step 3."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", str(partial_index), raising=False)
    result = unified_optimizer.run_optimization(
        method="dominating-set",
        candidates=list(_CANDIDATES),
        fg_prefixes=[str(partial_index / "fg")],
        fg_seq_lengths=[100000],
        target_size=3,
        verbose=False,
    )
    assert result is not None


def test_the_count_is_recorded_on_the_result_when_the_run_proceeds(partial_index, monkeypatch):
    """The refusal covers the pipeline path. On the paths that proceed the
    count still has to be auditable, so it reaches the run summary."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", str(partial_index), raising=False)
    result = unified_optimizer.run_optimization(
        method="dominating-set",
        candidates=list(_CANDIDATES),
        fg_prefixes=[str(partial_index / "fg")],
        fg_seq_lengths=[100000],
        target_size=3,
        verbose=False,
    )
    assert result.unindexed_candidates == 2
    assert result.to_dict()["unindexed_candidates"] == 2


def test_a_fully_indexed_pool_is_not_refused(partial_index, monkeypatch):
    """The guard must fire on the missing entries and nothing else."""
    from neoswga.core import parameter

    with h5py.File(str(partial_index / "fg_12mer_positions.h5"), "a") as fh:
        fh.create_dataset(_CANDIDATES[1], data=[300, 9000])
        fh.create_dataset(_CANDIDATES[2], data=[700, 15000])

    monkeypatch.setattr(parameter, "data_dir", str(partial_index), raising=False)
    result = unified_optimizer.run_optimization(
        method="dominating-set",
        fg_prefixes=[str(partial_index / "fg")],
        fg_seq_lengths=[100000],
        target_size=3,
        verbose=False,
    )
    assert result.unindexed_candidates == 0


def test_a_primer_absent_from_the_background_only_is_not_unindexed(partial_index):
    """A specific primer has no host sites and so no background entry. That is
    the design working, not a broken index, and must not trip the guard."""

    class _Cache:
        missing_primers = [("bg", _CANDIDATES[0]), ("bg", _CANDIDATES[1])]

    assert pipeline.unindexed_candidates(_Cache(), ["fg"]) == []


def test_a_primer_missing_from_one_of_two_targets_is_not_unindexed():
    """Multi-target designs legitimately carry primers that bind one target."""

    class _Cache:
        missing_primers = [("fgA", _CANDIDATES[0])]

    assert pipeline.unindexed_candidates(_Cache(), ["fgA", "fgB"]) == []
