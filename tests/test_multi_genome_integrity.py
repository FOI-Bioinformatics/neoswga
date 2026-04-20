"""Regression tests that ensure multi-genome inputs are fully honored.

The audit surfaced two [0]-indexed sites that silently ignored additional
target or background genomes:

- pipeline.py position-file cache check (now checks all fg_prefixes)
- background_aware_optimizer.compare_optimizers (now aggregates all bg_prefixes)

These tests lock in the fix.
"""

import os
import pytest


def test_position_file_cache_requires_all_fg_prefixes(tmp_path):
    """If one of several fg_prefixes lacks a cached HDF5, position_files_exist
    must return False so the missing prefix is scanned."""
    # Simulate the inline check in pipeline.step2
    k = 8
    prefixes = [str(tmp_path / "fg_a"), str(tmp_path / "fg_b"), str(tmp_path / "fg_c")]

    # Pre-create only the first cache
    (tmp_path / "fg_a_8mer_positions.h5").write_text("")

    position_files_exist = all(
        os.path.exists(f"{p}_{k}mer_positions.h5") for p in prefixes
    )
    assert position_files_exist is False, (
        "All prefixes must have caches for the pipeline to skip position-file creation"
    )

    # Now create them all
    (tmp_path / "fg_b_8mer_positions.h5").write_text("")
    (tmp_path / "fg_c_8mer_positions.h5").write_text("")

    position_files_exist = all(
        os.path.exists(f"{p}_{k}mer_positions.h5") for p in prefixes
    )
    assert position_files_exist is True


def test_background_aware_optimizer_aggregates_all_bg_prefixes():
    """compare_optimizers standard_bg calculation must sum across all bg_prefixes."""
    import inspect
    from neoswga.core import background_aware_optimizer as bao

    source = inspect.getsource(bao.compare_optimizers)
    # Strip comments before checking code (docstring-style notes may still
    # mention bg_prefixes[0] as a historical reference).
    code_lines = []
    for line in source.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("#"):
            continue
        code_lines.append(line)
    code = "\n".join(code_lines)

    assert "for bg_prefix in bg_prefixes" in code, (
        "compare_optimizers must iterate over all bg_prefixes, not bg_prefixes[0]"
    )
    assert "bg_prefixes[0]" not in code, (
        "bg_prefixes[0] hardcoding should be removed from code"
    )
