"""Regression tests that ensure multi-genome inputs are fully honored.

The audit surfaced two [0]-indexed sites that silently ignored additional
target or background genomes:

- the position-file cache check, which is now taken per foreground prefix in
  string_search rather than once for all of them in pipeline.py
- background_aware_optimizer.compare_optimizers (now aggregates all bg_prefixes)

These tests lock in the fix.
"""

import pytest


def test_position_file_cache_is_decided_per_fg_prefix(tmp_path):
    """A prefix without a cached HDF5 is scanned even when another prefix has one.

    `pipeline.step2` used to answer this with one `all(...)` over the prefixes at
    `parameter.min_k`, which was read by nothing but a log line. The decision now
    lives in `string_search._split_already_scanned` and is taken per prefix and
    per k, so this asserts it there.
    """
    import h5py

    from neoswga.core import string_search

    k = 8
    primers = ["ACGTACGT", "TTTTGGGG"]
    cached = str(tmp_path / "fg_a")
    uncached = str(tmp_path / "fg_b")

    with h5py.File(f"{cached}_{k}mer_positions.h5", "w") as handle:
        for primer in primers:
            handle.create_dataset(primer, data=[10])

    to_scan, reusable = string_search._split_already_scanned(primers, cached, k)
    assert to_scan == [] and reusable == primers

    to_scan, reusable = string_search._split_already_scanned(primers, uncached, k)
    assert (
        to_scan == primers and reusable == []
    ), "a prefix with no cached position file must still be scanned in full"


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

    assert (
        "for bg_prefix in bg_prefixes" in code
    ), "compare_optimizers must iterate over all bg_prefixes, not bg_prefixes[0]"
    assert "bg_prefixes[0]" not in code, "bg_prefixes[0] hardcoding should be removed from code"
