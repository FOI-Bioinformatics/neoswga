"""The BAM expansion helpers must not raise NameError on their first line.

`_alias_keys` and `_bam_gaps_for_expansion` in `neoswga/cli/iterate.py` both
read `parameter`, and neither imported it. Every other function in that module
imports it locally, so the omission is a slip rather than a convention. Calling
either raised NameError immediately, before doing anything.

Nothing caught it because these are the `--bam` paths, and this repository
contains no BAM or CRAM file: CLAUDE.md records that `calibrate-reach` has
never been run against measured sequencing depth. A path nothing exercises is
a path where an undefined name survives indefinitely.

These tests call the functions rather than inspecting them, because an import
that is present but wrong is still a NameError and only a call finds it.
"""

import pytest

from neoswga.cli import iterate
from neoswga.core.exceptions import ReferenceDataError


def test_alias_keys_runs(monkeypatch):
    """It reads parameter.fg_genomes, so the name must resolve.

    An empty prefix list still executes the parameter read, which is the line
    that raised.
    """
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "fg_genomes", [], raising=False)
    assert iterate._alias_keys([]) == []


def test_alias_keys_handles_a_prefix_without_a_genome(monkeypatch):
    """The truncation is deliberate per its own comment: this runs while
    explaining a failure, so a short genome list must not raise."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "fg_genomes", [], raising=False)
    assert isinstance(iterate._alias_keys(["some_prefix"]), list)


def test_bam_gaps_for_expansion_resolves_its_names(monkeypatch):
    """Reaches the parameter reads before any BAM work, which is where the
    NameError was.

    Only NameError and SystemExit are treated as outcomes of this test. A
    broad `except Exception` here hid a TypeError from calling the function
    with the wrong number of arguments, which made the test pass while
    exercising nothing.
    """
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "fg_genomes", [], raising=False)
    monkeypatch.setattr(parameter, "fg_circular", False, raising=False)

    class Args:
        bam = "/nonexistent.bam"
        min_depth = 1
        min_gap_size = 1
        reference = None
        contig_alias = None

    try:
        iterate._bam_gaps_for_expansion(Args(), [], [], True)
    except NameError as exc:  # the defect this pins
        pytest.fail(f"undefined name reached at runtime: {exc}")
    except SystemExit:
        pass  # pysam absent; the module exits rather than raising
    except ReferenceDataError:
        # Reaching this PROVES the names resolved: the parameter reads are
        # arguments to `bam_gaps`, so they are evaluated before the call that
        # fails on the missing BAM. CLAUDE.md records that this error is
        # deliberately not caught here, so it is the correct outcome.
        pass
    except (OSError, ValueError, RuntimeError):
        pass  # a missing BAM is not what this pins
