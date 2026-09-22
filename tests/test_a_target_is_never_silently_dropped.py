"""Pairing prefixes with lengths must not silently lose a target.

`bam_coverage` pairs `fg_prefixes` with `fg_seq_lengths` using `zip`, which
stops at the shorter sequence. The two lists come from separate params.json
keys and nothing on this path checks they agree, so a file naming three targets
and two lengths analyses two of them:

    match_contigs(refs, lens, ["chrA","chrB","chrC"], [1000, 2000, 3000])
        -> chrA, chrB, chrC
    match_contigs(refs, lens, ["chrA","chrB","chrC"], [1000, 2000])
        -> chrA, chrB

No error, no warning. The third target is not analysed, and every gap in it is
then reported as absent -- so `expand-primers`, whose job is to add oligos for
the gaps, would never add one for that target.

`optimizer_factory.create` already refuses a mismatch, which is why the
optimizer path is safe and this one is not: the coverage commands
(`calibrate-reach`, `analyze-coverage`) do not build an optimizer.

Found 2026-09-22 by reading the `zip()`-without-`strict` lint findings rather
than by a failing test, which is the argument for making that check blocking.
Of 64 in the package this is the group that matters: the rest pair values with
their own derived sequences, where a mismatch is impossible by construction.

This is the silent-zero family in a new shape. Not a scan that found nothing,
an integer that saturated, a cache asked for what it does not hold, or a
dictionary lookup whose fallback is the answer everyone wants -- a pairing that
truncated.
"""

import pytest

from neoswga.core.bam_coverage import match_contigs

BAM_REFS = ["chrA", "chrB", "chrC"]
BAM_LENGTHS = [1000, 2000, 3000]


def test_every_named_target_is_matched_when_the_lists_agree():
    """Guard the guard: a function that refuses everything passes the rest."""
    matched = match_contigs(BAM_REFS, BAM_LENGTHS, ["chrA", "chrB", "chrC"], [1000, 2000, 3000])

    assert sorted(matched) == ["chrA", "chrB", "chrC"]


def test_too_few_lengths_is_refused_rather_than_truncated():
    """The defect. Three targets and two lengths analysed two of them."""
    with pytest.raises(ValueError) as excinfo:
        match_contigs(BAM_REFS, BAM_LENGTHS, ["chrA", "chrB", "chrC"], [1000, 2000])

    message = str(excinfo.value)
    assert "fg_prefixes" in message and "fg_seq_lengths" in message
    assert "3" in message and "2" in message, "the message must name both counts"


def test_too_many_lengths_is_refused_too():
    """The other direction is just as wrong and was equally silent."""
    with pytest.raises(ValueError):
        match_contigs(BAM_REFS, BAM_LENGTHS, ["chrA", "chrB"], [1000, 2000, 3000])


def test_the_gap_finder_refuses_the_same_mismatch(tmp_path):
    """`bam_gaps` pairs the same two lists and had the same hole.

    It is the function `expand-primers` and `analyze-coverage` reach, so a
    dropped target there means no oligo is ever proposed for it.
    """
    from neoswga.core.bam_coverage import bam_gaps

    with pytest.raises(ValueError):
        bam_gaps(str(tmp_path / "absent.bam"), ["chrA", "chrB"], [1000])


def test_the_refusal_names_the_two_keys_a_user_would_edit():
    """A count mismatch is a params.json problem, so the message has to point
    at params.json rather than at a function signature."""
    with pytest.raises(ValueError) as excinfo:
        match_contigs(BAM_REFS, BAM_LENGTHS, ["chrA", "chrB", "chrC"], [1000])

    assert "params.json" in str(excinfo.value)
