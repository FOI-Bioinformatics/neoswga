"""Sequencing feedback must know which sequence it is reading.

Task 9 of the 2026-09-21 valid-design plan. Feedback drives redesign: low-depth
regions become targeted additions. So a depth profile read from the wrong
reference does not produce a wrong number a reader can see, it produces a
redesign aimed at gaps in something else.

Two contracts here.

**Coordinate identity.** A BAM contig is bound to a foreground reference by
name, by an explicit alias, or not at all. It is NOT bound because the two
happen to have the same length. Equal length is not identity: two plasmids,
two chromosomes from different assemblies, or a target and its own decoy can
share a length exactly, and when they do, every coordinate lines up and
nothing downstream can notice.

**No leakage.** A calibration fitted on a set of experiments cannot be
evaluated on one of them. `require_disjoint_experiments` is the guard, and it
exists because the failure is silent: an in-sample fit reports a good score
and the report says "held out".
"""

import pytest

from neoswga.core.bam_coverage import match_contigs


# ---------------------------------------------------------------------------
# Coordinate identity
# ---------------------------------------------------------------------------


def test_an_exact_name_match_binds():
    mapping = match_contigs(
        bam_refs=["pcDNA", "pLTR"],
        bam_ref_lengths=[5386, 5386],
        fg_prefixes=["pcDNA"],
        fg_seq_lengths=[5386],
    )

    assert mapping == {"pcDNA": "pcDNA"}


def test_a_basename_match_binds():
    """The prefix is a path; the BAM carries the record name."""
    mapping = match_contigs(
        bam_refs=["pcDNA"],
        bam_ref_lengths=[5386],
        fg_prefixes=["/data/run/pcDNA"],
        fg_seq_lengths=[5386],
    )

    assert mapping == {"/data/run/pcDNA": "pcDNA"}


def test_a_chr_prefix_difference_binds():
    """`chrI` and `I` are two spellings of one name, which is a convention."""
    mapping = match_contigs(
        bam_refs=["chrI"],
        bam_ref_lengths=[230218],
        fg_prefixes=["I"],
        fg_seq_lengths=[230218],
    )

    assert mapping == {"I": "chrI"}


def test_an_explicit_alias_binds():
    mapping = match_contigs(
        bam_refs=["NC_001133.9"],
        bam_ref_lengths=[230218],
        fg_prefixes=["target"],
        fg_seq_lengths=[230218],
        aliases={"target": "NC_001133.9"},
    )

    assert mapping == {"target": "NC_001133.9"}


def test_a_length_match_alone_does_not_bind():
    """The rule this file exists for.

    Names that share nothing, lengths that happen to agree. Binding them makes
    every coordinate line up against a sequence the design was not made for,
    and the depth profile then names gaps in the wrong genome. The user is
    told to pass an alias, which is a claim they can make and this code
    cannot.
    """
    mapping = match_contigs(
        bam_refs=["some_other_contig"],
        bam_ref_lengths=[5386],
        fg_prefixes=["pcDNA"],
        fg_seq_lengths=[5386],
    )

    assert mapping == {}


def test_two_references_of_the_same_length_are_both_unbound():
    """The case that makes length matching indefensible rather than merely weak.

    `pcDNA` and `pLTR` in this repository are both 5,386 bp. Under a
    length rule the only thing deciding which is which is whether the other
    one happened to be in the BAM.
    """
    mapping = match_contigs(
        bam_refs=["contig_a"],
        bam_ref_lengths=[5386],
        fg_prefixes=["pcDNA", "pLTR"],
        fg_seq_lengths=[5386, 5386],
    )

    assert mapping == {}


def test_an_unmatched_prefix_names_the_remedy(caplog):
    import logging

    with caplog.at_level(logging.WARNING):
        match_contigs(
            bam_refs=["unrelated"],
            bam_ref_lengths=[5386],
            fg_prefixes=["pcDNA"],
            fg_seq_lengths=[5386],
        )

    assert "contig-alias" in caplog.text


def test_an_alias_to_a_contig_the_bam_does_not_hold_is_not_honoured():
    """An alias is a claim about this BAM. It cannot invent a record."""
    mapping = match_contigs(
        bam_refs=["present"],
        bam_ref_lengths=[100],
        fg_prefixes=["target"],
        fg_seq_lengths=[100],
        aliases={"target": "absent"},
    )

    assert mapping == {}


def test_a_missing_contig_is_omitted_rather_than_guessed():
    """Partial feedback is legitimate; invented feedback is not."""
    mapping = match_contigs(
        bam_refs=["chrA"],
        bam_ref_lengths=[1000],
        fg_prefixes=["chrA", "chrB"],
        fg_seq_lengths=[1000, 2000],
    )

    assert mapping == {"chrA": "chrA"}


def test_a_name_match_wins_even_when_the_lengths_disagree():
    """A same-named record of a different length is the stale-reference case.

    Binding it is right -- the names agree, which is a claim someone made --
    but the disagreement is worth a warning, because it means the BAM was
    aligned against a different version of this sequence.
    """
    import logging

    mapping = match_contigs(
        bam_refs=["pcDNA"],
        bam_ref_lengths=[9999],
        fg_prefixes=["pcDNA"],
        fg_seq_lengths=[5386],
    )

    assert mapping == {"pcDNA": "pcDNA"}


def test_the_length_disagreement_on_a_name_match_is_reported(caplog):
    import logging

    with caplog.at_level(logging.WARNING):
        match_contigs(
            bam_refs=["pcDNA"],
            bam_ref_lengths=[9999],
            fg_prefixes=["pcDNA"],
            fg_seq_lengths=[5386],
        )

    assert "9999" in caplog.text and "5386" in caplog.text


# ---------------------------------------------------------------------------
# No leakage
# ---------------------------------------------------------------------------


def test_an_experiment_cannot_validate_its_own_fit():
    from neoswga.core.sequencing_feedback import require_disjoint_experiments

    with pytest.raises(ValueError, match="overlap"):
        require_disjoint_experiments(("run_a", "run_b"), ("run_b",))


def test_disjoint_experiments_are_accepted():
    from neoswga.core.sequencing_feedback import require_disjoint_experiments

    require_disjoint_experiments(("run_a", "run_b"), ("run_c",))


def test_the_overlapping_experiment_is_named():
    """A guard that says only "overlap" leaves the user to find which."""
    from neoswga.core.sequencing_feedback import require_disjoint_experiments

    with pytest.raises(ValueError) as caught:
        require_disjoint_experiments(("run_a", "run_b", "run_c"), ("run_b", "run_z"))

    assert "run_b" in str(caught.value)
    assert "run_z" not in str(caught.value), "only the overlapping one is at fault"


def test_an_empty_validation_set_is_refused():
    """Nothing held out is not the same as nothing overlapping.

    A caller passing an empty validation set would pass the disjointness check
    trivially and then report an out-of-sample result computed on no samples.
    """
    from neoswga.core.sequencing_feedback import require_disjoint_experiments

    with pytest.raises(ValueError, match="empty"):
        require_disjoint_experiments(("run_a",), ())


def test_an_empty_training_set_is_refused():
    from neoswga.core.sequencing_feedback import require_disjoint_experiments

    with pytest.raises(ValueError, match="empty"):
        require_disjoint_experiments((), ("run_a",))


def test_a_repeated_experiment_within_one_set_is_not_an_overlap():
    """Listing a run twice is untidy, not leakage."""
    from neoswga.core.sequencing_feedback import require_disjoint_experiments

    require_disjoint_experiments(("run_a", "run_a"), ("run_b",))
