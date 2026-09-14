"""The exclusion genome is the one thing a user explicitly asks to avoid.

`_filter_exclusion_genome` had no test of any kind, and two silent-zero shapes.

**Strand.** The k-mer tables are counted with jellyfish `-C`, so each k-mer and
its reverse complement share one canonical entry and only one of the pair
appears. The lookup compared the primer against the stored key by equality, so a
primer whose canonical form is its reverse complement matched nothing and was
recorded as binding the exclusion genome zero times. Demonstrated on the shipped
plasmid table: `AAACAAATTCCT` is rejected and `AGGAATTTGTTT`, which anneals to
exactly the same sites on the other strand, passes.

**A missing table.** A k-mer file that does not exist contributed no hits, so
every primer passed. That is indistinguishable from an exclusion genome nothing
binds, which is the failure mode Known Issues 5, 6 and 13 all share.
"""

import pytest

from neoswga.core.pipeline import _filter_exclusion_genome
from neoswga.core.thermodynamics import reverse_complement

PRIMER = "AAACAAATTCCT"
RC = reverse_complement(PRIMER)


@pytest.fixture
def excl(tmp_path):
    """A canonical table holding one strand of one k-mer, as jellyfish -C writes."""
    prefix = str(tmp_path / "excl")
    (tmp_path / "excl_12mer_all.txt").write_text(f"{PRIMER} 4\nGGGGGGGGGGGG 1\n")
    return prefix


def test_a_primer_present_in_the_exclusion_genome_is_rejected(excl):
    assert _filter_exclusion_genome([PRIMER], [excl], threshold=0) == [False]


def test_its_reverse_complement_is_rejected_too(excl):
    """Same duplex, other strand. The table stores only the canonical form."""
    assert _filter_exclusion_genome([RC], [excl], threshold=0) == [False]


def test_a_primer_absent_from_the_exclusion_genome_passes(excl):
    assert _filter_exclusion_genome(["ACACACACACAC"], [excl], threshold=0) == [True]


@pytest.mark.parametrize("threshold,expected", [(0, False), (3, False), (4, True), (10, True)])
def test_the_threshold_is_an_inclusive_maximum(excl, threshold, expected):
    """Four hits: rejected below four, allowed at four and above."""
    assert _filter_exclusion_genome([PRIMER], [excl], threshold=threshold) == [expected]


def test_hits_accumulate_across_several_exclusion_genomes(tmp_path):
    first, second = str(tmp_path / "a"), str(tmp_path / "b")
    (tmp_path / "a_12mer_all.txt").write_text(f"{PRIMER} 2\n")
    (tmp_path / "b_12mer_all.txt").write_text(f"{PRIMER} 2\n")

    assert _filter_exclusion_genome([PRIMER], [first, second], threshold=3) == [False]
    assert _filter_exclusion_genome([PRIMER], [first, second], threshold=4) == [True]


def test_a_missing_table_is_not_a_measurement_of_zero(tmp_path):
    """Asking to exclude a genome whose index is absent must not pass silently."""
    with pytest.raises(FileNotFoundError, match="exclusion"):
        _filter_exclusion_genome([PRIMER], [str(tmp_path / "never_counted")], threshold=0)
