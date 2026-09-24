"""The depth policy in the output is the one the measurement used.

`coverage_gaps.json` and `primer_expansion`'s record both wrote
`DepthPolicy().to_dict()`, constructed fresh at the write site. That is an
assertion about what ran, not a report of it, and the two agree only while
nothing can configure a policy. `core/depth_policy.py` exists precisely
because "a breadth figure means nothing without the rule that generated it",
so a record that cannot be wrong is also a record that cannot be right.

`bam_gaps` now takes a `policy`, defaults it ONCE, and forwards it to both
branches -- the record path through `bam_depth_profile` and the legacy prefix
path through `compute_bam_depth`. Both callers pass the object they then
record.

This is the Known Issue 8 shape in provenance form: the value was written and
the path that would make it true did not exist. So these assert the PATH, by
driving a real `bam_gaps` with a non-default policy and checking the delivered
gaps move, rather than checking that a `policy=` keyword appears in the source.
"""

import pytest

pysam = pytest.importorskip("pysam")

from neoswga.core.bam_coverage import bam_gaps  # noqa: E402
from neoswga.core.depth_policy import DepthPolicy  # noqa: E402

LENGTH = 1_000


@pytest.fixture
def supplementary_only(tmp_path):
    """A reference covered ONLY by supplementary alignments.

    The default policy excludes them, so the whole contig reads as a gap;
    `count_supplementary=True` counts them and the gap disappears. That makes
    the delivered result depend on the policy, which is the only way to show
    the policy arrived.
    """
    fasta = tmp_path / "g.fasta"
    fasta.write_text(">chr1\n" + "ACGT" * (LENGTH // 4) + "\n")
    path = str(tmp_path / "a.bam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": LENGTH}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for start in range(0, LENGTH - 100, 25):
            read = pysam.AlignedSegment()
            read.query_name = f"s{start}"
            read.query_sequence = "A" * 100
            read.flag = 0x800  # supplementary
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigartuples = [(0, 100)]
            read.query_qualities = pysam.qualitystring_to_array("I" * 100)
            out.write(read)
    pysam.index(path)
    return path, str(fasta)


def _gaps(path, fasta, policy, **kwargs):
    return bam_gaps(
        path,
        ["g"],
        [LENGTH],
        min_depth=1,
        min_gap_size=100,
        fg_genomes=[fasta],
        contig_aliases={"chr1": "chr1"},
        policy=policy,
        **kwargs,
    )


def test_the_policy_reaches_the_record_path(supplementary_only):
    """`fg_genomes` is set, so this is the branch every production caller
    takes: `_bam_gaps_by_record` -> `bam_depth_profile` -> `compute_bam_depth`."""
    path, fasta = supplementary_only

    assert _gaps(path, fasta, DepthPolicy()), "default excludes supplementary, so a gap"
    assert not _gaps(path, fasta, DepthPolicy(count_supplementary=True)), "policy did not arrive"


def test_the_policy_reaches_the_legacy_prefix_path(supplementary_only):
    """No `fg_genomes`, so `match_contigs` binds and `compute_bam_depth` is
    called directly. That call omitted `policy=` entirely."""
    path, _fasta = supplementary_only

    def gaps(policy):
        return bam_gaps(
            path,
            ["chr1"],
            [LENGTH],
            min_depth=1,
            min_gap_size=100,
            policy=policy,
        )

    assert gaps(DepthPolicy())
    assert not gaps(DepthPolicy(count_supplementary=True)), "policy did not arrive"


def test_no_site_records_a_freshly_built_default():
    """The defect itself, as a source check. A record built at the write site
    is right by luck: it reports the default rather than what ran, and nothing
    downstream can tell the two apart."""
    import pathlib

    package = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    offenders = [
        str(path.relative_to(package.parent))
        for path in package.rglob("*.py")
        if "DepthPolicy().to_dict()" in path.read_text()
    ]

    assert not offenders, (
        "record the policy object that was passed to the measurement, not a "
        f"fresh default built beside the record: {offenders}"
    )
