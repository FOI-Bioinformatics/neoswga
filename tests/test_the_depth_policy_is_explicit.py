"""Which reads count toward depth should be a stated choice, not a default.

Finding F8 of the 2026-09-16 pipeline audit, the read-selection half.
`compute_bam_depth` called `count_coverage` with `quality_threshold=0` and
pysam's `read_callback='all'`, and nothing recorded or printed what that meant.
Two of its consequences are wrong for this application and neither was visible.

**Duplicates were excluded.** pysam's `'all'` skips `BAM_FDUP`, and every
duplicate marker assumes that two fragments sharing a start coordinate are one
molecule sequenced twice. Hyperbranched amplification violates that assumption
by construction: multiple displacement amplification generates independent
priming events at the same position, so a duplicate marker removes real signal
and makes an amplified region look unamplified. That is the one reversal this
policy makes deliberately.

**Supplementary alignments were counted.** `'all'` does not skip
`BAM_FSUPPLEMENTARY`, so a chimeric read contributed at its primary locus and
again at every split. Branched products make these common, and counting them
counts one molecule at several places.

Both defaults now say so, and both are recorded rather than assumed.
"""

import pytest

from neoswga.core.depth_policy import DepthPolicy


class _Read:
    def __init__(
        self,
        mapping_quality=60,
        is_duplicate=False,
        is_supplementary=False,
        is_secondary=False,
        is_qcfail=False,
        is_unmapped=False,
    ):
        self.mapping_quality = mapping_quality
        self.is_duplicate = is_duplicate
        self.is_supplementary = is_supplementary
        self.is_secondary = is_secondary
        self.is_qcfail = is_qcfail
        self.is_unmapped = is_unmapped


class TestTheDefaultsAreChosenForThisApplication:
    def test_duplicates_are_counted(self):
        """The reversal. Every other tool's default is the opposite, and for
        MDA the opposite discards real coverage."""
        assert DepthPolicy().count_duplicates is True

    def test_supplementary_alignments_are_not(self):
        """They were, silently, because pysam's `'all'` does not skip them."""
        assert DepthPolicy().count_supplementary is False

    def test_secondary_and_qc_failed_reads_are_not(self):
        policy = DepthPolicy()

        assert policy.count_secondary is False
        assert policy.count_qcfail is False

    def test_no_mapping_quality_floor_by_default(self):
        """A floor would turn an unmappable repeat into a coverage gap, and a
        gap is what expansion then designs primers for. An ambiguously mapped
        read still came from somewhere; the region amplified."""
        assert DepthPolicy().min_mapping_quality == 0

    def test_no_base_quality_floor_by_default(self):
        """The question is whether a region amplified, not whether a base call
        is trustworthy. pysam's own default of 15 answers the other one."""
        assert DepthPolicy().min_base_quality == 0


class TestTheFilterIsWhatThePolicySays:
    def test_a_plain_read_is_counted(self):
        assert DepthPolicy().accepts(_Read()) is True

    def test_a_duplicate_is_counted_by_default(self):
        assert DepthPolicy().accepts(_Read(is_duplicate=True)) is True

    def test_a_duplicate_is_dropped_when_asked(self):
        assert DepthPolicy(count_duplicates=False).accepts(_Read(is_duplicate=True)) is False

    def test_a_supplementary_alignment_is_dropped(self):
        assert DepthPolicy().accepts(_Read(is_supplementary=True)) is False

    def test_a_supplementary_alignment_is_kept_when_asked(self):
        policy = DepthPolicy(count_supplementary=True)

        assert policy.accepts(_Read(is_supplementary=True)) is True

    @pytest.mark.parametrize("flag", ["is_secondary", "is_qcfail", "is_unmapped"])
    def test_the_standard_exclusions_hold(self, flag):
        assert DepthPolicy().accepts(_Read(**{flag: True})) is False

    def test_an_unmapped_read_is_never_counted(self):
        """Not configurable: an unmapped read has no position to count at."""
        assert "unmapped" not in DepthPolicy().to_dict()

    def test_a_mapping_quality_floor_is_applied_when_set(self):
        policy = DepthPolicy(min_mapping_quality=30)

        assert policy.accepts(_Read(mapping_quality=29)) is False
        assert policy.accepts(_Read(mapping_quality=30)) is True


class TestItIsRecordable:
    def test_every_field_reaches_the_dictionary(self):
        policy = DepthPolicy(min_mapping_quality=7, count_duplicates=False)
        payload = policy.to_dict()

        for field in policy.__dataclass_fields__:
            assert field in payload, f"{field} would not be recorded in an output file"

    def test_it_round_trips_as_json(self):
        import json

        payload = DepthPolicy().to_dict()

        assert json.loads(json.dumps(payload)) == payload

    def test_it_describes_itself_in_one_line(self):
        """Printed beside every breadth figure, so a reader knows what was
        counted without opening a file."""
        line = DepthPolicy().describe()

        assert "\n" not in line
        assert "duplicates" in line.lower()

    def test_the_description_names_a_non_default_choice(self):
        line = DepthPolicy(min_mapping_quality=30).describe()

        assert "30" in line


class TestWhatThisPolicyDeliberatelyDoesNotCover:
    """Two of the plan's knobs are not offered, each for its own reason.

    A field that records a choice nothing enforces is the defect this project
    keeps closing, so neither is present as an inert setting.
    """

    def test_overlapping_mates_are_not_a_setting(self):
        """`count_coverage` cannot deduplicate an overlapping pair; that needs
        the pileup API. Offering the knob here would record a policy that was
        never applied."""
        assert not [f for f in DepthPolicy.__dataclass_fields__ if "overlap" in f]

    def test_deletions_are_not_a_setting(self):
        """`count_coverage` tallies A/C/G/T from the read sequence, so a
        deleted base contributes nothing whatever a policy said. It is a
        property of the counter, not a choice."""
        assert not [f for f in DepthPolicy.__dataclass_fields__ if "deletion" in f]

    def test_the_docstring_says_so(self):
        """Otherwise the absence reads as an oversight."""
        text = DepthPolicy.__doc__ or ""

        assert "overlap" in text.lower()
        assert "deletion" in text.lower()


class TestTheDepthPathUsesIt:
    def test_compute_bam_depth_takes_a_policy(self):
        import inspect

        from neoswga.core.bam_coverage import compute_bam_depth

        assert "policy" in inspect.signature(compute_bam_depth).parameters

    def test_the_profile_carries_the_policy_it_used(self):
        import inspect

        from neoswga.core.bam_coverage import DepthProfile

        assert "policy" in DepthProfile.__dataclass_fields__
        assert "policy" in inspect.signature(
            __import__("neoswga.core.bam_coverage", fromlist=["x"]).bam_depth_profile
        ).parameters
