"""Which reads count toward depth, stated rather than inherited.

Finding F8 of the 2026-09-16 pipeline audit, the read-selection half.
`compute_bam_depth` called pysam's `count_coverage` with `quality_threshold=0`
and `read_callback='all'`, and nothing recorded or printed what that meant. Two
of its consequences are wrong for this application and neither was visible in
any output.

**Duplicates were excluded, and should not be.** pysam's `'all'` skips
`BAM_FDUP`. Every duplicate marker rests on the assumption that two fragments
sharing a start coordinate are one molecule sequenced twice. Multiple
displacement amplification violates that by construction: hyperbranched
priming generates independent events at the same position, so marking them
duplicate removes real signal and makes an amplified region look unamplified.
Since a region that looks unamplified is what BAM-guided expansion designs
primers for, this one mattered twice over. Counting them is the single
reversal this policy makes deliberately.

**Supplementary alignments were counted, and should not be.** `'all'` does not
skip `BAM_FSUPPLEMENTARY`, so a chimeric read contributed at its primary locus
and again at every split. Branched products make chimeras common, and counting
them counts one molecule in several places.

Two defaults are deliberately permissive, for a reason specific to this
question rather than to sequencing in general.

**No mapping-quality floor.** A floor turns an unmappable repeat into a
coverage gap, and a gap is what expansion then designs primers for. An
ambiguously mapped read still came from somewhere, and the region it came from
did amplify. Raise it when the question is variant calling rather than breadth.

**No base-quality floor.** The question is whether a region amplified, not
whether a base call is trustworthy. pysam's own default of 15 answers the
other one.

Two knobs the audit plan lists are NOT offered here, because nothing on this
path would enforce them and a recorded policy that was never applied is the
defect this project keeps closing:

- **Overlapping mates.** Deduplicating an overlapping pair needs the pileup
  API; `count_coverage` cannot do it. A fragment whose mates overlap is
  therefore counted twice in that overlap.
- **Deletions.** `count_coverage` tallies A/C/G/T from the read sequence, so a
  base deleted relative to the reference contributes nothing whatever a policy
  said. That is a property of the counter, not a choice.
"""

from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any, Dict


@dataclass(frozen=True)
class DepthPolicy:
    """The read-selection rules one depth measurement was taken under.

    Frozen, and carried beside the depth it produced: a breadth figure means
    nothing without the rule that generated it, and two runs under different
    rules are not comparable. See this module's docstring for why each default
    is what it is, and for the two knobs (overlapping mates, deletions) that
    are deliberately absent because nothing here would enforce them.
    """

    min_mapping_quality: int = 0
    min_base_quality: int = 0
    count_duplicates: bool = True
    count_supplementary: bool = False
    count_secondary: bool = False
    count_qcfail: bool = False

    def accepts(self, read: Any) -> bool:
        """Whether this read contributes to depth.

        Unmapped reads are never counted and that is not configurable: an
        unmapped read has no position to be counted at.
        """
        if getattr(read, "is_unmapped", False):
            return False
        if not self.count_secondary and getattr(read, "is_secondary", False):
            return False
        if not self.count_supplementary and getattr(read, "is_supplementary", False):
            return False
        if not self.count_qcfail and getattr(read, "is_qcfail", False):
            return False
        if not self.count_duplicates and getattr(read, "is_duplicate", False):
            return False
        if self.min_mapping_quality:
            if int(getattr(read, "mapping_quality", 0) or 0) < self.min_mapping_quality:
                return False
        return True

    def to_dict(self) -> Dict[str, Any]:
        """Plain data for `coverage_gaps.json`, the reach output and the manifest."""
        return {field.name: getattr(self, field.name) for field in fields(self)}

    def describe(self) -> str:
        """One line, printed beside every breadth figure."""
        counted = ["duplicates"] if self.count_duplicates else []
        counted += ["supplementary"] if self.count_supplementary else []
        counted += ["secondary"] if self.count_secondary else []
        counted += ["QC-failed"] if self.count_qcfail else []
        parts = [
            f"counting {', '.join(counted)}" if counted else "counting primary alignments only",
            f"MAPQ >= {self.min_mapping_quality}" if self.min_mapping_quality else "no MAPQ floor",
            (
                f"base quality >= {self.min_base_quality}"
                if self.min_base_quality
                else "no base-quality floor"
            ),
        ]
        return "Depth policy: " + "; ".join(parts) + "."
