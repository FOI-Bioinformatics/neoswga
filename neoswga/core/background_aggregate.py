"""What "the background" means when a design names several host genomes.

Today every host is pooled. `core/filter.py` divides a summed background count
by the summed background length, and `core/base_optimizer.py` sums site loads
across background prefixes. That is one reading of "the background" -- every
listed host present at once -- and it is the only one the pipeline can express.

It is not the only reading, and the difference was measured. Against the
*Wolbachia* reference with three hosts pooled (3.3 Mb, 144 Mb and 3.30 Gb, a
span of 985), hg38 is 95.7 percent of the pooled length, so the pooled gate's
ceiling is within 5 percent of hg38's own while the numerator still charges a
candidate for its *Drosophila* sites. Under that gate a candidate frequent on a
small host and clean on the large one passes, and the delivered panels differ:
pooling the three hosts reached only 73.6 to 92.4 percent of the
hg38-only design's occupancy-weighted density on hg38, and the smallest host
became the worst one
(`docs/validation/2026-10-02-diversity-baseline.md`, Q4 and Q5).

`background_aggregate` is that choice, and this module is the one place its
vocabulary is defined and its value resolved.

- `sum`: every listed host present at once. Today's behaviour, and the default,
  so no existing design moves.
- `worst-case`: one host out of the panel, each host judged on its own length.
  Reachable and deliberately NOT the default. Which reading is right for a
  panel of several hosts is the open question the Phase 7 measurement was
  written to answer, and nothing measured so far supports changing the default.

`mean` is not offered. `host_profile.aggregate_loads` accepts it, because that
function serves a panel of modelled profiles, but an average over hosts is
neither "all of them at once" nor "the worst one", and no use for it has been
shown. The vocabulary a user can write down is the one here.

**Wired, not yet acting.** Phase 7a plumbs the key only: the schema entry, the
`parameter` global, the resolved request and its hash, and the run manifest's
`effective_conditions`. No gate and no metric reads it yet, so with the key
absent -- and with it set to `sum` -- every figure is what it was. Phase 7b
applies it to the step-2 background gate and Phase 7c to the step-4 per-host
figures. Both read it through `resolve_background_aggregate` or through
`DesignRequest.background_aggregate`, and nowhere else.
"""

from typing import Any

__all__ = [
    "BACKGROUND_AGGREGATES",
    "SUM",
    "WORST_CASE",
    "resolve_background_aggregate",
]

#: Every host present at once: the summed counts over the summed length that
#: `filter` and `base_optimizer` already compute.
SUM = "sum"

#: One host out of the panel, each judged on its own length.
WORST_CASE = "worst-case"

#: The values the key accepts, in the order a refusal names them.
BACKGROUND_AGGREGATES = (SUM, WORST_CASE)


def resolve_background_aggregate(value: Any = None) -> str:
    """Which background aggregate applies, validated.

    This is the seam. There are two routes to the value and one definition of
    it:

    - code that holds a resolved request reads
      `design_request.DesignRequest.background_aggregate`, which this function
      produced at resolution time;
    - code running under the `parameter` globals, which is everything inside a
      step, calls this with no argument and gets the same answer.

    Nothing else reads the key. A stage reaching into `parameter` for it would
    carry its own fallback, which is the inert-option defect of Known Issue 8
    one layer down: the key would appear to work while one of its readers
    quietly disagreed.

    An absent key resolves to `sum`, which is the shipped behaviour, so an
    absent key and an explicit `sum` are one design. An unrecognised value
    raises rather than falling back: pooling hosts for a user who asked for a
    worst case would report a specificity the panel does not have.
    """
    from .exceptions import InvalidDesignRequest

    if value is None:
        from . import parameter

        value = getattr(parameter, "background_aggregate", None)
    if value is None:
        return SUM
    text = str(value)
    if text not in BACKGROUND_AGGREGATES:
        raise InvalidDesignRequest(
            "background_aggregate",
            "must be one of "
            + ", ".join(BACKGROUND_AGGREGATES)
            + " ('mean' is deliberately not offered: an average over hosts is "
            "neither every host at once nor the worst one)",
            text,
        )
    return text
