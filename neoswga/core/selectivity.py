"""Selectivity arithmetic: binding loads into a ratio, and into a density.

Extracted from `base_optimizer` on 2026-09-30, when adding one configuration
field pushed that module past its size budget. Nothing moved but the text: both
functions and both constants are re-exported from `base_optimizer`, so every
existing importer reaches the same objects.

They belong together because they answer the same question two ways, and the
difference between them is the subject of Known Issue 6: the count ratio moves
with how much background sequence the caller supplied, and the density ratio
does not.
"""

# Selectivity ratio that scores a full 1.0 on the normalised scale. Not a
# threshold for a usable design -- it is the point past which more selectivity
# stops earning score, chosen so the term discriminates across the range real
# designs occupy rather than saturating inside it.
SELECTIVITY_REFERENCE = 100.0

# Stands in for unbounded selectivity when no background binding is
# detected at all, so the value stays finite and serialisable. Same
# convention as multi_genome_filter.MAX_ENRICHMENT.
MAX_SELECTIVITY = 1e6


def selectivity_from_loads(fg_load: float, bg_load: float) -> float:
    """Selectivity from foreground and background binding loads.

    The integer version guarded with `max(bg, 1)`, which is right for counts --
    you cannot have a fraction of a site. Occupancy-weighted loads are floats,
    and carrying that guard over silently divided by 1.0 whenever the
    background load fell below it: a set with 0.3 effective background sites
    reported a third of its true selectivity, and one with 0.0 reported the
    foreground load as though it were a ratio.

    A genuinely zero background load is unbounded selectivity, not a large
    number. `MAX_SELECTIVITY` stands in so the value stays finite and
    JSON-serialisable, following the same convention and for the same reason as
    `multi_genome_filter.MAX_ENRICHMENT`: it means "no background binding was
    detected", not "measured this well".
    """
    if bg_load <= 0:
        return MAX_SELECTIVITY if fg_load > 0 else 0.0
    return fg_load / bg_load


def selectivity_density_from_loads(
    fg_load: float, fg_length: float, bg_load: float, bg_length: float
) -> float:
    """Binding load per base of target against per base of background.

    `selectivity_from_loads` divides two counts and carries no genome length,
    so its value moves with how much background sequence the caller supplied.
    Substituting a whole human genome (3.1 Gb) for the single chromosome often
    used as a stand-in (chr21, 46.7 Mb) multiplies the background load by
    roughly the length ratio and divides the reported selectivity by it, with
    nothing about the primers changed.

    Enrichment does not work that way. Priming is proportional to sites per
    genome copy, and a host contributes its sites spread over its own length,
    so the quantity that survives the substitution -- and the one comparable
    between two designs scored against different backgrounds -- is the ratio of
    densities.

    Both are reported. The count ratio remains what the objective is scored on:
    within a single run the background is fixed, so the two rank sets
    identically and selection is unaffected.
    """
    if fg_length <= 0 or bg_length <= 0:
        # No length is an absent measurement, not a clean background.
        return 0.0

    fg_density = fg_load / fg_length
    bg_density = bg_load / bg_length

    if bg_density <= 0:
        # Same convention as the count ratio: unbounded, not merely large.
        return MAX_SELECTIVITY if fg_density > 0 else 0.0
    return fg_density / bg_density


# The underscored names are what `base_optimizer` exposed and what six test
# files and two modules import. Kept as aliases rather than renaming the
# callers, so this extraction moves no other file.
_selectivity_from_loads = selectivity_from_loads
_selectivity_density_from_loads = selectivity_density_from_loads
