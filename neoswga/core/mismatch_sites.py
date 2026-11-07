"""Mismatch neighbours of a primer, kept individually instead of counted by class.

``mismatch_counts.mismatch_class_counts`` enumerates every Hamming neighbour of
a primer, looks up each one's count in a k-mer table, and then collapses the
result to one total per mismatch distance. Position and identity are discarded
before any thermodynamics runs, so no penalty applied downstream can depend on
them -- which is why ``occupancy.mismatch_tm`` can only be a uniform number per
mismatch, and why its own docstring records that the uniform form is wrong in
the direction that matters.

This module keeps the neighbour. The enumeration already happens and the count
lookup is already one per neighbour, so the counting cost is unchanged; what
changes is that the counts are not summed before anything weights them.

**Canonicalisation is not re-derived here.** ``mismatch_counts.canonical_kmer``
is the one definition of the form a table stores, and
``mismatch_counts._variants_at_distance`` the one definition of a neighbour.
Re-deriving either is how a neighbour comes to be counted twice, or half of a
background load to be dropped: roughly half of any primer's ``3k`` neighbours
are non-canonical, and looking them up as written would halve the denominator
of a selectivity ratio.

**One canonical form is one site group, returned once.** Two distinct
neighbours share a canonical form when one is the reverse complement of the
other, which happens for near-palindromic primers. Those are two physically
different duplexes and a site stored under that form could be either, so the
group carries both ``readings`` and the count once. Summing ``count`` over the
returned groups therefore cannot double count, whatever a scoring function does
with the readings.

**What this module does NOT decide.** It reports counts and geometry. It applies
no thermodynamics, no penalty and no 3'-end rule; those live in
``mismatch_model`` so that a change to the model cannot quietly change what was
counted.

**Depth.** Groups are deduplicated by canonical form exactly as
``mismatch_class_counts`` deduplicates -- within a distance AND across them --
so a load computed from them reduces to the uniform one at any depth. A form
reached at two distances is one physical duplex, so it yields ONE group, at the
nearest distance that reaches it, carrying every reading that does. The caller
sums ``count * weight`` over groups, so two groups would be two sites; and
``_group_weight`` takes the most stable reading rather than summing them, so
carrying the farther reading costs nothing and keeps the "a site binds in
whichever orientation binds better" rule intact.

Both halves are needed. Within a distance, a variant and its reverse complement
can both be neighbours of a near-palindromic primer. Across distances, the same
stored k-mer can be a 1-mismatch neighbour read one way and a 2-mismatch
neighbour read the other, because a count file holds only one of each k-mer and
its reverse complement. Until 2026-11-01 this module documented the second case
as a known limitation inherited from the uniform path; the uniform path was
fixed the same day, so leaving it would have made the two disagree -- which is
the "two implementations of one quantity" failure this project keeps paying for.
Reachable only above the shipped default depth of 1.
``docs/validation/2026-11-01-mismatch-classes-across-strands.md``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from neoswga.core.mismatch_counts import (
    _variants_at_distance,
    canonical_kmer,
    load_kmer_counts,
)

__all__ = ["NeighbourReading", "NeighbourSite", "neighbour_sites"]


@dataclass(frozen=True)
class NeighbourReading:
    """One way of reading a site group as a duplex with the primer.

    A group has two readings only when the primer is near palindromic enough
    that a neighbour and its reverse complement are both neighbours. The two
    are different duplexes with different mismatch geometry, and nothing here
    chooses between them.
    """

    #: The neighbour in the primer's own orientation and index order, so
    #: ``primer[i]`` and ``neighbour[i]`` are the two bases of one duplex
    #: position.
    neighbour: str
    #: Mismatch offsets measured from the 3' end, ascending. 0 is the terminal
    #: base. Empty for the exact match.
    offsets_from_three_prime: tuple[int, ...]
    #: ``(primer base, neighbour base)`` at each mismatch, in the same order as
    #: ``offsets_from_three_prime``.
    pairs: tuple[tuple[str, str], ...]

    @property
    def is_three_prime_terminal(self) -> bool:
        """Whether a mismatch sits on the 3'-terminal base itself."""
        return 0 in self.offsets_from_three_prime

    @property
    def is_at_a_duplex_terminus(self) -> bool:
        """Whether a mismatch sits on the first or last base of the duplex.

        Both ends, not only the 3' one. A nearest-neighbour INTERNAL mismatch
        parameter set does not describe either terminus: a terminal position has
        one flanking stack instead of two, and dangling-end and terminal-mismatch
        parameters are the applicable set. ``mismatch_model`` asks this so it can
        decline to apply an internal table out of its domain.
        """
        if not self.offsets_from_three_prime:
            return False
        last = len(self.neighbour) - 1
        return 0 in self.offsets_from_three_prime or last in self.offsets_from_three_prime


@dataclass(frozen=True)
class NeighbourSite:
    """A canonical k-mer the primer binds, and how many sites carry it.

    ``count`` is a measured number of sites: the table holds a count for every
    k-mer present in the genome, so a canonical form absent from it occurs zero
    times. A table that could not be read is a different thing and raises
    before any of these are built.
    """

    #: The form the table stores; the identity of the site group.
    canonical: str
    #: Sites carrying this canonical form, summed over the prefixes asked for.
    count: int
    #: Hamming distance from the primer. 0 is the exact match.
    distance: int
    #: Every duplex this group can form with the primer. At least one.
    readings: tuple[NeighbourReading, ...]


def _read(primer: str, neighbour: str) -> NeighbourReading:
    """Offsets and base pairs for one neighbour, measured from the 3' end.

    From the 3' end because that is the end a polymerase extends from, and the
    only index in which a 3'-end rule can be stated without also knowing the
    oligo length.
    """
    last = len(primer) - 1
    found: list[tuple[int, tuple[str, str]]] = []
    for index, (expected, actual) in enumerate(zip(primer, neighbour, strict=True)):
        if expected != actual:
            found.append((last - index, (expected, actual)))
    found.sort(key=lambda item: item[0])
    return NeighbourReading(
        neighbour=neighbour,
        offsets_from_three_prime=tuple(offset for offset, _ in found),
        pairs=tuple(pair for _, pair in found),
    )


def _count_over(canonical: str, tables: Sequence[Mapping[str, int]]) -> int:
    """Sites carrying ``canonical``, summed over every genome asked for.

    Additive across prefixes, matching ``mismatch_counts._count_of``: a primer
    binding two hosts carries both loads.
    """
    return sum(table.get(canonical, 0) for table in tables)


def neighbour_sites(
    primer: str,
    prefixes: Sequence[str],
    max_mismatches: int = 1,
    *,
    tables: Sequence[Mapping[str, int]] | None = None,
) -> tuple[NeighbourSite, ...]:
    """Every site group ``primer`` binds within ``max_mismatches``, kept apart.

    Args:
        primer: Primer sequence, 5' to 3'.
        prefixes: K-mer table prefixes to sum counts over.
        max_mismatches: Highest distance to enumerate. 0 gives the exact match
            alone, which is the reduction the model is checked against.
        tables: Already-loaded count tables, for a caller scoring many primers
            against one genome set. Loaded from ``prefixes`` when omitted.

    Returns:
        One ``NeighbourSite`` per distinct canonical form, ordered by distance
        and then by canonical form. Groups with a zero count are included: zero
        sites is a measurement, and dropping them would make the reported
        geometry depend on the genome in a way a caller cannot see.

    Raises:
        FileNotFoundError, OSError: propagated from ``load_kmer_counts`` when a
            table is missing. A neighbour whose count cannot be looked up makes
            the whole load unavailable, and an absent table read as zero sites
            is the strongest claim this model could make from a missing input.
    """
    if tables is None:
        tables = [load_kmer_counts(prefix, len(primer)) for prefix in prefixes]

    primer_canonical = canonical_kmer(primer)
    groups: list[NeighbourSite] = [
        NeighbourSite(
            canonical=primer_canonical,
            count=_count_over(primer_canonical, tables),
            distance=0,
            readings=(_read(primer, primer),),
        )
    ]

    # One group per canonical form, at the NEAREST distance that reaches it,
    # carrying every reading that does. Merged across distances and not only
    # within one, because a form reached at two distances is still one physical
    # duplex and the caller sums `count * weight` over groups.
    #
    # The primer's own form is dropped for the same reason: class 0 already
    # holds it.
    by_canonical: dict[str, tuple[int, list[str]]] = {}
    for distance in range(1, max(0, int(max_mismatches)) + 1):
        for variant in _variants_at_distance(primer, distance):
            canonical = canonical_kmer(variant)
            if canonical == primer_canonical:
                continue
            found = by_canonical.get(canonical)
            if found is None:
                by_canonical[canonical] = (distance, [variant])
            else:
                found[1].append(variant)

    for canonical, (distance, variants) in sorted(
        by_canonical.items(), key=lambda item: (item[1][0], item[0])
    ):
        groups.append(
            NeighbourSite(
                canonical=canonical,
                count=_count_over(canonical, tables),
                distance=distance,
                readings=tuple(_read(primer, variant) for variant in sorted(variants)),
            )
        )

    return tuple(groups)
