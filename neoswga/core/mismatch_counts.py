"""Background site counts grouped by how well each site matches.

Selectivity counts exact k-mer matches only, so the sites that stringency
actually discriminates against -- those carrying a mismatch or two -- are
invisible to the model. Occupancy then has nothing to weight, and an additive
has nothing to act on.

Counting them is nearly free. The jellyfish ``*_all.txt`` files already hold a
count for every k-mer in the genome, so a primer's 1-mismatch background load
is the sum over its ``3k`` neighbours: 30 lookups for a 10-mer, not a rescan of
the genome.

Two details decide whether the resulting numbers mean anything.

**Canonicalisation.** Jellyfish runs with ``-C``, so a count file holds only the
lexicographically smaller of each k-mer and its reverse complement, never both.
A primer and its reverse complement are the same duplex on double-stranded DNA
and have to count the same. Roughly half of any primer's ``3k`` neighbours are
non-canonical, so looking them up as written would drop half the background
load -- and halving the denominator of a selectivity ratio makes a design look
twice as specific as it is. ``filter._get_rate_for_one_file`` has the same
shape (it tests ``parts[0] in primer_set``) and is latent there only because the
default candidate list is itself read out of these files, and so is already
canonical.

**Double counting.** Sums run over distinct canonical forms. Two different
1-mismatch variants can share a canonical form when a primer is near
palindromic, and the exact match is excluded from the mismatch classes so a
site is never weighted at two different occupancies.

Memory: ``load_kmer_counts`` holds one dict per (prefix, k). That is the same
order as the exact path already uses and is fine to bacterial scale; a
host-sized background belongs on the Bloom plus sampled-index route instead,
which answers the same question without the dict.
"""

from collections.abc import Iterable
from functools import lru_cache

from neoswga.core.thermodynamics import reverse_complement

_BASES = ("A", "C", "G", "T")


def canonical_kmer(kmer: str) -> str:
    """The form jellyfish ``-C`` stores: the lexicographic minimum of a k-mer
    and its reverse complement.

    Both describe the same physical duplex, so both must resolve to one count.
    """
    rc = reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def one_mismatch_variants(seq: str) -> set[str]:
    """Every sequence differing from ``seq`` at exactly one position.

    ``3k`` of them, excluding ``seq`` itself. Shared with
    ``background_filter.BackgroundBloomFilter`` rather than reimplemented, so
    the two cannot disagree about what a mismatch neighbour is.
    """
    variants = set()
    for i, original in enumerate(seq):
        for base in _BASES:
            if base != original:
                variants.add(seq[:i] + base + seq[i + 1 :])
    return variants


def _variants_at_distance(seq: str, n_mismatches: int) -> set[str]:
    """Sequences differing from ``seq`` at exactly ``n_mismatches`` positions.

    Grown one mismatch at a time and filtered by actual Hamming distance,
    because expanding neighbours of neighbours also reaches sequences closer
    than ``n`` (changing a base back, or two edits landing on one position).
    """
    if n_mismatches <= 0:
        return set()

    frontier = {seq}
    for _ in range(n_mismatches):
        grown = set()
        for candidate in frontier:
            grown |= one_mismatch_variants(candidate)
        frontier = grown

    return {v for v in frontier if _hamming(seq, v) == n_mismatches}


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b, strict=True) if x != y)


@lru_cache(maxsize=32)
def load_kmer_counts(prefix: str, k: int) -> dict[str, int]:
    """Canonical k-mer -> count, from a jellyfish ``{prefix}_{k}mer_all.txt``.

    Cached per (prefix, k): mismatch counting makes tens of lookups per primer,
    which is the wrong shape for the line-scan-per-batch the exact path uses.

    Raises rather than returning an empty dict when the file is missing. A zero
    from an absent file is indistinguishable from zero background binding, and
    the second is a strong claim to make from a missing input.
    """
    # Through the table layer, so a KMC database answers as well as a text
    # dump. This still materialises the whole table as a dict, which is its own
    # ceiling on a host genome; asking only for each candidate's neighbours
    # (`kmer_tables.counts_for`) is the change that removes it, and it is a
    # restructuring of the callers rather than of this function.
    from neoswga.core import kmer_tables

    return dict(kmer_tables.iter_table(prefix, k))


def _count_of(kmers: Iterable[str], tables: list[dict[str, int]]) -> int:
    """Total count over distinct canonical forms, across all genomes."""
    total = 0
    for canonical in {canonical_kmer(kmer) for kmer in kmers}:
        for table in tables:
            total += table.get(canonical, 0)
    return total


def mismatch_class_counts(
    primer: str, prefixes: list[str], max_mismatches: int = 1
) -> dict[int, int]:
    """Binding-site counts for ``primer``, grouped by mismatch count.

    Args:
        primer: Primer sequence.
        prefixes: K-mer file prefixes to sum over. Multiple background genomes
            contribute additively -- a primer binding two hosts carries both
            loads.
        max_mismatches: Highest class to report. 0 reproduces exact-match
            counting, which is the reduction the occupancy model is checked
            against.

    Returns:
        ``{mismatch_count: sites}`` for 0..max_mismatches. Classes are
        disjoint: the exact match never appears in class 1.
    """
    k = len(primer)
    tables = [load_kmer_counts(prefix, k) for prefix in prefixes]

    counts = {0: _count_of([primer], tables)}
    for distance in range(1, max_mismatches + 1):
        variants = _variants_at_distance(primer, distance)
        # A variant can canonicalise onto the primer itself for a palindromic
        # sequence; dropping it keeps the classes disjoint.
        primer_canonical = canonical_kmer(primer)
        variants = {v for v in variants if canonical_kmer(v) != primer_canonical}
        counts[distance] = _count_of(variants, tables)

    return counts
