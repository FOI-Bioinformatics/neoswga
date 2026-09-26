"""Count a KNOWN k-mer set in a reference, without building the reference's table.

`kmer_tables.counts_for` answers the same question from a counted table. This
answers it from the FASTA, and the difference is what it costs to ask:

  counted table   one pass over the reference to build it, then set operations
                  per query batch. The table is kept: 138 MB for hg38 at k=12,
                  and bounded by the reference rather than by the k-mer space
                  above about k=15 -- hg38 at k=18 is roughly 78 to 84 GB as
                  text. Every later batch is then cheap.
  this module     one pass over the reference per BATCH, and nothing written
                  to disk at all.

**Counting is usually faster, and that is not what this is for.** Measured on
Drosophila at k=12 with 2,000 queries, jellyfish counts in 0.9 s and answers
the batch in 0.2 s, against 4.0 s to scan. The scan wins on what it does NOT
need: no counter installed, no table built, and no disk for one. That is the
whole case for it, and it is a real one only where the table is the problem --
a host at large k, where the table is tens of gigabytes, or a reference
queried once.

Memory is bounded by the largest RECORD and by the chunk, not by the query
set: the scan holds one record as bytes and again as codes, so Drosophila's
32 Mb chromosome dominates its 332 MB peak. A reference of many small records
costs little; one huge chromosome costs about 2 bytes per base. An earlier
draft of this file claimed memory proportional to the query set, which the
first measurement refuted.

Agreement is the point, so the rules that decide a count are the ones the rest
of this package already uses:

  canonical       a k-mer and its reverse complement are one quantity, counted
                  under the smaller 2-bit encoding. Both counters here run
                  canonical (`-C` for jellyfish; KMC's default).
  record-aware    a k-mer never spans the join between two FASTA records.
                  Records are processed one at a time, so a join cannot
                  produce one. This is the rule `string_search.spans_a_record_join`
                  enforces for binding positions, and the defect it was written
                  for -- fabricated sites at every contig join -- is the same
                  defect here.
  ambiguity       a window containing any base other than A, C, G or T is not
                  counted, rather than being counted as though the ambiguous
                  base were an A.

Verified against KMC on the whole human genome: at k=12, 921,263 of 921,553
queried wMel k-mers are present, 774,238,630 occurrences, with the present set
and every individual count agreeing exactly.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Sequence

import numpy as np

logger = logging.getLogger(__name__)

#: Positions scanned per pass. Chunks overlap by k-1 so a window spanning a
#: chunk edge is counted exactly once.
#:
#: Measured on Drosophila at k=12, where a larger chunk is both slower and
#: heavier -- 8,000,000 costs 6.6 s and 1,194 MB against 3.9 s and 351 MB
#: here, because the working set stops fitting in cache. Below about 60,000
#: the time stops improving and per-chunk overhead starts to show, so this
#: sits at the flat end of the measured range rather than at its edge.
DEFAULT_CHUNK = 125_000

#: Two bits per base caps the encodable length. 31 keeps `2 * k` at 62 bits,
#: inside uint64 with room for the shift, and the schema admits at most k=30.
MAX_ENCODABLE_K = 31

#: Maps a base to its 2-bit code; every other byte maps to this sentinel,
#: which marks a position no k-mer may span.
_AMBIGUOUS = 4
_BASE_CODE = np.full(256, _AMBIGUOUS, dtype=np.uint8)
for _code, _base in enumerate(b"ACGT"):
    _BASE_CODE[_base] = _code
    _BASE_CODE[_base + 32] = _code  # the same base written in lower case


def encode(kmer: str) -> int | None:
    """The 2-bit code of one k-mer, or None if it holds a base that is not ACGT.

    None rather than a substitute code: a k-mer with an ambiguous base occurs
    nowhere, and giving it the code of some other k-mer would report that
    other k-mer's count under its name.
    """
    code = 0
    for base in kmer.encode("ascii", "ignore"):
        value = int(_BASE_CODE[base])
        if value == _AMBIGUOUS:
            return None
        code = (code << 2) | value
    return code


def canonical_code(kmer: str) -> int | None:
    """The smaller of a k-mer's code and its reverse complement's."""
    forward = encode(kmer)
    if forward is None:
        return None
    k = len(kmer)
    reverse = 0
    remaining = forward
    for _ in range(k):
        reverse = (reverse << 2) | (3 - (remaining & 3))
        remaining >>= 2
    return min(forward, reverse)


def iter_records(path: str) -> Iterator[bytes]:
    """Each FASTA record's sequence, one at a time.

    Yielding per record is what keeps a k-mer from spanning a join: the
    concatenated form this package once scanned produced up to k-1 k-mers per
    join that occur in no record at all.
    """
    sequence = bytearray()
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if sequence:
                    yield bytes(sequence)
                # A new object rather than `sequence.clear()`: the buffer just
                # yielded is still the caller's, and clearing would empty it
                # underneath them.
                sequence = bytearray()
            else:
                sequence += line.strip()
    if sequence:
        yield bytes(sequence)


def _count_codes(
    genome: str, k: int, queries: np.ndarray, chunk: int = DEFAULT_CHUNK
) -> np.ndarray:
    """Occurrences of each code in `queries`, which must be sorted and unique."""
    counts = np.zeros(len(queries), dtype=np.int64)
    if len(queries) == 0:
        return counts
    mask = np.uint64((1 << (2 * k)) - 1) if k < 32 else np.uint64(0xFFFFFFFFFFFFFFFF)
    last = len(queries) - 1

    for record in iter_records(genome):
        bases = _BASE_CODE[np.frombuffer(record, dtype=np.uint8)]
        start = 0
        while start + k <= len(bases):
            window = bases[start : start + chunk + k - 1]
            n = len(window) - k + 1
            ambiguous = window == _AMBIGUOUS
            # An ambiguous base is encoded as A so the rolling shift stays
            # branch-free, then every window containing one is dropped below.
            defined = np.where(ambiguous, 0, window).astype(np.uint64)
            complement = np.uint64(3) - defined

            forward = np.zeros(n, dtype=np.uint64)
            reverse = np.zeros(n, dtype=np.uint64)
            for offset in range(k):
                forward = ((forward << np.uint64(2)) | defined[offset : offset + n]) & mask
                reverse = (
                    (reverse << np.uint64(2)) | complement[k - 1 - offset : k - 1 - offset + n]
                ) & mask

            spanned = np.concatenate(([0], np.cumsum(ambiguous, dtype=np.int64)))
            usable = (spanned[k : k + n] - spanned[:n]) == 0
            canonical = np.minimum(forward, reverse)[usable]

            slot = np.searchsorted(queries, canonical)
            within = slot <= last
            slot = np.minimum(slot, last)
            hit = within & (queries[slot] == canonical)
            counts += np.bincount(slot[hit], minlength=len(queries))
            start += chunk
    return counts


def count_kmers(
    genome: str, k: int, kmers: Sequence[str], chunk: int = DEFAULT_CHUNK
) -> dict[str, int]:
    """Occurrences of each k-mer in `kmers` within `genome`.

    Every requested k-mer appears in the result, with 0 for one the reference
    does not hold. That is `kmer_tables.counts_for`'s contract and it is the
    important half: a key left out would read as unknown, and an unknown
    background count clears the gate it was supposed to face.

    A k-mer and its reverse complement are the same quantity, so both spellings
    receive the canonical count. `counts_for` answers 0 for the spelling a
    canonical table does not store, which every caller in this package avoids
    by passing k-mers read from such a table. The two therefore agree on every
    query the pipeline makes, and this one is the better answer on the rest.
    """
    if not kmers:
        return {}
    if not 1 <= k <= MAX_ENCODABLE_K:
        raise ValueError(
            f"k={k} cannot be scanned: this encoding holds at most "
            f"{MAX_ENCODABLE_K} bases in a 64-bit code."
        )
    if not os.path.exists(genome):
        raise FileNotFoundError(f"Reference not found, so nothing can be counted: {genome}")

    codes: dict[str, int] = {}
    for kmer in kmers:
        if len(kmer) != k:
            continue
        code = canonical_code(kmer)
        if code is not None:
            codes[kmer] = code

    queries = np.array(sorted(set(codes.values())), dtype=np.uint64)
    counts = _count_codes(genome, k, queries, chunk=chunk)
    by_code = dict(zip(queries.tolist(), counts.tolist(), strict=True))
    return {kmer: by_code.get(codes.get(kmer, -1), 0) for kmer in kmers}


def describe_cost(genome: str, k: int, queries: int) -> str:
    """One line on what a scan of this reference will cost, for a log.

    Quoted from the two hg38 measurements rather than from a model: the pass is
    linear in the reference and nearly flat in the query count, so the figure
    that matters to someone waiting is the reference's size.
    """
    try:
        size = os.path.getsize(genome)
    except OSError:
        size = 0
    # Measured on Drosophila: 144 Mb of FASTA in 4.0 s, about 36 Mb/s. The
    # pass is linear in the reference and nearly flat in the query count, so
    # the reference's size is the figure that matters to someone waiting.
    seconds = size / 36_000_000 if size else 0
    return (
        f"Counting {queries:,} {k}-mers directly in {os.path.basename(genome)} "
        f"({size / 1e9:.2f} Gb), with no table to build and none kept. "
        f"Expect roughly {seconds / 60:.1f} min for this reference."
    )
