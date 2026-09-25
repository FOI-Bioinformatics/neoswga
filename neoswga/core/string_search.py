import bisect
import logging
import multiprocessing
import os

import h5py

from neoswga.core import kmer_tables, parameter, position_index
from neoswga.core.thermodynamics import reverse_complement

logger = logging.getLogger(__name__)

try:
    import ahocorasick

    AHOCORASICK_AVAILABLE = True
except ImportError:
    AHOCORASICK_AVAILABLE = False
    logger.info(
        "pyahocorasick not installed. Position scanning will use per-k "
        "sequential search. Install for 5-7x speedup: "
        "pip install neoswga[fast]"
    )

# =============================================================================
# Genome Sequence Cache (eliminates redundant file reads)
# =============================================================================

# Module-level cache for genome sequences
_genome_cache: dict[str, str] = {}

# Offsets in the concatenated sequence at which a new FASTA record begins,
# excluding 0. Records are joined without a separator -- the concatenation is
# the coordinate system every stored position is expressed in -- so a k-mer can
# straddle a join and be reported as present when neither record contains it.
# These offsets are what `get_all_positions_multi_k` rejects such a match with.
_record_boundary_cache: dict[str, list[int]] = {}


def get_cached_record_boundaries(seq_fname: str) -> list[int]:
    """Record-start offsets for a genome, loading it if it is not cached."""
    if seq_fname not in _record_boundary_cache:
        get_cached_genome_sequence(seq_fname)
    return _record_boundary_cache.get(seq_fname, [])


def get_cached_genome_sequence(seq_fname: str) -> str:
    """
    Get genome sequence with caching.

    Caches the full genome string to avoid re-reading for each k-mer length.
    Provides 5-10x speedup when processing multiple k-mer lengths.

    Args:
        seq_fname: Path to the FASTA file.

    Returns:
        The genome sequence as a string.
    """
    if seq_fname in _genome_cache:
        return _genome_cache[seq_fname]

    logger.info(f"Loading genome sequence from {seq_fname}...")

    # Read the genome into memory one record at a time.
    #
    # NOT `"".join(utility.read_fasta_file(...))`: that helper yields one
    # character per base, and `str.join` materialises its argument into a list
    # first, so the list holds one pointer per base -- about 26 GB for hg38
    # (3.3 Gbp) before a single base is concatenated. The symptom was a
    # SIGKILL with "Loading genome sequence from ..." as the last log line.
    #
    # NOT `"".join(loader.load_genome_streaming(...))` either, which replaced
    # it. `str.join` calls `PySequence_Fast` on its argument, so a generator is
    # materialised into a list of every record before any concatenation begins:
    # every chromosome string is alive at the moment the finished sequence is
    # allocated. For hg38 that is roughly 3.3 GB of records plus a 3.3 GB
    # result, which is the recorded 8.5 GB peak and the reason hg38 filter runs
    # had to be serialised.
    #
    # Appending in a loop is the one construction that never holds more than
    # one previously-yielded record alongside the accumulator, which is what
    # `tests/test_genome_cache_join_granularity.py` pins. CPython grows a
    # string in place when the accumulator's refcount is 1, so this is O(n) in
    # time and its peak is roughly the finished sequence plus one record.
    #
    # Measured on this machine (CPython 3.11.14, macOS/arm64) against the join
    # it replaces, on a synthetic 960 Mb genome in 24 records:
    #
    #                     time      tracemalloc peak   process max RSS
    #     join            7.68 s          1927 MB           2183 MB
    #     append loop     8.03 s          1168 MB           1464 MB
    #
    # About 40 percent less peak memory by both measures, for a time difference
    # inside run-to-run variation. Five replicates each, in separate processes:
    # every join run reported 1927.0 MB and every loop run 1167.6 MB, so the
    # allocation figures are deterministic even though the wall times are not.
    #
    # Holding the GENOME at 96 Mb and varying the record count, the loop's peak
    # falls from 123 MB at 24 records to 106 MB at 192, approaching the genome
    # size, while the join stays flat at 199 MB, which is two copies. That is
    # the signature of the in-place growth path being taken.
    #
    # `sequence` MUST stay a function-local. CPython's in-place append is
    # implemented for `STORE_FAST`/`LOAD_FAST`, so it does not fire for a
    # module-level name, which goes through the globals dict instead. Move this
    # accumulator to module scope and every `+=` becomes a full copy: the same
    # code turns quadratic, silently, with no test failing.
    #
    # Two independent attempts to measure this reached the opposite conclusion
    # by two different routes, which is why both cautions are recorded here. One
    # used `resource.getrusage`, which reports a high-water mark for the whole
    # process lifetime, so timing both constructions in one process makes
    # whichever runs second inherit the first one's peak and can reverse the
    # result. The other hand-wrote the loop at module scope in a benchmark
    # script, and so measured the quadratic path rather than the one that ships.
    # Measure them in separate
    # processes, and prefer `tracemalloc`, which attributes allocations rather
    # than reporting resident pages the allocator has not returned.
    #
    # `load_genome_streaming` already uppercases each record, which is why no
    # `.upper()` follows: that call was another full-length copy of the genome.
    from neoswga.core import genome_io

    loader = genome_io.GenomeLoader()
    sequence = ""
    boundaries: list[int] = []
    for record in loader.load_genome_streaming(seq_fname):
        if sequence:
            boundaries.append(len(sequence))
        sequence += record

    _genome_cache[seq_fname] = sequence
    _record_boundary_cache[seq_fname] = boundaries
    logger.info(f"Cached genome sequence: {len(sequence):,} bp")

    return sequence


def preload_genomes(seq_fnames: list[str]) -> None:
    """
    Preload multiple genome sequences into cache.

    Call this before processing to eliminate I/O during computation.

    Args:
        seq_fnames: List of FASTA file paths to preload.
    """
    for fname in seq_fnames:
        if fname not in _genome_cache:
            get_cached_genome_sequence(fname)
    logger.info(f"Preloaded {len(_genome_cache)} genome(s) into cache")


def clear_genome_cache() -> None:
    """Clear the genome sequence cache to free memory."""
    global _genome_cache, _record_boundary_cache
    cache_size = len(_genome_cache)
    _genome_cache = {}
    _record_boundary_cache = {}
    logger.info(f"Cleared genome cache ({cache_size} entries)")


def get_genome_cache_stats() -> dict[str, int]:
    """Return statistics about the genome cache."""
    total_bp = sum(len(seq) for seq in _genome_cache.values())
    return {
        "num_genomes": len(_genome_cache),
        "total_bp": total_bp,
        "memory_mb": total_bp // (1024 * 1024),
    }


# pyahocorasick's `Automaton.iter()` indexes with a 32-bit int. Handed a string
# longer than 2**31 - 1 its scan loop never executes: it yields nothing, raises
# nothing and warns nothing. Measured with the needle planted at offset 1000,
# length 2**31 - 10 finds it and 2**31 + 10 does not.
#
# Human is 3.1 Gb and mouse 2.7 Gb, so for the hosts this tool exists to
# discriminate against, every primer came back with an empty position list and
# `total_bg_sites` / `bg_coverage` reported 0 -- indistinguishable, downstream,
# from a perfectly specific primer set. Scanning below the limit is the fix.
MAX_SCAN_CHUNK = 2**30


def spans_a_record_join(start: int, length: int, boundaries) -> bool:
    """Whether a match at `start` exists only because two records were joined.

    Records are concatenated with no separator, so the last k-1 bases of one
    record and the first bases of the next form k-mers that occur in neither.
    On a draft assembly every join manufactures up to k-1 of them, and each one
    becomes a binding site for a primer that does not occur there: a foreground
    site inflates coverage, a background site deflates specificity, and nothing
    downstream can tell them from real ones.

    `bisect_right` gives the first boundary strictly after the match start. A
    match that begins exactly ON a boundary starts a record and is real; one
    that contains a boundary strictly inside it spans the join.

    Extracted on 2026-09-21 because the two scanners disagreed.
    `get_all_positions_multi_k` applied this rule and
    `get_all_positions_per_k` did not, so the same reference produced different
    site sets depending on whether pyahocorasick was installed. Two
    implementations of one quantity is how this codebase has produced
    disagreeing coverage numbers before.
    """
    if not boundaries:
        return False
    nxt = bisect.bisect_right(boundaries, start)
    return nxt < len(boundaries) and boundaries[nxt] < start + length


def get_all_positions_multi_k(primer_lists_by_k, seq_fname, circular, chunk_size=None):
    """Scan genome once for primers of all k-values using Aho-Corasick.

    Requires the optional ``pyahocorasick`` package.  When available this
    replaces per-k scanning and provides ~5-7x speedup for default k
    ranges (6-12).

    The genome is scanned in chunks of at most ``MAX_SCAN_CHUNK`` bases, for
    the 32-bit reason described above the constant. Consecutive chunks overlap
    by ``k - 1`` so that a match straddling a boundary is still seen whole, and
    each match is attributed to the chunk its START falls in, so the overlap
    cannot double-count it. ``chunk_size`` exists so tests can exercise the
    boundary logic without allocating two gigabytes.

    Args:
        primer_lists_by_k: Dict mapping k -> list of primers.
        seq_fname: Path to genome FASTA.
        circular: Whether the genome is circular.
        chunk_size: Override the scan chunk length. Tests only.

    Returns:
        Dict mapping primer -> list of start positions.
    """
    A = ahocorasick.Automaton()
    all_primers = {}
    for k, primers in primer_lists_by_k.items():
        for primer in primers:
            A.add_word(primer, primer)
            all_primers[primer] = []
    A.make_automaton()

    sequence = get_cached_genome_sequence(seq_fname)
    boundaries = get_cached_record_boundaries(seq_fname)
    seq_len = len(sequence)
    max_k = max(primer_lists_by_k.keys())

    if circular:
        search_seq = sequence + sequence[: max_k - 1]
    else:
        search_seq = sequence

    limit = chunk_size or MAX_SCAN_CHUNK
    overlap = max_k - 1
    total = len(search_seq)

    start = 0
    while start < total:
        stop = min(start + limit, total)
        window = search_seq[start : min(stop + overlap, total)]
        for end_pos, primer in A.iter(window):
            abs_start = start + end_pos - len(primer) + 1
            # Attribute to the chunk owning the start; the overlap is read but
            # never claimed twice.
            if not (start <= abs_start < stop):
                continue
            if abs_start >= seq_len:
                continue
            if spans_a_record_join(abs_start, len(primer), boundaries):
                continue
            all_primers[primer].append(abs_start)
        start = stop

    return all_primers


# #everything should be 5' to 3' written
def get_all_positions_per_k(kmer_list, seq_fname, circular, fname_prefix=None):
    """
    Gets all the positions of k-mers for one value of k. It assumes the first kmer in the list is the
    same length as all kmers in the list. Everything is written in the 5' to 3' direction.

    Now uses genome caching for 5-10x speedup when processing multiple k-mer lengths.

    Args:
        kmer_list: The list of all kmers where k is a specific and single value.
        seq_fname: The path to the fasta file of the genome.
        fname_prefix: This argument is only necessary for knowing which fasta file is currently being processed.

    Returns:
        kmer_dict: Dictionary of kmer to count in the fasta file seq_fname.

    """
    if len(kmer_list) == 0:
        return {}
    k = len(kmer_list[0])

    # Initialize result dictionary
    kmer_dict = {kmer: [] for kmer in kmer_list}

    # Use cached genome sequence (avoids re-reading file for each k)
    sequence = get_cached_genome_sequence(seq_fname)

    if fname_prefix is not None:
        print("Starting the search for " + fname_prefix + " " + str(k) + "mers...")

    # Optimized sliding window search using cached sequence
    seq_len = len(sequence)
    search_len = seq_len if not circular else seq_len + k - 1
    boundaries = get_cached_record_boundaries(seq_fname)

    for i in range(search_len - k + 1):
        if i < seq_len - k + 1:
            # Normal position within sequence
            current_kmer = sequence[i : i + k]
        else:
            # Circular wrap-around: combine end with beginning
            wrap_pos = i - (seq_len - k + 1)
            current_kmer = sequence[seq_len - k + 1 + wrap_pos :] + sequence[: wrap_pos + 1]

        if current_kmer not in kmer_dict:
            continue
        # The same rule the Aho-Corasick path applies. Without it this scanner
        # reported k-mers that occur in no record: on a two-record fixture
        # `ACGGTA` is absent from both and was stored at the join.
        if spans_a_record_join(i, k, boundaries):
            continue
        kmer_dict[current_kmer].append(i)

    return kmer_dict


# Dataset name holding the record-start offsets for a prefix. '#' cannot occur
# in a primer, so this can never collide with a k-mer key. Defined in
# `position_index`, which owns every name inside the file.
RECORD_STARTS_KEY = position_index.RECORD_STARTS_KEY

# Bumped when the on-disk geometry an index carries changes shape.
#
# Version 1 is the first to store record starts, so an index without this
# attribute predates record-aware scanning and its coverage windows crossed
# contig boundaries.
#
# Version 2 (2026-09-21) is the first whose SITES are record-aware.
# `get_all_positions_per_k` scanned the concatenated sequence without applying
# the join rule its Aho-Corasick counterpart already had, so on a multi-record
# reference it stored k-mers that occur in no record -- up to k-1 fabricated
# sites per join. Nothing in the file records which scanner wrote it, so a
# version 1 index cannot be vouched for and is refused rather than read.
# Single-record references are unaffected either way: with no joins the two
# scanners always agreed.
INDEX_FORMAT_VERSION = 2


def write_to_h5py(kmer_dict, fname_prefix, replace=False, record_starts=None, genome_fname=None):
    """
    Writes the kmer counts to an h5py file, which allows for efficient access in terms of looking up the
    frequency of a particular k-mer. If the kmer already exists in the dataset, the entry in the h5py file
    is overwritten with the new data.

    Args:
        kmer_dict: A dictionary of all the k-mers and their respective frequencies.
        fname_prefix: The file path prefix for the output h5py file. If k=6, for example, '_6mer_positions.h5'
        will be appended to the file name.
        replace: Discard whatever the file already holds instead of merging into
            it. The reuse gate sets this when a file's provenance record does not
            match the current run, so the datasets it holds for primers outside
            this run's pool -- which were scanned from something else -- do not
            survive. The old content is dropped here, at the moment the
            replacement is in hand, rather than when the mismatch is detected: an
            interrupted run then leaves the previous file intact rather than an
            empty one, which the step-4 prerequisite check catches on the
            foreground but not on the background.
    """
    if not kmer_dict:
        return  # Nothing to write

    k = len(next(iter(kmer_dict.keys())))
    h5_path = position_file_path(fname_prefix, k)
    attrs = {}
    if record_starts is not None:
        # Format version and reference identity, so a later design can tell
        # a modern index from a concatenation-era one and can tell which
        # reference it was built from. Stored as root attributes, which no
        # reader can mistake for a k-mer.
        attrs["index_format_version"] = INDEX_FORMAT_VERSION
        if genome_fname is not None:
            from neoswga.core import kmer_counter

            try:
                attrs["reference_digest"] = kmer_counter.genome_fingerprint(genome_fname)
            except OSError as exc:  # pragma: no cover - unreadable reference
                logger.debug("Could not fingerprint %s: %s", genome_fname, exc)
    # Written in the sorted-blocks layout (`core/position_index.py`). The
    # per-dataset layout spent most of its bytes on HDF5 object headers and
    # was rewritten entry by entry in place; this writes the whole index to
    # a new file and moves it into place, so an index converted from the old
    # layout, or merged with this run's entries, is never left half written.
    #
    # Record starts say where each FASTA record begins in the concatenated
    # coordinate system every stored position uses. Without them a coverage
    # window anchored near the end of one record extends into the next.
    position_index.write_entries(
        h5_path,
        kmer_dict,
        replace=replace,
        record_starts=None if record_starts is None else list(record_starts),
        attrs=attrs,
    )


def check_which_primers_absent_in_h5py(primer_list, fname_prefix):
    """
    This function checks if the primers in a given list are missing from an h5py file.

    Args:
        primer_list: The list of primers that need to be checked exist in the h5py file.
        fname_prefix: The prefix of the h5py file--basically the path minus '_6mer_positions.h5' where k = 6.

    Returns:
        filtered_primer_list: All the primers from the given list of primers that are missing from the h5py file.
    """
    if not primer_list:
        return primer_list

    k = len(primer_list[0])
    h5_path = fname_prefix + "_" + str(k) + "mer_positions.h5"

    if not os.path.exists(h5_path):
        # Create empty HDF5 file using context manager
        with h5py.File(h5_path, "a"):
            pass
        return primer_list

    with position_index.open_index(h5_path) as index:
        keys = set(index.keys())

    # Which of the not-yet-indexed primers occur in this genome at all. This
    # used to read EVERY k-mer in the genome into a Python set to answer a
    # question about a known, much smaller list -- 8.4 million entries for
    # hg38 at k=12, and billions at k=16 or above, where a host genome cannot
    # hold a text table in the first place. `counts_for` asks only about the
    # primers, which is a bounded query whatever the genome's size.
    unindexed = [primer for primer in primer_list if primer not in keys]
    if not unindexed:
        return []
    present = kmer_tables.counts_for(fname_prefix, k, unindexed)
    return [primer for primer in unindexed if present.get(primer, 0) > 0]


# Genome fingerprints already computed in this process, keyed by path, size and
# mtime. See `_position_provenance`.
_fingerprint_cache: dict[tuple, str] = {}


def position_file_path(fname_prefix, k):
    """Where the binding positions for one prefix and one k live."""
    return fname_prefix + "_" + str(k) + "mer_positions.h5"


def _position_provenance(genome_fname, k, circular):
    """What a position file has to have been produced under to be reusable.

    Genome identity plus the scan parameters that change the answer. `circular`
    belongs here because `get_all_positions_multi_k` appends `sequence[:max_k-1]`
    before scanning, so a primer straddling the origin is found only when it is
    true: a pool scanned once with `circular=False` and reused under
    `circular=True` silently loses those sites.

    The fingerprint is `kmer_counter.genome_fingerprint`, the same cheap
    size-plus-first-and-last-megabyte identifier the k-mer tables use, rather
    than a second mechanism. A full hash of hg38 costs seconds on an invocation
    that usually skips.

    Returns None when the genome cannot be read, which is treated as
    unverifiable and therefore not reusable.
    """
    from neoswga.core import kmer_counter

    try:
        stat = os.stat(genome_fname)
        key = (os.path.abspath(genome_fname), stat.st_size, stat.st_mtime_ns)
        # Memoised because the record is built once per k and again on the way
        # out, so a five-k run would otherwise read two megabytes of the genome
        # fourteen times per prefix. The size and mtime are part of the key, so
        # a file replaced under a running process is fingerprinted again.
        fingerprint = _fingerprint_cache.get(key)
        if fingerprint is None:
            fingerprint = kmer_counter.genome_fingerprint(genome_fname)
            _fingerprint_cache[key] = fingerprint
    except OSError as exc:
        logger.debug(f"Could not fingerprint {genome_fname} ({exc})")
        return None
    return {
        "genome": os.path.abspath(genome_fname),
        "fingerprint": fingerprint,
        "k": int(k),
        "circular": int(bool(circular)),
    }


def _provenance_matches(stored, expected):
    """Whether a stored record describes the scan this run is asking for.

    Compared field by field rather than as whole dicts, because the record also
    carries the genome path and the path is not part of the identity. It is
    `os.path.abspath`, which normalises `..` and the working directory but not
    symlinks, so the same file reached through a symlinked parent, an automount,
    or `/tmp` against `/private/tmp` on macOS spells differently on the two runs.
    Comparing whole dicts made that a mismatch: the file was truncated and fully
    rescanned, and since the mismatch was neither absent nor a fingerprint
    difference the diagnostic blamed `circular` and printed the same value twice.

    The fingerprint is what identifies the genome, which is the reason for
    having one. The path stays in the record for the message.
    """
    return all(stored[field] == expected[field] for field in ("fingerprint", "k", "circular"))


def _stored_provenance(handle):
    """Read the provenance record off an open position index, or None."""
    attrs = handle.attrs
    if "provenance_fingerprint" not in attrs:
        return None
    try:
        return {
            "genome": str(attrs["provenance_genome"]),
            "fingerprint": str(attrs["provenance_fingerprint"]),
            "k": int(attrs["provenance_k"]),
            "circular": int(attrs["provenance_circular"]),
        }
    except (KeyError, TypeError, ValueError):
        return None


def write_position_provenance(fname_prefix, genome_fname, k, circular):
    """Record what a position file was scanned from, for the next run's check."""
    record = _position_provenance(genome_fname, k, circular)
    if record is None:
        return
    with h5py.File(position_file_path(fname_prefix, k), "a") as handle:
        handle.attrs["provenance_genome"] = record["genome"]
        handle.attrs["provenance_fingerprint"] = record["fingerprint"]
        handle.attrs["provenance_k"] = record["k"]
        handle.attrs["provenance_circular"] = record["circular"]


def clear_position_provenance(fname_prefix, k):
    """Drop the provenance record from a position file.

    Called by the fallback path, which reuses datasets without checking any
    record. A file it has written therefore holds data of unverified origin, and
    an absent record makes the Aho-Corasick path rescan rather than trust it.
    """
    h5_path = position_file_path(fname_prefix, k)
    if not os.path.exists(h5_path):
        return
    try:
        with h5py.File(h5_path, "a") as handle:
            for name in (
                "provenance_genome",
                "provenance_fingerprint",
                "provenance_k",
                "provenance_circular",
            ):
                if name in handle.attrs:
                    del handle.attrs[name]
    except OSError as exc:
        logger.debug(f"Could not clear provenance on {h5_path} ({exc})")


def _reusable_positions(primers, fname_prefix, genome_fname, k, circular):
    """Return `(positions, to_scan)` for `primers` against the HDF5 file.

    `positions` holds every primer that can be reused, mapped to the positions
    recorded for it; `to_scan` holds the rest. Together they always cover
    `primers`, which is the contract rather than an optimisation:
    `get_positions` documents a complete map and
    `primer_attributes.get_gini_from_txt_for_one_k` reads it with
    `.get((prefix, primer), [])`, so a primer left out of it reads as a primer
    with no binding sites, which the Gini gate turns into NaN and drops. On a
    re-run that would empty the whole candidate pool.

    The provenance check and the read happen under one file open, so the key set
    and the datasets cannot come from different states of the file and the
    record cannot be verified against one state and the data taken from another.

    A primer counts as reusable when the index holds an entry for it. On the
    Aho-Corasick path that test is sound: `write_to_h5py` is handed
    `all_positions.get(p, [])` for every pattern in the automaton, so a primer
    that occurs nowhere still gets an empty entry, and entry presence
    means "this primer was scanned" rather than "this primer was found". That
    invariant is what the reuse rests on; a scan path that stopped writing empty
    entries would turn it into a silent wrong answer.

    Returns `(positions, to_scan, replace)`. `replace` says the existing file
    must be discarded rather than merged into when the rescan is written, and it
    is set whenever the record is absent or does not match. Merging would leave
    the datasets this run does not ask about in place, and the file would then be
    stamped with the current genome while still holding another genome's
    positions for those primers, which a later run with a different candidate
    pool would reuse. The discard happens at write time rather than here so that
    an interrupted run leaves the previous file intact rather than an empty one.

    Any read problem returns everything for scanning, which is the behaviour
    before reuse existed.
    """
    h5_path = position_file_path(fname_prefix, k)
    if not os.path.exists(h5_path):
        return {}, list(primers), False

    expected = _position_provenance(genome_fname, k, circular)
    if expected is None:
        return {}, list(primers), True

    try:
        with position_index.open_index(h5_path) as index:
            stored = _stored_provenance(index)
            if stored is not None and _provenance_matches(stored, expected):
                positions = {
                    primer: found.tolist() for primer, found in index.get_many(primers).items()
                }
                to_scan = [primer for primer in primers if primer not in positions]
                return positions, to_scan, False
    except (OSError, KeyError) as exc:
        logger.debug(f"Could not read {h5_path} ({exc}); scanning all primers")
        return {}, list(primers), True

    if stored is None:
        logger.warning(
            f"{h5_path} exists but has no provenance record, so it cannot be "
            f"matched to {genome_fname}. Rescanning. Position files written "
            f"before this check was added are each rescanned once."
        )
    elif stored["fingerprint"] != expected["fingerprint"]:
        logger.warning(
            f"{h5_path} was scanned from a different genome "
            f"({stored.get('genome', 'unrecorded')}), not from {genome_fname}. "
            f"Rescanning. Reusing it would have reported the other genome's "
            f"binding sites."
        )
    elif stored["circular"] != expected["circular"]:
        logger.warning(
            f"{h5_path} was scanned with circular={bool(stored['circular'])} and "
            f"this run is circular={bool(circular)}. Rescanning. Reusing it "
            f"would have missed or invented sites spanning the origin."
        )
    else:
        logger.warning(
            f"{h5_path} records k={stored['k']} and this run asked for "
            f"k={expected['k']}. Rescanning."
        )

    return {}, list(primers), True


def get_positions(
    primer_list,
    fname_prefixes,
    fname_genomes,
    circular,
    overwrite=False,
    no_all_primer_files=False,
    k_values=None,
):
    """
    Launches a multiprocessing pool to check if all primers exists in their relevant h5py file and modifies the file
    if frequencies for that k-mer is missing.

    When the optional ``pyahocorasick`` package is installed, all k-values
    are scanned in a single genome pass (~5-7x faster for the default
    k-mer range).

    Args:
        primer_list: List of k-mers to be checked that exist in the h5py files or to modify them if not.
        fname_prefixes: The path prefixes for the h5py files, basically the path minus '_6mer_positions.h5' where k = 6.
        fname_genomes: A list of paths to the fasta files.
        overwrite: Boolean which when set to true means overwrite the k-mer entries in the h5py file if it already exists.
        k_values: Explicit k-mer lengths to scan. Defaults to the lengths
            actually present in `primer_list`.

    Returns:
        Dictionary mapping (prefix, primer) -> list of positions when using
        Aho-Corasick path. None when using the fallback parallel path.
        Callers can use this to avoid re-reading HDF5 files for immediate
        downstream calculations (e.g. Gini index).
    """
    # Scan the lengths we were actually asked about.
    #
    # This used to iterate `range(parameter.min_k, parameter.max_k + 1)`, so a
    # primer whose length fell outside the configured window was silently
    # dropped -- no error, no warning, and an empty position list that is
    # indistinguishable from "binds nowhere". That defeats the bring-your-own-
    # oligo case this function is the entry point for: an externally designed
    # 14-mer set against a params.json configured for 10-mers returned nothing
    # for every primer.
    #
    # An explicit `k_values` still wins, for callers that want to constrain it.
    if k_values is not None:
        k_range = sorted(set(k_values))
    else:
        k_range = sorted({len(primer) for primer in primer_list})

        outside = sorted(
            {length for length in k_range if length < parameter.min_k or length > parameter.max_k}
        )
        if outside:
            logger.info(
                "Scanning k-mer length(s) %s, which fall outside the configured "
                "min_k-max_k window (%s-%s). They are included because they were "
                "asked for; position files are written per length.",
                outside,
                parameter.min_k,
                parameter.max_k,
            )

    # Collect in-memory positions when possible (avoids HDF5 round-trip)
    position_cache = {}

    if AHOCORASICK_AVAILABLE and not overwrite:
        # Single-pass multi-k Aho-Corasick search.
        #
        # Primers whose positions are already in the HDF5 file are read back
        # rather than rescanned. `pipeline.step2` has always logged "Reusing
        # existing position files (incremental update only)" when the files were
        # present, and until this branch consulted them the log line was the only
        # thing that reused anything.
        #
        # Reuse is gated on a provenance record naming the genome and the scan
        # parameters, because this branch previously rescanned unconditionally
        # and so was immune to a stale file. `_reusable_positions` carries that
        # gate; see its docstring for why a mismatch truncates the file.
        for i, fg_prefix in enumerate(fname_prefixes):
            primer_lists_by_k = {}
            replace_k = set()
            reused_count = 0
            for k in k_range:
                k_primers = [p for p in primer_list if len(p) == k]
                if not k_primers:
                    continue
                rc_primers = [reverse_complement(p) for p in k_primers]
                # Sorted rather than `list(set(...))`: the scan result does not
                # depend on the order, but a stable order makes the pattern
                # count reproducible across runs.
                wanted = sorted(set(k_primers + rc_primers))

                reused, to_scan, replace = _reusable_positions(
                    wanted, fg_prefix, fname_genomes[i], k, circular
                )
                if replace:
                    replace_k.add(k)
                for primer, positions in reused.items():
                    position_cache[(fg_prefix, primer)] = positions
                if reused:
                    reused_count += len(reused)
                    logger.debug(
                        f"Reusing {len(reused):,} cached {k}-mer position set(s) "
                        f"for {fg_prefix}; scanning {len(to_scan):,}"
                    )
                if to_scan:
                    primer_lists_by_k[k] = sorted(to_scan)

            if reused_count:
                # Reported with counts because the previous announcement of reuse
                # in `pipeline.step2` was not conditional on anything being
                # reused, and nothing was.
                to_scan_count = sum(len(v) for v in primer_lists_by_k.values())
                logger.info(
                    f"Reusing {reused_count:,} position set(s) already in the HDF5 "
                    f"file(s) for {fg_prefix}; scanning {to_scan_count:,}"
                )

            if not primer_lists_by_k:
                continue

            # Ensure HDF5 files exist for each k before writing
            for k in primer_lists_by_k:
                h5_path = position_file_path(fg_prefix, k)
                if not os.path.exists(h5_path):
                    with h5py.File(h5_path, "a"):
                        pass

            all_positions = get_all_positions_multi_k(primer_lists_by_k, fname_genomes[i], circular)

            # Store in cache and write to HDF5 for persistence
            for primer, positions in all_positions.items():
                position_cache[(fg_prefix, primer)] = positions

            # Write to HDF5 grouped by k, then record what the file was scanned
            # from so the next run can decide whether it may be reused.
            for k, primers in primer_lists_by_k.items():
                k_dict = {p: all_positions.get(p, []) for p in primers}
                if k_dict:
                    write_to_h5py(
                        k_dict,
                        fg_prefix,
                        replace=k in replace_k,
                        record_starts=[0] + get_cached_record_boundaries(fname_genomes[i]),
                        genome_fname=fname_genomes[i],
                    )
                write_position_provenance(fg_prefix, fname_genomes[i], k, circular)

        return position_cache
    else:
        # Fallback: per-k parallel approach
        tasks = []
        for i, fg_prefix in enumerate(fname_prefixes):
            for k in k_range:
                tasks.append(
                    (
                        [primer for primer in primer_list if len(primer) == k],
                        fg_prefix,
                        fname_genomes[i],
                        k,
                        circular,
                        overwrite,
                    )
                )

        # Use context manager to ensure proper pool cleanup (prevents resource leaks)
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            pool.map(append_positions_to_h5py_file, tasks)

        return None


def append_positions_to_h5py_file(task):
    """
    Check if the primer exists in the h5py file and modify all
    the file if frequencies for that k-mer is missing.

    Args:
        task: A tuple consisting of the following arguments:
            primer_list: List of k-mers to be checked that exist in the h5py files or to modify them if not.
            fname_prefix: The path prefix for the h5py files, basically the path minus '_6mer_positions.h5' where k = 6.
            fname_genome: A list of paths to the fasta files.
            k: The length of the k-mers.
            overwrite: Boolean which when set to true means overwrite the k-mer entries in the h5py file if it already exists.

    """
    primer_list, fname_prefix, fname_genome, k, circular, overwrite = task

    if len(primer_list) == 0:
        return

    new_list = set(primer_list)
    new_list.update([reverse_complement(primer) for primer in new_list])
    one_k_list = [primer for primer in sorted(list(new_list)) if len(primer) == k]

    if not overwrite:
        filtered_list = check_which_primers_absent_in_h5py(one_k_list, fname_prefix)
    else:
        filtered_list = one_k_list

    if len(filtered_list) > 0:
        kmer_dict = get_all_positions_per_k(
            list(set(filtered_list)), fname_genome, circular, fname_prefix
        )
        write_to_h5py(
            kmer_dict,
            fname_prefix,
            record_starts=[0] + get_cached_record_boundaries(fname_genome),
            genome_fname=fname_genome,
        )
        # This path decides what to rescan with `check_which_primers_absent_in_h5py`,
        # which consults no provenance record, so the file it leaves behind holds
        # datasets of unverified origin. Dropping the record stops the
        # Aho-Corasick path trusting them; it rescans instead.
        clear_position_provenance(fname_prefix, k)
