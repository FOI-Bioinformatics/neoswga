import logging
import multiprocessing
import os
from typing import Dict, List, Optional

import h5py

from neoswga.core import parameter
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
_genome_cache: Dict[str, str] = {}


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
    for record in loader.load_genome_streaming(seq_fname):
        sequence += record

    _genome_cache[seq_fname] = sequence
    logger.info(f"Cached genome sequence: {len(sequence):,} bp")

    return sequence


def preload_genomes(seq_fnames: List[str]) -> None:
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
    global _genome_cache
    cache_size = len(_genome_cache)
    _genome_cache = {}
    logger.info(f"Cleared genome cache ({cache_size} entries)")


def get_genome_cache_stats() -> Dict[str, int]:
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
            if abs_start < seq_len:
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

    for i in range(search_len - k + 1):
        if i < seq_len - k + 1:
            # Normal position within sequence
            current_kmer = sequence[i : i + k]
        else:
            # Circular wrap-around: combine end with beginning
            wrap_pos = i - (seq_len - k + 1)
            current_kmer = sequence[seq_len - k + 1 + wrap_pos :] + sequence[: wrap_pos + 1]

        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].append(i)

    return kmer_dict


def write_to_h5py(kmer_dict, fname_prefix):
    """
    Writes the kmer counts to an h5py file, which allows for efficient access in terms of looking up the
    frequency of a particular k-mer. If the kmer already exists in the dataset, the entry in the h5py file
    is overwritten with the new data.

    Args:
        kmer_dict: A dictionary of all the k-mers and their respective frequencies.
        fname_prefix: The file path prefix for the output h5py file. If k=6, for example, '_6mer_positions.h5'
        will be appended to the file name.
    """
    if not kmer_dict:
        return  # Nothing to write
    k = len(next(iter(kmer_dict.keys())))
    h5_path = fname_prefix + "_" + str(k) + "mer_positions.h5"
    with h5py.File(h5_path, "r+") as f:
        for kmer, positions in kmer_dict.items():
            if kmer not in f:
                f.create_dataset(kmer, data=positions)
            elif len(f[kmer]) == len(positions):
                # In-place overwrite avoids HDF5 file fragmentation
                f[kmer][...] = positions
            else:
                del f[kmer]
                f.create_dataset(kmer, data=positions)


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

    # Load all k-mers present in genome
    all_present_kmers_in_genome = set()
    txt_path = fname_prefix + "_" + str(k) + "mer_all.txt"
    with open(txt_path, "r") as txt_f:
        for line in txt_f:
            curr_kmer = line.split(" ")[0]
            all_present_kmers_in_genome.add(curr_kmer)

    # Get existing keys from HDF5 file using context manager
    with h5py.File(h5_path, "r") as f:
        keys = set(f.keys())

    filtered_primer_list = []
    for primer in primer_list:
        if primer not in keys:
            if primer in all_present_kmers_in_genome:
                filtered_primer_list.append(primer)
    return filtered_primer_list


# Genome fingerprints already computed in this process, keyed by path, size and
# mtime. See `_position_provenance`.
_fingerprint_cache: Dict[tuple, str] = {}


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


def _stored_provenance(handle):
    """Read the provenance record off an open position file, or None."""
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

    A primer counts as reusable when the file holds a dataset for it. On the
    Aho-Corasick path that test is sound: `write_to_h5py` is handed
    `all_positions.get(p, [])` for every pattern in the automaton, so a primer
    that occurs nowhere still gets a key with an empty dataset, and key presence
    means "this primer was scanned" rather than "this primer was found". That
    invariant is what the reuse rests on; a scan path that stopped writing empty
    entries would turn it into a silent wrong answer.

    A file whose record is absent or does not match is truncated rather than
    merely bypassed. Bypassing leaves the datasets that this run does not ask
    about in place, and the file is then stamped with the current genome while
    still holding another genome's positions for those primers, which a later
    run with a different candidate pool would reuse.

    Any read problem returns everything for scanning, which is the behaviour
    before reuse existed.
    """
    h5_path = position_file_path(fname_prefix, k)
    if not os.path.exists(h5_path):
        return {}, list(primers)

    expected = _position_provenance(genome_fname, k, circular)
    if expected is None:
        return {}, list(primers)

    try:
        with h5py.File(h5_path, "r") as handle:
            stored = _stored_provenance(handle)
            if stored == expected:
                positions = {
                    primer: handle[primer][:].tolist() for primer in primers if primer in handle
                }
                to_scan = [primer for primer in primers if primer not in positions]
                return positions, to_scan
    except (OSError, KeyError) as exc:
        logger.debug(f"Could not read {h5_path} ({exc}); scanning all primers")
        return {}, list(primers)

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
    else:
        logger.warning(
            f"{h5_path} was scanned with circular={bool(stored['circular'])} and "
            f"this run is circular={bool(circular)}. Rescanning. Reusing it "
            f"would have missed or invented sites spanning the origin."
        )

    try:
        with h5py.File(h5_path, "w"):
            pass
    except OSError as exc:
        logger.debug(f"Could not truncate {h5_path} ({exc})")
    return {}, list(primers)


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

                reused, to_scan = _reusable_positions(
                    wanted, fg_prefix, fname_genomes[i], k, circular
                )
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
                    write_to_h5py(k_dict, fg_prefix)
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
        write_to_h5py(kmer_dict, fname_prefix)
        # This path decides what to rescan with `check_which_primers_absent_in_h5py`,
        # which consults no provenance record, so the file it leaves behind holds
        # datasets of unverified origin. Dropping the record stops the
        # Aho-Corasick path trusting them; it rescans instead.
        clear_position_provenance(fname_prefix, k)
