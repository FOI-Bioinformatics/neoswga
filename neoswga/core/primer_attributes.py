import logging
import multiprocessing
import os

import h5py
import numpy as np

from neoswga.core import string_search as _string_search
from neoswga.core import utility as _utility
from neoswga.core.melting_temp import temp as _melting_temp
from neoswga.core.thermodynamics import reverse_complement

logger = logging.getLogger(__name__)


def get_melting_tm(primer):
    """Tm on the bundled model's fixed scale. NOT the Tm of your reaction.

    This returns `melting_temp.temp` at its defaults, 10 mM Na and 20 mM Mg.
    Those are the conditions the random forest's `melting_tm` feature is defined
    on, and they do not move with `na_conc`, `mg_conc`, the polymerase or any
    additive. On a run configured for 50 mM Na and 10 mM Mg the value is about
    7.9 C away from the reaction.

    It has no caller on the pipeline path, and is kept because it names the
    model's scale. For the Tm of a primer under the reaction actually being run,
    use `ReactionConditions.calculate_effective_tm` (what `filter.filter_extra`
    gates on) or `thermodynamics.calculate_tm_with_salt`.

    Args:
        primer: The sequence of the primer.

    Returns:
        Predicted melting temperature at 10 mM Na / 20 mM Mg.
    """
    return _melting_temp(primer)


def get_gini(primer, fname_prefixes):
    """
    Get the gini index of the gaps between all adjacent positions a primer may bind. This computed from positional gaps
    both the forward and reverse strand.

    Args:
        primer: The sequence of the primer.
        fname_prefixes: The list of path prefixes to the h5py files.

    Returns:
        gini: The computed gini index of the positional gaps.

    """
    k = len(primer)
    positions_diffs = []
    for i, fname_prefix in enumerate(fname_prefixes):
        h5_path = fname_prefix + "_" + str(k) + "mer_positions.h5"
        if os.path.exists(h5_path):
            with h5py.File(h5_path, "r") as db:
                if primer in db:
                    position_diffs_forward = _utility.get_positional_gap_lengths(db[primer])
                    positions_diffs.extend(position_diffs_forward)
                rc_primer = reverse_complement(primer)
                if rc_primer in db:
                    position_diffs_reverse = _utility.get_positional_gap_lengths(db[rc_primer])
                    positions_diffs.extend(position_diffs_reverse)
        else:
            logger.warning(f"Cannot find HDF5 file for prefix: {fname_prefix}")

    # `_utility.gini` does not exist -- the module exposes `gini_exact` -- so
    # every call to this function raised AttributeError. It has no caller inside
    # the package (`pipeline` uses the same-named helper in `filter`), which is
    # why that went unnoticed.
    #
    # Unlike get_gini_from_txt_for_one_k this treats the genome as linear: it
    # has no seq_length to wrap around with, so the origin-spanning gap is not
    # counted. Prefer that function when circularity matters.
    return _utility.gini_exact(positions_diffs)


def _load_positions_from_h5(primer_list, fname_prefix):
    """Try to load primer positions from existing HDF5 files.

    After step 2 creates position files, this avoids a redundant full
    genome scan by reading positions directly from HDF5.

    Args:
        primer_list: List of primer sequences (all same k).
        fname_prefix: Path prefix for HDF5 files.

    Returns:
        Dictionary mapping primer -> list of positions, or None if the
        HDF5 file does not exist or cannot be read.
    """
    if not primer_list:
        return None
    k = len(primer_list[0])
    h5_path = fname_prefix + "_" + str(k) + "mer_positions.h5"
    try:
        kmer_dict = {}
        with h5py.File(h5_path, "r") as f:
            for primer in primer_list:
                if primer in f:
                    kmer_dict[primer] = list(f[primer][:])
                else:
                    kmer_dict[primer] = []
        return kmer_dict
    except (FileNotFoundError, OSError):
        return None


# This probably belongs in a different file--its not a primer attribute.
# Below this many recorded binding sites, the Gini index of gap lengths is not
# a measurement of anything. One site gives no gap; two give a single gap,
# whose Gini is identically 0.0 -- the BEST score available. `filter.get_gini`
# keeps a primer when `gini.notna() & (gini < max_gini)`, so a 0.0 from an
# unmeasurable primer passed the evenness gate on an absence of evidence.
#
# Measured on the shipped pools: 86% of the 10,000-primer Prevotella-vs-chr21
# pool and 96.2% of the 500-row plasmid pool scored exactly 0.0, and every one
# of those had two or fewer foreground sites. The three whole-genome GC-tier
# pools have no rows at 0.0 at all, which is why this was invisible there.
#
# Three is the first count at which the number can vary. Sites are counted
# across both strands together; a per-strand rule would reject every primer
# that binds one strand only, which is a larger and separate change.
#
# This widens a narrower rule taken deliberately on 2026-08-01 in 0a40533,
# which held that only a TOTAL absence of positions is unmeasurable. That rule
# left a single-site primer scoring 0.0, the best value the gate can see.
#
# This is the DEFAULT, not a fixed rule. `min_gini_sites` in params.json and
# --min-gini-sites on `filter` override it. The value is threaded in as an
# argument rather than read from a module global here, because this function
# runs in a multiprocessing worker and macOS spawns rather than forks, so a
# global read inside the worker would see this default and not the configured
# value.
DEFAULT_MIN_GINI_SITES = 3


def get_gini_from_txt_for_one_k(
    primer_list,
    fname_prefix,
    fname_genome,
    seq_length,
    circular,
    position_cache=None,
    min_sites=None,
):
    """
    Measures the gini index of the gaps between all adjacent positions any primer in primer_list may bind to.

    Attempts to use the in-memory position cache first, then falls back
    to HDF5 files, then to a full genome scan.

    Args:
        primer_list: The list of primers to consider.
        fname_prefix: The path prefix to the h5py file.
        fname_genome: The path to the fasta file.
        seq_length: The length of the genome contained in fname_genome.
        circular: Whether the genome is circular.
        position_cache: Optional dict mapping (prefix, primer) -> positions
            from get_positions(). Avoids HDF5 round-trip when available.
        min_sites: Minimum combined forward-plus-reverse site count at which
            evenness is considered measurable. ``None`` means
            :data:`DEFAULT_MIN_GINI_SITES`. Passed as data rather than read
            from a global because this runs in a multiprocessing worker.

    Returns:
        primer_to_ginis: A dictionary of primers to a tuple of the gini indices (the first being computed from the
        forward strand, and the second being computed from the reverse strand).
    """
    rc_primer_list = [reverse_complement(primer) for primer in primer_list]
    all_primer_list = list(set(primer_list + rc_primer_list))

    kmer_dict = None

    # Try in-memory cache first (fastest, avoids all I/O)
    if position_cache is not None:
        kmer_dict = {}
        for primer in all_primer_list:
            kmer_dict[primer] = position_cache.get((fname_prefix, primer), [])

    # Try HDF5 files next
    if kmer_dict is None:
        kmer_dict = _load_positions_from_h5(all_primer_list, fname_prefix)

    if kmer_dict is None:
        # Fallback: scan genome (standalone usage without prior position creation)
        kmer_dict = _string_search.get_all_positions_per_k(
            kmer_list=all_primer_list,
            seq_fname=fname_genome,
            circular=circular,
            fname_prefix=fname_prefix,
        )

    # None means "the default". The caller resolves the configured value and
    # passes it in, because this function runs in a multiprocessing worker.
    min_sites = DEFAULT_MIN_GINI_SITES if min_sites is None else int(min_sites)

    ginis = []

    for primer in primer_list:
        # Evenness of spacing is only measurable from enough sites to have more
        # than one gap. See DEFAULT_MIN_GINI_SITES: below the threshold the
        # honest answer is "not measurable", and NaN is what makes
        # `filter.get_gini`'s existing `.notna()` guard fire. Returning 0.0
        # instead handed the gate the best score available on no evidence.
        forward_positions = kmer_dict.get(primer, [])
        reverse_positions = kmer_dict.get(reverse_complement(primer), [])

        if len(forward_positions) + len(reverse_positions) < min_sites:
            ginis.append((float("nan"), float("nan")))
            continue

        position_diffs_forward = _utility.get_positional_gap_lengths(
            forward_positions, circular, seq_length=seq_length
        )
        gini_forward = _utility.gini_exact(position_diffs_forward)
        position_diffs_reverse = _utility.get_positional_gap_lengths(
            reverse_positions, circular, seq_length=seq_length
        )
        gini_reverse = _utility.gini_exact(position_diffs_reverse)
        ginis.append((gini_forward, gini_reverse))

    primer_to_ginis = dict(zip(primer_list, ginis))
    return primer_to_ginis


def get_gini_from_txt_for_one_k_helper(args):
    """Multiprocessing wrapper for :func:`get_gini_from_txt_for_one_k`.

    Unpacks a positional-argument tuple and delegates to the underlying
    function. Accepts 5, 6 or 7 elements; the sixth is the position cache and
    the seventh is the minimum site count. The count is passed rather than read
    from the parameter module because this runs in a spawned worker, which
    would see the module's import-time default instead of the configured value.

    Args:
        args: Tuple of (primer_list, fname_prefix, fname_genome,
            seq_length, circular[, position_cache[, min_sites]]).

    Returns:
        dict mapping each primer to a (gini_forward, gini_reverse) tuple.
    """
    primer_list, fname_prefix, fname_genome, seq_length, circular = args[:5]
    position_cache = args[5] if len(args) > 5 else None
    min_sites = args[6] if len(args) > 6 else None
    return get_gini_from_txt_for_one_k(
        primer_list,
        fname_prefix,
        fname_genome,
        seq_length,
        circular,
        position_cache,
        min_sites=min_sites,
    )


def get_gini_from_txt(
    primer_list, fname_prefixes, fname_genomes, seq_lengths, circular, position_cache=None
):
    """
    This runs get_gini_from_txt_for_one_k in a multiprocessed fashion where the task is divided based on the length
    of the primers.

    Args:
        primer_list: The list of primers to consider.
        fname_prefixes: The list of path prefixes to the h5py file.
        fname_genomes: The list of paths to the fasta files.
        seq_lengths: The length of all the genomes in fname_genomes (in the same order!).
        position_cache: Optional dict mapping (prefix, primer) -> positions
            from get_positions(). When provided, Gini calculation uses
            in-memory data instead of reading back from HDF5.

    Returns:
        The average gini_index across all gini indices computed for each primer.
    """
    # Use parameter k-mer range instead of hardcoded [6-12]
    import neoswga.core.parameter as parameter

    k_range = range(parameter.min_k, parameter.max_k + 1)

    # Resolved here, in the parent, and passed as data. See the helper: a
    # spawned worker reading the module global would see the import-time
    # default rather than the configured value, and the parallel path would
    # then disagree with the serial one on the same primers.
    configured = getattr(parameter, "min_gini_sites", None)
    min_sites = (
        int(configured)
        if isinstance(configured, int) and not isinstance(configured, bool) and configured > 0
        else DEFAULT_MIN_GINI_SITES
    )

    # When we have an in-memory cache, compute directly (no multiprocessing
    # needed since the data is already loaded -- avoids pickle overhead)
    if position_cache is not None:
        results = []
        for i, fg_prefix in enumerate(fname_prefixes):
            for k in k_range:
                primer_list_a = [primer for primer in primer_list if len(primer) == k]
                if len(primer_list_a) > 0:
                    result = get_gini_from_txt_for_one_k(
                        primer_list_a,
                        fg_prefix,
                        fname_genomes[i],
                        seq_lengths[i],
                        circular,
                        position_cache,
                        min_sites=min_sites,
                    )
                    results.append(result)
    else:
        tasks = []
        for i, fg_prefix in enumerate(fname_prefixes):
            for k in k_range:
                primer_list_a = [primer for primer in primer_list if len(primer) == k]
                if len(primer_list_a) > 0:
                    # The explicit None in sixth position: the helper unpacks
                    # positionally, so the site count has to sit seventh.
                    tasks.append(
                        [
                            primer_list_a,
                            fg_prefix,
                            fname_genomes[i],
                            seq_lengths[i],
                            circular,
                            None,
                            min_sites,
                        ]
                    )

        # Use context manager to ensure proper pool cleanup (prevents resource leaks)
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.map(get_gini_from_txt_for_one_k_helper, tasks)

    primer_to_all_ginis = {}

    for primer_to_ginis in results:
        for primer, (gini_forward, gini_reverse) in primer_to_ginis.items():
            if primer not in primer_to_all_ginis:
                primer_to_all_ginis[primer] = []
            primer_to_all_ginis[primer].append(np.mean([gini_forward, gini_reverse]))
    gini_mean = [np.mean(primer_to_all_ginis[primer]) for primer in primer_list]
    return gini_mean


def get_rate_from_h5py(primer, fname_prefixes):
    """
    Gets the frequency of a primer in both the forward and reverse strand.

    Args:
        primer: The sequence of the primer.
        fname_prefixes: The list of path prefixes to all the relevant h5py files.

    Returns:
        count: The frequency of a primer in both the forward and reverse strand.
    """
    k = len(primer)
    count = 0
    for i, fname_prefix in enumerate(fname_prefixes):
        h5_path = fname_prefix + "_" + str(k) + "mer_positions.h5"
        if os.path.exists(h5_path):
            with h5py.File(h5_path, "r") as db:
                if primer in db:
                    count += len(db[primer])
                rc_primer = reverse_complement(primer)
                if rc_primer in db:
                    count += len(db[rc_primer])
        else:
            logger.warning(f"Cannot find HDF5 file for prefix: {fname_prefix}")
    return count
