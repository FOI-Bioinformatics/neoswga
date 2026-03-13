from neoswga.core import utility as _utility
import os
import h5py
from neoswga.core import string_search as _string_search
from neoswga.core.melting_temp import temp as _melting_temp
import multiprocessing
import numpy as np
import logging

from neoswga.core.thermodynamics import reverse_complement

logger = logging.getLogger(__name__)

def get_melting_tm(primer):
    """
    Gets the predicted melting temperature of a primer using the melt package (https://pypi.org/project/melt/).

    Args:
        primer: The sequence of the primer.

    Returns:
        The predicted melting temperature.
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
        h5_path = fname_prefix + '_' + str(k) + 'mer_positions.h5'
        if os.path.exists(h5_path):
            with h5py.File(h5_path, 'r') as db:
                if primer in db:
                    position_diffs_forward = _utility.get_positional_gap_lengths(db[primer])
                    positions_diffs.extend(position_diffs_forward)
                rc_primer = reverse_complement(primer)
                if rc_primer in db:
                    position_diffs_reverse = _utility.get_positional_gap_lengths(db[rc_primer])
                    positions_diffs.extend(position_diffs_reverse)
        else:
            logger.warning(f"Cannot find HDF5 file for prefix: {fname_prefix}")
    gini = _utility.gini(positions_diffs)
    return gini

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
    h5_path = fname_prefix + '_' + str(k) + 'mer_positions.h5'
    try:
        kmer_dict = {}
        with h5py.File(h5_path, 'r') as f:
            for primer in primer_list:
                if primer in f:
                    kmer_dict[primer] = list(f[primer][:])
                else:
                    kmer_dict[primer] = []
        return kmer_dict
    except (FileNotFoundError, OSError):
        return None


# This probably belongs in a different file--its not a primer attribute.
def get_gini_from_txt_for_one_k(primer_list, fname_prefix, fname_genome, seq_length, circular, position_cache=None):
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
            kmer_list=all_primer_list, seq_fname=fname_genome,
            circular=circular, fname_prefix=fname_prefix
        )

    ginis = []

    for primer in primer_list:
        position_diffs_forward = _utility.get_positional_gap_lengths(kmer_dict.get(primer, []), circular, seq_length=seq_length)
        gini_forward = _utility.gini_exact(position_diffs_forward)
        position_diffs_reverse = _utility.get_positional_gap_lengths(kmer_dict.get(reverse_complement(primer), []), circular, seq_length=seq_length)
        gini_reverse = _utility.gini_exact(position_diffs_reverse)
        ginis.append((gini_forward, gini_reverse))

    primer_to_ginis = dict(zip(primer_list, ginis))
    return primer_to_ginis

def get_gini_from_txt_for_one_k_helper(args):
    """Multiprocessing wrapper for :func:`get_gini_from_txt_for_one_k`.

    Unpacks a positional-argument tuple and delegates to the underlying
    function. Accepts either 5 or 6 elements; when 6 are provided the
    last element is used as the position cache.

    Args:
        args: Tuple of (primer_list, fname_prefix, fname_genome,
            seq_length, circular[, position_cache]).

    Returns:
        dict mapping each primer to a (gini_forward, gini_reverse) tuple.
    """
    if len(args) == 6:
        primer_list, fname_prefix, fname_genome, seq_length, circular, position_cache = args
    else:
        primer_list, fname_prefix, fname_genome, seq_length, circular = args
        position_cache = None
    return get_gini_from_txt_for_one_k(primer_list, fname_prefix, fname_genome, seq_length, circular, position_cache)

def get_gini_from_txt(primer_list, fname_prefixes, fname_genomes, seq_lengths, circular, position_cache=None):
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

    # When we have an in-memory cache, compute directly (no multiprocessing
    # needed since the data is already loaded -- avoids pickle overhead)
    if position_cache is not None:
        results = []
        for i, fg_prefix in enumerate(fname_prefixes):
            for k in k_range:
                primer_list_a = [primer for primer in primer_list if len(primer) == k]
                if len(primer_list_a) > 0:
                    result = get_gini_from_txt_for_one_k(
                        primer_list_a, fg_prefix, fname_genomes[i],
                        seq_lengths[i], circular, position_cache
                    )
                    results.append(result)
    else:
        tasks = []
        for i, fg_prefix in enumerate(fname_prefixes):
            for k in k_range:
                primer_list_a = [primer for primer in primer_list if len(primer) == k]
                if len(primer_list_a) > 0:
                    tasks.append([primer_list_a, fg_prefix, fname_genomes[i], seq_lengths[i], circular])

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
        h5_path = fname_prefix + '_' + str(k) + 'mer_positions.h5'
        if os.path.exists(h5_path):
            with h5py.File(h5_path, 'r') as db:
                if primer in db:
                    count += len(db[primer])
                rc_primer = reverse_complement(primer)
                if rc_primer in db:
                    count += len(db[rc_primer])
        else:
            logger.warning(f"Cannot find HDF5 file for prefix: {fname_prefix}")
    return count