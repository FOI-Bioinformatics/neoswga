import neoswga.core.utility
import neoswga.core.parameter as parameter
import multiprocessing
import h5py
import os
import logging
from typing import Dict, List, Optional

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

    # Read entire genome into memory
    seq_generator = neoswga.core.utility.read_fasta_file(seq_fname)
    sequence = ''.join(seq_generator).upper()

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
        'num_genomes': len(_genome_cache),
        'total_bp': total_bp,
        'memory_mb': total_bp // (1024 * 1024)
    }


def get_all_positions_multi_k(primer_lists_by_k, seq_fname, circular):
    """Scan genome once for primers of all k-values using Aho-Corasick.

    Requires the optional ``pyahocorasick`` package.  When available this
    replaces per-k scanning and provides ~5-7x speedup for default k
    ranges (6-12).

    Args:
        primer_lists_by_k: Dict mapping k -> list of primers.
        seq_fname: Path to genome FASTA.
        circular: Whether the genome is circular.

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

    if circular:
        max_k = max(primer_lists_by_k.keys())
        search_seq = sequence + sequence[:max_k - 1]
    else:
        search_seq = sequence

    for end_pos, primer in A.iter(search_seq):
        start_pos = end_pos - len(primer) + 1
        if start_pos < seq_len:
            all_primers[primer].append(start_pos)

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
        print("Starting the search for " + fname_prefix + ' ' + str(k) + 'mers...')

    # Optimized sliding window search using cached sequence
    seq_len = len(sequence)
    search_len = seq_len if not circular else seq_len + k - 1

    for i in range(search_len - k + 1):
        if i < seq_len - k + 1:
            # Normal position within sequence
            current_kmer = sequence[i:i+k]
        else:
            # Circular wrap-around: combine end with beginning
            wrap_pos = i - (seq_len - k + 1)
            current_kmer = sequence[seq_len - k + 1 + wrap_pos:] + sequence[:wrap_pos + 1]

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
    h5_path = fname_prefix + '_' + str(k) + 'mer_positions.h5'
    with h5py.File(h5_path, 'r+') as f:
        for kmer, positions in kmer_dict.items():
            if kmer not in f:
                f.create_dataset(kmer, data=positions)
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
    h5_path = fname_prefix + '_' + str(k) + 'mer_positions.h5'

    if not os.path.exists(h5_path):
        # Create empty HDF5 file using context manager
        with h5py.File(h5_path, 'a'):
            pass
        return primer_list

    # Load all k-mers present in genome
    all_present_kmers_in_genome = set()
    txt_path = fname_prefix + '_' + str(k) + 'mer_all.txt'
    with open(txt_path, 'r') as txt_f:
        for line in txt_f:
            curr_kmer = line.split(" ")[0]
            all_present_kmers_in_genome.add(curr_kmer)

    # Get existing keys from HDF5 file using context manager
    with h5py.File(h5_path, 'r') as f:
        keys = set(f.keys())

    filtered_primer_list = []
    for primer in primer_list:
        if primer not in keys:
            if primer in all_present_kmers_in_genome:
                filtered_primer_list.append(primer)
    return filtered_primer_list

def get_positions(primer_list, fname_prefixes, fname_genomes, circular, overwrite=False, no_all_primer_files=False):
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

    Returns:
        Dictionary mapping (prefix, primer) -> list of positions when using
        Aho-Corasick path. None when using the fallback parallel path.
        Callers can use this to avoid re-reading HDF5 files for immediate
        downstream calculations (e.g. Gini index).
    """
    # Use parameter k-mer range instead of hardcoded [6-12]
    k_range = range(parameter.min_k, parameter.max_k + 1)

    # Collect in-memory positions when possible (avoids HDF5 round-trip)
    position_cache = {}

    if AHOCORASICK_AVAILABLE and not overwrite:
        # Single-pass multi-k Aho-Corasick search
        for i, fg_prefix in enumerate(fname_prefixes):
            primer_lists_by_k = {}
            for k in k_range:
                k_primers = [p for p in primer_list if len(p) == k]
                if not k_primers:
                    continue
                rc_primers = [reverse_complement(p) for p in k_primers]
                primer_lists_by_k[k] = list(set(k_primers + rc_primers))

            if not primer_lists_by_k:
                continue

            # Ensure HDF5 files exist for each k before writing
            for k in primer_lists_by_k:
                h5_path = fg_prefix + '_' + str(k) + 'mer_positions.h5'
                if not os.path.exists(h5_path):
                    with h5py.File(h5_path, 'a'):
                        pass

            all_positions = get_all_positions_multi_k(
                primer_lists_by_k, fname_genomes[i], circular
            )

            # Store in cache and write to HDF5 for persistence
            for primer, positions in all_positions.items():
                position_cache[(fg_prefix, primer)] = positions

            # Write to HDF5 grouped by k
            for k, primers in primer_lists_by_k.items():
                k_dict = {p: all_positions.get(p, []) for p in primers}
                if k_dict:
                    write_to_h5py(k_dict, fg_prefix)

        return position_cache
    else:
        # Fallback: per-k parallel approach
        tasks = []
        for i, fg_prefix in enumerate(fname_prefixes):
            for k in k_range:
                tasks.append(([primer for primer in primer_list if len(primer) == k], fg_prefix, fname_genomes[i], k, circular, overwrite))

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
        kmer_dict = get_all_positions_per_k(list(set(filtered_list)), fname_genome, circular, fname_prefix)
        write_to_h5py(kmer_dict, fname_prefix)