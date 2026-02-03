"""
K-mer counting utilities using Jellyfish.

DEPRECATED: This module is deprecated. Use neoswga.core.kmer_counter instead,
which provides:
- MultiGenomeKmerCounter class for multi-genome support
- Better error handling and type hints
- Caching and performance optimizations

This module is kept for backwards compatibility but will be removed in a future version.
"""

import os
import subprocess
import logging
import warnings
import neoswga.core.parameter as parameter
import melting

logger = logging.getLogger(__name__)

# Issue deprecation warning at module import
warnings.warn(
    "neoswga.core.kmer is deprecated. Use neoswga.core.kmer_counter instead. "
    "See MultiGenomeKmerCounter for improved functionality.",
    DeprecationWarning,
    stacklevel=2
)


def run_jellyfish(genome_fname=None, output_prefix=None, min_k=6, max_k=12):
    """
    Run jellyfish to count k-mers and generate output files.

    DEPRECATED: Use neoswga.core.kmer_counter.MultiGenomeKmerCounter instead.

    Args:
        genome_fname: The fasta file used to count kmers.
        output_prefix: The output path prefix for the output files.
            Resulting output files will be suffixed by _kmer_all.txt
            for k from min_k to max_k inclusive.
        min_k: Minimum k-mer length (default: 6)
        max_k: Maximum k-mer length (default: 12, up to 18 with additives)

    Raises:
        FileNotFoundError: If genome_fname does not exist.
        subprocess.CalledProcessError: If jellyfish command fails.
    """
    warnings.warn(
        "run_jellyfish() is deprecated. Use MultiGenomeKmerCounter.count_kmers_jellyfish() instead.",
        DeprecationWarning,
        stacklevel=2
    )
    if genome_fname is None:
        raise ValueError("genome_fname is required")
    if output_prefix is None:
        raise ValueError("output_prefix is required")
    if not os.path.exists(genome_fname):
        raise FileNotFoundError(f"Genome file not found: {genome_fname}")

    for k in range(min_k, max_k + 1):
        jf_file = f"{output_prefix}_{k}mer_all.jf"
        txt_file = f"{output_prefix}_{k}mer_all.txt"

        if not os.path.exists(txt_file):
            # SECURE: Use list arguments to prevent command injection
            count_cmd = [
                'jellyfish', 'count',
                '-m', str(k),
                '-s', '1000000',
                '-t', str(parameter.cpus),
                genome_fname,
                '-o', jf_file
            ]
            logger.debug(f"Running: {' '.join(count_cmd)}")
            subprocess.run(count_cmd, check=True, capture_output=True)

            # SECURE: Use list arguments and file redirection via Python
            dump_cmd = ['jellyfish', 'dump', '-c', jf_file]
            logger.debug(f"Running: {' '.join(dump_cmd)}")
            with open(txt_file, 'w') as f_out:
                subprocess.run(dump_cmd, check=True, stdout=f_out)

        # Clean up intermediate .jf file
        if os.path.exists(jf_file):
            os.remove(jf_file)

def get_kmer_to_count_dict(f_in_name):
    """
    Computes the counts of all kmers in the 5' to 3' direction.

    DEPRECATED: Use neoswga.core.kmer_counter.MultiGenomeKmerCounter instead.

    Args:
        f_in_name: The path to the k-mer count file.

    Returns:
        primer_to_count_dict: Dictionary mapping kmer to count in the 5' to 3' direction.
    """
    warnings.warn(
        "get_kmer_to_count_dict() is deprecated. Use MultiGenomeKmerCounter instead.",
        DeprecationWarning,
        stacklevel=2
    )
    primer_to_count_dict = {}

    with open(f_in_name, 'r') as f_in:
        for line in f_in:
            primer = line.split(" ")[0]
            count = int(line.split(" ")[1].rstrip())
            primer_to_count_dict[primer] = count

    return primer_to_count_dict

def get_primer_list_from_kmers(prefixes, kmer_lengths=None):
    """
    Gets all the kmers from the jellyfish output files.

    DEPRECATED: Use neoswga.core.kmer_counter.MultiGenomeKmerCounter instead.

    Args:
        prefixes: The prefix path that all the jellyfish output files share.
        kmer_lengths: Optional list of k-mer lengths to include.

    Returns:
        primer_list: List of all the kmers that occur at least and satisfy the temperature conditions.
    """
    warnings.warn(
        "get_primer_list_from_kmers() is deprecated. Use MultiGenomeKmerCounter instead.",
        DeprecationWarning,
        stacklevel=2
    )
    primer_list = []

    if not kmer_lengths:
        kmer_lengths = range(6,13,1)

    for prefix in prefixes:
        for k in kmer_lengths:
            with open(prefix+'_'+str(k)+'mer_all.txt', 'r') as f_in:
                for line in f_in:
                    curr_kmer = line.split(" ")[0]
                    tm = melting.temp(curr_kmer)
                    if tm < parameter.max_tm and tm > parameter.min_tm:
                        primer_list.append(curr_kmer)
    return primer_list
