"""
Multi-genome k-mer counter for efficient k-mer frequency analysis.

Uses Jellyfish for fast k-mer counting across multiple genomes.
Jellyfish is a required dependency - the module will raise an error if not available.
"""

import logging
import os
import subprocess
import tempfile
import shutil
import sys
from typing import Dict, List, Optional
from collections import defaultdict
from pathlib import Path

from neoswga.core.thermodynamics import reverse_complement

logger = logging.getLogger(__name__)


def _print_k_progress(k: int, min_k: int, max_k: int):
    """Print inline progress for k-mer counting."""
    total = max_k - min_k + 1
    current = k - min_k + 1
    sys.stdout.write(f"\r  k={k}bp ({current}/{total})")
    sys.stdout.flush()
    if k == max_k:
        sys.stdout.write(" done\n")


def check_jellyfish_available() -> bool:
    """Check if jellyfish is available in PATH."""
    return shutil.which('jellyfish') is not None


def get_jellyfish_version() -> Optional[str]:
    """Parse the Jellyfish version from ``jellyfish --version``.

    Returns:
        Version string (e.g. '2.3.1'), or None if unavailable.
    """
    if not check_jellyfish_available():
        return None
    try:
        result = subprocess.run(
            ['jellyfish', '--version'],
            capture_output=True, text=True, timeout=10
        )
        # Output is typically "jellyfish 2.3.1" or just "2.3.1"
        text = result.stdout.strip() or result.stderr.strip()
        for token in text.split():
            if token[0].isdigit():
                return token
    except (subprocess.TimeoutExpired, OSError):
        pass
    return None


def require_jellyfish():
    """Raise error if Jellyfish is not available or is version 1.x."""
    if not check_jellyfish_available():
        raise RuntimeError(
            "Jellyfish is required but not found in PATH. "
            "Please install Jellyfish: https://github.com/gmarcais/Jellyfish"
        )
    version = get_jellyfish_version()
    if version:
        logger.info(f"Jellyfish version: {version}")
        major = version.split('.')[0]
        if major == '1':
            raise RuntimeError(
                f"Jellyfish version {version} detected. NeoSWGA requires "
                f"Jellyfish 2.x (1.x has an incompatible CLI). "
                f"Please upgrade: https://github.com/gmarcais/Jellyfish"
            )


class MultiGenomeKmerCounter:
    """
    Efficient k-mer counter for multiple genomes using Jellyfish.

    Jellyfish is required for this class - it provides 100x faster k-mer
    counting compared to pure Python implementations.
    """

    def __init__(self, cpus: int = 4, output_dir: Optional[str] = None):
        """
        Initialize counter.

        Args:
            cpus: Number of CPUs for Jellyfish
            output_dir: Directory for k-mer count files (temp if None)

        Raises:
            RuntimeError: If Jellyfish is not available
        """
        require_jellyfish()

        self.cpus = cpus
        self.output_dir = output_dir
        self._temp_dir = None

        # Cache
        self.genome_fastas: Dict[str, str] = {}  # name -> fasta_path
        self.genome_lengths: Dict[str, int] = {}
        self.kmer_counts: Dict[str, Dict[str, int]] = {}  # genome -> {kmer: count}
        self.kmer_files: Dict[str, Dict[int, str]] = {}  # genome -> {k: file_path}

    def _get_work_dir(self) -> str:
        """Get working directory for intermediate files."""
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
            return self.output_dir
        if self._temp_dir is None:
            self._temp_dir = tempfile.mkdtemp(prefix='neoswga_kmer_')
        return self._temp_dir

    def add_genome(self, name: str, fasta_path: str):
        """
        Add a genome for k-mer counting.

        Args:
            name: Genome identifier
            fasta_path: Path to FASTA file
        """
        if not os.path.exists(fasta_path):
            raise FileNotFoundError(f"Genome file not found: {fasta_path}")

        self.genome_fastas[name] = fasta_path

        # Calculate genome length
        length = 0
        with open(fasta_path) as f:
            for line in f:
                if not line.startswith('>'):
                    length += len(line.strip())
        self.genome_lengths[name] = length

        logger.debug(f"Added genome {name}: {length:,} bp from {fasta_path}")

    def count_kmers_jellyfish(self, genome_name: str, k: int) -> Dict[str, int]:
        """
        Count k-mers using Jellyfish.

        Args:
            genome_name: Name of genome to count
            k: K-mer length

        Returns:
            Dict mapping k-mer to count
        """
        if genome_name not in self.genome_fastas:
            raise ValueError(f"Unknown genome: {genome_name}")

        # Check cache
        cache_key = f"{genome_name}_{k}"
        if cache_key in self.kmer_counts:
            return self.kmer_counts[cache_key]

        fasta_path = self.genome_fastas[genome_name]
        work_dir = self._get_work_dir()

        # Run jellyfish count
        jf_file = os.path.join(work_dir, f"{genome_name}_{k}mer.jf")
        txt_file = os.path.join(work_dir, f"{genome_name}_{k}mer_all.txt")

        if not os.path.exists(txt_file):
            # Run jellyfish count
            count_cmd = [
                'jellyfish', 'count',
                '-m', str(k),
                '-s', '1000000',
                '-t', str(self.cpus),
                '-C',  # Canonical k-mers (both strands)
                fasta_path,
                '-o', jf_file
            ]
            subprocess.run(count_cmd, check=True, capture_output=True)

            # Dump to text
            dump_cmd = ['jellyfish', 'dump', '-c', jf_file]
            with open(txt_file, 'w') as f:
                subprocess.run(dump_cmd, check=True, stdout=f)

            # Clean up .jf file
            if os.path.exists(jf_file):
                os.remove(jf_file)

        # Parse results
        counts = {}
        with open(txt_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    kmer = parts[0]
                    count = int(parts[1])
                    counts[kmer] = count

        # Cache and return
        self.kmer_counts[cache_key] = counts
        return counts

    def count_kmers(self, genome_name: str, k: int) -> Dict[str, int]:
        """
        Count k-mers using Jellyfish.

        Args:
            genome_name: Name of genome to count
            k: K-mer length

        Returns:
            Dict mapping k-mer to count
        """
        return self.count_kmers_jellyfish(genome_name, k)

    def get_kmer_count(self, kmer: str, genome_name: str) -> int:
        """
        Get count of a specific k-mer in a genome.

        Args:
            kmer: K-mer sequence
            genome_name: Genome to query

        Returns:
            Count of k-mer (including reverse complement)
        """
        k = len(kmer)
        counts = self.count_kmers(genome_name, k)

        kmer = kmer.upper()
        kmer_rc = reverse_complement(kmer)
        canonical = min(kmer, kmer_rc)

        return counts.get(canonical, 0)

    def get_kmer_frequency(self, kmer: str, genome_name: str) -> float:
        """
        Get frequency of k-mer in a genome.

        Args:
            kmer: K-mer sequence
            genome_name: Genome to query

        Returns:
            Frequency (count / genome_length)
        """
        if genome_name not in self.genome_lengths:
            return 0.0

        count = self.get_kmer_count(kmer, genome_name)
        return count / self.genome_lengths[genome_name]

    def cleanup(self):
        """Remove temporary files."""
        if self._temp_dir and os.path.exists(self._temp_dir):
            shutil.rmtree(self._temp_dir)
            self._temp_dir = None


def count_kmers_in_sequence(sequence: str, k: int) -> Dict[str, int]:
    """
    Count all k-mers in a sequence (pure Python).

    Args:
        sequence: DNA sequence (will be uppercased)
        k: K-mer length

    Returns:
        Dict mapping k-mer to count
    """
    sequence = sequence.upper()
    counts = defaultdict(int)

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        # Only count canonical DNA k-mers
        if all(base in 'ACGT' for base in kmer):
            counts[kmer] += 1

    return dict(counts)


# =============================================================================
# Standalone Functions (replacements for deprecated kmer.py)
# =============================================================================

def run_jellyfish(genome_fname: str, output_prefix: str,
                  min_k: int = 6, max_k: int = 12, cpus: int = 4) -> None:
    """
    Run jellyfish to count k-mers and generate output files.

    This is a standalone function that wraps MultiGenomeKmerCounter
    for simple use cases.

    Args:
        genome_fname: Path to FASTA file
        output_prefix: Output path prefix (files will be suffixed by _Xmer_all.txt)
        min_k: Minimum k-mer length (default: 6)
        max_k: Maximum k-mer length (default: 12)
        cpus: Number of CPUs (default: 4)

    Raises:
        FileNotFoundError: If genome_fname does not exist
        RuntimeError: If Jellyfish is not available
    """
    require_jellyfish()

    if not os.path.exists(genome_fname):
        raise FileNotFoundError(f"Genome file not found: {genome_fname}")

    # Ensure output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for k in range(min_k, max_k + 1):
        _print_k_progress(k, min_k, max_k)

        jf_file = f"{output_prefix}_{k}mer_all.jf"
        txt_file = f"{output_prefix}_{k}mer_all.txt"

        if not os.path.exists(txt_file):
            # Run jellyfish count
            count_cmd = [
                'jellyfish', 'count',
                '-m', str(k),
                '-s', '1000000',
                '-t', str(cpus),
                genome_fname,
                '-o', jf_file
            ]
            logger.debug(f"Running: {' '.join(count_cmd)}")
            subprocess.run(count_cmd, check=True, capture_output=True)

            # Dump to text file
            dump_cmd = ['jellyfish', 'dump', '-c', jf_file]
            logger.debug(f"Running: {' '.join(dump_cmd)}")
            with open(txt_file, 'w') as f_out:
                subprocess.run(dump_cmd, check=True, stdout=f_out)

        # Clean up intermediate .jf file
        if os.path.exists(jf_file):
            os.remove(jf_file)


def get_kmer_to_count_dict(f_in_name: str) -> Dict[str, int]:
    """
    Read a jellyfish dump file and return k-mer counts as a dictionary.

    Args:
        f_in_name: Path to jellyfish dump file (space-separated: kmer count)

    Returns:
        Dictionary mapping k-mer sequence to count
    """
    kmer_to_count = {}

    with open(f_in_name, 'r') as f_in:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) >= 2:
                kmer = parts[0]
                count = int(parts[1])
                kmer_to_count[kmer] = count

    return kmer_to_count


def get_primer_list_from_kmers(prefixes: List[str],
                                kmer_lengths: Optional[range] = None,
                                min_tm: float = 15.0,
                                max_tm: float = 55.0) -> List[str]:
    """
    Get all k-mers from jellyfish output files, filtered by Tm.

    Args:
        prefixes: List of path prefixes for jellyfish output files
        kmer_lengths: Range of k-mer lengths (default: 6-12)
        min_tm: Minimum melting temperature (default: 15.0)
        max_tm: Maximum melting temperature (default: 55.0)

    Returns:
        List of k-mer sequences that pass Tm filter
    """
    import melting

    primer_list = []

    if kmer_lengths is None:
        kmer_lengths = range(6, 13)

    for prefix in prefixes:
        for k in kmer_lengths:
            fpath = f"{prefix}_{k}mer_all.txt"
            if not os.path.exists(fpath):
                logger.warning(f"K-mer file not found: {fpath}")
                continue

            with open(fpath, 'r') as f_in:
                for line in f_in:
                    parts = line.strip().split()
                    if parts:
                        curr_kmer = parts[0]
                        try:
                            tm = melting.temp(curr_kmer)
                            if min_tm < tm < max_tm:
                                primer_list.append(curr_kmer)
                        except (ValueError, TypeError, KeyError) as e:
                            # Skip k-mers with invalid sequences (e.g. ambiguous bases
                            # cause KeyError in the melting library's complement lookup)
                            logger.debug(f"Skipping k-mer {curr_kmer}: Tm calculation failed ({e})")

    return primer_list
