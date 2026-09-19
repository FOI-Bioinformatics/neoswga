"""
Multi-genome k-mer counter for efficient k-mer frequency analysis.

Uses Jellyfish for fast k-mer counting across multiple genomes.
Jellyfish is a required dependency - the module will raise an error if not available.
"""

import concurrent.futures
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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


def _adaptive_hash_size(genome_path: str) -> int:
    """Estimate Jellyfish hash size from genome file size.

    Avoids the default 1M hash, which causes Jellyfish to resize
    internally (2-3x slowdown) for large genomes.  For small genomes
    the minimum of 1M is retained.

    Args:
        genome_path: Path to the FASTA file.

    Returns:
        Hash table size for Jellyfish ``-s`` flag.
    """
    try:
        file_size = os.path.getsize(genome_path)
        return max(1_000_000, file_size // 10)
    except OSError:
        return 1_000_000


def check_jellyfish_available() -> bool:
    """Check if jellyfish is available in PATH."""
    return shutil.which("jellyfish") is not None


def get_jellyfish_version() -> Optional[str]:
    """Parse the Jellyfish version from ``jellyfish --version``.

    Returns:
        Version string (e.g. '2.3.1'), or None if unavailable.
    """
    if not check_jellyfish_available():
        return None
    try:
        result = subprocess.run(
            ["jellyfish", "--version"], capture_output=True, text=True, timeout=10
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
        major = version.split(".")[0]
        if major == "1":
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
            self._temp_dir = tempfile.mkdtemp(prefix="neoswga_kmer_")
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
                if not line.startswith(">"):
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
        output_prefix = os.path.join(work_dir, genome_name)
        jf_file = os.path.join(work_dir, f"{genome_name}_{k}mer.jf")
        txt_file = os.path.join(work_dir, f"{genome_name}_{k}mer_all.txt")

        if not _table_is_current(output_prefix, fasta_path, k):
            # Run jellyfish count
            count_cmd = [
                "jellyfish",
                "count",
                "-m",
                str(k),
                "-s",
                str(_adaptive_hash_size(fasta_path)),
                "-t",
                str(self.cpus),
                "-C",  # Canonical k-mers (both strands)
                fasta_path,
                "-o",
                jf_file,
            ]
            subprocess.run(count_cmd, check=True, capture_output=True)

            # Dump to text
            dump_cmd = ["jellyfish", "dump", "-c", jf_file]
            with open(txt_file, "w") as f:
                subprocess.run(dump_cmd, check=True, stdout=f)

            _write_table_provenance(output_prefix, fasta_path, k)

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
        kmer = sequence[i : i + k]
        # Only count canonical DNA k-mers
        if all(base in "ACGT" for base in kmer):
            counts[kmer] += 1

    return dict(counts)


# =============================================================================
# Standalone Functions (replacements for deprecated kmer.py)
# =============================================================================

_JELLYFISH_TIMEOUT = 3600  # 1 hour per k-value


def table_provenance_path(output_prefix: str, k: int) -> str:
    """Where the record of what a k-mer table was counted from lives."""
    return f"{output_prefix}_{k}mer_all.provenance.json"


# Named in every provenance record, so one written under the partial hash this
# replaced is recognisable as UNKNOWN rather than reported as a mismatch.
DIGEST_ALGORITHM = "sha256-full"

# path -> (size, mtime_ns, digest). Computed once per input per run, which is
# what makes a full digest affordable: `run_jellyfish` fingerprints once per k,
# so seven k values used to mean seven passes over the genome.
_DIGEST_CACHE: Dict[str, Tuple[int, int, str]] = {}


def genome_fingerprint(genome_fname: str) -> str:
    """SHA-256 of the whole genome file, computed once per input per run.

    This hashed the size, the first 1 MB and the last 1 MB, on the stated
    grounds that a full digest of hg38 "costs several seconds on every
    invocation". Finding F6: a substitution anywhere in the middle of a file
    over 2 MB left it unchanged, so a same-length consensus or a
    sample-specific assembly reused the previous genome's counts, index and
    inventory -- silently, and through the whole design. That is the case the
    sidecar exists to catch, missed in exactly the situation where the file
    does not change size.

    The cost objection does not survive measurement. SHA-256 runs at about
    2.5 GB/s, so hg38 is around a second, and the cache below makes it once per
    run rather than once per k. The counting step reads the whole file anyway.

    The cache is keyed on size and modification time as well as the path,
    because caching on the path alone would stop looking at the one file the
    guard exists to notice changing.
    """
    stat = os.stat(genome_fname)
    key = (stat.st_size, stat.st_mtime_ns)
    cached = _DIGEST_CACHE.get(genome_fname)
    if cached is not None and cached[:2] == key:
        return cached[2]

    digest = hashlib.sha256()
    with open(genome_fname, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    value = digest.hexdigest()
    _DIGEST_CACHE[genome_fname] = (key[0], key[1], value)
    return value


def table_provenance_is_comparable(output_prefix: str, k: int) -> bool:
    """Whether this table's record can be compared with today's digest.

    A record written under the partial hash carries a fingerprint that cannot
    be checked against a full one. That is UNKNOWN, not wrong, and the
    distinction matters: step 1 recounts an unknown table once, while step 2
    REFUSES a table it believes was counted from another genome. Treating the
    first as the second would refuse every working data directory on upgrade,
    over a fact nobody measured.
    """
    path = table_provenance_path(output_prefix, k)
    if not os.path.exists(path):
        return False
    try:
        with open(path) as fh:
            record = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return False
    return record.get("digest_algorithm") == DIGEST_ALGORITHM


def _table_is_current(output_prefix: str, genome_fname: str, k: int) -> bool:
    """Whether an existing table was counted from this genome.

    The existence check this replaces keyed on the output prefix, which comes
    from params.json and not from the genome. Repointing `fg_genomes` at a new
    assembly while leaving `fg_prefixes` alone therefore reused the previous
    organism's counts through the whole design, silently.
    """
    txt_file = f"{output_prefix}_{k}mer_all.txt"
    if not os.path.exists(txt_file):
        return False

    record_path = table_provenance_path(output_prefix, k)
    if not os.path.exists(record_path):
        logger.warning(
            f"{txt_file} exists but has no provenance record, so it cannot be "
            f"matched to {genome_fname}. Recounting. Tables written before this "
            f"check was added will each be recounted once."
        )
        return False

    try:
        with open(record_path) as fh:
            record = json.load(fh)
    except (OSError, ValueError) as exc:
        logger.warning(f"Could not read {record_path} ({exc}). Recounting.")
        return False

    if record.get("digest_algorithm") != DIGEST_ALGORITHM:
        # Written under the partial hash, whose value cannot be compared with
        # this one. Unknown rather than stale: recount once, as a table with no
        # sidecar at all already does.
        logger.info(
            "%s was fingerprinted with an older algorithm; recounting once so "
            "the record can be trusted.",
            os.path.basename(genome_fname),
        )
        return False
    if record.get("fingerprint") != genome_fingerprint(genome_fname):
        logger.warning(
            f"{txt_file} was counted from a different genome "
            f"({record.get('genome', 'unrecorded')}), not from {genome_fname}. "
            f"Recounting. Reusing it would have built this design from the "
            f"other genome's k-mer counts."
        )
        return False

    return True


def _write_table_provenance(output_prefix: str, genome_fname: str, k: int) -> None:
    """Record what this table was counted from, for the next run's check."""
    record = {
        "genome": os.path.abspath(genome_fname),
        "fingerprint": genome_fingerprint(genome_fname),
        "digest_algorithm": DIGEST_ALGORITHM,
        "k": k,
    }
    with open(table_provenance_path(output_prefix, k), "w") as fh:
        json.dump(record, fh, indent=2)


def _run_jellyfish_for_k(
    output_prefix: str, genome_fname: str, k: int, cpus: int, hash_size: int
) -> None:
    """Run Jellyfish count + dump for a single k-value.

    Args:
        output_prefix: Output path prefix for result files.
        genome_fname: Path to FASTA file.
        k: K-mer length.
        cpus: Number of threads for Jellyfish.
        hash_size: Jellyfish hash table size (``-s`` flag).

    Raises:
        RuntimeError: If jellyfish count or dump fails.
        subprocess.TimeoutExpired: If a command exceeds the timeout.
    """
    jf_file = f"{output_prefix}_{k}mer_all.jf"
    txt_file = f"{output_prefix}_{k}mer_all.txt"

    if not _table_is_current(output_prefix, genome_fname, k):
        count_cmd = [
            "jellyfish",
            "count",
            "-m",
            str(k),
            "-s",
            str(hash_size),
            "-t",
            str(cpus),
            "-C",  # Canonical k-mers (count both strands together)
            genome_fname,
            "-o",
            jf_file,
        ]
        logger.debug(f"Running: {' '.join(count_cmd)}")
        try:
            result = subprocess.run(
                count_cmd,
                check=True,
                capture_output=True,
                text=True,
                timeout=_JELLYFISH_TIMEOUT,
            )
            if result.stderr:
                logger.debug(f"jellyfish count stderr: {result.stderr.strip()}")
        except subprocess.CalledProcessError as e:
            stderr_msg = (e.stderr or "").strip()
            raise RuntimeError(
                f"jellyfish count failed for k={k}: {stderr_msg}\n"
                f"Command: {' '.join(count_cmd)}"
            ) from e
        except subprocess.TimeoutExpired:
            # Clean up partial output
            if os.path.exists(jf_file):
                os.remove(jf_file)
            raise

        dump_cmd = ["jellyfish", "dump", "-c", jf_file]
        logger.debug(f"Running: {' '.join(dump_cmd)}")
        try:
            with open(txt_file, "w") as f_out:
                subprocess.run(
                    dump_cmd,
                    check=True,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=_JELLYFISH_TIMEOUT,
                )
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            # Clean up partial txt file
            if os.path.exists(txt_file):
                os.remove(txt_file)
            raise

        _write_table_provenance(output_prefix, genome_fname, k)

    if os.path.exists(jf_file):
        os.remove(jf_file)


def run_jellyfish(
    genome_fname: str, output_prefix: str, min_k: int = 6, max_k: int = 12, cpus: int = 4
) -> None:
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

    hash_size = _adaptive_hash_size(genome_fname)
    num_k = max_k - min_k + 1
    max_workers = min(num_k, max(1, (os.cpu_count() or 1) // max(cpus, 1)))

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _run_jellyfish_for_k, output_prefix, genome_fname, k, cpus, hash_size
            ): k
            for k in range(min_k, max_k + 1)
        }
        errors = []
        for future in concurrent.futures.as_completed(futures):
            k = futures[future]
            try:
                future.result()
            except Exception as e:
                logger.error(f"Jellyfish failed for k={k}: {e}")
                errors.append((k, e))
                continue
            _print_k_progress(k, min_k, max_k)

    if errors:
        failed_ks = ", ".join(str(k) for k, _ in errors)
        raise RuntimeError(
            f"Jellyfish failed for k-mer lengths: {failed_ks}. " f"First error: {errors[0][1]}"
        )

    # Validate output files
    missing = []
    empty = []
    for k in range(min_k, max_k + 1):
        txt_file = f"{output_prefix}_{k}mer_all.txt"
        if not os.path.exists(txt_file):
            missing.append(txt_file)
        elif os.path.getsize(txt_file) == 0:
            empty.append(txt_file)

    if missing:
        raise RuntimeError(f"Jellyfish output files missing after counting: {missing}")
    if empty:
        logger.warning(
            f"Empty k-mer files (genome may be shorter than k): "
            f"{[os.path.basename(f) for f in empty]}"
        )


def get_kmer_to_count_dict(f_in_name: str) -> Dict[str, int]:
    """
    Read a jellyfish dump file and return k-mer counts as a dictionary.

    Args:
        f_in_name: Path to jellyfish dump file (space-separated: kmer count)

    Returns:
        Dictionary mapping k-mer sequence to count
    """
    kmer_to_count = {}

    with open(f_in_name, "r") as f_in:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) >= 2:
                kmer = parts[0]
                count = int(parts[1])
                kmer_to_count[kmer] = count

    return kmer_to_count


_UNAMBIGUOUS_BASES = frozenset("ACGT")


def _gc_content(seq: str) -> float:
    """Fast GC fraction for a DNA sequence (no validation)."""
    gc = 0
    for c in seq:
        if c == "G" or c == "C" or c == "g" or c == "c":
            gc += 1
    return gc / len(seq) if seq else 0.0


def get_primer_list_from_kmers(
    prefixes: List[str],
    kmer_lengths: Optional[range] = None,
    min_tm: float = 15.0,
    max_tm: float = 55.0,
    wide_tm_margin: float = 2.0,
    gc_min: float = 0.10,
    gc_max: float = 0.90,
    conditions=None,
) -> List[str]:
    """
    Get all k-mers from jellyfish output files, filtered by GC content and Tm.

    Applies a fast GC pre-filter before the Tm calculation.

    The Tm is `ReactionConditions.calculate_effective_tm`, which is the same
    call `filter.filter_extra` makes when it decides which candidates survive.
    It used to be `melting_temp.temp`, whose docstring records a deliberate bug
    kept "so that Tm values remain consistent with the RF model". That model was
    retired from the default path on 2026-09-05, and the bug outlived its
    purpose: measured on 20,000 random 12-mers under the E. coli tier's
    conditions the shim reads 10.12 C higher (sd 0.30 C).

    With a systematic +10 C offset and the old symmetric 15 C margin, the
    window this loader applied was `[min_tm - 25, max_tm + 5]` in true-Tm terms.
    That lost nothing on plain phi29 and grew with additives -- 9.6% of k=12
    candidates that pass the real window were discarded under DMSO 10% plus
    betaine 1.5 M, silently.

    The margin is now 2 C of stated headroom rather than a correction for
    disagreement between two estimators, because there is only one estimator.

    Args:
        prefixes: List of path prefixes for jellyfish output files
        kmer_lengths: Range of k-mer lengths (default: 6-12)
        min_tm: Minimum melting temperature (default: 15.0)
        max_tm: Maximum melting temperature (default: 55.0)
        wide_tm_margin: Headroom added to both sides of the Tm window
            (default: 2.0)
        gc_min: Minimum GC fraction for pre-filter (default: 0.10)
        gc_max: Maximum GC fraction for pre-filter (default: 0.90)
        conditions: ReactionConditions to measure Tm under. Defaults to a
            plain reaction, for library callers that reach this without a
            configured run.

    Returns:
        List of k-mer sequences that pass GC and Tm pre-filters
    """
    from neoswga.core.reaction_conditions import ReactionConditions

    if conditions is None:
        conditions = ReactionConditions()

    primer_list = []
    gc_rejected = 0
    ambiguous_rejected = 0

    if kmer_lengths is None:
        kmer_lengths = range(6, 13)

    wide_min = min_tm - wide_tm_margin
    wide_max = max_tm + wide_tm_margin

    for prefix in prefixes:
        for k in kmer_lengths:
            fpath = f"{prefix}_{k}mer_all.txt"
            if not os.path.exists(fpath):
                logger.warning(f"K-mer file not found: {fpath}")
                continue

            with open(fpath, "r") as f_in:
                for line in f_in:
                    parts = line.strip().split()
                    if not parts:
                        continue
                    curr_kmer = parts[0]
                    # Fast GC pre-filter (avoids the Tm calculation)
                    gc = _gc_content(curr_kmer)
                    if gc < gc_min or gc > gc_max:
                        gc_rejected += 1
                        continue
                    # `calculate_effective_tm` warns and substitutes penalty
                    # values for an unknown base rather than raising, so an
                    # ambiguous k-mer would otherwise be admitted with a
                    # meaningless number.
                    if not set(curr_kmer.upper()) <= _UNAMBIGUOUS_BASES:
                        ambiguous_rejected += 1
                        continue
                    try:
                        tm = conditions.calculate_effective_tm(curr_kmer)
                    except (ValueError, TypeError, KeyError) as e:
                        logger.debug(f"Skipping k-mer {curr_kmer}: Tm calculation failed ({e})")
                        continue
                    if wide_min < tm < wide_max:
                        primer_list.append(curr_kmer)

    if gc_rejected > 0:
        logger.info(
            f"GC pre-filter removed {gc_rejected} k-mers outside {gc_min:.0%}-{gc_max:.0%} range"
        )
    if ambiguous_rejected > 0:
        logger.info(f"Removed {ambiguous_rejected} k-mers containing ambiguous bases")

    return primer_list
