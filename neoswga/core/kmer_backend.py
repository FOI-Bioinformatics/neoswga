"""Which k-mer counter runs, and how it is invoked.

KMC3 is the default and jellyfish is selectable. The choice is recorded here
once so no caller builds a command line of its own.

**The parity flags are the substance of this module.** The two counters do not
compute the same quantity under their own defaults, and both KMC differences
are failures this repository already carries:

  `-ci1`        KMC excludes k-mers occurring once unless told otherwise.
                Measured on the shipped genomes at k=18 that default drops
                96.9% of wMel's k-mers and 95.6% of Drosophila's. A background
                counted that way reports a host as almost k-mer-free and every
                candidate as perfectly specific, which is the silent-zero shape
                of Known Issues 5, 6, 13 and 15.
  `-cs1000000`  KMC's counter ceiling is 255 by default, so a host k-mer
                occurring 10,000 times reads 255. That is Known Issue 7
                exactly: an integer ceiling hiding a primer's host load, which
                there changed the delivered panel once the true figure was
                visible.
  canonical     KMC is canonical unless `-b` is passed; jellyfish needs `-C`.

Measurements, including why KMC was chosen despite using more memory than
jellyfish: docs/validation/kmer_counter_comparison_2026-09-25.md
"""

import logging
import os
import shutil
import subprocess

logger = logging.getLogger(__name__)

#: KMC's counter ceiling. High enough that a host k-mer's true count survives;
#: the default of 255 does not.
_KMC_MAX_COUNT = 1_000_000

#: KMC refuses anything below this. Measured: `kmc -m1` exits with
#: "min memory must be at least 2GB".
_KMC_MIN_RAM_GB = 2


class KmerBackend:
    """A k-mer counter this project can drive."""

    name = "unnamed"
    #: What to tell a user who does not have it.
    install_hint = ""

    def binary(self, name):
        """Resolve a tool, honouring an explicit directory override."""
        override = os.environ.get(f"{self.name.upper()}_BIN")
        if override:
            candidate = os.path.join(override, name)
            return candidate if os.path.exists(candidate) else None
        return shutil.which(name)

    def available(self) -> bool:
        raise NotImplementedError

    def require_available(self) -> None:
        """Refuse with the install command rather than a bare OSError.

        The default counter changed, so a machine set up for the old default
        will arrive here, and "FileNotFoundError: kmc" does not tell anyone
        what to do about it.
        """
        if not self.available():
            raise RuntimeError(
                f"The {self.name} k-mer counter is not installed, and it is "
                f"the configured counter. {self.install_hint} Or select the "
                f'other counter with "kmer_counter" in params.json.'
            )

    def count(self, genome: str, k: int, out_prefix: str, threads: int = 4) -> None:
        raise NotImplementedError

    def text_table(self, prefix: str, k: int) -> str:
        """Where the text form of a table lives, for either backend."""
        return f"{prefix}_{k}mer_all.txt"


class KmcBackend(KmerBackend):
    name = "kmc"
    install_hint = "Install it with: conda install -c bioconda kmc"

    def available(self) -> bool:
        return self.binary("kmc") is not None and self.binary("kmc_dump") is not None

    def database(self, prefix: str, k: int) -> str:
        """KMC writes two files from this stem: .kmc_pre and .kmc_suf."""
        return f"{prefix}_{k}mer"

    def count(self, genome: str, k: int, out_prefix: str, threads: int = 4) -> None:
        self.require_available()
        db = self.database(out_prefix, k)
        os.makedirs(os.path.dirname(os.path.abspath(db)) or ".", exist_ok=True)
        workdir = f"{db}.tmp"
        os.makedirs(workdir, exist_ok=True)
        try:
            subprocess.run(
                [
                    self.binary("kmc"),
                    f"-k{k}",
                    "-ci1",  # count singletons; see the module docstring
                    f"-cs{_KMC_MAX_COUNT}",  # do not saturate at 255
                    "-fm",  # FASTA input
                    f"-t{threads}",
                    f"-m{_KMC_MIN_RAM_GB}",
                    genome,
                    db,
                    workdir,
                ],
                check=True,
                capture_output=True,
            )
        finally:
            shutil.rmtree(workdir, ignore_errors=True)


class JellyfishBackend(KmerBackend):
    name = "jellyfish"
    install_hint = (
        "Install it with: conda install -c bioconda kmer-jellyfish, or brew install jellyfish"
    )

    def available(self) -> bool:
        return self.binary("jellyfish") is not None

    def database(self, prefix: str, k: int) -> str:
        return f"{prefix}_{k}mer.jf"

    def _hash_size(self, genome: str) -> int:
        """Jellyfish sizes its hash upfront and resizes at 2-3x cost if the
        guess is low. This is the estimate `kmer_counter` has always used."""
        try:
            return max(1_000_000, os.path.getsize(genome) // 10)
        except OSError:
            return 1_000_000

    def count(self, genome: str, k: int, out_prefix: str, threads: int = 4) -> None:
        self.require_available()
        db = self.database(out_prefix, k)
        subprocess.run(
            [
                self.binary("jellyfish"),
                "count",
                "-m",
                str(k),
                "-s",
                str(self._hash_size(genome)),
                "-t",
                str(threads),
                "-C",  # canonical, matching KMC's default
                genome,
                "-o",
                db,
            ],
            check=True,
            capture_output=True,
        )


_BACKENDS = {"kmc": KmcBackend, "jellyfish": JellyfishBackend}

#: KMC3 is the default. It is faster at scale on the measured genomes; it also
#: uses more memory than jellyfish and will not run under 2 GB, which the
#: validation document states rather than leaving a reader to assume it won on
#: every axis.
DEFAULT_BACKEND = "kmc"


def select_backend(name: str | None = None) -> KmerBackend:
    """The counter to use, by name, falling back to the configured default."""
    if name is None:
        from neoswga.core import parameter

        name = getattr(parameter, "kmer_counter", None) or DEFAULT_BACKEND
    key = str(name).strip().lower()
    if key not in _BACKENDS:
        raise ValueError(
            f"Unknown k-mer counter {name!r}. Supported: " f"{', '.join(sorted(_BACKENDS))}."
        )
    return _BACKENDS[key]()
