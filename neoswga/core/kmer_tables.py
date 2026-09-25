"""One place every consumer asks about a k-mer table.

Callers name a prefix and a k, never a file. That is what lets the table live
as a binary database, a text dump, or both, and lets a question be answered by
whichever is cheapest without the caller knowing.

Two access patterns cover every consumer in the package:

  `counts_for`  the counts of a KNOWN list of k-mers. This is what `filter`
                asks, and it is a set operation rather than a scan. Measured
                on Drosophila at k=18 with 2,000 candidates: 2.5 s through
                `kmc_tools simple ... intersect` against 17.9 s to write a
                2.3 GB text table and stream it into Python.
  `iter_table`  every k-mer with its count, for the callers that genuinely
                need all of them.

Absence and failure are kept apart. A prefix that was never counted RAISES;
reporting zeros there would say every candidate is absent from the background,
which is the most permissive answer available and the silent-zero shape this
repository carries as Known Issues 5, 6, 13 and 15.

Measurements: docs/validation/kmer_counter_comparison_2026-09-25.md
"""

import logging
import os
import subprocess
import tempfile
from collections.abc import Iterator, Sequence

from neoswga.core.kmer_backend import JellyfishBackend, KmcBackend

logger = logging.getLogger(__name__)


def text_table_path(prefix: str, k: int) -> str:
    return f"{prefix}_{k}mer_all.txt"


def _kmc_database(prefix: str, k: int) -> str | None:
    """The KMC database stem for this prefix, if both of its files exist."""
    stem = KmcBackend().database(prefix, k)
    if os.path.exists(stem + ".kmc_pre") and os.path.exists(stem + ".kmc_suf"):
        return stem
    return None


def _jellyfish_database(prefix: str, k: int) -> str | None:
    path = JellyfishBackend().database(prefix, k)
    return path if os.path.exists(path) else None


def table_exists(prefix: str, k: int) -> bool:
    """Whether this prefix has been counted at this k, in any form."""
    return bool(
        os.path.exists(text_table_path(prefix, k))
        or _kmc_database(prefix, k)
        or _jellyfish_database(prefix, k)
    )


def _stream_pairs(cmd: list[str]) -> Iterator[tuple[str, int]]:
    """Yield (k-mer, count) from a dump command's stdout.

    Both counters are driven the same way, which is why this is one function:
    reusing a single `backend` name for two different classes is what mypy
    objected to when they were written out separately.

    A consumer that stops early leaves the writer with a closed pipe. That is
    the expected shape of an early exit, not an error, so the reader is closed
    and the child reaped without letting BrokenPipeError escape.
    """
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            parts = line.split()
            if len(parts) >= 2:
                yield parts[0], int(parts[1])
    finally:
        proc.stdout.close()
        proc.wait()


def discover_prefixes(directory: str, k: int) -> list[str]:
    """Prefixes in this directory that have a table at this k, in any form.

    Auto-discovery used to glob `*_{k}mer_all.txt` alone, so a directory
    counted into KMC databases looked empty and every reference in it went
    unregistered.
    """
    import glob as _glob

    patterns = {
        f"*_{k}mer_all.txt": f"_{k}mer_all.txt",
        f"*_{k}mer.kmc_pre": f"_{k}mer.kmc_pre",
        f"*_{k}mer.jf": f"_{k}mer.jf",
    }
    prefixes: list[str] = []
    for pattern, suffix in patterns.items():
        for path in _glob.glob(os.path.join(directory, pattern)):
            prefix = path[: -len(suffix)]
            if prefix not in prefixes and table_exists(prefix, k):
                prefixes.append(prefix)
    return sorted(prefixes)


def table_files(prefix: str, k: int) -> list[str]:
    """Every file making up this prefix's table at this k, in any form.

    A KMC table is TWO files. Linking or copying only one leaves a half a
    database, which `table_exists` correctly rejects but which would look like
    a corrupted count rather than an incomplete copy.
    """
    found = []
    stem = KmcBackend().database(prefix, k)
    if os.path.exists(stem + ".kmc_pre") and os.path.exists(stem + ".kmc_suf"):
        found.extend([stem + ".kmc_pre", stem + ".kmc_suf"])
    jf = JellyfishBackend().database(prefix, k)
    if os.path.exists(jf):
        found.append(jf)
    text = text_table_path(prefix, k)
    if os.path.exists(text):
        found.append(text)
    return found


def link_table(src_prefix: str, dst_prefix: str, k: int) -> bool:
    """Symlink a counted table from one prefix to another.

    Used when the genome library already holds a count. It links whatever form
    the source has, rather than a text file by name: linking
    `{src}_{k}mer_all.txt` when the source is a KMC database links nothing and
    the caller then recounts a genome it already has.

    Returns whether anything was linked. An existing destination is left
    alone.
    """
    linked = False
    for source in table_files(src_prefix, k):
        destination = source.replace(src_prefix, dst_prefix, 1)
        if source == destination or os.path.exists(destination):
            continue
        os.symlink(source, destination)
        linked = True
    return linked


def _counts_from_stream(prefix: str, k: int) -> Iterator[tuple[str, int]]:
    """Every (k-mer, count) pair, from whichever form is present and readable.

    A database is only readable with its counter's dump tool. A directory
    counted with KMC on one machine and read on another without KMC would
    otherwise fail inside `Popen` with a message about `None`. It falls back
    to a text dump when one is present, and otherwise says what to install.
    """
    database = _kmc_database(prefix, k)
    if database:
        dump = KmcBackend().binary("kmc_dump")
        if dump:
            yield from _stream_pairs([dump, database, "/dev/stdout"])
            return
        if not os.path.exists(text_table_path(prefix, k)):
            raise RuntimeError(
                f"The {k}-mer table for {prefix!r} is a KMC database, and KMC is "
                f"not installed here, so it cannot be read. Install it with: "
                f"conda install -c bioconda kmc"
            )

    jf = _jellyfish_database(prefix, k)
    if jf:
        tool = JellyfishBackend().binary("jellyfish")
        if tool:
            yield from _stream_pairs([tool, "dump", "-c", jf])
            return
        if not os.path.exists(text_table_path(prefix, k)):
            raise RuntimeError(
                f"The {k}-mer table for {prefix!r} is a jellyfish database, and "
                f"jellyfish is not installed here, so it cannot be read. Install "
                f"it with: conda install -c bioconda kmer-jellyfish"
            )

    with open(text_table_path(prefix, k)) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                yield parts[0], int(parts[1])


def iter_table(prefix: str, k: int) -> Iterator[tuple[str, int]]:
    """Every k-mer and its count. Raises if this prefix was never counted."""
    _require_counted(prefix, k)
    yield from _counts_from_stream(prefix, k)


def _require_counted(prefix: str, k: int) -> None:
    if table_exists(prefix, k):
        return
    raise FileNotFoundError(
        f"No {k}-mer table for {prefix!r}: neither a KMC database, a jellyfish "
        f"database, nor {text_table_path(prefix, k)!r}. Run 'neoswga "
        f"count-kmers -j params.json' to build one. Reporting zero counts here "
        f"would say every candidate is absent from this reference, which is "
        f"the most permissive answer available."
    )


def _counts_by_scanning(prefix: str, k: int, wanted: set[str]) -> dict[str, int]:
    """Stream the table and keep the requested k-mers.

    Stops once every one has been seen. That early exit is why this is not as
    slow as the full scan whenever the candidates are all present; when some
    are absent, which is the usual case for a background, it reads everything.
    """
    found: dict[str, int] = {}
    for kmer, count in _counts_from_stream(prefix, k):
        if kmer in wanted:
            found[kmer] = count
            if len(found) == len(wanted):
                break
    return found


def _counts_by_intersection(
    database: str, k: int, wanted: set[str], workdir: str
) -> dict[str, int]:
    """Ask KMC for the counts of exactly these k-mers.

    `-ocleft` keeps the counters of the FIRST database. Without it the output
    carries the candidate database's counts, which are all 1, and every
    background count would read 1 -- a wrong answer that looks plausible.
    """
    backend = KmcBackend()
    candidates_fasta = os.path.join(workdir, "candidates.fna")
    with open(candidates_fasta, "w") as fh:
        for index, kmer in enumerate(sorted(wanted)):
            fh.write(f">c{index}\n{kmer}\n")

    backend.count(candidates_fasta, k, os.path.join(workdir, "cand_prefix"))
    candidate_db = backend.database(os.path.join(workdir, "cand_prefix"), k)

    result_db = os.path.join(workdir, "hits")
    subprocess.run(
        [
            backend.binary("kmc_tools"),
            "-hp",
            "simple",
            database,
            candidate_db,
            "intersect",
            result_db,
            "-ocleft",
        ],
        check=True,
        capture_output=True,
    )

    dumped = os.path.join(workdir, "hits.txt")
    subprocess.run(
        [backend.binary("kmc_dump"), result_db, dumped],
        check=True,
        capture_output=True,
    )

    found: dict[str, int] = {}
    with open(dumped) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2:
                found[parts[0]] = int(parts[1])
    return found


def counts_for(prefix: str, k: int, kmers: Sequence[str]) -> dict[str, int]:
    """Counts of exactly these k-mers, zero for any the table does not hold.

    Every requested k-mer appears in the result. A caller cannot tell "absent"
    from "not asked" by inspecting the keys, which is the distinction the
    background frequency gate depends on.
    """
    if not kmers:
        return {}

    # A k-mer of the wrong length cannot be in a table of this k. Asking KMC
    # about it would fail the whole batch rather than answer the rest.
    usable = {kmer for kmer in kmers if len(kmer) == k}
    if not usable:
        _require_counted(prefix, k)
        return dict.fromkeys(kmers, 0)

    _require_counted(prefix, k)

    database = _kmc_database(prefix, k)
    if database and KmcBackend().available():
        with tempfile.TemporaryDirectory(prefix="neoswga_kmer_lookup_") as workdir:
            found = _counts_by_intersection(database, k, usable, workdir)
    else:
        found = _counts_by_scanning(prefix, k, usable)

    return {kmer: found.get(kmer, 0) for kmer in kmers}
