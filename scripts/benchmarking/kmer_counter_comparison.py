"""Jellyfish against KMC3 on the genomes this repository ships.

Measures wall time, peak resident memory and output size for each tool at each
k, and CHECKS THE TWO AGREE on the distinct canonical k-mer count. A speed
comparison between tools computing different quantities is worthless, and the
defaults do compute different quantities:

  KMC `-ci2`   excludes k-mers occurring once. Jellyfish counts them. Without
               `-ci1` KMC reports far fewer distinct k-mers, and on a genome
               where most k-mers are unique that is most of them.
  KMC `-cs255` saturates the counter at 255. A host k-mer occurring 10,000
               times reads 255, which is the silent-saturation family this
               repository already carries as Known Issue 7. `-cs` is raised
               here so counts are comparable.
  KMC default  IS canonical; `-b` turns it off. Jellyfish needs `-C` to match.

Jellyfish is invoked the way neoswga.core.kmer_counter invokes it, including
its hash size guess of file_size // 10, so the figures describe this project
rather than an idealised run.

Usage:
    python scripts/benchmarking/kmer_counter_comparison.py <genome.fna> <k> [threads]
"""

import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time

# macOS reports ru_maxrss in bytes, Linux in kilobytes.
_RSS_SCALE = 1 if sys.platform == "darwin" else 1024

def _tool(name, env_var):
    """Resolve a binary from PATH, or from an explicit override.

    KMC is not a dependency of this project and is usually installed in its
    own environment, so the path is not assumed. Set the override rather than
    editing this file:

        conda create -n kmc3 -c conda-forge -c bioconda kmc
        KMC_BIN=$CONDA_PREFIX/envs/kmc3/bin python scripts/benchmarking/...
    """
    override = os.environ.get(env_var)
    if override:
        candidate = os.path.join(override, name)
        if os.path.exists(candidate):
            return candidate
        raise SystemExit(f"{env_var} is set to {override!r} but {name} is not there")
    found = shutil.which(name)
    if found:
        return found
    raise SystemExit(
        f"{name} not found on PATH. Install it, or set {env_var} to the "
        f"directory holding it. For KMC: "
        f"conda create -n kmc3 -c conda-forge -c bioconda kmc"
    )


JELLYFISH = _tool("jellyfish", "JELLYFISH_BIN")
KMC = _tool("kmc", "KMC_BIN")
KMC_DUMP = _tool("kmc_dump", "KMC_BIN")


def _run(cmd, **kw):
    """Run a command, returning (seconds, peak child RSS in bytes)."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * _RSS_SCALE
    start = time.perf_counter()
    subprocess.run(cmd, check=True, capture_output=True, **kw)
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * _RSS_SCALE
    return elapsed, max(after, before)


def _run_callable(fn):
    """Same measurement as `_run`, for work driven through the package."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * _RSS_SCALE
    start = time.perf_counter()
    fn()
    elapsed = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * _RSS_SCALE
    return elapsed, max(after, before)


def _distinct(path):
    with open(path, "rb") as fh:
        return sum(1 for _ in fh)


def run_jellyfish(genome, k, threads, workdir):
    jf = os.path.join(workdir, "j.jf")
    txt = os.path.join(workdir, "j.txt")
    size = max(1_000_000, os.path.getsize(genome) // 10)  # neoswga's own guess
    t_count, rss = _run(
        [JELLYFISH, "count", "-m", str(k), "-s", str(size), "-t", str(threads),
         "-C", genome, "-o", jf]
    )
    start = time.perf_counter()
    with open(txt, "w") as fh:
        subprocess.run([JELLYFISH, "dump", "-c", jf], check=True, stdout=fh)
    t_dump = time.perf_counter() - start
    return {
        "count_s": t_count, "dump_s": t_dump, "total_s": t_count + t_dump,
        "peak_rss": rss, "distinct": _distinct(txt),
        "db_bytes": os.path.getsize(jf), "txt_bytes": os.path.getsize(txt),
        "hash_size": size,
    }


def run_kmc(genome, k, threads, workdir):
    """Counted through the SHIPPING backend, not a command written here.

    This used to build its own command line with `-m4`. KMC uses roughly the
    memory you allow it -- measured on hg38 at k=12, 3.36 GB at `-m4` against
    1.97 GB at the `-m2` the backend passes -- so the benchmark was reporting
    a memory figure no production run would produce. A benchmark that drifts
    from the code it measures is worse than none.
    """
    sys.path.insert(0, os.getcwd())
    from neoswga.core.kmer_backend import KmcBackend

    backend = KmcBackend()
    prefix = os.path.join(workdir, "kmcrun")
    db = backend.database(prefix, k)
    txt = os.path.join(workdir, "k.txt")
    t_count, rss = _run_callable(lambda: backend.count(genome, k, prefix, threads=threads))
    start = time.perf_counter()
    subprocess.run([backend.binary("kmc_dump"), db, txt], check=True, capture_output=True)
    t_dump = time.perf_counter() - start
    db_bytes = sum(
        os.path.getsize(db + suffix) for suffix in (".kmc_pre", ".kmc_suf")
    )
    return {
        "count_s": t_count, "dump_s": t_dump, "total_s": t_count + t_dump,
        "peak_rss": rss, "distinct": _distinct(txt),
        "db_bytes": db_bytes, "txt_bytes": os.path.getsize(txt),
        "hash_size": None,
    }


def main():
    genome, k = sys.argv[1], int(sys.argv[2])
    threads = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    name = os.path.basename(genome)

    results = {}
    for label, fn in (("jellyfish", run_jellyfish), ("kmc3", run_kmc)):
        with tempfile.TemporaryDirectory(prefix=f"{label}_") as workdir:
            results[label] = fn(genome, k, threads, workdir)

    j, m = results["jellyfish"], results["kmc3"]
    agree = j["distinct"] == m["distinct"]
    print(f"genome={name} k={k} threads={threads}")
    for label, r in results.items():
        print(
            f"  {label:10s} count={r['count_s']:7.2f}s dump={r['dump_s']:7.2f}s "
            f"total={r['total_s']:7.2f}s peakRSS={r['peak_rss']/1e6:8.1f}MB "
            f"db={r['db_bytes']/1e6:8.1f}MB txt={r['txt_bytes']/1e6:8.1f}MB "
            f"distinct={r['distinct']:,}"
        )
    print(f"  AGREE on distinct k-mers: {agree}"
          + ("" if agree else f"  jellyfish={j['distinct']:,} kmc={m['distinct']:,}"))
    if not agree:
        raise SystemExit("the two tools counted different things; the timings are not comparable")


if __name__ == "__main__":
    main()
