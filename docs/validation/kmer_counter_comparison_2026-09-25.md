# Jellyfish against KMC3 on the genomes this repository ships

Measured 25 September 2026 on wMel (1.2 MB, 1 record) and *Drosophila*
(139 MB, 1,870 records) from `examples/wolbachia_pool_design/input/`, at
k = 12, 16 and 18, with 4 threads. Jellyfish 2.3.1, KMC 3.2.4 from bioconda.
Script: [`scripts/benchmarking/kmer_counter_comparison.py`](../../scripts/benchmarking/kmer_counter_comparison.py).
KMC is not a dependency of this project, so the script takes its location
rather than assuming one:

```bash
conda create -n kmc3 -c conda-forge -c bioconda kmc
KMC_BIN=$HOME/miniforge3/envs/kmc3/bin \
  python scripts/benchmarking/kmer_counter_comparison.py \
  examples/wolbachia_pool_design/input/drosophila.fna 18 4
```

Prompted by a report that COATswga counts with KMC3.
[tool_comparison.md](tool_comparison.md) leaves COATswga's counter
unspecified and attributes DSK to swga 1.0, so there is no contradiction, but
**neither claim is verified here**: this document measures the two counters,
not what any other tool does with them.

Jellyfish is invoked exactly as `neoswga.core.kmer_counter` invokes it,
including its hash-size guess of `file_size // 10`, so the figures describe
this project rather than an idealised run.

## Result

| genome | k | tool | count s | dump s | total s | peak RSS | database | text dump | distinct k-mers |
|---|---|---|---|---|---|---|---|---|---|
| wMel | 12 | jellyfish | 0.07 | 0.18 | **0.25** | 12 MB | 6.5 MB | 13.8 MB | 921,553 |
| wMel | 12 | KMC3 | 0.23 | 0.51 | 0.74 | 433 MB | 4.2 MB | 13.8 MB | 921,553 |
| wMel | 16 | jellyfish | 0.11 | 0.22 | **0.33** | 15 MB | 9.2 MB | 21.9 MB | 1,152,648 |
| wMel | 16 | KMC3 | 1.01 | 0.11 | 1.12 | 57 MB | 9.0 MB | 21.9 MB | 1,152,648 |
| wMel | 18 | jellyfish | 0.10 | 0.22 | **0.32** | 18 MB | 10.4 MB | 24.4 MB | 1,159,838 |
| wMel | 18 | KMC3 | 0.99 | 0.11 | 1.11 | 56 MB | 9.2 MB | 24.4 MB | 1,159,838 |
| *Drosophila* | 12 | jellyfish | 9.62 | 1.85 | 11.47 | 77 MB | 58 MB | 129 MB | 8,300,897 |
| *Drosophila* | 12 | KMC3 | 0.85 | 0.54 | **1.39** | 753 MB | 34 MB | 129 MB | 8,300,897 |
| *Drosophila* | 16 | jellyfish | 15.87 | 22.30 | 38.17 | 526 MB | 845 MB | 2,008 MB | 105,662,529 |
| *Drosophila* | 16 | KMC3 | 3.16 | 6.80 | **9.96** | 1,199 MB | 636 MB | 2,008 MB | 105,662,529 |
| *Drosophila* | 18 | jellyfish | 17.25 | 22.67 | 39.92 | 631 MB | 1,050 MB | 2,451 MB | 116,702,442 |
| *Drosophila* | 18 | KMC3 | 3.26 | 7.93 | **11.19** | 1,201 MB | 818 MB | 2,451 MB | 116,702,442 |

KMC3 is 3.6x faster at *Drosophila* k=18 and 8.2x at k=12. On wMel it is 2 to
3x SLOWER, because its fixed startup dominates a 1.2 MB genome.

## The comparison is only valid because both counted the same thing

The script refuses to report timings unless the two tools agree on the
distinct k-mer count, and they agree exactly in all six rows. That check is not
ceremony: the defaults make the two tools count different things, and a speed
comparison between different quantities is worthless.

## The memory argument for KMC3 is wrong, and it was mine

Before measuring, this project's reasoning held that KMC3 is the
memory-friendly option because it is disk-based, and that it would therefore
relieve the pressure that
[the Bloom filter path](../../CLAUDE.md) exists to relieve. Measured, the
opposite is true.

| | jellyfish | KMC3 |
|---|---|---|
| *Drosophila* k=18 peak RSS | **631 MB** | 1,201 MB |
| wMel k=18 peak RSS | **18 MB** | 56 MB |
| minimum accepted `-m` | n/a | 2 GB |

`kmc -m1` refuses outright: "min memory must be at least 2GB". With `-m2 -sm`
(strict) it completes at 1,201 MB, and `-m4` and `-m12` also peak at about
1,200 MB, so it uses what it needs rather than filling the cap. But the floor
is real: on a machine with under 2 GB to spare KMC will not start, while
jellyfish counted wMel in 18 MB.

So switching counters would not have helped the problem that motivated the
Bloom filter. It would have made it worse.

One figure is unexplained: KMC peaks at 433 MB on wMel k=12 but 56 MB at
k=16 and k=18 on the same genome. Recorded rather than explained.

## KMC's defaults are two traps

**`-ci2` excludes k-mers occurring once**, which on these genomes is almost
everything:

| | `-ci1` | `-ci2` (default) | lost |
|---|---|---|---|
| wMel k=18 | 1,159,838 | 36,037 | 96.9% |
| *Drosophila* k=18 | 116,702,442 | 5,163,663 | 95.6% |

A background table built with the default would report a host as almost
k-mer-free, and every candidate as specific. That is the silent-zero shape
this repository carries as Known Issues 5, 6, 13 and 15, reached through a
counter's default.

**`-cs255` saturates the counter at 255.** A host k-mer occurring 10,000 times
reads 255. That is Known Issue 7 exactly: an integer ceiling hiding a primer's
true host load, which there changed the delivered panel once the true figure
was visible.

Jellyfish needs `-C` for canonical counting; KMC is canonical by default and
`-b` turns it off. NeoSWGA already passes `-C`.

## The text dump dominates, and that is this project's choice

| *Drosophila* k=18 | count | dump | dump share |
|---|---|---|---|
| jellyfish | 17.25 s | 22.67 s | 57% |
| KMC3 | 3.26 s | 7.93 s | 71% |

Both tools spend most of their wall time writing 2.45 GB of text, because
`kmer_counter` dumps to `*_{k}mer_all.txt` and the rest of the package parses
that text. Reading the binary database instead would save more than changing
counters would, for either tool.

## What this means for NeoSWGA

Nothing is changed on the strength of this. Recorded so the question does not
have to be re-argued from recollection.

- **At the shipped k = 12 defaults the choice barely matters.** *Drosophila*
  takes 11.5 s with jellyfish and 1.4 s with KMC3; both are noise beside the
  rest of a design run.
- **KMC3's speed advantage is real at long k on large references**, which is
  the Bst regime (k 15-25) against a host.
- **It is not a memory fix.** See above.
- **The change that the measurement actually points at is the text dump**, and
  after that `kmc_tools` set operations, which could express "foreground
  k-mers not in background" without materialising either table in Python.
  Neither requires abandoning jellyfish.

Switching would also cost every path that reads `*_{k}mer_all.txt`, the
provenance sidecar, `_table_is_current` and the Jellyfish version guard.

## What was not measured

- hg38. The largest reference here is 139 MB; a 3.3 GB human background is
  where both tools' behaviour would matter most and neither was run on it.
- Disk I/O was not isolated. KMC writes to a working directory; these runs
  used local SSD and a slower volume would change its figures and not
  jellyfish's.
- One machine, one run per cell. The sub-second wMel figures are noise: a
  repeat of wMel k=12 gave KMC 0.24 s against the 0.74 s in the table, a 3x
  swing on the same command. Read the wMel rows as "both are instant", not as
  a ranking. The *Drosophila* ratios, 3.6x at k=18 and 8.2x at k=12, are far
  outside that variance.
- Whether COATswga uses KMC3, and what it does with it.
