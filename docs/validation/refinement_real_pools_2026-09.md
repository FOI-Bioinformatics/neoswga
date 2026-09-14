# Refinement on saved Prevotella pools, 2026-09-12

Network refinement remains the default. The 10,000-evaluation swap budget made
no replacements in any tested case. At 100,000 evaluations, swaps increased
coverage in five of six cases and matched it in one, but increased exact host
binding in the two 10-mer cases. Runtime was not consistently lower.

## Inputs and protocol

The available inputs were three saved Prevotella candidate pools: 2,000
10-mers, 2,000 11-mers, and 1,999 12-mers. The target length was 3,168,282 bp;
background was human chromosome 21, 46,709,983 bp. These are three pools on
one target, not three independent species. The earlier GC-tier run directories
were absent. Candidate pools were consumed as saved, without regenerating
filtering or efficacy rankings.

All cases used background-aware hybrid optimization, phi29 at 30 C, the saved
buffer settings/defaults, 3,000 bp selection reach, and a strict 3 bp pairwise
dimer limit. Each case ran in a fresh process, sequentially, with seed 42 and
one OpenMP/OpenBLAS thread. Network/swap order alternated between the initial
repeats. Two later swap-only sweeps used 100,000 evaluations; their timings
were not interleaved with network repeats. Swap search had a 10-second limit;
each complete worker had a 60-second timeout.

There were 36 completed runs: six pool/size combinations, three configurations
and two repeats. All delivered the requested size and zero violating pairs.
Panels were identical across repeats of each configuration. All measured
metrics and selected sequences are in [measurements.json](refinement_comparison_2026-09/measurements.json).

## Coverage and background

Coverage below is independently recomputed raw union coverage, with circular
foreground geometry, rather than the binned selection objective. Background
counts are exact-match sites; zero does not mean absence of near-match binding
or experimental background amplification. Occupancy-weighted coverage,
selectivity density, maximum gaps and binned coverage are retained in the data.

| Pool | Primers | Network coverage | Swap 10k coverage | Swap 100k coverage | Network host sites | Swap 100k host sites |
|---|---:|---:|---:|---:|---:|---:|
| 10-mer | 6 | 47.890% | 47.890% | 48.243% | 40 | 47 |
| 10-mer | 12 | 58.309% | 58.309% | 59.638% | 76 | 81 |
| 11-mer | 6 | 8.539% | 8.539% | 8.953% | 0 | 0 |
| 11-mer | 12 | 13.219% | 13.219% | 13.292% | 0 | 0 |
| 12-mer | 6 | 6.007% | 5.547% | 6.007% | 0 | 0 |
| 12-mer | 12 | 9.605% | 9.525% | 9.645% | 0 | 0 |

The largest raw coverage gain over network was 1.329 percentage points for
12 selected 10-mers, accompanied by an increase from 76 to 81 exact host sites.
The six-primer 10-mer panel gained 0.353 points and increased host sites from
40 to 47. This follows the implemented objective: coverage takes precedence,
and background breaks ties. It is not a guarantee of lower host binding.

## Optimizer cost

Medians of two runs, in seconds. Optimizer time includes stage preprocessing
and final network statistics; it excludes initial cache loading and independent
panel evaluation. Peak RSS is the fresh worker process maximum, including
cache loading and evaluation. It is not memory attributable to refinement alone.

| Pool | Primers | Network seconds | Swap 10k seconds | Swap 100k seconds | Network peak MiB | Swap 100k peak MiB |
|---|---:|---:|---:|---:|---:|---:|
| 10-mer | 6 | 1.508 | 1.417 | 1.807 | 401.0 | 420.8 |
| 10-mer | 12 | 2.582 | 1.434 | 2.422 | 403.2 | 431.0 |
| 11-mer | 6 | 0.536 | 1.030 | 1.170 | 607.0 | 594.6 |
| 11-mer | 12 | 0.702 | 1.041 | 1.212 | 608.9 | 612.4 |
| 12-mer | 6 | 0.516 | 0.992 | 1.043 | 1258.5 | 1202.5 |
| 12-mer | 12 | 0.644 | 1.019 | 1.171 | 1262.4 | 1268.5 |

Process wall times and cache times are recorded separately. Two repeats are
insufficient for a reliable speedup estimate, and filesystem caches remained
warm across processes. This measures the library optimizer, not end-to-end CLI
candidate generation, automatic prefiltering or alternative-set collection.

## Search-budget finding

All twelve 10k-budget runs stopped at the evaluation limit with zero swaps.
With about 2,000 candidates, a full pass requires roughly 12,000 or 24,000
candidate/removal checks at panel sizes six and twelve, respectively. Even
dimer-incompatible attempts consume the current evaluation budget.

At 100k, five configurations reached a one-for-one local optimum in 23,916 to
47,712 evaluations. The 12-primer 10-mer configuration accepted four swaps and
still reached the 100,000-evaluation limit. Both repeats agreed. No search
hit its time limit. Local optimality here concerns only the tested one-for-one
binned-coverage/background objective, not global or experimental optimality.

## Decision and next implementation

Keep network refinement and the current opt-in defaults unchanged. For an
explicit larger search on similar pools, pass `--swap-max-evaluations 100000`
and inspect both coverage and host binding. These results do not justify a
general default change.

The next useful algorithm change is to precompute each incoming candidate's
conflicts with the panel and skip impossible removals cheaply, then make the
budget report distinguish compatibility checks from feasible objective
evaluations. A separately configurable host-load constraint would also let
background-aware users prevent the coverage-first trade-off measured here.
Rebenchmark those changes before selecting a default.

## Reproduction and evidence

Run from the repository root with locally available candidate CSVs and HDF5 indexes:

```bash
python scripts/benchmarking/compare_refinement.py \
  --params tests/validation/genomes/params.json \
  --pools tests/validation/genomes/out10/step3_df.csv \
          tests/validation/genomes/out11/step3_df.csv \
          tests/validation/genomes/out12/step3_df.csv \
  --sizes 6 12 --repeats 2 --background-aware --timeout 60 \
  --output /tmp/refinement-comparison-new
```

For larger-budget checks, add `--modes swap --evaluations 100000` and use a
fresh output directory. The harness refuses to overwrite existing evidence.
It writes only to that directory and reads the supplied pools and indexes.
Core source, candidate CSV, params and index SHA-256 hashes are in each
manifest. Logs and complete measurements are retained under
`refinement_comparison_2026-09/`. Input CSVs and indexes are not vendored in
this evidence directory, so external reproduction requires those inputs.
