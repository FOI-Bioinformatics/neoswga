# What an index without record geometry costs on a real assembly

Measured 21 September 2026 on the shipped Wolbachia design,
`examples/wolbachia_pool_design/work`: the delivered 12-oligo panel scored
against *Drosophila* (143,726,002 bp, 1,870 records) at 3 kb reach.

## Result

| Quantity | Value |
|---|---|
| host sites in the panel | 52 |
| `bg_coverage`, windows confined to their record | 0.00207162 |
| `bg_coverage`, windows unconfined | 0.00217080 |
| inflation | 14,255 bp, **+4.788% relative** |

Without record starts a coverage window anchored near the end of one contig
runs on into the next, crediting the panel with covering bases on a contig its
site is not on. On a 1,870-record assembly there are 1,869 places for that to
happen, and with only 52 sites several of them land near one.

The error is in the direction that makes a panel look WORSE on host coverage,
so it is not flattering. It is still wrong, and `max_host_coverage` is a
configurable panel limit, so a design held to one could be rejected for 4.8%
of coverage it does not have.

## Why the shipped example is now refused

Both indexes in that directory predate record geometry: neither carries
`#record_starts`, a format version, or a reference digest. `plan-pool` already
refused them before any of this work, because `require_record_metadata` has
always treated missing geometry as disqualifying. What changed on 2026-09-21
is that `optimize` applies the same standard.

The two references need different verdicts, which is why the check moved:

- **wMel** is a single record. There are no joins, so confining a window to
  its record and not confining it are the same operation. Its index is usable
  and is now accepted.
- **Drosophila** has 1,870 records and the figure above. Its index is refused.

`PositionCache` cannot draw that distinction, because which genome a prefix
belongs to is a relation only the resolved request knows; deciding it inside
the evaluator means reading `parameter.fg_genomes` and pairing it with the
prefixes the call was given, which is the defect that made a design refuse its
own index under `pytest -n 8`. So `reference_check.verify_index_geometry`
decides it against the manifest, beside the digest check.

Record counting reads header lines rather than loading the sequence. The
genome loader would hold 8.5 GB for hg38 to answer a question that counting
`>` settles in under a second.

## Remedy

Regenerate the affected index:

```
neoswga count-kmers -j params.json
neoswga filter -j params.json
```

Only references with more than one record need this. Every complete bacterial
genome in `tests/validation/genomes` is a single record and is unaffected;
*Drosophila*, hg38 (705 records) and Prevotella (2) are not.

## Where it does not matter: Prevotella

Measured the same way on the two-chromosome *Prevotella melaninogenica*
reference (3,168,282 bp, one join at 1,796,408) with the delivered 32-oligo
panel from `tests/validation/genomes`:

| Quantity | Value |
|---|---|
| target sites in the panel | 724 |
| `fg_coverage` confined | 0.52842455 |
| `fg_coverage` unconfined | 0.52842455 |
| inflation | **0 bp, 0.0000%** |

Exactly zero. With 724 sites in 3.17 Mb the region either side of the single
join is already covered from both directions, so confining the windows that
straddle it removes nothing that other windows do not supply.

Two data points, and together they say what the magnitude depends on: the
number of joins times the sparsity of the sites. One join with dense coverage
costs nothing; 1,869 joins with 52 sites costs 4.8%. Neither figure
generalises on its own.

## The pipeline can produce what the check demands

Worth stating because every index in this repository predates the feature, so
nothing here would have caught the opposite. A fresh `count-kmers`, `filter`,
`prepare-candidates`, `optimize` over a three-record reference writes format
version 2 with record starts `[0, 9000, 17000]` and all four steps complete.
`tests/integration/test_strict_design_pipeline.py` pins it, including that the
pipeline's own output satisfies `verify_index_geometry`.

Without that, the refusal would make multi-record references unusable rather
than merely requiring a recount.

## Caveat

Two panels, one reach each, two assemblies. The inflation scales with how many
sites fall within a reach of a join, so it rises with site count, with reach,
and with how fragmented the assembly is. A draft assembly of thousands of short
contigs would be worse than the Drosophila figure by a wide margin, and a
single-contig reference is unaffected by construction.
