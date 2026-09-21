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

## Caveat

One panel, one reach, one assembly. The inflation scales with how many sites
fall within a reach of a join, so it rises with site count, with reach, and
with how fragmented the assembly is. A draft assembly of thousands of short
contigs would be worse than this by a wide margin, and a single-contig
reference is unaffected by construction.
