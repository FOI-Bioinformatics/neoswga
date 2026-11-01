# One site, two mismatch classes, two occupancies

Found and fixed 1 November 2026, on `main` at `db3d2fc` plus the change this
document records, during a chemistry and design-flow audit. Reproduce with:

```
PYTHONPATH=. python docs/validation/mismatch_class_overlap_2026-11-01/probe.py
```

The script builds its own genome from a fixed seed, so it re-runs from the file
alone; `--prefix <work>/host --k 10` points it at a real table instead. The
regressions are the six tests under "Disjointness ACROSS classes" in
`tests/test_mismatch_counts.py`.

## The defect

`mismatch_counts.mismatch_class_counts` groups a primer's binding sites by how
many mismatches they carry, and `occupancy.weighted_site_load` weights each
class at its own occupancy — class *j* at `theta(dH, Tm - j*penalty, T)`. That
arithmetic is only meaningful if the classes are disjoint, which the module
docstring claimed: "a site is never weighted at two different occupancies."

Two kinds of collision can break it, and only one was handled.

**Within a class** two different variants can share a canonical form when a
primer is near-palindromic. `_count_of` summed over distinct canonical forms,
so that case was covered.

**Across classes** the same stored k-mer can be a 1-mismatch neighbour of the
primer read one way and a 2-mismatch neighbour read the other. Nothing stopped
it. A count file holds only the lexicographic minimum of each k-mer and its
reverse complement (jellyfish `-C`), and the two readings of one duplex are not
at the same Hamming distance from the primer.

Verified rather than argued, on `ACTGTATAGT` (reverse complement
`ACTATACAGT`):

| stored canonical form | reached at distance 1 as | reached at distance 2 as |
|---|---|---|
| `ACTATAAAGT` | `ACTTTATAGT` | `ACTATAAAGT` |
| `ACTATAGAGT` | `ACTCTATAGT` | `ACTATAGAGT` |
| `ACTCTACAGT` | `ACTGTAGAGT` | `ACTCTACAGT` |
| `ACTGTAAAGT` | `ACTGTAAAGT` | `ACTTTACAGT` |

Both classes then added that site's genome count, and the occupancy model
weighted it twice — once as a 1-mismatch duplex and once as a 2-mismatch one, a
whole mismatch penalty apart.

The existing code had the shape of the fix without its generality: it dropped
variants canonicalising onto the primer itself, which is the same collision
between class 0 and class 1. The case one class further out was missed.

## Scope: only above the default

`max_mismatches` is 1 by default and the schema's minimum, and at depth 1 class
1 is the only mismatch class — there is nothing above it to collide with, and
the class 0 / class 1 pair was already handled. **So no default run is
affected, and no delivered panel moves.** `test_the_default_depth_is_unchanged`
pins that, and passes against the pre-fix code as well as the fixed one, which
is what a regression guard should do.

Depth 2 is documented and supported (`CLAUDE.md`: "distance 2 costs about
19x"), and it is where this fires. It reaches `selectivity_ratio` and
`selectivity_density` — the figures a panel is ranked on and accepted against
`min_selectivity_density` — through `base_optimizer._effective_site_load`, and
the filter's occupancy candidate ranking through the same function.

## How much it moved

300 random 10-mers at `max_mismatches=2`, against a synthetic 200 kb genome and
against the prepared 4 Mb host table:

| | 200 kb, text table | 4 Mb, KMC database |
|---|---|---|
| primers with a shared form | 4 / 300 (1.3%) | 4 / 300 (1.3%) |
| class 1 + class 2 sites, corrected | 50,233 | 997,937 |
| class 1 + class 2 sites, pre-fix | 50,242 | 998,043 |
| double-counted | 9 (**0.02%**) | 106 (**0.01%**) |
| worst single primer | 1.021x | 1.010x |

**This is a correctness fix, not a materially different answer.** The incidence
is about 1.3% of primers and the pool-level load error is one part in five to
ten thousand; the worst single primer carried 1-2% more near-match load than it
should. The value is that the classes now mean what they are documented to
mean, and a site is weighted at one occupancy — not that any recorded figure was
wrong by a visible amount.

**Both genomes are synthetic**, so k-mer counts are near-uniform; no real
reference is reachable from this container (`tests/validation/genomes/` is
gitignored and absent, and NCBI is refused by the egress proxy). A shared form
that happened to be a repeat would carry a larger count than anything here, so
**the incidence generalises and the magnitude is for this site frequency
distribution.** The second column does exercise the real KMC read path and a
4 Mb table.

## The fix

The classes are walked nearest first, carrying the canonical forms already
counted, so a site falls in the class of its **best** match — which is what one
physical duplex means. `_count_of` takes the running set, skips a form already
claimed, and extends it with the forms it counts.

Totals do not depend on set iteration order: within one call each distinct
canonical form is visited once either way, and the only ordering that decides
which class claims a shared form is the caller's loop over distances, which is
`range(1, max_mismatches + 1)`.

This subsumes the in-loop "drop variants canonicalising onto the primer"
filter, because class 0 claims the primer's own form before class 1 runs.

## Scope of the claim

One k length, two genomes, one depth, 300 primers per genome. Nothing here is a
wet-lab outcome: the loads are modelled occupancy weightings, not measured
amplification. The linear per-mismatch Tm penalty the classes are weighted
through is itself an approximation that `occupancy.mismatch_tm` states plainly,
and this fix does not touch it.
