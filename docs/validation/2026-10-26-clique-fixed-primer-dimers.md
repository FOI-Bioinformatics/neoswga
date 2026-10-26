# The clique method broke its own dimer-free guarantee with fixed primers

Found and fixed 26 October 2026, on `main` at `ecb3553`, during an optimizer
audit. Reproduce the proof with the probe inline below; the regression is
`tests/test_clique_fixed_primers_and_defaults.py`.

## The guarantee, and the hole

`clique_optimizer` exists for one reason: it is "the only method that guarantees
no dimerising pair" (CLAUDE.md). Its module docstring promises "every dimer-free
primer set of the target size is enumerable". With fixed primers supplied, it
delivered sets that contained a dimerising pair, silently.

Two independent faults, both in `optimize`'s fixed-primer handling.

### 1. The heterodimer filter was dead

`optimize` removes the fixed primers from the candidate pool before building the
compatibility graph:

```python
candidates = [p for p in candidates if p not in fixed_set]   # fixed removed
...
G = build_compatibility_graph(candidates, ...)               # graph has no fixed primers
```

Then it tried to keep only candidates compatible with each fixed primer by
looking up the fixed primer's neighbours in that graph:

```python
for fp in fixed:
    if fp in G:                       # ALWAYS False -- fp was removed above
        compatible_nodes &= set(G.neighbors(fp))
    else:
        logger.warning(... "possible self-dimer")   # the branch always taken
```

Because the fixed primers are never nodes in `G`, `fp in G` is always False, the
`else` "proceed anyway" branch always runs, and `compatible_nodes` is never
intersected. A candidate that heterodimers with a fixed primer survives the
filter and is delivered in `fixed + primer_set`.

This is not the self-dimer edge case the audit first guessed; it fires for
**every** fixed primer. Proved at the exact filter step:

```python
from neoswga.core.clique_optimizer import build_compatibility_graph
from neoswga.core.dimer import is_dimer_fast
F, partner = "AAAAGGGGAAAA", "TTTTCCCCTTTT"          # heterodimer; F does not self-dimer
benign = ["ACACACATACAC","ATACACAACAAT","AATACAATACAA","ACAATACAAATA","TACAAATAACAA"]
candidates = [p for p in benign + [partner] if p != F]   # optimize() removes F
G = build_compatibility_graph(candidates, max_dimer_bp=3, max_self_dimer_bp=4,
                              check_self_dimers=True)
compatible = set(G.nodes())
for fp in [F]:
    if fp in G:                                      # False
        compatible &= set(G.neighbors(fp))
# partner survives, though is_dimer_fast(F, partner, 3) is True
assert partner in compatible and is_dimer_fast(F, partner, 3)
```

### 2. A self-dimering fixed primer was delivered, not refused

The method drops a self-dimering *candidate* from the graph and refuses a
mutually dimering *pair* of fixed primers, but a self-dimering *fixed* primer hit
the same "proceed anyway" branch and was returned. A self-dimer is a dimerising
pair (the primer with itself), so a set containing it is not dimer-free. The
guarantee cannot hold, so the method must refuse, not deliver.

## The fix

Both in `optimize`'s fixed-primer handling.

- **Refuse what cannot be dimer-free.** A new `_reject_undimerfree_fixed` holds
  fixed primers to the same standard as the rest: a self-dimering one (gated on
  `check_self_dimers`, as candidates are) or a mutually dimering pair is
  refused with a clear message, the same answer the existing mutual-pair check
  already gave.
- **Screen candidates against fixed primers directly**, not through the graph:

  ```python
  compatible_nodes = [
      node for node in G.nodes()
      if not any(is_dimer_fast(node, fp, max_dimer_bp) for fp in fixed)
  ]
  ```

  This enforces fixed-vs-candidate heterodimer compatibility regardless of
  whether the fixed primer is a graph node, which it never is.

After the fix, the probe's `partner` is excluded, and a self-dimering fixed
primer returns `OptimizationStatus.ERROR` with "forms a self-dimer".

## Scope and cost

Narrow in reach, central in meaning. `clique` is opt-in
(`--optimization-method clique`) and the fault needs fixed primers, which come
from iterative design (`improve-set`, `expand-primers`) or an explicit list. But
within that path it defeated the one property the method is chosen for, so the
severity is the guarantee, not the frequency.

No behaviour changes for a clique run without fixed primers: the graph and its
cliques are unchanged, and `_reject_undimerfree_fixed` returns None. Only the
fixed-primer path changes, and only to honour the guarantee the method already
claimed. `optimize` grew past its function-length budget with the added checks,
so the validation was extracted to `_reject_undimerfree_fixed`.

Pinned by `test_a_candidate_that_dimers_with_a_fixed_primer_is_never_delivered`
and `test_a_self_dimering_fixed_primer_is_refused`.

## Not fixed here (optimizer audit, same day)

The audit that found this surfaced other items, left for a decision:

- `objective_scan_width` (default 64) is read only on the `plan-pool` repair
  path; an `optimize` run's refinement never sees it, bounded only by
  `swap_max_evaluations`. Documented without that scoping. Wiring it into
  `optimize` needs coverage bins and a measurement; the honest alternative is to
  scope its documentation.
- `minimal_primer_selector.py` is dead (the real `--minimize-primers`
  post-process is `optimization_service.reduce_result`); CLAUDE.md and
  implementation-notes describe it as the post-process. A doc correction plus a
  decision on removing the module.
- Known Issue 14 (neither Stage-2 path carries both the shared objective as its
  primary ranking and the host term as a weighted axis) remains a deliberate,
  unmeasured modelling choice, not a defect.
