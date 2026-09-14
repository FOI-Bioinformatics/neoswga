# Initial swap-refinement check, 2026-09-12

This is a synthetic implementation check, not evidence of experimental
amplification performance. The fixture in
`tests/test_stage2_preserves_coverage.py` uses a 100,000 bp target, 40 candidate
10-mers, random seed 4242, and 3,000 bp coverage reach. Both methods used the
same position cache and default strict dimer threshold. Swaps were limited to
10,000 evaluations and 10 seconds of search.

| Requested primers | Network covered bases | Swap covered bases | Network seconds | Swap seconds |
|---|---:|---:|---:|---:|
| 4 | 49,750 | 49,750 | 0.2362 | 0.0065 |
| 6 | 62,500 | 67,000 | 0.0345 | 0.0060 |
| 10 | 83,500 | 88,000 | 0.0133 | 0.0085 |

All panels had the requested count. Coverage here is the union of the shared
optimizer bins, weighted by their base lengths. These single-run timings were
collected in one process, with network preceding swap at each size. They
include optimizer execution and final network statistics, but exclude initial
cache construction. Import and cache warming affect these figures; they do not
establish a general speedup.

The six- and ten-primer results each gained 4,500 covered bases on this fixture.
Network refinement remains the default pending measurements across real
candidate pools, background compositions, panel sizes, and repeated timings.

Regression tests cover full-pool search, fixed primers, pairwise compatibility,
base weights, background tie-breaking, time/evaluation limits, and the hybrid
execution path. Exact-model tests compare a conflicting pair under strict and
coverage-only formulations, enforce fixed-primer budgets, and check agreement
with the independently constructed benchmark model.
