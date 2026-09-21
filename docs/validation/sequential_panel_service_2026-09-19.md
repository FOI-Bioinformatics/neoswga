# Sequential panel-stage check, 19 September 2026

An exploratory comparison on the saved Wolbachia wMel 12-mer candidate pool,
using Drosophila as background. The run used 2,000 saved candidates, a 70%
effective-coverage target for reduction and an explicit selectivity-density
floor of 20. The reaction was resolved from the example params file (phi29,
30 C, no DMSO/betaine/trehalose). These are model-based design assumptions.

The first attempt found that the older benchmark directories no longer held
indexes and that the main example indexes lacked record geometry. Results from
those attempts are not used. The successful run regenerated exact positions
for the saved candidates and their reverse complements from the FASTAs, in a
separate directory. Both target and background indexes passed the current
record-metadata guard before optimization.

Each size starts from one dominating-set proposal, shared between the comparison
arms. The single-pass arm runs shared repair/refinement; the sequential arm also
reduces the pool and revisits refinement. Both use 1,000 evaluations per swap or
repair stage and a 30-second swap limit. Reduction has its own same-sized budget.
These are equal per-stage budgets, not equal total search budgets. Times exclude
initial proposal construction/indexing, and regression tests were running during
the measurement, so they are descriptive rather than isolated performance data.

| Starting size | Workflow | Delivered | Effective coverage | Density | Background sites | Stage time (s) |
|---:|---|---:|---:|---:|---:|---:|
| 6 | single_pass | 6 | 64.17% | 599.2 | 77 | 0.68 |
| 6 | sequential_reduction | 6 | 64.17% | 599.2 | 77 | 0.69 |
| 12 | single_pass | 12 | 71.80% | 474.3 | 136 | 0.70 |
| 12 | sequential_reduction | 10 | 70.78% | 448.2 | 129 | 1.55 |
| 24 | single_pass | 24 | 77.69% | 305.3 | 264 | 0.98 |
| 24 | sequential_reduction | 9 | 70.07% | 463.7 | 111 | 4.28 |

Starting from 24 oligos, the sequential search found a 9-oligo panel retaining
70.07% effective coverage and passing the configured density and dimer limits.
Starting from 12 found a different 10-oligo panel, which illustrates search-path
dependence. Starting from 6 did not reach the 70% target. These results do not
establish that nine is the minimum, that the density floor is appropriate in the
laboratory, or that sequencing will recover 70% of the genome.

The next evidence step is an equal-total-evaluation comparison across multiple
initial proposals and coverage targets, followed by evaluation against sequencing
measurements. The staged service still exposes per-stage budgets only.

Reproduction from the repository root:

```sh
python -m scripts.benchmarking.sequential_panel_search \
  --design examples/wolbachia_pool_design \
  --output examples/wolbachia_pool_design/sequential_service_2026-09-19/results.json \
  --sizes 6 12 24 --evaluations 1000 --rebuild-indexes
```

Machine-readable results and selected sequences are in
`examples/wolbachia_pool_design/sequential_service_2026-09-19/results.json`.
