# An external report on neoswga, checked claim by claim

Assessed 21 September 2026. A "Strategic Technical Evaluation and Documentation
Optimization Framework for NeoSWGA" arrived from outside this project. It reads
as competent domain writing and its code-checkable claims about this tool are
mostly wrong.

Every verdict below carries file:line evidence gathered from this repository at
the commit this document lands on. Where a claim concerns another tool it is
marked unverifiable here, because an in-repo statement about another project is
not verification either.

## Two recommendations that must not be followed

### "Explicitly set `min_amp_pred` to 10 or higher to eliminate low-efficiency primers"

**False, and silently so.** The amplification gate was retired from the default
path on 2026-09-05. `min_amp_pred` defaults to 10.0
(`core/parameter.py:352`), and the warning that the gate does nothing only
fires when the configured value DIFFERS from that default
(`core/pipeline.py:1402-1403`). So a user following this advice literally sets
the value to its default, gets no warning, and believes a filter ran that did
not.

The advice would be poor even with the gate enabled. It removed 7 of 1222
candidates on the S. aureus panel and none at all on E. coli (0 of 449) or
M. tuberculosis (0 of 319), because the scores cluster well above 10.0. And
every step-4 consumer reads only the primer column.

`--amp-model` restores the score and its gate for anyone who wants it.

### "Phi29 processivity is ~70 kbp. The optimization logic must place primers within this fragment range"

**The number is right and the prescription is wrong.** 70 kb is real
(`core/registry/polymerases.py:152`) and is the reach used for
amplification-network CONNECTIVITY. Selection and scoring use the realistic
per-primer reach of about 3 kb
(`core/coverage.py:367,374`, `coverage_metric='realistic'`).

Designing to 70 kb inflates reported foreground coverage by roughly 5 to 20
times, because a coverage window is twice the reach wide. A user who followed
this would read a coverage figure describing a genome their panel does not
cover. The two reaches are documented as distinct in CLAUDE.md for exactly this
reason.

## Wrong about the tool

| Claim | Verdict and evidence |
|---|---|
| "Restricted to fixed k-mer lengths (6-12 bp)" | **False on both halves.** `core/schema/params.schema.json` sets `min_k`/`max_k` to integers from 4 to 30. 6-12 is phi29's default window; bst is 15-25 and equiphi29 10-18. A mixed k 7-11 design is now demonstrated end to end in `tests/integration/test_variable_oligo_length.py`, and measured in [variable_oligo_length_2026-09-21.md](variable_oligo_length_2026-09-21.md). |
| "RF trained via Active Learning with iterative RCA experiments" | **False.** The bundled model is fitted to synthetic data whose labels come from a hand-written rule with Gaussian noise; `core/rf_preprocessing.py:33` states "No measured amplification outcome enters the fit." The RCA claim belongs to swga 2.0, as this repository's own `tool_comparison.md:237` records. |
| "A decisive departure from traditional greedy search" | **False.** The default `hybrid` method's Stage 1 IS a greedy set cover with an ln(n) approximation bound (`core/dominating_set_optimizer.py:10`). ILP exists only as a diagnostic. |
| "Target a pool of ~100 candidates for Stage 4" | **False.** `max_primer` defaults to 500 (`core/parameter.py`), and `candidate_retention` defaults to `all_qc`, which keeps every candidate clearing the hard gates addressable. |
| "Set `max_bg_freq` as low as possible, but not below 3 x 10^6" | **False and internally incoherent.** It is a dimensionless frequency bounded to 0-1 by the schema, default 5e-6. 3e6 is out of range by six orders of magnitude. Reading it as the plausible 3e-6 does not rescue the advice: that is STRICTER than the default, so "as low as possible but not below" contradicts itself. |
| "Binding propensity through dG_T, a smoother metric than exact-match counts" | **False.** No free-energy term appears in `base_optimizer`, `occupancy` or `dominating_set_optimizer`. Selectivity is exact-match counts, occupancy-weighted from enthalpy and effective Tm. Known Issue 17 separately measured that occupancy barely discriminates at phi29 30 C. |
| "Off-target binding to host DNA, the primary driver of sWGA failure" | **Contradicted by this repository's own evidence.** `published_primer_sets.md` finds mean binding distance does not predict success, that the two wet-lab datasets disagree about which statistic does, and that site ratio does not predict selectivity. |
| Stage 3 described as "Efficacy Prediction" | **Stale.** The stage is `prepare-candidates` since 2026-09-21 and predicts nothing by default. |

## Correct

The 3' GC clamp exists and is applied (`core/filter.py:92-94`, `:425`). It is
a band with a lower bound rather than a simple cap, it adapts to genome GC
content, and it is joined by a stricter rule on the last three bases. The
report describes the concept correctly.

Python 3.11 or later, and the `init`, `start` and `validate --quick` commands,
all exist as described. Jellyfish counting is lock-free.

## Unverifiable here

swga 2.0's reported 58-minute runtime, COATswga's KMC3 engine and its 98.9%
*Plasmodium falciparum* mapping figure, and the claim that interval tiling
"mimics phi29 more closely". None is refuted. None is checkable from this
repository, and none should be repeated as fact. The same caution applies to
this project's own `tool_comparison.md`, whose statements about other tools
were read from their source and not re-verified.

## What the report got right in spirit

Its criticism of variable oligo length pointed at something real even though it
was wrong about the tool. Nothing here required one k, but nothing here had
DEMONSTRATED a mixed-length design either, and every saved run used a single
length. That gap is now closed by a test and a measurement, and the measurement
found that mixing bought nothing on the one pair tried, which is an answer the
project did not previously have.

Its point about documentation also stands. No document joins the biology in
`SWGA_SCIENCE.md` to the four pipeline steps; `QUICK_START.md` gives commands
without rationale and `user-guide.md` is mechanical. That gap is real and is
not addressed by this document.

## Method, and its limits

Each claim was checked by reading the named code and, where behaviour was at
issue, by running it. The verdicts are about what this repository does today.
They are not a judgement of the report's domain reasoning, much of which
concerns wet-lab practice this project cannot evaluate.
