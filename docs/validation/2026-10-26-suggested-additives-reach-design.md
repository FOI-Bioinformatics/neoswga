# Do the additives the tool suggests reach a design decision?

Audited 26 October 2026, on `main` at `853c54e`. Reproduce with:

```
PYTHONPATH=. python docs/validation/suggested_additives_2026-10-26/probe.py
```

Continues [the chemistry audit](../reports/AUDIT_2026-10-26_chemistry_in_pool_design.md),
narrowing to one question: for every additive the tool suggests, is its effect
actually used in design? "Used in design" means one of two reachable paths, both
testable without a genome:

- **effective Tm moves** -> the `filter` Tm-window gate (which candidates are
  eligible) and the occupancy weighting (selectivity and effective coverage,
  the figures `optimize` scores on) both see it;
- **the design accepts it** -> `model_evidence.require_model_support` does not
  refuse the request before any index opens.

An additive that moves no Tm and is accepted silently reaches neither of the
default scoring paths. One that is recommended but refused cannot be used at all.

## Who suggests what

Three recommenders, and they do not agree with each other or with what the
design accepts:

| Recommender | Additives it outputs |
|---|---|
| `condition_suggester.suggest_conditions` | betaine, dmso, trehalose |
| `additive_optimizer.AdditiveOptimizer` | dmso, betaine, trehalose, mg; formamide and glycerol on opt-in; an **ssb** flag |
| `show-presets` | fixed preset recipes (not audited per-additive here) |

## The result

Effective-Tm delta measured at equiphi29 42 C over four oligos; design outcome
from driving each additive through `require_model_support`.

| Additive | Suggested by | dTm (C) | Design outcome | Evidence |
|---|---|---|---|---|
| dmso_percent | suggester + optimizer | 2.79 | accepted | estimated |
| betaine_m | suggester + optimizer | 1.98 | accepted | estimated |
| trehalose_m | suggester + optimizer | 1.82 | accepted | estimated |
| formamide_percent | optimizer | 1.32 | accepted | measured |
| **glycerol_percent** | **optimizer** | 0.00 | **refused** | absent |
| **ssb** | **optimizer (default-on)** | 0.00 | **accepted, never checked** | none |
| ethanol_percent | — | 1.22 | accepted | estimated |
| urea_m | — | 1.39 | accepted | estimated |
| tmac_m | — | 0.09 | accepted | estimated |
| propanediol_m | — | 5.47 | accepted | empirical |
| bsa_ug_ml | — | 0.00 | accepted | assumed |
| peg_percent | — | 0.00 | accepted | assumed |
| dtt_mm | — | 0.00 | accepted, never checked | none |

## The four core suggestions reach design

dmso, betaine, trehalose and formamide each move effective Tm by 1.3-2.8 C at a
protocol dose and are accepted. That Tm reaches two design decisions on a
default run: the `filter` Tm-window gate decides which candidates are eligible,
and occupancy weighting sets the `selectivity_density` and `effective_fg_coverage`
that `optimize` ranks on. So these four are used in design. How much they move a
delivered panel depends on the regime -- near 1x at phi29 30 C, up to 3x at
equiphi29 with k >= 10 -- which is
[additive_specificity.md](additive_specificity.md), not re-measured here. The
separate `stage2_tm_weight` selection axis is off by default
([the Tm-term record](2026-10-26-wiring-the-tm-term.md)); it is not needed for
these to reach design, because the gate and occupancy already carry the Tm.

## Two suggested additives do not reach design

### glycerol: recommended, then refused

`AdditiveOptimizer.optimize(include_glycerol=True)` adds glycerol at 5% or 10%
to its grid and can return a recommendation carrying `glycerol_percent`
(`additive_optimizer.py:313,424`). Any design that then sets glycerol is
**refused** by `require_model_support`: `tm_glycerol` is status `absent`, no
term computes glycerol's duplex effect, and the contract refuses rather than
report the same Tm with and without it. So the optimizer can recommend a recipe
the design will reject. The mechanistic model does carry glycerol enzyme-speed
terms, which is why the parameter is accepted at construction, but the design
refuses before any of that is reached.

Severity: `include_glycerol` is a public `optimize()` parameter, not wired to a
CLI flag, so the CLI `optimize-conditions` path does not currently recommend
glycerol. The contradiction is reachable programmatically and by a future flag,
not from the shipped CLI defaults. It is a latent inconsistency between what the
optimizer may propose and what the design accepts, and the two should be
reconciled: either give glycerol a duplex model (needs evidence, which the
registry says does not exist) or drop it from the optimizer's grid and
recommendation.

### ssb: recommended by default, inert on the default design

`AdditiveOptimizer.optimize` has `include_ssb=True` by default
(`additive_optimizer.py:274`), and the CLI calls `optimize` without overriding
it (`cli/commands.py:150`), so for any template with GC above 0.55 the
recommendation prints **"SSB: Yes (recommended)"** (`:153,338`). This is the one
suggested additive that fires from the shipped CLI by default.

SSB moves no effective Tm, so it reaches neither the Tm gate nor occupancy. It
is not in `_ADDITIVE_EVIDENCE`, so `require_model_support` does not even examine
it -- "accepted, never checked". Its one modeled consumer is the mechanistic
amplification term: `mechanistic_model` multiplies `kon` by `ssb_kon_multiplier`
when `conditions.ssb` (`mechanistic_model.py:613`). But that term is off by
default (`use_mechanistic_model=False`), and where enabled it is read only by
the `NetworkOptimizer` greedy selection the `network` method uses. On the
default `hybrid` method the mechanistic weight goes to the build-but-never-read
`NetworkOptimizer` (`build_unread_network_optimizer`, the same inert path as
`tm_weight`), so it reaches no delivered panel.

Net: a user who runs the default optimizer on a GC-rich target is told to add
SSB protein, and the tool's own design model then selects the identical panel
whether or not they do. Nothing in the output says the recommendation is one the
design cannot act on. SSB may well help a real reaction; the point is that this
tool does not model it on the path it designs with, so it should not present
SSB as a design recommendation without saying so.

`dtt_mm` is the same "accepted, never checked" shape (chemistry-audit gap 12)
but no recommender suggests it, so it is a smaller gap: a user only hits it by
setting it themselves.

## What this is and is not

It is a reachability audit: which suggested additives' effects are consumed by a
design decision. It is not a claim about whether any additive helps a real
reaction, and it does not re-measure the Tm coefficients (the chemistry audit
covers their evidence status). The Tm deltas are one reaction and four oligos,
enough to separate "moves Tm" from "does not"; the design outcomes are exact,
driven through the real refusal path.

## Suggested next, not done here

- Reconcile glycerol: drop it from the optimizer grid and recommendation, or
  give it a duplex model with evidence. A decision, because removing it changes
  what the optimizer proposes.
- ~~Make the SSB recommendation honest: either gate "SSB: Yes" behind a note
  that the design model does not use it on the default method, or stop
  defaulting `include_ssb` on.~~ **Done (the note).** The recommendation summary
  and rationale now state that the design model does not use SSB and that it
  changes a panel only under `--optimization-method=network` with a mechanistic
  weight; `tests/test_additive_optimizer.py` pins it. Whether to also stop
  defaulting `include_ssb` on is left open, since that changes what the
  optimizer recommends rather than only how it is explained.
- Both are instances of the chemistry audit's gap that the suggesters and the
  design accept different additive sets; this record is the per-additive
  evidence for that gap.
