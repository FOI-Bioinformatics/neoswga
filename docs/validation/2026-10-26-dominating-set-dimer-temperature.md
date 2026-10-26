# The dominating-set free-energy dimer screen ran at 37 C, not at the reaction

Found and fixed 26 October 2026, on `main` at `40af6e0`, by asking the
"does it reach the innermost decider?" question of the two methods the
scoring-engine ratchet does not cover. Reproduce with:

```
PYTHONPATH=. python docs/validation/dominating_set_dimer_temp_2026-10-26/probe.py
```

The regressions are `tests/test_optimizer_config_reaches_optimizers.py`
(`test_the_reaction_temperature_reaches_the_dominating_set_dimer_screen`,
`test_no_conditions_leaves_the_dominating_set_screen_temperature_alone`).

## The defect

`DominatingSetAdapter` forwards `max_dimer_dg` -- the free-energy floor on
duplex stability -- to the inner `DominatingSetOptimizer`, and did not forward
`dimer_temp`, the temperature that floor is evaluated **at**. The inner kept its
37 C constructor default on every run.

```python
# dominating_set_optimizer.py, the decision site
matrix = dimer_screen(
    pool, self.max_dimer_bp, max_dimer_dg=self.max_dimer_dg, temp=self.dimer_temp,
)
```

Both ends of the option existed and the path between them did not: Known Issue
8's usual shape, and the *second* instance of it in this same constructor call.
The comment one line above the fix records the first, for `max_dimer_bp`.

It is also an inter-method inconsistency. Of the four `dimer_screen` call sites,
three pass the reaction temperature and only this one did not:

| call site | temperature |
|---|---|
| `network_optimizer.py:872` | `reaction_temp` |
| `swap_refinement.py:520` | `reaction_temp` |
| `panel_refinement.py:80` | `conditions.temp` |
| `dominating_set_optimizer.py:489` | `self.dimer_temp`, never filled -> **37 C** |

## Why it matters: the error is in the permissive direction

Free energy is strongly temperature-dependent, so the screen's verdict changes
with the temperature it is evaluated at. Over 30,000 random 14-mer pairs,
keeping only those that PASS the base-pair screen at `max_dimer_bp=7` (those are
the only pairs the floor is ever consulted for), the verdict differs between the
phi29 default 30 C and the hardcoded 37 C:

| floor (kcal/mol) | pairs whose 30 C and 37 C verdicts differ |
|---|---|
| -1.0 | 2,416 |
| -2.0 | 1,335 |
| -3.0 | 1,316 |
| -4.0 | 723 |
| -5.0 | 409 |
| **-6.0** (Rychlik 1995 guideline, the documented default) | **171** |
| -7.0 | 59 |
| -8.0 | 29 |

**Every one of those flips is "dimer at 30 C, not a dimer at 37 C."** That is
the direction that matters: screening at 37 C under a 30 C reaction *under-calls*
dimers, so a pair that would dimerise in the actual reaction passes the screen
and can reach the delivered panel. This project has shipped an 11 bp heterodimer
against a configured 3 once already
([`2026-10-26-clique-fixed-primer-dimers.md`](2026-10-26-clique-fixed-primer-dimers.md)),
and under-calling is how that happens.

For a hot polymerase the error runs the other way and merely over-rejects: `bst`
at 60-65 C screened at 37 C discards pairs that would not survive the reaction,
shrinking the pool without saying so.

## The fix

The adapter reads the temperature off the attached `ReactionConditions`, the way
`panel_refinement` does, and forwards it:

```python
dimer_temp=float(getattr(conditions, "temp", 37.0) or 37.0),
```

## What does NOT change

- **Nothing moves on a default run.** `max_dimer_dg` is unset by default, and
  with no floor the screen never consults a temperature at all. The fix is only
  observable for a user who opted into a free-energy floor.
- **No temperature is invented.** With no conditions attached, the inner's own
  37 C default stands rather than a polymerase preset being substituted --
  pinned by `test_no_conditions_leaves_the_dominating_set_screen_temperature_alone`.
  This is the "unknown is not zero" rule applied to a default: the fix routes a
  measured temperature when there is one and fabricates nothing when there is not.

## Scope of the claim, and what this audit did not fix

The envelope above is a count of random pairs, not a delivered-panel difference:
it establishes that the temperature changes the screen's verdict, and in which
direction, not how many primers a real design gains or loses. No end-to-end
design was re-run, because the defect is unreachable without `max_dimer_dg`,
which no shipped preset sets.

The same question asked of `clique` surfaced a **different** gap, left unfixed
and deliberately not bundled here: `clique_optimizer` reads no `max_dimer_dg` at
all (zero occurrences in the module, and `CliqueOptimizerConfig` has no such
field). It screens on base-pair count only, via `is_dimer_fast` and
`heterodimer_matrix_fast`. So a user who sets a free-energy floor and selects
the method advertised as "the only method that guarantees no dimerising pair"
gets no free-energy screening. Honouring it there means giving clique a second
screen path, which changes what its enumeration guarantee covers, so it is a
design decision rather than a forwarding repair.
