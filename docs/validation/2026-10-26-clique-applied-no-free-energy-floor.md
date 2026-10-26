# `clique` applied no free-energy floor

Found 26 October 2026 while asking the "does it reach the innermost decider?"
question of the methods the scoring-engine ratchet does not cover (recorded as
left-open in
[`2026-10-26-dominating-set-dimer-temperature.md`](2026-10-26-dominating-set-dimer-temperature.md)),
and fixed the same day on `main` at `91ee82d`. Reproduce with:

```
PYTHONPATH=. python docs/validation/clique_dimer_dg_2026-10-26/probe.py
```

The regressions are in `tests/test_clique_fixed_primers_and_defaults.py`
(`test_the_fixture_isolates_the_free_energy_floor`,
`test_a_candidate_that_only_fails_the_free_energy_floor_is_never_delivered`,
`test_without_a_floor_the_same_candidate_is_delivered`,
`test_fixed_primers_that_only_fail_the_free_energy_floor_are_refused`,
`test_the_floor_only_ever_removes_edges`).

## The defect

`clique` is documented as "the only method that guarantees no dimerising pair".
It read `max_dimer_dg` **nowhere**: zero occurrences in `clique_optimizer.py`,
and no such field on `CliqueOptimizerConfig`. All three of its heterodimer
questions went through `is_dimer_fast` or `heterodimer_matrix_fast`, which count
complementary base pairs and know nothing about duplex stability:

| site | what it screened |
|---|---|
| `build_compatibility_graph` edges | base pairs only |
| `_reject_undimerfree_fixed` fixed pairs | base pairs only |
| `optimize`'s fixed-compatibility filter | base pairs only |

So a user who configured a free-energy floor and chose the method sold on a hard
dimer guarantee got base-pair screening alone, silently. This is Known Issue 8
in the shape where an option is declared on the base `OptimizerConfig`, honoured
by other methods, and one method never looks at it.

## That it mattered

The fixture is chosen so the floor is the **only** discriminator: across the
pool no pair exceeds `max_dimer_bp=3` and none self-dimers at 4, and exactly one
pair falls below a -1.0 kcal/mol floor at 37 C. Anything the floor excludes,
only the floor could have excluded.

With that candidate made the most attractive in the pool, so only a screen can
keep it out:

| `max_dimer_dg` | offending candidate delivered? |
|---|---|
| unset (default) | **True** -- the control: it clears every screen clique applies |
| -1.0 | **False** |

Before the fix the second row also read True. The delivered panel contained a
pair below the configured stability floor, from the method whose entire purpose
is to forbid dimerising pairs.

## The fix

The floor is consulted only for pairs the base-pair screen has **already
cleared**, which is what makes it purely additive -- and it is routed through
`lazy_dimer.dimer_screen`, the project's one door for dimer screening, so the
two tests compose here in the same order as everywhere else
(`LazyDimerCompatibility._pair` applies `is_dimer_fast` first and the floor only
to survivors). A floor forces the pairwise screen, whose object is
pool-independent.

One pair test, `CliqueOptimizer._heterodimer_test`, now answers every
heterodimer question in a run, so the fixed primers and the candidates are held
to the same screen. `_screen_inputs` resolves `(max_dimer_dg, temp)`.

**The temperature is the reaction's, not 37 C.** That is the lesson of the two
fixes this one followed: delta-G is strongly temperature-dependent, and a floor
evaluated at the wrong temperature is
[`2026-10-26-dominating-set-dimer-temperature.md`](2026-10-26-dominating-set-dimer-temperature.md).
It is read off the attached `ReactionConditions`. The same floor over 60 random
14-mers at `max_dimer_bp=7`:

| | 30 C | 65 C |
|---|---|---|
| floor -3.0 | -242 edges | -5 edges |
| floor -6.0 | -11 edges | 0 edges |

from 1,768 unscreened edges. The screen genuinely follows the reaction.

The self-dimer checks stay on base pairs: `max_self_dimer_bp` is its own knob and
`lazy_dimer` has no self-dimer free-energy notion to route through.

## What does NOT change

- **Nothing moves by default.** `max_dimer_dg` ships unset, and with no floor
  `build_compatibility_graph` builds exactly the graph it built before, from the
  same base-pair matrix. `test_without_a_floor_the_same_candidate_is_delivered`
  is the control that pins this: without that test, the behavioural test above
  could pass for an unrelated reason and the change would be moving the default
  rather than honouring an opted-in setting.
- **The floor can only tighten.** Pinned over two floors and two temperatures by
  `test_the_floor_only_ever_removes_edges`, which the fix keeps true by
  construction -- the floor never sees a pair the base-pair screen rejected.
- `build_compatibility_graph`'s new parameters default to `None` and 37.0, so a
  direct caller of the helper is unchanged, as is
  `_reject_undimerfree_fixed`'s new `dimerises` argument.

Clique's guarantee is now strictly stronger: no pair dimerising by base-pair
count **or**, when a floor is configured, by duplex stability at the reaction
temperature.

## Scope of the claim

Synthetic primers throughout; no delivered panel from a real design was
re-measured, because generating a candidate pool needs a k-mer counter and none
is installed in this container. What is established is that the option now
reaches all three screening sites, changes the delivered set in the direction it
claims, follows the reaction temperature, and only tightens. How many primers a
real clique design loses to a given floor is not determined here.

Full suite before and after: the same 60 failures, all environmental here (no
k-mer counter, no `pysam`). Zero added, zero removed.
