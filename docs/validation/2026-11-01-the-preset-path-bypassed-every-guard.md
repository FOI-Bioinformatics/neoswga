# `--preset` bypassed every guard, and overrode the recommendations it logged

Found and fixed 1 November 2026, auditing the rest of the preset path after
[`2026-11-01-a-preset-carried-an-additive-a-design-refuses.md`](2026-11-01-a-preset-carried-an-additive-a-design-refuses.md)
closed one instance of it. On `main` at `85ddbe5` plus the change this document
records.

Reproduce the guard matrix with:

```
PYTHONPATH=. python - <<'PY'
from neoswga.cli._common import load_preset_conditions
from neoswga.core.design_request import resolve_design_request
base = {"fg_prefixes": ["t"], "fg_seq_lengths": [1270000], "min_k": 6, "max_k": 12}
for name in ("standard_phi29", "bst", "klenow", "enhanced_equiphi29"):
    merged = dict(base); merged.update(load_preset_conditions(name))
    try:
        resolve_design_request(merged); print(name, "accepted")
    except Exception as exc: print(name, "REFUSED", exc)
PY
```

The regressions are the two sections added at the end of
`tests/test_preset_round_trip.py` (17 tests).

## One cause, two consequences

A preset reached the `parameter` module globals by `setattr` and nothing else.
It never entered the params file, and it never entered
`parameter._json_data`. Both omissions matter, because two different
mechanisms ask two different questions and neither could see a preset:

- `resolve_design_request` reads the params FILE, deliberately, so a request
  is a function of its input. **So none of the four guards in
  `require_model_support` ever saw a preset's chemistry.**
- `_apply_adaptive_additives` asks whether a key is PRESENT in `_json_data`,
  deliberately, to tell a user's explicit `0.0` from a key they never wrote.
  **So a preset's explicit `0.0` read as silence.**

## Consequence 1: the oligo-length guard, live on eight of eleven presets

A preset sets the **polymerase**; `min_k` and `max_k` come from the file.
Against the documented 6-12 window:

| preset | enzyme | enzyme's range | refused for |
|---|---|---|---|
| `enhanced_equiphi29`, `high_gc_genome`, `long_primers_15mer`, `q_solution`, `gc_melt`, `extreme_gc` | equiphi29 | 10-18 nt | k 6, 7, 8, 9 |
| `bst` | bst | 15-25 nt | k 6-12 |
| `klenow` | klenow | 8-15 nt | k 6, 7 |
| `standard_phi29`, `crude_sample`, `low_temp` | phi29 | 6-12 nt | — accepted |

Written into params.json, those eight are refused. Via `--preset` all eight
were accepted, and the computation ran. Measured on the plasmid example, a
fresh data directory per arm, with the no-preset arm confirmed deterministic
across two runs:

| arm | delivered candidates | oligo lengths |
|---|---|---|
| no preset | 37 | k 7-11 |
| `--preset bst` | **38** | k 7-**12** |
| `--preset klenow` | **34** | — |

So `--preset bst` computed Tm for 7-12mers at 63 C from a parameter set
covering 15-25 nt, and the shortlist changed. That is the defect
`_require_length_support` was written for, in its own words: *"A Bst design
filtered through phi29's 6-12 bp window is a recorded defect in this project.
The Tm computed for a 6-mer under Bst at 63 C is not wrong in a way anyone can
see: it is a number with no evidence behind it."*

The protection that would normally catch it is the GC-adaptive clamp into the
enzyme's own range -- and it is skipped exactly when the user sets
`min_k`/`max_k`, which the documented workflow does.

## Consequence 2: a preset overrode GC-adaptive, and the log said otherwise

`--preset standard_phi29` ought to be a no-op on a params.json already naming
phi29 at 30 C with the same salts. It was not. `effective_conditions` differed
in exactly one field of twenty-one -- `betaine_m`, 1.0 against 0.0 -- and the
run delivered 38 candidates against 37.

The preset writes fifteen explicit zeros, including every additive. Before the
fix, on one run:

```
INFO: GC-adaptive: Setting betaine to 1.0M     <- logged
INFO: Loading preset: standard_phi29           <- applied after, silently
effective_conditions betaine_m: 0.0            <- what ran
```

`_apply_adaptive_additives` carries a dated rationale for keying on presence:
*"A params.json excluding betaine was therefore run at the recommended
concentration, and every melting temperature in that design was computed
against a buffer the user had ruled out."* A preset's zero is the same
statement by the user, and the guard could not read it.

## The fix

**A preset is a user choice, so it joins the configuration where user choices
live.** `parameter._preset_overlay` is set by the CLI before `_initialize()`
and merged into `data_extra` by `get_params`, so the preset's keys land in
`_json_data` alongside the file's, winning over them -- which is the precedence
the CLI already applied by overwriting the globals afterwards. CLI flags are
applied later still and keep winning over both.

Every key-presence guard then works for presets at no further cost. GC-adaptive
now logs what the run uses:

```
INFO: GC-adaptive: Preserving user-specified betaine 0.0M
INFO: GC-adaptive: Preserving user-specified DMSO 0.0%
```

**No delivered value moves.** Before the fix GC-adaptive set a value and the
preset overwrote it; after, GC-adaptive preserves the preset's value. Same
final chemistry, so the four arms deliver what they delivered before -- 37, 38,
38, 34 for no preset, `standard_phi29`, `bst` and `klenow`. Only the record
stops contradicting the run.

**The guards are applied to the resolved configuration.**
`_common.validate_effective_chemistry` builds the effective mapping from
`_json_data` plus `_effective_conditions(parameter)` plus the resolved
`min_k`/`max_k`, and runs it through `resolve_design_request`. Called from
`run_step2` after the polymerase is settled and before anything is counted.

It runs **only when a preset was applied**, which bounds the change: without
one the effective configuration IS the file, which `optimize` already resolves,
so no run that works today starts failing. Extending it to every `filter` run
would be a further step, and a defensible one -- the project's own
`effective_design_request` argues "validation belongs before the position cache
is built; a typo should not cost a full index" -- but it would widen the
refusal surface beyond the gap this record is about.

The remedy the message names works, verified: `bst` with `min_k`/`max_k` of
15-25 or 15-18 is accepted, `equiphi29` at 10-18, `klenow` at 8-15.

## What this does not do

- **It does not close the class.** Known Issue 8 is not to be recorded as
  closed, and this is a new shape within it: a GUARD missing a door, rather
  than an option missing a reader. Two of the four guards remain structurally
  bypassable with no live instance -- every preset names a registered
  polymerase, and DMSO is the only additive carrying a temperature domain
  (30-45 C) with all seven DMSO-bearing presets inside it, `gc_melt` at 45 C on
  the boundary. Either would become live from one registry edit.
- **It does not change what any preset contains.** `long_primers_15mer` still
  sets no `min_k`/`max_k` despite its name, so it is refused against the
  default window rather than fixed. Making presets set the oligo window was
  considered and rejected: it would silently rewrite the user's filtering
  range.
- **It is not a measurement of chemistry.** The 37/38/34 figures are one 6 kb
  plasmid; they establish that the geometry changes the delivered shortlist,
  not by how much on a real genome.
