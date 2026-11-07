# One chemistry, two answers: `--preset` accepted the glycerol a design refuses

Found and fixed 1 November 2026, on `main` at `41d0eb5` plus the change this
document records, in the chemistry and design-flow audit that also produced
[`2026-11-01-mismatch-classes-across-strands.md`](2026-11-01-mismatch-classes-across-strands.md).

Reproduce the before-and-after with:

```
PYTHONPATH=. python - <<'PY'
from neoswga.cli._common import load_preset_conditions
from neoswga.core.design_request import resolve_design_request
base = {"fg_prefixes": ["t"], "fg_seq_lengths": [1270000], "min_k": 10, "max_k": 12}
for preset in ("q_solution", "crude_sample", "gc_melt"):
    params = dict(base); params.update(load_preset_conditions(preset))
    resolve_design_request(params)   # refused for the first two before the fix
PY
```

The regressions are `test_an_additive_a_design_refuses_is_absent_rather_than_zeroed`
and `test_at_least_one_preset_exercises_the_drop` in
`tests/test_preset_round_trip.py`.

## The two doors

`model_evidence.require_model_support` refuses an additive whose duplex effect
nothing computes — evidence status `absent`, which today is glycerol alone.
`CLAUDE.md` states it plainly: "A design that sets glycerol is refused: nothing
computes its Tm effect." That is true of `params.json`.

It was not true of `--preset`. `resolve_design_request` reads the params FILE,
deliberately: `design_request_for_run` documents that the `parameter` module is
"one resolution of it" and that reading the module would make an identical
request resolve differently depending on what ran before. A preset reaches the
`parameter` globals and never the file, so its additives never passed that
refusal.

Two of the eleven `--preset` choices on `neoswga filter` set glycerol:

| preset | glycerol | also sets |
|---|---|---|
| `q_solution` | 10.0% | betaine 1.5 M, BSA 200 µg/mL |
| `crude_sample` | 10.0% | betaine 1.0 M, BSA 400 µg/mL, PEG 5% |

So `neoswga filter --preset crude_sample` succeeded, and the value was
accepted, range-checked against 0–15%, and **written into the run manifest's
`effective_conditions`**:

```
manifest glycerol_percent: 10.0
```

while contributing exactly nothing to any Tm — no term computes it, which is
the whole reason a design refuses it. The record claimed a reaction the
computation did not use. The same number written into params.json was refused
outright.

Both halves are defensible on their own. Reading the file is what makes a
request a function of its input; refusing an unmodelled additive is what stops
a fabricated Tm. Together they left one chemistry with two answers depending on
which door it came through, and the silent door was the one that wrote the
record.

## The fix

`model_evidence.refuses_additive(field)` is the predicate, public and in one
place, because the question is now asked from three: the refusal itself, the
additive search's grid (which already had its own copy as
`AdditiveOptimizer._additive_is_usable`, now delegating), and the preset
loader.

`load_preset_conditions` leaves such a field **out** of the flattened dict and
names it:

```
WARNING: Preset 'crude_sample' sets glycerol_percent=10.0, which is left out:
no term models its melting-temperature effect (evidence status 'absent'), and a
design refuses it. Every other setting in the preset is applied.
```

After, on the same run:

```
manifest glycerol_percent: 0.0
manifest bsa_ug_ml      : 400.0
manifest peg_percent    : 5.0
manifest betaine_m      : 1.0
```

**Dropped, not zeroed.** Omitting the key leaves `parameter.glycerol_percent`
at its module default, so the preset simply does not speak about glycerol.
Writing 0.0 into the dict would be a different claim — "this reaction has no
glycerol" rather than "this preset has no opinion" — and would stop an explicit
params.json value being overridden by nothing.

**Dropped, not refused.** The rest of both presets is sound, and `crude_sample`
is the configuration the tool recommends for crude lysates, blood and soil.
Refusing the whole preset would withdraw a recommended setup over one
unmodelled component, which is a larger change than the evidence asks for.

Presets that set no refused additive are untouched: the drop is guarded on a
non-zero value, so `glycerol_percent: 0.0` still travels for the other nine.

## What this does not change

- **No delivered panel moves.** Glycerol entered no Tm before the fix and
  enters none after; the only difference is that the run no longer records a
  component it did not use, and the preset no longer disagrees with the file.
- **No threshold or coefficient moves.** Glycerol's evidence status is
  unchanged, and nothing here models it. `mechanistic_model` does read
  `glycerol_percent` for polymerase stability and speed, which is a different
  quantity from duplex stability and is not what the registry's `absent`
  status is about; that path is untouched.
- **It is not a general audit of the preset path.** One registry status on one
  additive, across eleven presets. BSA and PEG carry status `assumed` with no
  model — a deliberate absence, since neither acts on duplex stability — and
  are correctly still applied.
