# What making lint blocking would actually catch

Measured 22 September 2026, on `neoswga/` at the commit this document lands on.
Written because "make lint blocking" is a plan item that sounds obviously good
and is mostly not, and the difference matters for how much work to spend on it.

## The headline number is misleading

| scope | findings |
|---|---|
| `neoswga/` and `tests/` | 3,297 |
| `neoswga/` alone | 2,696 |
| of those, auto-fixable | 2,221 |

Almost all of it is one thing. Typing modernisation, `typing.List` to `list`
and `Optional[X]` to `X | None`, accounts for about 2,180:

| rule | count | what it is |
|---|---|---|
| UP006 | 1,352 | non-PEP585 annotation |
| UP045 | 511 | non-PEP604 optional |
| UP035 | 314 | deprecated import |
| UP007, UP015, UP037 | 34 | the same family |

None of it can change behaviour. Blocking on it would mean a 2,200-line
mechanical diff for no defect caught, and would bury the rules that do catch
defects underneath it.

## What the defect-prone rules actually found

| rule | count | verdict |
|---|---|---|
| B905 `zip()` without `strict` | 64 | **one real defect**, see below |
| F401 unused import | 187 | none; mostly re-exports |
| F811 redefinition | 7 | none; all redundant re-imports of `json`, `hashlib`, `time` |
| B023 loop-variable closure | 5 | none, but see below |
| F841 unused variable | 20 | none |
| B007 unused loop variable | 20 | none |
| E711, E712, F632, B006, B008, B017 | 0 | already clean |

### The one real defect

`bam_coverage` paired `fg_prefixes` with `fg_seq_lengths` using `zip`, which
truncates. A params file naming three targets and two lengths analysed two of
them, silently, and every gap in the third was then reported as absent. Fixed
separately.

That is one defect from 2,696 findings. It is also a consequential one, in the
family this project spends most of its effort on, and it was not caught by any
test.

### The five that look like defects and are not

`optimization_service._run_panel_stages` defines five closures inside a loop
that capture `result`, `previous` and `diagnostics`. That is the classic
late-binding shape, and `execute` calls `operation()` immediately, so every one
of them is correct today.

Worth recording rather than dismissing: if `execute` ever became deferred, for
instance to queue stage work, all five would silently use the wrong value, and
nothing in the suite would notice. The rule is a real guard against a change
nobody would think to check.

## What was done

Eleven rules that catch wrong behaviour and are **already at zero** are now
blocking in CI: `E711`, `E712`, `E713`, `E714`, `F632`, `F702`, `B006`, `B008`,
`B017`, `S102`, `S307`. That costs nothing today and refuses the next instance.
Verified load-bearing by introducing three violations and watching the check
fail.

The full run stays informational, because of the 2,180 cosmetic findings above.

## What is left, and why it was not done here

**B905 is the rule that earned its keep and is not yet blocking.** Sixty-four
sites need a judgement each: `strict=True` where the two sequences must
correspond, `strict=False` where one is derived from the other and a mismatch
is impossible by construction. Both are better than silence, and which applies
cannot be decided in bulk.

The 2,180 typing findings should be fixed eventually, in one labelled
mechanical commit with no other change in it, so the diff can be skimmed rather
than reviewed. They are not urgent and they are not defects.

## Caveat

This measures what lint finds TODAY, not what it prevents. A guard that catches
nothing on a clean codebase is still worth having; the argument for blocking is
about the next defect, not the current backlog. What this document rules out is
the idea that 2,696 findings represent 2,696 hidden problems.
