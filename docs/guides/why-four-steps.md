# Why the pipeline has four steps

This document explains what each pipeline step is *for*. It sits between
[SWGA_SCIENCE.md](../SWGA_SCIENCE.md), which explains the biology and
chemistry without reference to this tool, and
[QUICK_START.md](QUICK_START.md), which gives the commands without the
reasoning. Read it when you want to know why a step exists, what decision it
makes on your behalf, and which parameter changes that decision.

It is deliberately honest about limits. Several things this tool reports are
measurements it makes but does not act on, and a few are not modelled at all.
Those are marked, because a design guide that only lists capabilities is how
people come to believe a number means more than it does.

## The problem being solved

Selective whole-genome amplification amplifies a target genome out of a sample
dominated by something else, usually a host. A strand-displacing polymerase
extends from wherever a short primer happens to bind. If your primers bind the
target far more often than the host, the target is amplified preferentially.

So a design is a *set* of short oligos with three properties, and the four
steps exist because these are three different kinds of problem:

1. **Each oligo must be commoner in the target than in the host.** That is a
   counting problem over the two genomes.
2. **Each oligo must actually bind at the reaction temperature, and bind in
   enough places, spread out.** That is a per-oligo thermodynamics and geometry
   problem.
3. **The set together must cover the target without its members sticking to
   each other.** That is a combinatorial problem over sets, and it cannot be
   solved by picking the individually best oligos.

Steps 1 and 2 handle the first two. Step 4 handles the third. Step 3 prepares
the hand-off.

## Step 1, `count-kmers`: the census

At the lengths SWGA uses, specificity is **compositional**. A 12-mer is not
selective because it binds the target more tightly; it is selective because the
sequence occurs more often in the target. So before anything else you need to
know how often every k-mer occurs in each genome, and that is all this step
does.

It shells out to Jellyfish, counting **canonical** k-mers, meaning a sequence
and its reverse complement are counted together. That is correct here: a primer
binds whichever strand presents its complement, so both occurrences are real
binding sites.

One table is written per length in your `min_k`-`max_k` range, per genome.
A design mixing lengths therefore produces several tables and, later, several
position indexes, which is what makes mixed-length designs work at all.

Each table gets a provenance sidecar recording which genome it was counted
from, as a full SHA-256. This is not bookkeeping. Without it, repointing
`fg_genomes` at a new assembly and skipping this step would build your design
from the previous organism's counts, silently. The sidecar is what makes
step 2 refuse instead.

**What you control:** `min_k`, `max_k`, and the genome paths. Nothing else in
this step is a design decision.

## Step 2, `filter`: the gates, and why each one exists

This is where most candidates die, and where the position index is built. Each
gate answers a different question.

### The frequency gates

`min_fg_freq` is a floor on target occurrence and `max_bg_freq` a ceiling on
host occurrence. Together they are the compositional criterion from above, and
they do most of the work: on the bundled plasmid example the background gate
alone removes about half the candidates.

**The gotcha worth knowing.** Both thresholds are rescaled per primer length by
`4**(10-k)`, because a 15-mer has exponentially fewer exact matches than a
10-mer in the same genome. For the foreground floor that is plainly right —
without it every long primer would be rejected. For the background ceiling the
argument is weaker, and past k=13 the gate effectively means "zero exact host
matches" and stops responding to your parameter at all. The filter therefore
prints the effective threshold and the site count it implies for every length
it runs, and warns when a ceiling has fallen below one site. Read those lines;
they tell you what was actually enforced.

### The thermodynamic gate

`min_tm` and `max_tm` bound the effective melting temperature, corrected for
salts and any additives you configured. An oligo whose Tm is far below the
reaction temperature does not prime; one far above may bind the host as
readily as the target.

If you do not set these, they are resolved **from the polymerase**, not from a
fixed default. That matters: a Bst design run at 63 C through a window built
for phi29 would keep primers that cannot prime and reject ones that can.

### Composition and structure

GC content, homopolymer runs, a 3' GC clamp and self-dimer screening. These are
synthesis and mispriming concerns rather than selectivity ones. The GC window
adapts to your target's own GC content, which matters on extreme genomes: a
fixed lower bound would exclude exactly the zero-GC primers that published
AT-rich designs are built from.

### The evenness gate

`max_gini` rejects a primer whose binding sites are clustered. A primer with
fifty sites all in one region amplifies that region and nothing else.

Evenness is only measurable from enough sites. Below `min_gini_sites` (default
3) the index is treated as **unmeasured** rather than scored, because one site
gives no gap and two give a single gap whose Gini is identically 0.0 — the best
value available. Before that rule existed, unmeasurable primers ranked first.
On a small target where most primers bind once, lower this to 1 or 2 or you
will lose nearly the whole pool.

### What this step also builds

The **position index**: for every retained candidate, where it binds each
genome. Everything downstream measures against this rather than re-scanning,
and a candidate missing from it would score as covering nothing and binding the
host nowhere. The pipeline refuses rather than reporting a coverage figure that
describes only part of the pool.

## Step 3, `prepare-candidates`: ordering, not scoring

This step used to be called `score` and used to run a random forest to predict
amplification. **It no longer does**, and the rename has no alias precisely so
that old examples fail loudly rather than quietly doing something else.

The model was retired because it did not earn its place: its gate removed 7 of
1222 candidates on one panel and none at all on two others, every downstream
consumer read only the primer column, and taking the top half of a pool by its
score produced a worse design than three random halves. It was also fitted to
synthetic data from a hand-written rule, so it reproduced an opinion rather
than measured amplification. `--amp-model` brings it back if you want it.

What the step still does is write the candidate pool in a **deterministic
order**, led by step 2's own ranking. That order is not cosmetic: the optimizer
is order-sensitive, so this is what makes an unseeded run reproducible.

## Step 4, `optimize`: the set problem

Everything so far judged oligos one at a time. This step chooses a set, and the
best set is not the set of individually best oligos — two excellent primers
that bind the same region are worth less together than either is alone.

### Coverage, and the reach that defines it

A primer site "covers" the sequence a polymerase can reach from it. So coverage
depends entirely on what reach you assume, and **this tool uses two different
ones on purpose**:

| reach | value for phi29 | used for |
|---|---|---|
| realistic per-primer | ~3 kb | selecting and scoring coverage |
| single-molecule processivity | ~70 kb | amplification-network connectivity |

Selection and the reported coverage both use the realistic figure, so they
agree. If you see advice elsewhere to design to the 70 kb processivity, ignore
it: it inflates reported coverage by roughly five to twentyfold and describes a
genome your panel does not cover.

A hybrid run prints two coverage numbers, labelled `(estimated, binned)` and
`(measured)`. They disagree by construction. Read the measured one.

### The dimer constraint is hard, not tradeable

Two primers whose ends are complementary prime each other instead of the
genome. `max_dimer_bp` bounds the longest complementary run between any two
delivered primers, and selection **stops** rather than admitting a violating
pair.

This is why `num_primers` is a request and not a guarantee. A pool supports a
bounded panel size at a given threshold, and on real pools that bound is often
well below what people ask for. A short panel is reported with its reason.
`--allow-dimer-relaxation` trades the constraint for size and warns about every
admission by name.

### Choosing the size

Three tools bear on `num_primers` and they answer different questions.
`--auto-size` inverts a saturation curve and never looks at your candidate
pool; `--show-frontier` builds a real coverage-against-specificity trade-off
but stops at 20 primers; the marginal coverage table that every run prints has
no size limit and shows what each additional primer actually bought.

Coverage rises monotonically with size while selectivity peaks and then decays,
so the coverage curve alone will not tell you where to stop.

## What this tool does not do

Worth stating plainly, because each of these is a place where a number could be
over-read.

- **No amplification outcome is modelled from data.** Nothing here was fitted
  against measured sequencing yield. The retired random forest was fitted to
  synthetic labels.
- **Gap statistics are reported, never ranked.** `max_gap`, mean gap and gap
  evenness reach the summary and the report, and no optimizer scores them. That
  is deliberate: no threshold derived from the polymerase reach separates the
  18 published primer sets with known wet-lab outcomes, so this tool will not
  pick one. You can impose your own through the panel limits.
- **Background site *position* barely enters.** Selectivity is additive in
  per-primer counts, so two hosts with identical per-primer counts score alike
  whether their sites are clustered or dispersed — even though clustering is
  most of what drives off-target amplification. Only `bg_coverage` sees
  position, and nothing selects on it.
- **No background amplification network is built**, in contrast to the
  foreground one.
- **Occupancy barely discriminates on phi29 at 30 C.** Most 12-mers are bound
  essentially all the time, so thermodynamic discrimination is not available
  there. The lever that changes this is the reaction — additives, or a
  higher-temperature enzyme — not the candidate filter.
- **Some additives are accepted and change no melting temperature**, glycerol
  among them. Rather than invent a coefficient, a design that sets one is
  refused.

## Where to go next

- [SWGA_SCIENCE.md](../SWGA_SCIENCE.md) — the biology and chemistry behind the
  above, independent of this tool.
- [QUICK_START.md](QUICK_START.md) — the commands.
- [production-scenarios.md](../production-scenarios.md) — worked configurations
  by organism and situation. This is the organism-oriented parameter guide;
  start here if your question is "what should I set for my genome".
- [params-reference.md](../params-reference.md) — every parameter, generated
  from the schema.
- [TROUBLESHOOTING.md](TROUBLESHOOTING.md) — when a run fails or disappoints.
- `docs/validation/` — the measurements behind specific claims and defaults.
  Where this document says a default was chosen for a reason, that directory
  usually holds the measurement.
