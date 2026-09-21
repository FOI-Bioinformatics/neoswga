# `--application` does not change what `hybrid` selects

Measured 21 September 2026, on the shipped Wolbachia design and the plasmid
example. Found while checking an external report's claim that the melting-
temperature term "collapses a mixed-length panel by construction on every
default run".

## Result

All four application profiles deliver an **identical** panel from the default
`hybrid` method. Wolbachia wMel against *Drosophila*, 400 candidates, target
size 10, seed 1:

| `--application` | Delivered panel |
|---|---|
| clinical | 10 primers |
| discovery | identical |
| enrichment | identical |
| metagenomics | identical |

Directly on `HybridOptimizer`, on a 12-primer Wolbachia panel where Stage 2
certainly runs, and on the plasmid example:

| Weight | Effect on delivered panel |
|---|---|
| `tm_weight` 0.25 | none |
| `tm_weight` 1.0 | none |
| `uniformity_weight` 0.5 | none |
| both together | none |

## Mechanism

`HybridOptimizer.__init__` constructs a `NetworkOptimizer` and passes it both
weights. **Nothing in the package or the tests ever reads
`self.network_optimizer` again** — `grep -rn "\.network_optimizer\b"` returns
only the assignment and unrelated module imports.

Stage 2 is `HybridOptimizer._network_refine`, a method on a different object,
and it carries no Tm or uniformity term.

The comment beside the construction claimed the opposite. It read: "Omitting
these left the object that performs refinement at its own defaults, so
`--application` and a configured `max_dimer_bp` reached the adapter and stopped
there." So an earlier fix wired these parameters into an object that does not
perform refinement, and the comment recorded the intent as though it had
worked. `max_dimer_bp` is genuinely fine — it reaches selection through
`self.max_dimer_bp` on the class itself.

This is the Known Issue 8 class again, in the shape `attach_search_config`
records: both ends exist and the path does not.

## What the term would do if connected

The Tm score is a Gaussian on the effective melting temperature, peaked at
`reaction_temp + 5` and halving 3.7 °C away. One step in k moves an effective
Tm by 4–7 °C, so its output varies enormously with oligo length. Median
`tm_score` over the plasmid candidate pool at phi29 30 °C:

| k | median effective Tm | median `tm_score` |
|---|---|---|
| 7 | 20.1 | 0.000014 |
| 8 | 27.0 | 0.044 |
| 9 | 32.3 | 0.697 |
| 10 | 39.7 | 0.334 |
| 11 | 43.5 | 0.027 |

A 48,488-fold span between the best and worst length in one pool.

So the external report was right that this term would punish a mixed-length
panel severely, and wrong that it does: the term is not applied on the default
method. Its criticism pointed at something real by a mechanism that does not
operate.

It also sits in unresolved tension with `occupancy.site_occupancy`, which is
monotone increasing in Tm — longer is strictly better, 0.996 for the 11-mer
that the Gaussian scores 0.027. For a cycled reaction the peaked form is
defensible; for a strand-displacing isothermal one the case for penalising a
more stably bound primer is not made anywhere in this repository. The
justification is a one-line comment, `# Optimal annealing ~5C above reaction`.

## What was done

Nothing that moves a panel. The misleading comment is corrected, the
construction now warns when either weight is set, naming `network` as the
method where a Tm-weighted selection actually happens, and
`tests/test_selection_weights_reach_the_method.py` pins present behaviour
including an assertion that fails if anything starts reading
`self.network_optimizer`.

## The decision this leaves open

Wiring the weights into `_network_refine` would make `--application` mean
something on the default method, and would change every delivered panel. That
needs the same treatment the beam switch is getting: measure across several
targets first. It should not be done by reconnecting a wire, because the
choice between the peaked Tm model and the monotone occupancy model is
unresolved, and connecting the Gaussian would silently pick one.

## Caveat

Two designs, four profiles, two panel sizes. The weights are inert by
construction rather than by coincidence — the object is never read — so this
generalises in a way a coverage measurement would not. What does not
generalise is the 48,488-fold figure, which is specific to this pool, this
polymerase and this temperature.
