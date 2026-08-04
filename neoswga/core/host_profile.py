"""Background specificity without a background genome.

Metagenomic SWGA has no host genome to design against: the background varies
between samples, or is simply unknown. Host-free mode currently drops the
background entirely, which makes every design report the same selectivity
(`MAX_SELECTIVITY`, "no background binding detected") and leaves coverage as
the only thing the optimizer can rank on. That is not host-free design, it is
design with specificity switched off.

A compositional profile is enough to restore it. If a host's k-mer statistics
are approximately Markov, the expected count of any k-mer follows from a
transition table -- 16 numbers at order 1, 64 at order 2 -- rather than from
its genome. Measured against real human chr21 12-mer counts:

    profile fitted to                    order   Spearman rho
    chr21 itself (the actual host)         2        0.776
    chr21 itself                           1        0.763
    Prevotella (a wrong organism)          1        0.579
    GC content alone                       0        0.386

Two things follow. Even a **wrong-organism** profile ranks host abundance far
better than nothing, so generic designs are possible. And higher order helps
only when the profile is the real host: order 2 beats order 1 on the true host
(0.776 vs 0.763) and loses on a foreign one (0.520 vs 0.579), because
higher-order structure is organism-specific and transfers worse. `recommended_order`
encodes that rather than fixing one value.

What does NOT work, recorded so it is not rebuilt: fitting the null to the
*target* and treating over-representation as specificity. Correlation with real
host load was rho = 0.03-0.26 with the unhelpful sign. Over-representation
inside the target measures target repeats, not host absence.

Expected counts feed the same occupancy machinery as real ones, so a synthetic
background is weighted by mismatch class and reaction conditions exactly as a
real genome is, and the two are directly comparable.

**Use it for ranking, not for an absolute selectivity number.** A profile
predicts what a *typical* k-mer of a given composition does. A primer chosen
because it is unusually rare in the background is by definition atypical, and
the model cannot know that -- being rarer than your composition predicts is
precisely the residual only the real genome carries. Measured on a set selected
against real chr21 counts, the modelled background load came out 2.4-2.8x above
the measured one, and the selectivity correspondingly low (0.85-0.99 modelled
against 2.35 measured). The ranking is sound (rho 0.888 on an unbiased sample);
the level is a composition-based expectation and will be pessimistic for any
set selected for cleanliness. Metrics computed this way report
`selectivity_mode="modelled"` so the two are never confused.
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence

logger = logging.getLogger(__name__)

BASES = "ACGT"

# Order to use when the profile came from the actual host versus from some
# other organism. See the module docstring: higher order captures more of a
# genome's own structure and correspondingly less of anything else's.
ORDER_FOR_OWN_HOST = 2
ORDER_FOR_FOREIGN_HOST = 1


@dataclass(frozen=True)
class HostProfile:
    """A Markov description of a background's sequence composition.

    `provenance` is carried, not decorative: it records what the profile was
    fitted from, so a report can say "expected counts under a vertebrate
    profile" rather than presenting a modelled background as a measured one.
    """

    name: str
    order: int
    base_freq: Dict[str, float]
    transitions: Dict[str, Dict[str, float]]
    provenance: str = ""
    fitted_length: int = 0

    def log_probability(self, kmer: str) -> float:
        """log P(kmer) under this profile. Log space: a 12-mer's probability is
        ~1e-8 and products of 12 of those underflow float64 in longer designs."""
        kmer = kmer.upper()
        if self.order == 0:
            return sum(math.log(max(self.base_freq.get(b, 0.25), 1e-12)) for b in kmer)

        total = sum(math.log(max(self.base_freq.get(b, 0.25), 1e-12)) for b in kmer[: self.order])
        for i in range(len(kmer) - self.order):
            context = kmer[i : i + self.order]
            nxt = kmer[i + self.order]
            row = self.transitions.get(context)
            p = row.get(nxt, 0.25) if row else 0.25
            total += math.log(max(p, 1e-12))
        return total

    def expected_count(self, kmer: str, genome_length: int) -> float:
        """Expected occurrences of `kmer` in a genome of this composition.

        Counts both strands and returns the count for the *canonical* k-mer,
        matching how jellyfish `-C` stores real counts: a k-mer and its reverse
        complement are one entry. Adding the two probabilities rather than
        doubling one keeps palindromes from being counted twice.
        """
        from neoswga.core.thermodynamics import reverse_complement

        positions = max(0, genome_length - len(kmer) + 1)
        if positions == 0:
            return 0.0
        rc = reverse_complement(kmer.upper())
        p_forward = math.exp(self.log_probability(kmer))
        p_reverse = p_forward if rc == kmer.upper() else math.exp(self.log_probability(rc))
        return positions * (p_forward + p_reverse)


def recommended_order(is_actual_host: bool) -> int:
    """Model order for a profile, given whether it describes the real host.

    Not a tuning knob: fitting order 2 to a foreign organism measurably ranks
    worse than order 1 (rho 0.520 vs 0.579 predicting chr21 from a Prevotella
    profile), because what the extra order captures is the wrong genome's
    structure.
    """
    return ORDER_FOR_OWN_HOST if is_actual_host else ORDER_FOR_FOREIGN_HOST


def fit_profile(
    sequences: Iterable[str],
    order: int = 1,
    name: str = "fitted",
    provenance: str = "",
) -> HostProfile:
    """Fit a Markov profile from sequence.

    Takes an iterable of chunks so a multi-GB genome can be streamed rather
    than held in memory; contexts spanning a chunk boundary are handled by
    carrying the last `order` bases forward. Non-ACGT characters break the
    context, which is what should happen at an N-run: a transition observed
    across a gap is not a transition.
    """
    if order < 0:
        raise ValueError(f"Markov order must be >= 0, got {order}")

    base_counts = {b: 0 for b in BASES}
    transition_counts: Dict[str, Dict[str, int]] = {}
    carry = ""

    for chunk in sequences:
        chunk = chunk.upper()
        for b in chunk:
            if b in base_counts:
                base_counts[b] += 1
        if order == 0:
            continue
        window = carry + chunk
        for i in range(len(window) - order):
            context = window[i : i + order]
            nxt = window[i + order]
            if nxt not in base_counts or any(c not in base_counts for c in context):
                continue
            transition_counts.setdefault(context, {b: 0 for b in BASES})[nxt] += 1
        carry = window[-order:] if len(window) >= order else window

    total_bases = sum(base_counts.values())
    if total_bases == 0:
        raise ValueError("No ACGT bases found; cannot fit a profile.")
    base_freq = {b: base_counts[b] / total_bases for b in BASES}

    transitions: Dict[str, Dict[str, float]] = {}
    for context, row in transition_counts.items():
        total = sum(row.values())
        if total:
            transitions[context] = {b: row[b] / total for b in BASES}

    return HostProfile(
        name=name,
        order=order,
        base_freq=base_freq,
        transitions=transitions,
        provenance=provenance or f"fitted from {total_bases:,} bp at order {order}",
        fitted_length=total_bases,
    )


def fit_profile_from_fasta(path: str, order: int = 1, name: str = "") -> HostProfile:
    """Fit a profile from a FASTA file, streaming it line by line."""

    def _chunks():
        with open(path) as handle:
            for line in handle:
                if not line.startswith(">"):
                    yield line.strip()

    return fit_profile(
        _chunks(),
        order=order,
        name=name or path.rsplit("/", 1)[-1],
        provenance=f"fitted from {path} at order {order}",
    )


def expected_site_load(
    primers: Sequence[str],
    profile: HostProfile,
    genome_length: int,
    conditions,
    max_mismatches: int = 1,
) -> float:
    """Occupancy-weighted background load expected under a compositional profile.

    Deliberately mirrors `occupancy.weighted_site_load`, which computes the same
    quantity from real count files, so a modelled background and a real one are
    the same number measured two ways and can be compared directly. The mismatch
    classes are the primer's actual neighbours, each with its own expected count
    under the profile -- not the perfect match's expectation reused, which would
    ignore that a mismatched variant has a different composition and therefore a
    different abundance.
    """
    from neoswga.core.mismatch_counts import _variants_at_distance
    from neoswga.core.occupancy import default_mismatch_penalty, mismatch_tm, site_occupancy
    from neoswga.core.thermodynamics import calculate_enthalpy_entropy

    penalty = default_mismatch_penalty()
    total = 0.0
    for primer in primers:
        primer = primer.upper()
        dh, _ds = calculate_enthalpy_entropy(primer)
        tm = conditions.calculate_effective_tm(primer)
        for distance in range(0, max_mismatches + 1):
            variants = {primer} if distance == 0 else _variants_at_distance(primer, distance)
            expected = sum(profile.expected_count(v, genome_length) for v in variants)
            if expected:
                total += expected * site_occupancy(
                    dh, mismatch_tm(tm, distance, penalty), conditions.temp
                )
    return total


def aggregate_loads(loads: Sequence[float], mode: str = "worst-case") -> float:
    """Combine background loads from several candidate hosts.

    The mode encodes what "background" means, and the two answers differ:

    - `worst-case` (max): one unknown host drawn from this panel. A set that is
      clean against human and filthy against mouse is not a generic set, and
      averaging would hide that. This is the right default for metagenomics and
      for designs meant to survive an unspecified host.
    - `sum`: every listed background present at once, which is the co-infection
      or defined-community case.

    `occupancy.weighted_site_load` sums across prefixes, so the existing
    behaviour is `sum`. Choosing it silently for a panel of *alternative* hosts
    would penalise a primer for background it will never meet.
    """
    values = [v for v in loads if v is not None]
    if not values:
        return 0.0
    if mode == "worst-case":
        return max(values)
    if mode == "sum":
        return sum(values)
    if mode == "mean":
        return sum(values) / len(values)
    raise ValueError(f"Unknown aggregation mode {mode!r}; expected worst-case, sum or mean.")
