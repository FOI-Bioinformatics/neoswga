"""What a primer set costs to synthesise.

Set size is the main lever on coverage -- 8 to 32 primers moves effective
coverage 0.38 to 0.54 -- so a design decision that used to be "how many
primers can the optimizer justify" is really "how many oligos are you willing
to buy". Until the second question is answerable the frontier has a missing
axis, and the tool was answering it with `len(primers) * 5.0`, which ignores
both length and modifications.

Length matters because synthesis is priced per coupling. Modifications matter
more: SWGA primers are ordered with 3' phosphorothioate (PTO) bonds because
phi29 carries a 3'->5' exonuclease that degrades unprotected primers, and a PTO
coupling is charged at a large multiple of a standard base. For a 10-mer with
the standard two PTO bonds, the modification is a substantial fraction of the
per-oligo price rather than a rounding error -- which is precisely the regime
this tool designs in, since short primers make the fixed and modification costs
dominate.

The numbers below are order-of-magnitude list prices, not quotes. Oligo pricing
varies severalfold with vendor, scale, purification and contract, so these are
exposed as parameters and every figure is reported as an estimate. What the
model is for is comparing designs -- is 32 primers twice the price of 16, or
five times? -- and that comparison is robust to the absolute level being wrong.
"""

from dataclasses import dataclass, replace
from typing import Dict, Iterable, List, Optional, Sequence

# Order-of-magnitude list prices, USD, 25 nmol desalted synthesis.
#
# Deliberately conservative and deliberately round: a false precision here would
# invite the numbers being quoted as if they came from a vendor. Override with
# `--oligo-*` flags or an OligoCostModel built from a real quote.
DEFAULT_SETUP_COST = 3.00  # per-oligo fixed charge (handling, QC, shipping share)
DEFAULT_PER_BASE_COST = 0.25  # per standard coupling
DEFAULT_PER_PTO_COST = 2.50  # per phosphorothioate linkage -- ~10x a standard base
DEFAULT_FIVE_PRIME_BLOCK_COST = 25.00  # C18/C3 spacer, charged as a modification


@dataclass(frozen=True)
class OligoCostModel:
    """Per-oligo pricing. All values USD; see module docstring on precision."""

    setup_cost: float = DEFAULT_SETUP_COST
    per_base_cost: float = DEFAULT_PER_BASE_COST
    per_pto_cost: float = DEFAULT_PER_PTO_COST
    five_prime_block_cost: float = DEFAULT_FIVE_PRIME_BLOCK_COST

    @classmethod
    def from_overrides(cls, **overrides) -> "OligoCostModel":
        """Build from a partial set of overrides, ignoring None."""
        clean = {k: v for k, v in overrides.items() if v is not None}
        return replace(cls(), **clean) if clean else cls()

    def oligo_cost(self, primer: str, pto_bonds: int = 2, five_prime_block=None) -> float:
        """Cost of one oligo.

        PTO bonds are counted as linkages, and an n-mer has only n-1 of them, so
        a request for more bonds than the sequence can carry is clamped rather
        than silently charged for -- at the 6-8 bp end of the SWGA range that is
        a reachable mistake, not a hypothetical one.
        """
        length = len(primer)
        bonds = max(0, min(int(pto_bonds), max(0, length - 1)))
        cost = self.setup_cost + self.per_base_cost * length + self.per_pto_cost * bonds
        if five_prime_block:
            cost += self.five_prime_block_cost
        return cost


@dataclass(frozen=True)
class SetCost:
    """Cost of a primer set, broken down so the driver is visible."""

    total: float
    n_primers: int
    setup: float
    bases: float
    modifications: float
    per_primer_mean: float
    pto_bonds: int
    five_prime_block: Optional[str]

    @property
    def modification_fraction(self) -> float:
        """Share of the total spent on modifications rather than sequence.

        Worth surfacing: for short SWGA primers this is routinely the largest
        component, which is not the intuition carried over from PCR primers.
        """
        return self.modifications / self.total if self.total else 0.0

    def to_dict(self) -> Dict:
        return {
            "total_usd": round(self.total, 2),
            "n_primers": self.n_primers,
            "per_primer_mean_usd": round(self.per_primer_mean, 2),
            "setup_usd": round(self.setup, 2),
            "bases_usd": round(self.bases, 2),
            "modifications_usd": round(self.modifications, 2),
            "modification_fraction": round(self.modification_fraction, 3),
            "pto_bonds": self.pto_bonds,
            "five_prime_block": self.five_prime_block,
            "estimate": True,
        }


def set_cost(
    primers: Sequence[str],
    model: Optional[OligoCostModel] = None,
    pto_bonds: int = 2,
    five_prime_block: Optional[str] = None,
) -> SetCost:
    """Cost of synthesising a primer set, itemised."""
    model = model or OligoCostModel()
    primers = list(primers)
    n = len(primers)
    if n == 0:
        return SetCost(0.0, 0, 0.0, 0.0, 0.0, 0.0, pto_bonds, five_prime_block)

    setup = model.setup_cost * n
    bases = sum(model.per_base_cost * len(p) for p in primers)
    mods = 0.0
    for primer in primers:
        bonds = max(0, min(int(pto_bonds), max(0, len(primer) - 1)))
        mods += model.per_pto_cost * bonds
        if five_prime_block:
            mods += model.five_prime_block_cost
    total = setup + bases + mods
    return SetCost(
        total=total,
        n_primers=n,
        setup=setup,
        bases=bases,
        modifications=mods,
        per_primer_mean=total / n,
        pto_bonds=pto_bonds,
        five_prime_block=five_prime_block,
    )


def cost_from_modifications(
    primers: Sequence[str], modifications, model: Optional[OligoCostModel] = None
) -> SetCost:
    """Cost using an existing `export.PrimerModifications`.

    Keeps the ordering profile and the price of that profile in one place, so a
    design exported with `--modification-profile low-input` is costed with the
    C18 spacer it will actually be ordered with.
    """
    return set_cost(
        primers,
        model=model,
        pto_bonds=getattr(modifications, "pto_bonds", 2),
        five_prime_block=getattr(modifications, "five_prime_block", None),
    )


def format_cost(cost: SetCost) -> str:
    """Render the breakdown, not just the total.

    The split is the useful part. A set of short primers spends most of its
    budget on fixed and modification charges, so buying more coverage by adding
    primers is close to linear in count and barely sensitive to their length --
    the opposite of the intuition that longer primers cost more.
    """
    if cost.n_primers == 0:
        return "No primers; nothing to cost."
    lines = [
        f"Estimated synthesis cost: ${cost.total:,.2f} "
        f"({cost.n_primers} oligos, ${cost.per_primer_mean:,.2f} each)",
        f"  setup         ${cost.setup:>9,.2f}",
        f"  bases         ${cost.bases:>9,.2f}",
        f"  modifications ${cost.modifications:>9,.2f}"
        f"   ({cost.modification_fraction:.0%} of total, "
        f"{cost.pto_bonds} PTO bond(s)"
        f"{', 5-prime ' + cost.five_prime_block if cost.five_prime_block else ''})",
        "",
        "Order-of-magnitude list prices at 25 nmol desalted, not a quote; "
        "oligo pricing varies severalfold with vendor, scale and contract. "
        "Use for comparing designs rather than for budgeting.",
    ]
    return "\n".join(lines)


def cost_per_coverage_point(cost: SetCost, coverage: float) -> float:
    """Dollars per percentage point of coverage.

    The number that makes a frontier actionable: it says what the next
    percentage point costs, which is what decides whether to buy it. Returns
    infinity for zero coverage rather than raising, so a failed design sorts
    last instead of crashing a comparison.
    """
    if coverage <= 0:
        return float("inf")
    return cost.total / (coverage * 100.0)
