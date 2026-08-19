"""Six enzymes are accepted; four of them break something.

The registry, the schema and every `--polymerase` choice list agree on six
enzymes. The code that consumes the choice does not: `bst3.0` and `bsu` were
added to the registry and never exercised anywhere else, and `klenow` acquired a
correct-but-tiny processivity (40 bp, replacing a wrong 10,000) that arithmetic
elsewhere was not written for.

Three faults, all reachable from a plain command line:

- `neoswga suggest --polymerase bst` exits 1 with "Missing required parameter:
  'primer_length_range'". It reads that key off `get_polymerase_params`, which
  returns the *mechanistic enzyme* parameters and has never carried it. The
  primer-length range lives in the registry, which is where every other consumer
  reads it.

- `optimize --auto-size` gave bst and bst3.0 a 42 C default reaction temperature
  from a `30.0 if polymerase == "phi29" else 42.0` two-way branch. Both have a
  hard range starting at 50 C, so `ReactionConditions` refused it and the
  mechanistic model fell back to baseline effects behind a warning -- the
  set-size recommendation was then made with no chemistry in it at all.

- `results_interpreter` scales its enrichment estimate by
  `log2(processivity / 1000)`, which is **negative** for klenow (40 bp) and bsu
  (50 bp). A negative bonus makes the estimate negative, which the floor then
  clamps to 1.0 -- so both distributive enzymes report "no enrichment" whatever
  the primer set does. The floor hid the sign error rather than containing it.
"""

import math

import pytest

from neoswga.core.registry import POLYMERASES

ENZYMES = sorted(POLYMERASES)


# ----------------------------------------------------------------------
# suggest
# ----------------------------------------------------------------------


@pytest.mark.parametrize("polymerase", ENZYMES)
def test_suggest_resolves_a_default_primer_length(polymerase):
    """The crash, at its smallest."""
    from neoswga.cli.commands import _default_primer_length_for

    length = _default_primer_length_for(polymerase)
    low, high = POLYMERASES[polymerase].primer_length_range

    assert low <= length <= high


@pytest.mark.parametrize("polymerase", ENZYMES)
def test_suggest_runs_for_every_enzyme(polymerase, capsys):
    """End to end through the handler, which is where it exited 1."""
    from neoswga.cli.commands import run_suggest

    class Args:
        genome = None
        genome_gc = 0.5
        primer_length = None
        use_optimizer = False
        optimize_for = "amplification"

        def __init__(self, polymerase):
            self.polymerase = polymerase

        def __getattr__(self, name):
            return None

    run_suggest(Args(polymerase))

    assert capsys.readouterr().out.strip()


# ----------------------------------------------------------------------
# Reaction temperature defaults
# ----------------------------------------------------------------------


@pytest.mark.parametrize("polymerase", ENZYMES)
def test_the_default_temperature_is_inside_the_enzymes_hard_range(polymerase):
    """A two-way phi29/else branch cannot cover an enzyme running at 63 C."""
    from neoswga.core.parameter import default_reaction_temp
    from neoswga.core.reaction_conditions import ReactionConditions

    temp = default_reaction_temp(polymerase)
    low, high = POLYMERASES[polymerase].temp_hard_range

    assert low <= temp <= high
    ReactionConditions(temp=temp, polymerase=polymerase)  # must not raise


def test_no_two_way_polymerase_temperature_branch_remains():
    """Guards the shape, not one instance: `X if polymerase == 'phi29' else Y`
    is wrong for any enzyme count above two, and there are six.
    """
    import ast
    from pathlib import Path

    root = Path(__file__).parent.parent / "neoswga"
    offenders = []
    for path in root.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.IfExp):
                continue
            test = node.test
            if not isinstance(test, ast.Compare) or not isinstance(test.ops[0], ast.Eq):
                continue
            names = {n.id for n in ast.walk(test) if isinstance(n, ast.Name)}
            literals = {n.value for n in ast.walk(test) if isinstance(n, ast.Constant)}
            if "polymerase" in names and literals & set(POLYMERASES):
                offenders.append(f"{path.relative_to(root)}:{node.lineno}")

    assert not offenders, (
        "two-way polymerase branch; use the registry so all six enzymes are "
        f"covered: {offenders}"
    )


# ----------------------------------------------------------------------
# Distributive enzymes
# ----------------------------------------------------------------------


@pytest.mark.parametrize("polymerase", ENZYMES)
def test_the_processivity_bonus_is_never_negative(polymerase):
    """`log2(processivity / 1000)` goes negative below 1 kb, and klenow is 40 bp.

    A negative multiplier drove the enrichment estimate below zero, where the
    1.0 floor clamped it -- so every klenow and bsu design reported "no
    enrichment" regardless of its primer set.
    """
    from neoswga.core.results_interpreter import processivity_bonus

    assert processivity_bonus(polymerase) > 0


def test_a_longer_reach_still_earns_a_larger_bonus():
    """Fixing the sign must not flatten the term into a constant."""
    from neoswga.core.results_interpreter import processivity_bonus

    assert processivity_bonus("phi29") > processivity_bonus("bst")
    assert processivity_bonus("bst") > processivity_bonus("klenow")


def test_a_distributive_enzyme_can_still_report_enrichment():
    """The visible consequence: a klenow result that is not stuck at the floor."""
    from neoswga.core.results_interpreter import processivity_bonus

    mean_factor, n_primers = 5.0, 8
    estimate = mean_factor * n_primers * processivity_bonus("klenow")

    assert estimate > 1.0
