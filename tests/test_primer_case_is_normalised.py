"""One lowercase primer must not make a perfect duplex invisible.

`_validate_candidates` checks `all(c in "ATCG" for c in primer.upper())` and
then appends the primer UNCHANGED, so a lowercase sequence passes validation and
enters the pool as-is. The dimer scan translates with a table covering both
cases and then compares characters directly, so 'a' never matches 'A':

    pair casing        longest complementary run   is_dimer(thr=3)
    both uppercase                            10              True
    both lowercase                            10              True
    one lowercased                             0             False

A full-length duplex therefore reads as no dimer at all -- in the clique
optimizer's compatibility graph and in the network optimizer's penalty alike.
This is the one input that breaks the dimer-free guarantee outright rather than
weakening it, and it is reachable through fixed primers, `expand-primers` and
the library API. Sequences the pipeline generates itself are uppercase, which is
why it has stayed latent.
"""

import pytest

from neoswga.core.dimer import is_dimer_fast, max_complementary_run
from neoswga.core.thermodynamics import reverse_complement

#: Deliberately NOT self-complementary, so the pair is two distinct
#: oligos and de-duplication cannot hide the case question.
DUPLEX = "AAACCCGGGT"
PARTNER = reverse_complement(DUPLEX)


def test_the_pair_is_a_full_duplex_when_both_are_uppercase():
    """Guard the guard."""
    assert max_complementary_run(DUPLEX, PARTNER) == len(DUPLEX)
    assert is_dimer_fast(DUPLEX, PARTNER, 3)


def test_validation_uppercases_what_it_admits():
    """The pool the optimizers see must be case-normal.

    Fixing the comparison instead would leave every other consumer of the pool
    -- Tm, GC, k-mer lookup, the HDF5 keys -- to make the same allowance
    separately.
    """
    from neoswga.core.base_optimizer import BaseOptimizer

    validated = BaseOptimizer._validate_candidates(None, [DUPLEX.lower(), PARTNER])

    assert validated == [DUPLEX, PARTNER], (
        f"validation returned {validated}; a lowercase primer was admitted "
        "unchanged and will read as non-complementary against every uppercase "
        "partner in the pool"
    )


def test_a_mixed_case_pool_still_reports_its_dimer():
    """The consequence, stated end to end."""
    from neoswga.core.base_optimizer import BaseOptimizer
    from neoswga.core.unified_optimizer import worst_heterodimer

    pool = BaseOptimizer._validate_candidates(None, [DUPLEX.lower(), PARTNER])
    length, pair = worst_heterodimer(pool)

    assert length == len(DUPLEX), (
        f"the worst pair in a mixed-case pool measured {length} bp; the duplex "
        "is being missed because the two sequences differ only in case"
    )
    assert pair is not None


def test_invalid_sequences_are_still_rejected():
    """Normalising must not start admitting non-DNA."""
    from neoswga.core.base_optimizer import BaseOptimizer

    validated = BaseOptimizer._validate_candidates(None, [DUPLEX.lower(), "ACGUX", ""])
    assert validated == [DUPLEX]


def test_duplicates_that_differ_only_in_case_collapse():
    """Two spellings of one oligo are one oligo, and ordering both is waste."""
    from neoswga.core.base_optimizer import BaseOptimizer

    validated = BaseOptimizer._validate_candidates(None, [DUPLEX, DUPLEX.lower()])
    assert validated == [DUPLEX]
