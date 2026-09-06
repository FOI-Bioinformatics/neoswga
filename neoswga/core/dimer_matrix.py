"""The pairwise dimer relation for a pool of primers, as one boolean matrix.

Measured on the real 449-candidate E. coli pool (`runs/gc_tiers/mid_ecoli/step3_df.csv`)
on an unloaded machine, `dimer.is_dimer_fast` costs 12.9 us a pair against
`secondary_structure.check_heterodimer`'s 627.1 us a pair -- the substring test is
49x cheaper. At 12.9 us, `is_dimer_fast` is not itself the bottleneck; it is merely
too slow to run 4.5 million times (the pair count of a 3000-candidate pool) inside
an optimizer that calls it once per alternative set. The thermodynamic screen at
627 us a pair is the dominant cost this module removes: about 1.0 minute pairwise
for the substring test against about 47 minutes for the thermodynamic one, over
those 4.5 million pairs. Neither screen could run over a whole pool at that cost,
so the greedy that needed one did not have one.

The reformulation here is exact. `is_dimer_fast(a, b, m)` is True when the
longest common substring of `a` and the reverse complement of `b` is longer than
`m`. Writing `t = m + 1`, such a substring exists exactly when some `t`-mer of
`a` equals some `t`-mer of `revcomp(b)`, which is exactly when some `t`-mer of
`a` is the reverse complement of some `t`-mer of `b`.

So the relation is a product of two indicator matrices over the `4**t` possible
`t`-mers, and numpy computes the whole thing at once. At the default
`max_dimer_bp` of 3 that is 256 columns.
"""

import logging
from dataclasses import dataclass
from typing import Dict, Iterator, List, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

_BASE_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
_COMPLEMENT_CODE = {0: 3, 1: 2, 2: 1, 3: 0}

# Above this many t-mer codes the indicator matrices stop being the cheap
# option. 4**8 is 65536 columns, which is 196 MB of bool for a 3000-primer
# pool. Callers asking for a larger max_dimer_bp fall back to the pairwise
# function rather than allocating that.
MAX_CODES = 4**8


def _codes(sequence: str, t: int) -> List[int]:
    """Integer codes of every t-mer in `sequence`, base-4, A=0 C=1 G=2 T=3.

    Any t-mer containing a character outside ACGT is skipped. A primer with an
    ambiguity code cannot be said to dimerise at that position, and treating an
    N as a wildcard match would flag pairs that do not form.
    """
    out: List[int] = []
    n = len(sequence)
    for start in range(n - t + 1):
        code = 0
        ok = True
        for offset in range(t):
            base = _BASE_CODE.get(sequence[start + offset])
            if base is None:
                ok = False
                break
            code = code * 4 + base
        if ok:
            out.append(code)
    return out


def _revcomp_code(code: int, t: int) -> int:
    """The code of the reverse complement of the t-mer with this code."""
    out = 0
    for _ in range(t):
        out = out * 4 + _COMPLEMENT_CODE[code % 4]
        code //= 4
    return out


@dataclass
class DimerMatrix:
    """Which primers in a pool dimerise with which others.

    `pairs` is symmetric and its diagonal is meaningless: a primer against
    itself is a self-dimer, governed by `max_self_dimer_bp`, and is not this
    class's question.
    """

    pairs: np.ndarray
    index: Dict[str, int]
    max_dimer_bp: int

    def dimerises(self, primer: str, others: Sequence[str]) -> bool:
        """Whether `primer` dimerises with any of `others`."""
        i = self.index.get(primer.upper())
        if i is None or not others:
            return False
        cols = [self.index[o.upper()] for o in others if o.upper() in self.index]
        if not cols:
            return False
        return bool(self.pairs[i, cols].any())

    def flagged_pairs(self) -> Iterator[Tuple[int, int]]:
        """Every dimerising pair, as `i < j` index tuples."""
        rows, cols = np.nonzero(np.triu(self.pairs, k=1))
        for i, j in zip(rows.tolist(), cols.tolist()):
            yield i, j


def build(primers: Sequence[str], max_dimer_bp: int) -> DimerMatrix:
    """The dimer relation over `primers` at the configured threshold."""
    t = max_dimer_bp + 1
    n_codes = 4**t
    upper = [str(p).upper() for p in primers]
    n = len(upper)

    longest = max((len(s) for s in upper), default=0)
    if t > longest:
        # No primer is long enough to contain a run of t bases, so no pair can
        # dimerise at this threshold. Short-circuit rather than allocating the
        # code space: params.schema.json allows max_dimer_bp up to 15, which
        # would be 4**16 columns.
        empty = np.zeros((n, n), dtype=bool)
        return DimerMatrix(
            pairs=empty, index={s: i for i, s in enumerate(upper)}, max_dimer_bp=max_dimer_bp
        )

    if n_codes > MAX_CODES:
        raise ValueError(
            f"max_dimer_bp={max_dimer_bp} needs {n_codes} t-mer codes, above the "
            f"{MAX_CODES} this representation allocates. Use dimer.is_dimer_fast "
            f"pairwise for a threshold this loose."
        )

    contains = np.zeros((n, n_codes), dtype=bool)
    contains_rc = np.zeros((n, n_codes), dtype=bool)
    for i, sequence in enumerate(upper):
        for code in _codes(sequence, t):
            contains[i, code] = True
            contains_rc[i, _revcomp_code(code, t)] = True

    # One matmul over uint8 rather than a boolean product, because numpy
    # dispatches the integer form to BLAS and the boolean form to a Python-level
    # loop over the object dtype.
    pairs = (contains.astype(np.uint8) @ contains_rc.astype(np.uint8).T) > 0
    pairs = pairs | pairs.T
    np.fill_diagonal(pairs, False)

    index = {sequence: i for i, sequence in enumerate(upper)}
    logger.debug(
        "Dimer matrix: %d primers, t=%d, %d dimerising pairs",
        n,
        t,
        int(np.triu(pairs, k=1).sum()),
    )
    return DimerMatrix(pairs=pairs, index=index, max_dimer_bp=max_dimer_bp)
