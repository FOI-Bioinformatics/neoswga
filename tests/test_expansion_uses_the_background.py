"""expand-primers was handed a host genome and never looked at it.

Fixing the position cache to hold the background prefixes was necessary and not
sufficient. Measured on the plasmid example afterwards, a real `expand-primers`
run queried the foreground prefix 136 times and the background prefix ZERO
times: the data was available and nothing asked for it.

Three defects, all covered here.

1. `HybridOptimizer.background_pruning` defaults to False and the expansion path
   never set it, so the one stage that reads `bg_prefixes` never ran.
2. `PrimerExpander.expand` recognised only `hybrid`, `two-stage` and
   `dominating-set`. The other four methods the CLI advertises -- including
   `background-aware`, which is precisely what a user picks when they care
   about the host -- fell through to a branch that ran hybrid instead. The
   query count with `--optimization-method background-aware` was identical to
   the default: 136 and zero.
3. `_prune_background` has no notion of fixed primers, and `optimize` hands it
   `fixed_primers + newly_selected`. Enabling pruning without fixing that would
   let expansion DROP primers from the panel the user asked to expand, which is
   worse than ignoring the background.
4. Fixing 1-3 made the run READ the background without ACTING on it. Stage 1.5
   pruning is not the stage that chooses the panel; Stage 2 `_network_refine`
   is, and it ranked on amplification connectivity and unique coverage bins
   alone. Enabling the background stage therefore only shrank the pool Stage 2
   drew from, and a background-blind choice over a smaller pool is not a better
   choice: on a 40-candidate expansion over a 300 kb synthetic pair it moved
   delivered host binding from 32 sites to 45. Stage 2 now carries a host term,
   gated so that plain `hybrid` is untouched.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.position_cache import PositionCache
from neoswga.core.primer_expansion import PrimerExpander
from neoswga.core.utility import reverse_complement

# The target has to be large enough, and the candidate pool deep enough, that
# dropping one primer does not cost coverage nobody else supplies. A first
# version of this fixture used 10 primers with three sites each on 60 kb, where
# every primer was the sole coverer of its own bins; there the coverage term
# protects a host binder as firmly as a clean one, and the fixture cannot show
# a difference that is really there. 40 primers with six sites each on 300 kb
# gives the redundancy that makes preferring a clean primer cheap.
GENOME = 300_000
N_PRIMERS = 40
FINAL_COUNT = 8


@pytest.fixture(scope="module")
def planted():
    import random

    rng = random.Random(20260910)
    primers = []
    while len(primers) < N_PRIMERS:
        p = "".join(rng.choice("ACGT") for _ in range(10))
        if p not in primers:
            primers.append(p)

    fg = list("".join(rng.choice("ACGT") for _ in range(GENOME)))
    bg = list("".join(rng.choice("ACGT") for _ in range(GENOME)))

    for i, primer in enumerate(primers):
        for j in range(6):
            pos = (i * 700 + j * 4_900) % (GENOME - 20)
            fg[pos : pos + 10] = list(primer)

    # Half the primers also bind the host, densely. They are the ones a
    # background-aware run should avoid, and their foreground layout is the
    # same stride as the other half's, so avoiding them is close to free.
    half = N_PRIMERS // 2
    dirty = primers[:half]
    for i, primer in enumerate(dirty):
        for j in range(12):
            pos = (i * 1_100 + j * 2_300) % (GENOME - 20)
            bg[pos : pos + 10] = list(primer)

    return {
        "primers": primers,
        "dirty": dirty,
        "clean": primers[half:],
        "fg": "".join(fg),
        "bg": "".join(bg),
    }


@pytest.fixture
def cache_and_prefixes(tmp_path, planted):
    def write(prefix, sequence):
        with h5py.File(f"{prefix}_10mer_positions.h5", "w") as fh:
            for primer in planted["primers"]:
                for key in {primer, reverse_complement(primer)}:
                    hits, i = [], sequence.find(key)
                    while i != -1:
                        hits.append(i)
                        i = sequence.find(key, i + 1)
                    if hits:
                        fh.create_dataset(key, data=np.array(hits, dtype=np.int64))

    fg, bg = str(tmp_path / "fg"), str(tmp_path / "bg")
    write(fg, planted["fg"])
    write(bg, planted["bg"])
    return PositionCache([fg, bg], planted["primers"]), fg, bg


def _expander(cache, fg, bg):
    return PrimerExpander(
        position_cache=cache,
        fg_prefixes=[fg],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[bg],
        bg_seq_lengths=[GENOME],
    )


def test_the_fixture_really_has_a_dirty_and_a_clean_half(cache_and_prefixes, planted):
    """Guard the guard: without this, a background nobody binds would make
    every assertion below pass for the wrong reason.

    The clean half is not required to be at exactly zero. A 10-mer occurs by
    chance about 0.29 times per strand in 300 kb, so twenty of them collect
    roughly a dozen incidental host sites however carefully the sequence is
    generated. What the fixture needs is a wide separation, which is what is
    asserted: chance hits must not be mistakable for the planted ones.
    """
    cache, _fg, bg = cache_and_prefixes
    dirty = sum(len(cache.get_positions(bg, p, "both")) for p in planted["dirty"])
    clean = sum(len(cache.get_positions(bg, p, "both")) for p in planted["clean"])

    assert dirty > 0
    assert clean * 10 < dirty, f"planted {dirty} host sites, incidental {clean}"


def test_the_expansion_optimizer_reads_the_background(cache_and_prefixes, planted):
    """The regression, stated as a query count.

    Zero reads of the background prefix is what the whole finding was.
    """
    cache, fg, bg = cache_and_prefixes
    reads = {"bg": 0}
    real = type(cache).get_positions

    def spy(self, prefix, primer, strand="both"):
        if prefix == bg:
            reads["bg"] += 1
        return real(self, prefix, primer, strand)

    type(cache).get_positions = spy
    try:
        optimizer = _expander(cache, fg, bg)._build_hybrid_optimizer()
        optimizer.optimize(planted["primers"], 4, verbose=False)
    finally:
        type(cache).get_positions = real

    assert reads["bg"] > 0, "the background prefix was never queried"


def test_background_aware_is_honoured_rather_than_silently_substituted(cache_and_prefixes, planted):
    """`--optimization-method background-aware` ran plain hybrid.

    In this codebase background-aware IS hybrid with pruning switched on, so
    honouring it is a real behaviour, not a rename.
    """
    cache, fg, bg = cache_and_prefixes
    expander = _expander(cache, fg, bg)

    built = expander._build_hybrid_optimizer(background_pruning=True)
    assert built.background_pruning is True


def test_an_unsupported_method_is_refused_rather_than_replaced(cache_and_prefixes, planted):
    """Running something other than what was asked for, with a warning, is how
    a flag comes to mean nothing."""
    cache, fg, bg = cache_and_prefixes
    expander = _expander(cache, fg, bg)

    with pytest.raises(ValueError) as excinfo:
        expander.expand(
            candidates=planted["primers"],
            fixed_primers=planted["primers"][:2],
            target_new=2,
            optimization_method="clique",
            verbose=False,
        )

    message = str(excinfo.value)
    assert "clique" in message
    assert "hybrid" in message, "the message must name what IS supported"


def test_pruning_never_removes_a_fixed_primer(cache_and_prefixes, planted):
    """Expansion must not drop the panel it was asked to expand.

    `optimize` hands `_prune_background` the fixed primers together with the
    newly selected ones, and pruning ranks purely on background per unit of
    coverage. A fixed primer that binds the host heavily is exactly what it
    would remove first.
    """
    cache, fg, bg = cache_and_prefixes
    optimizer = _expander(cache, fg, bg)._build_hybrid_optimizer(background_pruning=True)

    # Fix the two dirtiest primers: the ones pruning most wants to drop.
    fixed = planted["dirty"][:2]
    kept, _coverage, _bg = optimizer._prune_background(
        list(planted["primers"]), target_size=4, verbose=False, fixed_primers=fixed
    )

    assert set(fixed) <= set(kept), f"pruning removed a fixed primer: {fixed} -> {kept}"


def test_pruning_still_removes_unfixed_background_binders(cache_and_prefixes, planted):
    """The guard must not turn pruning into a no-op."""
    cache, fg, bg = cache_and_prefixes
    optimizer = _expander(cache, fg, bg)._build_hybrid_optimizer(background_pruning=True)

    kept, _coverage, _bg = optimizer._prune_background(
        list(planted["primers"]), target_size=5, verbose=False, fixed_primers=[]
    )

    assert len(kept) < len(planted["primers"])


# ---------------------------------------------------------------------------
# The delivered panel, which is the only thing a user sees.
#
# Reading the background is not the same as acting on it. With the three fixes
# above in place, a 40-candidate expansion measured on a 300 kb synthetic pair
# gave 32 host sites with the background stage OFF and 45 with it ON: enabling
# background awareness made the delivered panel LESS specific.
#
# The cause is that `_prune_background` is Stage 1.5 and does not choose the
# panel. Stage 2, `_network_refine`, cuts the pruned set down to the requested
# size on amplification connectivity and unique coverage bins alone. Pruning
# therefore only shrinks the pool Stage 2 draws from, and a background-blind
# choice over a smaller pool can land anywhere -- including on more host
# binders than it would have picked from the larger one.
# ---------------------------------------------------------------------------


def _delivered_host_sites(cache, fg, bg, planted, background_pruning):
    optimizer = _expander(cache, fg, bg)._build_hybrid_optimizer(
        background_pruning=background_pruning
    )
    result = optimizer.optimize(list(planted["primers"]), final_count=FINAL_COUNT, verbose=False)
    panel = list(result.primers)
    sites = sum(len(cache.get_positions(bg, p, "both")) for p in panel)
    return panel, sites


def test_the_delivered_panel_binds_the_host_less_when_the_host_is_considered(
    cache_and_prefixes, planted
):
    """The property the command exists to preserve.

    Every clean primer in this fixture has the same number of foreground sites,
    laid out on the same stride, as every dirty one. Preferring a clean primer
    therefore costs no coverage: there is no trade-off to lose, only a choice
    to make. A run that cannot tell the two apart lands on host binders at
    chance.
    """
    cache, fg, bg = cache_and_prefixes

    blind_panel, blind_sites = _delivered_host_sites(cache, fg, bg, planted, False)
    aware_panel, aware_sites = _delivered_host_sites(cache, fg, bg, planted, True)

    assert (
        blind_sites > 0
    ), "the blind run picked no host binder, so this fixture cannot show a difference"
    assert aware_sites < blind_sites, (
        f"background awareness did not reduce host binding: "
        f"{blind_sites} sites blind ({blind_panel}) vs {aware_sites} aware ({aware_panel})"
    )


def test_a_background_blind_run_is_left_exactly_as_it_was(tmp_path, planted):
    """The background term must not leak into plain `hybrid`.

    `hybrid` is the default pipeline method and does not enable background
    pruning. Its panels must not move underneath users who did not ask for a
    background-aware run, so this states the gate as a property rather than as
    a repeat call: the same foreground scored against a host the primers hammer
    and against a host they barely touch has to give the identical panel.

    Asserting the two calls merely agree with each other would pass even if the
    gate were removed entirely.
    """
    import random

    rng = random.Random(4242)
    empty_host = "".join(rng.choice("ACGT") for _ in range(GENOME))

    def write(prefix, sequence):
        with h5py.File(f"{prefix}_10mer_positions.h5", "w") as fh:
            for primer in planted["primers"]:
                for key in {primer, reverse_complement(primer)}:
                    hits, i = [], sequence.find(key)
                    while i != -1:
                        hits.append(i)
                        i = sequence.find(key, i + 1)
                    if hits:
                        fh.create_dataset(key, data=np.array(hits, dtype=np.int64))

    fg = str(tmp_path / "fg")
    loud = str(tmp_path / "loud")
    quiet = str(tmp_path / "quiet")
    write(fg, planted["fg"])
    write(loud, planted["bg"])
    write(quiet, empty_host)

    panels = []
    for host in (loud, quiet):
        cache = PositionCache([fg, host], planted["primers"])
        panel, _sites = _delivered_host_sites(cache, fg, host, planted, False)
        panels.append(panel)

    assert panels[0] == panels[1]
