"""A dimer screen should be able to bound stability as well as run length.

`dimer.is_dimer_thermodynamic` has existed, tested, with no production caller.
Item 4 of `docs/validation/getting_ahead_on_spacing_2026-09-18.md`.

Measured before building, and the measurements decided the design
(`scripts/benchmarking/dimer_policy_comparison.py`):

- **Cost is not the obstacle.** 11.3 us per pair against 10.5 us for the run
  screen, a ratio of 1.1. The audit guessed this needed care; it does not.
- **At the shipped `max_dimer_bp` of 3 a free-energy floor decides nothing.**
  Every pair a -6 kcal/mol floor rejects, the run screen already rejects, at
  35%, 50% and 65% GC alike. So it cannot be sold as catching what the run
  screen misses at the default.
- **It is worth having because the run screen is what bounds panel size.** At
  run <= 3 a 200-primer pool supports a greedy panel of 14 to 20. At run <= 5
  with a -4 floor it supports 52 to 79, with a hard 5 bp length cap still in
  force. That is a principled alternative to `--allow-dimer-relaxation`, which
  admits arbitrary run lengths with a warning.
- **A floor ALONE is unsafe.** dG > -6 with no length cap admits complementary
  runs of 8 bp. This project delivered an 11 bp heterodimer against a
  configured 3 once already.

So the floor ADDS to `max_dimer_bp` and can never relax it, which is what
`TestTheFloorOnlyEverAdds` pins.
"""

import pytest

from neoswga.core.lazy_dimer import (
    LAZY_DIMER_POOL_THRESHOLD,
    LazyDimerCompatibility,
    dimer_screen,
)

# 8 bp of perfect complementarity: fails any run limit, and stable enough to
# fail any floor. Verified, not assumed.
STICKY_A = "GGGGCCCCAT"
STICKY_B = "GGGGCCCCTA"

# Inside a 5 bp run limit (4 bp run) and too unstable for a -6 floor at 30 C.
LOOSE_A = "ATTATTATTT"
LOOSE_B = "AAATTATTTT"

# Inside a 5 bp run limit (5 bp run) and MARGINAL against a -6 floor: a dimer
# at 20 C and not at 90 C. Found by search, because a hand-picked pair is a
# guess and two of mine were wrong.
MARGINAL_A = "CAAATCGCTG"
MARGINAL_B = "GCGATGGGAC"


def _screen(max_dimer_bp=3, **kw):
    return LazyDimerCompatibility(max_dimer_bp, **kw)


class TestNothingChangesWithoutAFloor:
    def test_the_default_carries_no_floor(self):
        assert _screen().max_dimer_dg is None

    def test_the_config_default_carries_no_floor(self):
        from neoswga.core.base_optimizer import OptimizerConfig

        assert OptimizerConfig().max_dimer_dg is None

    def test_the_parameter_global_is_unset(self):
        from neoswga.core import parameter

        assert parameter.max_dimer_dg is None

    def test_the_schema_declares_it(self):
        import json
        import pathlib

        # Anchored on this file, not the working directory: the session-scoped
        # example priming fixture chdirs, so a relative path resolves
        # differently depending on what ran first.
        root = pathlib.Path(__file__).resolve().parents[1]
        schema = json.loads(
            (root / "neoswga" / "core" / "schema" / "params.schema.json").read_text()
        )

        assert "max_dimer_dg" in schema["properties"]

    def test_a_small_pool_still_gets_the_dense_matrix(self):
        """The floor is what forces the lazy branch, so without one the dense
        matrix must still be chosen for a small pool."""
        screen = dimer_screen(["ACGTACGTAC", "TTGCATGCAT"], 3)

        assert not isinstance(screen, LazyDimerCompatibility)


class TestTheFloorBinds:
    def test_a_stable_pair_is_rejected_even_inside_the_run_limit(self):
        """The point of the floor: a pair short enough to pass the length
        screen but stable enough to prime off itself."""
        loose = _screen(max_dimer_bp=8)
        floored = _screen(max_dimer_bp=8, max_dimer_dg=-2.0)

        assert loose.dimerises(STICKY_A, [STICKY_B]) is False
        assert floored.dimerises(STICKY_A, [STICKY_B]) is True

    def test_an_unstable_pair_still_passes(self):
        floored = _screen(max_dimer_bp=5, max_dimer_dg=-6.0, temp=30.0)

        assert floored.dimerises(LOOSE_A, [LOOSE_B]) is False

    def test_the_reaction_temperature_reaches_the_model(self):
        """Free energy is dH - T*dS, so a hotter reaction destabilises a duplex
        and the same pair can pass. A screen that ignored the temperature would
        answer identically at 20 C and 90 C.

        This is the one property no published SWGA tool has: all three call the
        same melting routine with no arguments, so their dimer screens cannot
        see the reaction at all.
        """
        cold = _screen(max_dimer_bp=5, max_dimer_dg=-6.0, temp=20.0)
        hot = _screen(max_dimer_bp=5, max_dimer_dg=-6.0, temp=90.0)

        assert cold.dimerises(MARGINAL_A, [MARGINAL_B]) is True
        assert hot.dimerises(MARGINAL_A, [MARGINAL_B]) is False

    def test_a_hot_platform_does_not_raise(self):
        """`ReactionConditions` validates its temperature against the
        polymerase and refuses 63 C under the phi29 default, so the screen must
        not build one to pass a temperature through."""
        screen = _screen(max_dimer_bp=5, max_dimer_dg=-6.0, temp=63.0)

        screen.dimerises(MARGINAL_A, [MARGINAL_B])


class TestTheFloorOnlyEverAdds:
    """The safety property, and the reason a floor is not offered on its own.

    Measured: dG > -6 with no length cap admits 8 bp complementary runs. This
    project delivered an 11 bp heterodimer against a configured 3 once, and
    `max_dimer_bp` is what stops that.
    """

    def test_a_pair_over_the_run_limit_is_rejected_whatever_the_floor(self):
        for floor in (None, -0.001, -6.0, -100.0):
            screen = _screen(max_dimer_bp=3, max_dimer_dg=floor)
            assert screen.dimerises(STICKY_A, [STICKY_B]) is True, floor

    def test_a_lenient_floor_cannot_readmit_a_long_run(self):
        """A floor so lenient that nothing meets it must still leave the run
        limit in force, or the floor would be a relaxation."""
        screen = _screen(max_dimer_bp=3, max_dimer_dg=-999.0)

        assert screen.dimerises(STICKY_A, [STICKY_B]) is True


class TestTheScreenChoiceRespectsTheFloor:
    def test_a_floor_forces_the_pairwise_screen(self):
        """The dense matrix codes t-mers in a 4**8 space and cannot express free
        energy, so using it would silently ignore a configured floor. That is
        the shape of the defect `dimer_screen`'s docstring already records."""
        screen = dimer_screen(["ACGTACGTAC", "TTGCATGCAT"], 3, max_dimer_dg=-6.0)

        assert isinstance(screen, LazyDimerCompatibility)

    def test_the_floor_and_temperature_survive_the_choice(self):
        screen = dimer_screen(["ACGTACGTAC", "TTGCATGCAT"], 3, max_dimer_dg=-6.0, temp=42.0)

        assert screen.max_dimer_dg == -6.0
        assert screen.temp == 42.0

    def test_a_large_pool_without_a_floor_is_unaffected(self):
        pool = [f"ACGT{i:06d}"[:10] for i in range(LAZY_DIMER_POOL_THRESHOLD + 1)]

        screen = dimer_screen(pool, 3)

        assert isinstance(screen, LazyDimerCompatibility)
        assert screen.max_dimer_dg is None


class TestTheCacheStillWorks:
    def test_a_pair_is_computed_once_with_a_floor_in_force(self):
        screen = _screen(max_dimer_bp=3, max_dimer_dg=-6.0)

        screen.dimerises(STICKY_A, [STICKY_B])
        screen.dimerises(STICKY_B, [STICKY_A])

        assert screen.computations == 1


class TestItReachesTheSearches:
    @pytest.mark.parametrize(
        "module,attr",
        [
            ("dominating_set_optimizer", "max_dimer_dg"),
            ("network_optimizer", "max_dimer_dg"),
        ],
    )
    def test_the_optimizer_carries_the_floor(self, module, attr):
        import importlib

        mod = importlib.import_module(f"neoswga.core.{module}")
        source = mod.__file__
        with open(source) as handle:
            text = handle.read()

        assert f"self.{attr}" in text, f"{module} does not carry {attr}"

    def test_every_screen_call_site_forwards_the_floor(self):
        """Seven call sites build a screen. One that forwarded only the run
        limit would silently drop a configured floor, which is the class Known
        Issue 16 records."""
        import pathlib
        import re

        offenders = []
        for path in pathlib.Path("neoswga").rglob("*.py"):
            text = path.read_text()
            for match in re.finditer(r"(dimer_screen|LazyDimerCompatibility)\(", text):
                start = match.end()
                call = text[start : start + 220]
                depth, end = 1, 0
                for index, char in enumerate(call):
                    depth += (char == "(") - (char == ")")
                    if depth == 0:
                        end = index
                        break
                args = call[:end]
                if "max_dimer_bp: int" in args or "def " in args:
                    continue
                if "max_dimer_dg" not in args:
                    offenders.append(f"{path}: {match.group(1)}({args.strip()[:70]})")
        assert not offenders, "screen built without forwarding the floor:\n" + "\n".join(offenders)


class TestTheFloorReachesTheOptimizersThatBuildScreens:
    """Forwarding at the SCREEN sites is not enough, and I shipped that gap.

    The first version of this file asserted only that every `dimer_screen` call
    forwards the floor. It did, and a real run still ignored the floor
    entirely, because `HybridOptimizer` built its two inner optimizers without
    passing it: the screen sites forwarded an attribute that no constructor had
    set. That is Known Issue 16's shape one level up, and it survived a test
    written specifically to catch it.
    """

    PRODUCTION_SITES = {
        "neoswga/core/hybrid_optimizer.py",
        "neoswga/core/dominating_set_adapter.py",
        "neoswga/core/network_optimizer.py",
    }

    # Constructions inside functions no command reaches. Listed by the
    # enclosing function so the exemption is legible, and the list can only
    # shrink: a production path added here would defeat the test.
    EXEMPT_FUNCTIONS = {"benchmark_network_vs_ratio"}

    def test_every_production_optimizer_construction_forwards_the_floor(self):
        import ast
        import pathlib as _p

        offenders = []
        for rel in sorted(self.PRODUCTION_SITES):
            path = _p.Path(__file__).resolve().parents[1] / rel
            tree = ast.parse(path.read_text())
            exempt_lines = set()
            for parent in ast.walk(tree):
                if isinstance(parent, ast.FunctionDef) and parent.name in self.EXEMPT_FUNCTIONS:
                    exempt_lines.update(
                        range(parent.lineno, (parent.end_lineno or parent.lineno) + 1)
                    )
            for node in ast.walk(tree):
                if not isinstance(node, ast.Call):
                    continue
                name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
                if name not in {"DominatingSetOptimizer", "NetworkOptimizer"}:
                    continue
                if node.lineno in exempt_lines:
                    continue
                if not any(kw.arg == "max_dimer_dg" for kw in node.keywords):
                    offenders.append(f"{rel}:{node.lineno} {name}(...)")
        assert not offenders, (
            "an optimizer is constructed without the stability floor, so a "
            "configured floor reaches its screen as None:\n" + "\n".join(offenders)
        )

    def test_the_hybrid_forwards_its_own_floor_to_both_stages(self):
        """Both inner optimizers must carry the floor the hybrid was given.

        Driven through a real construction rather than read off
        `__init__`'s source. The source-text form of this test broke on
        2026-09-21 when the unread `NetworkOptimizer` construction moved to
        `core/selection_weights.py` -- the floor was still forwarded, and the
        test saw only that the call had left the function it was reading.
        That is the failure `attach_search_config` records from the other
        side: asserting where code sits rather than what it does.
        """
        from neoswga.core.hybrid_optimizer import HybridOptimizer

        floor = -7.25
        optimizer = HybridOptimizer(
            position_cache=None,
            fg_prefixes=["fg"],
            fg_seq_lengths=[10_000],
            bg_prefixes=[],
            bg_seq_lengths=[],
            max_dimer_dg=floor,
        )

        reached = {
            "DominatingSetOptimizer": optimizer.dominating_optimizer.max_dimer_dg,
            "NetworkOptimizer": optimizer.network_optimizer.max_dimer_dg,
        }

        assert reached == {
            "DominatingSetOptimizer": floor,
            "NetworkOptimizer": floor,
        }, f"the hybrid does not forward its own floor to both stages; got {reached}"
