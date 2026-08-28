"""Behavioural tests for the setup wizard (`neoswga init`).

`wizard.py` is the first thing a new user runs: it reads their target genome
and writes the params.json every other command depends on. Two classes of
failure matter here more than line coverage:

1. Inert recommendations -- a "GC-adaptive" wizard that gives the same advice
   for a 19%-GC genome and a 72%-GC genome is not adaptive at all.
2. A written config that downstream code rejects, or that silently disables a
   correctness fix elsewhere in the pipeline. This module has done both
   before: it once wrote `mg_conc: 0.0` (covered by
   tests/cli/test_setup_commands.py) and, discovered while writing these
   tests, it independently recomputed the primer GC window with a plain
   +/- clamp instead of calling `parameter.adaptive_gc_window()`. Because
   `get_params()` treats an explicit `gc_min`/`gc_max` in params.json as a
   user override, that silently suppressed the pipeline's own adaptive-GC
   release for every wizard-generated config on exactly the extreme AT-rich
   and GC-rich genomes it exists to help (see
   tests/test_adaptive_filters_behaviour.py for the release itself). It also
   hardcoded `schema_version: 1` while the pipeline's current schema is 2,
   so every fresh config triggered a "results will differ from a v1 run"
   warning about scientific corrections it never actually predates.

The tests below exercise the module's own decision logic (not the CLI
plumbing, which tests/cli/test_setup_commands.py already covers) against real
synthetic genomes of varying GC content, and check the config it produces
both validates and matches what the pipeline would compute on its own.
"""

import json
import random

import pytest

from neoswga.core import parameter
from neoswga.core.param_validator import ParamValidator, ValidationLevel
from neoswga.core.parameter import CURRENT_SCHEMA_VERSION, adaptive_gc_window, default_mg_conc
from neoswga.core.wizard import SetupWizard, format_bp, run_wizard


def _make_fasta(path, gc_fraction, length=60_000, lowercase=False, seed=0):
    """Write a synthetic single-contig FASTA with a controlled GC fraction."""
    rng = random.Random(seed)
    bases = []
    for _ in range(length):
        if rng.random() < gc_fraction:
            bases.append(rng.choice("GC"))
        else:
            bases.append(rng.choice("AT"))
    seq = "".join(bases)
    if lowercase:
        seq = seq.lower()
    lines = [">contig1"] + [seq[i : i + 70] for i in range(0, len(seq), 70)]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def at_rich_fasta(tmp_path):
    """~19% GC, like P. falciparum -- the extreme-AT case the docstring cites."""
    return _make_fasta(tmp_path / "at_rich.fasta", 0.19, seed=1)


@pytest.fixture
def gc_rich_fasta(tmp_path):
    """~72% GC, like Streptomyces -- the extreme-GC mirror case."""
    return _make_fasta(tmp_path / "gc_rich.fasta", 0.72, seed=2)


@pytest.fixture
def balanced_fasta(tmp_path):
    """~50% GC, like E. coli -- the ordinary case that should not be widened."""
    return _make_fasta(tmp_path / "balanced.fasta", 0.50, seed=3)


def _validator_errors(config):
    messages = ParamValidator().validate_params(config)
    return [m for m in messages if m.level == ValidationLevel.ERROR]


# ----------------------------------------------------------------------
# Recommendations respond to genome composition (not inert)
# ----------------------------------------------------------------------


def test_recommendations_differ_between_at_rich_and_gc_rich_genomes(at_rich_fasta, gc_rich_fasta):
    """A GC-adaptive wizard that gives identical advice for 19% and 72% GC
    genomes is not adaptive -- this is the inert-logic failure mode."""
    wizard_at = SetupWizard(interactive=False)
    wizard_at.analyze_genome(str(at_rich_fasta))

    wizard_gc = SetupWizard(interactive=False)
    wizard_gc.analyze_genome(str(gc_rich_fasta))

    assert wizard_at.recommended_polymerase != wizard_gc.recommended_polymerase
    assert wizard_at.recommended_kmer_range != wizard_gc.recommended_kmer_range
    assert wizard_at.recommended_temp != wizard_gc.recommended_temp
    assert wizard_at.gc_class != wizard_gc.gc_class


def test_gc_rich_genome_recommends_additives(gc_rich_fasta):
    """Extreme-GC targets need betaine/DMSO to denature secondary structure;
    a recommendation of 0M for both would leave the reaction non-functional."""
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(gc_rich_fasta))
    assert wizard.recommended_additives["betaine_m"] > 0
    assert wizard.recommended_additives["dmso_percent"] > 0


def test_at_rich_genome_flags_a_warning(at_rich_fasta):
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(at_rich_fasta))
    assert any("AT-rich" in w for w in wizard.warnings)


def test_gc_rich_genome_flags_a_warning(tmp_path):
    """The warning threshold (>75% GC) is stricter than the extreme_gc
    classification threshold (>70%), so this needs its own fixture rather
    than reusing gc_rich_fasta (72%)."""
    fasta = _make_fasta(tmp_path / "very_gc_rich.fasta", 0.80, seed=4)
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(fasta))
    assert any("GC-rich" in w for w in wizard.warnings)


# ----------------------------------------------------------------------
# GC window must match the pipeline's own adaptive computation
# ----------------------------------------------------------------------


def test_gc_window_releases_lower_bound_for_extreme_at_rich_genome(at_rich_fasta):
    """Regression: the wizard used to compute gc_min/gc_max itself with a
    plain clamp to [0, 1], instead of calling adaptive_gc_window(). Because an
    explicit gc_min/gc_max in params.json suppresses get_params()'s own
    adaptive call, that silently re-imposed a positive lower bound and would
    have rejected the zero-GC primers real AT-rich SWGA designs use."""
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(at_rich_fasta))
    config = wizard.generate_config()

    assert config["gc_min"] == 0.0, "extreme AT-rich genome should release the GC lower bound"


def test_gc_window_releases_upper_bound_for_extreme_gc_rich_genome(gc_rich_fasta):
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(gc_rich_fasta))
    config = wizard.generate_config()

    assert config["gc_max"] == 1.0, "extreme GC-rich genome should release the GC upper bound"


def test_gc_window_matches_pipeline_adaptive_window_across_gc_range(
    at_rich_fasta, gc_rich_fasta, balanced_fasta
):
    """The wizard's written gc_min/gc_max must agree with what get_params()
    would compute on its own for the same genome_gc, since an explicit value
    in the file takes precedence over the adaptive call at run time."""
    for fasta in (at_rich_fasta, gc_rich_fasta, balanced_fasta):
        wizard = SetupWizard(interactive=False)
        wizard.analyze_genome(str(fasta))
        config = wizard.generate_config()

        expected_min, expected_max = adaptive_gc_window(config["genome_gc"], config["gc_tolerance"])
        assert config["gc_min"] == round(expected_min, 2)
        assert config["gc_max"] == round(expected_max, 2)


def test_gc_window_stays_centred_for_a_balanced_genome(balanced_fasta):
    """The widening must be confined to compositionally extreme targets --
    a balanced genome should keep an ordinary +/- tolerance window, not also
    get released to [0, 1]."""
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(balanced_fasta))
    config = wizard.generate_config()

    assert config["gc_min"] > 0.0
    assert config["gc_max"] < 1.0
    assert config["gc_min"] <= config["genome_gc"] <= config["gc_max"]


# ----------------------------------------------------------------------
# Written config must be schema-current and pass the validator
# ----------------------------------------------------------------------


def test_written_schema_version_matches_current_pipeline_schema(balanced_fasta):
    """Regression: the wizard hardcoded schema_version: 1 while the pipeline's
    current schema is 2, so every freshly generated config triggered a
    'results will differ from a v1 run' warning about scientific corrections
    (Mg2+ defaults, Klenow processivity, Tm parameters) it never predates."""
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(balanced_fasta))
    config = wizard.generate_config()

    assert config["schema_version"] == CURRENT_SCHEMA_VERSION


@pytest.mark.parametrize("fixture_name", ["at_rich_fasta", "gc_rich_fasta", "balanced_fasta"])
def test_generated_config_passes_the_validator(fixture_name, request):
    fasta = request.getfixturevalue(fixture_name)
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(fasta))
    config = wizard.generate_config()

    errors = _validator_errors(config)
    assert not errors, f"wizard's own output failed validation: {errors}"


def test_magnesium_default_is_nonzero_and_polymerase_aware(at_rich_fasta, gc_rich_fasta):
    """Regression guard for the mg_conc: 0.0 bug fixed prior to this audit
    (see tests/cli/test_setup_commands.py) -- keep it covered at the module
    level too, across both polymerase choices the wizard can select."""
    for fasta in (at_rich_fasta, gc_rich_fasta):
        wizard = SetupWizard(interactive=False)
        wizard.analyze_genome(str(fasta))
        config = wizard.generate_config()

        assert config["mg_conc"] > 0
        assert config["mg_conc"] == default_mg_conc(config["polymerase"])


def test_optimization_method_default_is_valid(balanced_fasta):
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(balanced_fasta))
    config = wizard.generate_config()

    assert config["optimization_method"] in VALID_OPTIMIZATION_METHODS


# ----------------------------------------------------------------------
# User overrides actually take effect
# ----------------------------------------------------------------------


def test_polymerase_override_changes_reaction_temp_and_tm_bounds(gc_rich_fasta):
    """A user can override the recommended polymerase (e.g. force phi29 on a
    GC-rich genome for cost reasons). The reaction temperature and Tm bounds
    in the written config must follow the chosen polymerase, not the original
    recommendation -- otherwise the config would ask for a phi29 reaction run
    at equiphi29's 42C."""
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(gc_rich_fasta))
    assert wizard.recommended_polymerase == "equiphi29"

    wizard.user_overrides["polymerase"] = "phi29"
    config = wizard.generate_config()

    # Against the registry rather than repeated literals. The wizard used to
    # carry its own `10.0 if phi29 else 20.0` / `45.0 if phi29 else 55.0`, which
    # matched no enzyme's `primer_tm_range` -- and because an explicit key
    # suppresses the polymerase-aware default in `get_params`, writing those
    # numbers pinned a window nothing else in the tool agreed with.
    assert config["polymerase"] == "phi29"
    assert config["reaction_temp"] == parameter.default_reaction_temp("phi29")
    assert (config["min_tm"], config["max_tm"]) == parameter.default_tm_range("phi29")


def test_kmer_range_override_is_respected(balanced_fasta):
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(balanced_fasta))
    wizard.user_overrides["min_k"] = 9
    wizard.user_overrides["max_k"] = 16

    config = wizard.generate_config()
    assert config["min_k"] == 9
    assert config["max_k"] == 16


def test_num_primers_override_sets_both_size_fields(balanced_fasta):
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(balanced_fasta))
    wizard.user_overrides["num_primers"] = 12

    config = wizard.generate_config()
    assert config["num_primers"] == 12
    assert config["target_set_size"] == 12


# ----------------------------------------------------------------------
# Case sensitivity (soft-masked FASTA)
# ----------------------------------------------------------------------


def test_lowercase_soft_masked_genome_gives_the_same_recommendation(tmp_path):
    """Soft-masked genomes (repeat regions in lowercase) are common inputs.
    The GC classification feeding every downstream recommendation must not
    depend on case."""
    upper = _make_fasta(tmp_path / "upper.fasta", 0.72, seed=5, lowercase=False)
    lower = _make_fasta(tmp_path / "lower.fasta", 0.72, seed=5, lowercase=True)

    wizard_upper = SetupWizard(interactive=False)
    wizard_upper.analyze_genome(str(upper))

    wizard_lower = SetupWizard(interactive=False)
    wizard_lower.analyze_genome(str(lower))

    assert wizard_upper.gc_class == wizard_lower.gc_class
    assert wizard_upper.recommended_polymerase == wizard_lower.recommended_polymerase
    assert wizard_upper.target_stats["gc_content"] == pytest.approx(
        wizard_lower.target_stats["gc_content"]
    )


# ----------------------------------------------------------------------
# Background and blacklist genomes are recorded correctly
# ----------------------------------------------------------------------


def test_background_and_blacklist_paths_are_written_absolute(
    balanced_fasta, at_rich_fasta, gc_rich_fasta
):
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(str(balanced_fasta))
    wizard.analyze_background(str(at_rich_fasta))
    wizard.analyze_blacklist([str(gc_rich_fasta)])

    config = wizard.generate_config()

    assert config["bg_genomes"] == [str(at_rich_fasta.resolve())]
    assert config["bl_genomes"] == [str(gc_rich_fasta.resolve())]
    assert config["bl_penalty"] > 0


# ----------------------------------------------------------------------
# End-to-end via run_wizard (non-interactive)
# ----------------------------------------------------------------------


def test_run_wizard_non_interactive_writes_a_config_that_validates(tmp_path, balanced_fasta):
    out = tmp_path / "params.json"
    results_dir = tmp_path / "results"

    config = run_wizard(
        genome_path=str(balanced_fasta),
        output_path=str(out),
        output_dir=str(results_dir),
        interactive=False,
        auto_approve=True,
    )

    assert out.exists()
    written = json.loads(out.read_text())
    assert written == config
    assert not _validator_errors(written)


def test_run_wizard_creates_output_directory(tmp_path, balanced_fasta):
    out = tmp_path / "params.json"
    results_dir = tmp_path / "nested" / "results"

    run_wizard(
        genome_path=str(balanced_fasta),
        output_path=str(out),
        output_dir=str(results_dir),
        interactive=False,
        auto_approve=True,
    )

    assert results_dir.exists()


def test_run_wizard_raises_on_missing_genome(tmp_path):
    with pytest.raises(FileNotFoundError):
        run_wizard(
            genome_path=str(tmp_path / "does_not_exist.fasta"),
            output_path=str(tmp_path / "params.json"),
            interactive=False,
            auto_approve=True,
        )


def test_run_wizard_raises_on_missing_background(tmp_path, balanced_fasta):
    with pytest.raises(FileNotFoundError):
        run_wizard(
            genome_path=str(balanced_fasta),
            background_path=str(tmp_path / "no_such_background.fasta"),
            output_path=str(tmp_path / "params.json"),
            interactive=False,
            auto_approve=True,
        )


# ----------------------------------------------------------------------
# format_bp
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "length,expected",
    [(500, "500 bp"), (1500, "1.5 kb"), (2_500_000, "2.50 Mb")],
)
def test_format_bp(length, expected):
    assert format_bp(length) == expected


def test_extreme_gc_warning_fires_where_the_pipeline_actually_adapts(tmp_path):
    """The warning thresholds must not be narrower than the behaviour change.

    `adaptive_gc_window` releases the GC bound below EXTREME_AT_GENOME_GC (0.30)
    and above EXTREME_GC_GENOME_GC (0.70). The warnings used a separate,
    narrower pair (0.20 / 0.75), so a genome at 22% or 27% GC had its filtering
    silently widened with nothing said about it -- exactly the targets where a
    user most needs to know the tool is adapting.

    Four different "extreme GC" thresholds exist in this codebase. This pins the
    two that must agree: the one that changes behaviour and the one that
    explains it.
    """
    import random

    from neoswga.core.parameter import EXTREME_AT_GENOME_GC, EXTREME_GC_GENOME_GC
    from neoswga.core.wizard import SetupWizard

    def genome_at(target_gc, name):
        rng = random.Random(17)
        seq = "".join(
            rng.choice("GC") if rng.random() < target_gc else rng.choice("AT") for _ in range(4_000)
        )
        path = tmp_path / f"{name}.fasta"
        path.write_text(f">{name}\n{seq}\n")
        return str(path)

    # Just inside the adaptive band: 0.27 is AT-extreme to the pipeline but was
    # above the old 0.20 warning threshold.
    wizard = SetupWizard(interactive=False)
    wizard.analyze_genome(genome_at(0.27, "at_rich"))

    gc = wizard.target_stats["gc_content"]
    if gc < EXTREME_AT_GENOME_GC:
        assert any(
            "AT-rich" in w for w in wizard.warnings
        ), f"no warning at {gc:.1%} GC, where the pipeline releases the GC bound"
