"""Behavioural tests for the multi-genome capability.

`multi-genome` designs one primer set that amplifies several related targets
while avoiding backgrounds and blacklists -- pan-genome design. It sat at 24%
coverage and nothing constructed the pipeline, only inspected its result
dataclass.

Running it end to end found that it did not finish at all: an uncapped O(n^2)
heterodimer check ran over every candidate that passed the per-primer filters.
On a 15 kb toy genome that is ~15,000 candidates and ~10^8 pairs; on a 5 Mb
bacterial genome roughly 10^10. The command appeared to hang on any real input.

These tests use genomes small enough to run in seconds while still exercising
the real Jellyfish counting, filtering, scoring and selection path.
"""

import json
import math
import random
import shutil

import pytest

from neoswga.core.multi_genome_filter import (
    MAX_ENRICHMENT,
    GenomeRole,
    GenomeSet,
    MultiGenomeFilter,
)


def _write_fasta(path, name, sequence):
    path.write_text(
        f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    return str(path)


@pytest.fixture(scope="module")
def genomes(tmp_path_factory):
    """Two related targets sharing a core, plus an unrelated background."""
    rng = random.Random(7)
    core = "".join(rng.choice("ACGT") for _ in range(9_000))
    tail = lambda: "".join(rng.choice("ACGT") for _ in range(1_000))

    root = tmp_path_factory.mktemp("multi_genome")
    return {
        "t1": _write_fasta(root / "target1.fasta", "t1", core + tail()),
        "t2": _write_fasta(root / "target2.fasta", "t2", core + tail()),
        "bg": _write_fasta(
            root / "background.fasta",
            "bg",
            "".join(rng.choice("ACGT") for _ in range(10_000)),
        ),
        "root": root,
    }


@pytest.fixture
def genome_set(genomes):
    gs = GenomeSet()
    gs.add_genome("t1", genomes["t1"], "target")
    gs.add_genome("t2", genomes["t2"], "target")
    gs.add_genome("bg", genomes["bg"], "background")
    return gs


# ----------------------------------------------------------------------
# GenomeSet
# ----------------------------------------------------------------------


def test_role_accepts_the_enum_it_exports(genomes):
    """`GenomeRole` is exported next to `GenomeSet`, so passing it is natural.

    It used to raise `AttributeError: 'GenomeRole' object has no attribute
    'lower'` from inside the setup call, before any work started.
    """
    gs = GenomeSet()
    gs.add_genome("t1", genomes["t1"], GenomeRole.TARGET)
    gs.add_genome("bg", genomes["bg"], GenomeRole.BACKGROUND)

    assert [g.name for g in gs.targets] == ["t1"]
    assert [g.name for g in gs.backgrounds] == ["bg"]


def test_role_still_accepts_strings(genome_set):
    assert len(genome_set.targets) == 2
    assert len(genome_set.backgrounds) == 1


def test_unknown_role_is_rejected_with_the_valid_options(genomes):
    gs = GenomeSet()
    with pytest.raises(ValueError, match="target"):
        gs.add_genome("x", genomes["t1"], "decoy")


def test_blacklist_gets_a_heavier_penalty_than_background(genomes):
    """The whole point of three roles: blacklist avoidance is stronger."""
    gs = GenomeSet()
    gs.add_genome("t", genomes["t1"], GenomeRole.TARGET)
    gs.add_genome("bg", genomes["bg"], GenomeRole.BACKGROUND)
    gs.add_genome("bl", genomes["bg"], GenomeRole.BLACKLIST)

    assert gs.blacklists[0].penalty_weight > gs.backgrounds[0].penalty_weight


def test_a_set_without_a_target_is_invalid(genomes):
    gs = GenomeSet()
    gs.add_genome("bg", genomes["bg"], "background")
    with pytest.raises(ValueError, match="target"):
        gs.validate()


def test_duplicate_names_are_rejected(genomes):
    gs = GenomeSet()
    gs.add_genome("same", genomes["t1"], "target")
    gs.add_genome("same", genomes["t2"], "target")
    with pytest.raises(ValueError, match="[Dd]uplicate"):
        gs.validate()


# ----------------------------------------------------------------------
# Enrichment: no background binding is the best case, not an error
# ----------------------------------------------------------------------


def _filter_with_counts(genome_set, target_counts, background_counts):
    mg = MultiGenomeFilter(genome_set=genome_set)
    mg.load_genome_counts("t1", target_counts, 10_000)
    mg.load_genome_counts("t2", target_counts, 10_000)
    mg.load_genome_counts("bg", background_counts, 10_000)
    return mg


def test_zero_background_gives_a_finite_enrichment(genome_set):
    """It used to be `inf`, which is not valid JSON and broke ranking."""
    mg = _filter_with_counts(genome_set, {"ACGTACGTAC": 5}, {})
    score = mg.score_primer("ACGTACGTAC")

    assert math.isfinite(score.enrichment_score)
    assert score.enrichment_score == MAX_ENRICHMENT


def test_zero_background_still_clears_a_min_enrichment_threshold(genome_set):
    """The sentinel must behave like inf did for thresholding, or every
    perfectly selective primer would suddenly be filtered out."""
    mg = MultiGenomeFilter(genome_set=genome_set, min_enrichment=10.0)
    mg.load_genome_counts("t1", {"ACGTACGTAC": 5}, 10_000)
    mg.load_genome_counts("t2", {"ACGTACGTAC": 5}, 10_000)
    mg.load_genome_counts("bg", {}, 10_000)

    assert mg.score_primer("ACGTACGTAC").passes


def test_ranking_orders_equally_selective_primers_by_target_binding(genome_set):
    """With inf enrichment the composite was inf for every zero-background
    primer, so ordering among the best candidates was arbitrary.

    This is the bug that mattered: on a background sharing few k-mers with the
    target, most of the pool has zero background binding.
    """
    mg = _filter_with_counts(genome_set, {"AAACCCGGGT": 20, "TTTGGGCCCA": 2}, {})
    ranked = [p for p, _s in mg.rank_primers(["TTTGGGCCCA", "AAACCCGGGT"])]

    assert ranked[0] == "AAACCCGGGT", "the better-binding primer did not rank first"


def test_a_primer_binding_neither_scores_zero(genome_set):
    mg = _filter_with_counts(genome_set, {}, {})
    assert mg.score_primer("ACGTACGTAC").enrichment_score == 0.0


def test_background_binding_reduces_enrichment(genome_set):
    clean = _filter_with_counts(genome_set, {"ACGTACGTAC": 10}, {})
    dirty = _filter_with_counts(genome_set, {"ACGTACGTAC": 10}, {"ACGTACGTAC": 10})

    assert (
        dirty.score_primer("ACGTACGTAC").enrichment_score
        < clean.score_primer("ACGTACGTAC").enrichment_score
    )


# ----------------------------------------------------------------------
# The pipeline end to end
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def pipeline_result(genomes, tmp_path_factory):
    """Run the whole pipeline once. Skips without Jellyfish."""
    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available")

    from neoswga.core.multi_genome_pipeline import MultiGenomePipeline

    gs = GenomeSet()
    gs.add_genome("t1", genomes["t1"], GenomeRole.TARGET)
    gs.add_genome("t2", genomes["t2"], GenomeRole.TARGET)
    gs.add_genome("bg", genomes["bg"], GenomeRole.BACKGROUND)

    out = tmp_path_factory.mktemp("mg_out")
    pipeline = MultiGenomePipeline(
        genome_set=gs, output_dir=str(out), kmer_range=(10, 10), primer_count=5, cpus=1
    )
    result = pipeline.run(verbose=False)
    pipeline.save_results(result)
    return {"result": result, "dir": out, "pipeline": pipeline}


def test_pipeline_completes_and_selects_primers(pipeline_result):
    """It did not previously finish at all: the heterodimer check ran over
    ~15,000 candidates, which is ~10^8 pairs."""
    result = pipeline_result["result"]

    assert result.primers, "no primers selected"
    assert len(result.primers) <= 5
    assert all(set(p) <= set("ACGT") for p in result.primers)


def test_pipeline_runs_in_seconds_not_hours(pipeline_result):
    """Guards the O(n^2) regression specifically.

    The pairwise heterodimer check now runs on the ranked shortlist rather than
    the whole candidate pool. If it moves back, this genome takes tens of
    minutes instead of tens of seconds.
    """
    assert pipeline_result["result"].total_runtime < 120


def test_selected_primers_are_mutually_compatible(pipeline_result):
    """The heterodimer check was deferred, not dropped -- it now runs where it
    matters, among primers that will share a reaction."""
    from neoswga.core import dimer, parameter

    primers = pipeline_result["result"].primers
    max_bp = getattr(parameter, "max_dimer_bp", 3) or 3

    for i, a in enumerate(primers):
        for b in primers[i + 1 :]:
            assert not dimer.is_dimer_fast(a, b, max_bp), f"{a} and {b} dimerise"


def test_results_json_is_valid_json(pipeline_result):
    """`inf` enrichment made this unparseable by anything but Python."""

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    path = pipeline_result["dir"] / "results.json"
    assert path.is_file()
    payload = json.loads(path.read_text(), parse_constant=reject)
    assert payload["primers"]


def test_reported_enrichment_is_finite(pipeline_result):
    result = pipeline_result["result"]
    assert math.isfinite(result.mean_enrichment)
    assert math.isfinite(result.min_enrichment)


def test_pipeline_writes_the_documented_outputs(pipeline_result):
    written = {p.name for p in pipeline_result["dir"].iterdir() if p.is_file()}
    assert {"results.json", "primers.fasta", "protocol.txt", "summary.txt"} <= written


def test_primers_fasta_matches_the_selection(pipeline_result):
    fasta = (pipeline_result["dir"] / "primers.fasta").read_text()
    for primer in pipeline_result["result"].primers:
        assert primer in fasta
