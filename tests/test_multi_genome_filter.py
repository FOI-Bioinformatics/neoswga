"""Tests for multi-genome filtering module.

Validates GenomeRole, GenomeEntry, GenomeSet, MultiGenomeScore,
and MultiGenomeFilter data structures and scoring logic.
"""

import math
import pytest
from pathlib import Path
from unittest.mock import patch

from neoswga.core.multi_genome_filter import (
    GenomeRole,
    GenomeEntry,
    GenomeSet,
    MultiGenomeScore,
    MultiGenomeFilter,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@pytest.fixture
def tmp_fasta(tmp_path):
    """Create minimal FASTA files for testing."""
    def _make(name="genome.fna"):
        p = tmp_path / name
        p.write_text(">seq1\nATCGATCG\n")
        return p
    return _make


@pytest.fixture
def target_fasta(tmp_fasta):
    return tmp_fasta("target.fna")


@pytest.fixture
def background_fasta(tmp_fasta):
    return tmp_fasta("background.fna")


@pytest.fixture
def blacklist_fasta(tmp_fasta):
    return tmp_fasta("blacklist.fna")


@pytest.fixture
def genome_set_mixed(target_fasta, background_fasta, blacklist_fasta, tmp_fasta):
    """GenomeSet with two targets, one background, and one blacklist."""
    gs = GenomeSet()
    gs.add_genome("Target1", str(target_fasta), "target")
    gs.add_genome("Target2", str(tmp_fasta("target2.fna")), "target")
    gs.add_genome("Background1", str(background_fasta), "background")
    gs.add_genome("Blacklist1", str(blacklist_fasta), "blacklist")
    return gs


# ---------------------------------------------------------------------------
# GenomeRole
# ---------------------------------------------------------------------------

class TestGenomeRole:

    def test_enum_values(self):
        assert GenomeRole.TARGET.value == "target"
        assert GenomeRole.BACKGROUND.value == "background"
        assert GenomeRole.BLACKLIST.value == "blacklist"

    def test_from_string(self):
        assert GenomeRole("target") == GenomeRole.TARGET
        assert GenomeRole("background") == GenomeRole.BACKGROUND
        assert GenomeRole("blacklist") == GenomeRole.BLACKLIST

    def test_invalid_role_raises(self):
        with pytest.raises(ValueError):
            GenomeRole("invalid")


# ---------------------------------------------------------------------------
# GenomeEntry
# ---------------------------------------------------------------------------

class TestGenomeEntry:

    def test_creation_with_all_fields(self, target_fasta):
        entry = GenomeEntry(
            name="Borrelia",
            fasta_path=target_fasta,
            role=GenomeRole.TARGET,
            gc_content=0.28,
            size=1_500_000,
        )
        assert entry.name == "Borrelia"
        assert entry.fasta_path == target_fasta
        assert entry.role == GenomeRole.TARGET
        assert entry.gc_content == pytest.approx(0.28)
        assert entry.size == 1_500_000

    def test_target_penalty_weight_set_to_zero(self, target_fasta):
        entry = GenomeEntry(name="T", fasta_path=target_fasta, role=GenomeRole.TARGET)
        assert entry.penalty_weight == 0.0

    def test_background_default_penalty_weight(self, background_fasta):
        entry = GenomeEntry(
            name="BG", fasta_path=background_fasta, role=GenomeRole.BACKGROUND
        )
        assert entry.penalty_weight == 1.0

    def test_blacklist_default_penalty_weight(self, blacklist_fasta):
        entry = GenomeEntry(
            name="BL", fasta_path=blacklist_fasta, role=GenomeRole.BLACKLIST
        )
        assert entry.penalty_weight == 5.0

    def test_string_path_converted_to_path_object(self, target_fasta):
        entry = GenomeEntry(
            name="T", fasta_path=str(target_fasta), role=GenomeRole.TARGET
        )
        assert isinstance(entry.fasta_path, Path)

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            GenomeEntry(
                name="Missing",
                fasta_path=tmp_path / "nonexistent.fna",
                role=GenomeRole.TARGET,
            )


# ---------------------------------------------------------------------------
# GenomeSet
# ---------------------------------------------------------------------------

class TestGenomeSet:

    def test_add_genome_target(self, target_fasta):
        gs = GenomeSet()
        gs.add_genome("T1", str(target_fasta), "target")
        assert len(gs.targets) == 1
        assert gs.targets[0].name == "T1"

    def test_add_genome_background(self, background_fasta):
        gs = GenomeSet()
        gs.add_genome("BG", str(background_fasta), "background")
        assert len(gs.backgrounds) == 1

    def test_add_genome_blacklist(self, blacklist_fasta):
        gs = GenomeSet()
        gs.add_genome("BL", str(blacklist_fasta), "blacklist")
        assert len(gs.blacklists) == 1

    def test_mixed_roles(self, genome_set_mixed):
        gs = genome_set_mixed
        assert len(gs.targets) == 2
        assert len(gs.backgrounds) == 1
        assert len(gs.blacklists) == 1

    def test_get_all_genomes(self, genome_set_mixed):
        assert len(genome_set_mixed.get_all_genomes()) == 4

    def test_get_non_targets(self, genome_set_mixed):
        non_targets = genome_set_mixed.get_non_targets()
        assert len(non_targets) == 2
        roles = {g.role for g in non_targets}
        assert GenomeRole.TARGET not in roles

    def test_validate_no_targets_raises(self):
        gs = GenomeSet()
        with pytest.raises(ValueError, match="At least one target"):
            gs.validate()

    def test_validate_duplicate_names_raises(self, tmp_fasta):
        gs = GenomeSet()
        gs.add_genome("Dup", str(tmp_fasta("a.fna")), "target")
        gs.add_genome("Dup", str(tmp_fasta("b.fna")), "background")
        with pytest.raises(ValueError, match="Duplicate"):
            gs.validate()

    def test_validate_success(self, genome_set_mixed):
        assert genome_set_mixed.validate() is True

    def test_summary_contains_genome_names(self, genome_set_mixed):
        summary = genome_set_mixed.summary()
        assert "Target1" in summary
        assert "Background1" in summary
        assert "Blacklist1" in summary

    def test_empty_genome_set(self):
        gs = GenomeSet()
        assert len(gs.targets) == 0
        assert len(gs.backgrounds) == 0
        assert len(gs.blacklists) == 0
        assert gs.get_all_genomes() == []
        assert gs.get_non_targets() == []


# ---------------------------------------------------------------------------
# MultiGenomeScore
# ---------------------------------------------------------------------------

class TestMultiGenomeScore:

    def test_creation_and_fields(self):
        score = MultiGenomeScore(
            primer="ATCGATCG",
            target_frequency=1e-4,
            background_frequency=1e-5,
            blacklist_frequency=0.0,
            enrichment_score=10.0,
            penalty_score=1e-5,
            passes=True,
            details={"Target": 1e-4, "Background": 1e-5},
        )
        assert score.primer == "ATCGATCG"
        assert score.target_frequency == pytest.approx(1e-4)
        assert score.background_frequency == pytest.approx(1e-5)
        assert score.blacklist_frequency == pytest.approx(0.0)
        assert score.enrichment_score == pytest.approx(10.0)
        assert score.penalty_score == pytest.approx(1e-5)
        assert score.passes is True
        assert len(score.details) == 2

    def test_default_details_is_empty(self):
        score = MultiGenomeScore(
            primer="AAA",
            target_frequency=0.0,
            background_frequency=0.0,
            blacklist_frequency=0.0,
            enrichment_score=0.0,
            penalty_score=0.0,
            passes=False,
        )
        assert score.details == {}

    def test_str_representation(self):
        score = MultiGenomeScore(
            primer="ATCG",
            target_frequency=1e-4,
            background_frequency=1e-6,
            blacklist_frequency=0.0,
            enrichment_score=100.0,
            penalty_score=1e-6,
            passes=True,
        )
        text = str(score)
        assert "ATCG" in text
        assert "Passes: True" in text


# ---------------------------------------------------------------------------
# MultiGenomeFilter
# ---------------------------------------------------------------------------

class TestMultiGenomeFilter:

    def test_initialization(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        assert mgf.min_target_freq == pytest.approx(1e-5)
        assert mgf.max_background_freq == pytest.approx(1e-4)
        assert mgf.max_blacklist_freq == pytest.approx(1e-6)
        assert mgf.min_enrichment == pytest.approx(10.0)
        assert mgf.genome_set is genome_set_mixed

    def test_init_validates_genome_set(self):
        gs = GenomeSet()
        with pytest.raises(ValueError, match="At least one target"):
            MultiGenomeFilter(gs)

    def test_load_genome_counts(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        counts = {"ATCG": 100, "GCTA": 50}
        mgf.load_genome_counts("Target1", counts, 1_000_000)
        assert "Target1" in mgf.kmer_counts
        assert mgf.genome_sizes["Target1"] == 1_000_000

    def test_calculate_frequency(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        mgf.load_genome_counts("Target1", {"ATCG": 100}, 1_000_000)
        freq = mgf._calculate_frequency("ATCG", "Target1")
        assert freq == pytest.approx(1e-4)

    def test_calculate_frequency_missing_genome(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        freq = mgf._calculate_frequency("ATCG", "Unknown")
        assert freq == 0.0

    def test_calculate_frequency_missing_primer(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        mgf.load_genome_counts("Target1", {"ATCG": 100}, 1_000_000)
        freq = mgf._calculate_frequency("GGGG", "Target1")
        assert freq == 0.0


class TestMultiGenomeFilterScoring:
    """Test the scoring and filtering logic of MultiGenomeFilter."""

    @pytest.fixture
    def loaded_filter(self, genome_set_mixed):
        """Filter with pre-loaded k-mer counts for all genomes."""
        mgf = MultiGenomeFilter(
            genome_set_mixed,
            min_target_freq=1e-5,
            max_background_freq=1e-4,
            max_blacklist_freq=1e-6,
            min_enrichment=10.0,
        )
        # Target1: high binding
        mgf.load_genome_counts("Target1", {"ATCG": 200}, 1_000_000)
        # Target2: moderate binding
        mgf.load_genome_counts("Target2", {"ATCG": 100}, 1_000_000)
        # Background1: low binding
        mgf.load_genome_counts("Background1", {"ATCG": 5}, 1_000_000)
        # Blacklist1: no binding
        mgf.load_genome_counts("Blacklist1", {}, 1_000_000)
        return mgf

    def test_score_primer_target_frequency_is_mean(self, loaded_filter):
        """Target frequency is the mean across target genomes."""
        score = loaded_filter.score_primer("ATCG")
        # Mean of 200/1e6 and 100/1e6 = 1.5e-4
        assert score.target_frequency == pytest.approx(1.5e-4)

    def test_score_primer_background_frequency_is_max(self, loaded_filter):
        """Background frequency is the max across background genomes."""
        score = loaded_filter.score_primer("ATCG")
        assert score.background_frequency == pytest.approx(5e-6)

    def test_score_primer_blacklist_frequency_zero(self, loaded_filter):
        score = loaded_filter.score_primer("ATCG")
        assert score.blacklist_frequency == pytest.approx(0.0)

    def test_enrichment_ratio(self, loaded_filter):
        """Enrichment = target_freq / background_freq."""
        score = loaded_filter.score_primer("ATCG")
        expected = 1.5e-4 / 5e-6  # 30x
        assert score.enrichment_score == pytest.approx(expected)

    def test_enrichment_infinite_when_no_background_binding(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        mgf.load_genome_counts("Target1", {"ATCG": 100}, 1_000_000)
        mgf.load_genome_counts("Target2", {"ATCG": 100}, 1_000_000)
        mgf.load_genome_counts("Background1", {}, 1_000_000)
        mgf.load_genome_counts("Blacklist1", {}, 1_000_000)
        score = mgf.score_primer("ATCG")
        assert math.isinf(score.enrichment_score)

    def test_enrichment_zero_when_no_target_and_no_background(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed)
        mgf.load_genome_counts("Target1", {}, 1_000_000)
        mgf.load_genome_counts("Target2", {}, 1_000_000)
        mgf.load_genome_counts("Background1", {}, 1_000_000)
        mgf.load_genome_counts("Blacklist1", {}, 1_000_000)
        score = mgf.score_primer("ATCG")
        assert score.enrichment_score == pytest.approx(0.0)

    def test_blacklist_penalty_weight_applied(self, genome_set_mixed):
        """Blacklist penalty weight (default 5.0) amplifies blacklist frequency
        in the penalty score."""
        mgf = MultiGenomeFilter(genome_set_mixed)
        mgf.load_genome_counts("Target1", {"ATCG": 100}, 1_000_000)
        mgf.load_genome_counts("Target2", {"ATCG": 100}, 1_000_000)
        mgf.load_genome_counts("Background1", {"ATCG": 10}, 1_000_000)
        # Blacklist genome has some binding
        mgf.load_genome_counts("Blacklist1", {"ATCG": 10}, 1_000_000)

        score = mgf.score_primer("ATCG")
        bg_freq = 10 / 1_000_000
        bl_freq = 10 / 1_000_000

        # penalty = bg_freq * bg_weight + bl_freq * bl_weight
        # bg_weight for add_genome default = 1.0 (see __post_init__)
        # bl_weight for add_genome default blacklist = 5.0
        bg_entry = genome_set_mixed.backgrounds[0]
        bl_entry = genome_set_mixed.blacklists[0]
        expected_penalty = bg_freq * bg_entry.penalty_weight + bl_freq * bl_entry.penalty_weight
        assert score.penalty_score == pytest.approx(expected_penalty)
        # Blacklist contributes 5x more to penalty than background at equal frequency
        assert bl_entry.penalty_weight == pytest.approx(5.0)

    def test_primer_passes_when_criteria_met(self, loaded_filter):
        score = loaded_filter.score_primer("ATCG")
        assert score.passes

    def test_primer_fails_low_target_freq(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed, min_target_freq=1.0)
        mgf.load_genome_counts("Target1", {"ATCG": 1}, 1_000_000)
        mgf.load_genome_counts("Target2", {"ATCG": 1}, 1_000_000)
        mgf.load_genome_counts("Background1", {}, 1_000_000)
        mgf.load_genome_counts("Blacklist1", {}, 1_000_000)
        score = mgf.score_primer("ATCG")
        assert not score.passes

    def test_primer_fails_high_blacklist_freq(self, genome_set_mixed):
        mgf = MultiGenomeFilter(genome_set_mixed, max_blacklist_freq=1e-8)
        mgf.load_genome_counts("Target1", {"ATCG": 100}, 1_000_000)
        mgf.load_genome_counts("Target2", {"ATCG": 100}, 1_000_000)
        mgf.load_genome_counts("Background1", {}, 1_000_000)
        mgf.load_genome_counts("Blacklist1", {"ATCG": 10}, 1_000_000)
        score = mgf.score_primer("ATCG")
        assert not score.passes

    def test_filter_primers(self, loaded_filter):
        """filter_primers returns passing primers and all scores."""
        mgf = loaded_filter
        # Add a second primer that has no target binding (should fail)
        mgf.kmer_counts["Target1"]["GGGG"] = 0
        mgf.kmer_counts["Target2"]["GGGG"] = 0

        passing, scores = mgf.filter_primers(["ATCG", "GGGG"])
        assert "ATCG" in passing
        assert "GGGG" not in passing
        assert len(scores) == 2

    def test_details_contain_per_genome_frequencies(self, loaded_filter):
        score = loaded_filter.score_primer("ATCG")
        assert "Target1" in score.details
        assert "Target2" in score.details
        assert "Background1" in score.details


class TestMultiGenomeFilterTargetsOnly:
    """Edge case: genome set with only target genomes, no background or blacklist."""

    @pytest.fixture
    def targets_only_filter(self, tmp_fasta):
        gs = GenomeSet()
        gs.add_genome("T1", str(tmp_fasta("t1.fna")), "target")
        mgf = MultiGenomeFilter(gs)
        mgf.load_genome_counts("T1", {"ATCG": 100}, 1_000_000)
        return mgf

    def test_background_frequency_zero(self, targets_only_filter):
        score = targets_only_filter.score_primer("ATCG")
        assert score.background_frequency == pytest.approx(0.0)

    def test_blacklist_frequency_zero(self, targets_only_filter):
        score = targets_only_filter.score_primer("ATCG")
        assert score.blacklist_frequency == pytest.approx(0.0)

    def test_enrichment_infinite(self, targets_only_filter):
        score = targets_only_filter.score_primer("ATCG")
        assert math.isinf(score.enrichment_score)

    def test_penalty_zero(self, targets_only_filter):
        score = targets_only_filter.score_primer("ATCG")
        assert score.penalty_score == pytest.approx(0.0)

    def test_passes_with_sufficient_target_freq(self, targets_only_filter):
        score = targets_only_filter.score_primer("ATCG")
        assert score.passes


class TestMultiGenomeFilterRanking:
    """Test rank_primers composite scoring."""

    @pytest.fixture
    def ranking_filter(self, tmp_fasta):
        gs = GenomeSet()
        gs.add_genome("T1", str(tmp_fasta("t.fna")), "target")
        gs.add_genome("BG", str(tmp_fasta("bg.fna")), "background")
        mgf = MultiGenomeFilter(gs, min_enrichment=1.0)
        mgf.load_genome_counts("T1", {"AAA": 500, "BBB": 100}, 1_000_000)
        mgf.load_genome_counts("BG", {"AAA": 10, "BBB": 50}, 1_000_000)
        return mgf

    def test_rank_primers_order(self, ranking_filter):
        """Primer with higher target freq and lower background should rank first."""
        ranked = ranking_filter.rank_primers(["BBB", "AAA"])
        assert ranked[0][0] == "AAA"

    def test_rank_primers_top_n(self, ranking_filter):
        ranked = ranking_filter.rank_primers(["BBB", "AAA"], top_n=1)
        assert len(ranked) == 1
