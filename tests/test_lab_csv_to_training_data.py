"""Tests for scripts/lab_csv_to_training_data.py."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

import lab_csv_to_training_data as bridge  # noqa: E402

from neoswga.core.reaction_conditions import ReactionConditions  # noqa: E402


@pytest.fixture
def fg_fixture(tmp_path):
    p1 = "ATCGATCGATCGAT"
    p2 = "GCTAGCTAGCTAGC"
    p3 = "AAAGGGCCCTTTAA"
    spacer = "ACGT" * 25
    seq = spacer + p1 + spacer + p2 + spacer + p3 + spacer + p1 + spacer
    fasta = tmp_path / "fg.fna"
    fasta.write_text(f">chr1\n{seq}\n")
    return fasta, [p1, p2, p3]


@pytest.fixture
def conditions():
    return ReactionConditions(temp=42.0, polymerase="equiphi29", na_conc=50.0, mg_conc=10.0)


def test_per_primer_mode_aggregates_replicates(fg_fixture, conditions, tmp_path):
    fasta, (p1, p2, p3) = fg_fixture
    lab = tmp_path / "lab.csv"
    lab.write_text(
        "primer,enrichment_fold\n"
        f"{p1},10.0\n"
        f"{p1},20.0\n"
        f"{p2},5.0\n"
        f"{p3},2.0\n"
    )

    df = bridge.build_training_data(
        lab_csv=lab,
        fg_genome_fasta=fasta,
        conditions=conditions,
        mode="per-primer",
    )

    assert len(df) == 3
    assert set(df["primer"]) == {p1, p2, p3}

    row = df.set_index("primer").loc[p1]
    assert row["mean_enrichment"] == pytest.approx(15.0)
    assert row["log_enrichment"] == pytest.approx(np.log10(15.0))
    assert row["num_sets"] == 2

    # AdvancedFeatureEngineer contributed feature columns
    feature_cols = [c for c in df.columns if c not in bridge.LABEL_COLUMNS and c != "primer"]
    assert len(feature_cols) >= 50


def test_per_set_mode_distributes_label_to_each_primer(fg_fixture, conditions, tmp_path):
    fasta, (p1, p2, p3) = fg_fixture
    lab = tmp_path / "lab.csv"
    lab.write_text(
        "primers,enrichment_fold\n"
        f"{p1};{p2},10.0\n"
        f"{p2};{p3},20.0\n"
    )

    df = bridge.build_training_data(
        lab_csv=lab,
        fg_genome_fasta=fasta,
        conditions=conditions,
        mode="per-set",
    )

    by_primer = df.set_index("primer")
    assert by_primer.loc[p1, "mean_enrichment"] == pytest.approx(10.0)
    assert by_primer.loc[p2, "mean_enrichment"] == pytest.approx(15.0)
    assert by_primer.loc[p3, "mean_enrichment"] == pytest.approx(20.0)
    assert by_primer.loc[p2, "num_sets"] == 2


def test_log_enrichment_clipped_for_zero_enrichment(fg_fixture, conditions, tmp_path):
    fasta, (p1, _, _) = fg_fixture
    lab = tmp_path / "lab.csv"
    lab.write_text(
        "primer,enrichment_fold\n"
        f"{p1},0.0\n"
    )

    df = bridge.build_training_data(
        lab_csv=lab,
        fg_genome_fasta=fasta,
        conditions=conditions,
        mode="per-primer",
    )

    assert np.isfinite(df["log_enrichment"]).all()


def test_missing_required_columns_raises(tmp_path, fg_fixture, conditions):
    fasta, _ = fg_fixture
    lab = tmp_path / "lab.csv"
    lab.write_text("primer,whatever\nATCG,1.0\n")

    with pytest.raises(ValueError, match="enrichment_fold"):
        bridge.build_training_data(
            lab_csv=lab,
            fg_genome_fasta=fasta,
            conditions=conditions,
            mode="per-primer",
        )
