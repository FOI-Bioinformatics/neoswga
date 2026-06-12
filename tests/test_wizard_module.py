"""Coverage for neoswga.core.wizard in non-interactive mode."""

import json

import pytest

from neoswga.core import wizard as wz


def test_format_bp():
    assert wz.format_bp(500) == "500 bp"
    assert wz.format_bp(1500).endswith("kb")
    assert wz.format_bp(2_000_000).endswith("Mb")


def test_check_jellyfish_available_returns_bool():
    assert isinstance(wz.check_jellyfish_available(), bool)


def _write_fasta(path, seq):
    path.write_text(
        ">contig1\n" + "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n"
    )


@pytest.fixture
def genome(tmp_path):
    # ~3 kb balanced-GC sequence.
    seq = "ATCGATCGGC" * 300
    fa = tmp_path / "target.fasta"
    _write_fasta(fa, seq)
    return fa


def test_wizard_analyze_and_generate_config(genome, tmp_path):
    wizard = wz.SetupWizard(interactive=False)
    result = wizard.analyze_genome(str(genome))
    # analyze_genome returns recommendations; raw stats are nested under "stats".
    assert result["stats"]["length"] > 0
    assert 0.0 <= result["stats"]["gc_content"] <= 1.0
    assert result["polymerase"] is not None
    # recommendations populated on the wizard too
    assert wizard.recommended_polymerase is not None
    assert wizard.recommended_kmer_range is not None

    config = wizard.generate_config(output_dir=str(tmp_path / "results"))
    assert isinstance(config, dict)
    assert "fg_genomes" in config or "min_k" in config

    out = wizard.write_params(str(tmp_path / "params.json"), output_dir=str(tmp_path / "results"))
    written = json.loads(open(out).read())
    assert isinstance(written, dict) and written


def test_wizard_prompt_overrides_noop_when_non_interactive(genome):
    wizard = wz.SetupWizard(interactive=False)
    wizard.analyze_genome(str(genome))
    # In non-interactive mode this must not block on input().
    wizard.prompt_overrides()  # should return immediately


def test_analyze_genome_missing_file_raises():
    wizard = wz.SetupWizard(interactive=False)
    with pytest.raises(FileNotFoundError):
        wizard.analyze_genome("/nonexistent/genome.fasta")
