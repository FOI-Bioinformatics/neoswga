"""Multi-genome scenario: multiple targets + background + blacklist.

Reuses the plasmid example by duplicating the existing pcDNA.fasta target
into two "pan-target" copies so the pipeline exercises multi-fg code paths
without needing new jellyfish output. pLTR serves double duty as background
and blacklist to verify both paths.
"""

import json
import os
import shutil
import tempfile
from pathlib import Path

import pytest


EXAMPLE_DIR = Path(__file__).resolve().parent.parent.parent / "examples" / "plasmid_example"


def _reset_pipeline_state(params_file):
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter
    pipeline_mod._initialized = False
    pipeline_mod.fg_prefixes = None
    pipeline_mod.bg_prefixes = None
    pipeline_mod.fg_genomes = None
    pipeline_mod.bg_genomes = None
    pipeline_mod.fg_seq_lengths = None
    pipeline_mod.bg_seq_lengths = None
    pipeline_mod.fg_circular = None
    pipeline_mod.bg_circular = None
    parameter.json_file = params_file


@pytest.fixture
def multi_genome_workdir():
    """Lay out a tmpdir with two foreground targets (pcDNA copies) plus
    pLTR as background and as an aliased blacklist. Uses the pre-built
    k-mer count files so no jellyfish run is required."""
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("examples/plasmid_example not available")

    tmpdir = Path(tempfile.mkdtemp(prefix="neoswga_multigenome_"))

    # Copy the plasmid example files
    for fname in os.listdir(EXAMPLE_DIR):
        src = EXAMPLE_DIR / fname
        if src.is_file():
            shutil.copy2(src, tmpdir / fname)

    # Duplicate pcDNA as target_b so we have two foregrounds.
    # Symlink the k-mer files and positions.h5 rather than rerunning jellyfish.
    for k in range(6, 13):
        src_kmer = tmpdir / f"pcDNA_{k}mer_all.txt"
        dst_kmer = tmpdir / f"target_b_{k}mer_all.txt"
        if src_kmer.is_file() and not dst_kmer.exists():
            shutil.copy2(src_kmer, dst_kmer)

    # Copy positions.h5 if present
    for k in range(6, 13):
        src_h5 = tmpdir / f"pcDNA_{k}mer_positions.h5"
        dst_h5 = tmpdir / f"target_b_{k}mer_positions.h5"
        if src_h5.is_file() and not dst_h5.exists():
            shutil.copy2(src_h5, dst_h5)

    shutil.copy2(tmpdir / "pcDNA.fasta", tmpdir / "target_b.fasta")

    # Also set up pLTR as a blacklist (bl_pLTR) with k-mer files
    for k in range(6, 13):
        src = tmpdir / f"pLTR_{k}mer_all.txt"
        dst = tmpdir / f"bl_pLTR_{k}mer_all.txt"
        if src.is_file() and not dst.exists():
            shutil.copy2(src, dst)

    params_path = tmpdir / "params.json"
    with open(params_path) as fh:
        params = json.load(fh)

    # Two targets (same content, different paths), one background, one blacklist
    params["fg_genomes"] = [str(tmpdir / "pcDNA.fasta"), str(tmpdir / "target_b.fasta")]
    params["fg_prefixes"] = [str(tmpdir / "pcDNA"), str(tmpdir / "target_b")]
    params["fg_seq_lengths"] = [6157, 6157]
    params["bg_genomes"] = [str(tmpdir / "pLTR.fasta")]
    params["bg_prefixes"] = [str(tmpdir / "pLTR")]
    params["bg_seq_lengths"] = [6258]
    params["bl_genomes"] = [str(tmpdir / "pLTR.fasta")]
    params["bl_prefixes"] = [str(tmpdir / "bl_pLTR")]
    params["bl_seq_lengths"] = [6258]
    params["max_bl_freq"] = 1e-3  # Allow some tolerance
    params.setdefault("schema_version", 1)
    params.setdefault("min_amp_pred", 0.0)

    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)

    original_cwd = os.getcwd()
    os.chdir(tmpdir)
    yield tmpdir
    os.chdir(original_cwd)
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.mark.integration
@pytest.mark.slow
def test_multi_target_background_blacklist_filter_runs(multi_genome_workdir):
    """Filter step completes with 2 targets + 1 background + 1 blacklist."""
    params_file = multi_genome_workdir / "params.json"
    _reset_pipeline_state(str(params_file))

    from neoswga.core.pipeline import step2
    df = step2()

    assert (multi_genome_workdir / "step2_df.csv").is_file()
    # Primers must have survived both bg filtering and the blacklist filter
    assert len(df) > 0, "Multi-genome scenario produced no candidates"


@pytest.mark.integration
@pytest.mark.slow
def test_multi_target_parameter_propagation(multi_genome_workdir):
    """All foreground prefixes must be populated end-to-end in `parameter`."""
    from types import SimpleNamespace
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    attrs = {
        "json_file": str(multi_genome_workdir / "params.json"),
        "data_dir": None, "src_dir": None,
        "min_fg_freq": None, "max_bg_freq": None,
        "min_tm": None, "max_tm": None,
        "max_gini": None, "max_primer": None, "min_amp_pred": None,
        "cpus": None, "max_dimer_bp": None, "max_self_dimer_bp": None,
        "verbose": None, "drop_iterations": None, "iterations": None,
        "top_set_count": None, "retries": None, "max_sets": None,
        "selection_metric": None, "fg_circular": None, "bg_circular": None,
        "min_k": None, "max_k": None,
        "fasta_fore": None, "fasta_back": None,
        "kmer_fore": None, "kmer_back": None,
    }
    data = parameter.get_params(SimpleNamespace(**attrs))

    # Module globals
    assert len(parameter.fg_genomes) == 2, "2 foreground genomes expected"
    assert len(parameter.fg_prefixes) == 2
    assert len(parameter.bg_prefixes) == 1
    assert len(parameter.bl_prefixes) == 1
    assert len(parameter.bl_seq_lengths) == 1
    # fg_seq_lengths flows via the returned data dict (not module globals)
    assert len(data["fg_seq_lengths"]) == 2
    assert len(data["bg_seq_lengths"]) == 1
