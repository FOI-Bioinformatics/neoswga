"""Shared fixtures for integration tests that exercise the full
count-kmers / filter / score / optimize pipeline on example genomes.
"""

import pytest


@pytest.fixture
def reset_pipeline_state():
    """Reset module-level globals on neoswga.core.pipeline + neoswga.core.parameter
    so a test can drive the pipeline from a clean slate.

    The pipeline caches per-run state on the module objects (fg/bg
    prefixes, sequence lengths, circular flags, initialisation marker).
    Tests that run step2/step3/step4 sequentially or invoke the pipeline
    multiple times must reset this state between calls; otherwise each
    successive call inherits the previous test's parameters and reads
    stale k-mer files. The fixture returns a callable so tests can
    re-reset between sequential stages with the same params file.
    """
    def _reset(params_file: str) -> None:
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

    return _reset
