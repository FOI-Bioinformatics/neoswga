"""Shared fixtures for integration tests that exercise the full
count-kmers / filter / score / optimize pipeline on example genomes.
"""

import glob
import os

import pytest

_EXAMPLE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..", "examples", "plasmid_example"
)


@pytest.fixture(scope="session", autouse=True)
def _ensure_example_kmers():
    """Generate the plasmid-example k-mer count files once if they are missing.

    Most integration tests copy ``examples/plasmid_example/`` into a tmpdir and
    then run ``filter``/``score`` (step2/step3), which require the
    ``*_Nmer_all.txt`` k-mer files produced by ``count-kmers`` (step1). Those
    files are generated artifacts (gitignored), so on a clean checkout / CI they
    do not exist and the tests fail with a STEP 2 PREREQUISITE ERROR.

    This session-scoped autouse fixture runs ``count-kmers`` once in the example
    directory (the default 6-12 k range, which covers every range the
    integration tests use). It requires jellyfish; if jellyfish is unavailable
    the files are left absent and the dependent integration tests will fail with
    the pipeline's own prerequisite error (jellyfish is a hard dependency for the
    pipeline, so the integration suite genuinely needs it).
    """
    if not os.path.isdir(_EXAMPLE_DIR):
        return  # tests that need it will skip on missing EXAMPLE_DIR

    if glob.glob(os.path.join(_EXAMPLE_DIR, "*mer_all.txt")):
        return  # already present (generated earlier or committed)

    try:
        from neoswga.core.kmer_counter import check_jellyfish_available
    except Exception:
        return
    if not check_jellyfish_available():
        return  # leave absent; dependent tests surface the missing dependency

    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core import parameter

    cwd = os.getcwd()
    try:
        os.chdir(_EXAMPLE_DIR)
        pipeline_mod._initialized = False
        parameter.json_file = "params.json"
        pipeline_mod._initialize()
        pipeline_mod.step1()  # count-kmers
    except Exception:
        # Best-effort; if generation fails the dependent tests will report it.
        pass
    finally:
        # Reset pipeline global state so it does not leak into the tests that
        # drive the pipeline themselves.
        pipeline_mod._initialized = False
        os.chdir(cwd)


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
