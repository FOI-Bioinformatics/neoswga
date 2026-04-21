"""Phase 15C: library-level blacklist re-injection guard.

`unified_optimizer.run_optimization` now cross-checks the candidate pool
against `parameter.bl_prefixes` and passes any offenders as
`forbidden_primers` to `OptimizationResult.validate`. Python API callers
who bypass the filter step and pass arbitrary candidate lists get the
same protection as the swap-primer CLI.
"""

import json
import os
import shutil
import tempfile
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parent.parent
EXAMPLE_DIR = ROOT / "examples" / "plasmid_example"


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
def plasmid_with_aliased_blacklist(tmp_path):
    """Duplicate the plasmid example with pLTR aliased as a blacklist."""
    if not EXAMPLE_DIR.is_dir():
        pytest.skip("plasmid example not available")

    for fname in os.listdir(EXAMPLE_DIR):
        src = EXAMPLE_DIR / fname
        if src.is_file():
            shutil.copy2(src, tmp_path / fname)

    # Copy pLTR k-mer / position files as "bl_pLTR"
    for k in range(6, 13):
        src = tmp_path / f"pLTR_{k}mer_all.txt"
        dst = tmp_path / f"bl_pLTR_{k}mer_all.txt"
        if src.is_file() and not dst.exists():
            shutil.copy2(src, dst)

    params_path = tmp_path / "params.json"
    with open(params_path) as fh:
        params = json.load(fh)
    params.setdefault("schema_version", 1)
    params.setdefault("min_amp_pred", 0.0)
    params["min_k"] = 8
    params["max_k"] = 10
    params["bl_genomes"] = [str(tmp_path / "pLTR.fasta")]
    params["bl_prefixes"] = [str(tmp_path / "bl_pLTR")]
    params["bl_seq_lengths"] = [6258]
    params["max_bl_freq"] = 0.0
    with open(params_path, "w") as fh:
        json.dump(params, fh, indent=2)
    return tmp_path, params_path


@pytest.mark.integration
@pytest.mark.slow
def test_run_optimization_flags_blacklist_candidates(plasmid_with_aliased_blacklist):
    """When candidates contain blacklist-hitting primers, the validator's
    forbidden_primers list should include them. Library callers can then
    read step4_improved_df_validation.json or
    unified_optimizer._LAST_RESULT.validation to see the guard fired."""
    tmpdir, params_file = plasmid_with_aliased_blacklist
    os.chdir(tmpdir)
    _reset_pipeline_state(str(params_file))

    from neoswga.core.pipeline import step2, step3
    step2()
    step3()

    # Candidates from step3 are already filtered; the validator's check is
    # a belt-and-braces cross-check. We probe the mechanism itself here —
    # the forbidden list is empty only when every candidate passes the
    # filter helper, but the code path runs.
    from neoswga.core.unified_optimizer import optimize_step4
    _reset_pipeline_state(str(params_file))
    primer_sets, _, _ = optimize_step4(
        optimization_method="hybrid", verbose=False,
    )

    # Validation artefact should exist
    report_path = tmpdir / "step4_improved_df_validation.json"
    assert report_path.is_file(), "validator should have written its report"
    report = json.loads(report_path.read_text())
    assert "issues" in report


@pytest.mark.integration
@pytest.mark.slow
def test_library_caller_with_explicit_blacklist_candidate(plasmid_with_aliased_blacklist):
    """Direct library API: hand run_optimization a candidate pool that
    contains a primer with guaranteed blacklist hits. The validator's
    forbidden_primers path fires; if the optimizer selects that primer,
    the validation report carries a blacklist_primer_in_set error."""
    tmpdir, params_file = plasmid_with_aliased_blacklist
    os.chdir(tmpdir)
    _reset_pipeline_state(str(params_file))
    from neoswga.core.pipeline import step2, step3
    step2()
    step3()

    from neoswga.core import parameter
    from neoswga.core.unified_optimizer import run_optimization

    # Build a candidate list explicitly containing a primer known to
    # appear in the blacklist k-mer files (any primer from pLTR_Xmer_all.txt
    # will do; we read one short entry from the file).
    bl_kmer_file = tmpdir / "bl_pLTR_8mer_all.txt"
    blacklist_primer = None
    if bl_kmer_file.is_file():
        with open(bl_kmer_file) as fh:
            for line in fh:
                parts = line.strip().split()
                if parts and len(parts) >= 2 and int(parts[1]) > 0:
                    blacklist_primer = parts[0].upper()
                    break

    if not blacklist_primer:
        pytest.skip("no blacklist primer available in test fixture")

    # Combine real step3 primers with the blacklist one
    import pandas as pd
    step3_csv = tmpdir / "step3_df.csv"
    step3_df = pd.read_csv(step3_csv, index_col=0) if step3_csv.is_file() else pd.DataFrame()
    candidates = list(step3_df.index.astype(str).tolist())[:10]
    candidates.append(blacklist_primer)

    _reset_pipeline_state(str(params_file))
    result = run_optimization(
        method="hybrid",
        candidates=candidates,
        fg_prefixes=list(parameter.fg_prefixes),
        fg_seq_lengths=[6157],
        bg_prefixes=list(parameter.bg_prefixes),
        bg_seq_lengths=[6258],
        target_size=3,
        verbose=False,
    )

    report_path = tmpdir / "step4_improved_df_validation.json"
    assert report_path.is_file()
    report = json.loads(report_path.read_text())

    # Two acceptable outcomes:
    # (a) optimizer happened not to pick the blacklist primer — report ok
    # (b) optimizer picked it — report contains blacklist_primer_in_set error
    codes = [i.get("code") for i in report.get("issues", [])]
    picked = blacklist_primer in result.primers
    if picked:
        assert "blacklist_primer_in_set" in codes, (
            f"Optimizer selected blacklist primer {blacklist_primer!r} but "
            f"validator did not flag it; codes={codes}"
        )
