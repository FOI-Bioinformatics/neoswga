"""End-to-end test running all four pipeline steps on the plasmid example."""

import os
import shutil
import subprocess

import pytest

from neoswga.core.kmer_counter import check_jellyfish_available

PLASMID_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, "examples", "plasmid_example"
)

requires_jellyfish = pytest.mark.skipif(
    not check_jellyfish_available(), reason="Jellyfish not installed"
)


@pytest.mark.integration
@requires_jellyfish
class TestPlasmidE2E:
    """Run the full pipeline on the plasmid example and verify outputs."""

    @pytest.fixture(autouse=True)
    def setup_workdir(self, tmp_path):
        """Copy plasmid example to a temporary directory with relaxed params."""
        import json

        self.workdir = tmp_path / "plasmid"
        shutil.copytree(PLASMID_DIR, self.workdir)
        self.params = str(self.workdir / "params.json")

        # Relax min_amp_pred for small plasmid genomes where
        # the RF model scores all primers below the default threshold
        with open(self.params) as f:
            params = json.load(f)
        params["min_amp_pred"] = 0
        with open(self.params, "w") as f:
            json.dump(params, f)

    def _run(self, command):
        """Run a neoswga CLI command and assert success."""
        result = subprocess.run(
            ["neoswga", command, "-j", self.params],
            capture_output=True,
            text=True,
            cwd=str(self.workdir),
            timeout=300,
        )
        assert result.returncode == 0, (
            f"'{command}' failed (rc={result.returncode}):\n"
            f"stdout: {result.stdout[-500:]}\n"
            f"stderr: {result.stderr[-500:]}"
        )
        return result

    def test_full_pipeline(self):
        """Run count-kmers, filter, score, optimize and verify outputs."""
        # Step 1: count-kmers
        self._run("count-kmers")
        kmer_files = [
            f for f in os.listdir(self.workdir)
            if f.endswith("_all.txt")
        ]
        assert len(kmer_files) > 0, "No k-mer output files produced"

        # Step 2: filter
        self._run("filter")
        step2 = self.workdir / "step2_df.csv"
        assert step2.exists(), "step2_df.csv not produced"
        assert step2.stat().st_size > 100, "step2_df.csv is too small"

        # Step 3: score
        self._run("score")
        step3 = self.workdir / "step3_df.csv"
        assert step3.exists(), "step3_df.csv not produced"
        assert step3.stat().st_size > 100, "step3_df.csv is too small"

        # Step 4: optimize
        self._run("optimize")
        step4 = self.workdir / "step4_improved_df.csv"
        assert step4.exists(), "step4_improved_df.csv not produced"
        assert step4.stat().st_size > 50, "step4_improved_df.csv is too small"
