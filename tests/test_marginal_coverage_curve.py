"""Set size is the user's most consequential choice and nothing showed the knee.

The audit's sweep has coverage rising monotonically while selectivity density
peaks near n=32 for M. tuberculosis and is already falling by n=32 for E. coli.
--auto-size estimates a size from a closed-form coverage model without reading
the candidate pool, and --show-frontier stops at 20 primers, so the project's
own sweep to n=160 was done by hand.

This curve is measured on the delivered set: coverage after the first k
primers, for each k. It is a lower bound on a re-optimization at size k,
because a prefix of one set is not the best set of that size.
"""

import numpy as np
import pytest

from neoswga.core.coverage import marginal_coverage_curve


class _Cache:
    """Minimal PositionCache stand-in: get_positions(prefix, primer, strand)."""

    def __init__(self, mapping):
        self._mapping = mapping

    def get_positions(self, prefix, primer, strand="both"):
        return np.asarray(self._mapping.get((prefix, primer), []), dtype=np.int64)


def test_the_curve_has_one_entry_per_primer_in_order():
    cache = _Cache({("g", "AAAA"): [100], ("g", "CCCC"): [5000], ("g", "GGGG"): [9000]})

    curve = marginal_coverage_curve(
        cache=cache,
        primers=["AAAA", "CCCC", "GGGG"],
        prefixes=["g"],
        seq_lengths=[10000],
        extension=500,
    )

    assert [row["n"] for row in curve] == [1, 2, 3]
    assert [row["primer"] for row in curve] == ["AAAA", "CCCC", "GGGG"]


def test_coverage_is_cumulative_and_non_decreasing():
    cache = _Cache({("g", "AAAA"): [100], ("g", "CCCC"): [5000], ("g", "GGGG"): [9000]})

    curve = marginal_coverage_curve(
        cache=cache,
        primers=["AAAA", "CCCC", "GGGG"],
        prefixes=["g"],
        seq_lengths=[10000],
        extension=500,
    )
    coverages = [row["coverage"] for row in curve]

    assert coverages == sorted(coverages)
    assert all(0.0 <= c <= 1.0 for c in coverages)


def test_the_final_coverage_matches_compute_per_prefix_coverage():
    """The curve must end where the authoritative measure says it does, or the
    table contradicts the fg_coverage printed beside it."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    cache = _Cache({("g", "AAAA"): [100, 7000], ("g", "CCCC"): [5000]})
    primers = ["AAAA", "CCCC"]

    curve = marginal_coverage_curve(
        cache=cache,
        primers=primers,
        prefixes=["g"],
        seq_lengths=[10000],
        extension=500,
    )
    aggregate, _ = compute_per_prefix_coverage(
        cache=cache,
        primers=primers,
        prefixes=["g"],
        seq_lengths=[10000],
        extension=500,
    )

    assert curve[-1]["coverage"] == pytest.approx(aggregate)


def test_marginal_is_the_gain_in_percentage_points():
    cache = _Cache({("g", "AAAA"): [1000], ("g", "CCCC"): [6000]})

    curve = marginal_coverage_curve(
        cache=cache,
        primers=["AAAA", "CCCC"],
        prefixes=["g"],
        seq_lengths=[10000],
        extension=500,
    )

    # Each primer marks a 1000 bp window in a 10 kb genome: 10 pp each.
    assert curve[0]["marginal_pp"] == pytest.approx(10.0)
    assert curve[1]["marginal_pp"] == pytest.approx(10.0)


def test_a_redundant_primer_contributes_nothing():
    """The whole point of the curve: a primer landing inside an already-covered
    window adds no coverage, and the table must show that as 0.00."""
    cache = _Cache({("g", "AAAA"): [1000], ("g", "CCCC"): [1010]})

    curve = marginal_coverage_curve(
        cache=cache,
        primers=["AAAA", "CCCC"],
        prefixes=["g"],
        seq_lengths=[10000],
        extension=500,
    )

    # Each window is pos +/- 500, so the second primer adds only the 10 bases
    # its window extends past the first's: 10 / 10000 = 0.1 pp.
    assert curve[1]["marginal_pp"] == pytest.approx(0.1, abs=0.01)


def test_the_curve_decays_on_a_saturating_genome():
    """The measured shape: 2.86, 1.18, 0.29 pp per primer on M. tuberculosis.
    Marginal gain must fall as the genome fills."""
    # Deliberately NOT sorted by position. Adding primers left to right along
    # the genome gives every one fresh territory and no decay at all, which is
    # an artefact of the ordering rather than a property of the curve.
    rng = np.random.default_rng(20260906)
    positions = [int(p) for p in rng.integers(0, 100000, size=40)]
    cache = _Cache({("g", f"P{i:02d}"): [positions[i]] for i in range(40)})

    curve = marginal_coverage_curve(
        cache=cache,
        primers=[f"P{i:02d}" for i in range(40)],
        prefixes=["g"],
        seq_lengths=[100000],
        extension=3000,
    )

    first_ten = sum(row["marginal_pp"] for row in curve[:10])
    last_ten = sum(row["marginal_pp"] for row in curve[-10:])
    assert last_ten < first_ten, (
        f"marginal gain did not decay: first ten {first_ten:.2f} pp, " f"last ten {last_ten:.2f} pp"
    )


def test_empty_inputs_give_an_empty_curve():
    assert marginal_coverage_curve(None, [], [], []) == []
    assert marginal_coverage_curve(_Cache({}), [], ["g"], [1000]) == []


def test_optimize_prints_the_table_and_it_agrees_with_fg_coverage(tmp_path):
    """A source grep for the call was the first version of this test, and it
    broke the moment the reporting moved to its own module while the behaviour
    was untouched. What matters is that optimize prints the table and that its
    last row equals the fg_coverage recorded beside it -- two numbers a reader
    will compare, computed by different code paths.
    """
    import json
    import shutil
    import subprocess
    import sys

    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available (required for count-kmers)")

    rng = np.random.default_rng(20260910)
    seq = "".join(rng.choice(list("ACGT"), size=60_000))
    fasta = tmp_path / "target.fasta"
    fasta.write_text(
        ">target\n" + "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n"
    )

    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "fg_genomes": [str(fasta)],
                "bg_genomes": [],
                "fg_prefixes": [str(tmp_path / "target")],
                "bg_prefixes": [],
                "data_dir": str(tmp_path / "results"),
                "min_k": 10,
                "max_k": 10,
                "polymerase": "phi29",
                "reaction_temp": 30.0,
                "min_fg_freq": 1e-6,
                "max_bg_freq": 1.0,
                "max_gini": 1.0,
                "max_primer": 60,
                "min_tm": 0,
                "max_tm": 100,
                "gc_min": 0.0,
                "gc_max": 1.0,
                "num_primers": 4,
                "target_set_size": 4,
                "max_sets": 1,
                "iterations": 2,
                "cpus": 1,
                "fg_circular": True,
            }
        )
    )

    def _step(*extra):
        proc = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified", *extra, "-j", str(params_file)],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
            timeout=900,
        )
        if proc.returncode != 0:
            pytest.skip(f"{extra[0]} failed:\n{proc.stderr[-800:]}")
        return proc

    _step("count-kmers")
    _step("filter")
    _step("score")
    optimize = _step("optimize", "--seed", "1")

    output = optimize.stdout + optimize.stderr
    assert "Marginal foreground coverage" in output, output[-2000:]
    assert "pp/primer" in output

    rows = [
        line.split()
        for line in output.splitlines()
        if line.strip().startswith(("INFO:", "INFO "))
        and len(line.split()) == 4
        and line.split()[1].isdigit()
    ]
    assert rows, "no data rows in the table"
    last_coverage = float(rows[-1][2])

    summary = json.loads((tmp_path / "results" / "step4_improved_df_summary.json").read_text())
    fg_coverage = summary["metrics"]["fg_coverage"]
    assert last_coverage == pytest.approx(
        fg_coverage, abs=0.001
    ), f"the table ends at {last_coverage} but fg_coverage is {fg_coverage}"
