#!/usr/bin/env python3
"""Where does the panel assessment disagree with the post-optimization validator?

Both verdicts are written for every run. `validation["ok"]` is the validator's
and `validation["assessment"]["qualified"]` is the assessment's. They check
different things, so they are expected to differ; this records on which runs and
WHY, before either is allowed to gate more than it already does.

The `why` is the point. "The assessment is stricter about panel size" and "the
assessment catches a limit the validator misses" are different findings with
different consequences, and a bare disagreement count conflates them.

Run:
    python scripts/benchmarking/assessment_vs_validator.py
"""

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

VALIDATION = "step4_improved_df_validation.json"

#: Each case is (label, source directory, panel sizes, extra params).
#: The plasmid is primed by the test suite and runs in about 4 s.
CASES = [
    # n=1 is the control. The plasmid pool supports one primer at this
    # dimer threshold, so it is the only size where the panel MEETS the
    # request -- and therefore the only row that isolates the size rule
    # from everything else the two records check.
    ("plasmid", Path("examples/plasmid_example"), [1, 2, 4, 6], {}),
    ("plasmid+limit", Path("examples/plasmid_example"), [6], {"min_selectivity_density": 1e9}),
    # Wolbachia is deliberately absent: its shipped indexes carry no
    # reference digest, so `optimize` refuses before any panel exists. The
    # measurement document records the attempt and the exact refusal.
]

#: Which side each violation phrase comes from. The assessment states its
#: violations in prose, so classifying them here is what turns a count into a
#: finding. An unrecognised phrase is reported as such rather than bucketed.
REASON_KINDS = (
    ("size", ("is below the requested",)),
    ("dimer", ("heterodimer", "complementary run")),
    ("configured limit", ("below minimum", "above maximum", "exceeds")),
)


def classify(violation):
    lowered = violation.lower()
    for kind, needles in REASON_KINDS:
        if any(needle in lowered for needle in needles):
            return kind
    return f"unclassified: {violation}"


def run_one(source, size, overrides, workdir):
    shutil.copytree(source, workdir)
    params_path = workdir / "params.json"
    if not params_path.exists():
        return None, ["no params.json in the source directory"]
    params = json.loads(params_path.read_text())
    params["num_primers"] = size
    params["target_set_size"] = size
    params.update(overrides)
    params_path.write_text(json.dumps(params, indent=2))

    completed = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "optimize", "-j", "params.json"],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    report = workdir / VALIDATION
    if not report.exists():
        tail = [line for line in completed.stderr.strip().splitlines() if line.strip()]
        return None, tail[-1:] or ["no validation artifact and no message"]
    return json.loads(report.read_text()), []


def main():
    rows = []
    for label, source, sizes, overrides in CASES:
        if not (source / "params.json").exists():
            print(f"skipping {label}: no params.json under {source}")
            continue
        for size in sizes:
            with tempfile.TemporaryDirectory() as tmp:
                payload, why = run_one(source, size, overrides, Path(tmp) / "run")
            if payload is None:
                rows.append((label, size, "no artifact", "", "; ".join(why)))
                continue
            assessment = payload.get("assessment") or {}
            violations = assessment.get("violations") or []
            kinds = sorted({classify(v) for v in violations})
            rows.append(
                (
                    label,
                    size,
                    payload.get("ok"),
                    assessment.get("qualified"),
                    ", ".join(kinds),
                )
            )

    print()
    print(f"{'case':<16}{'n':>4}  {'validator ok':<14}{'qualified':<12}assessment faulted")
    print("-" * 88)
    for label, size, ok, qualified, kinds in rows:
        print(f"{label:<16}{size:>4}  {str(ok):<14}{str(qualified):<12}{kinds}")

    real = [r for r in rows if r[2] != "no artifact"]
    disagree = [r for r in real if r[2] != r[3]]
    print()
    print(f"runs with an artifact: {len(real)} of {len(rows)}")
    print(f"disagreements: {len(disagree)}")
    if disagree:
        kinds = sorted({k for r in disagree for k in r[4].split(", ") if k})
        print(f"reasons behind them: {', '.join(kinds)}")


if __name__ == "__main__":
    main()
