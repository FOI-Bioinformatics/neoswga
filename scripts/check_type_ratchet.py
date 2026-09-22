#!/usr/bin/env python3
"""Type errors may shrink, never grow, and a new module must be clean.

A repository-wide type gate is not affordable here: 164 errors across 46 files,
most in modules nobody is editing. Turning it on would mean either a large
mechanical diff or a permanently red check, and a permanently red check is
worse than none -- it teaches people to merge past a failing gate, which is
what the doubly-suppressed dependency audit was teaching.

A ratchet is affordable, and it is the shape this repository already uses for
module size and for search loops the budget cannot see.

## Why this is a script and not a test

Two reasons, both found by trying the other way.

**The counts depend on the environment.** Measured 2026-09-22: the same commit
and the same mypy version reported 200 errors in 53 files on a developer
machine and 164 in 46 in CI's lint job, which installs the linters and nothing
else. mypy sees more when it can resolve an import, so the counts move with
which optional dependencies happen to be installed. Running this in the
six-cell test matrix compared one environment's measurement against another's
baseline and failed for a reason that was not a defect. It runs in the lint job
only, with mypy pinned exactly, for the same reason black and isort are pinned
there.

**`tests/conftest.py` imports numpy**, so pytest cannot even start in that
environment.

`tests/test_the_type_errors_can_only_shrink.py` tests `compare` against
synthetic counts, which needs no mypy and no environment at all.

## Usage

    python scripts/check_type_ratchet.py            # exit 1 on a regression
    python scripts/check_type_ratchet.py --record   # rewrite the baseline

Reproduce the measuring environment locally:

    python3.11 -m venv /tmp/lintenv && /tmp/lintenv/bin/pip install mypy==1.19.1
    /tmp/lintenv/bin/python scripts/check_type_ratchet.py
"""

import argparse
import json
import pathlib
import re
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
BASELINE = ROOT / "tests" / "type_error_baseline.json"

#: Exactly the invocation `.github/workflows/ci.yml` makes in its informational
#: mypy step. If the two drift, the gate measures something CI does not print.
MYPY_ARGS = ["neoswga", "--ignore-missing-imports", "--no-strict-optional", "--no-color-output"]

#: The version the baseline was recorded under, pinned in the lint job. An
#: upgrade moves the counts with nothing in this repository changing -- the
#: same drift that turned every model-loading test red when skops 0.15
#: narrowed its default trust list.
BASELINE_MYPY_VERSION = "1.19.1"

_ERROR = re.compile(r"^neoswga/(?P<path>[^:]+):\d+: error:")


def measure():
    """Per-file error counts, keyed by path relative to the package."""
    completed = subprocess.run(
        [shutil.which("mypy") or "mypy", *MYPY_ARGS],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=900,
    )
    counts: dict = {}
    for line in completed.stdout.splitlines():
        found = _ERROR.match(line)
        if found:
            counts[found.group("path")] = counts.get(found.group("path"), 0) + 1
    return counts


def compare(counts, baseline):
    """Every way the measurement fails the ratchet, as (kind, detail) pairs.

    Three rules, and the second is the one a repository-wide count cannot give
    you. New code is where a type error is cheapest to fix and likeliest to be
    a real defect, and it is exactly the code a total protects least.
    """
    problems = []

    for path, count in sorted(counts.items()):
        allowed = baseline.get(path)
        if allowed is None:
            problems.append(
                (
                    "unlisted",
                    f"{path}: {count} error(s), and it is not in the baseline. A new "
                    f"or renamed module starts clean; adding it to the baseline is "
                    f"not the fix.",
                )
            )
        elif count > allowed:
            problems.append(("worse", f"{path}: {count} now against {allowed} in the baseline"))

    for path in sorted(set(baseline) - set(counts)):
        problems.append(("stale", f"{path} is now clean, so its baseline entry can go"))

    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description="Type-error ratchet")
    parser.add_argument("--record", action="store_true", help="rewrite the baseline")
    args = parser.parse_args()

    if shutil.which("mypy") is None:
        print("mypy is not installed; this gate cannot run", file=sys.stderr)
        return 2

    version = subprocess.run(
        [shutil.which("mypy"), "--version"], capture_output=True, text=True
    ).stdout.strip()
    counts = measure()

    if args.record:
        BASELINE.write_text(json.dumps(dict(sorted(counts.items())), indent=2) + "\n")
        print(f"recorded {sum(counts.values())} errors in {len(counts)} files " f"under {version}")
        return 0

    problems = compare(counts, json.loads(BASELINE.read_text()))
    if not problems:
        print(f"type-error ratchet: {sum(counts.values())} errors in {len(counts)} files, held")
        return 0

    print(f"type-error ratchet FAILED ({version}):", file=sys.stderr)
    for kind, detail in problems:
        print(f"  [{kind}] {detail}", file=sys.stderr)
    if BASELINE_MYPY_VERSION not in version:
        print(
            f"\nThe baseline was recorded under mypy {BASELINE_MYPY_VERSION} and this "
            f"is {version}. A version difference moves the counts with no change to "
            f"this repository: re-record with --record in the same commit as the "
            f"pin bump.",
            file=sys.stderr,
        )
    else:
        print(
            "\nThe counts also depend on which optional dependencies are "
            "importable. Reproduce the measuring environment with a venv holding "
            f"mypy=={BASELINE_MYPY_VERSION} and nothing else; see this file's "
            "docstring.",
            file=sys.stderr,
        )
    return 1


if __name__ == "__main__":
    sys.exit(main())
