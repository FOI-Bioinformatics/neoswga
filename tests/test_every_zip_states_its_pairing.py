"""Every `zip()` in the package says whether its inputs must be the same length.

`zip()` truncates to the shorter input and says nothing. In a tool that pairs
one list of prefixes with another of genome lengths, of genomes, of profiles or
of frequencies, that silence is the defect this repository keeps meeting: a
favourable answer standing in for an unmeasured one.

It has bitten here already. `compute_per_prefix_coverage` paired
`fg_prefixes` with `fg_seq_lengths`, and a params.json whose two lists
disagreed dropped the last target from coverage entirely -- reporting a higher
figure computed over fewer genomes, with nothing to say a target was missing.

So `strict=` is required on every call, and the interesting half is not
`strict=True`. It is that the three calls which pass `strict=False` had to
write down WHY the lengths legitimately differ, which is a judgement no
automatic fix could have made:

- `reach_calibration._held_out_error` pairs `edges` with `edges[1:]`, which
  differ by one by construction.
- `secondary_structure` compares two terminal slices, and a primer shorter than
  the requested terminal length yields fewer bases. Short oligos are what this
  tool designs.
- `primer_expansion` builds a reference manifest, where a prefix with no genome
  must be ABSENT rather than paired with a guess -- the same decision
  `DesignRequest.reference_manifest` documents.

This test is the ratchet, not the lint rule: CI runs ruff against the package,
and this runs against the package and the scripts that ship beside it, and
fails with the reason rather than the code.
"""

import ast
import pathlib

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "neoswga"


def zip_calls_without_strict(tree):
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
        if name != "zip":
            continue
        if not any(kw.arg == "strict" for kw in node.keywords):
            found.append(node.lineno)
    return found


def test_no_zip_in_the_package_leaves_its_pairing_unstated():
    offenders = []
    for path in sorted(PACKAGE.rglob("*.py")):
        lines = zip_calls_without_strict(ast.parse(path.read_text(encoding="utf-8")))
        offenders.extend(f"{path.relative_to(PACKAGE.parent)}:{line}" for line in lines)

    assert not offenders, (
        "these zip() calls do not say whether their inputs must be the same "
        "length, so an input that is short truncates the pairing silently:\n  "
        + "\n  ".join(offenders)
    )


def test_every_deliberate_truncation_is_explained():
    """A `strict=False` is a claim that unequal lengths are correct here.

    Requiring the word nearby is a weak check and it is the right strength: it
    cannot verify the reason, but it does stop `strict=False` becoming the
    quiet way to silence the lint rule, which is the only failure mode that
    would put the original defect back.
    """
    unexplained = []
    for path in sorted(PACKAGE.rglob("*.py")):
        lines = path.read_text(encoding="utf-8").splitlines()
        for index, line in enumerate(lines):
            if "strict=False" not in line:
                continue
            window = "\n".join(lines[max(0, index - 10) : index + 1]).lower()
            if "deliberate" not in window:
                unexplained.append(f"{path.relative_to(PACKAGE.parent)}:{index + 1}")

    assert not unexplained, (
        "strict=False asserts that unequal lengths are correct at this site. "
        "Say why within the ten lines above it:\n  " + "\n  ".join(unexplained)
    )


def test_the_check_would_catch_a_new_one():
    """Load-bearing, verified in place rather than by mutating the package."""
    tree = ast.parse("a = zip(xs, ys)\nb = zip(xs, ys, strict=True)\n")

    assert zip_calls_without_strict(tree) == [1]
