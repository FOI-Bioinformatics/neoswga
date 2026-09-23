"""Pairing a reference prefix with its length must refuse a mismatch.

`fg_prefixes` and `fg_seq_lengths` are two separate params.json keys that must
correspond one to one. Twenty-nine places pair them with `zip`, which stops at
the shorter sequence, so a file naming three targets and two lengths processed
two of them and said nothing.

That is not hypothetical. It was a real defect in `bam_coverage`, merged
2026-09-22: every gap in the dropped target was then reported as absent, so
`expand-primers` -- whose job is to add oligos for the gaps -- would never add
one for it. The entry points there now carry a named guard with a message
pointing at params.json.

This file covers the other twenty-odd sites, and the reason they get
`strict=True` rather than a guard each: a named guard is worth writing where a
user's mistake arrives, and `strict=True` is the right answer everywhere the
lists are already deep inside a computation. Both are loud. Silence is what
was wrong.

The AST check is the durable part. Adding a thirtieth pairing is easy and
forgetting `strict=` is easier, and nothing else in the suite would notice.
"""

import ast
import pathlib

import pytest

PACKAGE = pathlib.Path(__file__).resolve().parent.parent / "neoswga"

#: Name fragments that mark one side of the pairing this file is about.
PREFIX_HINTS = ("prefix", "prefixes")
PARTNER_HINTS = ("length", "lengths", "seq_lengths", "genome", "genomes")


def _name_of(node):
    """The readable name of a zip argument, attribute access included."""
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        return node.attr
    return ""


def _prefix_length_zips():
    """Every `zip(<...prefix...>, <...length/genome...>)` in the package."""
    for path in sorted(PACKAGE.rglob("*.py")):
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except SyntaxError:  # pragma: no cover - defensive
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if getattr(node.func, "id", None) != "zip":
                continue
            if len(node.args) < 2:
                continue
            first = _name_of(node.args[0]).lower()
            second = _name_of(node.args[1]).lower()
            if not any(h in first for h in PREFIX_HINTS):
                continue
            if not any(h in second for h in PARTNER_HINTS):
                continue
            strict = any(kw.arg == "strict" for kw in node.keywords)
            yield f"{path.relative_to(PACKAGE)}:{node.lineno}", strict


def test_the_check_finds_the_pairings_it_is_about():
    """Guard the guard: a broken matcher would make the next test vacuous."""
    found = list(_prefix_length_zips())

    assert len(found) >= 20, f"only {len(found)} pairings found; the matcher is broken"


def test_no_prefix_length_pairing_can_truncate_silently():
    """The rule. `strict=` must be given, either value, deliberately."""
    silent = [where for where, strict in _prefix_length_zips() if not strict]

    assert not silent, (
        "these pair a reference prefix with its length and will silently drop "
        "a target if the two params.json lists disagree; pass strict=True:\n  "
        + "\n  ".join(silent)
    )


# ---------------------------------------------------------------------------
# The behaviour, not just the spelling
# ---------------------------------------------------------------------------


def test_coverage_refuses_a_mismatch_rather_than_covering_fewer_targets():
    """A reachable function, so this is not only an AST assertion.

    Silently covering two of three targets reports a coverage figure for a
    design the user did not describe, which is worse than failing.
    """
    from neoswga.core.coverage import compute_per_prefix_coverage

    class Cache:
        def get_positions(self, prefix, primer, strand):
            return []

    with pytest.raises(ValueError):
        compute_per_prefix_coverage(Cache(), ["ACGTACGTACGT"], ["a", "b", "c"], [1000, 2000], 3000)


def test_coverage_still_works_when_the_lists_agree():
    """Guard the guard again: refusing everything would pass the test above."""
    from neoswga.core.coverage import compute_per_prefix_coverage

    class Cache:
        def get_positions(self, prefix, primer, strand):
            return []

    result = compute_per_prefix_coverage(Cache(), ["ACGTACGTACGT"], ["a", "b"], [1000, 2000], 3000)

    assert result is not None
