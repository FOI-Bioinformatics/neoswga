"""A genome must be checked before the pipeline spends an hour on it.

`GenomeLoader.validate_genome` computes exactly the right things -- length, N
fraction, empty sequences, implausible GC -- and its only non-test caller was
inside an `if __name__ == "__main__"` demo block. Nothing on the pipeline path
called it, so whatever a download left behind went straight into jellyfish.

Measured before this check existed, on a set of malformed inputs:

    input                        count-kmers   filter        what the user saw
    HTML 404 page as .fna        fails         --            a C++ crash from
                                                             jellyfish, naming
                                                             neither file nor cause
    header only, 0 bp            PASSES        fails         Step 2 failed: 'primer'
    all N                        PASSES        fails         Step 2 failed: 'primer'
    truncated to a third         PASSES        50 candidates nothing at all

The last row is the dangerous one. A partial download is still valid FASTA, so
it passes every check the pipeline makes and yields a panel designed against a
third of the target, silently.

WHY LENGTH IS NOT A HARD GATE. `validate_genome` defaults `min_length` to
100,000 bp, and the shipped plasmid example is 6,157 bp. Plasmids, viruses and
organelles are legitimate targets, so refusing a short genome would reject real
work to catch a truncation the tool cannot distinguish from it. The rule adopted
here is narrower: a genome that carries NO sequence at all is refused, because
nothing legitimate is zero bases, and everything else is reported -- the
measured length is logged for every input, so a truncation is visible to the
person who knows what the genome should weigh.
"""

import pytest

from neoswga.core.genome_io import GenomeLoader


def _write(tmp_path, name, body):
    p = tmp_path / name
    p.write_text(body)
    return str(p)


def _genome(tmp_path, name, bases, line=70):
    seq = (bases * ((60_000 // max(len(bases), 1)) + 1))[:60_000]
    return _write(
        tmp_path,
        name,
        ">chr1\n" + "\n".join(seq[i : i + line] for i in range(0, len(seq), line)) + "\n",
    )


# ---------------------------------------------------------------------------
# What the check refuses
# ---------------------------------------------------------------------------


def test_a_genome_with_no_sequence_is_refused(tmp_path):
    """The truncated-download case that used to reach jellyfish."""
    from neoswga.core.pipeline import check_genome_inputs

    path = _write(tmp_path, "header_only.fna", ">chr1 truncated\n")
    problems = check_genome_inputs([path])

    assert problems, "a 0 bp genome was accepted"
    assert "header_only.fna" in problems[0], f"the file is not named: {problems[0]}"
    assert "0" in problems[0]


def test_a_missing_file_is_refused_by_name(tmp_path):
    from neoswga.core.pipeline import check_genome_inputs

    problems = check_genome_inputs([str(tmp_path / "absent.fna")])
    assert problems and "absent.fna" in problems[0]


def test_a_non_fasta_file_is_refused(tmp_path):
    """An HTML error page saved with a genome's name is the classic artifact of
    a download that failed without anyone noticing."""
    from neoswga.core.pipeline import check_genome_inputs

    path = _write(
        tmp_path,
        "html_404.fna",
        "<!DOCTYPE html>\n<html><head><title>404 Not Found</title></head></html>\n",
    )
    problems = check_genome_inputs([path])

    assert problems, "an HTML error page was accepted as a genome"
    assert "html_404.fna" in problems[0]


# ---------------------------------------------------------------------------
# What the check must NOT refuse
# ---------------------------------------------------------------------------


def test_a_plasmid_sized_genome_is_accepted(tmp_path):
    """`validate_genome` calls anything under 100 kb "too short". The shipped
    plasmid example is 6,157 bp, so length cannot be a hard gate without
    rejecting the tool's own example and every plasmid design."""
    from neoswga.core.pipeline import check_genome_inputs

    seq = "ACGTACGTAC" * 600  # 6,000 bp
    path = _write(tmp_path, "plasmid.fna", ">p\n" + seq + "\n")

    assert check_genome_inputs([path]) == []


def test_an_ordinary_genome_is_accepted(tmp_path):
    from neoswga.core.pipeline import check_genome_inputs

    assert check_genome_inputs([_genome(tmp_path, "good.fna", "ACGTGGATCCAT")]) == []


def test_a_partly_masked_genome_is_reported_but_not_refused(tmp_path):
    """Real assemblies carry N runs. Worth saying, not worth refusing."""
    from neoswga.core.pipeline import check_genome_inputs

    # Half N, half sequence: gappy, but primers can still bind.
    path = _genome(tmp_path, "gappy.fna", "ACGTNNNNNN")
    assert check_genome_inputs([path]) == []


def test_an_all_n_genome_is_refused(tmp_path):
    """Not one A, C, G or T is as degenerate as no bases at all: no primer can
    bind anywhere. Refusing it here replaces a `Step 2 failed: 'primer'` two
    steps downstream."""
    from neoswga.core.pipeline import check_genome_inputs

    problems = check_genome_inputs([_genome(tmp_path, "all_n.fna", "N")])

    assert problems, "a genome of pure N was accepted"
    assert "all_n.fna" in problems[0]
    assert "N" in problems[0]


def test_an_empty_list_is_not_an_error(tmp_path):
    """A run with no background genomes is normal."""
    from neoswga.core.pipeline import check_genome_inputs

    assert check_genome_inputs([]) == []
    assert check_genome_inputs(None) == []


# ---------------------------------------------------------------------------
# What it reports
# ---------------------------------------------------------------------------


def test_every_genome_has_its_measured_length_logged(tmp_path, caplog):
    """The signal that catches a truncation without false-positiving on a
    plasmid: state the size and let the person who knows the organism judge it."""
    import logging

    from neoswga.core.pipeline import check_genome_inputs

    path = _genome(tmp_path, "good.fna", "ACGTGGATCCAT")
    with caplog.at_level(logging.INFO):
        check_genome_inputs([path])

    logged = caplog.text
    assert "good.fna" in logged
    assert "60,000" in logged or "60000" in logged, (
        f"the measured length is not reported; a truncated genome would be "
        f"invisible. Logged: {logged!r}"
    )


# ---------------------------------------------------------------------------
# The raise path, which the helper tests above do not reach
# ---------------------------------------------------------------------------


def test_step1_refuses_to_start_on_a_bad_genome(tmp_path, monkeypatch):
    """Exercises the `raise`, not just the check.

    `StepPrerequisiteError` takes a `StepValidationResult`, not a string, and
    its `__init__` reads `validation.error_message`. Passing a bare string
    raises `AttributeError` from inside the error handler -- a failure while
    reporting a failure. The helper tests above all call `check_genome_inputs`
    directly and would not have noticed.
    """
    import neoswga.core.pipeline as pipeline_mod
    from neoswga.core.pipeline import StepPrerequisiteError

    path = _write(tmp_path, "header_only.fna", ">chr1 truncated\n")

    monkeypatch.setattr(pipeline_mod, "_initialize", lambda: None)
    monkeypatch.setattr(pipeline_mod, "fg_genomes", [path], raising=False)
    monkeypatch.setattr(pipeline_mod, "bg_genomes", [], raising=False)
    monkeypatch.setattr(pipeline_mod, "fg_prefixes", ["fg"], raising=False)
    monkeypatch.setattr(pipeline_mod, "bg_prefixes", [], raising=False)

    with pytest.raises(StepPrerequisiteError) as excinfo:
        pipeline_mod.step1()

    message = str(excinfo.value)
    assert "header_only.fna" in message, f"the offending file is not named: {message}"
    assert "download" in message.lower(), "the remediation does not mention the likely cause"
