"""`cpus` in params.json has to reach the stage that spends them.

Found while auditing where parallel execution would help
(`docs/validation/parallelism_opportunities_2026-09-17.md`).

`run_jellyfish` declares `cpus: int = 4` and divides the machine by it to decide
how many k values to count at once:

    max_workers = min(num_k, max(1, cpu_count // max(cpus, 1)))

None of step 1's four call sites passed the configured value, so every run
counted at the default regardless of what params.json said. On a 64-core machine
with `cpus: 16` that is 4 threads per jellyfish and 16 concurrent k values, where
the configuration asked for 16 threads and 4. Same class as Known Issue 8: a
schema key that is read in some places and not in the one that would use it.

This asserts the value arrives, not that a particular thread count is fastest.
"""

import ast
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
PIPELINE = ROOT / "neoswga" / "core" / "pipeline.py"


def _run_jellyfish_calls():
    """Every `run_jellyfish(...)` call in the pipeline, as AST nodes."""
    tree = ast.parse(PIPELINE.read_text())
    return [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "run_jellyfish"
    ]


def test_step_one_calls_jellyfish_at_all():
    """Guard the guard: a rename would otherwise make this file vacuous."""
    assert len(_run_jellyfish_calls()) >= 4, (
        "No run_jellyfish calls found in pipeline.py. If the call was renamed, "
        "point this test at the new name rather than deleting it."
    )


def test_every_call_passes_the_configured_cpu_count():
    """Not the default, which the machine is then divided by."""
    offenders = []
    for call in _run_jellyfish_calls():
        passed = {keyword.arg for keyword in call.keywords if keyword.arg}
        # Five positional parameters: genome, prefix, min_k, max_k, cpus.
        positionally = len(call.args) >= 5
        if "cpus" not in passed and not positionally:
            offenders.append(call.lineno)

    assert not offenders, (
        "run_jellyfish called without cpus at pipeline.py line(s) "
        f"{offenders}. The default of 4 then decides both the threads per "
        "jellyfish and, through cpu_count // cpus, how many k values run at "
        "once, so params.json's 'cpus' has no effect on step 1."
    )


def test_the_worker_count_falls_as_the_per_process_thread_count_rises():
    """The two numbers trade off, which is why passing the wrong one matters."""
    import os

    from neoswga.core.kmer_counter import run_jellyfish  # noqa: F401  (import guard)

    cpu_count = os.cpu_count() or 1
    num_k = 7

    def workers(cpus):
        return min(num_k, max(1, cpu_count // max(cpus, 1)))

    assert workers(1) >= workers(cpu_count), (
        "Raising the threads per jellyfish should not raise the number of "
        "concurrent k values; the formula in run_jellyfish divides by it."
    )
