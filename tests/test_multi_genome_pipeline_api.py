"""Guard: MultiGenomePipeline must match the MultiGenomeKmerCounter API.

`multi_genome_pipeline` had drifted away from `kmer_counter` at four call sites and
nothing caught it, because the existing multi-genome tests only inspect the
`MultiGenomePipelineResult` dataclass fields and never construct the pipeline. The
drift:

  * ``MultiGenomeKmerCounter(use_parallel=True)`` -- the real signature is
    ``(cpus, output_dir)``.
  * ``self.kmer_counter.genome_sequences`` -- the attribute is ``genome_fastas``.
  * ``add_genome(name, sequence)`` -- it takes a FASTA *path*, not a sequence.
  * ``count_candidates_all_genomes(...)`` -- no such method existed.

These tests are deliberately import/introspection based so they run without
Jellyfish; the end-to-end counting path is covered by the integration suite.
"""

import inspect

import pytest


def test_kmer_counter_constructor_signature_is_what_the_pipeline_uses():
    from neoswga.core.kmer_counter import MultiGenomeKmerCounter

    params = inspect.signature(MultiGenomeKmerCounter.__init__).parameters
    assert "cpus" in params
    assert "output_dir" in params
    assert "use_parallel" not in params, (
        "if MultiGenomeKmerCounter gained use_parallel, update multi_genome_pipeline"
    )


def test_pipeline_constructs_counter_with_real_kwargs():
    import neoswga.core.multi_genome_pipeline as mgp

    src = inspect.getsource(mgp.MultiGenomePipeline.__init__)
    assert "use_parallel" not in src, "pipeline still passes the non-existent use_parallel kwarg"
    assert "cpus=" in src and "output_dir=" in src


def test_pipeline_uses_attributes_that_exist_on_the_counter():
    from neoswga.core.kmer_counter import MultiGenomeKmerCounter
    import neoswga.core.multi_genome_pipeline as mgp

    counter = MultiGenomeKmerCounter.__init__
    src_counter = inspect.getsource(counter)
    assert "self.genome_fastas" in src_counter
    assert "self.genome_lengths" in src_counter

    src_pipeline = inspect.getsource(mgp)
    assert "kmer_counter.genome_sequences" not in src_pipeline, (
        "pipeline references genome_sequences, which does not exist on the counter"
    )


def test_pipeline_only_calls_methods_the_counter_defines():
    from neoswga.core.kmer_counter import MultiGenomeKmerCounter
    import neoswga.core.multi_genome_pipeline as mgp
    import re

    # Class-level members (methods, properties) ...
    available = {name for name, _ in inspect.getmembers(MultiGenomeKmerCounter)}
    # ... plus instance attributes, which only appear as `self.X =` in __init__.
    available |= set(
        re.findall(
            r"self\.([A-Za-z_][A-Za-z0-9_]*)\s*[:=]",
            inspect.getsource(MultiGenomeKmerCounter.__init__),
        )
    )

    src = inspect.getsource(mgp)
    called = set(re.findall(r"kmer_counter\.([A-Za-z_][A-Za-z0-9_]*)", src))
    unknown = {c for c in called if c not in available}
    assert not unknown, (
        f"multi_genome_pipeline calls kmer_counter.{{{', '.join(sorted(unknown))}}} "
        f"which MultiGenomeKmerCounter does not define"
    )


def test_add_genome_takes_a_path_not_a_sequence(tmp_path):
    """The pipeline used to hand add_genome a loaded sequence string."""
    from neoswga.core.kmer_counter import MultiGenomeKmerCounter

    params = list(inspect.signature(MultiGenomeKmerCounter.add_genome).parameters)
    assert params[1:] == ["name", "fasta_path"]


def test_count_candidates_helper_exists_and_is_used():
    import neoswga.core.multi_genome_pipeline as mgp

    assert hasattr(MGP := mgp.MultiGenomePipeline, "_count_candidates_all_genomes")
    src = inspect.getsource(MGP._apply_multi_genome_filter) if hasattr(
        MGP, "_apply_multi_genome_filter"
    ) else inspect.getsource(mgp)
    assert "count_candidates_all_genomes" in src


def test_count_candidates_groups_by_length_and_canonicalises():
    """Counts must be looked up under the canonical (Jellyfish) form."""
    import neoswga.core.multi_genome_pipeline as mgp

    src = inspect.getsource(mgp.MultiGenomePipeline._count_candidates_all_genomes)
    assert "reverse_complement" in src, "candidate lookup must canonicalise like Jellyfish"
    assert "count_kmers(" in src


def test_pipeline_module_imports_cleanly():
    """A NameError-level regression would surface here."""
    import importlib

    import neoswga.core.multi_genome_pipeline as mgp

    importlib.reload(mgp)
    assert mgp.MultiGenomePipeline is not None
