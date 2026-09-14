"""Read-only numerical audit probes; run with PYTHONPATH=. python this_file."""

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

from neoswga.core.additives import AdditiveConcentrations
from neoswga.core.base_optimizer import BaseOptimizer, OptimizerConfig
from neoswga.core.coverage import compute_per_prefix_coverage
from neoswga.core.mechanistic_model import MechanisticModel
from neoswga.core.reaction_conditions import ReactionConditions, build_reaction_conditions
from neoswga.core.string_search import get_cached_genome_sequence


def main():
    primer = "ATCGATCGATCG"
    evidence = {}
    cases = []
    for conc in (0.05e-6, 5e-6):
        conditions = build_reaction_conditions(
            SimpleNamespace(
                polymerase="phi29",
                reaction_temp=30,
                primer_conc=conc,
                na_conc=50,
                mg_conc=10,
            )
        )
        cases.append(
            dict(
                configured_concentration=conc,
                default_call_tm=conditions.calculate_effective_tm(primer),
                explicit_call_tm=conditions.calculate_effective_tm(primer, conc),
            )
        )
    evidence["primer_concentration"] = cases
    a = AdditiveConcentrations(propanediol_m=1)
    c = ReactionConditions.from_additives(a)
    evidence["propanediol_conversion"] = dict(input_m=a.propanediol_m, output_m=c.propanediol_m)
    base = ReactionConditions()
    model = MechanisticModel(base)
    evidence["ignored_tm_modifier"] = dict(
        correction_minus_1=model._calculate_effective_tm(primer, 0.5, -1),
        correction_minus_10=model._calculate_effective_tm(primer, 0.5, -10),
    )
    evidence["dmso_correction_by_temperature"] = {
        temp: AdditiveConcentrations(dmso_percent=5).calculate_tm_correction(
            0.5, 12, reaction_temp_celsius=temp
        )
        for temp in (30, 37, 42)
    }

    class EmptyCache:
        def get_positions(self, *args):
            return []

    evidence["circular_without_sites"] = compute_per_prefix_coverage(
        EmptyCache(), [primer], ["target"], [1000], extension=1000, circular=True
    )[0]
    fake = SimpleNamespace(
        config=OptimizerConfig(extension_reach=100, fg_circular=False),
        conditions=SimpleNamespace(temp=30, calculate_effective_tm=lambda p: 30),
    )
    evidence["occupancy_grouping"] = {
        "one_site": BaseOptimizer._compute_effective_coverage(fake, {primer: [500]}, 1000),
        "two_same_primer_sites": BaseOptimizer._compute_effective_coverage(
            fake, {primer: [500, 510]}, 1000
        ),
        "two_distinct_primers": BaseOptimizer._compute_effective_coverage(
            fake, {primer: [500], "ATCGATCGATCA": [510]}, 1000
        ),
    }
    with TemporaryDirectory() as directory:
        fasta = Path(directory) / "two_contigs.fna"
        fasta.write_text(">one\nAAAAAACCCCCC\n>two\nGGGGGGTTTTTT\n")
        sequence = get_cached_genome_sequence(str(fasta))
        evidence["contig_boundary"] = dict(
            concatenated_sequence=sequence, artificial_join_match=sequence.find("CCCCCCGGGGGG")
        )
    root = Path(__file__).resolve().parents[3]
    saved = root / "examples/wolbachia_pool_design/results/pool_plan.json"
    if saved.exists():
        from neoswga.core.position_cache import PositionCache

        plan = json.loads(saved.read_text())
        row = next(r for r in plan["rows"] if r["requested_size"] == 26)
        prefix = plan["inputs"]["foreground_prefixes"][0]
        if Path(f"{prefix}_12mer_positions.h5").exists():
            cache = PositionCache([prefix], row["primers"], on_missing="error")
            positions = {p: list(cache.get_positions(prefix, p, "both")) for p in row["primers"]}
            conditions = build_reaction_conditions(SimpleNamespace(**plan["design_parameters"]))
            sensitivity = {}
            for reach in (1000, 3000, 5000, 10000):
                subject = SimpleNamespace(
                    config=OptimizerConfig(extension_reach=reach, fg_circular=True),
                    conditions=conditions,
                )
                sensitivity[reach] = BaseOptimizer._compute_effective_coverage(
                    subject, positions, 1267782
                )
            evidence["wmel_same_26_oligos_reach_sensitivity"] = sensitivity
    destination = Path(__file__).with_name("probe_results.json")
    destination.write_text(json.dumps(evidence, indent=2) + "\n")
    print(destination.read_text())


if __name__ == "__main__":
    main()
