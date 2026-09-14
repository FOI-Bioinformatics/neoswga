"""Read-only pipeline audit; writes the adjacent numerical evidence JSON."""

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from neoswga.core.network_optimizer import NetworkOptimizer
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamic_filter import ThermodynamicCriteria, ThermodynamicFilter

primer = "ATCGATCGATCG"
conditions = ReactionConditions(temp=30, na_conc=50, mg_conc=10, dmso_percent=5, betaine_m=1)
network = SimpleNamespace(conditions=conditions, _tm_cache={})
screen = ThermodynamicFilter(ThermodynamicCriteria(na_conc=50, mg_conc=10))
evidence = {
    "diagnostic_primer": primer,
    "conditions": conditions.to_dict(),
    "tm_celsius": {
        "canonical": conditions.calculate_effective_tm(primer),
        "network": NetworkOptimizer._get_primer_tm(network, primer),
        "secondary_screen": screen.analyze_primer(primer).tm,
    },
}
root = Path(__file__).resolve().parents[3]
work = root / "examples/wolbachia_pool_design/work"
if (work / "step3_df.csv").exists():
    step2 = pd.read_csv(work / "step2_df.csv")["primer"]
    step3 = pd.read_csv(work / "step3_df.csv")["primer"]
    evidence["wmel_candidates"] = {
        "step2": len(step2),
        "step3": len(step3),
        "same_sequences": set(step2) == set(step3),
        "funnel": json.loads((work / "filter_stats.json").read_text()),
    }
Path(__file__).with_name("pipeline_probe_results.json").write_text(
    json.dumps(evidence, indent=2) + "\n"
)
print(json.dumps(evidence, indent=2))
