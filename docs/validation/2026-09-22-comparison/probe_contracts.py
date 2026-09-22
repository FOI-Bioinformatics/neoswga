"""Read-only probes for the 22 September comparative review.

Run from the repository with its Python environment:
    python docs/validation/2026-09-22-comparison/probe_contracts.py

These describe observed behavior; they are not regression tests approving it.
No genome processing, external code execution or optimization is performed.
"""

import ast
import json
import subprocess
import tempfile
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from neoswga.core.base_optimizer import PrimerSetMetrics
from neoswga.core.design_request import design_request_for_run, resolve_design_request
from neoswga.core.panel_evaluation import evaluate_panel


def main():
    root = Path(__file__).resolve().parents[3]
    base = dict(
        fg_prefixes=["target"], fg_seq_lengths=[100000], bg_prefixes=[],
        bg_seq_lengths=[], polymerase="phi29", reaction_temp=30.0,
        min_k=12, max_k=12, num_primers=12, data_dir="results",
    )
    request = resolve_design_request(base)
    initial_hash = request.request_hash
    request.conditions.temp = 31.0
    evidence = {
        "commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "conditions_mutation_accepted": request.conditions.temp == 31.0,
        "mutation_changes_request_hash": request.request_hash != initial_hash,
        "same_hash_after_setting_change": {},
    }
    for key, value in (
        ("bg_circular", True), ("stage1_objective_width", 64),
        ("swap_max_evaluations", 123),
    ):
        evidence["same_hash_after_setting_change"][key] = (
            resolve_design_request(base).request_hash
            == resolve_design_request(dict(base, **{key: value})).request_hash
        )
    request = resolve_design_request(
        dict(base, concentration_mode="fixed_total", total_primer_molar=12e-6)
    )
    evidence["concentration"] = {
        "declared_two_oligo_allocation_molar": request.concentrations_molar(
            ("ACGGACGGACGG", "AGGAGGAGGAGG")
        ),
        "reaction_primer_conc_molar": request.conditions.primer_conc,
    }
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "params.json"
        path.write_text(json.dumps(base))
        resolved = design_request_for_run(
            SimpleNamespace(json_file=str(path), reaction_temp=42.0, seed=17), None
        )
        evidence["request_with_cli_overrides"] = {
            "requested_temperature": 42.0, "resolved_temperature": resolved.conditions.temp,
            "requested_seed": 17, "resolved_seed": resolved.seed,
        }
    request = resolve_design_request(base)
    metrics = replace(PrimerSetMetrics.empty(), mean_gap=0.0, max_gap=0.0)
    assessment = evaluate_panel(request, ("AAAAAAAAAAAA",), metrics)
    evidence["assessment_without_objective"] = {
        "qualified": assessment.qualified,
        "requested_size": request.target_size,
        "assessed_size": len(assessment.primers),
        "has_effective_coverage": "effective_fg_coverage" in assessment.metrics,
    }
    evidence["explicit_zero_target_resolves_to"] = resolve_design_request(
        dict(base, target_set_size=0)
    ).target_size
    names = ("evaluate_panel", "concentrations_molar", "require_disjoint_experiments")
    calls = {name: [] for name in names}
    for path in sorted((root / "neoswga").rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Call):
                name = getattr(node.func, "id", getattr(node.func, "attr", None))
                if name in calls:
                    calls[name].append(f"{path.relative_to(root)}:{node.lineno}")
    evidence["named_production_calls"] = calls
    print(json.dumps(evidence, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
