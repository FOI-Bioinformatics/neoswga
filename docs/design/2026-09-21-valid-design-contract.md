# Valid oligo-pool design contract

Status: proposed requirements, 21 September 2026. This specification supersedes
earlier requirements to preserve legacy command behavior or delivered panels.

## Goal

Design a small oligo pool that satisfies explicit target-coverage and specificity
requirements under a declared model, and report the evidence and search limits
needed to interpret it. Wolbachia wMel against its supplied host reference is the
initial application fixture; the implementation is organism-independent.

## Requirements

1. Backward compatibility is not a requirement. Remove obsolete parameters,
   aliases, adapters and artifact readers when they weaken the design contract.
   Diagnose unsupported artifacts and give explicit regeneration instructions.
2. Missing data, stale provenance, unsupported requested chemistry, non-finite
   calculations and software failures must not be replaced with another model,
   zero binding, an empty index or a successful result.
3. All design commands use one immutable resolved request, one evaluator and one
   final acceptance procedure. No stage changes chemistry or thresholds.
4. Count, QC, candidate preparation and optimization retain distinct contracts.
   Every QC-passing candidate remains available to search. Ranking and bounded
   search can postpone examination but must not silently delete candidates.
5. The count/index coordinate system records reference identity, record lengths,
   strand semantics and per-record circularity. Missing answers differ from
   verified zero occurrences.
6. Chemistry models declare their supported sequence lengths, units, salts,
   temperature ranges, additives, combinations and oligo modifications. Numerical
   support and empirical validation are separate properties. An absent effect
   model must not be represented as a known zero effect.
7. Geometric coverage, occupancy-weighted coverage and measured sequencing breadth
   are different quantities. Use the requested quantity consistently in selection,
   reduction, acceptance and the report. Do not present a proxy as measured or
   validated recovery.
8. Pool concentration is explicit: fixed per-oligo concentrations or a fixed total
   concentration with declared allocation. Re-evaluate affected quantities when
   composition or size changes. Concentration-dependent QC cannot permanently
   discard a candidate using an unrelated provisional concentration.
9. Hard constraints include fixed/excluded oligos, candidate eligibility, pool
   size, configured dimer restrictions, target-specific coverage and configured
   background limits. None can be traded away by a composite score.
10. Report the smallest qualifying pool found and coverage-versus-count tradeoffs.
    Claim a global minimum only with a valid certificate for the full stated
    candidate universe and objective. Exhausting candidates does not exhaust
    all candidate subsets.
11. Budget exhaustion is a recorded search termination, not a model failure or
    proof of infeasibility. Preserve verified incumbents; any interrupted or
    failed run retains its own status even when a diagnostic incumbent exists.
12. Sequencing feedback must verify reference and experiment provenance, report
    observed breadth at declared depth, and separate fitting, selection and
    held-out evaluation. A fitted reach is not automatically a physical extension
    length or an estimate of each oligo's individual effect.
13. Documentation uses modest scientific language. Do not use Unicode in Nextflow
    files. Keep existing module/function size ceilings; extract code instead of
    raising them.

## Supported result claims

| Analysis | Permitted claim | Required evidence |
|---|---|---|
| Geometric design | Fraction within declared binding-site windows | Verified positions, geometry and explicit reach convention |
| Chemistry-weighted design | Occupancy-weighted coverage proxy | Supported model domain, recorded assumptions and resolved chemistry |
| Calibrated prediction | Predicted sequencing breadth for a declared domain and depth | Frozen calibration artifact and held-out evaluation for that domain |
| Sequencing assessment | Observed breadth and uniformity | Verified alignments, explicit depth/mapping filters and reference denominator |

A supported numerical proxy can be useful before prospective validation. Its
result remains a proxy; software correctness does not upgrade its evidence level.
The plan does not assume literature constants are valid for short SWGA oligos or
for combinations of additives. That must be established parameter by parameter.

## Delivery policy

Only a finished run with a qualifying panel can write a recommended oligo pool.
A budget-limited finished search may recommend its verified incumbent, with the
budget limit prominent and no completeness claim. A failed run writes a failure
artifact and may retain a separately marked diagnostic checkpoint. Failure while
performing required final validation prevents recommendation export.

Exploratory geometry-only analysis and analysis without a background are explicit
request types with restricted claims. They are never selected after a chemistry
or background computation fails. Optional presentation failures may leave a valid
machine-readable result intact, but must be reported as rendering failures.
