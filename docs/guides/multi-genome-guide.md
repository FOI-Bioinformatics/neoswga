# Multi-Genome SWGA Guide

## Overview

NeoSWGA now supports **multi-genome primer selection** with differential penalty weights for:

1. **Target genomes**: Organisms you want to amplify (maximize binding)
2. **Background genomes**: Host DNA you want to avoid but can tolerate (e.g. a human, mouse or arthropod host)
3. **Blacklist genomes**: Organisms you must strongly avoid (e.g., co-existing bacteria, contaminating species)

This enables real-world applications like detecting pathogens in host tissue while avoiding both host DNA and other co-existing organisms.

---

## Key Innovation: Differential Penalties

Traditional SWGA treats all non-target genomes equally. NeoSWGA's multi-genome system uses **differential penalty weights**:

| Genome Type | Penalty Weight | Interpretation |
|-------------|---------------|----------------|
| Target | 0.0x | No penalty (maximize binding) |
| Background | 1.0x | Standard avoidance (tolerable) |
| Blacklist | 5.0-10.0x | Strong avoidance (must minimize) |

This allows you to:
- **Tolerate** some binding to host DNA (which is inevitable in real samples)
- **Strictly avoid** binding to blacklisted organisms (which would be false positives)

---

## Use Cases

### 1. A vector-borne target, with a co-occurring organism blacklisted

**Challenge**: amplify a target carried by an arthropod vector while avoiding:
- Vector DNA (host, far more abundant than the target)
- A second organism that occurs in the same vector, which would read as a false
  positive if its DNA amplified too

**Solution**:
```python
from neoswga.core.multi_genome_pipeline import MultiGenomePipeline, GenomeSet

genome_set = GenomeSet()

# Target: what we want to amplify
genome_set.add_genome(
    name="target_species",
    fasta_path="target.fasta",
    role="target"
)

# Background: the vector's own genome. Tolerate some binding; it is
# unavoidable in a real sample.
genome_set.add_genome(
    name="vector_host",
    fasta_path="vector.fasta",
    role="background",
    penalty_weight=1.0  # Standard avoidance
)

# Blacklist: the co-occurring organism. Avoid completely.
genome_set.add_genome(
    name="cooccurring_species",
    fasta_path="cooccurring.fasta",
    role="blacklist",
    penalty_weight=10.0  # Strong avoidance
)

pipeline = MultiGenomePipeline(genome_set, output_dir="results")
result = pipeline.run()
```

**Result**: primers that bind frequently to the target, less to the vector, and
rarely or never to the blacklisted organism.

---

### 2. Distinguishing one species from its close relatives

**Challenge**: amplify one species while avoiding:
- A large host genome (3 Gbp is common, and extremely abundant)
- Sibling species in the same genus, which would confuse identification

This is the pattern to use whenever the answer has to name a species rather
than a genus. The blacklist carries the relatives; the background carries the
host.

**Solution**:
```python
genome_set = GenomeSet()

genome_set.add_genome(
    name="target_species",
    fasta_path="target.fasta",
    role="target"
)

# Background: the host. A single chromosome is enough while testing.
genome_set.add_genome(
    name="host_chr1",
    fasta_path="host_chr1.fasta",
    role="background",
    penalty_weight=1.0
)

# Blacklist: every sibling species the assay must not call.
for relative in ["relative_a", "relative_b"]:
    genome_set.add_genome(
        name=relative,
        fasta_path=f"{relative}.fasta",
        role="blacklist",
        penalty_weight=10.0
    )

pipeline = MultiGenomePipeline(genome_set, output_dir="results")
result = pipeline.run()
```

**Result**: species-level discrimination rather than genus-level.

---

### 3. A clinical specimen with near relatives present

**Challenge**: amplify a target from a clinical specimen while avoiding:
- Host DNA
- Environmental relatives of the target that are commonly present in the same
  specimen type and would read as false positives

The difference from case 2 is the weighting: these relatives are a nuisance
rather than a competing diagnosis, so they take a lower penalty than a sibling
species would.

**Solution**:
```python
genome_set = GenomeSet()

genome_set.add_genome(
    name="target_species",
    fasta_path="target.fasta",
    role="target"
)

genome_set.add_genome(
    name="host",
    fasta_path="host.fasta",
    role="background"
)

# Blacklist: several related species at once.
for relative in ["relative_a", "relative_b", "relative_c"]:
    genome_set.add_genome(
        name=relative,
        fasta_path=f"{relative}.fasta",
        role="blacklist",
        penalty_weight=8.0
    )

pipeline = MultiGenomePipeline(genome_set, output_dir="results")
result = pipeline.run()
```

---

## API Reference

### GenomeSet

Container for organizing genomes by role.

```python
from neoswga.core.multi_genome_filter import GenomeSet, GenomeRole

genome_set = GenomeSet()

# Add genomes
genome_set.add_genome(
    name="Descriptive_name",
    fasta_path="/path/to/genome.fasta",
    role="target",  # or "background" or "blacklist"
    penalty_weight=1.0  # Optional, defaults based on role
)

# Validate
genome_set.validate()  # Raises error if invalid

# Summary
print(genome_set.summary())
```

**Roles**:
- `"target"`: Maximize binding (penalty = 0.0)
- `"background"`: Minimize but tolerate (penalty = 1.0 default)
- `"blacklist"`: Strongly avoid (penalty = 5.0 default)

The `GenomeRole` enum is accepted wherever a role string is, so
`role=GenomeRole.TARGET` and `role="target"` are equivalent.

### Reading the enrichment figure

`enrichment_score` is the target/background binding-frequency ratio. A primer
with no detectable background binding cannot have a ratio, and is reported as
`MAX_ENRICHMENT` (1e6) rather than infinity -- printed as "no background binding
detected". Infinity would make `results.json` invalid JSON, and it would flatten
ranking, since every perfectly selective primer would tie. With a finite
sentinel the composite score stays proportional to target binding frequency, so
equally selective primers are ordered by how well they bind the target.

### MultiGenomePipeline

Automated pipeline with multi-genome support.

```python
from neoswga.core.multi_genome_pipeline import MultiGenomePipeline

pipeline = MultiGenomePipeline(
    genome_set=genome_set,
    output_dir="results/",
    kmer_range=(8, 12),  # Optional, auto-selected if None
    preferred_polymerase="equiphi29",  # Optional, auto-selected if None
    primer_count=12,  # Target number of primers
    validate_with_simulation=False  # Optional Gillespie validation
)

# Run pipeline
result = pipeline.run(verbose=True)

# Save results
pipeline.save_results(result)
# Creates: primers.fasta, results.json, protocol.txt, summary.txt
```

### MultiGenomePipelineResult

Complete result with multi-genome metrics.

```python
# Result attributes
result.primers  # List of primer sequences
result.primer_count  # Number of primers

# Target genome info
result.target_genome_names  # List of target names
result.mean_target_gc  # Mean GC content
result.genome_classification  # "at_rich", "balanced", or "gc_rich"

# Non-target genomes
result.background_genome_names
result.blacklist_genome_names

# Strategy
result.polymerase  # Selected polymerase
result.reaction_temp  # Temperature (°C)
result.betaine_concentration  # Betaine (M)

# Multi-genome metrics
result.mean_target_frequency  # Average binding in targets
result.mean_background_frequency  # Average binding in backgrounds
result.mean_blacklist_frequency  # Average binding in blacklists
result.mean_enrichment  # Target/background ratio
result.min_enrichment  # Worst-case enrichment

# Protocol
result.protocol  # Ready-to-use experimental protocol
```

---

## Filtering Strategy

### Stage 1: Thermodynamic Filtering

All primers must pass thermodynamic criteria:
- Melting temperature (Tm) appropriate for SWGA
- GC content constraints
- No strong homodimers (ΔG threshold)
- No strong heterodimers between primers
- No stable hairpins

### Stage 2: Multi-Genome Frequency Filtering

Primers are filtered based on binding frequency with differential penalties:

```python
# For each primer:
target_score = target_frequency * 1.0  # Maximize
background_penalty = background_frequency * 1.0  # Standard penalty
blacklist_penalty = blacklist_frequency * 5.0  # High penalty

# Composite score (maximize)
score = target_score / (background_penalty + blacklist_penalty + ε)
```

**Thresholds** (configurable):
- Min target frequency: 1×10⁻⁵ (must bind to target)
- Max background frequency: 1×10⁻⁴ (tolerate some binding)
- Max blacklist frequency: 1×10⁻⁶ (10x stricter than background)
- Min enrichment: 10x (target/background ratio)

### Stage 3: Ranking and Selection

Primers are ranked by composite score and top N selected.

---

## Output Files

### 1. primers.fasta

Standard FASTA format:
```
>Primer_1
ATGCATGC
>Primer_2
GCTAGCTA
...
```

### 2. results.json

Complete metrics in JSON:
```json
{
  "primers": ["ATGCATGC", "GCTAGCTA", ...],
  "primer_count": 12,
  "target_genome_names": ["target_species"],
  "background_genome_names": ["vector_host"],
  "blacklist_genome_names": ["cooccurring_species"],
  "mean_enrichment": 45.2,
  "min_enrichment": 12.3,
  "mean_target_frequency": 2.3e-05,
  "mean_background_frequency": 5.1e-07,
  "mean_blacklist_frequency": 1.2e-08,
  ...
}
```

### 3. protocol.txt

Ready-to-use experimental protocol:
```
================================================================================
MULTI-GENOME SWGA EXPERIMENTAL PROTOCOL
================================================================================

GENOME CONFIGURATION
  Target genome(s): target_species
  Background genome(s): vector_host
  Blacklist genome(s): cooccurring_species

PRIMER SEQUENCES (12 primers)
  1. ATGCATGC
     Target frequency:     2.3e-05
     Background frequency: 5.1e-07
     Blacklist frequency:  1.2e-08
     Enrichment:           45.1x
  ...

REACTION CONDITIONS
  Polymerase: PHI29
  Temperature: 30°C
  Betaine: 0.0M

EXPECTED PERFORMANCE
  Mean enrichment: 45.2x
  Selectivity: Target-specific amplification
```

### 4. summary.txt

Quick summary:
```
Multi-Genome SWGA Summary

Targets: target_species
Backgrounds: vector_host
Blacklists: cooccurring_species

Primers: 12
Polymerase: phi29
Temperature: 30.0°C
Mean enrichment: 45.2x
```

---

## Customizing Penalty Weights

Default penalties:
- Background: 1.0x
- Blacklist: 5.0x

You can customize based on your specific needs:

```python
# Mild blacklist avoidance (3x)
genome_set.add_genome(
    name="Mild_avoid",
    fasta_path="genome.fasta",
    role="blacklist",
    penalty_weight=3.0
)

# Extreme blacklist avoidance (20x)
genome_set.add_genome(
    name="Must_avoid",
    fasta_path="genome.fasta",
    role="blacklist",
    penalty_weight=20.0
)

# Low-penalty background (accept more binding)
genome_set.add_genome(
    name="Tolerant_background",
    fasta_path="genome.fasta",
    role="background",
    penalty_weight=0.5
)

# High-penalty background (minimize binding)
genome_set.add_genome(
    name="Strict_background",
    fasta_path="genome.fasta",
    role="background",
    penalty_weight=2.0
)
```

**Guidelines**:
- **Background**: 0.5-2.0x (tolerable non-target DNA)
- **Blacklist**: 3.0-20.0x (organisms that must be avoided)
- Higher penalty = stronger avoidance

---

## Multiple Target Genomes

You can amplify multiple organisms simultaneously:

```python
genome_set = GenomeSet()

# Multiple targets (amplify all)
genome_set.add_genome("Pathogen_A", "pathA.fasta", "target")
genome_set.add_genome("Pathogen_B", "pathB.fasta", "target")
genome_set.add_genome("Pathogen_C", "pathC.fasta", "target")

# Common background
genome_set.add_genome("Host", "host.fasta", "background")

pipeline = MultiGenomePipeline(genome_set, output_dir="multi_target")
result = pipeline.run()
```

**Result**: Primers that bind to ALL target genomes while avoiding background.

**Use case**: detection at genus rather than species level (any member of a genus)

---

## Performance Considerations

### Where the time goes

The pipeline is dominated by two costs: Jellyfish k-mer counting (linear in
genome size x number of genomes) and per-candidate thermodynamics (linear in the
number of candidates).

Pairwise primer-primer compatibility is deliberately **not** in that list. The
heterodimer check is O(n^2), and it runs on the ranked shortlist rather than the
full candidate pool. That matters more than it sounds: a 15 kb genome yields
around 15,000 candidates after per-primer filtering, which is ~10^8 pairs, and a
5 Mb bacterial genome roughly 10^10. Running the check on the pool made the
command appear to hang on any real input. It is also the wrong place for it --
heterodimer compatibility only matters among primers that will share a reaction.

### Genome Size

- **Small genomes (<5 Mbp)**: Fast (<2 minutes)
- **Bacterial genomes (5-10 Mbp)**: Moderate (2-5 minutes)
- **Large genomes (>100 Mbp)**: Slower (10-30 minutes)
- **Very large (>1 Gbp)**: Consider using chromosome subset or Bloom filters

### Number of Genomes

Processing time scales linearly with number of genomes.

**Recommendations**:
- Targets: 1-3 (more is fine, but diminishing returns)
- Backgrounds: 1-2 (most common hosts)
- Blacklists: 1-5 (most likely confounders)

### K-mer Range

Smaller ranges = faster:
- Testing: `kmer_range=(8, 9)` (very fast, ~30s)
- Production: `kmer_range=(8, 12)` (optimal, ~5 min)
- Comprehensive: `kmer_range=(8, 15)` (slower, ~15 min)

---

## Best Practices

### 1. Start Simple

Begin with one target and one background:
```python
genome_set = GenomeSet()
genome_set.add_genome("Target", "target.fasta", "target")
genome_set.add_genome("Background", "background.fasta", "background")
```

Validate this works before adding blacklists.

### 2. Use Representative Backgrounds

For human samples, you don't need the entire human genome:
```python
# Use a representative chromosome
genome_set.add_genome(
    "Human_chr1",
    "human_chr1.fasta",
    "background"
)
```

Or combine several chromosomes.

### 3. Test with Small K-mer Range First

```python
pipeline = MultiGenomePipeline(
    genome_set,
    kmer_range=(8, 9),  # Fast testing
    primer_count=5
)
```

Once validated, expand to full range.

### 4. Validate Enrichment

Check `result.min_enrichment` - this is your worst-case scenario:
```python
if result.min_enrichment < 10.0:
    print("Warning: Low enrichment, may have false positives")
```

### 5. Inspect Per-Genome Frequencies

Look at the protocol to see binding in each genome:
```
Primer 1: ATGCATGC
  target_species:      2.3e-05  <- Good target binding
  vector_host:         5.1e-07  <- Low background
  cooccurring_species: 1.2e-08  <- Minimal blacklist
```

---

## Troubleshooting

### No primers pass filters

**Problem**: All primers filtered out

**Solutions**:
1. Relax thermodynamic criteria (check GC constraints)
2. Reduce min_enrichment threshold
3. Expand k-mer range
4. Check genome quality (Ns, low complexity)

### Low enrichment

**Problem**: `mean_enrichment < 10x`

**Solutions**:
1. Increase blacklist penalty weights
2. Reduce max_background_frequency
3. Try different k-mer range
4. Verify genomes are sufficiently different

### Too few primers

**Problem**: `primer_count < requested`

**Solutions**:
1. Reduce stringency (thermodynamic thresholds)
2. Expand k-mer range
3. Accept lower enrichment temporarily
4. Check if genomes are too similar

### Background frequency too high

**Problem**: Primers bind too much to host

**Solutions**:
1. Increase background penalty weight
2. Reduce max_background_frequency
3. Use different k-mer range
4. Consider pre-filtering (remove repetitive regions)

---

## Advanced: Custom Filtering Thresholds

```python
from neoswga.core.multi_genome_filter import MultiGenomeFilter

mg_filter = MultiGenomeFilter(
    genome_set=genome_set,
    min_target_freq=1e-5,      # Must bind target (adjust for genome size)
    max_background_freq=1e-4,   # Tolerate some background
    max_blacklist_freq=1e-6,    # Strict blacklist (10x lower)
    min_enrichment=20.0         # Higher enrichment requirement
)
```

**Tradeoffs**:
- Stricter thresholds → fewer primers, higher specificity
- Relaxed thresholds → more primers, may sacrifice specificity

---

## Future Enhancements

Planned features:
1. **MCP server integration**: Automatic download of reference genomes
2. **Batch processing**: Process multiple experiments in parallel
3. **Interactive mode**: Adjust parameters and re-filter in real-time
4. **Wet-lab validation feedback**: Update penalties based on experimental results
5. **Bloom filter optimization**: Handle 100+ Gbp backgrounds efficiently

---

## CLI Usage

The multi-genome pipeline is also available from the command line:

```bash
neoswga multi-genome \
  --genomes target1.fna target2.fna \
  --output results/
```

## Citation

If you use multi-genome SWGA in your research, please cite:

```bibtex
@article{dwivedi2023fast,
  title={A fast machine-learning-guided primer design pipeline for selective whole genome amplification},
  author={Dwivedi-Yu, Jane A and Oppler, Zachary J and Mitchell, Matthew W and Song, Yun S and Brisson, Dustin},
  journal={PLOS Computational Biology},
  volume={19},
  number={4},
  pages={e1010137},
  year={2023}
}
```
