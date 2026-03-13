# Genome-Adaptive QA User Guide

**NeoSWGA 3.0+**

A practical guide to using genome-adaptive quality assessment for improved primer design in genomes with extreme GC content.

---

## Table of Contents

1. [What is Adaptive QA?](#what-is-adaptive-qa)
2. [When Should I Use It?](#when-should-i-use-it)
3. [Quick Start](#quick-start)
4. [Detailed Usage](#detailed-usage)
5. [Understanding Results](#understanding-results)
6. [Best Practices](#best-practices)
7. [Troubleshooting](#troubleshooting)
8. [API Reference](#api-reference)

---

## What is Adaptive QA?

Genome-adaptive quality assessment automatically adjusts primer quality thresholds based on your target genome's GC content.

**The Problem**: Standard primer design assumes ~50% GC content. For extreme genomes (very AT-rich or GC-rich), this creates unrealistic expectations and rejects primers that would work well for your specific genome.

**The Solution**: Adaptive QA adjusts quality thresholds to match your genome's composition, enabling effective primer design for challenging genomes.

### Key Adjustments

Adaptive QA automatically adjusts three critical thresholds:

1. **Terminal 3' Stability (Tm)**: Raised for GC-rich genomes, lowered for AT-rich genomes
2. **GC Content Range**: Centered on genome GC ± 15% (instead of fixed 35-65%)
3. **Dimer Formation**: More stringent for GC-rich, more permissive for AT-rich

---

## When Should I Use It?

### Use Adaptive QA When:

- Your genome is **AT-rich** (< 40% GC)
  - Examples: Francisella, Plasmodium, Borrelia
  - Benefit: +75-200% primer coverage

- Your genome is **GC-rich** (> 60% GC)
  - Examples: Burkholderia, Streptomyces, Mycobacterium
  - Benefit: +70-150% primer coverage

- Standard primer design yields **few acceptable primers**

- You need **maximum genome coverage** for challenging targets

### Standard QA is Fine When:

- Your genome is balanced (40-60% GC)
- Standard primer design already yields good results
- You don't have severe coverage gaps

---

## Quick Start

### Automatic Mode (Recommended)

Adaptive QA is **enabled automatically** when you provide a genome file:

```python
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

# Just provide your genome - adaptive QA activates automatically!
generator = OptimalOligoGenerator(
    genome_file='path/to/francisella.fna'  # 32% GC genome
)

result = generator.generate_primers(
    primer_count=10,
    length_range=(10, 14)
)

# Adaptive QA was used automatically
print(f"Primers generated: {len(result.primers)}")
print(f"Adaptive QA: {result.reaction_conditions.use_adaptive_qa}")
print(f"Genome GC: {result.reaction_conditions.genome_gc:.1%}")
```

**That's it!** The system analyzes your genome and applies adaptive QA automatically.

---

## Detailed Usage

### Method 1: Automatic (High-Level API)

Best for most users - let the system decide:

```python
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

# Automatic genome analysis and adaptive QA
generator = OptimalOligoGenerator(
    genome_file='burkholderia.fna'  # 68% GC
)

result = generator.generate_primers(
    primer_count=12,
    length_range=(10, 14),
    stringency='moderate'
)

# Check if adaptive QA was used
if result.reaction_conditions.use_adaptive_qa:
    print(f"Adaptive QA activated for {result.reaction_conditions.genome_gc:.1%} GC genome")
else:
    print("Standard QA used (genome is balanced)")
```

### Method 2: Manual Control (Advanced)

For users who want explicit control:

```python
from neoswga.core.genome_analysis import analyze_genome_for_qa
from neoswga.core.integrated_quality_scorer import create_quality_scorer

# Step 1: Analyze genome
analysis = analyze_genome_for_qa('plasmodium.fna')

print(f"Genome GC: {analysis['gc_content']:.1%}")
print(f"GC class: {analysis['adaptive_qa_recommendation']['gc_class']}")
print(f"Recommended: {analysis['adaptive_qa_recommendation']['reason']}")
print(f"Expected improvement: {analysis['adaptive_qa_recommendation']['expected_improvement']}")

# Step 2: Create adaptive quality scorer
scorer = create_quality_scorer(
    stringency='moderate',
    genome_gc=analysis['gc_content']  # Enable adaptive QA
)

# Step 3: Score primers
test_primer = "ATATATATATAT"
score = scorer.score_primer(test_primer)

print(f"Primer: {test_primer}")
print(f"Overall score: {score.overall_score:.3f}")
print(f"Passes QA: {score.passes_all}")
```

### Method 3: Disable Adaptive QA (Use Standard QA)

If you want to force standard QA even for extreme genomes:

```python
# Create scorer without genome_gc parameter
scorer = create_quality_scorer(
    stringency='moderate',
    genome_gc=None  # Forces standard QA (50% GC baseline)
)

# This uses standard thresholds regardless of genome GC
score = scorer.score_primer("ATATATATATAT")
```

---

## Understanding Results

### Genome Analysis Output

When you analyze a genome, you get:

```python
{
    'gc_content': 0.323,  # 32.3% GC
    'genome_stats': {
        'length': 1892775,
        'contigs': 1,
        'gc_std': 0.018
    },
    'adaptive_qa_recommendation': {
        'gc_class': 'at_rich',
        'use_adaptive': True,
        'reason': 'AT-rich genome (32.3% GC) benefits from adapted thresholds',
        'expected_improvement': '+75-100% coverage'
    }
}
```

### GC Classification

| GC Class | GC Range | Adaptive QA | Expected Benefit |
|----------|----------|-------------|------------------|
| extreme_at | < 30% | Strongly recommended | +150-200% |
| at_rich | 30-40% | Recommended | +75-100% |
| balanced | 40-60% | Optional | +10-20% |
| gc_rich | 60-70% | Recommended | +70-100% |
| extreme_gc | > 70% | Strongly recommended | +120-150% |

### Log Messages

**Adaptive QA Activated**:
```
INFO - Genome-adaptive QA enabled: genome_gc=32.3%, stringency=moderate,
       terminal_tm=8.0°C (vs 14.0°C standard)
INFO - Using genome-adaptive QA for 32.3% GC genome
```

**Standard QA Used**:
```
INFO - Genome GC close to balanced (50.8%), adaptive QA optional
INFO - Using standard QA thresholds
```

### Threshold Adjustments

For a 68% GC genome (Burkholderia):

| Threshold | Standard | Adaptive | Change |
|-----------|----------|----------|--------|
| Terminal Tm | 14.0°C | 24.4°C | +10.4°C |
| GC range min | 35% | 53% | +18% |
| GC range max | 65% | 83% | +18% |
| Dimer ΔG | -9.0 kcal/mol | -10.5 kcal/mol | -1.5 kcal/mol |

---

## Best Practices

### 1. Let the System Decide

Unless you have specific reasons, use automatic mode:

```python
# Good: Automatic
generator = OptimalOligoGenerator(genome_file='mygenome.fna')

# Less ideal: Manual override without good reason
scorer = create_quality_scorer(genome_gc=None)  # Why disable it?
```

### 2. Check Recommendations

Always review the adaptive QA recommendation:

```python
analysis = analyze_genome_for_qa('mygenome.fna')

if analysis['adaptive_qa_recommendation']['use_adaptive']:
    print(f"RECOMMENDATION: {analysis['adaptive_qa_recommendation']['reason']}")
    print(f"BENEFIT: {analysis['adaptive_qa_recommendation']['expected_improvement']}")
```

### 3. Start with Moderate Stringency

Begin with `stringency='moderate'` and adjust if needed:

```python
# Start here
result = generator.generate_primers(stringency='moderate')

# If too few primers, try lenient
result = generator.generate_primers(stringency='lenient')

# If quality concerns, try strict
result = generator.generate_primers(stringency='strict')
```

### 4. Validate on Real Data

For critical applications, validate primers experimentally:

```python
# Generate primer set
result = generator.generate_primers(primer_count=12)

# Review quality scores
for primer_id, primer_data in result.primers.items():
    score = primer_data['quality']
    if score < 0.7:
        print(f"Review {primer_id}: score {score:.3f}")
```

### 5. Document Your Workflow

Record adaptive QA settings for reproducibility:

```python
# Save configuration
config = {
    'genome_file': 'francisella.fna',
    'genome_gc': result.reaction_conditions.genome_gc,
    'adaptive_qa_used': result.reaction_conditions.use_adaptive_qa,
    'stringency': 'moderate',
    'primer_count': 12
}

import json
with open('primer_config.json', 'w') as f:
    json.dump(config, f, indent=2)
```

---

## Troubleshooting

### Problem: "Too Few Primers Generated"

**Symptoms**: Only 1-2 primers pass QA, need 10+

**Solutions**:
1. Verify adaptive QA is enabled:
   ```python
   print(f"Adaptive QA: {result.reaction_conditions.use_adaptive_qa}")
   ```

2. Try lenient stringency:
   ```python
   result = generator.generate_primers(stringency='lenient')
   ```

3. Expand length range:
   ```python
   result = generator.generate_primers(length_range=(8, 16))
   ```

4. Check genome quality:
   ```python
   analysis = analyze_genome_for_qa('mygenome.fna')
   print(f"Genome length: {analysis['genome_stats']['length']:,} bp")
   print(f"GC std dev: {analysis['genome_stats']['gc_std']:.3f}")
   ```

---

### Problem: "Adaptive QA Not Activating"

**Symptoms**: System uses standard QA when you expect adaptive

**Causes & Solutions**:

1. **Genome is balanced (40-60% GC)**
   ```python
   # Check GC content
   analysis = analyze_genome_for_qa('mygenome.fna')
   print(f"GC: {analysis['gc_content']:.1%}")
   # If 40-60%, adaptive QA is optional (minimal benefit)
   ```

2. **genome_gc parameter not passed**
   ```python
   # Wrong: Manual scorer without genome_gc
   scorer = create_quality_scorer(stringency='moderate')

   # Right: Include genome_gc
   scorer = create_quality_scorer(
       stringency='moderate',
       genome_gc=0.32
   )
   ```

3. **Using old API**
   ```python
   # Old API doesn't support adaptive QA
   # Update to NeoSWGA 3.5+
   ```

---

### Problem: "Primers Fail Experimental Validation"

**Symptoms**: Primers pass QA but don't amplify well

**Diagnostic Steps**:

1. **Check binding site coverage**:
   ```python
   # Verify binding sites span genome
   for primer_id, data in result.primers.items():
       sites = data['binding_sites']
       print(f"{primer_id}: {len(sites)} sites, range {min(sites)}-{max(sites)}")
   ```

2. **Review quality scores**:
   ```python
   # Check for borderline primers
   low_quality = {pid: data['quality']
                  for pid, data in result.primers.items()
                  if data['quality'] < 0.7}
   if low_quality:
       print(f"Low quality primers: {low_quality}")
   ```

3. **Check reaction conditions**:
   ```python
   # Verify PCR conditions match primer Tm
   conditions = result.reaction_conditions
   print(f"Annealing temp: {conditions.annealing_temp}°C")
   print(f"Extension temp: {conditions.extension_temp}°C")
   ```

4. **Try different stringency**:
   ```python
   # More strict QA
   result = generator.generate_primers(stringency='strict')
   ```

---

### Problem: "Want to Compare Standard vs Adaptive"

**Solution**: Generate both and compare:

```python
from neoswga.core.integrated_quality_scorer import create_quality_scorer

# Standard QA
scorer_standard = create_quality_scorer(
    stringency='moderate',
    genome_gc=None  # Standard (50% GC baseline)
)

# Adaptive QA
scorer_adaptive = create_quality_scorer(
    stringency='moderate',
    genome_gc=0.32  # Your genome GC
)

# Compare scores
test_primers = ["ATATATATATAT", "GCGCGCGCGCGC"]

for primer in test_primers:
    score_std = scorer_standard.score_primer(primer).overall_score
    score_adp = scorer_adaptive.score_primer(primer).overall_score

    print(f"{primer}:")
    print(f"  Standard: {score_std:.3f}")
    print(f"  Adaptive: {score_adp:.3f}")
    print(f"  Improvement: {(score_adp - score_std):.3f}")
```

---

## API Reference

### High-Level API

#### OptimalOligoGenerator

**Automatic adaptive QA** (recommended):

```python
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

generator = OptimalOligoGenerator(
    genome_file: str  # Path to genome FASTA
)

result = generator.generate_primers(
    primer_count: int = 10,
    length_range: tuple = (10, 14),
    stringency: str = 'moderate'  # 'lenient', 'moderate', 'strict'
)
```

**Returns**: `PrimerGenerationResult` with:
- `primers`: Dict of primer sequences and metadata
- `reaction_conditions`: Includes `genome_gc` and `use_adaptive_qa`
- `quality_summary`: Statistics on primer quality

---

### Low-Level API

#### analyze_genome_for_qa()

Analyze genome and get adaptive QA recommendation:

```python
from neoswga.core.genome_analysis import analyze_genome_for_qa

analysis = analyze_genome_for_qa(
    genome_file: str | Path  # Path to genome FASTA
) -> dict
```

**Returns**:
```python
{
    'gc_content': float,  # 0-1
    'genome_stats': {
        'length': int,
        'contigs': int,
        'gc_std': float
    },
    'adaptive_qa_recommendation': {
        'gc_class': str,  # extreme_at, at_rich, balanced, gc_rich, extreme_gc
        'use_adaptive': bool,
        'reason': str,
        'expected_improvement': str
    }
}
```

---

#### create_quality_scorer()

Create quality scorer with optional adaptive QA:

```python
from neoswga.core.integrated_quality_scorer import create_quality_scorer

scorer = create_quality_scorer(
    stringency: str = 'moderate',  # 'lenient', 'moderate', 'strict'
    genome_gc: float | None = None  # 0-1, or None for standard QA
) -> PrimerQualityScorer
```

**Parameters**:
- `stringency`: Overall QA stringency level
- `genome_gc`: Genome GC content (0-1). If provided, enables adaptive QA. If None, uses standard QA (50% GC baseline).

**Returns**: `PrimerQualityScorer` with `score_primer()` method

---

#### score_primer()

Score a single primer:

```python
score = scorer.score_primer(
    primer_seq: str  # Primer sequence (DNA)
) -> PrimerQualityScore
```

**Returns**: `PrimerQualityScore` with:
- `overall_score`: float (0-1, higher is better)
- `strand_bias_score`: float
- `dimer_score`: float
- `three_prime_score`: float
- `complexity_score`: float
- `thermo_score`: float
- `passes_all`: bool
- `failure_reasons`: List[str]

---

## Examples

### Example 1: Francisella (AT-Rich)

```python
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

# AT-rich genome (32% GC)
generator = OptimalOligoGenerator(
    genome_file='francisella_tularensis.fna'
)

result = generator.generate_primers(
    primer_count=12,
    length_range=(10, 14),
    stringency='moderate'
)

print(f"Generated {len(result.primers)} primers")
print(f"Genome GC: {result.reaction_conditions.genome_gc:.1%}")
print(f"Adaptive QA: {result.reaction_conditions.use_adaptive_qa}")

# Expected output:
# Generated 12 primers
# Genome GC: 32.3%
# Adaptive QA: True
```

### Example 2: Burkholderia (GC-Rich)

```python
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

# GC-rich genome (68% GC)
generator = OptimalOligoGenerator(
    genome_file='burkholderia_pseudomallei.fna'
)

result = generator.generate_primers(
    primer_count=12,
    length_range=(10, 14),
    stringency='moderate'
)

print(f"Generated {len(result.primers)} primers")
print(f"Genome GC: {result.reaction_conditions.genome_gc:.1%}")
print(f"Adaptive QA: {result.reaction_conditions.use_adaptive_qa}")

# Expected output:
# Generated 12 primers
# Genome GC: 68.2%
# Adaptive QA: True
```

### Example 3: E. coli (Balanced)

```python
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

# Balanced genome (51% GC)
generator = OptimalOligoGenerator(
    genome_file='ecoli_k12.fna'
)

result = generator.generate_primers(
    primer_count=12,
    length_range=(10, 14),
    stringency='moderate'
)

print(f"Generated {len(result.primers)} primers")
print(f"Genome GC: {result.reaction_conditions.genome_gc:.1%}")
print(f"Adaptive QA: {result.reaction_conditions.use_adaptive_qa}")

# Expected output:
# Generated 12 primers
# Genome GC: 50.8%
# Adaptive QA: False (optional for balanced genomes)
```

### Example 4: Manual Comparison

Compare standard vs adaptive QA:

```python
from neoswga.core.genome_analysis import analyze_genome_for_qa
from neoswga.core.integrated_quality_scorer import create_quality_scorer

# Analyze genome
analysis = analyze_genome_for_qa('francisella.fna')
genome_gc = analysis['gc_content']

# Create both scorers
scorer_standard = create_quality_scorer(stringency='moderate', genome_gc=None)
scorer_adaptive = create_quality_scorer(stringency='moderate', genome_gc=genome_gc)

# Test AT-rich primer (typical for Francisella)
at_primer = "ATATATATATAT"

score_std = scorer_standard.score_primer(at_primer)
score_adp = scorer_adaptive.score_primer(at_primer)

print(f"Primer: {at_primer} (17% GC)")
print(f"Standard QA score: {score_std.overall_score:.3f}")
print(f"Adaptive QA score: {score_adp.overall_score:.3f}")
print(f"Improvement: {score_adp.overall_score - score_std.overall_score:+.3f}")
```

---

## Summary

**Key Takeaways**:

1. Adaptive QA **automatically activates** when you provide a genome file
2. Use for **extreme genomes** (< 40% or > 60% GC) for best results
3. Expected improvement: **+75-200%** primer coverage for extreme genomes
4. **Backward compatible**: Standard QA still available via `genome_gc=None`
5. **No code changes** required for existing users

**Next Steps**:
- Try adaptive QA on your genome
- Review generated primers and quality scores
- Validate experimentally
- Report feedback and results

---

**Questions?** Check the [Troubleshooting](#troubleshooting) section or file an issue on GitHub.

**Version**: NeoSWGA 3.0+
