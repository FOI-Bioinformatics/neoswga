# NeoSWGA Quick Start Guide

Get started with NeoSWGA primer design in minutes.

## Prerequisites

- Python 3.11 or later
- Jellyfish k-mer counter (v2.0+)
- 8+ GB RAM

## Installation

```bash
# Clone and install
git clone https://github.com/FOI-Bioinformatics/neoswga.git
cd neoswga
pip install -e .

# Verify installation
neoswga validate --quick
```

## Your First Primer Design

### 1. Prepare Input Files

You need:
- **Target genome**: FASTA file of the organism you want to amplify
- **Background genome**: FASTA file of contaminating DNA to avoid

### 2. Create Configuration

Create `params.json`:

```json
{
  "fg_genomes": ["target.fasta"],
  "bg_genomes": ["background.fasta"],
  "fg_prefixes": ["data/target"],
  "bg_prefixes": ["data/background"],
  "data_dir": "results/",
  "min_k": 6,
  "max_k": 12,
  "num_primers": 6
}
```

### 3. Run the Pipeline

The simplest approach is to run all four steps with a single command:

```bash
neoswga design -j params.json
```

Or run each step individually for more control:

```bash
# Step 1: Count k-mers (5-30 min depending on genome size)
neoswga count-kmers -j params.json

# Step 2: Filter candidates (1-5 min)
neoswga filter -j params.json

# Step 3: Score primers (1-5 min)
neoswga score -j params.json

# Step 4: Optimize selection (1-10 min)
neoswga optimize -j params.json
```

### 4. Get Results

Your optimized primers are in `results/step4_improved_df.csv`:

```csv
set_id,primers,coverage,enrichment
1,"ATCGATCG,GCTAGCTA,TGCATGCA,...",0.85,125.4
```

### 5. Export and Interpret

```bash
# Interpret quality and predicted enrichment
neoswga interpret -d results/

# FASTA for primer ordering
neoswga export -d results/ --format fasta -o primers.fasta

# BED file for genome browser visualization
neoswga export -d results/ --format bed -o primers.bed

# BedGraph for coverage depth visualization
neoswga export -d results/ --format bedgraph -o primers.bedgraph
```

## Common Workflows

### Standard Phi29 SWGA (6-12mer primers)

```json
{
  "polymerase": "phi29",
  "reaction_temp": 30.0,
  "min_k": 6,
  "max_k": 12,
  "min_tm": 15,
  "max_tm": 45
}
```

### EquiPhi29 with Longer Primers (12-18mer)

```json
{
  "polymerase": "equiphi29",
  "reaction_temp": 42.0,
  "dmso_percent": 5.0,
  "betaine_m": 1.0,
  "min_k": 12,
  "max_k": 18,
  "min_tm": 30,
  "max_tm": 65
}
```

### High-GC Genome

```json
{
  "polymerase": "equiphi29",
  "reaction_temp": 45.0,
  "dmso_percent": 10.0,
  "betaine_m": 2.0,
  "min_k": 12,
  "max_k": 18
}
```

### Host-Free Optimization (no background genome)

```bash
neoswga optimize -j params.json --no-background
```

### Clinical Application (minimize background)

```bash
neoswga optimize -j params.json --optimization-method=background-aware
```

### Fast Optimization (large primer pools)

```bash
neoswga optimize -j params.json --optimization-method=dominating-set
```

### Condition Sweep (find optimal reaction conditions)

```bash
neoswga suggest --genome target.fasta --sweep --output conditions.csv
```

## Choosing Optimization Methods

| Method | Speed | Best For |
|--------|-------|----------|
| `hybrid` | Medium | General use (default) |
| `dominating-set` | Fast | Large primer pools, quick results |
| `background-aware` | Slow | Clinical samples, low background |
| `genetic` | Slow | Complex optimization, exploration |

## Troubleshooting

### "Jellyfish not found"

Install Jellyfish and ensure it's in your PATH:
```bash
jellyfish --version
```

### "No primers pass filtering"

Try relaxing constraints:
```json
{
  "max_bg_freq": 1e-4,
  "max_gini": 0.7
}
```

### Out of memory

Use Bloom filter for large backgrounds:
```bash
neoswga build-filter human_genome.fasta ./filters/
# Then add to params.json:
# "background_bloom_path": "./filters/background_bloom.pkl"
```

## Next Steps

- [User Guide](user-guide.md) - Detailed documentation
- [API Reference](API_REFERENCE.md) - Python API
- [Architecture](ARCHITECTURE_DIAGRAMS.md) - System design

## Example Data

Test with provided example:
```bash
cd tests/integration/phi29_baseline
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```
