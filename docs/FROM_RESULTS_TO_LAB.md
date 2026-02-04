# From Results to Lab: Complete Workflow Guide

This guide covers the complete workflow from running NeoSWGA to having primers ready for your SWGA experiment in the lab.

## Overview

After running the NeoSWGA pipeline, you have optimized primer sequences. This guide covers:

1. Understanding your results
2. Exporting primers for ordering
3. Ordering from synthesis vendors
4. Receiving and preparing primers
5. Setting up your SWGA reaction

---

## 1. Understanding Your Results

### Pipeline Output Files

After running the four pipeline steps (`count-kmers`, `filter`, `score`, `optimize`), your results directory contains:

| File | Description | Key Columns |
|------|-------------|-------------|
| `step4_improved_df.csv` | Final optimized primer set(s) | seq, set_score, fg_freq, bg_freq, Tm |
| `step3_df.csv` | Scored candidates (before optimization) | seq, amp_pred, Tm, gini, fg_freq |
| `step2_df.csv` | Filtered candidates | seq, fg_freq, bg_freq, gini, Tm |
| `step1_*_Xmer_all.txt` | K-mer counts (X = primer length) | primer, count |
| `params.json` | Parameters used for the run | All pipeline settings |

### Key Metrics to Check

Before ordering primers, verify these metrics in `step4_improved_df.csv`:

**Quality indicators:**
- `fg_freq`: Foreground frequency (higher is better, typically 1e-5 to 1e-3)
- `bg_freq`: Background frequency (lower is better, should be < max_bg_freq)
- `gini`: Binding evenness (lower is better, < 0.6 ideal)
- `Tm`: Melting temperature (should match your reaction temperature)
- `set_score`: Overall set quality score (higher is better)

**Coverage indicators:**
- `coverage_fraction`: Genome coverage (target 70-90% depending on application)
- `enrichment`: Foreground/background ratio (higher is better for specificity)

### Checking Quality

Use the built-in quality assessment tools:

```bash
# Quick quality assessment with recommendations
neoswga interpret -d ./results/

# Generate detailed HTML report
neoswga report -d ./results/

# Full technical report
neoswga report -d ./results/ --level full
```

**Interpretation guidelines:**

| Grade | Coverage | Specificity | Recommendation |
|-------|----------|-------------|----------------|
| A | >80% | >90% | Excellent, proceed with confidence |
| B | 70-80% | 80-90% | Good, suitable for most applications |
| C | 60-70% | 70-80% | Acceptable, consider re-optimization |
| D | 50-60% | 60-70% | Marginal, likely needs adjustment |
| F | <50% | <60% | Poor, re-run with different parameters |

### If Results Are Suboptimal

**Poor coverage (<70%):**
```bash
# Increase primer set size
neoswga optimize -j params.json --target-set-size 10

# Try background-aware optimizer for better specificity
neoswga optimize -j params.json --optimization-method background-aware
```

**Poor specificity (high background):**
```bash
# Use stricter background filtering
# Edit params.json: "max_bg_freq": 1e-5 (was 5e-5)

# Or use background-aware optimizer
neoswga optimize -j params.json --optimization-method background-aware
```

**Inconsistent Tm values:**
```bash
# Tighten Tm range in params.json
# "min_tm": 55, "max_tm": 65 (instead of 50-70)
```

---

## 2. Exporting Primers for Ordering

### Quick Export (Recommended)

Export all formats at once with a single command:

```bash
# Export to new directory with project name
neoswga export -d ./results/ -o ./order/ --project MyPathogen

# This creates:
#   order/MyPathogen_primers.fasta         (FASTA with metadata)
#   order/MyPathogen_order_idt.csv         (IDT bulk upload)
#   order/MyPathogen_protocol.md           (Lab protocol)
```

### Export Options

**Different vendor formats:**
```bash
# IDT format (default)
neoswga export -d ./results/ -o ./order/ --vendor idt

# Twist Bioscience format
neoswga export -d ./results/ -o ./order/ --vendor twist

# Sigma-Aldrich format
neoswga export -d ./results/ -o ./order/ --vendor sigma

# Generic CSV format
neoswga export -d ./results/ -o ./order/ --vendor generic
```

**Single format export:**
```bash
# FASTA only
neoswga export -d ./results/ -o ./order/ --format fasta

# Vendor CSV only
neoswga export -d ./results/ -o ./order/ --format csv --vendor idt

# Protocol only
neoswga export -d ./results/ -o ./order/ --format protocol
```

**With custom parameters file:**
```bash
# Use params.json for reaction conditions in protocol
neoswga export -d ./results/ -o ./order/ -j params.json --project MyPathogen
```

### Output File Formats

**IDT Format (order_idt.csv):**
```
Name,Sequence,Scale,Purification
MyPathogen_001,ATCGATCGATCG,25nm,STD
MyPathogen_002,GCTAGCTAGCTA,25nm,STD
```
- Ready for IDT bulk input
- Upload directly to cart

**Twist Format (order_twist.csv):**
```
Name,Sequence
MyPathogen_001,ATCGATCGATCG
MyPathogen_002,GCTAGCTAGCTA
```
- Compatible with Twist Bioscience ordering

**Generic CSV Format (order_generic.csv):**
```
name,sequence,length,tm,gc
MyPathogen_001,ATCGATCGATCG,12,36.0,50.0%
MyPathogen_002,GCTAGCTAGCTA,12,36.0,50.0%
```
- Includes Tm and GC content
- Use for labs with custom ordering systems

**FASTA Format (primers.fasta):**
```
>MyPathogen_001 Tm=36.0C GC=50.0% len=12
ATCGATCGATCG
>MyPathogen_002 Tm=36.0C GC=50.0% len=12
GCTAGCTAGCTA
```
- Standard sequence format with metadata
- For archival and analysis

---

## 3. Ordering from Synthesis Vendors

### Recommended Vendors

| Vendor | Turnaround | Price Range | Best For |
|--------|------------|-------------|----------|
| **IDT (Integrated DNA Technologies)** | 1-2 days | $$ | Fast turnaround, reliable quality |
| **Twist Bioscience** | 5-7 days | $ | Bulk orders, primer pools |
| **Sigma-Aldrich** | 3-5 days | $ | Budget option, global availability |
| **Eurofins Genomics** | 2-3 days | $$ | European labs |

### Order Specifications

**Standard SWGA primers (6-12 bp):**

| Parameter | Recommendation | Notes |
|-----------|----------------|-------|
| **Scale** | 25 nmol | Sufficient for 100+ reactions |
| **Purification** | Desalt (standard) | Adequate for most applications |
| **Format** | Dry or 100 uM in TE | Dry for long-term storage |
| **Modifications** | None | No 5' or 3' modifications needed |

**Long primers (15-18 bp, EquiPhi29):**

| Parameter | Recommendation | Notes |
|-----------|----------------|-------|
| **Scale** | 25-50 nmol | Larger scale recommended |
| **Purification** | HPLC | Higher quality for longer oligos |
| **Format** | 100 uM in TE | Ready-to-use |

### Vendor-Specific Instructions

#### IDT (Integrated DNA Technologies)

**Web ordering:**
1. Go to https://www.idtdna.com/
2. Select "Bulk Input"
3. Upload `*_primers_idt.csv` file
4. Review cart, select shipping
5. Typical turnaround: 1-2 business days (standard desalt)

**Tips:**
- Order before 2 PM EST for same-day processing
- Use "Standard Desalt" for 6-12 bp primers
- Upgrade to "HPLC" only for primers >15 bp or critical applications
- Request "dry" format for long-term storage

#### Twist Bioscience

**Web ordering:**
1. Go to https://www.twistbioscience.com/
2. Select "Oligo Pools" for primer sets
3. Upload `*_primers_twist.csv` file
4. Specify concentration (recommend 100 uM per primer)
5. Typical turnaround: 5-7 business days

**Tips:**
- Best for ordering 10+ primers at once
- Can request primers pre-mixed at equimolar concentration
- Ideal for multi-genome or large primer sets

#### Sigma-Aldrich

**Web ordering:**
1. Go to https://www.sigmaaldrich.com/
2. Select "Custom DNA Oligos"
3. Enter primers individually or upload CSV
4. Select "Desalt" purification
5. Typical turnaround: 3-5 business days

**Tips:**
- Good global availability
- Budget-friendly option
- Reliable quality for short primers

### Cost Estimation

**Approximate costs (USD, as of 2024):**

| Scale | Purification | Price per Primer | 8-primer Set Cost |
|-------|--------------|------------------|-------------------|
| 25 nmol | Desalt | $5-10 | $40-80 |
| 100 nmol | Desalt | $10-15 | $80-120 |
| 25 nmol | HPLC | $20-30 | $160-240 |

**Budget considerations:**
- Start with 25 nmol desalt for initial testing
- Order 100 nmol if you'll be running many experiments
- HPLC purification usually not necessary for SWGA

---

## 4. Receiving and Preparing Primers

### Upon Arrival

**Inspection:**
1. Check that all primers arrived (compare to order sheet)
2. Verify lot numbers match between tubes and documentation
3. Note any shipping delays or temperature excursions

**Initial handling:**
1. Centrifuge tubes briefly (5 sec, 1000g) before opening
2. Keep tubes closed until ready to resuspend
3. Work in a clean area to avoid contamination

### Resuspension (Dry Oligos)

**Standard protocol (100 uM stock):**

1. Calculate volume needed:
   ```
   Volume (uL) = nmoles ordered * 10

   Example: 25 nmol primer -> 250 uL TE buffer
            100 nmol primer -> 1000 uL TE buffer
   ```

2. Add calculated volume of TE buffer (10 mM Tris pH 8.0, 0.1 mM EDTA)

3. Vortex thoroughly (30 sec) to dissolve

4. Centrifuge briefly to collect solution

5. Let stand 5 min at room temperature

6. Vortex again and store at -20C

**Quality check:**
```bash
# Measure concentration using NanoDrop or similar
# Expected: 100 uM (100 pmol/uL)
# A260/A280 ratio: 1.8-2.0 (pure DNA)
```

### Creating Working Stocks

**Individual primer working stocks (10 uM):**
```
For each primer:
- Mix 10 uL of 100 uM stock
- Add 90 uL of nuclease-free water
- Label as "10 uM working stock"
- Store at -20C
```

**Equimolar primer mix (10 uM each primer):**
```
For 8 primers in set:
- Take 10 uL of each 100 uM stock (8 primers = 80 uL total)
- Add 720 uL nuclease-free water
- Final: 800 uL at 10 uM per primer
- Label: "Primer Mix - 10 uM each"
- Aliquot into 100 uL portions
- Store at -20C
```

**Why equimolar mix?**
- Simplifies reaction setup
- Ensures balanced primer concentrations
- Reduces pipetting errors
- More reproducible results

### Storage Guidelines

| Stock Type | Concentration | Storage Temp | Stability |
|------------|---------------|--------------|-----------|
| Stock | 100 uM | -20C | >2 years |
| Working stock | 10 uM | -20C | 1 year |
| Primer mix | 10 uM each | -20C | 6 months |
| Thawed aliquot | Any | 4C | 1 week |

**Best practices:**
- Avoid repeated freeze-thaw (max 3 cycles)
- Use small aliquots (50-100 uL) of working stocks
- Keep track of freeze-thaw cycles on tube labels
- Store in frost-free areas of freezer

### Record Keeping

Update your order sheet (`*_order_sheet.csv`) with:
- Date received
- Lot numbers from vendor
- Concentration after resuspension (if measured)
- Date of first use
- Freeze-thaw cycles

---

## 5. Setting Up Your SWGA Reaction

### Reaction Components

NeoSWGA designs primers for specific polymerase and temperature conditions. Use the conditions matching your pipeline run.

### Standard Protocol (Phi29, 30C)

**For primers designed with:**
```json
{
  "polymerase": "phi29",
  "reaction_temp": 30
}
```

**25 uL reaction:**

| Component | Volume | Final Concentration | Source |
|-----------|--------|---------------------|--------|
| 10X Phi29 Buffer | 2.5 uL | 1X | NEB or Qiagen |
| dNTPs (10 mM each) | 1.0 uL | 0.4 mM each | Standard mix |
| Primer mix (10 uM each) | 3.0 uL | 200 nM each primer | Your mix |
| Phi29 Polymerase (10 U/uL) | 0.5 uL | 10 U total | NEB M0269 |
| Template DNA | 1.0 uL | 1-10 ng total | Your sample |
| Nuclease-free water | 17.0 uL | - | - |
| **Total** | **25.0 uL** | | |

**Thermal cycling:**
```
30C for 16-24 hours (overnight)
65C for 10 min (enzyme inactivation)
4C hold
```

**Expected yield:** 5-20 ug DNA from 10 ng input

### Enhanced Protocol (EquiPhi29, 42-45C)

**For primers designed with:**
```json
{
  "polymerase": "equiphi29",
  "reaction_temp": 42,
  "dmso_percent": 5.0,
  "betaine_m": 1.0
}
```

**25 uL reaction:**

| Component | Volume | Final Concentration |
|-----------|--------|---------------------|
| 10X EquiPhi29 Buffer | 2.5 uL | 1X |
| dNTPs (10 mM each) | 1.0 uL | 0.4 mM each |
| Primer mix (10 uM each) | 3.0 uL | 200 nM each |
| DMSO (100%) | 1.25 uL | 5% |
| Betaine (5 M) | 5.0 uL | 1 M |
| EquiPhi29 Polymerase | 0.5 uL | 10 U |
| Template DNA | 1.0 uL | 1-10 ng |
| Water | 10.75 uL | - |
| **Total** | **25.0 uL** | |

**Thermal cycling:**
```
42C for 16-24 hours
65C for 10 min
4C hold
```

**Benefits of enhanced protocol:**
- Higher specificity (longer primers, higher temperature)
- Better for GC-rich targets
- Reduced non-specific amplification

### Reaction Setup Tips

**Order of addition:**
1. Make master mix without template (for multiple reactions)
2. Aliquot master mix into PCR tubes
3. Add template DNA last (varies by sample)
4. Mix gently by pipetting
5. Quick spin to collect

**Master mix example (for 10 reactions + 10% excess):**
```
Component                 1X      11X (with excess)
10X Buffer               2.5 uL   27.5 uL
dNTPs                    1.0 uL   11.0 uL
Primer mix               3.0 uL   33.0 uL
DMSO (if using)          1.25 uL  13.75 uL
Betaine (if using)       5.0 uL   55.0 uL
Polymerase               0.5 uL   5.5 uL
Water                    10.75 uL 118.25 uL
Total per reaction:      24 uL
Add per tube:            1 uL template
```

### Controls

**Essential controls:**

| Control | Template | Primers | Purpose |
|---------|----------|---------|---------|
| Positive | Foreground DNA | Full set | Verify amplification |
| Negative | Water | Full set | Check for contamination |
| Background | Background DNA | Full set | Assess specificity |

**Optional controls:**
- No-primer control (template + polymerase only)
- Individual primer tests (one primer at a time)

### Validation

**After amplification, check:**

1. **Yield (NanoDrop or Qubit):**
   - Expected: 200-1000 ng/uL
   - Indicates successful amplification

2. **Size distribution (gel or Bioanalyzer):**
   - Expected: Smear from 500 bp to 10 kb
   - Indicates random priming and amplification

3. **Enrichment (qPCR or sequencing):**
   - Compare foreground vs background ratio
   - Expected: 10-1000X enrichment

**Example validation:**
```bash
# If you have access to qPCR
# Design qPCR assays for foreground and background loci
# Compare Ct values: Lower Ct = higher abundance

Foreground locus Ct: 20 (abundant)
Background locus Ct: 30 (depleted)
Enrichment: 2^(30-20) = 1000X
```

### Optimization

**If amplification is weak:**
- Increase primer concentration (try 400 nM each)
- Extend incubation time (24-48 hours)
- Increase template input (up to 50 ng)
- Add more polymerase (20 U instead of 10 U)

**If background is high:**
- Decrease primer concentration (try 100 nM each)
- Decrease template input (try 1 ng)
- Use tighter temperature control (avoid fluctuations)
- Consider re-running NeoSWGA with stricter parameters

**If primers form dimers:**
- Reduce primer concentration
- Increase temperature (if using EquiPhi29)
- Check for predicted dimers in results:
  ```bash
  neoswga analyze-dimers --primers-file ./results/step4_improved_df.csv
  ```

---

## Troubleshooting

### Common Problems and Solutions

| Problem | Possible Cause | Solution |
|---------|----------------|----------|
| No amplification | Poor primer quality | Reorder primers, verify sequences |
| | Wrong polymerase | Check params.json, match conditions |
| | Template degraded | Use fresh template, check quality |
| | Incorrect temperature | Verify thermocycler calibration |
| Low yield | Insufficient primers | Increase primer concentration |
| | Short incubation | Extend to 24 hours |
| | Inactive polymerase | Use fresh enzyme, check storage |
| Poor enrichment | Non-specific priming | Use EquiPhi29 + additives |
| | Background contamination | Change reagents, clean workspace |
| | Suboptimal primer set | Re-run optimization |
| High background | Too many primers | Reduce concentration |
| | Low temperature | Increase to 42-45C (EquiPhi29) |
| | Re-annealing | Minimize freeze-thaw cycles |
| Primer dimers | Self-complementary primers | Check with analyze-dimers |
| | High concentration | Reduce to 100 nM each |
| Batch variation | Inconsistent mixing | Use single master mix |
| | Temperature variation | Calibrate thermocycler |
| | Different primer lots | Re-validate each lot |

### Diagnostic Steps

**Step 1: Verify primer quality**
```bash
# Check primer concentrations
# Run on agarose gel (should see single band at expected size)
```

**Step 2: Test individual primers**
```bash
# Set up reactions with one primer at a time
# Identify which primer(s) are problematic
```

**Step 3: Check for dimers**
```bash
neoswga analyze-dimers \
  --primers-file ./results/step4_improved_df.csv \
  --output ./dimers/ \
  --visualize
```

**Step 4: Re-optimize if needed**
```bash
# Try stricter parameters
neoswga filter -j params.json  # Edit params.json first
neoswga score -j params.json
neoswga optimize -j params.json --optimization-method background-aware
```

### Getting Help

**Check documentation:**
```bash
# View reaction condition details
cat $NEOSWGA_DIR/docs/SWGA_SCIENCE.md

# View optimization strategies
neoswga optimize --help
```

**Analyze your results:**
```bash
# Generate detailed quality report
neoswga report -d ./results/ --level full

# Interpret with recommendations
neoswga interpret -d ./results/
```

**Community resources:**
- GitHub Issues: Report bugs or ask questions
- Discussions: Share experiences and protocols
- Examples: See `tests/integration/` for working examples

---

## Appendix A: Reagent Sources

### Polymerases

| Enzyme | Vendor | Catalog Number | Notes |
|--------|--------|----------------|-------|
| Phi29 DNA Polymerase | NEB | M0269 | Most common, 10 U/uL |
| EquiPhi29 Polymerase | ThermoFisher | A39391 | Higher temperature, 10 U/uL |
| Bst DNA Polymerase | NEB | M0275 | For high-temp applications |

### Additives

| Reagent | Source | Catalog | Notes |
|---------|--------|---------|-------|
| DMSO (molecular biology grade) | Sigma | D8418 | Use 100% stock |
| Betaine (5 M solution) | Sigma | B0300 | Ready-to-use |
| Trehalose | Sigma | T9449 | Make 1 M stock |

### Buffers and dNTPs

| Reagent | Source | Catalog |
|---------|--------|---------|
| 10X Phi29 Buffer | NEB | Included with enzyme |
| dNTP Set (100 mM each) | NEB | N0447 | Dilute to 10 mM |
| TE Buffer (pH 8.0) | Various | Standard recipe |

---

## Appendix B: Quick Reference Cards

### Export Cheat Sheet

```bash
# Standard export (all formats)
neoswga export -d ./results/ -o ./order/ --project MyProject

# Specific vendor
neoswga export -d ./results/ -o ./order/ --vendor twist --project MyProject

# Single format
neoswga export -d ./results/ -o ./order/ --format fasta --project MyProject

# With params file for protocol conditions
neoswga export -d ./results/ -o ./order/ -j params.json --project MyProject
```

### Reaction Setup Card (Phi29, 30C)

```
25 uL Reaction:
  10X Buffer:     2.5 uL
  dNTPs (10 mM):  1.0 uL
  Primers (10 uM):3.0 uL
  Phi29 (10 U/uL):0.5 uL
  Template:       1.0 uL
  Water:          17.0 uL

Incubate: 30C, 16-24h
Inactivate: 65C, 10 min
```

### Reaction Setup Card (EquiPhi29, 42C)

```
25 uL Reaction:
  10X Buffer:     2.5 uL
  dNTPs (10 mM):  1.0 uL
  Primers (10 uM):3.0 uL
  DMSO (100%):    1.25 uL
  Betaine (5 M):  5.0 uL
  EquiPhi29:      0.5 uL
  Template:       1.0 uL
  Water:          10.75 uL

Incubate: 42C, 16-24h
Inactivate: 65C, 10 min
```

---

## Related Documentation

- **[SWGA Science](SWGA_SCIENCE.md)**: Thermodynamics, mechanism, and theory
- **[README](../README.md)**: Quick start and installation
- **[CLAUDE.md](../CLAUDE.md)**: Complete technical reference
- **[Examples](../tests/integration/)**: Working examples for different conditions

---

**Last updated:** 2026-02-04
**Version:** 1.0.0
