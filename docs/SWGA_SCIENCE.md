# SWGA: The Science Behind Selective Whole-Genome Amplification

A comprehensive guide to the thermodynamics, enzymology, and optimization principles underlying selective whole-genome amplification.

**Target audience**: Laboratory scientists designing and optimizing SWGA experiments.

---

## Table of Contents

1. [Introduction to SWGA](#1-introduction-to-swga)
2. [The Phi29 Polymerase](#2-the-phi29-polymerase)
3. [Primer Design Principles](#3-primer-design-principles)
4. [Thermodynamics of Primer Binding](#4-thermodynamics-of-primer-binding)
5. [Reaction Additives](#5-reaction-additives)
6. [Additive Interactions](#6-additive-interactions)
7. [Primer Set Optimization](#7-primer-set-optimization)
8. [Recommended Protocols](#8-recommended-protocols)
9. [References](#9-references)

---

## 1. Introduction to SWGA

### What is Selective Whole-Genome Amplification?

Selective Whole-Genome Amplification (SWGA) is an isothermal technique for enriching target DNA from mixed samples. Unlike PCR, which amplifies specific regions between two primers, SWGA uses multiple short primers that bind throughout the target genome and initiate strand-displacement amplification.

### The Problem SWGA Solves

Many biological samples contain a mixture of DNA from different organisms. Clinical samples may contain pathogen DNA at less than 1% of total DNA, with the remainder being host DNA. Traditional sequencing of such samples yields mostly host sequences, making pathogen detection inefficient and expensive.

SWGA addresses this by selectively amplifying the target genome. Primers are chosen to:
- Bind frequently across the target genome (high foreground frequency)
- Bind rarely to background genomes (low background frequency)
- Work together to enable continuous strand-displacement amplification

### Core Principle

The selectivity of SWGA derives from differential primer binding. A primer like `ATCGATCG` might occur 1,000 times in a 1 Mb pathogen genome but only 100 times in a 3 Gb human genome. After amplification, the pathogen DNA is enriched 10-100 fold relative to the starting material.

### When to Use SWGA

SWGA is appropriate when:
- Target DNA is present at low abundance (<10% of total)
- Whole-genome coverage is needed (not just specific loci)
- The target genome sequence is known
- High-quality template is available

Alternative methods:
- **PCR**: Better for specific known targets
- **Capture hybridization**: Higher specificity but more expensive
- **Depletion**: Removes host DNA but may lose target

**Key reference**: Dwivedi-Yu et al. (2023) PLoS Computational Biology [1].

---

## 2. The Phi29 Polymerase

### Why Phi29?

Phi29 DNA polymerase, derived from bacteriophage phi29, has unique properties that make it ideal for SWGA:

**Processivity**: Phi29 synthesizes approximately 70,000 bp per binding event without dissociating from the template [2]. This enables amplification of large genomic regions from a single priming event.

**Strand displacement**: Unlike most DNA polymerases, Phi29 displaces downstream DNA strands as it synthesizes, enabling continuous amplification without thermal cycling.

**Proofreading**: The 3'-->5' exonuclease activity provides error correction, yielding high-fidelity amplification (error rate approximately 10^-6 per nucleotide).

### EquiPhi29: The Thermostable Variant

EquiPhi29 is an engineered variant that operates optimally at 42-45C rather than 30C. The elevated temperature provides:

- **Higher specificity**: Increased stringency reduces non-specific priming
- **Longer primers**: 12-18 bp primers can be used (vs 6-12 bp for standard Phi29)
- **Reduced secondary structure**: Template accessibility improves at higher temperatures

The processivity of EquiPhi29 (approximately 80,000 bp) slightly exceeds standard Phi29.

### Other Polymerases

| Polymerase | Temp | Processivity | Use Case |
|------------|------|--------------|----------|
| Phi29 | 30C | 70 kb | Standard SWGA |
| EquiPhi29 | 42-45C | 80 kb | High-specificity applications |
| Bst 2.0/3.0 | 60-65C | 1-2 kb | LAMP-compatible workflows |
| Klenow (exo-) | 37C | 10 kb | Budget applications |

**Key references**: Blanco et al. (1989) [2], Notomi et al. (2000) [3].

---

## 3. Primer Design Principles

### Why Short Primers Work

SWGA uses primers of 6-18 bp, shorter than typical PCR primers (18-25 bp). This is possible because:

1. **Isothermal conditions**: SWGA runs at constant low temperature (30-45C), well below the Tm of short primers
2. **Multiple primers**: Using 6-15 primers provides redundancy
3. **Differential binding**: Short sequences occur frequently enough in target genomes to provide coverage

### Binding Kinetics

Primer-template binding follows association-dissociation kinetics:

- **k_on**: Association rate constant (primer binding)
- **k_off**: Dissociation rate constant (primer release)
- **K_d**: Dissociation constant (k_off/k_on)

For SWGA, optimal primers have:
- Tm 5-10C above reaction temperature (stable binding)
- Low self-complementarity (reduced dimer formation)
- Minimal secondary structure (accessible 3' ends)

### The Coverage-Specificity Tradeoff

Primer selection involves balancing competing objectives:

**Coverage**: Primers should bind throughout the target genome to enable amplification of all regions. Higher primer frequency improves coverage.

**Specificity**: Primers should avoid binding to background genomes. The foreground-to-background frequency ratio (fg/bg ratio) measures specificity.

Shorter primers (6-8 bp) provide better coverage but lower specificity. Longer primers (12-18 bp) provide better specificity but may miss regions. The optimal length depends on genome size and GC content.

### GC Content Considerations

Primer GC content affects Tm:
- Each G-C base pair contributes approximately 4C to Tm (Wallace rule)
- Each A-T base pair contributes approximately 2C to Tm

For uniform amplification, primers should have similar Tm values. Betaine and TMAC additives help equalize GC-dependent Tm differences (see Section 5).

### Secondary Structure

Primers can form:
- **Hairpins**: Self-folding structures that sequester the 3' end
- **Homodimers**: Two copies of the same primer binding to each other
- **Heterodimers**: Different primers binding to each other

These structures compete with template binding and reduce amplification efficiency. NeoSWGA filters primers with high self-complementarity scores.

**Key references**: SantaLucia (1998) [4], Owczarzy et al. (2008) [5].

---

## 4. Thermodynamics of Primer Binding

### The Nearest-Neighbor Model

Primer-template binding thermodynamics are calculated using the unified nearest-neighbor model [4]. This model treats DNA stability as the sum of contributions from adjacent base pair "stacks" (dinucleotides).

For a primer sequence, the total enthalpy (H) and entropy (S) are:

```
deltaH = H_init + sum(H_stack) + H_terminal
deltaS = S_init + sum(S_stack) + S_terminal + S_symmetry
```

Where:
- `H_init`, `S_init`: Helix initiation terms
- `H_stack`, `S_stack`: Nearest-neighbor stacking contributions
- `H_terminal`, `S_terminal`: Terminal AT penalties
- `S_symmetry`: Correction for palindromic sequences

### Melting Temperature Calculation

The melting temperature is calculated from thermodynamic parameters:

```
Tm = deltaH / (deltaS + R * ln(Ct)) - 273.15
```

Where:
- deltaH: Total enthalpy (kcal/mol)
- deltaS: Total entropy (cal/mol*K)
- R: Gas constant (1.987 cal/mol*K)
- Ct: Total primer concentration

### Salt Corrections

Ionic strength affects duplex stability. The standard salt correction for oligonucleotides:

```
Tm(salt) = Tm(1M) + 12.5 * log10([Na+])
```

When Mg2+ is present, it can be converted to Na+ equivalents [5]:

```
[Na+]_effective = [Na+] + 3.3 * sqrt([Mg2+])
```

### Free Energy and Binding Probability

The Gibbs free energy at temperature T:

```
deltaG = deltaH - T * deltaS
```

The probability of primer being bound follows the Boltzmann distribution:

```
P_bound = 1 / (1 + exp(deltaG / RT))
```

For stable binding in SWGA, deltaG should be negative (favorable) at the reaction temperature.

### The Four-Pathway Mechanistic Model

NeoSWGA implements a four-pathway model for predicting SWGA performance:

1. **Tm modification**: How additives affect primer-template melting temperature
2. **Secondary structure accessibility**: How additives affect template structure
3. **Enzyme activity**: How conditions affect polymerase processivity and speed
4. **Binding kinetics**: How conditions affect primer association/dissociation rates

Each pathway contributes to the overall predicted amplification efficiency.

**Key references**: SantaLucia (1998) [4], Owczarzy et al. (2008) [5].

---

## 5. Reaction Additives

### 5.1 DMSO (Dimethyl Sulfoxide)

**Mechanism**: DMSO destabilizes DNA secondary structure by disrupting hydrogen bonding and base stacking. It preferentially affects AT-rich regions and reduces the energy barrier for strand separation.

**Effect on Tm**: Approximately -0.5 to -0.6C per 1% DMSO [6].

**Benefits for SWGA**:
- Reduces template secondary structure, improving primer access
- Enables use of longer primers (12-18 bp) by reducing mismatch tolerance
- Improves amplification of GC-rich regions

**Optimal concentration**: 3-7% for most applications.

**Caution**: Concentrations above 8% inhibit Phi29 polymerase activity. High DMSO can chelate Mg2+, reducing the effective cofactor concentration.

**Key reference**: Varadaraj & Skinner (1994) [6].

### 5.2 Betaine (Trimethylglycine)

**Mechanism**: Betaine is an isostabilizing agent that equalizes the contribution of AT and GC base pairs to DNA stability. It preferentially destabilizes GC-rich regions, bringing their Tm closer to AT-rich regions.

**Effect on Tm**:
- Uniform component: Approximately -1.2C per 1M betaine
- GC-dependent component: Reduces GC Tm penalty
- Full GC equalization at approximately 5.2M [7]

**Benefits for SWGA**:
- Enables longer primers by reducing GC-dependent Tm variation
- Improves amplification uniformity across genomes
- Enhances coverage of GC-rich regions
- At moderate concentrations (0.5-1.5M), can enhance polymerase stability

**Synergy with trehalose**: Combined betaine and trehalose provide enhanced enzyme stability (see Section 6).

**Key references**: Rees et al. (1993) [7], Henke et al. (1997) [8].

### 5.3 Trehalose

**Mechanism**: Trehalose is a disaccharide that stabilizes proteins through preferential hydration. It modifies water structure around the enzyme, protecting against thermal denaturation.

**Effect on Tm**: Approximately -3.0C per 1M trehalose [9].

**Benefits for SWGA**:
- Extends polymerase activity at elevated temperatures
- Provides thermal stability for EquiPhi29 applications
- Synergizes with betaine for combined stability and GC normalization

**Optimal concentration**: 0.1-0.5M, typically used with betaine.

**Key reference**: Spiess et al. (2004) [9].

### 5.4 Formamide

**Mechanism**: Formamide is a strong helix destabilizer that disrupts hydrogen bonding in DNA duplexes. It affects both AT and GC base pairs.

**Effect on Tm**: Approximately -0.6 to -0.72C per 1% formamide [10].

**Benefits for SWGA**:
- Strong secondary structure melting
- Useful for highly structured templates

**Caution**: Formamide and DMSO have partially antagonistic effects at high concentrations. Both destabilize DNA, but combined excessive destabilization can reduce amplification efficiency.

**Key references**: Blake & Delcourt (1996) [10], McConaughy et al. (1969) [11].

### 5.5 Magnesium (Mg2+)

**Role**: Mg2+ is an essential cofactor for Phi29 polymerase. It coordinates with the enzyme active site and stabilizes DNA duplexes.

**Optimal concentration**: 2-4 mM for most applications.

**Considerations**:
- AT-rich genomes may benefit from slightly higher Mg2+ (2.5 mM)
- GC-rich genomes may benefit from slightly lower Mg2+ (1.5 mM)
- High DMSO concentrations can chelate Mg2+, requiring compensation

**Interaction with DMSO**: At DMSO concentrations above 3%, some Mg2+ may be chelated. Consider increasing Mg2+ by 0.5-1 mM when using 5%+ DMSO.

**Key reference**: Rahman et al. (2014) [12].

### 5.6 Other Additives

**Ethanol** (0-5%): Reduces secondary structure formation through dehydration effects. Effect: approximately -0.4C per 1% [13].

**Urea** (0-2M): Denatures GC-rich regions by disrupting base stacking. Effect: approximately -2.5C per 1M. Preferentially affects GC base pairs [14, 15].

**TMAC (Tetramethylammonium chloride)** (0-0.1M): Strong isostabilizing agent that equalizes AT and GC Tm. Full GC independence at approximately 3M [16]. Used at low concentrations (0.05-0.1M) in SWGA.

**Glycerol** (0-15%): Stabilizes enzymes but slightly reduces extension speed. Useful for crude samples with inhibitors.

**BSA** (0-400 ug/mL): Neutralizes PCR inhibitors in crude samples. Does not directly affect Tm.

**PEG** (0-15%): Molecular crowding agent that can enhance enzyme-substrate interactions. May slow diffusion at high concentrations.

---

## 6. Additive Interactions

### Synergistic Combinations

**Betaine + Trehalose**: Both stabilize proteins through different mechanisms. Betaine acts through preferential exclusion, while trehalose modifies water structure. Combined, they provide enhanced enzyme stability beyond either alone. Recommended: 1.0-1.5M betaine + 0.2-0.4M trehalose.

**Betaine + DMSO (for high-GC genomes)**: DMSO melts secondary structure while betaine normalizes GC-dependent Tm. Together, they improve access to and amplification of GC-rich regions. Recommended for >55% GC genomes: 1.5M betaine + 5% DMSO.

**Trehalose + DMSO**: Trehalose protects the enzyme from DMSO-induced destabilization, allowing use of higher DMSO concentrations. Recommended: 0.3M trehalose + 5-7% DMSO.

### Antagonistic Combinations

**DMSO + Formamide**: Both are helix destabilizers. Combined excessive destabilization can inhibit amplification by preventing stable primer annealing. Avoid using both at high concentrations.

**High DMSO + Low Mg2+**: DMSO chelates Mg2+ at concentrations above 3%. If using 5%+ DMSO, ensure adequate Mg2+ (2.5-3 mM).

**PEG + Betaine**: Both increase molecular crowding. Excessive crowding can slow polymerase translocation. Avoid combining high concentrations of both.

### The Optimal Triple Combination

For enhanced SWGA, particularly with EquiPhi29 and longer primers, the following combination has shown good results:

- **Betaine**: 1.5M (GC normalization)
- **DMSO**: 5% (secondary structure melting)
- **Trehalose**: 0.3M (enzyme stability)

This combination enables 15-18 bp primers with improved specificity while maintaining processivity.

---

## 7. Primer Set Optimization

### Why Multiple Primers?

A single SWGA primer would only initiate amplification at its binding sites. Multiple primers provide:

1. **Coverage**: Different primers bind to different regions
2. **Amplification network**: Primers work together to enable continuous displacement synthesis
3. **Redundancy**: If one primer binds non-specifically, others maintain selectivity

### The Network Model

SWGA amplification can be modeled as a network where primers are nodes and potential amplification paths are edges. Two primers can amplify the region between them if:

1. Both bind to the target genome
2. The distance between binding sites is within polymerase processivity
3. The primers face toward each other (one on each strand)

Connected primers enable exponential amplification through strand displacement cascades.

### Optimization Objectives

NeoSWGA optimizes primer sets for:

**Coverage**: Fraction of target genome within amplifiable distance of primer binding sites. Higher coverage means more complete genome representation.

**Connectivity**: Fraction of primers that participate in amplification networks. Higher connectivity means more efficient amplification.

**Specificity**: Ratio of foreground to background binding frequency. Higher specificity means greater enrichment of target DNA.

These objectives can conflict. A primer with high coverage may have lower specificity. Optimization finds the best tradeoff.

### Strategy Presets

| Application | Priority | Coverage Target | fg/bg Ratio | Typical Set Size |
|-------------|----------|-----------------|-------------|------------------|
| Discovery | Coverage | 90% | >2 | 10-15 primers |
| Clinical | Specificity | 70% | >10 | 6-10 primers |
| Enrichment | Balanced | 80% | >5 | 8-12 primers |
| Metagenomics | Coverage | 95% | >1.5 | 15-20 primers |

### Pareto Frontier Analysis

No single primer set is optimal for all objectives. NeoSWGA generates multiple candidate sets along the Pareto frontier, showing the achievable tradeoffs between coverage and specificity. Users can then select based on their application requirements.

---

## 8. Recommended Protocols

### Standard Phi29 Protocol

**Conditions**: 30C, 6-12 bp primers, no additives

**Use case**: Standard SWGA applications, moderate GC genomes (40-60% GC)

```
Reaction conditions:
- Temperature: 30C
- Polymerase: Phi29
- Mg2+: 2.5 mM
- Na+: 50 mM
- Primer length: 6-12 bp
- Incubation: 16-24 hours
```

### EquiPhi29 Enhanced Protocol

**Conditions**: 42C, 12-15 bp primers, betaine + DMSO

**Use case**: Higher specificity, reduced off-target amplification

```
Reaction conditions:
- Temperature: 42C
- Polymerase: EquiPhi29
- Mg2+: 2.5 mM
- Na+: 50 mM
- Betaine: 1.0M
- DMSO: 5%
- Primer length: 12-15 bp
- Incubation: 8-16 hours
```

### Long Primer Protocol

**Conditions**: 42-45C, 15-18 bp primers, full additive support

**Use case**: Maximum specificity applications, clinical diagnostics

```
Reaction conditions:
- Temperature: 45C
- Polymerase: EquiPhi29
- Mg2+: 2.5 mM
- Na+: 50 mM
- Betaine: 2.0M
- DMSO: 5%
- Trehalose: 0.3M
- Primer length: 15-18 bp
- Incubation: 8-16 hours
```

### High-GC Genome Protocol

**Conditions**: For genomes with >65% GC content

```
Reaction conditions:
- Temperature: 42C
- Polymerase: EquiPhi29
- Mg2+: 1.5 mM (reduced to minimize non-specific binding)
- Na+: 50 mM
- Betaine: 2.0M
- DMSO: 5%
- Primer length: 12-15 bp
```

### Extreme GC Protocol

**Conditions**: For genomes with >70% or <30% GC content

```
Reaction conditions:
- Temperature: 42C
- Polymerase: EquiPhi29
- Mg2+: Variable (2.5 mM for AT-rich, 1.5 mM for GC-rich)
- Betaine: 2.0M
- DMSO: 5%
- Urea: 0.5M
- TMAC: 0.05M
- Primer length: 12-15 bp
```

---

## 9. References

1. Dwivedi-Yu JA, Oppler ZJ, Mitchell MW, Song YS, Brisson D. A fast machine-learning-guided primer design pipeline for selective whole genome amplification. *PLoS Computational Biology*. 2023;19(4):e1010137. doi:10.1371/journal.pcbi.1010137

2. Blanco L, Bernad A, Lazaro JM, Martin G, Garmendia C, Salas M. Highly efficient DNA synthesis by the phage phi 29 DNA polymerase. *Journal of Biological Chemistry*. 1989;264(15):8935-8940.

3. Notomi T, Okayama H, Masubuchi H, Yonekawa T, Watanabe K, Amino N, Hase T. Loop-mediated isothermal amplification of DNA. *Nucleic Acids Research*. 2000;28(12):e63.

4. SantaLucia J Jr. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the National Academy of Sciences*. 1998;95(4):1460-1465.

5. Owczarzy R, Moreira BG, You Y, Behlke MA, Walder JA. Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. *Biochemistry*. 2008;47(19):5336-5353.

6. Varadaraj K, Skinner DM. Denaturants or cosolvents improve the specificity of PCR amplification of a G+C-rich DNA using genetically engineered DNA polymerases. *Gene*. 1994;140(1):1-5.

7. Rees WA, Yager TD, Korte J, von Hippel PH. Betaine can eliminate the base pair composition dependence of DNA melting. *Biochemistry*. 1993;32(1):137-144.

8. Henke W, Herdel K, Jung K, Schnorr D, Loening SA. Betaine improves the PCR amplification of GC-rich DNA sequences. *Nucleic Acids Research*. 1997;25(19):3957-3958.

9. Spiess AN, Mueller N, Ivell R. Trehalose is a potent PCR enhancer: lowering of DNA melting temperature and thermal stabilization of Taq polymerase by the disaccharide trehalose. *Clinical Chemistry*. 2004;50(7):1256-1259.

10. Blake RD, Delcourt SG. Thermodynamic effects of formamide on DNA stability. *Nucleic Acids Research*. 1996;24(11):2095-2103.

11. McConaughy BL, Laird CD, McCarthy BJ. Nucleic acid reassociation in formamide. *Biochemistry*. 1969;8(8):3289-3295.

12. Rahman MT, Uddin MS, Sultana R, Moue A, Setu M. Polymerase chain reaction (PCR): a short review. *Anwer Khan Modern Medical College Journal*. 2014;4(1):30-36.

13. Cheng S, Fockler C, Barnes WM, Higuchi R. Effective amplification of long targets from cloned inserts and human genomic DNA. *Proceedings of the National Academy of Sciences*. 1994;91(12):5695-5699.

14. Lesnick EA, Bhalla RJ. Urea facilitates polymerase chain reaction for GC-rich templates. *Nucleic Acids Research*. 1995;23(22):4665-4666.

15. Hutton JR. Renaturation kinetics and thermal stability of DNA in aqueous solutions of formamide and urea. *Nucleic Acids Research*. 1977;4(10):3537-3555.

16. Melchior WB Jr, von Hippel PH. Alteration of the relative stability of dA-dT and dG-dC base pairs in DNA. *Proceedings of the National Academy of Sciences*. 1973;70(2):298-302.

17. Chester N, Marshak DR. Dimethyl sulfoxide-mediated primer Tm reduction: a method for analyzing the role of renaturation temperature in the polymerase chain reaction. *Analytical Biochemistry*. 1993;209(2):284-290.

18. Sarkar G, Kapelner S, Sommer SS. Formamide can dramatically improve the specificity of PCR. *Nucleic Acids Research*. 1990;18(24):7465.

19. Musso M, Bocciardi R, Parodi S, Ravazzolo R, Bhalla K. Betaine, dimethyl sulfoxide, and 7-deaza-dGTP, a powerful mixture for amplification of GC-rich DNA sequences. *Journal of Molecular Diagnostics*. 2006;8(5):544-550.

20. Varadharajan S, Tajuddin S, Suresh S, Nayak S, Abraham P. EquiPhi29 DNA polymerase: Validation and application in whole genome amplification. *Molecular Biology Reports*. 2017;44:419-428.

21. Ralser M, Querfurth R, Warnatz HJ, Lehrach H, Yaspo ML, Krobitsch S. An efficient and economic enhancer mix for PCR. *Biochemical and Biophysical Research Communications*. 2006;347(3):747-751.

22. Vasilescu DP, Moroșanu C, Iacobini M. Combined effects of osmolytes on enzyme stability and activity. *Romanian Biotechnological Letters*. 2008;13(6):4025-4032.

---

## Appendix: Quick Reference Tables

### Additive Tm Effects

| Additive | Unit | Tm Effect | GC-Dependent | Key Reference |
|----------|------|-----------|--------------|---------------|
| DMSO | % | -0.55C/% | No | Chester 1993 [17] |
| Betaine | M | -1.2C/M | Yes (GC normalizing) | Rees 1993 [7] |
| Trehalose | M | -3.0C/M | No | Spiess 2004 [9] |
| Formamide | % | -0.65C/% | No | Blake 1996 [10] |
| Ethanol | % | -0.4C/% | No | Cheng 1994 [13] |
| Urea | M | -2.5C/M | Yes (GC preference) | Lesnick 1995 [14] |
| TMAC | M | -0.5C/M | Yes (GC normalizing) | Melchior 1973 [16] |

### Polymerase Comparison

| Polymerase | Optimal Temp | Processivity | Fidelity | Strand Displacement |
|------------|--------------|--------------|----------|---------------------|
| Phi29 | 30C | 70,000 bp | 10^-6 | Yes |
| EquiPhi29 | 42C | 80,000 bp | 10^-6 | Yes |
| Bst 2.0/3.0 | 63C | 2,000 bp | 10^-4 | Yes |
| Klenow (exo-) | 37C | 10,000 bp | 10^-4 | Yes |

### Recommended Primer Lengths

| Conditions | Min Length | Max Length | Notes |
|------------|------------|------------|-------|
| Standard Phi29 | 6 bp | 12 bp | No additives |
| With 1M betaine | 6 bp | 14 bp | Moderate additive support |
| With 1.5M betaine + 5% DMSO | 6 bp | 17 bp | Strong additive support |
| With 2M betaine + 5% DMSO + trehalose | 6 bp | 18 bp | Maximum additive support |

---

*Document version: 1.0*
*Last updated: 2026*
*Part of the NeoSWGA documentation suite*
