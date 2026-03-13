# Migration Guide: NeoSWGA 3.5 (Genome-Adaptive QA)

**Upgrading from NeoSWGA < 3.5 to 3.5+**

This guide helps existing Neo SWGA users smoothly migrate to version 3.5, which introduces genome-adaptive quality assessment for improved primer design on extreme genomes.

---

## Table of Contents

1. [What Changed](#what-changed)
2. [Backward Compatibility](#backward-compatibility)
3. [Should I Upgrade?](#should-i-upgrade)
4. [Migration Steps](#migration-steps)
5. [Testing After Upgrade](#testing-after-upgrade)
6. [Expected Behavior Changes](#expected-behavior-changes)
7. [Rollback Instructions](#rollback-instructions)
8. [FAQ](#faq)

---

## What Changed

### New Features (Version 3.5)

**1. Genome-Adaptive Quality Assessment**
- Automatically adjusts QA thresholds based on genome GC content
- Improves primer coverage for AT-rich (< 40% GC) and GC-rich (> 60% GC) genomes
- Expected improvement: +75-200% primer coverage for extreme genomes

**2. Enhanced Genome Analysis**
- New `analyze_genome_for_qa()` function analyzes genome and provides recommendations
- Automatic GC classification (extreme_at, at_rich, balanced, gc_rich, extreme_gc)
- Expected improvement estimates for each genome type

**3. Updated Quality Scorer**
- `create_quality_scorer()` now accepts optional `genome_gc` parameter
- Automatically applies adaptive thresholds when `genome_gc` provided
- Maintains standard behavior when `genome_gc=None`

**4. New Primer Selection Strategies**
- Minimal set selection (dominating set algorithm)
- Cooperative binding network analysis
- Hybrid core+supplemental primer sets

### What Didn't Change

- All existing APIs remain functional
- Default behavior unchanged (unless you provide a genome file)
- Standard QA still available (genome_gc=None)
- File formats unchanged
- Configuration parameters unchanged

---

## Backward Compatibility

### 100% Backward Compatible

NeoSWGA 3.5 is **fully backward compatible** with previous versions:

**Existing Code Works Unchanged**:
```python
# This code from v3.4 works identically in v3.5
from neoswga.core.integrated_quality_scorer import create_quality_scorer

scorer = create_quality_scorer(stringency='moderate')
score = scorer.score_primer("ATGCATGCATGC")

# Behavior: Uses standard QA (50% GC baseline) - same as v3.4
```

**No Breaking Changes**:
- All function signatures preserved
- All return types unchanged
- All parameters optional (genome_gc defaults to None)
- Standard QA thresholds unchanged

**Opt-In, Not Forced**:
- Adaptive QA only activates when you provide genome_gc
- Existing workflows continue using standard QA
- You choose when to migrate

---

## Should I Upgrade?

### Upgrade If:

✓ You work with **AT-rich genomes** (< 40% GC)
  - Francisella, Plasmodium, Borrelia, etc.
  - Benefit: +75-200% primer coverage

✓ You work with **GC-rich genomes** (> 60% GC)
  - Burkholderia, Streptomyces, Mycobacterium, etc.
  - Benefit: +70-150% primer coverage

✓ Standard primer design yields **too few primers**
  - Struggling to get enough primers passing QA
  - Need better genome coverage

✓ You want **improved primer selection algorithms**
  - Minimal set selection
  - Cooperative binding networks
  - Hybrid strategies

### Stay on Current Version If:

- You only work with **balanced genomes** (40-60% GC)
- Current primer design works well
- You prefer stability over new features
- You need time to test before upgrading production systems

---

## Migration Steps

### Step 1: Backup Current Environment

```bash
# Save current version
pip freeze > requirements_backup.txt

# Backup existing primers and results
cp -r /path/to/results /path/to/results_backup
```

### Step 2: Update NeoSWGA

```bash
# Update to version 3.5
pip install --upgrade neoswga

# Verify version
python -c "import neoswga; print(neoswga.__version__)"
# Should print: 3.5.0 or higher
```

### Step 3: Test Basic Functionality

```python
# test_upgrade.py - Verify standard QA still works
from neoswga.core.integrated_quality_scorer import create_quality_scorer

# Test 1: Standard QA (should work identically to v3.4)
scorer = create_quality_scorer(stringency='moderate')
test_primer = "ATGCATGCATGC"
score = scorer.score_primer(test_primer)

print(f"Standard QA test:")
print(f"  Primer: {test_primer}")
print(f"  Score: {score.overall_score:.3f}")
print(f"  Passes: {score.passes_all}")

# Expected: Same results as v3.4
```

### Step 4: Test Adaptive QA (New Feature)

```python
# test_adaptive_qa.py - Try new adaptive QA
from neoswga.core.genome_analysis import analyze_genome_for_qa
from neoswga.core.integrated_quality_scorer import create_quality_scorer

# Analyze your genome
analysis = analyze_genome_for_qa('path/to/your/genome.fna')

print(f"Genome Analysis:")
print(f"  GC content: {analysis['gc_content']:.1%}")
print(f"  GC class: {analysis['adaptive_qa_recommendation']['gc_class']}")
print(f"  Use adaptive: {analysis['adaptive_qa_recommendation']['use_adaptive']}")
print(f"  Expected benefit: {analysis['adaptive_qa_recommendation']['expected_improvement']}")

# Create adaptive scorer
scorer = create_quality_scorer(
    stringency='moderate',
    genome_gc=analysis['gc_content']
)

# Test scoring
score = scorer.score_primer("ATATATATATAT")
print(f"\nAdaptive QA test:")
print(f"  Score: {score.overall_score:.3f}")
```

### Step 5: Run Full Pipeline Test

```python
# test_full_pipeline.py - End-to-end test
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

# Test with your genome
generator = OptimalOligoGenerator(
    genome_file='path/to/your/genome.fna'
)

result = generator.generate_primers(
    primer_count=10,
    length_range=(10, 14),
    stringency='moderate'
)

print(f"Full Pipeline Test:")
print(f"  Primers generated: {len(result.primers)}")
print(f"  Genome GC: {result.reaction_conditions.genome_gc:.1%}")
print(f"  Adaptive QA used: {result.reaction_conditions.use_adaptive_qa}")

# Review primers
for primer_id, data in list(result.primers.items())[:5]:
    print(f"  {primer_id}: quality={data['quality']:.3f}, sites={len(data['binding_sites'])}")
```

### Step 6: Compare Results (Optional)

Compare primers from v3.4 (standard) vs v3.5 (adaptive):

```python
# compare_versions.py
from neoswga.core.integrated_quality_scorer import create_quality_scorer

# Your genome GC (e.g., 32% for Francisella)
genome_gc = 0.32

# Standard QA (v3.4 behavior)
scorer_standard = create_quality_scorer(stringency='moderate', genome_gc=None)

# Adaptive QA (v3.5 feature)
scorer_adaptive = create_quality_scorer(stringency='moderate', genome_gc=genome_gc)

# Test primers from your v3.4 results
old_primers = ["ATATATATATAT", "GCGCGCGCGCGC"]  # Your actual primers

for primer in old_primers:
    score_old = scorer_standard.score_primer(primer).overall_score
    score_new = scorer_adaptive.score_primer(primer).overall_score

    print(f"{primer}:")
    print(f"  v3.4 (standard): {score_old:.3f}")
    print(f"  v3.5 (adaptive): {score_new:.3f}")
    print(f"  Change: {score_new - score_old:+.3f}")
```

---

## Testing After Upgrade

### Recommended Test Suite

**1. Backward Compatibility Test**
```bash
# Run your existing test suite (should pass identically)
pytest tests/  # or your test command
```

**2. Standard QA Test**
```python
# Verify standard QA unchanged
from neoswga.core.integrated_quality_scorer import create_quality_scorer

scorer = create_quality_scorer(stringency='moderate')  # No genome_gc
score = scorer.score_primer("ATGCATGCATGC")

assert 0.8 < score.overall_score < 0.9  # Typical range
```

**3. Adaptive QA Test**
```python
# Test adaptive QA on extreme genome
from neoswga.core.genome_analysis import analyze_genome_for_qa

# Use a genome you work with
analysis = analyze_genome_for_qa('your_genome.fna')

# Verify classification
assert analysis['gc_content'] >= 0 and analysis['gc_content'] <= 1
assert analysis['adaptive_qa_recommendation']['gc_class'] in [
    'extreme_at', 'at_rich', 'balanced', 'gc_rich', 'extreme_gc'
]
```

**4. Primer Generation Test**
```python
# Test full pipeline
from neoswga.core.optimal_oligo_generator import OptimalOligoGenerator

generator = OptimalOligoGenerator(genome_file='your_genome.fna')
result = generator.generate_primers(primer_count=10)

# Verify results
assert len(result.primers) > 0
assert result.reaction_conditions.genome_gc is not None
```

---

## Expected Behavior Changes

### What You'll Notice

**1. More Primers for Extreme Genomes**

Before (v3.4, standard QA):
```
Francisella (32% GC): 3 primers pass QA
```

After (v3.5, adaptive QA):
```
Francisella (32% GC): 8 primers pass QA (+167%)
```

**2. New Log Messages**

You'll see informative messages like:
```
INFO - Genome-adaptive QA enabled: genome_gc=32.3%, stringency=moderate,
       terminal_tm=8.0°C (vs 14.0°C standard)
INFO - Using genome-adaptive QA for 32.3% GC genome
```

**3. Enhanced Result Metadata**

Result objects now include:
```python
result.reaction_conditions.genome_gc  # 0.323
result.reaction_conditions.use_adaptive_qa  # True
```

### What Won't Change

- Standard QA behavior (genome_gc=None)
- Balanced genome results (40-60% GC)
- File formats
- Configuration files
- API function signatures
- Error handling

---

## Rollback Instructions

If you need to rollback to the previous version:

### Step 1: Reinstall Previous Version

```bash
# Check your backup
cat requirements_backup.txt | grep neoswga

# Reinstall specific version (replace X.X.X with your version)
pip install neoswga==3.4.0  # Example

# Verify version
python -c "import neoswga; print(neoswga.__version__)"
```

### Step 2: Restore Results (If Needed)

```bash
# If you need to restore old results
cp -r /path/to/results_backup/* /path/to/results/
```

### Step 3: Test Functionality

```bash
# Run your test suite
pytest tests/

# Or run basic functionality test
python test_basic.py
```

---

## FAQ

### Q: Do I need to change my existing code?

**A:** No. All existing code works unchanged. Adaptive QA only activates when you provide a genome file or explicitly pass `genome_gc` parameter.

```python
# Old code (v3.4) - still works identically
scorer = create_quality_scorer(stringency='moderate')

# New feature (v3.5) - opt-in only
scorer = create_quality_scorer(stringency='moderate', genome_gc=0.32)
```

---

### Q: Will my primers change?

**A:** Only if:
1. You provide a genome file (enables automatic adaptive QA), AND
2. Your genome is extreme (< 40% or > 60% GC)

For balanced genomes (40-60% GC), results will be nearly identical.

---

### Q: What if I want to keep using standard QA?

**A:** Simply don't provide the `genome_gc` parameter:

```python
# Forces standard QA (v3.4 behavior)
scorer = create_quality_scorer(stringency='moderate', genome_gc=None)
```

---

### Q: How do I know if adaptive QA is active?

**A:** Check the result object or logs:

```python
# Check result
if result.reaction_conditions.use_adaptive_qa:
    print("Adaptive QA is active")
    print(f"Genome GC: {result.reaction_conditions.genome_gc:.1%}")
```

Or watch for log messages:
```
INFO - Genome-adaptive QA enabled: genome_gc=68.2%...
```

---

### Q: Can I use adaptive QA for balanced genomes?

**A:** Yes, but benefit is minimal (+10-20%). Adaptive QA is most beneficial for extreme genomes.

---

### Q: Will this break my production pipeline?

**A:** No. The upgrade is designed to be transparent:
- Existing scripts run unchanged
- No configuration changes required
- Backward compatibility guaranteed
- You choose when to enable new features

---

### Q: How do I report issues?

**A:**
1. Check this migration guide and troubleshooting sections
2. Review the [User Guide](user-guide.md)
3. File an issue on GitHub with:
   - NeoSWGA version
   - Your genome GC content
   - Expected vs actual behavior
   - Minimal reproducible example

---

### Q: What if my genome is balanced (50% GC)?

**A:** Adaptive QA is optional for balanced genomes. The system will:
1. Analyze your genome
2. Note it's balanced
3. Optionally use adaptive QA (minimal difference from standard)
4. Work exactly as before

---

### Q: Can I disable adaptive QA globally?

**A:** Yes, by not providing genome files or genome_gc parameter. Your code will use standard QA as in v3.4.

---

## Summary Checklist

Before deploying to production:

- [ ] Backup current environment (`pip freeze`)
- [ ] Backup existing results
- [ ] Update to NeoSWGA 3.5 (`pip install --upgrade`)
- [ ] Test standard QA (backward compatibility)
- [ ] Test adaptive QA (new feature)
- [ ] Run full pipeline test
- [ ] Review log messages
- [ ] Compare results (optional)
- [ ] Update documentation/workflows if needed
- [ ] Train team on new features

---

## Support

**Documentation**:
- [User Guide](user-guide.md) - Usage tutorials and adaptive QA
- [API Reference](API_REFERENCE.md) - Complete API documentation
- [Module Reference](MODULE_REFERENCE.md) - All core modules

**Help**:
- GitHub Issues: Report bugs and request features
- Email: Contact development team
- Documentation: Check FAQ and troubleshooting sections

---

## Version History

**Version 3.5.0** (November 2025)
- Added genome-adaptive quality assessment
- New genome analysis functions
- Enhanced primer selection strategies
- 100% backward compatible with v3.4

**Previous Versions**
- v3.4.x: Standard QA only
- v3.3.x: Legacy system

---

**Migration Guide Version**: 1.0
**Last Updated**: November 2025
**Target Audience**: Existing NeoSWGA users (< v3.5)
