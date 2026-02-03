#!/usr/bin/env python3
"""
Retrain random forest model for sklearn 1.7 compatibility.

This script creates a new random forest model compatible with sklearn 1.7+
by generating synthetic training data based on primer characteristics and
known thermodynamic principles from SWGA literature.

The model predicts amplification scores based on:
- Sequence composition features (GC content, base frequencies)
- Melting temperature (Tm)
- 3' end stability (GC clamp)
- Homopolymer runs (repeat regions)
- Thermodynamic binding energy histogram

Usage:
    python scripts/retrain_rf_model.py --output neoswga/core/models/random_forest_filter.p
"""

import argparse
import os
import pickle
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor


# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

# Feature names matching rf_preprocessing.py
BASE_FEATURES = [
    'molarity', 'sequence.length', 'number.of.A', 'proportion.of.A',
    'number.of.T', 'proportion.of.T', 'number.of.G', 'proportion.of.G',
    'number.of.C', 'proportion.of.C', 'GC.content', 'melting_tm',
    'GC.clamp', 'longest.A.repeat', 'longest.T.repeat', 'longest.G.repeat',
    'longest.C.repeat', 'AA repeat', 'CC repeat', 'TT repeat', 'GG repeat',
    '3.end.first.base', '3.end.second.base', '3.end.third.base',
    '3.end.fourth.base', '3.end.fifth.base'
]

# Thermodynamic histogram bins
BINS = [-20, -18, -16, -14, -12, -10, -9, -8, -7, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
DELTA_G_FEATURES = ['on_target_' + str(b) for b in BINS[1:]]

ALL_FEATURES = BASE_FEATURES + DELTA_G_FEATURES


def generate_primer_sequence(length=12, gc_target=None):
    """Generate a random primer sequence with optional GC content target."""
    bases = 'ATGC'
    if gc_target is None:
        gc_target = random.uniform(0.3, 0.7)

    gc_count = int(length * gc_target)
    at_count = length - gc_count

    seq = ['G'] * (gc_count // 2) + ['C'] * (gc_count - gc_count // 2)
    seq += ['A'] * (at_count // 2) + ['T'] * (at_count - at_count // 2)
    random.shuffle(seq)
    return ''.join(seq)


def compute_melting_temp(seq):
    """Simple Tm calculation using basic nearest-neighbor approximation."""
    # Wallace rule for short primers
    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')

    if len(seq) < 14:
        # Wallace rule
        tm = 2 * (a_count + t_count) + 4 * (g_count + c_count)
    else:
        # Marmur-Doty approximation
        gc_content = (g_count + c_count) / len(seq)
        tm = 64.9 + 41 * (gc_content - 0.16)

    return tm


def extract_features(seq, molarity=2.5):
    """Extract features from a primer sequence."""
    features = {}

    # Basic features
    features['molarity'] = molarity
    features['sequence.length'] = len(seq)

    # Base counts and proportions
    for base in 'ATGC':
        count = seq.count(base)
        features[f'number.of.{base}'] = count
        features[f'proportion.of.{base}'] = count / len(seq)

    # GC content
    features['GC.content'] = (seq.count('G') + seq.count('C')) / len(seq)

    # Melting temperature
    features['melting_tm'] = compute_melting_temp(seq)

    # GC clamp (last 5 bases)
    last_five = seq[-5:]
    features['GC.clamp'] = last_five.count('G') + last_five.count('C')

    # Longest repeats
    for base in 'ATGC':
        max_run = 0
        current_run = 0
        for c in seq:
            if c == base:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        features[f'longest.{base}.repeat'] = max_run

    # Dinucleotide repeats
    for dinuc in ['AA', 'CC', 'TT', 'GG']:
        features[f'{dinuc} repeat'] = seq.count(dinuc)

    # 3' end bases (encoded as integers)
    base_to_int = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    for i, name in enumerate(['first', 'second', 'third', 'fourth', 'fifth']):
        if len(seq) > i:
            features[f'3.end.{name}.base'] = base_to_int.get(seq[-(i+1)], 0)
        else:
            features[f'3.end.{name}.base'] = 0

    # Generate synthetic thermodynamic histogram
    # Primers with good Tm and GC content should have more strong binding events
    total_sites = random.randint(100, 10000)
    gc = features['GC.content']
    tm = features['melting_tm']

    # Quality factor based on primer properties
    quality = 0.5
    if 0.4 <= gc <= 0.6:
        quality += 0.2
    if 30 <= tm <= 42:
        quality += 0.2
    if features['GC.clamp'] >= 2:
        quality += 0.1

    # Distribute binding sites across energy bins
    # Better primers have more in low-energy (strong binding) bins
    weights = np.array([quality ** (i * 0.3) for i in range(len(DELTA_G_FEATURES))])
    weights = weights / weights.sum()
    site_distribution = np.random.multinomial(total_sites, weights)

    for feat, count in zip(DELTA_G_FEATURES, site_distribution):
        features[feat] = count

    return features


def compute_target_score(features):
    """
    Compute target amplification score based on primer features.

    This follows known SWGA principles:
    - Optimal Tm for MDA/EquiPhi29: 30-45C
    - Optimal GC content: 40-60%
    - GC clamp at 3' end improves binding
    - Long homopolymer runs reduce specificity
    - Strong binding (low delta-G) improves amplification
    """
    score = 10.0  # baseline

    # Melting temperature effect
    tm = features['melting_tm']
    if 30 <= tm <= 42:
        score += 3.0
    elif 25 <= tm < 30 or 42 < tm <= 50:
        score += 1.5
    elif tm < 20 or tm > 55:
        score -= 3.0

    # GC content effect
    gc = features['GC.content']
    if 0.4 <= gc <= 0.6:
        score += 2.5
    elif 0.35 <= gc < 0.4 or 0.6 < gc <= 0.65:
        score += 1.0
    elif gc < 0.25 or gc > 0.75:
        score -= 3.0

    # GC clamp effect
    gc_clamp = features['GC.clamp']
    if 2 <= gc_clamp <= 3:
        score += 2.0
    elif gc_clamp == 1 or gc_clamp == 4:
        score += 0.5
    elif gc_clamp > 4:
        score -= 1.0

    # Homopolymer penalty
    max_repeat = max(
        features['longest.A.repeat'],
        features['longest.T.repeat'],
        features['longest.G.repeat'],
        features['longest.C.repeat']
    )
    if max_repeat >= 5:
        score -= 4.0
    elif max_repeat >= 4:
        score -= 2.0
    elif max_repeat >= 3:
        score -= 1.0

    # 3' end stability (G or C at 3' end)
    end_base = features['3.end.first.base']
    if end_base in [2, 3]:  # G or C
        score += 1.0

    # Thermodynamic binding quality
    strong_binding = sum(features.get(f, 0) for f in DELTA_G_FEATURES[:8])
    weak_binding = sum(features.get(f, 0) for f in DELTA_G_FEATURES[-8:])
    total_binding = sum(features.get(f, 0) for f in DELTA_G_FEATURES)

    if total_binding > 0:
        strong_ratio = strong_binding / total_binding
        score += strong_ratio * 4.0

        # Penalty for too much weak binding
        weak_ratio = weak_binding / total_binding
        score -= weak_ratio * 2.0

    # Add some noise
    score += random.gauss(0, 0.5)

    # Clamp to reasonable range [0, 20]
    return max(0.0, min(20.0, score))


def generate_training_data(n_samples=5000):
    """Generate synthetic training data."""
    print(f"Generating {n_samples} training samples...")

    data = []
    for i in range(n_samples):
        # Vary primer length (6-18 bp)
        length = random.choice([6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18])

        # Generate primer with varying GC content
        gc_target = random.uniform(0.25, 0.75)
        seq = generate_primer_sequence(length, gc_target)

        # Extract features
        features = extract_features(seq)

        # Compute target score
        score = compute_target_score(features)

        features['sequence'] = seq
        features['target_score'] = score
        data.append(features)

        if (i + 1) % 1000 == 0:
            print(f"  Generated {i + 1}/{n_samples} samples")

    return pd.DataFrame(data)


def train_model(df, output_path):
    """Train and save the random forest model."""
    print("\nTraining random forest model...")

    # Prepare features and target
    X = df[ALL_FEATURES]
    y = df['target_score']

    # Handle any NaN values
    X = X.fillna(0)

    # Train random forest
    rf = RandomForestRegressor(
        n_estimators=100,
        max_depth=15,
        min_samples_split=5,
        min_samples_leaf=2,
        n_jobs=-1,
        random_state=42
    )

    rf.fit(X, y)

    # Evaluate on training data
    y_pred = rf.predict(X)
    mse = np.mean((y - y_pred) ** 2)
    r2 = 1 - mse / np.var(y)

    print(f"  Training MSE: {mse:.4f}")
    print(f"  Training R2: {r2:.4f}")

    # Save model
    print(f"\nSaving model to {output_path}...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, 'wb') as f:
        pickle.dump(rf, f)

    print("Model saved successfully!")

    # Verify model can be loaded
    with open(output_path, 'rb') as f:
        loaded_model = pickle.load(f)

    test_pred = loaded_model.predict(X.head(5))
    print(f"\nVerification: Model loaded and predicted {len(test_pred)} samples")

    return rf


def main():
    parser = argparse.ArgumentParser(
        description='Retrain random forest model for sklearn 1.7 compatibility'
    )
    parser.add_argument('--output', '-o', type=str,
                        default='neoswga/core/models/random_forest_filter.p',
                        help='Output model path')
    parser.add_argument('--samples', '-n', type=int, default=5000,
                        help='Number of training samples')
    parser.add_argument('--backup', action='store_true',
                        help='Backup existing model before overwriting')

    args = parser.parse_args()

    # Resolve output path
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
    output_path = project_dir / args.output

    # Backup existing model if requested
    if args.backup and output_path.exists():
        backup_path = output_path.with_suffix('.p.bak')
        print(f"Backing up existing model to {backup_path}")
        import shutil
        shutil.copy(output_path, backup_path)

    # Generate training data
    df = generate_training_data(args.samples)

    # Train and save model
    train_model(df, str(output_path))

    print("\nDone! The new model is compatible with sklearn 1.7+")


if __name__ == '__main__':
    main()
