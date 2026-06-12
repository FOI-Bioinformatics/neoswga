#!/usr/bin/env python3
"""Convert the bundled random-forest pickle to the version-tolerant skops format.

One-time migration helper (kept for reproducibility). Loads the legacy
``random_forest_filter.p`` and writes ``random_forest_filter.skops``, verifies
the predictions are identical, prints the new SHA-256 for ``checksums.json``,
and confirms the skops archive contains no untrusted types (so it loads with
skops' secure default ``trusted=False``).

Usage:
    python scripts/convert_model_to_skops.py
"""

import hashlib
import os
import sys

import numpy as np

MODELS = os.path.join(os.path.dirname(__file__), "..", "neoswga", "core", "models")
PKL = os.path.abspath(os.path.join(MODELS, "random_forest_filter.p"))
SKOPS = os.path.abspath(os.path.join(MODELS, "random_forest_filter.skops"))


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    if not os.path.exists(PKL):
        sys.exit(f"Legacy pickle not found: {PKL}")

    from skops.io import dump, get_untrusted_types, load

    from neoswga.core.rf_preprocessing import load_model_safely

    clf = load_model_safely(PKL, verify_hash=False)
    dump(clf, SKOPS)

    untrusted = get_untrusted_types(file=SKOPS)
    if untrusted:
        sys.exit(f"Refusing to ship: skops archive has untrusted types {untrusted}")

    clf2 = load(SKOPS)  # secure default load
    rng = np.random.RandomState(0)
    X = rng.rand(50, clf.n_features_in_)
    if not np.allclose(clf.predict(X), clf2.predict(X)):
        sys.exit("Predictions differ between pickle and skops model!")

    print(f"Wrote {SKOPS} ({os.path.getsize(SKOPS) // 1024} KB)")
    print("Predictions identical; no untrusted types.")
    print(f'  "random_forest_filter.skops": "{_sha256(SKOPS)}"')


if __name__ == "__main__":
    main()
