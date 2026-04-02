"""
Restricted pickle loading to mitigate arbitrary code execution risk.

Standard pickle.load() can execute arbitrary Python code embedded in
pickle files. This module provides a RestrictedUnpickler that only
allows classes from an explicit allowlist.

Usage:
    from neoswga.core.safe_pickle import safe_load

    model = safe_load(path, context='sklearn_model')
    bloom = safe_load(path, context='bloom_filter')
"""

import io
import pickle
import logging

logger = logging.getLogger(__name__)


# Allowlists per context: (module, class_name) tuples
_ALLOWED_CLASSES = {
    'sklearn_model': {
        ('sklearn.ensemble._forest', 'RandomForestClassifier'),
        ('sklearn.ensemble._forest', 'RandomForestRegressor'),
        ('sklearn.ensemble.forest', 'RandomForestClassifier'),
        ('sklearn.ensemble.forest', 'RandomForestRegressor'),
        ('sklearn.tree._tree', 'Tree'),
        ('sklearn.tree.tree', 'DecisionTreeClassifier'),
        ('sklearn.tree._classes', 'DecisionTreeClassifier'),
        ('sklearn.tree._classes', 'DecisionTreeRegressor'),
        ('numpy', 'ndarray'),
        ('numpy', 'dtype'),
        ('numpy.core.multiarray', 'scalar'),
        ('numpy.core.multiarray', '_reconstruct'),
        ('numpy._core.multiarray', 'scalar'),
        ('numpy._core.multiarray', '_reconstruct'),
        ('numpy', 'int64'),
        ('numpy', 'float64'),
        ('builtins', 'slice'),
        ('collections', 'defaultdict'),
    },
    'bloom_filter': {
        ('pybloom_live', 'ScalableBloomFilter'),
        ('pybloom_live', 'BloomFilter'),
        ('pybloom_live.pybloom', 'ScalableBloomFilter'),
        ('pybloom_live.pybloom', 'BloomFilter'),
        ('collections', 'defaultdict'),
        ('builtins', 'set'),
        ('numpy', 'ndarray'),
        ('numpy', 'dtype'),
        ('numpy.core.multiarray', '_reconstruct'),
    },
}


class RestrictedUnpickler(pickle.Unpickler):
    """Unpickler that only allows classes from an explicit allowlist."""

    def __init__(self, file, allowed_classes):
        super().__init__(file)
        self._allowed = allowed_classes

    def find_class(self, module, name):
        if (module, name) in self._allowed:
            return super().find_class(module, name)
        # Allow standard builtins needed by pickle protocol
        if module == 'builtins' and name in ('dict', 'list', 'tuple', 'set',
                                              'frozenset', 'int', 'float',
                                              'str', 'bytes', 'bool', 'None',
                                              'complex', 'type'):
            return super().find_class(module, name)
        # Allow copy_reg/_reconstructor (used by many objects)
        if module == 'copy_reg' and name == '_reconstructor':
            return super().find_class(module, name)
        if module == 'copyreg' and name == '_reconstructor':
            return super().find_class(module, name)
        raise pickle.UnpicklingError(
            f"Blocked unpickling of {module}.{name}. "
            f"If this class is expected, add it to the allowlist in safe_pickle.py."
        )


def safe_load(path, context='sklearn_model'):
    """Load a pickle file with restricted class allowlist.

    Args:
        path: Path to pickle file.
        context: One of 'sklearn_model' or 'bloom_filter'. Determines
            which classes are permitted during deserialization.

    Returns:
        The deserialized object.

    Raises:
        pickle.UnpicklingError: If the file tries to instantiate a
            class not in the allowlist.
        ValueError: If context is unknown.
    """
    if context not in _ALLOWED_CLASSES:
        raise ValueError(
            f"Unknown context '{context}'. "
            f"Available: {list(_ALLOWED_CLASSES.keys())}"
        )

    allowed = _ALLOWED_CLASSES[context]
    with open(path, 'rb') as f:
        return RestrictedUnpickler(f, allowed).load()


def unsafe_load_with_warning(path, purpose=''):
    """Load pickle with a logged warning. Fallback when safe_load fails.

    Use this only when RestrictedUnpickler cannot handle the file
    (e.g., complex sklearn internals). Logs a warning about the risk.
    """
    msg = f"Loading pickle file {path} without restriction"
    if purpose:
        msg += f" ({purpose})"
    msg += ". Only load files from trusted sources."
    logger.warning(msg)
    with open(path, 'rb') as f:
        return pickle.load(f)
