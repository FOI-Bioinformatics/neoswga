"""Small shared I/O helpers for pipeline dataframes.

The core pipeline writes its primer column as ``'primer'`` (step2/step3/step4
CSVs), while the optional QA-integration path historically expected ``'seq'``.
`primer_column` resolves the difference in one place so call sites do not
re-implement the fallback.
"""

from __future__ import annotations


def primer_column(df) -> str:
    """Return the name of the primer-sequence column in a pipeline dataframe.

    Prefers ``'primer'`` (the canonical core-pipeline name), falling back to
    ``'seq'`` (the QA-integration name).

    Args:
        df: a pandas DataFrame.

    Returns:
        ``'primer'`` or ``'seq'``.

    Raises:
        KeyError: if neither column is present.
    """
    if "primer" in df.columns:
        return "primer"
    if "seq" in df.columns:
        return "seq"
    raise KeyError(
        "DataFrame has neither a 'primer' nor a 'seq' column; "
        f"columns present: {list(df.columns)}"
    )
