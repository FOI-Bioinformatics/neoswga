"""`step2_df.csv` and `step3_df.csv` survive an interrupted write.

Both were written with a plain `to_csv`, so a run killed mid-write left a short
file that still parses. Step 4 checks only that the file is non-empty, so it
accepted the truncated pool (2026-11-23 sweep, tier 2).
"""

import pandas as pd
import pytest

from neoswga.core import pipeline


def _interrupting_to_csv(self, path_or_buf=None, *args, **kwargs):
    handle = path_or_buf
    if isinstance(handle, str):
        handle = open(handle, "w")
    handle.write("primer\nPARTIAL\n")
    handle.flush()
    raise KeyboardInterrupt


@pytest.mark.parametrize(
    "writer,name",
    [(pipeline._write_step2_table, "step2_df.csv"), (pipeline._write_step3_table, "step3_df.csv")],
)
def test_an_interrupted_write_leaves_the_previous_table(tmp_path, monkeypatch, writer, name):
    target = tmp_path / name
    target.write_text("primer\nACGTACGTAC\nTTGACCATGA\n")
    monkeypatch.setattr(pd.DataFrame, "to_csv", _interrupting_to_csv)
    with pytest.raises(KeyboardInterrupt):
        writer(pd.DataFrame({"primer": ["GGCATTACGA"]}), str(tmp_path))
    assert target.read_text() == "primer\nACGTACGTAC\nTTGACCATGA\n"
