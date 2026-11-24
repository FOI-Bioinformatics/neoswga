"""`--no-modifications` with `--pto-bonds` refuses.

`run_export` took `--no-modifications` and dropped `--pto-bonds` without a
word, so an order asked for with three phosphorothioate bonds went out with
none (2026-11-23 sweep, tier 3).
"""

import sys

import pytest

from neoswga import cli_unified


def test_the_two_flags_together_refuse_before_writing(tmp_path, monkeypatch, capsys, caplog):
    out = tmp_path / "order"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "neoswga",
            "export",
            "-d",
            str(tmp_path),
            "-o",
            str(out),
            "--no-modifications",
            "--pto-bonds",
            "3",
        ],
    )
    with pytest.raises(SystemExit) as info:
        cli_unified.main()
    assert info.value.code != 0
    captured = capsys.readouterr()
    text = captured.out + captured.err + caplog.text
    assert "--no-modifications" in text and "--pto-bonds" in text
    assert not out.exists()
