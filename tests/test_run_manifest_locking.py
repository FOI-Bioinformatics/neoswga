"""The manifest append is read-modify-write, so it needs a lock.

Two runs sharing a data_dir -- an ensemble sweep, a batch of designs against
one background, two terminals -- both read the same list and both write their
own append, and one entry is lost. fcntl.flock is already used correctly in
experimental_tracker.py.
"""

import json
import os
import subprocess
import sys

from neoswga.core import run_manifest as rm


def test_fcntl_is_available_and_used():
    assert rm.HAS_FCNTL is True or os.name == "nt", (
        "fcntl should be available on this platform; the manifest append is "
        "a read-modify-write and needs it"
    )


def test_concurrent_appends_keep_every_entry(tmp_path):
    """Twelve writers, four processes. Without a lock this loses entries."""
    writer = tmp_path / "writer.py"
    writer.write_text(
        "import sys\n"
        "from neoswga.core import run_manifest as rm\n"
        "data_dir, tag = sys.argv[1], sys.argv[2]\n"
        "for i in range(3):\n"
        "    rm.write_manifest(step=f'{tag}-{i}', data_dir=data_dir)\n"
    )

    procs = [
        subprocess.Popen(
            [sys.executable, str(writer), str(tmp_path), f"w{n}"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        for n in range(4)
    ]
    for p in procs:
        out, err = p.communicate(timeout=180)
        assert p.returncode == 0, err.decode()

    entries = json.loads((tmp_path / rm.MANIFEST_FILENAME).read_text())["steps"]
    names = sorted(e["step"] for e in entries)
    expected = sorted(f"w{n}-{i}" for n in range(4) for i in range(3))
    assert (
        names == expected
    ), f"lost or duplicated entries: got {len(names)}, expected {len(expected)}"


def test_a_single_write_still_produces_valid_json(tmp_path):
    rm.write_manifest(step="filter", data_dir=str(tmp_path))
    rm.write_manifest(step="score", data_dir=str(tmp_path))

    data = json.loads((tmp_path / rm.MANIFEST_FILENAME).read_text())
    assert [e["step"] for e in data["steps"]] == ["filter", "score"]


def test_a_corrupt_manifest_is_replaced_not_appended_to(tmp_path):
    """The locked path reads the file it holds open rather than re-opening it.
    A truncated or hand-edited manifest must still leave a valid file behind."""
    (tmp_path / rm.MANIFEST_FILENAME).write_text('{"steps": [{"step": "filt')

    rm.write_manifest(step="score", data_dir=str(tmp_path))

    data = json.loads((tmp_path / rm.MANIFEST_FILENAME).read_text())
    assert [e["step"] for e in data["steps"]] == ["score"]
