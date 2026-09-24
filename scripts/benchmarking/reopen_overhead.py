"""What reopening the BAM once per record costs.

`bam_depth_profile` calls `compute_bam_depth` inside its loop over bound
records, and each call opens the file and loads the index. This measures the
toll at the *Drosophila* record count against one shared open.

The records here are tiny, so the ratio is inflated by construction and only
the ABSOLUTE overhead transfers: on real records the counting work dominates
and the ratio collapses while the seconds stay.

    python reopen_overhead.py [n_records] [record_length]
"""
import os
import sys
import tempfile
import time

import pysam

from neoswga.core.depth_policy import DepthPolicy


def main(n_records: int = 1870, record_length: int = 1000) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "many.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": f"ctg{i}", "LN": record_length} for i in range(n_records)]}
        with pysam.AlignmentFile(path, "wb", header=header) as out:
            for i in range(n_records):
                read = pysam.AlignedSegment()
                read.query_name = f"r{i}"
                read.flag = 0
                read.reference_id = i
                read.reference_start = 10
                read.mapping_quality = 60
                read.cigarstring = "50M"
                read.query_sequence = "A" * 50
                read.query_qualities = pysam.qualitystring_to_array("I" * 50)
                out.write(read)
        pysam.index(path)

        accepts = DepthPolicy().accepts
        names = [f"ctg{i}" for i in range(n_records)]

        start = time.perf_counter()
        for name in names:
            with pysam.AlignmentFile(path, "rb") as handle:
                handle.count_coverage(name, 0, record_length,
                                      quality_threshold=0, read_callback=accepts)
        per_record = time.perf_counter() - start

        start = time.perf_counter()
        with pysam.AlignmentFile(path, "rb") as handle:
            for name in names:
                handle.count_coverage(name, 0, record_length,
                                      quality_threshold=0, read_callback=accepts)
        shared = time.perf_counter() - start

        print(f"{n_records} records, {record_length} bp each")
        print(f"  open per record : {per_record:7.3f} s")
        print(f"  one shared open : {shared:7.3f} s")
        print(f"  overhead        : {per_record - shared:7.3f} s "
              f"({(per_record - shared) / n_records * 1000:.2f} ms per open)")


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]]
    sys.exit(main(*args))
