"""Peak RSS of `compute_bam_depth`, one contig length per process.

ru_maxrss is a high-water mark for the PROCESS, so it never falls and a second
size measured in the same process reports the first size's peak. Each point
must therefore be its own process; this script measures one and prints a row.

    for L in 5000000 20000000 40000000; do python count_coverage_rss.py $L; done
"""
import gc
import os
import resource
import sys
import tempfile

import pysam


def measure(length: int) -> None:
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "one.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "big", "LN": length}]}
        with pysam.AlignmentFile(path, "wb", header=header) as out:
            read = pysam.AlignedSegment()
            read.query_name = "r"
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 10
            read.mapping_quality = 60
            read.cigarstring = "50M"
            read.query_sequence = "A" * 50
            read.query_qualities = pysam.qualitystring_to_array("I" * 50)
            out.write(read)
        pysam.index(path)

        from neoswga.core.bam_coverage import compute_bam_depth

        gc.collect()
        base = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        depth = compute_bam_depth(path, "big", length)
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

        delta = peak - base  # bytes on macOS; KiB on Linux -- scale there
        print(f"{length:>12,} bp   peak delta {delta:>15,} B   "
              f"{delta / length:>6.1f} B/base   result {depth.nbytes / length:.1f} B/base")


if __name__ == "__main__":
    measure(int(sys.argv[1]))
