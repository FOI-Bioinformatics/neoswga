"""Small packaged genomes for `neoswga validate --smoke`.

A 6.2 kB plasmid as the target and a 6.3 kB one as the background, copied from
`examples/plasmid_example`. They live inside the package rather than being read
from `examples/` because that directory is not package data, so a smoke mode
reading from there would work in a source checkout and fail for every user who
installed from a wheel.
"""

from pathlib import Path

SMOKE_DIR = Path(__file__).parent
FOREGROUND = SMOKE_DIR / "pcDNA.fasta"
BACKGROUND = SMOKE_DIR / "pLTR.fasta"
