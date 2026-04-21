#!/usr/bin/env sh
# Fetch E. coli K-12 MG1655 (NC_000913.3) from NCBI as a single FASTA record.
# Requires: curl. No API key needed for single-record efetch at this scale.

set -eu

ACCESSION="NC_000913.3"
OUT="ecoli_K12_MG1655.fasta"
URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${ACCESSION}&rettype=fasta&retmode=text"

if [ -f "$OUT" ]; then
  echo "$OUT already exists; skipping download."
  exit 0
fi

echo "Fetching ${ACCESSION} from NCBI..."
curl -sS --fail "$URL" -o "$OUT"

# Basic sanity check: expect roughly 4.7 MB for the ~4.64 Mb chromosome plus headers.
SIZE=$(wc -c < "$OUT")
if [ "$SIZE" -lt 4000000 ]; then
  echo "Warning: downloaded file is only ${SIZE} bytes; expected ~4.7 MB. Check network or accession." >&2
  exit 1
fi

echo "Downloaded $OUT (${SIZE} bytes)"
