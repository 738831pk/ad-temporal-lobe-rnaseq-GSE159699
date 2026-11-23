#!/usr/bin/env bash
# Downloads GSE159699 files from NCBI GEO.
# curl with a browser-like user agent is used because R/GEOquery requests get HTTP 403.
set -euo pipefail
mkdir -p data/GSE159699
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159699"
curl -fL -A "Mozilla/5.0" -o data/GSE159699/GSE159699_summary_count.star.txt.gz "$BASE/suppl/GSE159699_summary_count.star.txt.gz"
curl -fL -A "Mozilla/5.0" -o data/GSE159699/GSE159699_series_matrix.txt.gz "$BASE/matrix/GSE159699_series_matrix.txt.gz"
ls -lh data/GSE159699
