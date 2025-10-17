#!/usr/bin/env bash
set -euo pipefail

# Run any number of arbitray files; assumes some input parameters (see below)

# Usage:
#   ./metal_meta.sh results/meta_all data/AoUv8.txt.gz data/UKB.txt.gz data/MGB.txt.gz
# Produces: results/meta_all.TBL, results/meta_all.INFO, etc.

out_prefix="${1:?First arg must be output prefix (e.g., results/meta_all)}"
shift
files=("$@")
(( ${#files[@]} > 0 )) || { echo "Give at least one input file."; exit 1; }

mkdir -p "$(dirname "$out_prefix")"

# If METAL isn't on PATH, set the full path here:
METAL_BIN="${METAL_BIN:-METAL/build/bin/metal}"

# Build the METAL script and pipe it to METAL
{
  cat <<'HEAD'
SCHEME STDERR
AVERAGEFREQ ON
# Uncomment one of the next lines if needed based on your files:
# SEPARATOR TAB
# SEPARATOR WHITESPACE

CUSTOMVARIABLE CHROM
CUSTOMVARIABLE GENPOS
CUSTOMVARIABLE N
CUSTOMVARIABLE N_CASES
CUSTOMVARIABLE N_CONTROLS
CUSTOMVARIABLE Ncarriers

MARKER ID
ALLELE ALLELE1 ALLELE0
EFFECT BETA
STDERR SE
PVALUE P
FREQLABEL A1FREQ

LABEL CHROM as CHROM
LABEL GENPOS as GENPOS
LABEL N as N
LABEL N_CASES as N_CASES
LABEL N_CONTROLS as N_CONTROLS
LABEL Ncarriers as Ncarriers
HEAD

  # Output prefix (let METAL add .TBL, .INFO, etc.)
  printf 'OUTFILE %s .\n' "$out_prefix"

  # One PROCESS line per file
  for f in "${files[@]}"; do
    # %q safely quotes paths (handles spaces)
    printf 'PROCESS %q\n' "$f"
  done

  cat <<'TAIL'
ANALYZE
# Optionally:
# ANALYZE HETEROGENEITY
QUIT
TAIL
} | "$METAL_BIN"
