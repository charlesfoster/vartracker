#!/usr/bin/env bash
set -euo pipefail

# Reproduces the PAO1 bacterial-genome validation described in the vartracker
# README ("Bacterial genomes" section under Limitations) and manuscript: a
# synthetic 80-variant scenario and a synthetic 1000-variant scenario, each
# spanning 6 longitudinal timepoints, simulated against Pseudomonas aeruginosa
# PAO1 (NC_002516.2) and analysed with `vartracker vcf`.
#
# Usage:
#   scripts/validation/pao1/run_validation.sh [outdir]
#
# Requires on PATH: vartracker, bcftools, bgzip, python3 (with numpy and
# pandas). Optional, to reproduce the reported timing/peak-memory figures:
# GNU time - `gtime` on macOS (brew install gnu-time) or `/usr/bin/time -v`
# on Linux. Without either, the analysis still runs, just unbenchmarked.

OUTDIR="${1:-pao1_validation}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$OUTDIR"

echo "== Step 1: build the PAO1 reference bundle (NC_002516.2) =="
vartracker prepare reference \
  --accessions NC_002516.2 \
  --outdir "$OUTDIR/refs/pao1" \
  --prefix pao1_ref

echo "== Step 2: simulate scenario_a (80 variants, uniform gene weighting) and scenario_b (1000 variants, Zipf-biased gene weighting) =="
python3 "$SCRIPT_DIR/simulate_pao1.py" \
  --fasta "$OUTDIR/refs/pao1/pao1_ref.fa" \
  --gff3 "$OUTDIR/refs/pao1/pao1_ref.gff3" \
  --base-outdir "$OUTDIR"

if command -v gtime >/dev/null 2>&1; then
  TIME_CMD="gtime -v"
elif /usr/bin/time -v true >/dev/null 2>&1; then
  TIME_CMD="/usr/bin/time -v"
else
  echo "No GNU time found (gtime / time -v); continuing without wall-clock/memory benchmarking." >&2
  TIME_CMD=""
fi

echo "== Step 3: run vartracker vcf on scenario_a =="
$TIME_CMD vartracker vcf "$OUTDIR/scenario_a/scenario_a_input.csv" \
  --reference "$OUTDIR/refs/pao1/pao1_ref.fa" \
  --gff3 "$OUTDIR/refs/pao1/pao1_ref.gff3" \
  --outdir "$OUTDIR/results/pao1_80_variants"

echo "== Step 4: run vartracker vcf on scenario_b =="
$TIME_CMD vartracker vcf "$OUTDIR/scenario_b/scenario_b_input.csv" \
  --reference "$OUTDIR/refs/pao1/pao1_ref.fa" \
  --gff3 "$OUTDIR/refs/pao1/pao1_ref.gff3" \
  --outdir "$OUTDIR/results/pao1_1000_variants"

echo
echo "Done."
echo "Cross-check $OUTDIR/results/pao1_*_variants/results.csv against"
echo "$OUTDIR/scenario_{a,b}/scenario_{a,b}_ground_truth.csv, and (if GNU time"
echo "ran above) compare Elapsed / Maximum resident set size against the"
echo "figures quoted in the README's Bacterial genomes section."
