#!/usr/bin/env bash
################################################################################
# 09_03_motif_prediction.sh
#
# Non-interactive wrapper for HOMER's findMotifsGenome.pl
# Usage:
#   ./09_03_motif_prediction.sh \<peaks.bed> \<background.bed> \<output_dir> [options]
#
# Required arguments:
#   PEAKS_BED      BED file of target peak regions
#   BG_BED         BED file of background regions
#   OUT_DIR        Directory to write motif results
#
# Optional environment vars / flags:
#   GENOME         Genome prefix (default: mm10)
#   THREADS        Number of threads for motif scanning (default: 15)
#   LENGTHS        Comma-separated motif lengths (default: 6,8,10,12)
#   SIZE_OPTION    Size window around peak centers (default: -250,250)
################################################################################
set -euo pipefail

# --- Defaults ---
: "${GENOME:=mm10}"
: "${THREADS:=15}"
: "${LENGTHS:=6,8,10,12}"
: "${SIZE_OPTION:=-250,250}"

usage() {
  cat <<EOF
Usage: $0 <peaks.bed> <background.bed> <out_dir>

Environment variables:
  GENOME       genome prefix (default: \$GENOME)
  THREADS      threads for motif discovery (default: \$THREADS)
  LENGTHS      motif lengths, comma-separated (default: \$LENGTHS)
  SIZE_OPTION  size window around peaks (default: \$SIZE_OPTION)
EOF
  exit 1
}

# --- Parse args ---
if [[ $# -ne 3 ]]; then
  usage
fi
PEAKS_BED="$1"
BG_BED="$2"
OUT_DIR="$3"

# --- Prepare output directory ---
mkdir -p "$OUT_DIR"

echo "Running motif discovery with HOMER:"
echo "  Peaks:      $PEAKS_BED"
echo "  Background: $BG_BED"
echo "  Genome:     $GENOME"
echo "  Threads:    $THREADS"
echo "  Lengths:    $LENGTHS"
echo "  Size:       $SIZE_OPTION"
echo "  Output:     $OUT_DIR"

# --- Execute findMotifsGenome.pl ---
findMotifsGenome.pl \
  "$PEAKS_BED" \
  "$GENOME" \
  "$OUT_DIR" \
  -bg "$BG_BED" \
  -S "$THREADS" \
  -len "$LENGTHS" \
  -size "$SIZE_OPTION"

# --- Completion message ---
echo "Motif discovery complete. Results in $OUT_DIR."