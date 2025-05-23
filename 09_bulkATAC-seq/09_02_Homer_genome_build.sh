#!/usr/bin/env bash
################################################################################
# 09_02_Homer_genome_build.sh
#
# Wrapper to install a HOMER genome if not already present.
# Usage:
#   ./09_02_Homer_genome_build.sh [GENOME]
#
# Defaults:
#   GENOME = mm10
#
# Requires:
#   - HOMER installed and in $PATH (configureHomer.pl available)
################################################################################
set -euo pipefail

# --- Usage ---
usage() {
  cat <<EOF
Usage: $0 [GENOME]

Installs a HOMER genome index if not already present.
GENOME: Genome prefix (default: mm10)
Examples:
  $0 mm10   # installs mm10
  $0 hg38   # installs hg38
EOF
  exit 1
}

# --- Parse args ---
if [[ $# -gt 1 ]]; then
  usage
fi
GENOME=${1:-mm10}

# --- Check for configureHomer.pl ---
if ! command -v configureHomer.pl &>/dev/null; then
  echo "Error: configureHomer.pl not found in PATH. Please install HOMER and retry." >&2
  exit 1
fi

# --- Install genome if missing ---
# HOMER genomes are installed under $(homer_path)/data/genomes/
HOMER_BASE=$(dirname "$(command -v configureHomer.pl)")
GENOME_DIR="$HOMER_BASE/../data/genomes/$GENOME"

if [[ -d "$GENOME_DIR" ]]; then
  echo "Genome '$GENOME' is already installed at $GENOME_DIR. Skipping."
else
  echo "Installing HOMER genome: $GENOME"
  configureHomer.pl -install "$GENOME"
  echo "Installation complete."
fi
