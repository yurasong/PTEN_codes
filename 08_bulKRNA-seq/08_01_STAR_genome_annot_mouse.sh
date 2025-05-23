#!/usr/bin/env bash

#######################################################################################
# 08_01_STAR_genome_annot_mouse.sh
# Generate a STAR genome index for mouse (GRCm38)
# Run in a bash environment.
# Only needed once per genome build.
# Reference FASTA and GTF from Ensembl (GRCm38.87).
#######################################################################################

set -euo pipefail

# === Configuration ===
THREADS=4
GENOME_DIR="reference_genomes/RNA-seq/GRCm38"
GENOME_FASTA="$GENOME_DIR/GRCm38.fa"
GTF_FILE="$GENOME_DIR/Mus_musculus.GRCm38.87.gtf"

# === Run STAR genomeGenerate ===
STAR \
  --runThreadN      $THREADS \
  --runMode         genomeGenerate \
  --genomeDir       "$GENOME_DIR" \
  --genomeFastaFiles "$GENOME_FASTA" \
  --sjdbGTFfile     "$GTF_FILE" \
  --sjdbOverhang    $SJDB_OVERHANG
