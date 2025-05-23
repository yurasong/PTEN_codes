#!/usr/bin/env bash
################################################################################
# 09_01_ATAC_alignment.sh
#
# Non-interactive ATAC-seq alignment pipeline (paired-end) to fixed mm10 genome.
# Steps:
#   1) Trim adapters (Trimmomatic)
#   2) Align to mm10 (bowtie2)
#   3) Filter mito & unwanted contigs
#   4) Tn5 shift correction
#   5) Mark duplicates (Picard)
#   6) Index & flagstat (samtools)
#   7) Peak calling (MACS2)
#   8) Remove blacklist regions (bedtools)
#
# Usage:
#   ./09_01_ATAC_alignment.sh <FASTQ_R1> <FASTQ_R2> <SAMPLE_NAME>
################################################################################

set -euo pipefail

# --- Parse arguments ---
if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <FASTQ_R1(.fq|.fq.gz)> <FASTQ_R2(.fq|.fq.gz)> <SAMPLE_NAME>"
  exit 1
fi
FASTQ_R1="$1"
FASTQ_R2="$2"
SAMPLE="$3"

# --- Configuration ---
THREADS=4
MAX_INSERT=2000
MINLEN=36
HEADCROP=8
SLIDINGWINDOW="4:15"
LEADING=3
TRAILING=3

# Fixed genome prefix and paths
GENOME="mm10"
GENOME_DIR="reference_genomes/\$GENOME"
BWT2_INDEX="\$GENOME_DIR/\$GENOME"
ADAPTERS="reference_genomes/adapters/Nextera1.fa"
BLACKLIST_BED="\$GENOME_DIR/\${GENOME}-blacklist.bed"
PICARD_JAR="\${PICARD_JAR:-/usr/local/bin/picard.jar}"

# --- Helper: decompress if gzipped ---
decompress_if_gz() {
  local fq="\$1"
  if [[ "\$fq" == *.gz ]]; then
    gunzip -c "\$fq"
  else
    cat "\$fq"
  fi
}

# --- 1) Trim adapters ---
TRIM1="\${SAMPLE}_trim_R1.fq"
TRIM2="\${SAMPLE}_trim_R2.fq"
UNTR1="\${SAMPLE}_untrim_R1.fq"
UNTR2="\${SAMPLE}_untrim_R2.fq"

echo "[1/8] Trimming adapters..."
decompress_if_gz "\$FASTQ_R1" > tmp_R1.fq
decompress_if_gz "\$FASTQ_R2" > tmp_R2.fq
TrimmomaticPE \
  -threads \$THREADS -phred33 \
  tmp_R1.fq tmp_R2.fq \
  "\$TRIM1" "\$UNTR1" \
  "\$TRIM2" "\$UNTR2" \
  ILLUMINACLIP:"\$ADAPTERS":2:30:10 \
  LEADING:\$LEADING TRAILING:\$TRAILING \
  SLIDINGWINDOW:\$SLIDINGWINDOW MINLEN:\$MINLEN \
  HEADCROP:\$HEADCROP
rm tmp_R1.fq tmp_R2.fq

# --- 2) Align with bowtie2 ---
SAM_RAW="\${SAMPLE}.sam"
echo "[2/8] Aligning to \$GENOME..."
bowtie2 \
  -x "\$BWT2_INDEX" \
  -1 "\$TRIM1" -2 "\$TRIM2" \
  -p \$THREADS \
  -X \$MAX_INSERT \
  --very-sensitive --no-discordant --no-mixed --no-unal \
  -S "\$SAM_RAW"

# --- 3) Filter mito & unwanted contigs ---
SAM_NO_MT="\${SAMPLE}_noMT.sam"
echo "[3/8] Filtering mitochondrial & unwanted contigs..."
sed '/^@/b; /MT\|chrM\|random\|Un\|hap\|KI\|GL/d' "\$SAM_RAW" > "\$SAM_NO_MT"
rm "\$SAM_RAW"

# --- 4) Convert to BAM & sort ---
BAM_NO_MT="\${SAMPLE}_noMT.bam"
echo "[4/8] Converting to BAM & sorting..."
samtools view -@ \$THREADS -b "\$SAM_NO_MT" \
  | samtools sort -@ \$THREADS -m 2G -o "\$BAM_NO_MT"
rm "\$SAM_NO_MT"

# --- 5) Tn5 shift correction ---
BAM_SHIFTED="\${SAMPLE}_shifted.bam"
echo "[5/8] Applying Tn5 shift correction..."
samtools view -h "\$BAM_NO_MT" \
  | awk 'BEGIN{OFS="\t"} /^@/ {print; next} { if(and(\$2,16)) {\$4=\$4-5} else {\$4=\$4+4} }1' \
  | samtools view -@ \$THREADS -b - \
  | samtools sort -@ \$THREADS -m 2G -o "\$BAM_SHIFTED"
rm "\$BAM_NO_MT"

# --- 6) Mark duplicates ---
BAM_NODUP="\${SAMPLE}_nodup.bam"
METRICS="\${SAMPLE}_dup_metrics.txt"
echo "[6/8] Marking duplicates..."
java -jar "\$PICARD_JAR" MarkDuplicates \
  I="\$BAM_SHIFTED" O="\$BAM_NODUP" M="\$METRICS" \
  REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
rm "\$BAM_SHIFTED"

# --- 7) Index & flagstat ---
echo "[7/8] Indexing & collecting stats..."
samtools index -@ \$THREADS "\$BAM_NODUP"
samtools flagstat -@ \$THREADS "\$BAM_NODUP" > "\${SAMPLE}_flagstat.txt"

# --- 8) Peak calling & blacklist removal ---
PEAK_DIR="\${SAMPLE}_peaks"
echo "[8/8] Calling peaks & filtering blacklist..."
mkdir -p "\$PEAK_DIR"
macs2 callpeak -t "\$BAM_NODUP" -f BAMPE -g mm -n "\$SAMPLE" -q 0.01 --nomodel --shift 0 --outdir "\$PEAK_DIR"
bedtools intersect -v -a "\$PEAK_DIR/\${SAMPLE}_peaks.narrowPeak" -b "\$BLACKLIST_BED" > "\$PEAK_DIR/\${SAMPLE}_final_peaks.bed"

# Done
cat <<EOF
Pipeline complete.
Outputs:
  - Trimmed FASTQ: $TRIM1, $TRIM2
  - BAM (no dup): $BAM_NODUP
  - Dup metrics:  $METRICS
  - Flagstat:     ${SAMPLE}_flagstat.txt
  - Peaks dir:    $PEAK_DIR
EOF
