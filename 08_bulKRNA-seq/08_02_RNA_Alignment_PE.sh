#!/usr/bin/env bash
################################################################################
# 08_02_RNA_Alignment_PE.sh
#
# Paired‐end RNA-seq alignment pipeline with STAR + Trimmomatic + HTSeq.
# Usage: ./08_02_RNA_Alignment_PE.sh <FASTQ_R1> <FASTQ_R2> <SAMPLE>
#
# Requires:
#   - reference_genomes/RNA-seq/GRCm38/ (STAR index built, GRCm38.fa, GTF)
#   - reference_genomes/adapters/TruSeq3-PE-2.fa
#   - reference_genomes/features/Mus_musculus.GRCm38.87.gtf
################################################################################
set -euo pipefail

usage() {
  echo "Usage: $0 <FASTQ_R1(.fq/.fq.gz)> <FASTQ_R2(.fq/.fq.gz)> <SAMPLE_NAME>"
  exit 1
}

#— Parse args
if [[ $# -ne 3 ]]; then
  usage
fi
FASTQ_R1="$1"
FASTQ_R2="$2"
SAMPLE="$3"

#— Config
THREADS=6 # Please adjust it
GENOME_DIR="reference_genomes/RNA-seq/GRCm38"
ADAPTERS="reference_genomes/adapters/TruSeq3-PE-2.fa"
GTF_FILE="reference_genomes/features/Mus_musculus.GRCm38.87.gtf"
SJDB_OVERHANG=98  # e.g. readLength-1 for 99bp reads

# Derived names
BASE1="${FASTQ_R1%%.gz}" && BASE1="${BASE1%%.fq}"
BASE2="${FASTQ_R2%%.gz}" && BASE2="${BASE2%%.fq}"
TRIM1="${BASE1}.trim.fq"
UNTRIM1="${BASE1}.untrim.fq"
TRIM2="${BASE2}.trim.fq"
UNTRIM2="${BASE2}.untrim.fq"

#— Decompress if needed
for fq in "$FASTQ_R1" "$FASTQ_R2"; do
  if [[ "$fq" == *.gz ]]; then
    echo "Decompressing $fq..."
    gunzip -c "$fq" > "${fq%.gz}"
  fi
done

#— 1) Trim adapters & low-quality bases
echo "Running Trimmomatic..."
java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" \
  PE -threads $THREADS -phred33 \
  "${FASTQ_R1%.gz}" "${FASTQ_R2%.gz}" \
  "$TRIM1" "$UNTRIM1" \
  "$TRIM2" "$UNTRIM2" \
  ILLUMINACLIP:"$ADAPTERS":2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:8

#— 2) STAR alignment
echo "Running STAR..."
STAR \
  --runThreadN      $THREADS \
  --genomeDir       "$GENOME_DIR" \
  --sjdbOverhang    $SJDB_OVERHANG \
  --readFilesIn     "$TRIM1" "$TRIM2" \
  --outFileNamePrefix "${SAMPLE}." \
  --outFilterIntronMotifs RemoveNoncanonicalUnannotated

#— 3) Convert SAM → BAM
echo "Converting SAM to BAM..."
sambamba view \
  -t $THREADS \
  -f bam \
  -S \
  -o "${SAMPLE}.bam" \
  "${SAMPLE}.Aligned.out.sam"

#— 4) Sort by coordinate
echo "Sorting BAM by coordinate..."
sambamba sort \
  -t $THREADS \
  -m 3GB \
  -p \
  -o "sorted_${SAMPLE}.bam" \
  "${SAMPLE}.bam"

#— 5) Remove duplicates
echo "Removing duplicates..."
java -jar "$EBROOTPICARD/picard.jar" MarkDuplicates \
  I="sorted_${SAMPLE}.bam" \
  O="dedup_${SAMPLE}.bam" \
  M="${SAMPLE}.rmdup_metrics.txt" \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=true

#— 6) Sort by read name (for HTSeq)
echo "Sorting BAM by name..."
sambamba sort \
  -t $THREADS \
  -m 3GB \
  -p \
  -n \
  -o "name_sorted_dedup_${SAMPLE}.bam" \
  "dedup_${SAMPLE}.bam"

#— 7) Count reads per gene
echo "Counting with HTSeq..."
htseq-count \
  -f bam \
  -r name \
  -s no \
  --nonunique all \
  "name_sorted_dedup_${SAMPLE}.bam" \
  "$GTF_FILE" \
  > "${SAMPLE}.counts.txt"

echo "Done! Outputs:"
echo "  Trimmed:    $TRIM1, $TRIM2"
echo "  Alignment:  ${SAMPLE}.Aligned.out.sam"
echo "  BAMs:       ${SAMPLE}.bam, sorted_${SAMPLE}.bam, dedup_${SAMPLE}.bam, name_sorted_dedup_${SAMPLE}.bam"
echo "  Counts:     ${SAMPLE}.counts.txt"
