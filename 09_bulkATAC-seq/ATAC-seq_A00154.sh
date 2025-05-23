#!/bin/bash
darkred='\e[0;31m'
white='\e[1;37m'
lightblue='\e[1;34m'
darkgreen='\e[0;32m'
pink='\e[1;35m'

##################################################################################################
#INSTRUCTIONS
echo -e "${darkred}INSTRUCTIONS OF THIS SCRIPT :"
echo -e "Input file format should be fastq format."
echo -e "The raw data should be paired-end data."
echo -e "The adapter sequences used for library preparation is required."
echo -e "The bowtie2 index of reference genome is required. The index file format is .ebwt."
echo -e "${darkgreen}To interrupt this script, type Ctrl+C"
######################################################################################################

###################################################################################################
# MANAGE GENOME VERSION ON WHICH THE ALIGMENT WILL BE PERFORMED
echo -e "${darkred}Enter the reference genome${pink}"
read genome
echo -e "${darkred}Is reference genome ${lightblue} $genome ${darkred}? If not, press Ctrl+C${white}"
if [[ "$genome" =~ ^mm ]]
then
size=`echo "mm"`
fi
if [[ "$genome" =~ ^Grcm38 ]]
then
size=`echo "mm"`
fi
if [[ "$genome" =~ ^hg ]]
then
size=`echo "hs"`
fi
if [[ "$genome" =~ ^GRCh ]]
then
size=`echo "hs"`
fi
###################################################################################################

###################################################################################################
# File with the list of all adapters used in the libraries preparation
echo -e "${darkred}Enter the location and name of adapter sequence file${pink}"
read adapter
####################################################################################################

####################################################################################################
#Files with the reads
fastqR1=()
fastqR2=()
fastqcutR1=()
fastqcutR2=()
echo -e "${darkred}Please inform the the absolute path and the name of raw fastq file. Since ATAC-seq is in paired-end, you need to check there are two fastq files called R1 and R2.${lightblue}"

	echo -e "${darkred}Indicate the information of R1 file ${pink}"
	read fastqR1
	echo -e "${darkred}Is your file name ${lightblue} $fastqR1 ${darkred}? If not, press Ctrl+C${white}"
	fastqR1=$fastqR1
	fastqcutR1=`awk -F'/' '{print $(NF)}' <<< $fastqR1`
	echo -e "${darkred}Indicate the information of R2 file ${pink}"
	read fastqR2
	echo -e "${darkred}Is your file name ${lightblue} $fastqR2 ${darkred}? If not, press Ctrl+C${white}"
	fastqR2=$fastqR2
	fastqcutR2=`awk -F'/' '{print $(NF)}' <<< $fastqR2`

########################################################################################################

#sample name retrieval
samplename=K5PTEN-6W-VP1

# Remove adapter sequences
echo -e "${darkgreen}Remove adapter sequences :${lightblue}"

trim1=`echo "$fastqR1.trim.fq"`
trim2=`echo "$fastqR2.trim.fq"`
untrim1=`echo "$fastqR1.untrim.fq"`
untrim2=`echo "$fastqR2.untrim.fq"`

TrimmomaticPE -threads 4 -phred33 $fastqR1 $fastqR2 $trim1 $untrim1 $trim2 $untrim2 ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#ALIGNMENT
echo -e "${darkgreen}Alignment to the reference genome with bowtie2:${pink}"
align=`echo "alignment_$samplename.sam"`

bowtie2 -X 2000 -q -p 4 -t --fr --very-sensitive --no-discordant --no-unal --no-mixed -x /home/audrey/reference_genomes/$genome -1 $trim1 -2 $trim2 -S $align

#--no-discordant: turn off the default setting, which looks for discordant alignments if bowtie cannot find any concordant alignments (unpaired reads suppression)
#--no-unal: suppress SAM record if the reads are failed to be aligned (non-mapped suppression)
#--no-mixed: turn off the default setting, when bowtie2 could not find a concordant or discordant alignment for a pair, then it tries to find alignment for individual mates
#--non-deterministic: Not necessarily report the same alignment for two individual reads, appropriate in situation where the input consists of many identical reads (really need it?)

#Mitochondreal reads and unmapped contigs filtration
echo -e "${darkgreen}Post-alignment quality control :${lightblue}"

align_noMT=`echo "alignment_noMT_$samplename.sam"`
alignbam_noMT=`echo "alignment_noMT_$samplename.bam"`

echo -e "${darkgreen}Mitochondreal reads and unmapped contigs filtration: ${lightblue}"

sed '/MT/d;/random/d;/Un/d;/gl0/d;/hap/d;/KI/d;/GL/d;/chrM/d;' < $align > $align_noMT #MT and chrM both should be put. SAM file generated from bowtie2 has MT as a chromosome information.
samtools view -bSh $align_noMT > $alignbam_noMT

#Read correction
#Adjusted the reads start sites to represent the center of the transposon binding event. As the Tn5 transposase show that the transposon binds as a dimer and inserts two adaptors separated by 9bp, reads aligning to the + strand are offset by +4bp and reads aligning to the - strand were offset -5bp (flag=16)

echo -e "${darkgreen}Adjusting the reads start sites to represent the center of the transposon binding event: ${lightblue}"

align_noMT_correct=`echo "alignment_noMT_corrected_$samplename.sam"`
alignbam_noMT_correct=`echo "alignment_noMT_corrected_$samplename.bam"`
alignbam_noMT_correct_sort=`echo "alignment_noMT_corrected_sort_$samplename.bam"`

awk '{if($2^VN || $2^SN || $2^ID){ print $_} else if($2=83 || $2==115 || $2==147 || $2==179){ $4=$4-5; print $_;} else { $4=$4+4; print $_;}}' $align_noMT > $align_noMT_correct
samtools view -bSh $align_noMT_correct > $alignbam_noMT_correct

rm $align_noMT
rm $align_noMT_correct #Deleting all the sam files since the further steps are always done with BAM files.

samtools sort -O BAM -T tempsort -o $alignbam_noMT_correct_sort $alignbam_noMT_correct
samtools index -b $alignbam_noMT_correct_sort $alignbam_noMT_correct_sort".bai" 

samtools flagstat $alignbam_noMT_correct_sort > $alignbam_noMT_correct_sort".stat"

##Reads with only one unique alignment
echo -e "${darkgreen}Reads with only one unique alignment: ${lightblue}"

alignbam_noMT_correct_sort_unique=`echo "alignment_noMT_corrected_sort_unique_$samplename.bam"`

samtools view -H $alignbam_noMT_correct_sort > header.sam 
samtools view -F 4 $alignbam_noMT_correct_sort | grep -v "XS:" | cat header.sam - | samtools view -b - > $alignbam_noMT_correct_sort_unique #putting headers back
#XS tag: Suboptimal alignment score

##Mapping quality at least 20
echo -e "${darkgreen}Mapping quality > 20: ${lightblue}"

alignbam_noMT_correct_sort_unique_q20=`echo "alignment_noMT_corrected_sort_unique_q20_$samplename.bam"`
samtools view -b -q 20 -o $alignbam_noMT_correct_sort_unique_q20 $alignbam_noMT_correct_sort_unique
samtools index -b $alignbam_noMT_correct_sort_unique_q20 $alignbam_noMT_correct_sort_unique_q20".bai"

##Retrieve only the reads who are properly mapped and not duplicate
echo -e "${darkgreen}Retrieving properly mapped reads: ${lightblue}"

alignbam_noMT_correct_sort_unique_filt=`echo "alignment_noMT_corrected_sort_unique_filtered_$samplename.bam"`
samtools view -h -b -F 1804 -f 2 $alignbam_noMT_correct_sort_unique_q20 > $alignbam_noMT_correct_sort_unique_filt

samtools flagstat $alignbam_noMT_correct_sort_unique_filt > $alignbam_noMT_correct_sort_unique_filt".stat"

# samtools view -F 1804 : reads or mates unmapped, non primary alignment, reads fails to quality check, PCR duplicates remove (from ENCODE)
# Remove reads unmapped, mate unmapped, not primary alignment, reads failing platform, duplicates (-F 1804)
# Retain properly paired reads -f 2

#Suppress duplicated reads by picard
echo -e "${darkgreen}Duplicated reads suppression by picard: ${lightblue}"

alignbam_noMT_correct_sort_unique_filt_nodup=`echo "alignment_noMT_corrected_sort_unique_filtered_nodup_$samplename.bam"`

picard.jar MarkDuplicates I=$alignbam_noMT_correct_sort_unique_filt O=$alignbam_noMT_correct_sort_unique_filt_nodup M=alignment_unique_nodup.metrix REMOVE_DUPLICATES=true
samtools index -b $alignbam_noMT_correct_sort_unique_filt_nodup $alignbam_noMT_correct_sort_unique_filt_nodup".bai"

samtools flagstat $alignbam_noMT_correct_sort_unique_filt_nodup > $alignment_noMT_corrected_sort_unique_filtered_nodup".stat" 

#Calculate matrix after quality control
after=`echo "Metrics_afterQuality"`
mkdir $after

picard.jar CollectMultipleMetrics I=$alignbam_noMT_correct_sort_unique_filt_nodup O=$after/"multiple_metrics" R="/home/audrey/reference_genomes/"$genome".fa" PROGRAM=null PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics

#ATAC-seq peak calling
echo -e "${darkgreen}Peak calling with MACS2: ${lightblue}"

peak=`echo "Peak_calling"`
mkdir $peak

macs2 callpeak -t $alignbam_noMT_correct_sort_unique_filt_nodup --outdir $peak -f BAMPE -g $size -n "$samplename" -q 0.01 --nomodel --shift 0 

#Blacklisted region removed
cd Peak_calling
sed '/^$/d' *.xls > peaks.bed
sed '/#/d' peaks.bed > peaks1.bed
sed '1d' peaks1.bed > peaks.bed
awk '{ print $1"\t"$2"\t"$3"\t"$10"\t"$7}' peaks.bed > peaks1.bed
mv peaks1.bed peaks.bed
bedtools intersect -v -a peaks.bed -b "/home/audrey/reference_genomes/blacklisted/"$genome-blacklist.bed > blacklist_removed_$samplename.bed
rm peaks.bed

cd ../

#count number of reads following the conditions
echo -e "${lightblue}The number of reads will be counted ${white}"

echo -e "${lightblue}The number of reads aligned: ${white}"
grep -c A00154 $align

echo -e "${lightblue}The number of reads after suppressing mitochondreal reads: ${white}"
samtools view -c $alignbam_noMT

echo -e "${lightblue}The number of reads uniquely mapped, no mitochondreal reads: ${white}"
samtools view -c $alignbam_noMT_correct_sort_unique

echo -e "${lightblue}The number of reads uniquely mapped, no mitochondreal reads and quality score > 20: ${white}"
samtools view -c $alignbam_noMT_correct_sort_unique_q20

echo -e "${lightblue}The number of reads uniquely and properly mapped, no mitochondreal reads and quality score > 20: ${white}"
samtools view -c $alignbam_noMT_correct_sort_unique_filt

echo -e "${lightblue}The number of reads uniquely and properly mapped, no mitochondreal reads, quality score > 20 and duplicated read removed: ${white}"
samtools view -c $alignbam_noMT_correct_sort_unique_filt_nodup

#DONE
echo -e "${darkgreen}All procedures have been done!"
