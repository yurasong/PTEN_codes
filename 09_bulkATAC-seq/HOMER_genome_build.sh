#!/bin/bash
# Submission script for Hercules2
#SBATCH --job-name=motif_genome
#SBATCH --time=01:30:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5120 # 5GB
#SBATCH --partition=batch
#
#SBATCH --mail-user=yura.song@ulb.ac.be
#SBATCH --mail-type=ALL

echo --- Start ---

ml load releases/2020a
ml load SAMtools
ml load Perl

perl /CECI/home/ulb/iribhm/ysong/HOMER/configureHomer.pl -install mm10
perl /CECI/home/ulb/iribhm/ysong/HOMER/configureHomer.pl -install hg38

module purge

echo --- DONE ---
