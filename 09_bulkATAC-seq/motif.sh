#!/bin/bash
# Submission script for Hercules2
#SBATCH --job-name=motif_test #please change it when you run your data
#SBATCH --time=05:00:00 # hh:mm:ss
#
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5120 # 5GB
#SBATCH --partition=batch
#
#SBATCH --mail-user=yura.song@ulb.be
#SBATCH --mail-type=ALL

echo --- Start ---

ml load releases/2020a
ml load SAMtools
ml load Perl

findMotifsGenome.pl sampleB_peaks.bed mm10 test_chr1_nobg -S 15 -len 6,8,10,12 -size -250,250 

#findMotifsGenome.pl sampleB_peaks.bed mm10 test_chr1_nobg -preparse -S 15 -len 6,8,10,12 -size -250,250 
# -preparse option should be added when you run HOMER first time in your cluster.


echo --- DONE ---
