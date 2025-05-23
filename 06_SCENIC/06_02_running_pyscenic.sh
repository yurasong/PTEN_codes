#!/usr/bin/env bash

#######################################################################################
# 06_02_running_pyscenic.sh
# Please run this on bash environment.
# This script is used for running pySCENIC on seurat object.
# Database are retrieved from cisTarget database (https://resources.aertslab.org/cistarget/databases/).
# Database is provided by Stein Aerts lab.
# Conda environment is included in 06_00_conda_env.txt
# To adjust stochastic bias, we run the pipeline for 10 times with same parameters.
#######################################################################################

conda activate pyscenic

# Number of runs (adjust if needed)

NUM_RUNS=10

# Paths (adjust if needed)

EXPR="exprMat.tsv"
TF_LIST="mm_mgi_tfs.txt"
CISTARGET_DB="/home/pySCENIC/cisTarget_databases/mm10/mm10__refseq-r80*"
M2T="/home/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
 
for i in $(seq 1 $NUM_RUNS); do
  echo "=== Run #$i ==="

  # 1) Infer GRN
  OUT_ADJ="adj_${i}.csv"
  pyscenic grn \
    --num_workers 10 \
    -o $OUT_ADJ \
    $EXPR \
    $TF_LIST

  # 2) Prune with cisTarget
  OUT_REG="regulons_${i}.csv"
  pyscenic ctx \
    $OUT_ADJ \
    $CISTARGET_DB \
    --annotations_fname $M2T \
    --expression_mtx_fname $EXPR \
    --mode dask_multiprocessing \
    --output $OUT_REG \
    --num_workers 10

  # 3) Score regulons (AUC)
  OUT_AUC="auc_mtx_${i}.csv"
  pyscenic aucell \
    $EXPR \
    $OUT_REG \
    -o $OUT_AUC \
    --num_workers 10

  echo "Outputs: $OUT_ADJ, $OUT_REG, $OUT_AUC"
done