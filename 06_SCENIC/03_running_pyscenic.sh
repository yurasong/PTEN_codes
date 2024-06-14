## Directly work on the expression data extracted from Seurat object. Not going through the loom file creating and building annotation files. 
## Database are retrieved from SCENIC database, provided by Stein Aerts lab.

conda activate pyscenic

pyscenic grn --num_workers 10 -o adj.csv exprMat.tsv resources/mm_mgi_tfs.txt 

pyscenic ctx adj.csv /home/pySCENIC/cisTarget_databases/mm10/mm10__refseq-r80* --annotations_fname /home/pySCENIC/cisTarget_databases/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname exprMat.tsv --mode "dask_multiprocessing" --output regulons.csv --num_workers 10

pyscenic aucell exprMat.tsv regulons.csv -o auc_mtx.csv --num_workers 10

 
