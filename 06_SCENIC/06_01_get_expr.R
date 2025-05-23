#######################################################################################
# 06_01_get_expr.R
# This script will extract gene expression matrix from proceeded Seurat object.
# Output will be used for running of pySCENIC.
#######################################################################################

exprMat <- t(as.matrix(seurat_object@assays$RNA@data))
write.table(data.frame("cell_id"=rownames(exprMat), exprMat), "exprMat.tsv", row.names=F, col.names=T, sep="\t", quote=F)
