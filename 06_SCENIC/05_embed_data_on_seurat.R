# Preparations

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(Matrix)
library(dplyr)
library(patchwork)
library(clustree)
library(tidyverse)
library(pheatmap)

# Put the datasets

## Average result

scenic_results <- read.table("AUC_Average_Collagen.csv", header=T, sep=",", row.names=1)
colnames(scenic_results) <- sub("...$","",colnames(scenic_results))

## Seurat object

seuset <- readRDS("seurat_object.rds")

### set labels as cell_type metadata (order levels by frequency)
seuset <- RenameIdents(seuset, `0` = "LC ER+ 1", `1` = "BC", `2` = "LC ER- 1") # This is still an example: need to adjust based on your data

### set manual labels as cell_type metadata (order levels by frequency)
seuset$cell_type <- fct_infreq(seuset@active.ident)

DimPlot(seuset, reduction = "umap", group.by = "cell_type", label = TRUE, label.size = 3) + theme_void() + theme(aspect.ratio = 1)

# Plotting

## UMAP - on R

d <- merge.data.frame(d,as.data.frame(scenic_results),by=0)
row.names(d) <- d$Row.names

# Embedding on Seurat

# AUC matrix as the assay of seurat object
## To put the AUC matrix as additional assay on seurat object

seuset_auc <- seuset 

auc_mtx <- t(scenic_results)
aucs <- CreateAssayObject(counts = auc_mtx)
seuset_auc[["AUC"]] <- aucs

DefaultAssay(seuset_auc) <- "AUC"

# rm(seuset)