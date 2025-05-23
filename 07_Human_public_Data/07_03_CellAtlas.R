#######################################################################################
# 07_03_CellAtlas.R
# This script is used for geneset enrichment analysis on public human prostate tumour scRNA-seq data.
# Reference article: https://www.sciencedirect.com/science/article/pii/S2211124721016284?via%3Dihub
# Public data source is mentioned in Data Availability section of manuscript.
#######################################################################################

# Library

library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)
library(SingleCellExperiment)
library(viridis)

# Load data

seuset <- readRDS("prostate_portal_300921_seuratobject.RDS")

seuset <- seuset[, seuset$group == "tumor"] # Take tumour cells only
seuset <- subset(seuset, idents=c("BE", "HE", "CE", "LE-KLK4", "LE-KLK3")) # Take only epithelial cells

# Load marker genes
prox <- readRDS("07_00_reprogramming_genes/Proximal.rds")
prox <- toupper(prox)

hy <- readRDS("07_00_reprogramming_genes/HY_BC_Prox.rds")
hy <- toupper(hy)

hillock <- readRDS("07_00_reprogramming_genes/BC_Hillock.rds")
hillock <- toupper(hillock)

wholeset <- c(prox, hy, hillock)

# Definition of signature

signature <- list(prox, hy, hillock, wholeset)
names(signature) <- c("Proximal", "HY_BC_Prox", "Hillock", "All")

# Perform enrichment analysis

seuset_mod <- AddModuleScore(
  object = seuset,
  features = signature,
  ctrl = 5,
  name = names(signature)
)

results <- seuset_mod@meta.data[, 12:15] # Take the enrichment score

# Plotting with colour scale

umap_coord <- seuset_mod@reductions$umap@cell.embeddings
head(umap_coord)

result_merged <- merge(umap_coord, results, by=0)
rownames(result_merged) <- result_merged$Row.names
result_merged <- result_merged[, -1]
head(result_merged)

myPalette <- viridis(n = 10, option = "C", direction = -1) # Palette

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- result_merged[,regulon]
  regulon_plot <- ggplot(data=result_merged) + 
    geom_point(aes(x=UMAP_1,y=UMAP_2,colour=regulon_d), size=1)
  regulon_plot <- regulon_plot + labs(colour=regulon)
  regulon_plot <- regulon_plot + sc
  regulon_plot <- regulon_plot + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
          #legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  regulon_plot
  
}

umap_color_scaled("Proximal1") + theme(legend.position = "none")
umap_color_scaled("HY_BC_Prox2") + theme(legend.position = "none")
umap_color_scaled("Hillock3") + theme(legend.position = "none")
umap_color_scaled("All4") + theme(legend.position = "none")

