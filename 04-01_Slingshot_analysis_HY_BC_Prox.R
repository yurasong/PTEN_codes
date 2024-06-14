## If destiny is not installed:
#install.packages("remotes")
#remotes::install_github("theislab/destiny")

# Library

library(Seurat)
library(SingleCellExperiment)
library(scater)
library(slingshot)
library(destiny) # Manual installation required if your bioconductor version > 3.10
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(GGally)
library(RColorBrewer)

## Data

seuset <- readRDS("PTEN_6W.rds")

seuset <- subset(seuset, 
                 idents=c("BC_p63high", "BC_p63low", "BC_Hillock", "Proximal", "HY_BC_Prox")) 

# Set figures and results paths

fig_path <- file.path("figures")
dir.create(fig_path, recursive = TRUE)
res_path <- file.path("results")
dir.create(res_path, recursive = TRUE)

# Load scripts

source("04_Slingshot_scripts/pseudotemporalDynamics.R")
source("04_Slingshot_scripts/plots.R")

# Slingshot on PCA

## Input data

PCAPlot(seuset, group.by = "cell_type")
ggsave(filename = file.path(fig_path, "Slingshot_input-PCA.pdf"))

## define figures path for slingshot results
sling_fig_path <- file.path(fig_path, "slingshot-PCA")
dir.create(sling_fig_path)

## Data extraction and convert to SCE object

exprs_data <- GetAssayData(seuset, slot = "data")
meta <- seuset@meta.data[c("cell_type")]
pca <- Embeddings(seuset, reduction.type = "pca")[, 1:30]

sce <- SingleCellExperiment(
  assays = exprs_data,
  colData = meta
)
logcounts(sce) <- exprs_data
reducedDim(sce, "PCA") <- pca

### Plotting the gene expression on SCE object

p1 <- plotPCA(sce, colour_by = "cell_type")
p2 <- plotPCA(sce, colour_by = "Krt14")
p3 <- plotPCA(sce, colour_by = "Tacstd2")
p4 <- plotPCA(sce, colour_by = "Krt4")
CombinePlots(list(p1, p2, p3, p4), ncol = 2)

# Apply Slingshot

## Fixed topology

sling_fixed <- slingshot(sce, clusterLabels = "cell_type", reducedDim = "PCA", 
                         start.clus = c("BC_p63high", "BC_p63low")) 
SlingshotDataSet(sling_fixed)

## Add slingshot pseudotimelines to seurat object and plot ordering of cells

n_lineages <- SlingshotDataSet(sling_fixed)@lineages %>% length()

# remove old slingshot results
seuset@meta.data[which(
  str_detect(names(seuset@meta.data), 
             "^slingPseudotime"))] <- NULL

# add new results
sling_names <- names(sling_fixed@colData) %>%
  .[str_detect(., "^slingPseudotime")]

seuset@meta.data[sling_names] <- as.data.frame(sling_fixed@colData[sling_names])

## Pseudotemporal gene dynamics

library(gam, quietly = TRUE)

# set pseudotime variables
t1 <- seuset$slingPseudotime_1
t2 <- seuset$slingPseudotime_2
t3 <- seuset$slingPseudotime_3
t4 <- seuset$slingPseudotime_4

# extract gene expression data for var.genes
Y <- FetchData(seuset, vars = VariableFeatures(seuset), slot = "data")

gam_fdr_t1 <- fitPseudotimeGAM_nopal(t1, Y) 
gam_fdr_t2 <- fitPseudotimeGAM_nopal(t2, Y)
gam_fdr_t3 <- fitPseudotimeGAM_nopal(t3, Y) 
gam_fdr_t4 <- fitPseudotimeGAM_nopal(t4, Y) 

seuset@misc$PT_DE_results$slingshot <- list(
  TI_method = "Slingshot",
  DE_method = "GAM",
  results = list(lineage1 = gam_fdr_t1,
                 lineage2 = gam_fdr_t2,
                 lineage3 = gam_fdr_t3,
                 lineage4 = gam_fdr_t4)
)

