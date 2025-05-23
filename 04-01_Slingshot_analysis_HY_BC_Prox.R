#######################################################################################
# 04-01_Slingshot_analysis_HY_BC_Prox.R
# This script will perform slingshot analysis.
# This will return the result for Fig.4f.
#######################################################################################

# Library

library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tidyverse)
library(ggpubr)
library(ggthemes)
library(GGally)
library(RColorBrewer)

## Data

pten6w <- readRDS("PTEN_6W.rds")

pten6w <- subset(pten6w, 
                 idents=c("BC p63high", "BC p63low", "BC Hillock", "Proximal", "HY BC Prox")) 

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

PCAPlot(pten6w, group.by = "cell_type")
ggsave(filename = file.path(fig_path, "Slingshot_input-PCA.pdf"))

## define figures path for slingshot results
sling_fig_path <- file.path(fig_path, "slingshot-PCA")
dir.create(sling_fig_path)

## Data extraction and convert to SCE object

exprs_data <- GetAssayData(pten6w, slot = "data")
meta <- pten6w@meta.data[c("cell_type")]
pca <- Embeddings(pten6w, reduction.type = "pca")[, 1:30]

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
                         start.clus = c("BC p63high", "BC p63low")) 
SlingshotDataSet(sling_fixed)

## Add slingshot pseudotimelines to seurat object and plot ordering of cells

n_lineages <- SlingshotDataSet(sling_fixed)@lineages %>% length()

# remove old slingshot results
pten6w@meta.data[which(
  str_detect(names(pten6w@meta.data), 
             "^slingPseudotime"))] <- NULL

# add new results
sling_names <- names(sling_fixed@colData) %>%
  .[str_detect(., "^slingPseudotime")]

pten6w@meta.data[sling_names] <- as.data.frame(sling_fixed@colData[sling_names])

# Visualisation on input UMAP plot

sling_names <- names(seuset@meta.data) %>%
  .[str_detect(., "^slingPseudotime")]

sling_umaps <- FeaturePlot(seuset,
                           reduction = "umap",
                           features = sling_names, combine = FALSE
) %>%
  map(~ . + scale_color_viridis_c(option = "inferno", na.value = "light grey") +
        theme_void() + theme(aspect.ratio = 1))

CombinePlots(plots = sling_umaps, ncol = 2) %>% 
  annotate_figure(
    top = text_grob("Slingshot with fixed topology",
                    face = "bold", size = 16)
  )

# Visualisation on input PCA plot

sling_names <- names(seuset@meta.data) %>%
  .[str_detect(., "^slingPseudotime")]

sling_umaps <- FeaturePlot(seuset,
                           reduction = "pca",
                           features = sling_names, combine = FALSE
) %>%
  map(~ . + scale_color_viridis_c(option = "inferno", na.value = "light grey") +
        theme_void() + theme(aspect.ratio = 1))

CombinePlots(plots = sling_umaps, ncol = 2) %>% 
  annotate_figure(
    top = text_grob("Slingshot with fixed topology",
                    face = "bold", size = 16)
  )

# Pseudotemporal gene dynamics

library(gam, quietly = TRUE)

# set pseudotime variables
t1 <- pten6w$slingPseudotime_1
t2 <- pten6w$slingPseudotime_2

# extract gene expression data for var.genes
Y <- FetchData(pten6w, vars = VariableFeatures(pten6w), slot = "data")

gam_fdr_t1 <- fitPseudotimeGAM_nopal(t1, Y) 
gam_fdr_t2 <- fitPseudotimeGAM_nopal(t2, Y)

pten6w@misc$PT_DE_results$slingshot <- list(
  TI_method = "Slingshot",
  DE_method = "GAM",
  results = list(lineage1 = gam_fdr_t1,
                 lineage2 = gam_fdr_t2)
)

## Heatmap of top100 DE

library(pheatmap)

# make heatmaps of top 100 differential gene along each lineage
sling_gam_results <- seuset@misc$PT_DE_results$slingshot$results
lineage_pt <- list(
  "lineage1" = seuset$slingPseudotime_1,
  "lineage2" = seuset$slingPseudotime_2)

# make heatmap for each lineage
pt_heatmaps <- pmap(
  list(lineage_pt, sling_gam_results, names(lineage_pt)),
  ~ makePseudotimeHeatmap(
    pt = ..1,
    gam_result = ..2,
    response_data = FetchData(seuset,
                              vars = VariableFeatures(seuset),
                              slot = "data"
    ),
    n = 100,
    anno_data = data.frame(
      "updated_identity" = droplevels(seuset$updated_identity),
      "Pseudotime" = .x
    ),
    anno_colors = list(
      "updated_identity" = ct_cols,
      "Pseudotime" = colorRampPalette(brewer.pal(11, "RdYlBu"))(50)
    ),
    show_rownames = TRUE,
    scaled = TRUE,
    main_title = paste("Slingshot", ..3, ": top 100 significant DE genes")
  )
)
