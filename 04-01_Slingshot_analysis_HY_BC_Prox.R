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

p1 <- seuset@meta.data %>%
  dplyr::mutate(cell_type = fct_reorder(cell_type, slingPseudotime_1)) %>%
  ggplot(aes(slingPseudotime_1, cell_type, col = cell_type)) +
  geom_jitter() + ggtitle("Slingshot cell ordering", subtitle = "Lineage 1")

ggsave(filename = file.path(sling_fig_path, "slingshot_ordering.pdf"), plot = p1, width = 16, height = 12)
p1

# Visualise on UMAP

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
ggsave(file.path(sling_fig_path, "sling_lineages_on_UMAP.pdf"), height=15, width=15, units="cm")

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

### Heatmaps of top 100 DE genes along pseudotime lineages

library(pheatmap)

# make heatmaps of top 100 differential gene along each lineage
sling_gam_results <- seuset@misc$PT_DE_results$slingshot$results
lineage_pt <- list(
  "lineage1" = seuset$slingPseudotime_1,
  "lineage2" = seuset$slingPseudotime_2,
  "lineage3" = seuset$slingPseudotime_3,
  "lineage4" = seuset$slingPseudotime_4)


# set celltype colors
#ct_colors <- brewer.pal(n = length(levels(droplevels(seuset$cell_type))), "Paired")
#names(ct_colors) <- levels(droplevels(seuset$cell_type))

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
      "cell_type" = droplevels(seuset$cell_type),
      "Pseudotime" = .x
    ),
    anno_colors = list(
      "cell_type" = ct_cols,
      "Pseudotime" = colorRampPalette(brewer.pal(11, "RdYlBu"))(50)
    ),
    show_rownames = TRUE,
    scaled = TRUE,
    main_title = paste("Slingshot", ..3, ": top 100 significant DE genes")
  )
)

# export plots
save_pheatmap_pdf <- function(x, filename, width=16, height=12) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

iwalk(
  pt_heatmaps,
  ~ save_pheatmap_pdf(
    .x,
    filename = file.path(sling_fig_path,
                         paste0(.y, "-top100-DE_heatmap.pdf"))
  )
)

