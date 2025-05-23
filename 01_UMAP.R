#######################################################################################
# 01_marker_genes.R
# This script will generate the UMAP of proceeded data, with custom colours.
#######################################################################################

# Library

library(Seurat)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)

# Proceeded Data

wt1 <- readRDS("WT_PTEN.rds")
pten6w <- readRDS("PTEN_6W.rds")
pten10m <- readRDS("PTEN_10M.rds")

# UMAP for control sample (Fig. 4a)

DimPlot(wt1, reduction = "umap", label = TRUE, 
        label.size = 4, repel = T,
        cols=c("BC p63 low" = "#F8766D",
               "BC Nkx3.1" = "#00C1A7",
               "AD" = "#53B400", 
               "Ventral" = "#00C094",
               "Lateral" = "#00B6EB",
               "Proximal" = "#FF63B6",
               "BC Pr" = "#B385FF",
               "Intersparced LC" = "#696969")) + 
  NoLegend()

# UMAP for PTEN 6W (Fig. 4b)

DimPlot(pten6w, reduction = "umap", label = TRUE, 
        label.size = 4, repel = T,
        cols=c("BC p63 high" = "#C49A00",
               "BC Nkx3.1" = "#00C1A7",
               "BC p63 low" = "#F8766D",
               "HY Nkx3.1" = "#FF4500", 
               "Ventral" = "#00C094",
               "AD" = "#53B400", 
               "BC Hillock" = "#DB8E00",
               "Proximal" = "#FF63B6",
               "HY BC Prox" = "#00A6FF",
               "LC meta" = "#EF67EB",
               "Prolif" = "#B385FF")) + 
  NoLegend()

# UMAP for PTEN 10M (Fig. 6a)

DimPlot(pten10m, reduction = "umap", label = TRUE, 
        label.size = 4, repel = T,
        cols=c("BC" = "#F8766D",
               "BC Interferon" = "#00C094",
               "HY BC Prox" = "#00A6FF",
               "HY Interferon" = "#C71585",
               "Prox-LC Interferon" = "#53B400",
               "Prox-LC ClassII Antigen" = "#FF63B6",
               "Prox-LC Chemokine" = "#F4A460",
               "Prox-LC Cx3cr1" = "#708090",
               "Pr" = "#B385FF")) + 
  NoLegend() +
  xlim(-5, 5) +
  ylim(-10, 8)