# Library

library(Seurat)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(scales)

# Load the data

wt_1 <- readRDS("WT_PTEN.rds")
pten_6w <- readRDS("PTEN_6W.rds")
pten_10m <- readRDS("PTEN_10M_maincluster_only.rds")

# UMAP for control sample

DimPlot(wt_1, reduction = "umap", label = TRUE, 
        label.size = 4, repel = T,
        cols=c("BC_p63low" = "#F8766D",
               "BC_p63high" = "#00C1A7",
               "AD" = "#53B400", 
               "Ventral" = "#00C094",
               "Lateral" = "#00B6EB",
               "Proximal" = "#FF63B6",
               "BC_Pr" = "#B385FF",
               "Intersparced_LC" = "#696969")) + 
  NoLegend()

# UMAP for PTEN 6W

DimPlot(seuset, reduction = "umap", label = TRUE, 
        label.size = 4, repel = T,
        cols=c("BC_p63high" = "#C49A00",
               "BC_Nkx3-1" = "#00C1A7",# need to change
               "BC_p63low" = "#F8766D",
               "HY_Nkx3-1" = "#FF4500", # need to change
               "Ventral" = "#00C094",
               "AD" = "#53B400", 
               "BC_Hillock" = "#DB8E00",
               "Proximal" = "#FF63B6",
               "HY_BC_Prox" = "#00A6FF",
               "LC_meta" = "#EF67EB",
               "Prolif" = "#B385FF")) + 
  NoLegend()

# UMAP for PTEN 10M

DimPlot(seuset, reduction = "umap", label = TRUE, 
        label.size = 4, repel = T,
        cols=c("BCs" = "#F8766D",
               "BC_Interferon" = "#00C094",
               "HY_BC_Prox" = "#00A6FF",
               "HY_Interferon" = "#C71585",
               "LC_Interferon" = "#53B400",
               "LC_ClassII_Antigen" = "#FF63B6",
               "LC_Chemokine" = "#F4A460",
               "LC_Cx3cr1/Tnfrsf9" = "#708090",
               "Prolif" = "#B385FF")) + 
  NoLegend() +
  xlim(-5, 5) +
  ylim(-10, 8)