#######################################################################################
# 05_DEG analysis.R
# This script is used for DEG analysis and plotting of volcano plot.
# This code will generate plots in Extended Data Fig. 5g-j.
# Labelling of each data point is done manually.
#######################################################################################

# Library

library(Seurat)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

# Load data

pten6w <- readRDS("PTEN_6W.rds")

# Extended Data Fig.5g: HY BC Prox vs BC

up_hy_bc <- FindMarkers(pten6w, ident.1="HY BC Prox", ident.2="BC_p63high", 
                        min.pct=0.25, logfc.threshold=0.25, only.pos=T)

up_hy_bc$sig <- "no"
up_hy_bc$sig[up_hy_bc$p_val_adj < 0.01] <- "Significant"
up_hy_bc$sig[up_hy_bc$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_bc, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()

# Extended Data Fig.5h: HY BC Prox vs Proximal

up_hy_prox <- FindMarkers(pten6w, ident.1="HY BC Prox", ident.2="Proximal", 
                        min.pct=0.25, logfc.threshold=0.25, only.pos=T)

up_hy_prox$sig <- "no"
up_hy_prox$sig[up_hy_prox$p_val_adj < 0.01] <- "Significant"
up_hy_prox$sig[up_hy_prox$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_prox, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()

# Extended Data Fig.5i: HY Nkx3.1 vs BC Nkx3.1

up_hy_nkx <- FindMarkers(pten6w, ident.1="HY Nkx3.1", ident.2="BC Nkx3.1", 
                          min.pct=0.25, logfc.threshold=0.25, only.pos=T)

up_hy_nkx$sig <- "no"
up_hy_nkx$sig[up_hy_nkx$p_val_adj < 0.01] <- "Significant"
up_hy_nkx$sig[up_hy_nkx$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_nkx, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()

# Extended Data Fig.5j: HY Nkx3.1 vs Ventral

up_hy_ventral <- FindMarkers(pten6w, ident.1="HY Nkx3.1", ident.2="Ventral", 
                         min.pct=0.25, logfc.threshold=0.25, only.pos=T)

up_hy_ventral$sig <- "no"
up_hy_ventral$sig[up_hy_ventral$p_val_adj < 0.01] <- "Significant"
up_hy_ventral$sig[up_hy_ventral$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_ventral, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()
