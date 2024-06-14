# Library

library(Seurat)
library(patchwork)
library(tidyverse)
library(data.table)
library(magrittr)

# Load data

seuset <- readRDS("PTEN_6W.rds")

# DEG calculation

## DEG up-regulated on HY BC Prox

up_hy_vs_BC <- FindMarkers(seuset, 
                           ident.1 = "HY_BC_Prox",
                           ident.2 = c("BC_p63high", "BC_p63low"),
                           min.pct = 0.25, logfc.threshold = 1, 
                           only.pos = T)

### Visualisation as volcano

up_hy_vs_BC$sig <- "no"
up_hy_vs_BC$sig[up_hy_vs_BC$p_val_adj < 0.01] <- "Significant"
up_hy_vs_BC$sig[up_hy_vs_BC$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_vs_BC, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()


up_hy_vs_prox <- FindMarkers(seuset, 
                           ident.1 = "HY_BC_Prox",
                           ident.2 = "Proximal",
                           min.pct = 0.25, logfc.threshold = 1, 
                           only.pos = T)

### Visualisation as volcano

up_hy_vs_prox$sig <- "no"
up_hy_vs_prox$sig[up_hy_vs_prox$p_val_adj < 0.01] <- "Significant"
up_hy_vs_prox$sig[up_hy_vs_prox$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_vs_prox, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()


## DEG up-regulated on HY Nkx3.1

up_hyn_vs_BC <- FindMarkers(seuset, 
                           ident.1 = "HY_Nkx3-1",
                           ident.2 = "BC_Nkx3-1",
                           min.pct = 0.25, logfc.threshold = 1, 
                           only.pos = T)

### Visualisation as volcano

up_hyn_vs_BC$sig <- "no"
up_hyn_vs_BC$sig[up_hyn_vs_BC$p_val_adj < 0.01] <- "Significant"
up_hyn_vs_BC$sig[up_hyn_vs_BC$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hy_vs_prox, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()


up_hyn_vs_ventral <- FindMarkers(seuset, 
                             ident.1 = "HY_Nkx3-1",
                             ident.2 = "Ventral",
                             min.pct = 0.25, logfc.threshold = 1, 
                             only.pos = T)

### Visualisation as volcano

up_hyn_vs_ventral$sig <- "no"
up_hyn_vs_ventral$sig[up_hyn_vs_ventral$p_val_adj < 0.01] <- "Significant"
up_hyn_vs_ventral$sig[up_hyn_vs_ventral$p_val_adj >= 0.01] <- "Not_significant"

ggplot(data=up_hyn_vs_ventral, aes(x=avg_log2FC, y=-log(p_val_adj, 10), col=sig)) + 
  geom_point() +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal()

