# Library

library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(patchwork)
library(clustree)
library(tidyverse)
library(scales)
library(MASS)

# Load data

seuset <- readRDS("PTEN_6W.rds")

## Get information from meta.data
meta <- seuset@meta.data[,c(3, ncol(seuset@meta.data))] # Should return nFeature_RNA and cell_type information

## Load marker genes

bc <- readRDS("02_Marker_Genes/BC_markers.rds")
lc <- readRDS("02_Marker_Genes/LC_markers.rds")

# Get a count table and proceed the data

c_seuset <- as.matrix(seuset@assays[["RNA"]]@counts) # Count table for whole genes

## Fraction calculation on LC markers

lc_seuset <- c_seuset[rownames(c_seuset) %in% lc, ]
fraction_lc <- colSums(lc_seuset != 0)/nrow(lc_seuset)
lc_seuset_up <- rbind(lc_seuset, fraction_lc)

## Fraction calculation on BC markers

bc_seuset <- c_seuset[rownames(c_seuset) %in% bc, ]
fraction_bc <- colSums(bc_seuset != 0)/nrow(bc_seuset)
bc_seuset_up <- rbind(bc_seuset, fraction_bc)

## Take ratio of LC-BC and build the table

bc_frac <- as.matrix(bc_seuset_up[nrow(bc_seuset_up),])
colnames(bc_frac) <- "fraction_bc"

lc_frac <- as.matrix(lc_seuset_up[nrow(lc_seuset_up),])
colnames(lc_frac) <- "fraction_lc"

bc_lc <- merge(bc_frac, lc_frac, by="row.names")
rownames(bc_lc) <- bc_lc$Row.names
bc_lc <- bc_lc[, 2:3]

bc_lc_mer <- merge(bc_lc, meta, by="row.names")

rownames(bc_lc_mer) <- bc_lc_mer$Row.names
bc_lc_mer <- bc_lc_mer[, 2:ncol(bc_lc_mer)]
head(bc_lc_mer)

## Ratio correction by nFeature_RNA

bc_l <- rlm(fraction_bc ~ nFeature_RNA, data=bc_lc_mer)
lc_l <- rlm(fraction_lc ~ nFeature_RNA, data=bc_lc_mer)

bc_fitted <- mean(bc_lc_mer$fraction_bc) + bc_l$residuals[rownames(bc_lc_mer)]
lc_fitted <- mean(bc_lc_mer$fraction_lc) + lc_l$residuals[rownames(bc_lc_mer)]

bc_lc_mer[rownames(bc_lc_mer), "adjusted_bc_prop"] <- bc_fitted
bc_lc_mer[rownames(bc_lc_mer), "adjusted_lc_prop"] <- lc_fitted

# Plotting

## For all clusters

ggplot(data=bc_lc_mer, aes(x=adjusted_bc_prop, y=adjusted_lc_prop, color = cell_type)) +
  geom_point() +
  xlab("Adjusted BC proportion") +
  ylab("Adjusted LC proportion") +
  geom_abline(slope = 1) +
  xlim(0, 1.2) +
  ylim(0, 1.2) +
  scale_color_manual(values=c("#F8766D", "#00C1A7", "#DB8E00", "#FF63B6", "#C49A00", "#B385FF", "#00A6FF",
                              "#00C094", "#53B400", "#FF4500", "#EF67EB")) +
  theme(legend.position = "none")

## For hybrid clusters

hybrid <- bc_lc_mer[bc_lc_mer$cell_type == "HY BC Prox",]

ggplot(data=hybrid, aes(x=adjusted_bc_prop, y=adjusted_lc_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#00A6FF")) +
  xlab("Adjusted BC proportion") +
  ylab("Adjusted LC proportion") +
  geom_abline(slope = 1) +
  xlim(0.15, 1.2) +
  ylim(0.15, 1.2) +
  theme(legend.position = "none")

hybrid <- bc_lc_mer[bc_lc_mer$cell_type == "HY Nkx3.1",]

ggplot(data=hybrid, aes(x=adjusted_bc_prop, y=adjusted_lc_prop, color = cell_type)) +
  geom_point() +
  scale_color_manual(values=c("#FF4500")) +
  xlab("Adjusted BC proportion") +
  ylab("Adjusted LC proportion") +
  geom_abline(slope = 1) +
  xlim(0.15, 1.2) +
  ylim(0.15, 1.2) +
  theme(legend.position = "none")
