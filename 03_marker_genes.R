#######################################################################################
# 03_marker_genes.R
# This script will generate the feature plots which shows gene expression.
# Each related figures are mentioned as title.
#######################################################################################

# Library

library(Seurat)

library(RColorBrewer)
library(dplyr)
library(plotly)
library(grid)
library(data.table)
library(tidyverse)

library(viridis)


# Set-up the colour scale

myPalette <- viridis(n = 10, option = "C", direction = -1)
sc <- scale_colour_gradientn(colours = rev(myPalette))

# Definition of plotting function

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- seuset_gene[,regulon]
  regulon_plot <- ggplot(data=seuset_gene) + geom_point(aes(x=UMAP_1,y=UMAP_2,colour=regulon_d), size=1)
  regulon_plot <- regulon_plot + labs(colour=regulon)
  regulon_plot <- regulon_plot + sc
  regulon_plot <- regulon_plot + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
    
  regulon_plot
}

# Load proceeded data

wt1 <- readRDS("WT_PTEN.rds")
pten6w <- readRDS("PTEN_6W.rds")
pten10m <- readRDS("PTEN_10M.rds")

# Gene expression on WT

umap_coord <- as.data.frame(wt1@reductions$umap@cell.embeddings)

seuset_gene <- merge.data.frame(umap_coord,t(as.data.frame(wt1@assays$RNA@data)),by=0)
row.names(seuset_gene) <- seuset_gene$Row.names

marker_genes <- c("Krt13", "Aqp3", "Krt4", "Nkx3-1", # Genes in Fig.4
                  "Krt14", "Krt5", "Trp63", "Psca", "Clu", "Wfdc2", "Ppp1r1b",
                  "Sbp", "Sbpl", "Spink1", "Msmb", "Cldn10", "C1rb", "Tgm4", "Gsmda",
                  "Mki67", "Cenpa") # Genes in Extended Fig. 3b-g

lapply(marker_genes, umap_color_scaled) # For printing out each files on R

for (i in 1:length(marker_genes)){
  umap_color_scaled(marker_genes[i])
  ggsave(paste(marker_genes[i],"_Exp_WT.pdf",sep=""), width = 10, height = 10, units = "cm") # To save
}

# Gene expression on PTEN 6W

umap_coord <- as.data.frame(pten6w@reductions$umap@cell.embeddings)

seuset_gene <- merge.data.frame(umap_coord,t(as.data.frame(pten6w@assays$RNA@data)),by=0)
row.names(seuset_gene) <- seuset_gene$Row.names

marker_genes <- c("Krt13", "Aqp3", "Krt4", "Nkx3-1", # Genes in Fig.4
                  "Irf6", "Irf7", "Irf9", "Stat1", "Stat2", "Nfkb2", "Relb", # Genes in Fig.5
                  "Krt14", "Krt5", "Trp63", "Psca", "Clu", "Wfdc2", "Ppp1r1b", "Ltf", "Pigr",
                  "Krt6a", "Ly6d", "Sbp", "Sbpl", "Spink1", 
                  "Tgm4", "Gsmda", "Mki67", "Cenpa", # Genes in Extended Data Fig. 3i-o
                  "Elf3", "Grhl3", "Creb5") # Genes in Extended Data Fig.5k-m

lapply(marker_genes, umap_color_scaled) # For printing out each files on R

for (i in 1:length(marker_genes)){
  umap_color_scaled(marker_genes[i])
  ggsave(paste(marker_genes[i],"_Exp_PTEN_6W.pdf",sep=""), width = 10, height = 10, units = "cm") # To save
}

# Gene expression on PTEN 10M

umap_coord <- as.data.frame(pten10m@reductions$umap@cell.embeddings)

seuset_gene <- merge.data.frame(umap_coord,t(as.data.frame(pten10m@assays$RNA@data)),by=0)
row.names(seuset_gene) <- seuset_gene$Row.names

marker_genes <- c("Irf6", "Irf7", "Stat1", "Stat2", "Nfkb2", "Relb", # Genes in Fig.6b,
                  "Krt14", "Krt13", "Aqp3", "Krt4", "Clu", "Wfdc2", "Pigr", "Ppp1r1b", # Genes in Fig. 6c-e,
                  "Tacstd2", "Ltf", "Top2a", "Mki67", "Cenpa", "Ly6d", "Krt6a",
                  "Nkx3-1", "Tgm4", "Cxcl1", "Cxcl2", "Cxcl5", "Cx3cr1",
                  "Ifit1", "Ifitm2", "Ifi202b", "Cd74", "H2-Aa", "H2-Ab1") # Genes in Extended Data Fig.8b-h

lapply(marker_genes, umap_color_scaled) # For printing out each files on R

for (i in 1:length(marker_genes)){
  umap_color_scaled(marker_genes[i])
  ggsave(paste(marker_genes[i],"_Exp_PTEN_10M.pdf",sep=""), width = 10, height = 10, units = "cm") # To save
}
