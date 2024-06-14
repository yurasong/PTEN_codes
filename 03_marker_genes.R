# Library

library(Seurat)

library(RColorBrewer)
library(dplyr)
library(plotly)
library(grid)
library(data.table)
library(tidyverse)

library(viridis)

# Load data

seuset <- readRDS("WT_PTEN.rds")

# Take UMAP coordinate and gene expresion

umap_coord <- as.data.frame(seuset@reductions$umap@cell.embeddings)

seuset_gene <- merge.data.frame(umap_coord,t(as.data.frame(seuset@assays$RNA@data)),by=0)
row.names(seuset_gene) <- seuset_gene$Row.names

# Set-up the palette``

myPalette <- viridis(n = 10, option = "C", direction = -1)
sc <- scale_colour_gradientn(colours = rev(myPalette))

# Plotting function

umap_color_scaled <- function(x)
{
  regulon <- x
  regulon_d <- seuset_gene[,regulon]
  regulon_plot <- ggplot(data=seuset_gene) + geom_point(aes(x=UMAP_1,y=UMAP_2,colour=regulon_d), size=1)
  regulon_plot <- regulon_plot + labs(colour=regulon)
  regulon_plot <- regulon_plot + sc
  regulon_plot <- regulon_plot + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  regulon_plot
}

# If there are several list of genes

selected_regulons <- c("Mki67", "Top2a", "Cenpa", "Krt14", "Krt5", "Trp63",
                       "Jun", "Fos", "Fosb", "Krt13",
                       "Hoxb13", "Abo", "Spink1", "Nkx3-1", "Mmp7", "Fgl1", "Pbsn",
                       "Dpp4", "Krt4", "Psca", "Foxi1", "Atp6v1g3", "Gsdma", "Tgm4",
                       "Cldn10", "Msmb", "Cldn10", "Trpv6", 
                       "Clu", "Ppp1r1b", "Tacstd2", "Wfdc2", "Krt7",
                       "Cx3cr1", "Cd74", "H2-Aa", "Cxcl1", "Cxcl2", "Ifit1", "Irf6", "Ifi202b") # Genes we used for this paper

for (i in 1:length(selected_regulons)){
  umap_color_scaled(selected_regulons[i])
  ggsave(paste(selected_regulons[i],"_Exp.pdf",sep=""), width = 10, height = 10, units = "cm")
}