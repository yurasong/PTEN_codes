###############################################################################

# load umap base theme gist
# devtools::source_gist("https://gist.github.com/milanmlft/e509140797bd0b88de2d7354c8bda579")


#' Make FeaturePlots with Viridis color scale
#' -----------------------------------------------
#' 
#' Parameters
#' ----------
#'   seuset :    Seurat object.
#'   markers :   markers to plot.
#'   reduction : DR embedding to use.
#'   viridis_scale : Viridis color scale to use (see viridis::scale_color_viridis)
#'                    Default: "inferno"
#'   ... :       arguments passed on to Seurat::FeaturePlot()
#'
#' Returns
#' -------
#' A list of ggplot2 objects for each marker.
#' If only one marker given, returns a ggplot2 object.

FeaturePlotViridis <- function(
  seuset,
  markers,
  reduction,
  viridis_scale = "inferno",
  ...) {

  # construct plot list with Seurat::FeaturePlot
  plots <- FeaturePlot(
    object = seuset,
    features = markers,
    reduction = reduction,
    combine = FALSE,
    ...
  ) %>%
    # add viridis color scale
    map(~ . + scale_color_viridis_c(option = "inferno") +
      # fix order of legends (useful when shape & color scales are used)
      guides(shape = guide_legend(order = 1), color = guide_colorbar(order = 2)))

  if (reduction == "dm") {
    # fix axes when using Diffusion Map
    plots <- plots %>%
      map(~. + scale_x_continuous() + scale_y_continuous())
  }
  
  # name plots
  names(plots) <- plots %>% map(~.$labels$title)
  
  if (length(markers) == 1) {
    plots <- plots[[1]]
  }
  
  return(plots)
}



#' Make pairs plot of reduced dimension components
#' -----------------------------------------------
#' 
#' Plots combinations of Principal or Diffusion Components on a scatter matrix.
#' Diagonal shows densities of components.
#' 
#' Parameters
#' ----------
#'   seuset :    Seurat object.
#'   reduction : DR embedding to use ("pca" or "dm").
#'   n_components : how many components to plot
#'   group.by :  name of one or more metadata columns to group (color) cells by
#'
#' Returns
#' -------
#' A ggmatrix object
plotDimRedPairs <- function(seuset, reduction, n_components, group.by = "ident") {
  meta <- FetchData(seuset, vars = group.by) %>% 
    rownames_to_column(var = "cell")
  
  Embeddings(seuset, reduction = reduction) %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>% 
    as_tibble() %>% 
    select(seq_len(1 + n_components)) %>% 
    left_join(meta, by = "cell") %>% 
    select(-cell) %>% 
    
    ggscatmat(columns = seq_len(n_components),
              color = group.by, alpha = 0.6) +
    theme(legend.title = element_text(face = "bold"))
}
