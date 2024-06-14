###############################################################################
#
# PSEUDOTEMPORAL DYNAMICS: generating results and making figures
# --------------------------------------------------------------
###############################################################################


# Fit GAM to pseudotemporal data
# ------------------------------
#
# Fit a Generalized Additive Model (GAM) to pseudotemporal data. Smooths the
# pseudotime variable with LOESS. Returns FDR-adjusted p-values for each 
# response variable to allow discovery of significant changes along the
# pseudotime trajectory.
#
# Depends on the gam package to fit the model.
# Depends on the parallel package to allow parallelized computation.
#
# Parameters
# ----------
#   pt :               pseudotime variable. Wil be LOESS-smoothed.
#   response_data :   matrix of response variables to fit the curve to.
#                     Variables should be in COLUMNS.
# 
# Returns
# -------
# FDR-adjusted p-values for the fit of each input variable to the smoothed 
# pseudotime.

fitPseudotimeGAM <- function(pt, response_data) {
  # initialize cluster
  #no_cores <- detectCores() - 3
  #cl <- makeCluster(spec=no_cores, type = "PSOCK")
  #cl <- makeCluster(no_cores, type = "FORK")
  cl <- parallel::makeCluster(10)
  
  # fit GAM with a loess term for pseudotime to each response variable
  gam_pvals <- parApply(cl, response_data, 2, function(z){
    d <- data.frame(z = z, t = pt)
    tmp <- gam(z ~ lo(t), data = d)
    # extract p-values
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  
  stopCluster(cl)
  
  # Return FDR-adjusted p-values (Benjamini-Hochberg)
  p.adjust(gam_pvals, method = "BH")
}

# Use this function if you cannot use the parallel function

fitPseudotimeGAM_nopal <- function(pt, response_data) {
  
  # fit GAM with a loess term for pseudotime to each response variable - should not use parapply
  gam_pvals <- lapply(response_data, function(z){
    d <- data.frame(z = z, t = pt)
    tmp <- gam(z ~ lo(t), data = d)
    # extract p-values
    p <- summary(tmp)[4][[1]][1,5]
    p
  })
  
  # Return FDR-adjusted p-values (Benjamini-Hochberg)
  p.adjust(gam_pvals, method = "BH")
}

# Make heatmap for pseudotemporal data
# ------------------------------------
# 
# Description.
#
# Parameters
# ----------
#   pt :            Named vector of pseudotime values with the corresponding 
#                   cell IDs as names.
#   gam_result :    named vector of p-values for pseudotemporal change of 
#                   response data (see fitPseudotimeGAM)
#   response_data : dataframe of data to display. Variables should be in COLUMNS.
#                   Rownames should correspond to names of the pt vector.
#   n :             number of variables to display.
#   anno_data :     dataframe of meta-data to annotate the heatmap
#   anno_colors:    list for specifying annotation_row and annotation_col track colors manually. 
#                   It is possible to define the colors for only some of the features. 
#   show_rownames:  boolean specifying if row names are to be shown. Default is TRUE.
#   scaled :        logical. Should data be scaled? Default is TRUE.
#   reversed :      logical. Should gene order in heatmap be reversed?
#                   Defaults to FALSE.
#   main_title :    character string. Main title to display on the heatmap.
#                   Default is no title.
#   file_name :     location to save the heatmap. If not given, just prints the
#                   heatmap as output.
#
# Returns
# -------
# If filename given, saves the heatmap at the given file location.
# Otherwise prints out heatmap

makePseudotimeHeatmap <- function(
  pt, 
  gam_result, 
  response_data,
  n,
  anno_data,
  anno_colors,
  show_rownames = TRUE,
  scaled = TRUE,
  main_title = NULL,
  file_name = NA) {
  
  # drop NA values from GAM result
  gam_result <- gam_result[!is.na(gam_result)]
  
  # if n > number of variables, use all variables
  if (n > length(gam_result)) {
    n <- length(gam_result)
  }
  
  # select n most signficant variables
  top_names <- names(sort(gam_result, decreasing = FALSE))[1:n]
  
  # extract response data for top_names and sort according to pt
  heatdata <- response_data[names(pt[order(pt, na.last = NA)]), top_names]
  
  if (scaled) {
    # scale data per gene for visualization
    heatdata <- scale(heatdata)
    
    # trim Z-scores
    heatdata[heatdata > 1.7] <- 1.7
    heatdata[heatdata < -1.7] <- -1.7
  }
  
  # transpose heatdata to get genes as rows in the heatmap
  heatdata <-  t(heatdata)
  
  # # set heatmap annotation parameters
  # # order meta-data according to pseudotime and add it as a variable
  # ann_columns <- list(
  #   "Cell_type" = anno_data$cell_type[order(pt, na.last = NA)],
  #   "Pseudotime" = sort(pt, na.last = NA))
  # ann_colors <- list(
  #   "Cell_type" = brewer.pal(9, "Set1"),
  #   "Pseudotime" = colorRampPalette(brewer.pal(9, "Greys"))(50)
  #   )
  
  # create heatmap
  pheatmap(
    heatdata,
    scale = "none",
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    annotation_legend = TRUE,
    annotation_col = anno_data,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    show_rownames = show_rownames,
    main = main_title,
    fontsize = 14,
    fontsize_row = 10,
    labAnn = TRUE,
    annLegend = TRUE,
    filename = file_name,
    width = 16,
    height = 12
  )
}





# Plot regulon activity and expression of its target genes along PT trajectory
# ----------------------------------------------------------------------------
#
# Description.
#
# Uses ggplot2 for plotting, dplyr for data manipulation and RColorBrewer for
# colors.
#
# Parameters
# ----------
#   tf :      character string. Transcription factor for which the regulon
#             should be plotted
#   scenic :  data frame. SCENIC results containing activity AUC values for
#             each regulon.
#   seuset :  Seurat object containing the expression data and pseudotime
#             trajectories.
# 
# Returns
# -------
# ggplot2 object.

plotRegulonTargetPT <- function(tf, scenic, seuset) {
  # get targets belonging to this regulon
  targets <- seuset@misc$regulon_targets[[tf]]
  
  # get regulon activity from scenic_results
  pt_regulon_targets <- scenic %>% 
    # extract TF activity from scenic results
    select(Cell, tf) %>% 
    # rename column containing the TF AUC values
    # to avoid confusion later, since the TF itself can be a target gene of its
    # regulon!
    dplyr::rename(REGULON = tf) %>%
    # extract pseudotime values and expression values of the target genes from
    # the seuset object
    left_join(as_tibble(FetchData(
      seuset, 
      vars.all = c("Sling_km_Pseudotime1", "Sling_km_Pseudotime2", targets)), 
      rownames = "Cell"), by = "Cell") %>%
    dplyr::rename(
      Trajectory_1 = Sling_km_Pseudotime1, 
      Trajectory_2 = Sling_km_Pseudotime2) %>%
    gather(key = "name",
           value = "value",
           -Cell,
           -Trajectory_1,
           -Trajectory_2) %>% 
    # specify whether the "name" and "value" columns denote the regulon and its
    # activity or a target gene and its expression
    mutate(type = ifelse(name == "REGULON", "Regulon_AUC", "Target_Expression")) %>%
    gather(
      key = "trajectory",
      value = "pseudotime",
      -Cell,
      -name,
      -value,
      -type,
      na.rm = TRUE
    ) %>%
    # convert target and trajectory columns to factors
    mutate(
      name = factor(name),
      type = factor(type),
      trajectory = factor(trajectory)
    )
  
  # set colors depending on the number of target genes
  # add 1 for the regulon itself!
  n <- length(targets) + 1
  if (n > 20) {
    cols <- colorRampPalette(tableau_color_pal("Tableau 20")(20))(n)
  } else if (n > 10) {
    cols <- tableau_color_pal("Tableau 20")(n)
  } else {
    cols <- tableau_color_pal("Tableau 10")(n)
  }
  
  # set span parameter for geom_smooth
  ggplot(pt_regulon_targets, aes(pseudotime, value, col = name)) +
    geom_point(alpha = 0.3) + 
    geom_smooth(method = "loess", se = FALSE, span = 1.5) +
    scale_color_manual(name = NULL, values = cols) +
    ggtitle(paste("Regulon:", tf)) +
    facet_grid(type ~ trajectory, scales = "free")
}
