#' Extract pseudotime data from Slingshot result
#' -----------------------------------------------
#' 
#' Parameters
#' ----------
#'   
#'
#' Returns
#' -------
#' 
get_sling_traject_data <- function(sling_set) {
  # get DR data and cluster labels
  d <- as_tibble(reducedDim(sling_set), rownames = "cell")
  d <- as_tibble(sling_set@colData$slingshot@elementMetadata$clusterLabels, rownames = "cell") %>% 
    right_join(d, by = "cell")

  # extract curve data
  curve_data <- sling_set@colData$slingshot@metadata$curves  
  d <- imap(curve_data, function(curve, curve_name) {
    # get DR coordinates of curve
    c <- as_tibble(curve$s, rownames = "cell") #here is the problem occurrs, *cell name is not retreived*
    # get pseudotime values of curve
    c$sling_ord <- curve$ord
    # set curve name
    c$curve <- curve_name
    c
  }) %>%
    bind_rows() %>% 
    right_join(d, by = "cell", suffix = c(".curve", ".cell"))

  return(d)
}
