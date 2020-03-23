#' This function takes in the H coefficient matrix object from and NMF object and returns a matrix object with the topic profile for each cell type
#'
#' @param h Object of class LDA_Gibbs.
#' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @return This function returns a list where the first element is a matrix with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
#' @export
#' @examples
#'

topic_profile_per_cluster_nmf <- function(H,
                                          train_cell_clust,
                                          clust_vr) {
  
  # Check variables
  if (!is(h, "matrix")) stop("ERROR: h must be a matric object!")
  if (! is(train_cell_clust, "vector")) stop("ERROR: train_cell_clust must be a vector/list object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  
  # Loading libraries
  suppressMessages(require(tibble))
  suppressMessages(require(dplyr))

  h_ds <- data.frame(t(H))
  h_ds[, clust_vr] <- allen_reference_down@meta.data[, clust_vr]
  
  ct_topic_profiles <- h_ds %>%
    dplyr::group_by(!!! syms(clust_vr)) %>%
    dplyr::summarise_all(list(median)) %>%
    tibble::column_to_rownames(clust_vr) %>% 
    as.matrix()
  
  ct_topic_profiles_t <- t(ct_topic_profiles)
  colnames(ct_topic_profiles_t) <- gsub("[\\+|\\ |\\/]", ".", colnames(ct_topic_profiles_t))
  
  return(ct_topic_profiles_t)
}
