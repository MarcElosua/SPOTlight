#' This function takes in the H coefficient matrix object from and NMF object and returns a matrix object with the topic profile for each cell type
#'
#' @param h Object of class matrix, coefficient matrix from NMF model.
#' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
#' @return This function returns a list where the first element is a matrix with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
#' @export
#' @examples
#'

topic_profile_per_cluster_nmf <- function(h,
                                          train_cell_clust) {

  # Check variables
  if (!is(h, "matrix")) stop("ERROR: h must be a matrix object!")
  if (! is(train_cell_clust, "vector")) stop("ERROR: train_cell_clust must be a vector/list object!")

  # Loading libraries
  suppressMessages(require(tibble))
  suppressMessages(require(dplyr))

  h_ds <- data.frame(t(h))
  h_ds[, "clust_vr"] <- train_cell_clust

  ct_topic_profiles <- h_ds %>%
    dplyr::group_by(clust_vr) %>%
    dplyr::summarise_all(list(median)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "clust_vr") %>%
    as.matrix()

  ct_topic_profiles_t <- t(ct_topic_profiles)
  colnames(ct_topic_profiles_t) <- gsub("[[:punct:]]|[[:blank:]]",
                                        ".",
                                        colnames(ct_topic_profiles_t))

  return(ct_topic_profiles_t)
}
