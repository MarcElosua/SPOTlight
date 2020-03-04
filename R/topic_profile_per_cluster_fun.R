#' This function takes in a seurat object and an LDA model generated from topic models and returns a list with a topic profile per cluster.
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param train_cell_clust Object of class vector with cluster of the cells used to train the model.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param verbose Object of class Logical determining if progress should be reported or not (TRUE by default).
#' @return This function returns a list where the first element is a matrix with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
#' @export
#' @examples
#'

topic_profile_per_cluster <- function(lda_mod, train_cell_clust, clust_vr) {

  # Check variables
  if (!is(lda_mod, "LDA_Gibbs")) stop("ERROR: lda_mod must be an LDA_Gibbs object!")
  if (! is(train_cell_clust, "vector")) stop("ERROR: se_obj must be a vector/list object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")

  # Loading libraries
  suppressMessages(require(tibble))
  suppressMessages(require(dplyr))
  suppressMessages(require(data.table))
  suppressMessages(require(dtplyr))
  suppressMessages(require(BiocGenerics))

  # Extract topics from LDA
  g_mtrx <- lda_mod@gamma # n of cells X n of topics
  g_mtrx <- as_tibble(lda_mod@gamma)

  # extract metadata from Seurat object
  colnames(g_mtrx) <- paste("topic_", seq_len(ncol(g_mtrx)), sep = "")

  # Make sure train_cell_clust is a vector
  train_cell_clust <- unlist(train_cell_clust)

  if (length(train_cell_clust) != nrow(g_mtrx)) stop("se_obj@meta.data and lda_mod@gamma don't have the same number of rows (cells). Please make sure you're passing the same Seurat object you passes to the function 'train_LDA_fun()'")
  # generate clust_profiles
  clust_profiles <- cbind("cluster" = train_cell_clust, round(g_mtrx, 4)) %>%
    data.frame() %>%
    # tibble::as_tibble() %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise_all(list(median)) %>%
    tibble::column_to_rownames("cluster")

  colnames(clust_profiles) <- seq_len(ncol(clust_profiles))

  return(clust_profiles)
}
