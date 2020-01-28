#' This function takes in a seurat object and an LDA model generated from topic models and returns a list with a topic profile per cluster.
#'
#' @param lda_mod Object of class LDA_Gibbs.
#' @param se_obj Object of class Seurat.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param verbose Object of class Logical determining if progress should be reported or not (TRUE by default).
#' @return This function returns a list where the first element is a matric with the topic profiles of all possible combinations and the 2nd element is the cell composition of each spot.
#' @export
#' @examples
#'

topic_profile_per_cluster <- function(lda_mod, se_obj, clust_vr){

  # Check variables
  if(is(lda_mod)[[1]] != "LDA_Gibbs") {stop("ERROR: lda_mod must be a LDA_Gibbs object!")}
  if(is(se_obj) != "Seurat") {stop("ERROR: se_obj must be a Seurat object!")}
  if(!is.character(clust_vr)){stop("ERROR: clust_vr must be a character string!")}

  # cat("Loading packages...", sep="\n")
  suppressMessages(require(tibble))
  suppressMessages(require(dplyr))
  suppressMessages(require(data.table))
  suppressMessages(require(dtplyr))
  suppressMessages(require(BiocGenerics))

  # cluster assignment
  se_obj$seurat_clusters <-  droplevels(factor(factor(se_obj@meta.data[,clust_vr])))

  # Extract topics from LDA
  g_mtrx <- lda_mod@gamma # n of cells X n of topics

  # extract metadata from Seurat object
  colnames(g_mtrx) <- paste('topic_',1:ncol(g_mtrx),sep='')
  se_meta <- se_obj@meta.data

  # generate clust_profiles
  clust_profiles <- cbind(se_meta,g_mtrx) %>%
    # dtplyr::lazy_dt() %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::select(seurat_clusters, paste('topic_',1:nlevels(se_obj$seurat_clusters),sep='')) %>%
    dplyr::summarise_all(list(median)) %>%
    BiocGenerics::as.data.frame(clust_profiles) %>%
    tibble::column_to_rownames('seurat_clusters')

  colnames(clust_profiles) <- 1:ncol(clust_profiles)

  return(clust_profiles)
}
