#' This function downsamples the number of cells and genes used to train the model
#'
#' @param se_obj Object of class Seurat with the data of interest
#' @param clust_vr Name of the variable containing the cell clustering
#' @param cluster_markers_all Object of class dataframe obtained from the function Seurat::FindAllMarkers()
#' @param cl_n Object of integer indicating how many cells to keep from each cluster. If a cluster has n < cl_n then all cells will be selected, if it has more then cl_n will be sampled randomly, 100 by default.
#' @return A downsampled Seurat object from the original
#' @export
#' @examples

downsample_se_obj <- function(se_obj,
                              clust_vr,
                              cluster_markers_all,
                              cl_n = 10,
                              hvg = 1000) {

  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.data.frame(cluster_markers_all)) stop("ERROR: cluster_markers_all must be a data frame object returned from Seurat::FindAllMarkers()!")
  if (!is.numeric(cl_n)) stop("ERROR: cl_n must be an object of class integer!")

  # load required packages
  suppressMessages(require(Seurat))
  suppressMessages(require(purrr))
  suppressMessages(require(dplyr))
  suppressMessages(require(tibble))

  # se_obj$seurat_clusters <- droplevels(factor(se_obj@meta.data[, clust_vr]))

  if (hvg > 0) {
    se_obj <- Seurat::FindVariableFeatures(object = se_obj, nfeatures = hvg)
  }

  #### Combine marker genes and highest variable genes and subset genes ####
  keep_genes <- unique(c(VariableFeatures(se_obj), cluster_markers_all$gene))

  #### Get cell IDs to subset by cluster ####
  keep_ids <- lapply(split(se_obj@meta.data, se_obj@meta.data[, clust_vr]), function(subdf) {
    # Determine n_sample, if the size of the group is < cl_n use all the group, if not just use cl_n
    n_sample <- if_else(nrow(subdf) < cl_n, as.numeric(nrow(subdf)), as.numeric(cl_n))
    # Subset a random selection of that group and get the identifiers
    tmp_ds <- subdf[sample(seq_len(nrow(subdf)), n_sample), ] %>%
      tibble::rownames_to_column("ID") %>%
      dplyr::pull(ID)
    return(tmp_ds)
  }) %>%
    purrr::flatten_chr() # flatten the list into a vector

  #### Subset seurat object ####
  se_obj <- se_obj[keep_genes, keep_ids]

  return(se_obj)
}
