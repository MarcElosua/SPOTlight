#' This functions carries out the NMF and returns and NMF object.
#'
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers()
#' @param se_obj Object of class Seurat with the data of interest
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), NULL (no transformation applied). By default CPM.
#' @param ntop Object of class "numeric"; number of unique markers per cluster used to seed the model, by default 100. If NULL it uses all of them.
#' @param clust_vr Object of class character; Name of the variable containing the cell clustering
#' @param method Object of class character; Type of method to us to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @return This function returns a list with the initialized matrices H and W.
#' @export
#' @examples
#'

train_nmf <- function(cluster_markers,
                      se_obj,
                      clust_vr,
                      ntop = 100,
                      transf = "cpm",
                      method = "nsNMF") {

  if (transf == "cpm") {
    counts <- as.matrix(se_obj@assays$RNA@counts)
    count_mtrx <- edgeR::cpm(counts,
                             normalized.lib.sizes = FALSE)

  } else if (transf == "uv") {
    counts <- as.matrix(se_obj@assays$RNA@counts)
    count_mtrx <- scale(t(counts),
                        center = FALSE,
                        scale = apply(t(counts), 2, sd, na.rm = TRUE))
    count_mtrx <- t(count_mtrx)

  } else if (transf == "sct") {
    # Can't use scale.data since it has negative values
    count_mtrx <- as.matrix(se_obj@assays$SCT@counts)

  } else if (is.null(transf)) {
    count_mtrx <- as.matrix(se_obj@assays$RNA@counts)

  }

  # Rank of the model equals the number of cell types
  k <- length(unique(se_obj@meta.data[, clust_vr]))

  #################################
  ###### Initialize the model #####
  #################################
  # Define initial seeding model and set the right type
  if (method == "nsNMF") mod <- "NMFns" else mod <- "NMFstd"

  # Get init seeding matrices
  init_mtrx <- seed_init_mtrx_nmf(cluster_markers = cluster_markers_filt,
                                  se_obj = se_obj,
                                  ntop = 100)

  # Initialize the matrix with the seeded matrices
  nmf_init <- nmfModel(W = init_mtrx[["W"]],
                       H = init_mtrx[["H"]],
                       model = mod)

  # Train NMF model
  start_t <- Sys.time()

  nmf_mod <- nmf(x = count_mtrx,
                 rank = k,
                 seed = nmf_init,
                 method = method)

  total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 2)
  print(sprintf("Time to train NMF model was %smins", total_t))

  return(list(nmf_mod, se_obj@meta.data[, clust_vr]))
}

