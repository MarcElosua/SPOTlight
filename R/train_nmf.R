#' This functions seeds initialization matrices H and W to perform NMF
#'
#' @param cluster_markers 
#' @param se_obj 
#' @param ntop
#' @param k
#' @param method 
#' @param seed_init 
#' @return This function returns a list with the initialized matrices H and W.
#' @export
#' @examples
#'

train_nmf <- function(cluster_markers, se_obj, ntop, transf, ntop, k, method, seed_init) {
  
  init_mtrx <- seed_init_mtrx_nmf(cluster_markers = cluster_markers_filt, 
                                  se_obj = allen_reference_down, 
                                  ntop = 100)
  
  if (transf == "cpm") {
    count_mtrx <- edgeR::cpm(counts, normalized.lib.sizes = FALSE)
  }
  
  if (transf == "uv") {
    count_mtrx <- counts_uv <- scale(t(counts), 
                                   center = FALSE, 
                                   scale = apply(t(counts), 2, sd, na.rm = TRUE))
  }
    
  if (transf == "sct") {
    # Can't use scale.data since it has negative values
    count_mtrx <- as.matrix(allen_reference_down@assays$RNA@counts)
  }
  
  if (method == "nsNMF") mod <- "NMFns" else mod <- "NMFstd"
  
  # Initialize the matrix with the seeded matrices
  nmfns_init <- nmfModel(W = init_mtrx[["W"]], 
                         H = init_mtrx[["H"]],
                         model = mod)
  
  start_t <- Sys.time()
  nmf_mod <- nmf(x = counts_cpm, 
                 rank = k, 
                 seed = nmf_init)
  total_t <- round(difftime(Sys.time(), start_t, units = "mins"), 2)
  print(sprintf("Time to train NMF model was %smins", total_t))
  
  return(nmf_mod)
}