#' This function takes in a seurat object with several tuning parameters and it assesses its performance on synthetic test spots
#'
#' @param se_obj Object of class Seurat.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param verbose Object of class Logical determining if progress should be reported or not (TRUE by default).
#' @param iter Object of class "integer"; number of Gibbs iterations, by default equals 5000.
#' @param nstart Object of class "integer". Number of repeated random starts.
#' @param keep Object of class "integer"; if a positive integer, the log-likelihood is saved every keep iterations, by default 100.
#' @param top_dist Object of class integer, on how many top euclidean distance are we going to calculate the JSD.
#' @param top_jsd Object of class integer, how many of the top spots according JSD distance are we going to use to determine the composition.
#' @param clust_vr Name of the variable containing the cell clustering
#' @return This function returns a list where the first element is the lda model trained, the second is a list with test spot counts + metadata and the third element are the raw_statistics.
#' @export
#' @examples
#'

spatial_decon_syn_assessment_fun <- function(se_obj, clust_vr, verbose = TRUE, iter = 7500, nstart = 1, keep = 100, top_dist = 1000, top_jsd = 15, cl_n = 100) {

  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.logical(verbose))stop("ERROR: verbose must be a logical object!")
  if (!is.numeric(iter)) stop("ERROR: iter must be of class integer!")
  if (!is.numeric(nstart)) stop("ERROR: nstart must be of class integer!")
  if (!is.numeric(nstart)) stop("ERROR: nstart must be of class integer!")
  if (!is.numeric(keep)) stop("ERROR: keep must be of class integer!")
  if (!is.numeric(top_dist)) stop("ERROR: top_dist must be an integer!")
  if (!is.numeric(top_jsd)) stop("ERROR: top_jsd must be an integer!")


  # Set Identities as the cluster variables
  Seurat::Idents(object = se_obj) <- se_obj@meta.data[, clust_vr]

  # Get marker genes for all the clusters
  cluster_markers_all <- Seurat::FindAllMarkers(object = se_obj, verbose = verbose, only.pos = T)

  # Filter marker genes by p value and logFC
  cluster_markers_all <- cluster_markers_all[cluster_markers_all$p_val_adj < 0.01 &
                                               cluster_markers_all$avg_logFC > 1, ]

  # Downsample seurat object to reduce n cells and n genes
  se_obj <- downsample_se_obj(se_obj = se_obj, clust_vr = clust_vr, cluster_markers_all = cluster_markers_all, cl_n = cl_n)

  #### Train LDA model ####
  set.seed(1000)
  start_time <- Sys.time()

  lda_mod_ls <- train_lda(se_obj = se_obj, clust_vr = clust_vr,
                          cluster_markers_all = cluster_markers_all, al = 0.01,
                          verbose = keep, iter = iter, burnin = 0,
                          best = TRUE, keep = keep, nstart = nstart)

  print(sprintf("Time to run LDA model for %s is %s",
                tech, round(difftime(Sys.time(), start_time, units = "mins"), 2)))

  # Select the best model
  lda_mod <- lda_mod_ls[[1]]

  # Generate test spots synthetically
  test_spots_ls <- test_spot_fun(se_obj = se_obj, clust_vr = clust_vr, n = 1000, verbose = verbose)

  test_spots_counts <- test_spots_ls[[1]]

  # Transpose spot_counts so its SPOTxGENES
  test_spots_counts <- BiocGenerics::t(test_spots_counts)
  test_spots_metadata <- test_spots_ls[[2]]
  test_spots_metadata <- as.matrix(test_spots_metadata[, which(colnames(test_spots_metadata) != "name")])

  # Perform spot deconvolution
  decon_mtrx <- spot_deconvolution(lda_mod = lda_mod, se_obj = se_obj,
                                   clust_vr = clust_vr,  spot_counts = test_spots_counts,
                                   verbose = verbose, ncores = 5, parallelize = TRUE,
                                   top_dist = top_dist, top_jsd = top_jsd)

  # Assess deconvolution performance
  raw_statistics_ls <- test_synthetic_performance(test_spots_metadata_mtrx = test_spots_metadata,
                                                  spot_composition_mtrx = decon_mtrx)

  return(list(lda_mod, test_spots_ls, raw_statistics_ls))
}
