#' This function takes in a seurat object with several tuning parameters and it assesses its performance on synthetic test spots
#'
#' @param se_sc Object of class Seurat with the scRNAseq data.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers().
#' @param n_syn_mixt Object of class integer specifying how many synthetic mixtures to generate.
#' @param cl_n Object of integer indicating how many cells to keep from each cluster. If a cluster has n < cl_n then all cells will be selected, if it has more then cl_n will be sampled randomly, 100 by default.
#' @param hvg Object of class numeric or "uns"; Number of highly variable genes to use on top of the marker genes, if "uns" then it is completely unsupervised and use top 3000 HVG.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), raw (no transformation applied). By default CPM.
#' @param method Object of class character; Type of method to us to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.01, by default, means that we will accept those cell types whose weight coefficient is at least 1\% of the total.
#' @return This function returns a list where the first element is a list with the NMF model trained and the cell labels, the second is a list with the raw_statistics.
#' @export
#' @examples
#'

spatial_decon_syn_assessment_nmf_fun <- function(se_sc,
                                                 clust_vr,
                                                 cluster_markers,
                                                 n_syn_mixt = 1000,
                                                 cl_n = 100,
                                                 hvg = 3000,
                                                 ntop = NULL,
                                                 transf = "uv",
                                                 method = "nsNMF",
                                                 min_cont = 0.01,
                                                 assay = "RNA",
                                                 slot = "counts") {

  ##########################
  ### Generate test data ###
  ##########################
  print(sprintf("Generating %s synthetic test mixtures", n_syn_mixt))
  test_spot_ls <- test_spot_fun(se_obj = se_sc,
                                clust_vr = clust_vr,
                                n = n_syn_mixt)

  test_spot_counts <- as.matrix(test_spot_ls[[1]])
  colnames(test_spot_counts) <- paste("mixt", 1:ncol(test_spot_counts), sep = "_")
  test_spot_metadata <- test_spot_ls[[2]]

  ####################
  ### Downsampling ###
  ####################
  # Downsample number of genes and number of samples
  print("Downsampling genes and cells")

  # Downsample scRNAseq to select gene set and number of cells to train the model
  se_sc_down <- downsample_se_obj(se_obj = se_sc,
                                  clust_vr = clust_vr,
                                  cluster_markers = cluster_markers,
                                  cl_n = cl_n,
                                  hvg = hvg)

  ###################
  #### Train NMF ####
  ###################
  print("Train NMF")
  # Train the NMF model
  nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers,
                          se_sc = se_sc_down,
                          mtrx_spatial = test_spot_counts,
                          ntop = ntop,
                          transf = transf,
                          clust_vr = clust_vr,
                          method = method,
                          assay = assay,
                          slot = slot)

  #################################
  #### Get mixture composition ####
  #################################
  # Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit.
  # Get cell type specific topif profiles
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod_ls[[1]]),
                                                     train_cell_clust = nmf_mod_ls[[2]])
  print("Deconvolute synthetic spots")
  # Perform deconvolution of the capture location mixtures
  pred_comp <- mixture_deconvolution_nmf(nmf_mod = nmf_mod_ls[[1]],
                                         mixture_transcriptome = test_spot_counts,
                                         transf = transf,
                                         reference_profiles = ct_topic_profiles,
                                         min_cont = min_cont)

  ################################
  #### Performance statistics ####
  ################################
  ct_cols <- colnames(pred_comp)[which(colnames(pred_comp) != "res_ss")]
  raw_statistics_ls <- test_synthetic_performance(test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
                                                  spot_composition_mtrx = pred_comp[, ct_cols])

  return(list("nmf_mod" = nmf_mod_ls, "stats" = raw_statistics_ls))
}
