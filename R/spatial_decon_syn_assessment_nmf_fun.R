#' This function takes in a seurat object with several tuning parameters and it assesses its performance on synthetic test spots
#'
#' @param se_obj Object of class Seurat.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param verbose Object of class Logical determining if progress should be reported or not (TRUE by default).
#' @param cl_n Object of class integer, how many cells to grab from each cluster, by default 10.
#' @param hvg Object of class integer, how many HVG to pass to the LDA model besides the cluster markers, by default 0
#' @param ntop Object of class "numeric"; number of unique markers per cluster used to seed the model, by default 100. If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: cpm (Counts per million), uv (unit variance), sct (Seurat::SCTransform), NULL (no transformation applied). By default CPM.
#' @param method Object of class character; Type of method to us to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @param logFC Object of class numeric; logFC to filter marker genes by, by default 1.
#' @param pct1 Object of class numeric; pct1 to filter marker genes by, by default 0.9.
#' @return This function returns a list where the first element is the lda model trained, the second is a list with test spot counts + metadata and the third element are the raw_statistics.
#' @export
#' @examples
#'

spatial_decon_syn_assessment_nmf_fun <- function(se_obj,
                                             clust_vr,
                                             verbose = TRUE,
                                             cl_n = 10,
                                             hvg = 0,
                                             ntop = 100,
                                             transf = "uv",
                                             method = "nsNMF",
                                             logFC = 1,
                                             pct1 = 0.9) {

  # Check variables
  if (is(se_obj) != "Seurat") stop("ERROR: se_obj must be a Seurat object!")
  if (!is.character(clust_vr)) stop("ERROR: clust_vr must be a character string!")
  if (!is.logical(verbose))stop("ERROR: verbose must be a logical object!")
  if (!is.numeric(cl_n)) stop("ERROR: ntop must be of class numeric!")
  if (!(is.numeric(hvg) | hvg == "uns")) stop("ERROR: hvg must be an object of class integer or 'uns'!")
  if (!is.numeric(ntop)) stop("ERROR: ntop must be of class numeric!")
  if (!is.character(transf)) stop("ERROR: ntop must be of class character!")
  if (!is.character(method)) stop("ERROR: ntop must be of class character!")

  ##########################
  ### Generate test data ###
  ##########################
  test_spot_ls <- test_spot_fun(se_obj = se_obj,
                                clust_vr = clust_vr,
                                n = 1000)

  test_spot_counts <- test_spot_ls[[1]]
  colnames(test_spot_counts) <- paste("mixt", 1:ncol(test_spot_counts), sep = "_")
  test_spot_metadata <- test_spot_ls[[2]]

  ############################
  #### Preprocessing data ####
  ############################
  ### Marker genes
  #### Extract the top marker genes from each cluster ####
  Seurat::Idents(object = se_obj) <- se_obj@meta.data[, clust_vr]
  # cluster_markers_all <- Seurat::FindAllMarkers(object = se_obj,
  #                                               verbose = TRUE,
  #                                               only.pos = TRUE,
  #                                               assay = "SCT",
  #                                               slot = "data",
  #                                               min.pct = 0.9,
  #                                               max.cells.per.ident = 100)

  # saveRDS(object = cluster_markers_all,
  #         file = sprintf("%s/%s/cluster_markers_%s.RDS",
  #                        an_mouse, robj_dir, id_comp))
  print("Loading markers")
  cluster_markers_all <- readRDS(file = sprintf("%s/%s/cluster_markers_%s.RDS",
                                                an_mouse, robj_dir, id_comp))

  cluster_markers_filt <- cluster_markers_all %>%
    filter(avg_logFC > logFC & pct.1 > pct1)

  ####################
  ### Downsampling ###
  ####################
  # Downsample number of genes and number of samples
  print("Downsampling genes and cells")
  se_obj_down <- downsample_se_obj(se_obj = se_obj,
                                   clust_vr = clust_vr,
                                   cluster_markers = cluster_markers_filt,
                                   cl_n = cl_n,
                                   hvg = hvg)
  ###################
  #### Train NMF ####
  ###################
  print("Train NMF")
  nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers_filt,
                          se_sc = se_obj_down,
                          mtrx_spatial = test_spot_counts,
                          ntop = ntop,
                          transf = transf,
                          clust_vr = clust_vr,
                          method = method)
  nmf_mod <- nmf_mod_ls[[1]]

  # get matrix H
  h <- coef(nmf_mod)

  #################################
  #### Get mixture composition ####
  #################################
  # Run test spots through the basis to get the pertinent coefficients. To do this for every spot we are going to set up a system of linear equations where we need to find the coefficient, we will use non-negative least squares to determine the best coefficient fit.
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = h,
                                                     train_cell_clust = nmf_mod_ls[[2]],
                                                     clust_vr = clust_vr)
  print("Deconvolute synthetic spots")
  pred_comp <- mixture_deconvolution_nmf(nmf_mod = nmf_mod,
                                         mixture_transcriptome = test_spot_counts,
                                         transf = transf,
                                         reference_profiles = ct_topic_profiles)

  ################################
  #### Performance statistics ####
  ################################
  ct_cols <- colnames(pred_comp)[which(colnames(pred_comp) != "res_ss")]
  raw_statistics_ls <- test_synthetic_performance(test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
                                                  spot_composition_mtrx = pred_comp[, ct_cols])

  return(list("nmf_mod" = nmf_mod_ls, "stats" = raw_statistics_ls))
}
