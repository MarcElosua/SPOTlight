#' This functions takes in a Seurat object with the scRNAseq data and the count matrix of the spatial transcriptomics data and returns the deconvoluted spots.
#'
#' @param se_sc Object of class Seurat with the scRNAseq data.
#' @param counts_spatial Object of class Matrix of shape GENESxSPOT, se_obj\@assays$Spatial@counts.
#' @param clust_vr Object of class character. Name of the variable containing the cell clustering.
#' @param cluster_markers Object of class dataframe obtained from the function Seurat::FindAllMarkers()
#' @param cl_n Object of class integer indicating how many cells to keep from each cluster. If a cluster has n < cl_n then all cells will be selected, if it has more then cl_n will be sampled randomly, 100 by default.
#' @param hvg Object of class numeric or NULL; Number of highly variable genes to use on top of the marker genes, if NULL then it is completely unsupervised and use top 3000 HVG.
#' @param ntop Object of class "numeric" or NULL; number of unique markers per cluster used to seed the model, by default NULL. If NULL it uses all of them.
#' @param transf Transformation to normalize the count matrix: uv (unit variance), raw (no transformation applied). By default UV.
#' @param method Object of class character; Type of method to us to find W and H. Look at NMF package for the options and specifications, by default nsNMF.
#' @param min_cont Object of class numeric; Indicates the minimum contribution we expect from a cell in that spot. Since we're working with proportions by setting 0.01, by default, means that we will accept those cell types whose weight coefficient is at least 1\% of the total.
#' @param assay Object of class character; From which assay to grab the expression data to train the model, by default "RNA".
#' @param slot Object of class character; From which slot to grab the expression data to train the model, by default "counts".
#' @return This function returns a matrix with the coefficients of the spatial mixtures.
#' @export
#' @examples
#'


spotlight_deconvolution <- function(se_sc,
                                    counts_spatial,
                                    clust_vr,
                                    cluster_markers,
                                    cl_n = 100,
                                    hvg = 3000,
                                    ntop = NULL,
                                    transf = "uv",
                                    method = "nsNMF",
                                    min_cont = 0.01,
                                    assay = "RNA",
                                    slot = "counts") {

  # Downsample scRNAseq to select gene set and number of cells to train the model
  se_sc_down <- downsample_se_obj(se_obj = se_sc,
                                  clust_vr = clust_vr,
                                  cluster_markers = cluster_markers,
                                  cl_n = cl_n,
                                  hvg = hvg)

  # Train the NMF model
  nmf_mod_ls <- train_nmf(cluster_markers = cluster_markers,
                          se_sc = se_sc_down,
                          mtrx_spatial = counts_spatial,
                          ntop = ntop,
                          transf = transf,
                          clust_vr = clust_vr,
                          method = method,
                          assay = assay,
                          slot = slot)

  # Get cell type specific topif profiles
  ct_topic_profiles <- topic_profile_per_cluster_nmf(h = coef(nmf_mod_ls[[1]]),
                                                     train_cell_clust = nmf_mod_ls[[2]])

  # Perform deconvolution of the capture location mixtures
  decon_mtrx <- mixture_deconvolution_nmf(nmf_mod = nmf_mod_ls[[1]],
                                          mixture_transcriptome = counts_spatial,
                                          transf = transf,
                                          reference_profiles = ct_topic_profiles,
                                          min_cont = min_cont)

  return(list(nmf_mod_ls, decon_mtrx))
}
