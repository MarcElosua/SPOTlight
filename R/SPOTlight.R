#' @name SPOTlight
#' @title Deconvolution of mixture using single-cell data
#'
#' @description This is the backbone function which takes in single cell
#'   expression data to deconvolute spatial transcriptomics spots.
#'
#' @param x,y single-cell and mixture dataset, respectively. Can be a
#'   numeric matrix, \code{SingleCellExperiment} or \code{SeuratObjecy}.
#' @param groups character vector of group labels for cells in \code{x}.
#'   When \code{x} is a \code{SingleCellExperiment} or \code{SeuratObject},
#'   defaults to \code{colLabels(x)} and \code{Idents(x)}, respectively.
#'   Make sure groups is not a Factor.
#' @param mgs \code{data.frame} or \code{DataFrame} of marker genes.
#'   Must contain columns holding gene identifiers, group labels and
#'   the weight (e.g., logFC, -log(p-value) a feature has in a given group.
#' @param hvg character vector containing hvg to include in the model.
#'   By default NULL.
#' @param gene_id,group_id,weight_id character specifying the column
#'   in \code{mgs} containing gene identifiers, group labels and weights,
#'   respectively.
#' @param scale logical specifying whether to scale single-cell counts to unit
#'   variance. This gives the user the option to normalize the data beforehand
#'   as you see fit (CPM, FPKM, ...) when passing a matrix or specifying the
#'   slot from where to extract the count data.
#' @param n_top integer scalar specifying the number of markers to select per
#'  group. By default NULL uses all the marker genes to initialize the model.
#' @param model character string indicating which model to use when running NMF.
#' Either "ns" (default) or "std".
#' @param assay_sc,assay_sp if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'   By default "RNA" and "Spatial".
#' @param slot_sc,slot_sp if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SingleCellExperiment} indicates matrix to use.
#'   By default "counts".
#' @param L1_nmf LASSO penalty in the range (0, 1] for NMF,
#'   larger values increase sparsity of each factor
#' @param L2_nmf RUDGE penalty >0 for NMF,
#'  larger values increase angle between factors and thus sparsity.
#' @param tol tolerance of the NMF model at convergence, the Pearson correlation 
#'   distance between models across consecutive iterations (1e-5 is publication quality)
#' @param maxit maximum number of NMF iterations for fitting
#' @param threads number of threads to use, default 0 (all threads)
#' @param verbose logical. Should information on progress be reported?
#' @param min_prop scalar in [0,1] setting the minimum contribution
#'   expected from a cell type in \code{x} to observations in \code{y}.
#'   By default 0.
#' @param L1_nnls_topics,L1_nnls_prop LASSO penalty in the range (0, 1] for NNLS
#'   when computing cell type topic profiles and cell type proportions
#'   respectively. Larger values remove "noisy" contributions more aggressively.
#' @param L2_nnls_topics,L2_nnls_prop RIDGE penalty >0 for NNLS when computing
#'   cell type topic profiles and cell type proportions respectively.
#'   Larger values remove "noisy" contributions more aggressively.
#' @param ... additional parameters.
#'
#' @return a numeric matrix with rows corresponding to samples
#'   and columns to groups
#'
#' @author Marc Elosua Bayes, Zach DeBruine, and Helena L Crowell
#'
#' @details SPOTlight uses a Non-Negative Matrix Factorization approach to learn
#'   which genes are important for each cell type. In order to drive the
#'   factorization and give more importance to cell type marker genes we
#'   previously compute them and use them to initialize the basis matrix. This
#'   initialized matrices will then be used to carry out the factorization with
#'   the single cell expression data. Once the model has learn the topic
#'   profiles for each cell type we use non-negative least squares (NNLS) to
#'   obtain the topic contributions to each spot. Lastly, NNLS is again used to
#'   obtain the proportion of each cell type for each spot by finding the
#'   fitting the single-cell topic profiles to the spots topic contributions.
#'
#' @examples
#' library(scater)
#' library(scran)
#' 
#' # Use Mock data
#' # Refer to the vignette for a full workflow
#' sce <- mockSC(ng = 200, nc = 10, nt = 3)
#' spe <- mockSP(sce)
#' mgs <- getMGS(sce)
#' 
#' res <- SPOTlight(
#'     x = counts(sce),
#'     y = counts(spe),
#'     groups = as.character(sce$type),
#'     mgs = mgs,
#'     hvg = NULL,
#'     weight_id = "weight",
#'     group_id = "type",
#'     gene_id = "gene")
NULL

#' @rdname SPOTlight
#' @export
#' @useDynLib SPOTlight, .registration = TRUE
#' 
SPOTlight <- function(
    x,
    y,
    groups = NULL,
    mgs,
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    hvg = NULL,
    scale = TRUE,
    model = c("ns", "std"),
    min_prop = 0.01,
    verbose = TRUE,
    assay_sc = "RNA",
    slot_sc = "counts",
    assay_sp = "Spatial",
    slot_sp = "counts",
    L1_nmf = 0,
    L2_nmf = 0,
    maxit = 100,
    threads = 0,
    tol = 1e-5,
    L1_nnls_topics = 0,
    L2_nnls_topics = 0,
    L1_nnls_prop = 0,
    L2_nnls_prop = 0,
  ...) {
    
    # train NMF model
    mod_ls <- trainNMF(
        x = x,
        y = rownames(y),
        groups = groups,
        mgs = mgs,
        n_top = n_top,
        gene_id = gene_id, 
        group_id = group_id, 
        weight_id = weight_id, 
        hvg = hvg, 
        model = model, 
        scale = scale,
        verbose = verbose,
        assay_sc = assay_sc,
        slot_sc = slot_sc,
        L1_nmf = L1_nmf,
        L2_nmf = L2_nmf,
        tol = tol,
        threads = threads,
        maxit = maxit,
        ...)
    
    # perform deconvolution
    res <- runDeconvolution(
        x = y,
        mod = mod_ls[["mod"]], 
        ref = mod_ls[["topic"]], 
        scale = scale, 
        min_prop = min_prop, 
        verbose = verbose,
        assay = assay_sp,
        slot = slot_sp,
        L1_nnls_topics = L1_nnls_topics,
        L2_nnls_topics = L2_nnls_topics,
        L1_nnls_prop = L1_nnls_prop,
        L2_nnls_prop = L2_nnls_prop)
    
    # return list of NMF model & deconvolution matrix
    list(
        "mat" = res[["mat"]],
        "res_ss" = res[["res_ss"]],
        "NMF" = mod_ls[["mod"]])
    
}
