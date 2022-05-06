#' @name SPOTlight
#' @title Deconvolution of mixture using single-cell data
#'
#' @description This is the backbone function which takes in single cell
#'   expression data to deconvolute spatial transcriptomics spots.
#'
#' @param x,y single-cell and mixture dataset, respectively. Can be a
#'   numeric matrix, \code{SingleCellExperiment} or \code{SeuratObjecy}.
#' @param groups vector of group labels for cells in \code{x}.
#'   When \code{x} is a \code{SingleCellExperiment} or \code{SeuratObject},
#'   defaults to \code{colLabels} and \code{Idents(x)}, respectively.
#' @param mgs \code{data.frame} or \code{DataFrame} of marker genes.
#'   Must contain columns holding gene identifiers, group labels and
#'   the weight (e.g., logFC, -log(p-value) a feature has in a given group.
#' @param pnmf character vector specifying which from which package to grab
#'   \code{nmf()}. It can be either \code{RcppML} (default) or \code{NMF}.
#' @param hvg character vector containing hvg to include in the model.
#'   By default NULL.
#' @param gene_id,group_id,weight_id character specifying the column
#'   in \code{mgs} containing gene identifiers, group labels and weights,
#'   respectively.
#' @param scale logical specifying whether to scale single-cell counts to unit
#'   variance. This gives the user the option to normalize the data beforehand
#'   as you see fit (CPM, FPKM, ...) when passing a matrix or specifying the
#'   slot from where to extract the count data.
#' @param min_prop scalar in [0,1] setting the minimum contribution
#'   expected from a cell type in \code{x} to observations in \code{y}.
#'   By default 0.
#' @param assay if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'   By default "RNA".
#' @param slot if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SingleCellExperiment} indicates matrix to use.
#'   By default "counts".
#' @param n_top integer scalar specifying the number of markers to select per
#'  group. By default NULL uses all the marker genes to initialize the model.
#' @param model character string indicating which model to use when running NMF.
#' Either "ns" (default) or "std".
#' @param L1 LASSO penalties in the range (0, 1], single value or array of
#'   length two for c(w, h). See ?RcppML::nmf() for more info.
#' @param L2 Ridge penalties greater than zero, single value or array of length
#'   two for c(w, h). See ?RcppML::nmf() for more info.
#' @param tol tolerance of the fit ?RcppML::nmf() for more info.
#' @param verbose logical. Should information on progress be reported?
#' @param ... additional parameters.
#'
#' @return a numeric matrix with rows corresponding to samples
#'   and columns to groups
#'
#' @author Marc Elosua-Bayes & Helena L. Crowell
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
SPOTlight <- function(
    x,
    y,
    groups = NULL,
    # markers
    mgs,
    pnmf = c("RcppML", "NMF"),
    n_top = NULL,
    gene_id = "gene",
    group_id = "cluster",
    weight_id = "weight",
    hvg = NULL,
    # NMF
    scale = TRUE,
    model = c("ns", "std"),
    # deconvolution
    min_prop = 0.01,
    # other
    verbose = TRUE,
    assay = "RNA",
    slot = "counts",
    L1 = 0.5,
    L2 = 0,
    tol = 1e-5,
    ...) {
    
    # train NMF model
    mod_ls <- trainNMF(x, y, groups, mgs, pnmf, n_top, gene_id, group_id,
        weight_id, hvg, model, scale, verbose, assay, slot, L1, L2, tol...)
    
    # perform deconvolution
    res <- runDeconvolution(y, mod_ls[["mod"]], mod_ls[["topic"]],
        scale, min_prop, verbose)

    # return list of NMF model & deconvolution matrix
    list(
        "mat" = res[["mat"]],
        "res_ss" = res[["res_ss"]],
        "NMF" = mod_ls[["mod"]])
    }
