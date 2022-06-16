#' @name runDeconvolution
#' @rdname runDeconvolution
#' @title Run Deconvolution using NNLS model
#'
#' @aliases runDeconvolution
#'
#' @description This function takes in the mixture data, the trained model & the
#'   topic profiles and returns the proportion of each cell type within each
#'    mixture
#'  
#' @param x mixture dataset. Can be a numeric matrix,
#'   \code{SingleCellExperiment}, \code{SpatialExperiment} or 
#'   \code{SeuratObjecy}.
#' @param mod object of class NMFfit as obtained from trainNMF. 
#' @param ref bject of class matrix containing the topic profiles for each cell
#'  type as obtained from trainNMF. 
#' @param scale logical specifying whether to scale single-cell counts to unit
#'   variance. This gives the user the option to normalize the data beforehand
#'   as you see fit (CPM, FPKM, ...) when passing a matrix or specifying the
#'   slot from where to extract the count data.
#' @param min_prop scalar in [0,1] setting the minimum contribution
#'   expected from a cell type in \code{x} to observations in \code{y}.
#'   By default 0.
#' @param assay if the object is of Class \code{Seurat}, character string
#'   specifying the assay from which to extract the expression matrix.
#'   By default "RNA". Ignore for the rest of x input classes.
#' @param slot if the object is of Class \code{Seurat}, character string
#'   specifying the slot from which to extract the expression matrix. If the
#'   object is of class \code{SpatialExperiment} indicates matrix to use.
#'   By default "counts".
#' @param L1_nnls L1/LASSO penalty to be subtracted from b. See ?RcppML::nnls()
#'   for more info.
#' @param L2_nnls Ridge penalty to be added to diagonal of a. See ?RcppML::nmf()
#'   for more info.
#' @param verbose logical. Should information on progress be reported?
#'
#'
#' @return base a list where the first element is an \code{NMFfit} object and
#'   the second is a matrix contatining the topic profiles learnt.
#'
#' @author Marc Elosua Bayes & Helena L Crowell
#'
#' @examples
#' set.seed(321)
#' # mock up some single-cell, mixture & marker data
#' sce <- mockSC(ng = 200, nc = 10, nt = 3)
#' spe <- mockSP(sce)
#' mgs <- getMGS(sce)
#' 
#' res <- trainNMF(
#'     x = sce,
#'     y = spe,
#'     groups = sce$type,
#'     mgs = mgs,
#'     weight_id = "weight",
#'     group_id = "type",
#'     gene_id = "gene")
#' # Run deconvolution
#' decon <- runDeconvolution(
#'     x = spe,
#'     mod = res[["mod"]],
#'     ref = res[["topic"]])
NULL

#' @rdname runDeconvolution
#' @importFrom RcppML predict.nmf
#' @export
runDeconvolution <- function(
    x,
    mod,
    ref,
    scale = TRUE,
    min_prop = 0.01,
    verbose = TRUE,
    assay = "RNA",
    slot = "counts",
    L1_nnls = 0,
    L2_nnls = 0) {
    
    # Class checks
    stopifnot(
        # Check x inputs
        is.matrix(x) | is(x, "DelayedMatrix") | is(x, "dgCMatrix") |
            is(x, "Seurat") | is(x, "SingleCellExperiment") |
            is(x, "SpatialExperiment"),
        # Check mod inputs
        is.list(mod),
        # check ref
        is.matrix(ref),
        # Check assay name
        is.character(assay), length(assay) == 1,
        # Check slot name
        is.character(slot), length(slot) == 1,
        # Check scale and verbose
        is.logical(scale), length(scale) == 1,
        is.logical(verbose), length(verbose) == 1,
        # Check min_prop numeric
        is.numeric(min_prop), length(min_prop) == 1,
        min_prop >= 0, min_prop <= 1
    )
    
    # Extract expression matrix
    if (!is.matrix(x))
        x <- .extract_counts(x, assay, slot)
    
    # Get topic profiles for mixtures
    mat <- .pred_prop(x, mod, scale)
    
    if (verbose) message("Deconvoluting mixture data...")
    if (verbose) message(paste0("Using L1 reg: ", L1_nnls))
    if (verbose) message(paste0("Using L2 reg: ", L2_nnls))
    ref_scale <- t(t(ref) / colSums(ref))
    pred <- RcppML::predict.nmf(w = ref_scale, data = mat, L1 = L1_nnls, L2 = L2_nnls)
    rownames(pred) <- rownames(ref_scale)
    
    # Proportions within each cell type
    res <- prop.table(pred, 2)
    
    # TODO Compute residuals
    ss <- colSums(mat^2)
    err <- rep(0, ncol(res))
    names(err) <- colnames(res)
    
    # res <- vapply(seq_len(ncol(mat)), function(i) {
    #     pred <- nnls::nnls(ref, mat[, i])
    #     # scale reference
    #     # ref_scale <- t(ref) / colSums(ref)
    #     # prop <- prop.table(pred$x)
    #     # drop groups that fall below 'min_prop' & update
    #     # prop[prop < min_prop] <- 0
    #     # prop <- prop.table(prop)
    #     # compute residual sum of squares
    #     ss <- sum(mat[, i]^2)
    #     # compute percentage of unexplained residuals
    #     err <- pred$deviance / ss
    #     c(prop, err)
    # }, numeric(ncol(ref) + 1))
    
    # set dimension names
    # rownames come from the reference
    # rownames(res) <- c(rownames(ref), "res_ss")
    # colnames(res) <- colnames(mat)
    
    # Separate residuals from proportions
    # Extract residuals
    # err <- res["res_ss", ]
    # Extract only deconvolution matrices
    # res <- res[-nrow(res), ]
    
    return(list("mat" = t(res), "res_ss" = err))
}

